"""
04_ai_biomarker_engine.py
AI 기반 바이오마커 예측 엔진

전략:
  1. 전이학습 (pan-disease → periodontitis 특화)
  2. 그래프 신경망 (GNN) — PPI + 오믹스 통합
  3. 메타분석 기반 앙상블 랭킹
  4. 검증 프레임워크 (내부 + 외부 + 범질환)
"""

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

CURATED = Path("../data/curated")
OUTDIR  = Path("../results/ai_engine")
OUTDIR.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
# 1. 데이터셋 클래스
# ═══════════════════════════════════════════════════════════════
class MultiOmicsDataset(Dataset):
    """PyTorch Dataset — 멀티오믹스 데이터 로더"""

    def __init__(self,
                 mrna: np.ndarray,
                 prot: np.ndarray,
                 metab: np.ndarray,
                 labels: np.ndarray):
        self.mrna   = torch.FloatTensor(mrna)
        self.prot   = torch.FloatTensor(prot)
        self.metab  = torch.FloatTensor(metab)
        self.labels = torch.LongTensor(labels)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return {
            'mrna':   self.mrna[idx],
            'prot':   self.prot[idx],
            'metab':  self.metab[idx],
            'label':  self.labels[idx]
        }


# ═══════════════════════════════════════════════════════════════
# 2. 전이학습 모델 (Pan-disease → Periodontitis)
# ═══════════════════════════════════════════════════════════════
class TransferLearningModel(nn.Module):
    """
    전이학습 구조:
    Step 1: TCGA/GEO 다양한 염증 질환 데이터로 사전학습
            (범질환 오믹스 패턴 학습)
    Step 2: 치주염 데이터로 fine-tuning
            (질환 특이적 바이오마커 특화)
    """

    def __init__(self,
                 mrna_dim:  int = 2000,
                 prot_dim:  int = 1000,
                 metab_dim: int = 500,
                 hidden_dim: int = 256,
                 latent_dim: int = 64,
                 n_classes:  int = 2):
        super().__init__()

        # 오믹스별 인코더 (공유 구조)
        self.mrna_encoder = nn.Sequential(
            nn.Linear(mrna_dim, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, hidden_dim),
            nn.ReLU()
        )

        self.prot_encoder = nn.Sequential(
            nn.Linear(prot_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, hidden_dim),
            nn.ReLU()
        )

        self.metab_encoder = nn.Sequential(
            nn.Linear(metab_dim, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, hidden_dim),
            nn.ReLU()
        )

        # Attention 기반 오믹스 통합 (각 오믹스 가중치 자동 학습)
        self.attention = nn.MultiheadAttention(
            embed_dim=hidden_dim,
            num_heads=4,
            dropout=0.1,
            batch_first=True
        )

        # 통합 표현 → 잠재 공간
        self.fusion = nn.Sequential(
            nn.Linear(hidden_dim * 3, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, latent_dim),
            nn.ReLU()
        )

        # 분류 헤드 (fine-tuning 시 교체)
        self.classifier = nn.Linear(latent_dim, n_classes)

        # 오믹스별 기여도 추적
        self.omics_weights = nn.Parameter(torch.ones(3) / 3)

    def forward(self, mrna, prot, metab):
        # 각 오믹스 인코딩
        z_mrna  = self.mrna_encoder(mrna)
        z_prot  = self.prot_encoder(prot)
        z_metab = self.metab_encoder(metab)

        # Attention 기반 통합
        # (batch, 3 views, hidden_dim)
        omics_stack = torch.stack([z_mrna, z_prot, z_metab], dim=1)
        attended, weights = self.attention(
            omics_stack, omics_stack, omics_stack
        )

        # Flatten + Fusion
        fused = attended.reshape(attended.size(0), -1)
        latent = self.fusion(fused)

        # 분류
        logits = self.classifier(latent)
        return logits, latent, weights

    def freeze_encoders(self):
        """Fine-tuning 시 인코더 동결"""
        for param in self.mrna_encoder.parameters():
            param.requires_grad = False
        for param in self.prot_encoder.parameters():
            param.requires_grad = False
        for param in self.metab_encoder.parameters():
            param.requires_grad = False

    def unfreeze_all(self):
        for param in self.parameters():
            param.requires_grad = True


# ═══════════════════════════════════════════════════════════════
# 3. 그래프 신경망 (PPI + 오믹스 통합)
# ═══════════════════════════════════════════════════════════════
class GNNBiomarkerModel(nn.Module):
    """
    Graph Neural Network — PPI 구조를 활용한 바이오마커 발굴
    노드: 유전자/단백질
    엣지: STRING PPI 상호작용
    노드 특성: 발현값 (mRNA, Protein 동시)
    """

    def __init__(self,
                 node_feat_dim: int = 2,   # [log2FC_mrna, log2FC_prot]
                 hidden_dim: int = 64,
                 n_classes: int = 2):
        super().__init__()

        try:
            from torch_geometric.nn import GCNConv, GATConv, global_mean_pool

            # GAT (Graph Attention Network) — 노드별 중요도 학습
            self.conv1 = GATConv(node_feat_dim, hidden_dim, heads=4, dropout=0.2)
            self.conv2 = GATConv(hidden_dim * 4, hidden_dim, heads=1, dropout=0.2)
            self.conv3 = GATConv(hidden_dim, hidden_dim // 2, heads=1)

            self.classifier = nn.Sequential(
                nn.Linear(hidden_dim // 2, 32),
                nn.ReLU(),
                nn.Linear(32, n_classes)
            )
            self.gnn_available = True

        except ImportError:
            print("torch_geometric 미설치: pip install torch_geometric")
            self.gnn_available = False
            # 대체: 단순 MLP
            self.classifier = nn.Sequential(
                nn.Linear(node_feat_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, n_classes)
            )

    def forward(self, x, edge_index=None, batch=None):
        if self.gnn_available and edge_index is not None:
            from torch_geometric.nn import global_mean_pool
            x = F.elu(self.conv1(x, edge_index))
            x = F.dropout(x, p=0.2, training=self.training)
            x = F.elu(self.conv2(x, edge_index))
            x = self.conv3(x, edge_index)
            if batch is not None:
                x = global_mean_pool(x, batch)
        return self.classifier(x)


# ═══════════════════════════════════════════════════════════════
# 4. 학습 + 검증 프레임워크
# ═══════════════════════════════════════════════════════════════
class BiomarkerTrainer:

    def __init__(self, model: nn.Module,
                 device: str = "cuda" if torch.cuda.is_available() else "cpu"):
        self.model  = model.to(device)
        self.device = device
        print(f"  디바이스: {device}")

    def pretrain(self, pan_disease_loader: DataLoader,
                  epochs: int = 50) -> list:
        """
        Pan-disease 사전학습 (TCGA, GEO 다양한 염증 데이터)
        목표: 범질환 오믹스 패턴 학습
        """
        optimizer = torch.optim.Adam(self.model.parameters(),
                                      lr=1e-3, weight_decay=1e-4)
        criterion = nn.CrossEntropyLoss()
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=epochs
        )

        losses = []
        for epoch in range(epochs):
            self.model.train()
            epoch_loss = 0
            for batch in pan_disease_loader:
                mrna   = batch['mrna'].to(self.device)
                prot   = batch['prot'].to(self.device)
                metab  = batch['metab'].to(self.device)
                labels = batch['label'].to(self.device)

                optimizer.zero_grad()
                logits, _, _ = self.model(mrna, prot, metab)
                loss = criterion(logits, labels)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()

            scheduler.step()
            losses.append(epoch_loss / len(pan_disease_loader))
            if epoch % 10 == 0:
                print(f"  Epoch {epoch}: loss={losses[-1]:.4f}")

        return losses

    def fine_tune(self, periodontitis_dataset: MultiOmicsDataset,
                   epochs: int = 100,
                   freeze_epochs: int = 20) -> dict:
        """
        치주염 데이터 Fine-tuning
        Phase 1: 인코더 동결 → 분류 헤드만 학습
        Phase 2: 전체 모델 학습
        """
        self.model.freeze_encoders()
        optimizer = torch.optim.Adam(
            filter(lambda p: p.requires_grad, self.model.parameters()),
            lr=1e-3
        )
        criterion = nn.CrossEntropyLoss()

        # LOOCV fine-tuning (n=60 소규모)
        n = len(periodontitis_dataset)
        all_probs  = np.zeros(n)
        all_labels = np.zeros(n)

        for test_idx in range(n):
            train_indices = [i for i in range(n) if i != test_idx]
            train_subset = torch.utils.data.Subset(
                periodontitis_dataset, train_indices
            )
            train_loader = DataLoader(train_subset, batch_size=16, shuffle=True)
            test_loader  = DataLoader(
                torch.utils.data.Subset(periodontitis_dataset, [test_idx]),
                batch_size=1
            )

            # 모델 리셋
            self.model.unfreeze_all()
            self.model.freeze_encoders()

            for phase_epoch in range(epochs):
                if phase_epoch == freeze_epochs:
                    self.model.unfreeze_all()
                    optimizer = torch.optim.Adam(
                        self.model.parameters(), lr=1e-4
                    )

                self.model.train()
                for batch in train_loader:
                    mrna   = batch['mrna'].to(self.device)
                    prot   = batch['prot'].to(self.device)
                    metab  = batch['metab'].to(self.device)
                    labels = batch['label'].to(self.device)
                    optimizer.zero_grad()
                    logits, _, _ = self.model(mrna, prot, metab)
                    loss = criterion(logits, labels)
                    loss.backward()
                    optimizer.step()

            # 테스트
            self.model.eval()
            with torch.no_grad():
                for batch in test_loader:
                    mrna   = batch['mrna'].to(self.device)
                    prot   = batch['prot'].to(self.device)
                    metab  = batch['metab'].to(self.device)
                    label  = batch['label'].item()
                    logits, _, _ = self.model(mrna, prot, metab)
                    prob = F.softmax(logits, dim=1)[0, 1].item()
                    all_probs[test_idx]  = prob
                    all_labels[test_idx] = label

            if test_idx % 10 == 0:
                print(f"  LOOCV {test_idx+1}/{n}...")

        auc = roc_auc_score(all_labels, all_probs)
        print(f"\n  Fine-tuning LOOCV AUC: {auc:.3f}")
        return {"auc": auc, "probs": all_probs, "labels": all_labels}


# ═══════════════════════════════════════════════════════════════
# 5. 바이오마커 랭킹 엔진 (통합 점수)
# ═══════════════════════════════════════════════════════════════
class BiomarkerRankingEngine:
    """
    모든 증거를 통합한 최종 바이오마커 랭킹

    증거 원천:
      A. 자체 연구 (n=60)
      B. GEO 공개 코호트 (메타분석)
      C. PRIDE 공개 프로테오믹스
      D. MetaboLights 공개 대사체
      E. AI 예측 점수
      F. PPI 네트워크 중심성
      G. 문헌 빈도 (PubMed 언급 횟수)
    """

    def calculate_evidence_score(self, feature: str,
                                  evidence: dict) -> float:
        """
        통합 증거 점수 계산 (0~100)

        점수 구성:
          통계적 유의성    (30점): 메타분석 p-value 기반
          재현성          (25점): 유의한 연구 수
          오믹스 레이어   (20점): 몇 개 오믹스에서 확인
          네트워크 중심성 (10점): PPI hub 여부
          AI 예측         (10점): 딥러닝 기여도
          문헌 근거        (5점): PubMed 언급 빈도
        """
        score = 0.0

        # 1. 통계적 유의성 (30점)
        meta_p = evidence.get("meta_pvalue", 1.0)
        if meta_p < 0.001:
            score += 30
        elif meta_p < 0.01:
            score += 25
        elif meta_p < 0.05:
            score += 20
        elif meta_p < 0.1:
            score += 10

        # 2. 재현성 (25점)
        n_sig = evidence.get("n_studies_sig", 0)
        score += min(25, n_sig * 8)  # 최대 3개 연구 = 25점

        # 3. 오믹스 레이어 (20점)
        n_omics = evidence.get("n_omics_layers", 0)
        score += n_omics * 6.7  # 3개 오믹스 = 20점

        # 4. PPI Hub (10점)
        if evidence.get("is_ppi_hub", False):
            hub_score = evidence.get("hub_centrality", 0.5)
            score += hub_score * 10

        # 5. AI 예측 (10점)
        ai_score = evidence.get("ai_confidence", 0.5)
        score += ai_score * 10

        # 6. 문헌 (5점)
        pubmed_count = min(evidence.get("pubmed_count", 0), 10)
        score += pubmed_count * 0.5

        return min(100, score)

    def rank_all_candidates(self,
                             own_study: pd.DataFrame,
                             geo_meta: pd.DataFrame,
                             pride_meta: pd.DataFrame = None,
                             ppi_hubs: list = None,
                             ai_scores: dict = None) -> pd.DataFrame:
        """전체 바이오마커 후보 통합 랭킹"""

        all_features = set(own_study['feature_name'].tolist())
        if geo_meta is not None:
            all_features |= set(geo_meta['feature_name'].tolist())

        records = []
        for feat in all_features:
            evidence = {}

            # 자체 연구 데이터
            own = own_study[own_study['feature_name'] == feat]
            if not own.empty:
                evidence['meta_pvalue']    = own.iloc[0].get('padj', 1.0)
                evidence['n_studies_sig']  = 1
                evidence['n_omics_layers'] = own.iloc[0].get('n_omics', 1)

            # GEO 메타분석
            if geo_meta is not None:
                geo = geo_meta[geo_meta['feature_name'] == feat]
                if not geo.empty:
                    evidence['meta_pvalue'] = min(
                        evidence.get('meta_pvalue', 1.0),
                        geo.iloc[0].get('pooled_pvalue', 1.0)
                    )
                    evidence['n_studies_sig'] = evidence.get('n_studies_sig',0) + \
                        int(geo.iloc[0].get('n_studies', 0))

            # PPI Hub
            evidence['is_ppi_hub'] = feat in (ppi_hubs or [])

            # AI 예측
            if ai_scores:
                evidence['ai_confidence'] = ai_scores.get(feat, 0.5)

            score = self.calculate_evidence_score(feat, evidence)
            confidence = ("High" if score >= 70
                         else "Medium" if score >= 40
                         else "Low")

            records.append({
                'feature_name':   feat,
                'biomarker_score': score,
                'confidence':     confidence,
                **evidence
            })

        df = pd.DataFrame(records).sort_values('biomarker_score', ascending=False)
        df.to_csv(OUTDIR / "integrated_biomarker_ranking.csv", index=False)
        print(f"\n✅ 통합 랭킹: {len(df)}개 후보")
        print(f"  High confidence: {(df['confidence']=='High').sum()}개")
        print(f"  Medium confidence: {(df['confidence']=='Medium').sum()}개")
        return df


# ═══════════════════════════════════════════════════════════════
# 6. 검증 프레임워크
# ═══════════════════════════════════════════════════════════════
class ValidationFramework:
    """
    다층 검증:
    L1: 내부 LOOCV (자체 n=60)
    L2: GEO 외부 코호트
    L3: PRIDE/MetaboLights 독립 코호트
    L4: TCGA 범질환 검증 (염증 공통성)
    """

    def validate_on_geo_cohort(self,
                                panel_features: list,
                                geo_expr: pd.DataFrame,
                                geo_labels: pd.Series,
                                method: str = "logistic") -> dict:
        """GEO 공개 코호트 외부 검증"""
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import StratifiedKFold

        # 공통 피처만 사용
        common = [f for f in panel_features if f in geo_expr.index]
        if len(common) < 3:
            return {"error": "공통 피처 부족", "n_common": len(common)}

        X = geo_expr.loc[common].T.values
        y = geo_labels.values

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # 5-fold CV (외부 코호트는 n이 클 수 있음)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        aucs = []
        for tr, te in cv.split(X_scaled, y):
            clf = LogisticRegression(C=0.1, random_state=42, max_iter=1000)
            clf.fit(X_scaled[tr], y[tr])
            prob = clf.predict_proba(X_scaled[te])[:, 1]
            if len(np.unique(y[te])) > 1:
                aucs.append(roc_auc_score(y[te], prob))

        result = {
            "cohort":      "GEO_external",
            "n_samples":   len(y),
            "n_features":  len(common),
            "mean_auc":    np.mean(aucs),
            "std_auc":     np.std(aucs),
            "n_folds":     len(aucs)
        }
        print(f"  외부 검증 AUC: {result['mean_auc']:.3f} ± {result['std_auc']:.3f}")
        return result

    def cross_disease_validation(self,
                                  panel_features: list,
                                  disease_datasets: dict) -> pd.DataFrame:
        """
        범질환 검증 (IBD, RA, SLE 데이터에서 같은 패턴 확인)
        → 범염증 마커 vs 치주염 특이 마커 구분
        """
        results = []
        for disease, (expr, labels) in disease_datasets.items():
            common = [f for f in panel_features if f in expr.index]
            if len(common) < 3:
                continue

            X = expr.loc[common].T.values
            y = labels.values

            from sklearn.linear_model import LogisticRegression
            from sklearn.preprocessing import StandardScaler
            from sklearn.model_selection import StratifiedKFold

            scaler = StandardScaler()
            X_sc = scaler.fit_transform(X)

            clf = LogisticRegression(C=0.1, random_state=42, max_iter=1000)
            cv  = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            probs_all = np.zeros(len(y))

            for tr, te in cv.split(X_sc, y):
                clf.fit(X_sc[tr], y[tr])
                probs_all[te] = clf.predict_proba(X_sc[te])[:, 1]

            auc = roc_auc_score(y, probs_all) if len(np.unique(y)) > 1 else 0.5

            results.append({
                "disease":    disease,
                "auc":        auc,
                "n_features": len(common),
                "interpretation": (
                    "범염증 마커 (전이 가능)" if auc > 0.7
                    else "치주염 특이 마커" if auc < 0.6
                    else "부분 전이"
                )
            })

        df = pd.DataFrame(results)
        df.to_csv(OUTDIR / "cross_disease_validation.csv", index=False)
        return df


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("AI Biomarker Engine — 사용 예시")
    print()
    print("# 1. 전이학습 모델 생성")
    print("model = TransferLearningModel(mrna_dim=2000, prot_dim=500)")
    print()
    print("# 2. Pan-disease 사전학습 (TCGA/GEO 데이터)")
    print("trainer = BiomarkerTrainer(model)")
    print("trainer.pretrain(pan_disease_loader, epochs=50)")
    print()
    print("# 3. 치주염 fine-tuning")
    print("results = trainer.fine_tune(periodontitis_dataset, epochs=100)")
    print()
    print("# 4. 통합 랭킹")
    print("ranker = BiomarkerRankingEngine()")
    print("ranking = ranker.rank_all_candidates(own_study, geo_meta)")
    print()
    print("# 5. 외부 검증")
    print("validator = ValidationFramework()")
    print("ext_result = validator.validate_on_geo_cohort(panel, geo_expr, geo_labels)")
