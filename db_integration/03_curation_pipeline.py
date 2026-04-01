"""
03_curation_pipeline.py
공개 DB 데이터 큐레이션 파이프라인
표준화 → 품질관리 → 배치보정 → 정규화 → 통합 DB 적재
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional
import warnings
warnings.filterwarnings('ignore')

RAW_DIR    = Path("../data/public_db")
CURATED    = Path("../data/curated")
CURATED.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
# 1. ID 표준화
# ═══════════════════════════════════════════════════════════════
class IDHarmonizer:
    """
    다양한 식별자 → 표준 ID 변환
    Gene: Probe/Ensembl/Entrez → HGNC symbol
    Protein: PD accession → UniProt Swiss-Prot
    Metabolite: name → InChIKey / HMDB ID
    """

    def __init__(self):
        self._gene_map   = {}
        self._prot_map   = {}
        self._metab_map  = {}

    def harmonize_genes(self, df: pd.DataFrame,
                         id_col: str = "gene",
                         id_type: str = "symbol") -> pd.DataFrame:
        """
        유전자 ID 표준화 → HGNC symbol
        id_type: 'symbol', 'ensembl', 'entrez', 'probe'
        """
        try:
            import mygene
            mg = mygene.MyGeneInfo()

            gene_ids = df[id_col].dropna().unique().tolist()
            print(f"  유전자 ID 변환: {len(gene_ids)}개 ({id_type} → HGNC)")

            scope = {
                'symbol':  'symbol',
                'ensembl': 'ensembl.gene',
                'entrez':  'entrezgene',
                'probe':   'reporter.HG-U133_Plus_2'
            }.get(id_type, 'symbol')

            result = mg.querymany(
                gene_ids, scopes=scope,
                fields='symbol,ensembl.gene,entrezgene',
                species='human', returnall=True
            )

            id_map = {}
            for r in result.get('out', []):
                query = r.get('query', '')
                symbol = r.get('symbol', '')
                if symbol:
                    id_map[query] = symbol

            df[id_col + '_std'] = df[id_col].map(id_map).fillna(df[id_col])
            mapped_pct = sum(df[id_col + '_std'] != df[id_col]) / len(df) * 100
            print(f"  변환 성공률: {mapped_pct:.1f}%")
            return df

        except ImportError:
            print("  mygene 미설치: pip install mygene")
            df[id_col + '_std'] = df[id_col]
            return df

    def harmonize_proteins(self, df: pd.DataFrame,
                            id_col: str = "protein") -> pd.DataFrame:
        """
        단백질 ID → UniProt Swiss-Prot
        PD output의 Accession 컬럼 처리
        """
        def extract_uniprot(acc_string):
            # PD 출력: "P02741|CRP_HUMAN" 형태 처리
            if not isinstance(acc_string, str):
                return acc_string
            parts = acc_string.split("|")
            if len(parts) >= 2:
                return parts[0]  # UniProt ID 반환
            return acc_string

        df[id_col + '_uniprot'] = df[id_col].apply(extract_uniprot)
        return df

    def harmonize_metabolites(self, df: pd.DataFrame,
                               name_col: str = "metabolite") -> pd.DataFrame:
        """
        대사체 이름 → InChIKey / HMDB ID
        (MetaboAnalyst API 활용)
        """
        import requests

        names = df[name_col].dropna().unique().tolist()
        print(f"  대사체 ID 변환: {len(names)}개")

        # MetaboAnalyst name mapping API
        url = "https://rest.xialab.ca/api/mapcompounds"
        batch_size = 20
        id_map = {}

        for i in range(0, len(names), batch_size):
            batch = names[i:i+batch_size]
            try:
                r = requests.post(url, json={"compounds": batch}, timeout=30)
                if r.status_code == 200:
                    results = r.json()
                    for item in results.get("CompoundResults", []):
                        query = item.get("Query", "")
                        hmdb  = item.get("HMDB", "")
                        kegg  = item.get("KEGG", "")
                        if query:
                            id_map[query] = {"hmdb": hmdb, "kegg": kegg}
            except Exception:
                pass

        df['hmdb_id'] = df[name_col].map(
            lambda x: id_map.get(x, {}).get('hmdb', '')
        )
        df['kegg_id'] = df[name_col].map(
            lambda x: id_map.get(x, {}).get('kegg', '')
        )
        return df


# ═══════════════════════════════════════════════════════════════
# 2. 품질 필터링
# ═══════════════════════════════════════════════════════════════
class QualityFilter:

    def filter_mrna(self, count_mat: pd.DataFrame,
                     min_cpm: float = 1.0,
                     min_sample_pct: float = 0.2) -> pd.DataFrame:
        cpm = count_mat.divide(count_mat.sum()) * 1e6
        keep = (cpm > min_cpm).mean(axis=1) >= min_sample_pct
        filtered = count_mat[keep]
        print(f"  mRNA 필터: {count_mat.shape[0]} → {filtered.shape[0]}개")
        return filtered

    def filter_proteomics(self, prot_mat: pd.DataFrame,
                           max_missing_pct: float = 0.5) -> pd.DataFrame:
        missing = prot_mat.isna().mean(axis=1)
        keep = missing < max_missing_pct
        filtered = prot_mat[keep]
        print(f"  단백체 필터: {prot_mat.shape[0]} → {filtered.shape[0]}개")
        return filtered

    def filter_metabolomics(self, metab_mat: pd.DataFrame,
                             max_missing_pct: float = 0.5,
                             max_rsd: float = 0.3,
                             qc_cols: list = None) -> pd.DataFrame:
        missing = metab_mat.isna().mean(axis=1)
        keep = missing < max_missing_pct

        if qc_cols:
            qc = metab_mat[qc_cols]
            rsd = qc.std(axis=1) / qc.mean(axis=1)
            keep = keep & (rsd < max_rsd)

        filtered = metab_mat[keep]
        print(f"  대사체 필터: {metab_mat.shape[0]} → {filtered.shape[0]}개")
        return filtered

    def detect_outlier_samples(self, expr_mat: pd.DataFrame,
                                n_sd: float = 3.0) -> list:
        """IQR 기반 이상치 샘플 탐지"""
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        X = StandardScaler().fit_transform(expr_mat.T)
        pca = PCA(n_components=2)
        scores = pca.fit_transform(X)

        # Mahalanobis distance
        center = scores.mean(axis=0)
        dists  = np.sqrt(((scores - center) ** 2).sum(axis=1))
        threshold = dists.mean() + n_sd * dists.std()
        outliers = [expr_mat.columns[i] for i, d in enumerate(dists)
                    if d > threshold]

        if outliers:
            print(f"  ⚠️  이상치 샘플 {len(outliers)}개 탐지: {outliers}")
        return outliers


# ═══════════════════════════════════════════════════════════════
# 3. 정규화 + 배치 보정
# ═══════════════════════════════════════════════════════════════
class Normalizer:

    def normalize_mrna(self, count_mat: pd.DataFrame,
                        method: str = "vst") -> pd.DataFrame:
        if method == "log2cpm":
            cpm = count_mat.divide(count_mat.sum()) * 1e6
            return np.log2(cpm + 1)
        elif method == "vst":
            print("  VST: R DESeq2 사용 권장. log2CPM으로 대체.")
            cpm = count_mat.divide(count_mat.sum()) * 1e6
            return np.log2(cpm + 1)
        return count_mat

    def normalize_proteomics(self, prot_mat: pd.DataFrame) -> pd.DataFrame:
        log2 = np.log2(prot_mat.replace(0, np.nan))
        median_per_sample = log2.median()
        global_median     = log2.stack().median()
        return log2 - (median_per_sample - global_median)

    def normalize_metabolomics(self, metab_mat: pd.DataFrame,
                                method: str = "pqn") -> pd.DataFrame:
        log2 = np.log2(metab_mat.replace(0, np.nan) + 1)
        if method == "pqn":
            reference = log2.median(axis=1)
            quotients = log2.divide(reference, axis=0)
            correction = quotients.median()
            return log2.divide(correction)
        return log2

    def batch_correct(self, expr_mat: pd.DataFrame,
                       batch_labels: pd.Series) -> pd.DataFrame:
        """ComBat 배치 보정 (Python 구현)"""
        try:
            from inmoose.pycombat import pycombat_norm
            corrected = pycombat_norm(expr_mat, batch_labels)
            print(f"  ComBat 배치 보정 완료 ({len(batch_labels.unique())}개 배치)")
            return pd.DataFrame(corrected,
                                 index=expr_mat.index,
                                 columns=expr_mat.columns)
        except ImportError:
            print("  inmoose 미설치: pip install inmoose")
            return expr_mat


# ═══════════════════════════════════════════════════════════════
# 4. 메타분석 (여러 코호트 결합)
# ═══════════════════════════════════════════════════════════════
class MetaAnalysis:
    """
    여러 공개 코호트 + 자체 연구 결과를 메타분석으로 통합
    Random Effects Model (DerSimonian-Laird)
    """

    def run_meta_analysis(self,
                           study_effects: list,
                           study_se: list,
                           study_names: list = None) -> dict:
        """
        study_effects: 각 연구의 log2FC 목록
        study_se:      각 연구의 표준오차 목록
        """
        k = len(study_effects)
        effects = np.array(study_effects)
        ses     = np.array(study_se)
        weights = 1 / (ses ** 2)

        # Fixed effects
        fe_effect = np.sum(weights * effects) / np.sum(weights)
        fe_se     = np.sqrt(1 / np.sum(weights))
        q_stat    = np.sum(weights * (effects - fe_effect) ** 2)

        # Cochran's Q test
        q_pvalue = 1 - __import__('scipy').stats.chi2.cdf(q_stat, df=k-1)

        # I² heterogeneity
        i_squared = max(0, (q_stat - (k-1)) / q_stat * 100) if q_stat > 0 else 0

        # Random effects (DerSimonian-Laird)
        c = np.sum(weights) - np.sum(weights**2) / np.sum(weights)
        tau_sq = max(0, (q_stat - (k-1)) / c)
        re_weights = 1 / (ses**2 + tau_sq)
        re_effect  = np.sum(re_weights * effects) / np.sum(re_weights)
        re_se      = np.sqrt(1 / np.sum(re_weights))

        from scipy import stats
        z = re_effect / re_se
        re_pvalue = 2 * (1 - stats.norm.cdf(abs(z)))

        return {
            "pooled_log2fc": re_effect,
            "pooled_se":     re_se,
            "pooled_pvalue": re_pvalue,
            "ci_lower":      re_effect - 1.96 * re_se,
            "ci_upper":      re_effect + 1.96 * re_se,
            "i_squared":     i_squared,
            "q_stat":        q_stat,
            "q_pvalue":      q_pvalue,
            "n_studies":     k,
            "method":        "random_effects_DL"
        }

    def forest_plot(self, results_df: pd.DataFrame,
                    feature_col: str = "feature",
                    save_path: Path = None):
        """Forest plot 생성"""
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        df = results_df.sort_values("pooled_log2fc")
        n  = len(df)

        fig, ax = plt.subplots(figsize=(10, max(6, n * 0.35)))

        for i, (_, row) in enumerate(df.iterrows()):
            color = '#F44336' if row['pooled_log2fc'] > 0 else '#2196F3'
            ax.errorbar(
                row['pooled_log2fc'], i,
                xerr=[[row['pooled_log2fc'] - row['ci_lower']],
                       [row['ci_upper'] - row['pooled_log2fc']]],
                fmt='s', color=color, capsize=4,
                markersize=7 * (1 - row.get('pooled_pvalue', 0.5)),
                linewidth=1.5
            )

        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels(df[feature_col], fontsize=8)
        ax.set_xlabel("Pooled log2FC [95% CI]", fontsize=11)
        ax.set_title("Meta-Analysis Forest Plot\n(Periodontitis vs Control)",
                     fontsize=12, fontweight='bold')

        red = mpatches.Patch(color='#F44336', label='Upregulated')
        blue = mpatches.Patch(color='#2196F3', label='Downregulated')
        ax.legend(handles=[red, blue], fontsize=9)
        ax.grid(True, alpha=0.3, axis='x')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"  Forest plot 저장: {save_path}")
        plt.close()


# ═══════════════════════════════════════════════════════════════
# 5. DB 적재 (PostgreSQL)
# ═══════════════════════════════════════════════════════════════
class DBLoader:
    """통합 PostgreSQL DB에 큐레이션된 데이터 적재"""

    def __init__(self, db_url: str = None):
        self.db_url = db_url or os.environ.get(
            "OMICS_DB_URL",
            "postgresql://user:password@localhost:5432/omics_db"
        )
        self.engine = None

    def connect(self):
        try:
            from sqlalchemy import create_engine
            self.engine = create_engine(self.db_url)
            print(f"  DB 연결 성공: {self.db_url[:30]}...")
            return True
        except Exception as e:
            print(f"  DB 연결 실패: {e}")
            return False

    def load_source_metadata(self, meta_df: pd.DataFrame) -> int:
        """데이터 소스 메타데이터 적재"""
        if self.engine is None:
            return 0
        n = meta_df.to_sql('data_sources', self.engine,
                            if_exists='append', index=False,
                            method='multi', chunksize=100)
        return n or 0

    def load_differential_results(self,
                                    diff_df: pd.DataFrame,
                                    source_id: int,
                                    omics_type: str) -> int:
        """DEG/DEP/대사체 통계 결과 적재"""
        if self.engine is None:
            return 0
        diff_df = diff_df.copy()
        diff_df['source_id']  = source_id
        diff_df['omics_type'] = omics_type

        n = diff_df.to_sql('differential_analysis', self.engine,
                            if_exists='append', index=False,
                            method='multi', chunksize=500)
        return n or 0

    def update_biomarker_scores(self, candidates_df: pd.DataFrame):
        """바이오마커 후보 점수 업데이트"""
        if self.engine is None:
            return
        candidates_df.to_sql('biomarker_candidates', self.engine,
                               if_exists='append', index=False,
                               method='multi', chunksize=100)


# ═══════════════════════════════════════════════════════════════
# 6. Parquet 대용량 행렬 저장
# ═══════════════════════════════════════════════════════════════
class ParquetManager:
    """대용량 발현 행렬을 Parquet 파일로 효율 저장"""

    def __init__(self, base_dir: Path = CURATED):
        self.base = base_dir

    def save_matrix(self, df: pd.DataFrame,
                     name: str,
                     partition_col: str = None):
        """발현 행렬 → Parquet 저장"""
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq

            outpath = self.base / f"{name}.parquet"
            df.to_parquet(outpath, compression='snappy', index=True)
            size_mb = outpath.stat().st_size / 1e6
            print(f"  Parquet 저장: {outpath} ({size_mb:.1f} MB)")
            return outpath

        except ImportError:
            # fallback to feather
            outpath = self.base / f"{name}.feather"
            df.to_feather(outpath)
            print(f"  Feather 저장: {outpath}")
            return outpath

    def load_matrix(self, name: str,
                     features: list = None) -> pd.DataFrame:
        """Parquet 로드 (선택적 피처 로드로 메모리 절약)"""
        import pyarrow.parquet as pq

        path = self.base / f"{name}.parquet"
        if not path.exists():
            raise FileNotFoundError(f"{path}")

        if features:
            table = pq.read_table(path, columns=features)
        else:
            table = pq.read_table(path)

        return table.to_pandas()

    def merge_cohorts(self, cohort_names: list,
                       omics_type: str) -> pd.DataFrame:
        """
        여러 코호트 발현 행렬 병합 (배치 보정 전)
        """
        dfs = []
        batch_labels = []

        for cohort in cohort_names:
            try:
                df = self.load_matrix(f"{omics_type}_{cohort}")
                dfs.append(df)
                batch_labels.extend([cohort] * df.shape[1])
            except FileNotFoundError:
                print(f"  ⚠️  {cohort} 없음")

        if not dfs:
            return pd.DataFrame()

        # 공통 피처만 유지
        common_features = list(set.intersection(*[set(df.index) for df in dfs]))
        merged = pd.concat([df.loc[common_features] for df in dfs], axis=1)
        print(f"  {len(dfs)}개 코호트 병합: "
              f"{len(common_features)}개 피처 × {merged.shape[1]}개 샘플")

        return merged, pd.Series(batch_labels, index=merged.columns)


# ═══════════════════════════════════════════════════════════════
# MAIN: 큐레이션 파이프라인 실행
# ═══════════════════════════════════════════════════════════════
def run_curation_pipeline(
    source_name: str,
    omics_type: str,
    raw_data_path: Path,
    metadata_path: Path,
    output_name: str
):
    print(f"\n{'='*60}")
    print(f"큐레이션: {source_name} — {omics_type}")
    print(f"{'='*60}")

    harmonizer = IDHarmonizer()
    qfilter    = QualityFilter()
    normalizer = Normalizer()
    parquet    = ParquetManager()

    # 1. 원데이터 로드
    df = pd.read_csv(raw_data_path, index_col=0)
    meta = pd.read_csv(metadata_path, index_col=0)
    print(f"  로드: {df.shape}")

    # 2. ID 표준화
    if omics_type == "mRNA":
        # df index → HGNC symbol
        pass  # harmonizer.harmonize_genes(...)

    # 3. 품질 필터
    if omics_type == "mRNA":
        df = qfilter.filter_mrna(df)
    elif omics_type == "Proteomics":
        df = qfilter.filter_proteomics(df)
    elif omics_type == "Metabolomics":
        df = qfilter.filter_metabolomics(df)

    # 이상치 샘플 탐지
    outliers = qfilter.detect_outlier_samples(df)
    if outliers:
        df = df.drop(columns=outliers)

    # 4. 정규화
    if omics_type == "mRNA":
        df = normalizer.normalize_mrna(df, method="log2cpm")
    elif omics_type == "Proteomics":
        df = normalizer.normalize_proteomics(df)
    elif omics_type == "Metabolomics":
        df = normalizer.normalize_metabolomics(df)

    # 5. Parquet 저장
    outfile = parquet.save_matrix(df, output_name)

    print(f"\n✅ 큐레이션 완료: {outfile}")
    return df


if __name__ == "__main__":
    import argparse, os
    parser = argparse.ArgumentParser()
    parser.add_argument("--source",  default="GEO_GSE173078")
    parser.add_argument("--omics",   default="mRNA")
    parser.add_argument("--raw",     default="../data/public_db/geo/GSE173078/matrix.txt")
    parser.add_argument("--meta",    default="../data/public_db/metadata/GSE173078_meta.csv")
    parser.add_argument("--out",     default="mrna_GSE173078")
    args = parser.parse_args()

    run_curation_pipeline(
        args.source, args.omics,
        Path(args.raw), Path(args.meta),
        args.out
    )
