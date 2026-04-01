# 논문 스켈레톤 — 치주염 멀티오믹스 바이오마커 발굴

**제목 (초안)**
> "Integrated Multi-Omics Analysis of Transcriptomics, Proteomics, and Metabolomics Reveals Novel Biomarker Panels for Periodontitis Diagnosis"

**투고 타깃 저널 (예시)**
- Journal of Clinical Periodontology (IF ~7)
- Journal of Proteome Research (IF ~4.5)
- Frontiers in Cellular and Infection Microbiology
- Journal of Periodontal Research

---

## Abstract (구조)

**Background:** Periodontitis is a chronic inflammatory disease...
**Objective:** We aimed to identify multi-omics biomarker panels...
**Methods:** Gingival tissue samples from 30 periodontitis patients and 30 healthy controls were subjected to bulk RNA-seq, quantitative DDA proteomics, and untargeted metabolomics. Single-omics analyses (DESeq2, limma, PLS-DA) were followed by multi-omics integration using MOFA+ and DIABLO. A stacking ensemble machine learning model (LASSO + Random Forest + XGBoost) with LOOCV validation was used for biomarker discovery. The inflammatome framework was applied to distinguish disease-specific from general inflammatory signatures.
**Results:** We identified X differentially expressed genes (DEGs), Y differentially expressed proteins (DEPs), and Z differential metabolites. MOFA+ integration revealed N disease-associated latent factors, while DIABLO identified correlated cross-omics signatures. The integrated biomarker panel (N_mRNA mRNA + N_prot protein + N_metab metabolite markers) achieved AUC = X.XX [CI], outperforming single-omics models (mRNA: X.XX, protein: X.XX, metabolite: X.XX). Inflammatome analysis revealed X% of DEGs were periodontitis-specific, with purine metabolism and NF-κB signaling as key shared pathways.
**Conclusions:** This study presents the first matched-sample tri-omics integration in periodontitis, providing a validated multi-omics biomarker panel with high diagnostic accuracy.

**Keywords:** Periodontitis, Multi-omics, Transcriptomics, Proteomics, Metabolomics, Biomarker, Machine learning

---

## 1. Introduction

### 1.1 Periodontitis: 임상적 중요성
- 전세계 유병률 (11%, 743 million명)
- 치아 손실의 주요 원인
- 전신 질환과의 연관성 (당뇨, 심혈관, 임신 합병증)
- 현재 진단의 한계: 임상/방사선학적 평가에만 의존

### 1.2 단일 오믹스 연구의 한계
- 기존 전사체 연구 (DESeq2, GEO 공개 데이터 활용)
- 단백체 / 대사체 단독 연구
- 비매칭 샘플 설계의 문제
- **gap**: 동일 샘플에서 3개 오믹스 통합 연구 부재

### 1.3 멀티오믹스 통합의 필요성
- 각 오믹스 레이어의 상호보완성
- MOFA+, DIABLO 등 최신 통합 방법론
- ML 기반 바이오마커 발굴 트렌드

### 1.4 연구 목적
> 본 연구는 치주염 환자 30명과 정상 대조군 30명의 **동일 치은 조직**에서 mRNA, 단백질, 대사체를 동시 측정하고, 최신 멀티오믹스 통합 방법론을 적용하여 진단 바이오마커 패널을 발굴하고 검증하는 것을 목적으로 한다.

---

## 2. Materials and Methods

### 2.1 Study Design and Sample Collection
```
- IRB 승인 번호 (XXX)
- 포함 기준: 2017 새 분류 기준 치주염 (Stage III/IV)
- 제외 기준: 항생제 투여, 흡연, 전신 질환
- 정상 대조군: 건강한 치은 조직
- 샘플: 치은 생검 (30 Control + 30 Periodontitis)
- 샘플 처리: 동일 샘플을 3등분 → 각 오믹스 실험
```

### 2.2 RNA Sequencing
```
- RNA 추출: TRIzol + RNeasy
- 라이브러리: poly-A enrichment, paired-end 150bp
- 시퀀서: Illumina NovaSeq 6000
- 매핑: STAR v2.7 → hg38
- 정량: Salmon (quasi-mapping)
- Count matrix: tximport
```

### 2.3 Quantitative DDA Proteomics
```
- 단백질 추출: RIPA buffer
- 소화: Trypsin (FASP)
- LC-MS/MS: Thermo Q Exactive
- 데이터베이스: UniProt Human (SwissProt)
- 소프트웨어: Proteome Discoverer 2.4
- 정량: Normalized Abundance (TMT or LFQ)
- FDR: 1% (protein level)
```

### 2.4 Untargeted Metabolomics
```
- 추출: methanol/water (80:20)
- 분석: UHPLC-MS/MS (Positive + Negative mode)
- 기기: Thermo Q Exactive HF
- Peak picking: XCMS or MZmine3
- 어노테이션: HMDB, KEGG, mzCloud (MS2)
- 정규화: PQN (Probabilistic Quotient)
```

### 2.5 mRNA Differential Expression Analysis
```
- 도구: DESeq2 (v1.42)
- 정규화: VST (Variance Stabilizing Transformation)
- 임계값: padj < 0.05, |log2FC| > 1.0
- 배치 보정: ComBat-seq (필요시)
- 모듈 분석: WGCNA (signed hybrid, power=X, minModuleSize=20)
- 경로 분석: GSEA (gseKEGG, gseGO), ORA (enrichKEGG)
- 면역세포 추정: CIBERSORT (absolute mode) + xCell
```

### 2.6 Proteomics Differential Expression Analysis
```
- 결측치 처리: MinProb (MNAR) + KNN (MAR), cutoff 50%
- 정규화: Median centering (재확인)
- DEP 분석: limma + BH correction
- 임계값: adj.P < 0.05, |log2FC| > 0.58 (1.5×)
- 모듈 분석: WGCNA (bicor, signed hybrid)
- mRNA-Protein concordance: Spearman 상관, 4-quadrant 분류
```

### 2.7 Metabolomics Analysis
```
- 전처리: Log2 변환, Pareto scaling
- 결측치: KNN imputation
- 다변량: PLS-DA + OPLS-DA (ropls 패키지)
- 유의성: VIP > 1.0 AND padj < 0.05 AND |FC| > 1.5
- 검증: Permutation test (200회, Q2 > 0)
- 경로: MetaboAnalyst 6.0
```

### 2.8 Multi-Omics Integration
```
2.8.1 MOFA+ (비지도)
- 패키지: mofapy2 (Python)
- Factors: 15 초기 → auto-pruning
- Likelihood: Gaussian (3 views)
- convergence: "medium", seed=42
- Factor-phenotype: Spearman 상관

2.8.2 DIABLO (지도)
- 패키지: mixOmics (R)
- Design matrix: 0.1 (off-diagonal)
- keepX: tune.block.splsda() LOO-CV, BER measure
- ncomp: 2
- Circos plot: |r| > 0.7
```

### 2.9 PPI Network Analysis
```
- 데이터베이스: STRING v12.0
- 신뢰도 임계값: medium confidence (≥ 0.4)
- 네트워크 분석: NetworkX
- Hub gene 기준: 통합 중심성 점수 (degree 40% + betweenness 40% + closeness 20%)
```

### 2.10 Inflammatome Analysis
```
- 유전자 세트: Inflammatome (Díaz-Pinés Cort et al., 2025/2026)
  - Inflammation signature: top 100 genes
  - Inflammatome: top 2,000 genes
- 분류: Concordant (DEG ∩ inflammatome) vs Periodontitis-specific (DEG - inflammatome)
- Inflammation score: ssGSEA (per-sample)
```

### 2.11 Machine Learning Biomarker Discovery
```
- Feature pool: 유의 DEG + DEP + 대사체 합집합
- Feature selection: LASSO (1-SE rule, 10-fold CV)
- Base models: LASSO-LR + Random Forest (n=500) + XGBoost
- Meta-learner: Logistic Regression (Stacking)
- 검증: LOOCV (n=60 최적)
- 해석: SHAP (TreeExplainer, beeswarm + per-omics contribution)
- 신뢰구간: Bootstrap AUC CI (n=1,000)
- 최소 패널: SHAP top 10 피처
```

### 2.12 External Validation
```
- GEO 데이터셋: GSE173078, GSE152042 (치주염 RNA-seq)
- MetaboLights: MTBLS8357 (치주염 gingival metabolomics)
- mRNA 바이오마커 적용 → AUC 계산
```

### 2.13 Joint Pathway Analysis
```
- 도구: MetaboAnalyst 6.0 Joint Pathway Analysis
- 입력: DEG 목록 + 유의 대사체 목록
- 데이터베이스: KEGG (Homo sapiens)
- 통계: Hypergeometric test, degree centrality topology
- 유의 기준: p.combine < 0.05, impact > 0.1
```

### 2.14 Statistical Analysis
```
- 모든 통계: R v4.4 또는 Python 3.11
- 다중 검정 보정: BH (Benjamini-Hochberg)
- 그룹 비교: Mann-Whitney U test (비모수)
- 유의 수준: α = 0.05
- 재현성: seed=42, 모든 코드 GitHub 공개
```

---

## 3. Results

### 3.1 Study Cohort and Quality Control
```
- 샘플 특성 테이블 (Table 1)
- 각 오믹스 QC 통과율
- PCA — 그룹 분리 확인 (Fig 1A-C)
- Outlier 없음 확인
```

### 3.2 Single-Omics Differential Analysis

#### 3.2.1 Transcriptomics
```
- DEG 수: X개 상향, Y개 하향 (Fig 2A: Volcano)
- Top DEG: [유전자명] 포함 언급
- GSEA 결과: NF-κB, IL signaling, ECM remodeling (Fig 2B: Dotplot)
- WGCNA: Z개 모듈, 치주염 연관 모듈 (Fig 2C: Heatmap)
```

#### 3.2.2 Proteomics
```
- DEP 수: X개 상향, Y개 하향 (Fig 3A: Volcano)
- MMP-8, S100A8/A9 등 언급
- limma pathway 결과 (Fig 3B)
- mRNA-Protein concordance: N개 일치, M개 불일치 (Fig 3C: 4분면)
```

#### 3.2.3 Metabolomics
```
- PLS-DA 그룹 분리: R²=X.XX, Q²=X.XX, permutation p<0.05 (Fig 4A)
- 유의 대사체: N개 (VIP>1, padj<0.05, |FC|>1.5) (Fig 4B: Volcano)
- 주요 대사체: deoxyinosine, arachidonic acid 등
- Pathway: 퓨린 대사, 아미노산 대사 (Fig 4C)
```

### 3.3 Immune Microenvironment
```
- CIBERSORT + xCell 일치 결과
- Neutrophil↑↑, M1 Macro↑, Treg↓, NK↓ (Fig 5: Violin + stacked bar)
- Hub gene — 면역세포 상관 (Fig 5B: Heatmap)
```

### 3.4 Multi-Omics Integration

#### 3.4.1 MOFA+ 비지도 통합
```
- N개 유의 Factor (|r| > 0.5 with periodontitis)
- Factor별 오믹스 설명력 (Fig 6A: R² decomposition)
- Factor score 분포 (Fig 6B: scatter)
- 각 Factor의 top loading 피처 (Fig 6C)
```

#### 3.4.2 DIABLO 지도 통합
```
- 오믹스 간 상관 구조 (Fig 7A: Circos plot)
- Sample plot (Fig 7B: PLS-DA space)
- Arrow plot — 오믹스 일치도 (Fig 7C)
- 선택된 피처: N_mRNA + N_prot + N_metab (Table 2)
```

### 3.5 PPI Network and Hub Genes
```
- PPI 네트워크 통계: N nodes, M edges
- Top hub gene: [유전자명] (degree, betweenness)
- 다중 증거 hub (DEG∩DEP∩PPI): N개 (Fig 8)
```

### 3.6 Inflammatome Analysis
```
- DEG 중 X%가 inflammatome 포함 (범염증)
- 치주염 특이 DEG: N개 상향, M개 하향 (Fig 9A: Venn)
- Inflammation score: 환자군에서 유의 상승 (Fig 9B)
- 특이 마커: [유전자명 3-5개] 언급
```

### 3.7 Biomarker Discovery and Validation

#### 3.7.1 Feature Selection
```
- LASSO 1-SE rule: N개 피처 선택 (전체 X개 중)
- 오믹스별 구성: N_mRNA + N_prot + N_metab
```

#### 3.7.2 Model Performance
```
- 단독 AUC:
  mRNA only:       AUC = X.XX [CI]
  Protein only:    AUC = X.XX [CI]
  Metabolite only: AUC = X.XX [CI]
- 통합 패널:       AUC = X.XX [CI] ← 가장 높음
- 비교 Figure (Fig 10A: ROC, 10B: AUC bar)
```

#### 3.7.3 SHAP 해석
```
- Top 20 피처 SHAP beeswarm (Fig 11A)
- 오믹스별 기여도: mRNA X%, Protein Y%, Metabolite Z% (Fig 11B: pie)
```

#### 3.7.4 최소 임상 패널
```
- SHAP top 10 → 최소 패널
- AUC = X.XX, Sensitivity = X.XX, Specificity = X.XX (Fig 12)
- 패널 구성: Table 3
```

#### 3.7.5 외부 검증
```
- GSE173078: AUC = X.XX (mRNA 마커)
- GSE152042: AUC = X.XX
- MTBLS8357: 대사체 마커 검증
```

### 3.8 Joint Pathway Analysis
```
- mRNA + 대사체 동시 유의 경로: N개
- 최상위 경로: 퓨린 대사 (p.combine=X.XX, impact=X.XX)
- Bubble chart: impact vs significance (Fig 13)
```

---

## 4. Discussion

### 4.1 멀티오믹스 통합의 가치
```
- 단독 AUC vs 통합 AUC 비교 → 통합의 당위성
- mRNA-Protein concordance: Concordant 마커의 강점
- Discordant 마커: post-translational regulation 시사
```

### 4.2 핵심 바이오마커 해석
```
- mRNA 상위 마커: 기존 문헌과 비교
- DEP 상위 마커: MMP-8, S100A8/A9 등 기존 문헌 일치
- 대사체 마커: 퓨린 대사 → Luo 2023, Chu 2024 일치
```

### 4.3 면역 microenvironment
```
- Neutrophil 상승 해석 → 조직 파괴 기전
- Treg 감소 → 면역 억제 기능 소실
- Hub gene과 면역세포 상관 → 분자 기전 연결
```

### 4.4 Inflammatome 적용의 의의
```
- 치주염 특이 마커 vs 범염증 마커 구분
- Disease-specific 마커의 진단적 가치
```

### 4.5 연구의 강점
```
✅ 동일 샘플 3개 오믹스 매칭 설계
✅ 적절한 샘플 크기 (n=60)
✅ LOOCV + Bootstrap CI 엄격한 검증
✅ 외부 코호트 검증
✅ 공개 코드 + 데이터 (FAIR)
```

### 4.6 연구의 한계
```
- 단일 기관 코호트 → 외부 검증 필요
- 치주염 중증도 스펙트럼 미반영
- 기능적 검증 없음 (in vitro/vivo 추가 필요)
- 표본 크기 (n=60) 제한
```

### 4.7 결론
> 본 연구는 치주염에서 최초의 매칭 3-오믹스 통합 분석을 수행하여 X개의 통합 바이오마커 패널을 발굴하고 LOOCV + 외부 코호트로 검증하였다. 특히 inflammatome 프레임워크를 통해 치주염 특이적 분자 서명을 분리하였으며, 이는 정밀 치주 진단의 토대를 제공한다.

---

## 5. References

### 필수 인용 논문 (최소 30편)

**치주염 오믹스:**
1. Luo et al. (2023). Multi-omics study... Arch Oral Biol. DOI:10.1016/j.archoralbio.2023.105680
2. Chu et al. (2024). Untargeted Metabolomics... J Proteome Res. DOI:10.1021/acs.jproteome.3c00105
3. Front. Med. (2025). Multi-omics immune biomarkers. DOI:10.3389/fmed.2025.1640961
4. Sci. Rep. (2025). Ensemble learning biomarkers. DOI:10.1038/s41598-025-18017-7

**방법론:**
5. DESeq2: Love et al. (2014). Genome Biol.
6. limma: Ritchie et al. (2015). Nucleic Acids Res.
7. WGCNA: Langfelder & Horvath (2008). BMC Bioinform.
8. MOFA+: Argelaguet et al. (2020). Genome Biol.
9. DIABLO/mixOmics: Singh et al. (2019). PLOS Comput Biol.
10. CIBERSORT: Newman et al. (2019). Nat Methods.
11. PathIntegrate: Wieder et al. (2024). PLoS Comput Biol.
12. Inflammatome: Díaz-Pinés Cort et al. (2025). Cell Reports.
13. MetaboAnalyst 6.0: Pang et al. (2024). Nucleic Acids Res.
14. SHAP guide: PMC11513550 (2024). Clin Transl Sci.

**유사 질환 참고:**
15. He et al. (2023). SLE Proteomics+Metabolomics. Clin Immunol.
16. IBD multi-omics. PMC11792892 (2025).
17. RA multi-omics. Arthritis Res Ther (2024).

---

## Tables (예정)

| Table | 내용 |
|-------|------|
| Table 1 | Clinical characteristics (Control vs Periodontitis) |
| Table 2 | DIABLO selected features (N_mRNA + N_prot + N_metab) |
| Table 3 | Final minimal biomarker panel (name, omics, SHAP, AUC) |
| Table S1 | All DEGs (DESeq2) |
| Table S2 | All DEPs (limma) |
| Table S3 | All differential metabolites |
| Table S4 | GSEA/ORA pathway results |
| Table S5 | MOFA+ factor-phenotype correlation |
| Table S6 | PPI hub gene centrality scores |
