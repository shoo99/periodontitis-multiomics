# 🦷 Periodontitis Multi-Omics Analysis Pipeline

치주염(Periodontitis) 환자 vs 정상 대조군의 **mRNA + Proteomics + Metabolomics** 통합 분석을 위한 자동화 파이프라인.

## 연구 개요

| 항목 | 내용 |
|------|------|
| 샘플 | 정상 30 vs 치주염 30 (총 60명) |
| 오믹스 ① | mRNA (STAR+Salmon count matrix) |
| 오믹스 ② | Quant DDA Proteomics (Proteome Discoverer normalized abundance) |
| 오믹스 ③ | Untargeted Metabolomics (intensity table) |
| 목표 | 바이오마커 발굴 + 통합 분석 보고서 자동 생성 |

---

## 📂 파이프라인 구조

```
pipeline/
├── 00_run_all.sh                  # 전체 실행 마스터 스크립트
│
├── 01_mrna_preprocessing.R        # mRNA: DESeq2 DEG + WGCNA + GSEA + 면역세포 디컨볼루션
├── 02_proteomics_preprocessing.R  # Proteomics: limma DEP + WGCNA(bicor) + mRNA-Protein concordance
├── 03_metabolomics_preprocessing.py  # Metabolomics: PLS-DA + VIP + 경로 분석
│
├── 04_mofa_integration.py         # MOFA+ 비지도 멀티오믹스 통합
├── 05_diablo_integration.R        # DIABLO 지도 멀티오믹스 통합 + Circos plot
├── 06_immune_deconvolution.R      # CIBERSORT + xCell 면역세포 추정
│
├── 07_stacking_biomarker_ml.py    # Stacking Ensemble + SHAP 바이오마커 발굴
├── 08_inflammatome_analysis.py    # Inflammatome 분류 (범염증 vs 치주염 특이)
├── 09_ppi_network.py              # STRING API PPI 네트워크 + Hub gene
└── 10_pathintegrate_pathway.py    # PathIntegrate + Joint Pathway 분석

references/
├── inflammatome_top100.txt        # 범염증 100 유전자 (Cell Reports 2025)
├── inflammatome_top2000.txt       # 범염증 2000 유전자
└── inflammatome_full_ranked.tsv   # 전체 ranke 목록 (ENSG ID + gene name)
```

---

## 🔬 분석 흐름

```
[mRNA count]  [PD Abundance]  [Metabolite intensity]
      ↓              ↓                 ↓
   DESeq2          limma           PLS-DA
   WGCNA           WGCNA(bicor)    VIP+t-test
   GSEA            mRNA↔Prot       Pathway(gseapy)
      ↓              ↓                 ↓
   ──────────────────────────────────────
         MOFA+ (비지도 통합)
         DIABLO (지도 통합)
         Circos (Cross-omics 상관)
         면역세포 디컨볼루션
   ──────────────────────────────────────
                    ↓
       LASSO(1-SE) + RF + XGBoost
       Stacking Ensemble (LOOCV)
       SHAP 해석 + 오믹스 기여도
                    ↓
       Inflammatome 분류
       PPI Hub gene
       Joint Pathway (MetaboAnalyst)
                    ↓
       최소 임상 패널 (≤10개) + AUC
       자동 보고서 생성
```

---

## 🛠️ 환경 설정

### Python 환경
```bash
conda create -n multiomics python=3.11
conda activate multiomics
pip install -r requirements_python.txt
```

### R 패키지
```r
# Bioconductor
BiocManager::install(c(
  "DESeq2", "WGCNA", "limma", "clusterProfiler",
  "org.Hs.eg.db", "GSVA", "mixOmics",
  "immunedeconv", "impute", "pcaMethods"
))
# CRAN
install.packages(c("ggplot2","pheatmap","RColorBrewer","patchwork",
                   "dplyr","tibble","ggrepel","biomaRt"))
```

---

## 🚀 실행

```bash
# 1. 데이터 준비 (아래 경로에 파일 위치)
#    data/processed/mrna_counts.csv
#    data/processed/proteomics_abundance.csv
#    data/processed/metabolomics_intensity.csv
#    data/processed/sample_metadata.csv  (columns: group)

# 2. 전체 파이프라인 실행
cd pipeline
bash 00_run_all.sh

# 3. 개별 스텝 실행
Rscript 01_mrna_preprocessing.R
python  03_metabolomics_preprocessing.py
```

---

## 📊 출력 결과

| 디렉토리 | 내용 |
|----------|------|
| `results/mrna/` | DEG, VST matrix, WGCNA hub, GSEA |
| `results/proteomics/` | DEP, log2 matrix, WGCNA hub, concordance |
| `results/metabolomics/` | 유의 대사체, PLS-DA, 경로 |
| `results/mofa/` | MOFA+ factor scores, 표현형 상관 |
| `results/diablo/` | DIABLO 선택 피처, BER |
| `results/immune/` | 면역세포 비율, Hub-면역 상관 |
| `results/ml/` | 바이오마커 패널, SHAP, AUC 비교 |
| `results/ppi/` | PPI 상호작용, Hub centrality |
| `figures/` | 모든 Figure (PDF/PNG) |

---

## 📚 주요 참고 문헌

| # | 논문 | 저널 | 연도 | 역할 |
|---|------|------|------|------|
| 1 | Luo et al. | Arch Oral Biol | 2023 | 치주염 transcriptomics+metabolomics |
| 2 | Chu et al. | J Proteome Res | 2024 | 치주염 gingival metabolomics |
| 3 | Díaz-Pinés Cort et al. | Cell Reports | 2025 | Inflammatome gene sets |
| 4 | Wieder et al. | PLoS Comput Biol | 2024 | PathIntegrate |
| 5 | Front. Med. | 2025 | 치주염 면역 멀티오믹스 |
| 6 | He et al. | Clin. Immunol. | 2023 | SLE proteomics+metabolomics RF |
| 7 | IBD PMC11792892 | 2025 | MOFA+ML 환자 서브그룹 |

---

## 👤 Contact

- **Disease**: Periodontitis (치주염)
- **Omics**: mRNA + DDA Proteomics + Untargeted Metabolomics
- **Samples**: n=60 (30 Control, 30 Periodontitis)
- **Pipeline version**: v1.0 (2026-04)
