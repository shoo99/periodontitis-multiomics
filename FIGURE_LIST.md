# Figure 목록 최종 확정 — 치주염 멀티오믹스 논문

## 📊 Main Figures (12개)

---

### Figure 1 — Study Overview & QC (3-panel)
| Panel | 내용 | 생성 스크립트 | 데이터 소스 |
|-------|------|-------------|------------|
| A | 연구 디자인 workflow 다이어그램 | 수동 제작 (BioRender) | — |
| B | 3개 오믹스 PCA (3×1 배열) | 01, 02, 03 | VST, prot_norm, metab_pareto |
| C | 샘플별 오믹스 커버리지 bar | 03 | 각 omic 검출률 |

**파일명:** `Fig1_Overview_QC.pdf`
**해상도:** 300 DPI, 180mm width

---

### Figure 2 — mRNA DEG Analysis (3-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | Volcano plot (DEG, n=X) | 01 |
| B | Heatmap — Top 50 DEG (Z-score) | 01 |
| C | GSEA KEGG dotplot (top 15 pathways) | 01 |

**파일명:** `Fig2_mRNA_DEG.pdf`
**주요 강조:** 상향/하향 유전자 색상 구분, 주요 유전자 레이블

---

### Figure 3 — Proteomics DEP Analysis (3-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | Volcano plot (DEP) | 02 |
| B | Heatmap — Top 40 DEP | 02 |
| C | mRNA-Protein Concordance 4분면 산점도 | 02 |

**파일명:** `Fig3_Proteomics_DEP.pdf`
**Panel C 특이점:** Concordant(보라)/Discordant(주황)/mRNA-only(파랑)/Prot-only(초록) 4색 구분

---

### Figure 4 — Metabolomics Analysis (3-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | PLS-DA score plot (R², Q², permutation) | 03 |
| B | Volcano plot (대사체, VIP>1 강조) | 03 |
| C | Top 12 대사체 boxplot (2×6 배열) | 03 |

**파일명:** `Fig4_Metabolomics.pdf`
**Panel A 필수 정보:** R²=X.XX, Q²=X.XX, permutation p=X.XX

---

### Figure 5 — Immune Microenvironment (2-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | Stacked bar plot (CIBERSORT 22 cells × 60 samples) | 06 |
| B | Violin plot — 유의 세포 유형 비교 (8종) | 06 |

**파일명:** `Fig5_Immune_Deconv.pdf`
**주요 예상 결과:** Neutrophil↑, M1 Macro↑, Treg↓, NK↓

---

### Figure 6 — MOFA+ Integration (3-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | Variance explained per factor per view (heatmap) | 04 |
| B | Factor score scatter (Factor1 vs Factor2, 색=그룹) | 04 |
| C | Factor-Periodontitis 상관 막대 | 04 |

**파일명:** `Fig6_MOFA_Integration.pdf`
**Panel A:** 오믹스별 R² decomposition — 각 view 기여도 시각화

---

### Figure 7 — DIABLO Integration (3-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | Circos plot (Cross-omics correlation, |r|>0.7) | 05 |
| B | Sample plot (Comp1 vs Comp2, 3 views) | 05 |
| C | Loading plot — 각 오믹스 top 피처 | 05 |

**파일명:** `Fig7_DIABLO_Integration.pdf`
**Panel A:** mRNA-Protein/mRNA-Metabolite/Protein-Metabolite 연결선 시각화

---

### Figure 8 — PPI Network (1-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| — | PPI 네트워크 (노드 크기=degree, 색=증거 레벨) | 09 |

**파일명:** `Fig8_PPI_Network.pdf`
**색 코딩:**
- 보라: Hub + DEG + DEP (최고 증거)
- 빨강: Hub + DEG
- 주황: Hub + DEP
- 파랑: Hub only

---

### Figure 9 — Inflammatome Analysis (2-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | Venn diagram (DEG vs inflammatome top100) | 08 |
| B | Inflammation score violin (Control vs Periodontitis) | 08 |

**파일명:** `Fig9_Inflammatome.pdf`
**Panel A:** 치주염 특이 vs 범염증 DEG 비율 — 논문 novelty 강조

---

### Figure 10 — Biomarker Model Comparison (2-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | ROC curve (통합 패널, LOOCV) | 07 |
| B | AUC bar + 95% CI (단독 vs 통합 비교) | 07 |

**파일명:** `Fig10_ROC_AUC_Comparison.pdf`
**핵심 메시지:** 통합 AUC > 단독 AUC → 멀티오믹스의 당위성

---

### Figure 11 — SHAP Interpretation (2-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | SHAP beeswarm plot (top 20 features) | 07 |
| B | 오믹스별 SHAP 기여도 파이차트 | 07 |

**파일명:** `Fig11_SHAP_Interpretation.pdf`
**Panel A:** 각 피처의 방향성(빨강=상승기여, 파랑=억제기여) 포함

---

### Figure 12 — Minimal Biomarker Panel (2-panel)
| Panel | 내용 | 생성 스크립트 |
|-------|------|-------------|
| A | 최소 패널 ROC curve (AUC + 95% CI) | 07 |
| B | Confusion matrix (Sensitivity, Specificity) | 07 |

**파일명:** `Fig12_Minimal_Panel_ROC.pdf`
**임상 적용성:** ≤10개 마커로 높은 AUC 달성 강조

---

## 📎 Supplementary Figures (6개)

| Figure | 내용 | 스크립트 |
|--------|------|---------|
| Fig S1 | WGCNA soft-power 선택 그래프 (mRNA + Prot) | 01, 02 |
| Fig S2 | WGCNA Module-Trait heatmap (mRNA + Prot) | 01, 02 |
| Fig S3 | Metabolomics top 20 VIP bar + OPLS-DA permutation | 03 |
| Fig S4 | Inflammatome — Specific vs Shared heatmap (top 30 DEG) | 08 |
| Fig S5 | Joint Pathway bubble chart (mRNA + Metabolomics) | 10 |
| Fig S6 | External validation ROC (GSE173078, GSE152042) | 별도 |

---

## 📐 Figure 제작 규격

| 항목 | 규격 |
|------|------|
| 해상도 | 300 DPI (최소) |
| 파일 형식 | PDF (벡터) + PNG (300 DPI) |
| 너비 | Full page: 180mm / Half page: 85mm |
| 폰트 | Arial 또는 Helvetica (8-10pt 권장) |
| 색상 모드 | RGB (온라인) / CMYK (인쇄) |
| 패널 레이블 | A, B, C ... (bold, 12pt) |
| 통계 표시 | * p<0.05, ** p<0.01, *** p<0.001, ns |

---

## 🎨 색상 팔레트 (전 Figure 공통)

```python
# 그룹 색상
CTRL_COLOR   = "#4CAF50"   # 초록 (Control)
PERIO_COLOR  = "#F44336"   # 빨강 (Periodontitis)

# 오믹스 색상
MRNA_COLOR   = "#2196F3"   # 파랑 (mRNA)
PROT_COLOR   = "#FF5722"   # 주황 (Protein)
METAB_COLOR  = "#4CAF50"   # 초록 (Metabolite)
INTEGR_COLOR = "#9C27B0"   # 보라 (Integrated)

# DEG 방향
UP_COLOR     = "#F44336"   # 상향
DOWN_COLOR   = "#2196F3"   # 하향
NS_COLOR     = "#9E9E9E"   # NS
```

---

## ✅ Figure 완성도 체크리스트

```
□ Fig 1  — Study Overview           (수동 + 자동)
□ Fig 2  — mRNA DEG                 (01_mrna_preprocessing.R)
□ Fig 3  — Proteomics DEP           (02_proteomics_preprocessing.R)
□ Fig 4  — Metabolomics             (03_metabolomics_preprocessing.py)
□ Fig 5  — Immune Deconv            (06_immune_deconvolution.R)
□ Fig 6  — MOFA+                    (04_mofa_integration.py)
□ Fig 7  — DIABLO Circos            (05_diablo_integration.R)
□ Fig 8  — PPI Network              (09_ppi_network.py)
□ Fig 9  — Inflammatome             (08_inflammatome_analysis.py)
□ Fig 10 — ROC Comparison           (07_stacking_biomarker_ml.py)
□ Fig 11 — SHAP                     (07_stacking_biomarker_ml.py)
□ Fig 12 — Minimal Panel            (07_stacking_biomarker_ml.py)
□ Fig S1 — WGCNA Power              (01, 02)
□ Fig S2 — WGCNA Module-Trait       (01, 02)
□ Fig S3 — VIP + Permutation        (03)
□ Fig S4 — Inflammatome Heatmap     (08)
□ Fig S5 — Joint Pathway            (10)
□ Fig S6 — External Validation      (별도 스크립트)
```

**총 Figure: Main 12개 + Supplementary 6개 = 18개**
