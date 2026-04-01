#!/bin/bash
# =============================================================
# 00_run_all.sh
# 치주염 멀티오믹스 분석 전체 파이프라인 실행
# 실행 순서: 01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10
# =============================================================

set -e   # 오류 발생 시 즉시 종료
LOGDIR="../logs"
mkdir -p "$LOGDIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "╔══════════════════════════════════════════════════════════╗"
echo "║  치주염 멀티오믹스 분석 파이프라인 시작                   ║"
echo "║  $(date)                                  ║"
echo "╚══════════════════════════════════════════════════════════╝"

# ── conda 환경 활성화 ──────────────────────────────────────
conda activate multiomics 2>/dev/null || true

# ── 실행 함수 ──────────────────────────────────────────────
run_step() {
  local step=$1
  local cmd=$2
  local desc=$3
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "  STEP ${step}: ${desc}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  eval "$cmd" 2>&1 | tee "${LOGDIR}/step${step}_${TIMESTAMP}.log"
  echo "  ✅ STEP ${step} 완료 ($(date +"%H:%M:%S"))"
}

# ── 파이프라인 실행 ────────────────────────────────────────
run_step "01" "Rscript 01_mrna_preprocessing.R" \
  "mRNA 전처리 (DESeq2 + WGCNA)"

run_step "02" "Rscript 02_proteomics_preprocessing.R" \
  "Proteomics 전처리 (limma + WGCNA)"

run_step "03" "python 03_metabolomics_preprocessing.py" \
  "Metabolomics 전처리 (PLS-DA + 통계)"

run_step "04" "python 04_mofa_integration.py" \
  "MOFA+ 비지도 통합"

run_step "05" "Rscript 05_diablo_integration.R" \
  "DIABLO 지도 통합"

run_step "06" "Rscript 06_immune_deconvolution.R" \
  "면역세포 디컨볼루션"

run_step "07" "python 07_stacking_biomarker_ml.py" \
  "Stacking Ensemble + SHAP 바이오마커 발굴"

run_step "08" "python 08_inflammatome_analysis.py" \
  "Inflammatome 분류 분석"

run_step "09" "python 09_ppi_network.py" \
  "PPI 네트워크 + Hub gene"

run_step "10" "python 10_pathintegrate_pathway.py" \
  "PathIntegrate + Joint Pathway 분석"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  ✅ 전체 파이프라인 완료!                                  ║"
echo "║  $(date)                                  ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "결과 위치:"
echo "  분석 결과: ../results/"
echo "  Figure:    ../figures/"
echo "  로그:      ${LOGDIR}/"
