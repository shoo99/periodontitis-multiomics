#!/usr/bin/env Rscript
# =============================================================
# 05_diablo_integration.R
# DIABLO (mixOmics) — 지도 멀티오믹스 통합 + 바이오마커 선별
# =============================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
})

set.seed(42)

OUTDIR <- "../results/diablo"
FIGDIR <- "../figures"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

cat("=" , rep("=",58), "\n", sep="")
cat("05. DIABLO Multi-Omics Integration\n")
cat("=" , rep("=",58), "\n\n", sep="")


# ══════════════════════════════════════════════════════════════
# STEP 1: 데이터 로드
# ══════════════════════════════════════════════════════════════
cat("[1/5] 데이터 로드...\n")

mrna_vst  <- as.matrix(read.csv("../results/mrna/mrna_vst_matrix.csv", row.names=1))
prot_log2 <- as.matrix(read.csv("../results/proteomics/proteomics_log2_normalized.csv", row.names=1))
metab_log2<- as.matrix(read.csv("../results/metabolomics/metabolomics_log2.csv", row.names=1))
meta      <- read.csv("../data/processed/sample_metadata.csv", row.names=1)

# 공통 샘플
common_s <- Reduce(intersect,
                   list(colnames(mrna_vst), colnames(prot_log2),
                        colnames(metab_log2), rownames(meta)))
meta <- meta[common_s, , drop=FALSE]
Y    <- factor(meta$group, levels=c("Control","Periodontitis"))

# 상위 변동 피처 선택 (DIABLO 속도 최적화)
top_var <- function(mat, n) {
  mads <- apply(mat, 1, mad)
  mat[names(sort(mads, decreasing=TRUE))[1:min(n, nrow(mat))], ]
}

X <- list(
  mRNA       = t(top_var(mrna_vst[, common_s], 2000)),
  Protein    = t(top_var(prot_log2[, common_s], 1000)),
  Metabolite = t(metab_log2[, common_s])
)
cat(sprintf("  mRNA: %d, Protein: %d, Metabolite: %d 피처\n",
            ncol(X$mRNA), ncol(X$Protein), ncol(X$Metabolite)))


# ══════════════════════════════════════════════════════════════
# STEP 2: Design Matrix 설정 + 초기 모델
# ══════════════════════════════════════════════════════════════
cat("\n[2/5] Design matrix 설정...\n")

# Design = 0.1: 분류력 최대화 (바이오마커 발굴 목적)
design <- matrix(0.1, nrow=3, ncol=3,
                 dimnames=list(c("mRNA","Protein","Metabolite"),
                               c("mRNA","Protein","Metabolite")))
diag(design) <- 0
cat("  Design matrix (off-diagonal = 0.1):\n")
print(design)

# 초기 모델 (keepX 튜닝 전)
diablo_init <- block.splsda(X, Y, ncomp=2, design=design)
cat("  초기 모델 완료\n")


# ══════════════════════════════════════════════════════════════
# STEP 3: keepX 튜닝 (LOO-CV)
# ══════════════════════════════════════════════════════════════
cat("\n[3/5] keepX 튜닝 (LOO-CV)...\n")
cat("  (n=60 소규모 → Leave-One-Out CV 사용)\n")

test_keepX <- list(
  mRNA       = c(5, 10, 15, 20, 25, 30),
  Protein    = c(5, 10, 15, 20, 25),
  Metabolite = c(5, 10, 15, 20)
)

tryCatch({
  tune_res <- tune.block.splsda(
    X            = X,
    Y            = Y,
    ncomp        = 2,
    test.keepX   = test_keepX,
    design       = design,
    validation   = "loo",         # LOO-CV
    measure      = "BER",         # Balanced Error Rate
    nrepeat      = 1,
    BPPARAM      = BiocParallel::SerialParam(),
    verbose      = FALSE
  )

  optimal_keepX <- tune_res$choice.keepX
  cat("  최적 keepX:\n")
  print(optimal_keepX)

  write.csv(as.data.frame(tune_res$error.rate),
            file.path(OUTDIR, "diablo_tune_BER.csv"))

}, error = function(e) {
  cat(sprintf("  튜닝 오류: %s\n", e$message))
  cat("  기본값 사용\n")
  optimal_keepX <<- list(mRNA=c(15,10), Protein=c(10,8), Metabolite=c(10,8))
})


# ══════════════════════════════════════════════════════════════
# STEP 4: 최종 DIABLO 모델
# ══════════════════════════════════════════════════════════════
cat("\n[4/5] 최종 DIABLO 모델 학습...\n")

diablo_final <- block.splsda(
  X       = X,
  Y       = Y,
  ncomp   = 2,
  keepX   = optimal_keepX,
  design  = design
)

# 선택된 피처 추출
selected_features <- lapply(names(X), function(omics) {
  feats <- lapply(1:2, function(comp) {
    selectVar(diablo_final, block=omics, comp=comp)$name
  })
  unique(unlist(feats))
})
names(selected_features) <- names(X)

cat("  선택된 피처:\n")
for (nm in names(selected_features)) {
  cat(sprintf("    %s: %d개\n", nm, length(selected_features[[nm]])))
}

# 저장
for (nm in names(selected_features)) {
  write.csv(
    data.frame(feature=selected_features[[nm]], omics=nm),
    file.path(OUTDIR, paste0("diablo_selected_", nm, ".csv")),
    row.names=FALSE
  )
}

# 전체 selected features 합치기
all_selected <- do.call(rbind, lapply(names(selected_features), function(nm) {
  data.frame(feature=selected_features[[nm]], omics=nm)
}))
write.csv(all_selected, file.path(OUTDIR, "diablo_all_selected_features.csv"),
          row.names=FALSE)


# ══════════════════════════════════════════════════════════════
# STEP 5: 시각화
# ══════════════════════════════════════════════════════════════
cat("\n[5/5] 시각화...\n")

# Sample plot
pdf(file.path(FIGDIR, "Fig_DIABLO_sample.pdf"), width=12, height=5)
par(mfrow=c(1,3))
for (nm in names(X)) {
  plotIndiv(diablo_final,
            ind.names=FALSE,
            group=Y,
            blocks=nm,
            comp=c(1,2),
            legend=TRUE,
            title=paste("DIABLO —", nm),
            col.per.group=c("Control"="#4CAF50","Periodontitis"="#F44336"))
}
dev.off()

# Circos plot (오믹스 간 상관관계)
# correlation cutoff = 0.7 (문헌 기준)
pdf(file.path(FIGDIR, "Fig_DIABLO_circos.pdf"), width=10, height=10)
tryCatch({
  circosPlot(diablo_final,
             cutoff=0.7,
             line=FALSE,
             color.blocks=c("#2196F3","#FF5722","#4CAF50"),
             color.cor=c("darkred","darkblue"),
             size.labels=1.2)
  title("DIABLO Circos Plot\n(Cross-Omics Correlation, |r| > 0.7)")
}, error=function(e) {
  cat(sprintf("  circosPlot 오류: %s\n", e$message))
  plot.new()
  text(0.5, 0.5, "circosPlot 실행 실패\n(피처 수 부족 가능)", cex=1.2)
})
dev.off()

# Variable importance plot (각 오믹스별)
pdf(file.path(FIGDIR, "Fig_DIABLO_variables.pdf"), width=14, height=5)
par(mfrow=c(1,3))
for (nm in names(X)) {
  plotLoadings(diablo_final,
               comp=1,
               block=nm,
               method="median",
               contrib="max",
               title=paste("Loadings —", nm, "(Comp 1)"))
}
dev.off()

# Arrow plot (샘플별 오믹스 일치도)
pdf(file.path(FIGDIR, "Fig_DIABLO_arrow.pdf"), width=8, height=7)
tryCatch({
  plotArrow(diablo_final,
            ind.names=FALSE,
            group=Y,
            col.per.group=c("Control"="#4CAF50","Periodontitis"="#F44336"),
            title="DIABLO Arrow Plot\n(Omics Consistency per Sample)")
}, error=function(e) {
  cat(sprintf("  Arrow plot 오류: %s\n", e$message))
})
dev.off()

# Performance 평가 (LOO-CV)
cat("\n  성능 평가 (LOO-CV)...\n")
tryCatch({
  perf_diablo <- perf(
    diablo_final,
    validation = "loo",
    nrepeat    = 1,
    dist       = c("centroids.dist", "max.dist")
  )

  cat("  BER per component:\n")
  print(perf_diablo$MajorityVote.error.rate)

  write.csv(
    as.data.frame(perf_diablo$MajorityVote.error.rate),
    file.path(OUTDIR, "diablo_loocv_BER.csv")
  )
}, error=function(e) {
  cat(sprintf("  성능 평가 오류: %s\n", e$message))
})

# 결과 요약 저장
summary_diablo <- data.frame(
  omics          = names(selected_features),
  n_features_sel = sapply(selected_features, length)
)
write.csv(summary_diablo, file.path(OUTDIR, "diablo_summary.csv"), row.names=FALSE)

cat("\n✅ DIABLO 완료!\n")
