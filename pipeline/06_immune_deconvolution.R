#!/usr/bin/env Rscript
# =============================================================
# 06_immune_deconvolution.R
# 면역세포 디컨볼루션 + 바이오마커 연계 분석
# CIBERSORT + xCell → 치주염 면역 microenvironment
# =============================================================

suppressPackageStartupMessages({
  library(immunedeconv)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(RColorBrewer)
  library(patchwork)
  library(pheatmap)
})

set.seed(42)

OUTDIR <- "../results/immune"
FIGDIR <- "../figures"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

cat("=" , rep("=",58), "\n", sep="")
cat("06. Immune Cell Deconvolution Analysis\n")
cat("=" , rep("=",58), "\n\n", sep="")


# ══════════════════════════════════════════════════════════════
# STEP 1: TPM 데이터 준비
# ══════════════════════════════════════════════════════════════
cat("[1/5] 데이터 로드...\n")

meta <- read.csv("../data/processed/sample_metadata.csv", row.names=1)
meta$group <- factor(meta$group, levels=c("Control","Periodontitis"))

# TPM 우선, 없으면 CPM 사용
tpm_path <- "../data/processed/mrna_tpm.csv"
cpm_path <- "../results/mrna/mrna_cpm.csv"

if (file.exists(tpm_path)) {
  expr_mat <- as.matrix(read.csv(tpm_path, row.names=1, check.names=FALSE))
  cat("  TPM 행렬 사용\n")
} else if (file.exists(cpm_path)) {
  expr_mat <- as.matrix(read.csv(cpm_path, row.names=1, check.names=FALSE))
  cat("  CPM 행렬 사용 (TPM 없음)\n")
} else {
  # count matrix → CPM 직접 계산
  counts <- read.csv("../data/processed/mrna_counts.csv", row.names=1)
  expr_mat <- sweep(as.matrix(counts), 2, colSums(counts), "/") * 1e6
  cat("  Count → CPM 변환\n")
}

# 샘플 맞추기
common_s  <- intersect(colnames(expr_mat), rownames(meta))
expr_mat  <- expr_mat[, common_s]
meta      <- meta[common_s, , drop=FALSE]

cat(sprintf("  유전자: %d개, 샘플: %d개\n", nrow(expr_mat), ncol(expr_mat)))


# ══════════════════════════════════════════════════════════════
# STEP 2: 디컨볼루션 실행
# ══════════════════════════════════════════════════════════════
cat("\n[2/5] 면역세포 디컨볼루션 실행...\n")

# 방법별 실행 함수
run_deconv <- function(method, mat) {
  tryCatch({
    cat(sprintf("  %s 실행 중...\n", method))
    res <- deconvolute(mat, method)
    cat(sprintf("  %s 완료: %d개 세포 유형\n", method, nrow(res)))
    return(res)
  }, error = function(e) {
    cat(sprintf("  %s 오류: %s\n", method, e$message))
    return(NULL)
  })
}

res_cibersort <- run_deconv("cibersort_abs", expr_mat)
res_xcell     <- run_deconv("xcell", expr_mat)
res_quantiseq <- run_deconv("quantiseq", expr_mat)   # 오픈소스 대안

# 결과 저장
results_all <- list()
if (!is.null(res_cibersort)) {
  write.csv(res_cibersort, file.path(OUTDIR, "cibersort_abs_results.csv"), row.names=FALSE)
  results_all[["CIBERSORT"]] <- res_cibersort
}
if (!is.null(res_xcell)) {
  write.csv(res_xcell, file.path(OUTDIR, "xcell_results.csv"), row.names=FALSE)
  results_all[["xCell"]] <- res_xcell
}
if (!is.null(res_quantiseq)) {
  write.csv(res_quantiseq, file.path(OUTDIR, "quantiseq_results.csv"), row.names=FALSE)
  results_all[["quanTIseq"]] <- res_quantiseq
}


# ══════════════════════════════════════════════════════════════
# STEP 3: 통계 분석 (그룹 간 비교)
# ══════════════════════════════════════════════════════════════
cat("\n[3/5] 그룹 간 비교 (Wilcoxon test)...\n")

analyze_deconv_results <- function(res_df, method_name) {
  if (is.null(res_df)) return(NULL)

  long_df <- res_df %>%
    pivot_longer(-cell_type, names_to="sample", values_to="fraction") %>%
    left_join(meta %>% rownames_to_column("sample"), by="sample") %>%
    filter(!is.na(fraction))

  # 각 세포 유형별 Wilcoxon test
  stat_results <- long_df %>%
    group_by(cell_type) %>%
    summarise(
      ctrl_mean  = mean(fraction[group=="Control"],      na.rm=TRUE),
      perio_mean = mean(fraction[group=="Periodontitis"], na.rm=TRUE),
      p_value    = tryCatch(
        wilcox.test(
          fraction[group=="Control"],
          fraction[group=="Periodontitis"]
        )$p.value,
        error=function(e) NA_real_
      )
    ) %>%
    mutate(
      log2FC = log2((perio_mean + 1e-6) / (ctrl_mean + 1e-6)),
      direction = ifelse(log2FC > 0, "Up in Periodontitis", "Down in Periodontitis")
    ) %>%
    arrange(p_value)

  stat_results$padj <- p.adjust(stat_results$p_value, method="BH")

  sig_cells <- stat_results %>% filter(padj < 0.05)
  cat(sprintf("  [%s] 유의한 세포 유형: %d개\n", method_name, nrow(sig_cells)))

  write.csv(long_df, file.path(OUTDIR, paste0(method_name, "_long.csv")), row.names=FALSE)
  write.csv(stat_results, file.path(OUTDIR, paste0(method_name, "_stats.csv")), row.names=FALSE)

  return(list(long=long_df, stats=stat_results))
}

analysis_results <- lapply(names(results_all), function(nm) {
  analyze_deconv_results(results_all[[nm]], nm)
})
names(analysis_results) <- names(results_all)


# ══════════════════════════════════════════════════════════════
# STEP 4: 시각화
# ══════════════════════════════════════════════════════════════
cat("\n[4/5] 시각화...\n")

# 주요 면역세포 목록 (치주염 관련)
key_cells_patterns <- c(
  "T cell", "B cell", "NK", "Macrophage", "Neutrophil",
  "Dendritic", "Monocyte", "Treg", "CD4", "CD8"
)

plot_deconv_violin <- function(long_df, stats_df, method_name) {
  if (is.null(long_df)) return(NULL)

  # 주요 세포만 필터
  key_filter <- grepl(
    paste(key_cells_patterns, collapse="|"),
    long_df$cell_type, ignore.case=TRUE
  )
  plot_df <- long_df[key_filter, ]

  if (nrow(plot_df) == 0) return(NULL)

  # 유의성 표시
  if (!is.null(stats_df)) {
    sig_stars <- stats_df %>%
      filter(cell_type %in% unique(plot_df$cell_type)) %>%
      mutate(
        star = case_when(
          padj < 0.001 ~ "***",
          padj < 0.01  ~ "**",
          padj < 0.05  ~ "*",
          TRUE         ~ "ns"
        )
      )
  }

  p <- ggplot(plot_df, aes(x=cell_type, y=fraction+1e-6, fill=group)) +
    geom_violin(alpha=0.7, position=position_dodge(0.8), trim=TRUE) +
    geom_boxplot(width=0.15, position=position_dodge(0.8),
                 outlier.shape=NA, alpha=0.9) +
    geom_jitter(aes(group=group),
                position=position_jitterdodge(jitter.width=0.1, dodge.width=0.8),
                size=1.5, alpha=0.5) +
    scale_fill_manual(values=c("Control"="#4CAF50","Periodontitis"="#F44336")) +
    scale_y_log10() +
    labs(title=paste("Immune Cell Composition —", method_name),
         subtitle="Control vs Periodontitis (log10 scale)",
         x="", y="Estimated Fraction (log10)") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
          legend.position="top")

  return(p)
}

for (nm in names(analysis_results)) {
  if (!is.null(analysis_results[[nm]])) {
    p <- plot_deconv_violin(
      analysis_results[[nm]]$long,
      analysis_results[[nm]]$stats,
      nm
    )
    if (!is.null(p)) {
      pdf(file.path(FIGDIR, paste0("Fig_Immune_", nm, "_violin.pdf")),
          width=14, height=7)
      print(p)
      dev.off()
    }
  }
}

# Stacked bar plot (구성 비율)
if (!is.null(res_cibersort)) {
  long_ciber <- analysis_results[["CIBERSORT"]]$long

  # 상위 10개 세포 유형만
  top_cells <- long_ciber %>%
    group_by(cell_type) %>%
    summarise(mean_frac=mean(fraction, na.rm=TRUE)) %>%
    top_n(10, mean_frac) %>%
    pull(cell_type)

  stack_df <- long_ciber %>%
    filter(cell_type %in% top_cells) %>%
    group_by(sample, cell_type, group) %>%
    summarise(fraction=mean(fraction, na.rm=TRUE)) %>%
    arrange(group, sample)

  p_stack <- ggplot(stack_df, aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat="identity") +
    facet_grid(~group, scales="free_x", space="free_x") +
    scale_fill_brewer(palette="Set3") +
    labs(title="Immune Cell Composition (CIBERSORT)",
         x="Sample", y="Estimated Fraction", fill="Cell Type") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=90, hjust=1, size=7),
          strip.background=element_rect(fill="lightgray"))

  pdf(file.path(FIGDIR, "Fig_Immune_stacked_bar.pdf"), width=14, height=7)
  print(p_stack)
  dev.off()
}


# ══════════════════════════════════════════════════════════════
# STEP 5: 면역세포 비율 vs 바이오마커 상관 분석
# ══════════════════════════════════════════════════════════════
cat("\n[5/5] 면역세포 vs 바이오마커 상관...\n")

# 이전 결과에서 바이오마커 로드
biomarker_path <- "../results/ml/final_biomarker_panel.csv"
hub_gene_path  <- "../results/mrna/WGCNA_hub_genes.csv"

if (!is.null(res_cibersort) &&
    (file.exists(biomarker_path) || file.exists(hub_gene_path))) {

  # DEG 발현값
  vst_mat <- read.csv("../results/mrna/mrna_vst_matrix.csv", row.names=1)

  # CIBERSORT 결과를 wide format으로
  ciber_wide <- res_cibersort %>%
    column_to_rownames("cell_type") %>%
    t() %>% as.data.frame()

  common_s2 <- intersect(rownames(ciber_wide), colnames(vst_mat))
  if (length(common_s2) >= 20) {

    # 유의 면역세포 목록
    sig_cells_list <- analysis_results[["CIBERSORT"]]$stats %>%
      filter(padj < 0.05) %>%
      pull(cell_type)

    if (length(sig_cells_list) > 0) {
      # 상관 행렬 계산
      hub_genes <- if (file.exists(hub_gene_path)) {
        read.csv(hub_gene_path)$gene[1:min(20, nrow(read.csv(hub_gene_path)))]
      } else {
        character(0)
      }

      valid_hubs <- intersect(hub_genes, rownames(vst_mat))
      if (length(valid_hubs) >= 5 && length(sig_cells_list) >= 3) {

        gene_expr   <- t(vst_mat[valid_hubs, common_s2])
        immune_frac <- ciber_wide[common_s2, sig_cells_list[1:min(8, length(sig_cells_list))]]

        corr_mat <- cor(gene_expr, immune_frac, use="pairwise.complete.obs",
                        method="spearman")

        pdf(file.path(FIGDIR, "Fig_Gene_ImmuneCell_Corr.pdf"), width=10, height=8)
        pheatmap(
          corr_mat,
          color       = colorRampPalette(c("#2196F3","white","#F44336"))(100),
          breaks      = seq(-1,1,length.out=101),
          main        = "Hub Gene – Immune Cell Correlation\n(Spearman r)",
          fontsize    = 10,
          display_numbers=TRUE,
          number_format="%.2f",
          number_color="black"
        )
        dev.off()

        write.csv(corr_mat,
                  file.path(OUTDIR, "hub_gene_immune_correlation.csv"))
        cat(sprintf("  상관 분석 완료: %d genes × %d cell types\n",
                    ncol(gene_expr), ncol(immune_frac)))
      }
    }
  }
}

cat("\n✅ Immune Deconvolution 완료!\n")
