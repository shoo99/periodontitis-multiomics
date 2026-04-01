#!/usr/bin/env Rscript
# =============================================================
# 02_proteomics_preprocessing.R
# Proteome Discoverer 결과 전처리 + limma DEP + WGCNA
# 입력: PD normalized abundance (protein × sample)
# 출력: DEP 결과, log2 matrix, WGCNA hub protein
# =============================================================

suppressPackageStartupMessages({
  library(limma)
  library(WGCNA)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(impute)       # KNN imputation
  library(pcaMethods)   # QRILC imputation
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(biomaRt)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)
set.seed(42)

INDIR    <- "../data/processed"
OUTDIR   <- "../results/proteomics"
FIGDIR   <- "../figures"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

cat("=" , rep("=",58), "\n", sep="")
cat("02. Proteomics Preprocessing & Analysis\n")
cat("=" , rep("=",58), "\n\n", sep="")


# ══════════════════════════════════════════════════════════════
# STEP 1: 데이터 로드
# ══════════════════════════════════════════════════════════════
cat("[1/5] 데이터 로드...\n")

prot <- read.csv(file.path(INDIR, "proteomics_abundance.csv"),
                 row.names=1, check.names=FALSE)
meta <- read.csv(file.path(INDIR, "sample_metadata.csv"), row.names=1)
meta <- meta[colnames(prot), , drop=FALSE]
meta$group <- factor(meta$group, levels=c("Control","Periodontitis"))

cat(sprintf("  단백질: %d개, 샘플: %d개\n", nrow(prot), ncol(prot)))

# 값이 0인 것 → NA 처리
prot[prot == 0] <- NA
prot_mat <- as.matrix(prot)


# ══════════════════════════════════════════════════════════════
# STEP 2: 결측치 처리 + 정규화
# ══════════════════════════════════════════════════════════════
cat("\n[2/5] 결측치 처리 + 정규화...\n")

# 결측 비율 계산
missing_ratio <- rowMeans(is.na(prot_mat))
cat(sprintf("  전체 결측률: %.1f%%\n", mean(missing_ratio)*100))

# 필터: 그룹별 50% 이상 결측 단백질 제거
ctrl_idx  <- which(meta$group == "Control")
perio_idx <- which(meta$group == "Periodontitis")

keep_ctrl  <- rowMeans(is.na(prot_mat[, ctrl_idx]))  < 0.5
keep_perio <- rowMeans(is.na(prot_mat[, perio_idx])) < 0.5
keep_prot  <- keep_ctrl | keep_perio   # 어느 한 그룹에서라도 충분하면 유지
prot_filt  <- prot_mat[keep_prot, ]
cat(sprintf("  필터 후 단백질: %d개\n", nrow(prot_filt)))

# Log2 변환
prot_log2 <- log2(prot_filt)

# 결측치 패턴 분류 및 대체
# MNAR (Missing Not At Random): MinProb imputation
# MAR  (Missing At Random)    : KNN imputation
# 단순 전략: 결측 비율 < 20% → KNN, ≥ 20% → MinProb
missing_per_protein <- rowMeans(is.na(prot_log2))

# KNN imputation (MAR)
prot_knn <- impute.knn(prot_log2, k=5)$data

# MinProb imputation (MNAR) — 하위 1% 값 기반 정규분포에서 샘플링
min_prob_impute <- function(mat, q=0.01, tune_sigma=0.3) {
  mat_imp <- mat
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    nas <- is.na(col)
    if (sum(!nas) < 3) next
    mu    <- quantile(col[!nas], q, na.rm=TRUE)
    sigma <- sd(col[!nas], na.rm=TRUE) * tune_sigma
    mat_imp[nas, j] <- rnorm(sum(nas), mean=mu, sd=sigma)
  }
  return(mat_imp)
}

# 통합 대체: 결측 > 20%인 단백질은 MinProb, 나머지는 KNN
prot_imputed <- prot_knn
high_missing <- missing_per_protein >= 0.2
if (any(high_missing)) {
  prot_minprob <- min_prob_impute(prot_log2[high_missing, ])
  prot_imputed[high_missing, ] <- prot_minprob
}

# Median centering (배치 간 보정)
median_per_sample <- apply(prot_imputed, 2, median, na.rm=TRUE)
global_median      <- median(prot_imputed, na.rm=TRUE)
prot_norm <- sweep(prot_imputed, 2, median_per_sample - global_median, "-")

cat(sprintf("  Imputation + normalization 완료\n"))
write.csv(prot_norm, file.path(OUTDIR, "proteomics_log2_normalized.csv"))

# PCA
pca_res <- prcomp(t(prot_norm), scale.=TRUE)
pca_df  <- data.frame(
  PC1   = pca_res$x[,1],
  PC2   = pca_res$x[,2],
  group = meta$group,
  label = rownames(meta)
)
pca_var <- round(100 * summary(pca_res)$importance[2, 1:2])

pdf(file.path(FIGDIR, "Fig_Prot_PCA.pdf"), width=7, height=6)
ggplot(pca_df, aes(PC1, PC2, color=group, label=label)) +
  geom_point(size=4, alpha=0.8) +
  geom_text_repel(size=3, max.overlaps=10) +
  scale_color_manual(values=c("Control"="#4CAF50","Periodontitis"="#FF5722")) +
  labs(title="Proteomics PCA (log2 normalized)",
       x=paste0("PC1: ",pca_var[1],"%"),
       y=paste0("PC2: ",pca_var[2],"%")) +
  theme_bw(base_size=13)
dev.off()


# ══════════════════════════════════════════════════════════════
# STEP 3: limma DEP 분석
# ══════════════════════════════════════════════════════════════
cat("\n[3/5] limma DEP 분석...\n")

design  <- model.matrix(~ group, data=meta)
fit     <- lmFit(prot_norm, design)
fit     <- eBayes(fit)

# 모든 단백질 결과
dep_all <- topTable(fit, coef="groupPeriodontitis",
                    n=Inf, sort.by="P") %>%
  rownames_to_column("protein")

# 유의 DEP
dep_sig <- dep_all %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.58) %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down"))

cat(sprintf("  유의 DEP: %d개 (상향 %d, 하향 %d)\n",
            nrow(dep_sig),
            sum(dep_sig$direction=="Up"),
            sum(dep_sig$direction=="Down")))

write.csv(dep_all, file.path(OUTDIR, "limma_all_results.csv"), row.names=FALSE)
write.csv(dep_sig, file.path(OUTDIR, "limma_significant_DEP.csv"), row.names=FALSE)

# Volcano
dep_plot <- dep_all %>%
  mutate(
    sig   = case_when(
      adj.P.Val < 0.05 & logFC >  0.58 ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.58 ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(sig!="NS" & abs(logFC)>1.5 & -log10(adj.P.Val)>5, protein, "")
  )

pdf(file.path(FIGDIR, "Fig_Prot_Volcano.pdf"), width=8, height=7)
ggplot(dep_plot, aes(logFC, -log10(adj.P.Val), color=sig, label=label)) +
  geom_point(alpha=0.5, size=1.5) +
  geom_text_repel(size=3, max.overlaps=15, color="black") +
  scale_color_manual(values=c("Up"="#FF5722","Down"="#FF9800","NS"="#9E9E9E")) +
  geom_vline(xintercept=c(-0.58,0.58), linetype="dashed", color="gray40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray40") +
  labs(title="Volcano Plot — Proteomics (Periodontitis vs Control)",
       x="log2 Fold Change", y="-log10(adjusted p-value)") +
  theme_bw(base_size=13)
dev.off()

# Heatmap (Top 40 DEP)
top_dep <- dep_sig %>% arrange(adj.P.Val) %>% head(40) %>% pull(protein)
if (length(top_dep) >= 10) {
  heat_prot <- prot_norm[top_dep, ]
  heat_prot_z <- t(scale(t(heat_prot)))
  ann_col <- data.frame(Group=meta$group, row.names=colnames(heat_prot))
  ann_colors <- list(Group=c(Control="#4CAF50",Periodontitis="#FF5722"))

  pdf(file.path(FIGDIR, "Fig_Prot_Heatmap_top40.pdf"), width=12, height=9)
  pheatmap(heat_prot_z, annotation_col=ann_col, annotation_colors=ann_colors,
           show_colnames=FALSE,
           color=colorRampPalette(c("#FF9800","white","#FF5722"))(100),
           breaks=seq(-3,3,length.out=101),
           main="Top 40 DEP Heatmap (Z-score)")
  dev.off()
}


# ══════════════════════════════════════════════════════════════
# STEP 4: WGCNA — 단백질 공동발현 네트워크
# ══════════════════════════════════════════════════════════════
cat("\n[4/5] WGCNA (단백체)...\n")

# 상위 변동 단백질
mad_prot <- apply(prot_norm, 1, mad)
top_prots <- names(sort(mad_prot, decreasing=TRUE))[
  1:min(2000, nrow(prot_norm))]

datExpr_prot <- t(prot_norm[top_prots, ])

# 품질 확인
gsg_prot <- goodSamplesGenes(datExpr_prot, verbose=0)
if (!gsg_prot$allOK) {
  datExpr_prot <- datExpr_prot[gsg_prot$goodSamples, gsg_prot$goodGenes]
}

# Soft threshold (단백체는 낮게)
powers_p <- c(1:15)
sft_prot <- pickSoftThreshold(datExpr_prot, powerVector=powers_p,
                               networkType="signed hybrid",
                               RsquaredCut=0.80, verbose=0)
power_prot <- sft_prot$powerEstimate
if (is.na(power_prot)) {
  power_prot <- sft_prot$fitIndices[which.max(sft_prot$fitIndices[,"SFT.R.sq"]), "Power"]
}
cat(sprintf("  단백체 soft power: %d\n", power_prot))

# 네트워크 구축 (단백체 전용 파라미터)
net_prot <- blockwiseModules(
  datExpr_prot,
  power          = power_prot,
  networkType    = "signed hybrid",
  TOMType        = "signed",
  minModuleSize  = 15,       # 단백 수 적으므로 낮게
  mergeCutHeight = 0.25,
  deepSplit      = 2,
  corType        = "bicor",  # ⭐ 단백체: biweight midcorrelation
  maxPOutliers   = 0.05,
  verbose        = 0
)

# Module-Trait 상관 (치주염)
MEs_prot    <- orderMEs(net_prot$MEs)
trait_prot  <- data.frame(
  Periodontitis = as.numeric(meta$group == "Periodontitis"),
  row.names     = rownames(datExpr_prot)
)
mtcor_prot <- cor(MEs_prot, trait_prot, use="p")
mtpval_prot <- corPvalueStudent(mtcor_prot, nrow(datExpr_prot))

sig_mod_prot <- rownames(mtcor_prot)[
  abs(mtcor_prot[, "Periodontitis"]) > 0.5 &
  mtpval_prot[, "Periodontitis"]     < 0.05
]
cat(sprintf("  단백체 치주염 관련 모듈: %s\n",
            paste(sig_mod_prot, collapse=", ")))

# Hub protein 추출
hub_prots_all <- data.frame()
for (mod in sig_mod_prot) {
  mod_name  <- sub("ME", "", mod)
  mod_prots <- names(net_prot$colors)[net_prot$colors == mod_name]
  kME_p <- cor(datExpr_prot[, mod_prots], MEs_prot[, mod], use="p")
  GS_p  <- cor(datExpr_prot[, mod_prots], trait_prot$Periodontitis, use="p")

  hub_df <- data.frame(
    protein = mod_prots,
    module  = mod_name,
    kME     = kME_p[,1],
    GS      = GS_p[,1]
  ) %>% filter(abs(kME) > 0.7, abs(GS) > 0.3) %>%  # 단백체: 0.7
    arrange(desc(abs(kME)))

  hub_prots_all <- rbind(hub_prots_all, hub_df)
}
cat(sprintf("  Hub proteins: %d개\n", nrow(hub_prots_all)))

write.csv(hub_prots_all, file.path(OUTDIR, "WGCNA_hub_proteins.csv"), row.names=FALSE)
write.csv(as.data.frame(mtcor_prot), file.path(OUTDIR, "WGCNA_moduleTraitCor_prot.csv"))


# ══════════════════════════════════════════════════════════════
# STEP 5: mRNA-Protein Concordance 분석
# ══════════════════════════════════════════════════════════════
cat("\n[5/5] mRNA-Protein Concordance 분석...\n")

# DEG & DEP 결과 합치기 (gene symbol 기준)
mrna_res_path <- file.path("../results/mrna", "DESeq2_all_results.csv")

if (file.exists(mrna_res_path)) {
  mrna_res <- read.csv(mrna_res_path)

  # Protein → Gene symbol 매핑 (UniProt accession 제거, gene symbol만 사용 가정)
  # PD output에서 protein 이름이 gene symbol인 경우
  common_genes <- intersect(mrna_res$gene, dep_all$protein)
  cat(sprintf("  공통 유전자/단백질: %d개\n", length(common_genes)))

  if (length(common_genes) >= 10) {
    conc_df <- data.frame(
      gene       = common_genes,
      mrna_lfc   = mrna_res$log2FoldChange[match(common_genes, mrna_res$gene)],
      prot_lfc   = dep_all$logFC[match(common_genes, dep_all$protein)],
      mrna_padj  = mrna_res$padj[match(common_genes, mrna_res$gene)],
      prot_padj  = dep_all$adj.P.Val[match(common_genes, dep_all$protein)]
    ) %>%
      mutate(
        type = case_when(
          mrna_padj < 0.05 & prot_padj < 0.05 &
            sign(mrna_lfc) == sign(prot_lfc) ~ "Concordant",
          mrna_padj < 0.05 & prot_padj < 0.05 &
            sign(mrna_lfc) != sign(prot_lfc) ~ "Discordant",
          mrna_padj < 0.05 & prot_padj >= 0.05 ~ "mRNA only",
          mrna_padj >= 0.05 & prot_padj < 0.05 ~ "Protein only",
          TRUE ~ "NS"
        )
      )

    write.csv(conc_df, file.path(OUTDIR, "mRNA_Protein_concordance.csv"), row.names=FALSE)
    cat(sprintf("  Concordant: %d, Discordant: %d, mRNA-only: %d, Prot-only: %d\n",
                sum(conc_df$type=="Concordant"),
                sum(conc_df$type=="Discordant"),
                sum(conc_df$type=="mRNA only"),
                sum(conc_df$type=="Protein only")))

    # 4분면 산점도
    conc_plot <- conc_df %>% filter(type != "NS") %>%
      mutate(label = ifelse(type == "Concordant" & abs(mrna_lfc) > 1.5, gene, ""))

    pdf(file.path(FIGDIR, "Fig_mRNA_Protein_Concordance.pdf"), width=8, height=7)
    p <- ggplot(conc_plot, aes(mrna_lfc, prot_lfc, color=type, label=label)) +
      geom_point(alpha=0.7, size=2.5) +
      geom_text_repel(size=3, max.overlaps=12, color="black") +
      scale_color_manual(values=c(
        "Concordant" ="#9C27B0",
        "Discordant" ="#FF5722",
        "mRNA only"  ="#2196F3",
        "Protein only"="#4CAF50"
      )) +
      geom_hline(yintercept=0, linetype="dashed", color="gray50") +
      geom_vline(xintercept=0, linetype="dashed", color="gray50") +
      labs(title="mRNA vs Protein log2FC (Periodontitis vs Control)",
           x="mRNA log2FC", y="Protein log2FC", color="Type") +
      theme_bw(base_size=13) +
      annotate("text", x=-Inf, y=Inf,
               label=paste0("Concordant: ", sum(conc_df$type=="Concordant"),
                             "\nDiscordant: ", sum(conc_df$type=="Discordant")),
               hjust=0, vjust=1, size=4)
    print(p)
    dev.off()
  }
}

cat("\n✅ Proteomics Analysis 완료!\n")
cat(sprintf("  DEP: %d개, Hub proteins: %d개\n",
            nrow(dep_sig), nrow(hub_prots_all)))
