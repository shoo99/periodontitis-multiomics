#!/usr/bin/env Rscript
# =============================================================
# 01_mrna_preprocessing.R
# mRNA 전처리: QC → DESeq2 DEG → WGCNA 모듈 분석
# 입력: count matrix (gene × sample), sample metadata
# 출력: DEG 결과, VST matrix, WGCNA 모듈/hub gene
# =============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(GSVA)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)
set.seed(42)

# ── 경로 설정 ──────────────────────────────────────────────
INDIR    <- "../data/processed"
OUTDIR   <- "../results/mrna"
FIGDIR   <- "../figures"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("=" , rep("=",58), "\n", sep="")
cat("01. mRNA Preprocessing & Analysis\n")
cat("=" , rep("=",58), "\n\n", sep="")


# ══════════════════════════════════════════════════════════════
# STEP 1: 데이터 로드 및 QC
# ══════════════════════════════════════════════════════════════
cat("[1/6] 데이터 로드...\n")

counts <- read.csv(file.path(INDIR, "mrna_counts.csv"),
                   row.names=1, check.names=FALSE)
meta   <- read.csv(file.path(INDIR, "sample_metadata.csv"),
                   row.names=1)

# 샘플 순서 맞추기
meta <- meta[colnames(counts), , drop=FALSE]
meta$group <- factor(meta$group, levels=c("Control", "Periodontitis"))

cat(sprintf("  유전자: %d개, 샘플: %d개\n", nrow(counts), ncol(counts)))
cat(sprintf("  Control: %d, Periodontitis: %d\n",
            sum(meta$group=="Control"), sum(meta$group=="Periodontitis")))

# 저발현 유전자 제거 (CPM > 1 in ≥ 20% 샘플)
cpm_mat  <- sweep(counts, 2, colSums(counts), "/") * 1e6
keep_genes <- rowSums(cpm_mat > 1) >= (0.2 * ncol(counts))
counts_filt <- counts[keep_genes, ]
cat(sprintf("  필터 후 유전자: %d개\n", nrow(counts_filt)))


# ══════════════════════════════════════════════════════════════
# STEP 2: DESeq2 DEG 분석
# ══════════════════════════════════════════════════════════════
cat("\n[2/6] DESeq2 DEG 분석...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = meta,
  design    = ~ group
)

# VST 정규화 (시각화 + WGCNA 입력용)
vsd <- vst(dds, blind=FALSE)
vst_mat <- assay(vsd)

# PCA — Outlier 확인
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"))

pdf(file.path(FIGDIR, "Fig_mRNA_PCA.pdf"), width=7, height=6)
ggplot(pca_data, aes(PC1, PC2, color=group, label=name)) +
  geom_point(size=4, alpha=0.8) +
  geom_text_repel(size=3, max.overlaps=10) +
  scale_color_manual(values=c("Control"="#4CAF50","Periodontitis"="#F44336")) +
  labs(title="mRNA PCA (VST)",
       x=paste0("PC1: ",pca_var[1],"% variance"),
       y=paste0("PC2: ",pca_var[2],"% variance")) +
  theme_bw(base_size=13)
dev.off()

# DESeq2 실행
dds <- DESeq(dds)

# 결과 추출 (Periodontitis vs Control)
res <- results(dds,
               contrast=c("group","Periodontitis","Control"),
               alpha=0.05)
res <- lfcShrink(dds, coef="group_Periodontitis_vs_Control",
                 type="apeglm", res=res)

# 유의 DEG
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

deg_sig <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1.0) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down"))

cat(sprintf("  유의 DEG: %d개 (상향 %d, 하향 %d)\n",
            nrow(deg_sig),
            sum(deg_sig$direction=="Up"),
            sum(deg_sig$direction=="Down")))

write.csv(res_df,  file.path(OUTDIR, "DESeq2_all_results.csv"), row.names=FALSE)
write.csv(deg_sig, file.path(OUTDIR, "DESeq2_significant_DEG.csv"), row.names=FALSE)
write.csv(vst_mat, file.path(OUTDIR, "mrna_vst_matrix.csv"))


# ── Volcano Plot ─────────────────────────────────────────────
res_plot <- res_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1.0 ~ "Up",
      padj < 0.05 & log2FoldChange < -1.0 ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(sig != "NS" & abs(log2FoldChange) > 2.5 &
                   -log10(padj) > 10, gene, "")
  )

pdf(file.path(FIGDIR, "Fig_mRNA_Volcano.pdf"), width=8, height=7)
ggplot(res_plot, aes(log2FoldChange, -log10(padj), color=sig, label=label)) +
  geom_point(alpha=0.5, size=1.5) +
  geom_text_repel(size=3, max.overlaps=15, color="black") +
  scale_color_manual(values=c("Up"="#F44336","Down"="#2196F3","NS"="#9E9E9E")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="gray40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray40") +
  labs(title="Volcano Plot — mRNA (Periodontitis vs Control)",
       x="log2 Fold Change", y="-log10(adjusted p-value)",
       color="Direction") +
  theme_bw(base_size=13)
dev.off()

# ── Heatmap (Top 50 DEG) ─────────────────────────────────────
top_deg <- deg_sig %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

if (length(top_deg) >= 10) {
  heat_mat <- vst_mat[top_deg, ]
  heat_mat_z <- t(scale(t(heat_mat)))  # Z-score per gene

  ann_col <- data.frame(
    Group = meta$group,
    row.names = colnames(heat_mat_z)
  )
  ann_colors <- list(Group=c(Control="#4CAF50", Periodontitis="#F44336"))

  pdf(file.path(FIGDIR, "Fig_mRNA_Heatmap_top50.pdf"), width=12, height=10)
  pheatmap(heat_mat_z,
           annotation_col=ann_col,
           annotation_colors=ann_colors,
           show_colnames=FALSE,
           color=colorRampPalette(c("#2196F3","white","#F44336"))(100),
           breaks=seq(-3,3,length.out=101),
           main="Top 50 DEG Heatmap (Z-score)")
  dev.off()
}


# ══════════════════════════════════════════════════════════════
# STEP 3: GSEA/ORA Pathway 분석
# ══════════════════════════════════════════════════════════════
cat("\n[3/6] Pathway 분석 (GSEA + ORA)...\n")

# ENTREZ ID 변환
gene_entrez <- bitr(res_df$gene, fromType="SYMBOL",
                    toType="ENTREZID", OrgDb=org.Hs.eg.db)
res_entrez <- merge(res_df, gene_entrez, by.x="gene", by.y="SYMBOL")

# GSEA
gsea_input <- res_entrez %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(ENTREZID, log2FoldChange) %>%
  distinct(ENTREZID, .keep_all=TRUE)

gene_list <- setNames(gsea_input$log2FoldChange, gsea_input$ENTREZID)

# KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList     = gene_list,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH",
  seed         = 42
)

# GO GSEA
gsea_go <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH",
  seed         = 42
)

# ORA (DEG 상향만)
deg_up_entrez <- res_entrez %>%
  filter(padj < 0.05, log2FoldChange > 1.0) %>%
  pull(ENTREZID)

ora_kegg <- enrichKEGG(
  gene         = deg_up_entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH"
)

# 저장
write.csv(as.data.frame(gsea_kegg), file.path(OUTDIR, "GSEA_KEGG_results.csv"), row.names=FALSE)
write.csv(as.data.frame(gsea_go),   file.path(OUTDIR, "GSEA_GO_BP_results.csv"), row.names=FALSE)
write.csv(as.data.frame(ora_kegg),  file.path(OUTDIR, "ORA_KEGG_upDEG.csv"), row.names=FALSE)

# Dotplot
if (nrow(as.data.frame(gsea_kegg)) > 0) {
  pdf(file.path(FIGDIR, "Fig_mRNA_GSEA_KEGG.pdf"), width=10, height=8)
  print(dotplot(gsea_kegg, showCategory=20, split=".sign") +
          facet_grid(.~.sign) +
          ggtitle("KEGG GSEA — mRNA"))
  dev.off()
}


# ══════════════════════════════════════════════════════════════
# STEP 4: WGCNA — 유전자 공동발현 네트워크
# ══════════════════════════════════════════════════════════════
cat("\n[4/6] WGCNA 분석...\n")

# 입력: VST matrix (상위 변동 유전자 5000개)
mad_vals  <- apply(vst_mat, 1, mad)
top_genes <- names(sort(mad_vals, decreasing=TRUE))[1:min(5000, nrow(vst_mat))]
datExpr   <- t(vst_mat[top_genes, ])  # sample × gene

# 샘플 품질 확인
gsg <- goodSamplesGenes(datExpr, verbose=0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("  샘플/유전자 제거 후: %d × %d\n", nrow(datExpr), ncol(datExpr)))
}

# Soft threshold power 선택
powers <- c(1:20)
sft    <- pickSoftThreshold(datExpr, powerVector=powers,
                             networkType="signed hybrid",
                             RsquaredCut=0.85, verbose=0)

optimal_power <- sft$powerEstimate
if (is.na(optimal_power)) {
  # R² 기준 미충족 시 최대값 power 선택
  rsq <- sft$fitIndices[, "SFT.R.sq"]
  optimal_power <- sft$fitIndices[which.max(rsq), "Power"]
}
cat(sprintf("  선택된 soft power: %d\n", optimal_power))

# Soft threshold plot 저장
pdf(file.path(FIGDIR, "Fig_WGCNA_softpower.pdf"), width=9, height=5)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold", ylab="Scale Free Topology R²",
     main="Scale Independence", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.85, col="red", lty=2)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold", ylab="Mean Connectivity",
     main="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

# 네트워크 구축
net <- blockwiseModules(
  datExpr,
  power              = optimal_power,
  networkType        = "signed hybrid",
  TOMType            = "signed",
  minModuleSize      = 20,
  mergeCutHeight     = 0.25,
  deepSplit          = 2,
  corType            = "pearson",
  numericLabels      = FALSE,
  pamRespectsDendro  = FALSE,
  saveTOMs           = FALSE,
  verbose            = 0
)

module_colors <- net$colors
module_table  <- table(module_colors)
cat(sprintf("  발견된 모듈: %d개\n", length(unique(module_colors))-1))
cat("  모듈 크기:\n")
print(sort(module_table, decreasing=TRUE))

# Module-Trait 상관
MEs <- orderMEs(net$MEs)
trait_df <- data.frame(
  Periodontitis = as.numeric(meta$group == "Periodontitis"),
  row.names     = rownames(datExpr)
)

moduleTraitCor <- cor(MEs, trait_df, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 유의 모듈 (|r| > 0.5, p < 0.05)
sig_modules <- rownames(moduleTraitCor)[
  abs(moduleTraitCor[, "Periodontitis"]) > 0.5 &
  moduleTraitPvalue[, "Periodontitis"] < 0.05
]
cat(sprintf("\n  치주염 연관 유의 모듈: %s\n", paste(sig_modules, collapse=", ")))

# Module-Trait heatmap
pdf(file.path(FIGDIR, "Fig_WGCNA_moduleTrait.pdf"), width=7, height=9)
textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")"
)
par(mar=c(6,8.5,3,3))
labeledHeatmap(
  Matrix       = moduleTraitCor,
  xLabels      = colnames(trait_df),
  yLabels      = rownames(moduleTraitCor),
  ySymbols     = rownames(moduleTraitCor),
  colorLabels  = FALSE,
  colors       = blueWhiteRed(50),
  textMatrix   = textMatrix,
  setStdMargins= FALSE,
  cex.text     = 0.6,
  zlim         = c(-1,1),
  main         = "Module-Trait Correlation (mRNA)"
)
dev.off()

# Hub gene 추출 (유의 모듈에서)
hub_genes_all <- data.frame()
for (mod in sig_modules) {
  mod_name <- sub("ME", "", mod)
  mod_genes <- names(module_colors)[module_colors == mod_name]

  # kME (module membership)
  kME <- cor(datExpr[, mod_genes], MEs[, mod], use="p")
  # GS (gene significance with trait)
  GS  <- cor(datExpr[, mod_genes], trait_df$Periodontitis, use="p")

  hub_df <- data.frame(
    gene   = mod_genes,
    module = mod_name,
    kME    = kME[,1],
    GS     = GS[,1]
  ) %>% filter(abs(kME) > 0.8, abs(GS) > 0.3) %>%
    arrange(desc(abs(kME)))

  hub_genes_all <- rbind(hub_genes_all, hub_df)
}
cat(sprintf("  Hub genes (kME>0.8, GS>0.3): %d개\n", nrow(hub_genes_all)))

write.csv(hub_genes_all, file.path(OUTDIR, "WGCNA_hub_genes.csv"), row.names=FALSE)
write.csv(as.data.frame(moduleTraitCor),
          file.path(OUTDIR, "WGCNA_moduleTraitCor.csv"))


# ══════════════════════════════════════════════════════════════
# STEP 5: 면역세포 디컨볼루션 (immunedeconv)
# ══════════════════════════════════════════════════════════════
cat("\n[5/6] 면역세포 디컨볼루션...\n")

# TPM 매트릭스 필요 (VST가 아닌 TPM)
# TPM이 없으면 CPM 근사 사용
tpm_path <- file.path(INDIR, "mrna_tpm.csv")
if (file.exists(tpm_path)) {
  tpm_mat <- read.csv(tpm_path, row.names=1, check.names=FALSE)
  tpm_mat <- as.matrix(tpm_mat)
} else {
  # CPM으로 대체
  cat("  TPM 없음 → CPM 사용\n")
  tpm_mat <- sweep(as.matrix(counts_filt), 2, colSums(counts_filt), "/") * 1e6
}

tryCatch({
  library(immunedeconv)

  # CIBERSORT (절대 정량)
  cat("  CIBERSORT 실행 중...\n")
  res_cibersort <- deconvolute(tpm_mat, "cibersort_abs")

  # xCell
  cat("  xCell 실행 중...\n")
  res_xcell <- deconvolute(tpm_mat, "xcell")

  # 결과 저장
  write.csv(res_cibersort, file.path(OUTDIR, "immunedeconv_CIBERSORT.csv"), row.names=FALSE)
  write.csv(res_xcell,     file.path(OUTDIR, "immunedeconv_xCell.csv"),     row.names=FALSE)

  # 그룹별 비교 (Wilcoxon test)
  # 주요 세포 비교
  key_cells <- c("T cell CD4+", "T cell CD8+", "B cell", "NK cell",
                 "Macrophage", "Neutrophil", "Dendritic cell")

  cib_long <- res_cibersort %>%
    tidyr::pivot_longer(-cell_type, names_to="sample", values_to="fraction") %>%
    filter(cell_type %in% key_cells) %>%
    left_join(meta %>% rownames_to_column("sample"), by="sample")

  write.csv(cib_long, file.path(OUTDIR, "immunedeconv_CIBERSORT_long.csv"), row.names=FALSE)

  # 시각화
  pdf(file.path(FIGDIR, "Fig_ImmuneDeconv_violin.pdf"), width=12, height=8)
  p <- ggplot(cib_long, aes(x=cell_type, y=fraction, fill=group)) +
    geom_violin(alpha=0.7, position=position_dodge(0.8)) +
    geom_boxplot(width=0.1, position=position_dodge(0.8), outlier.shape=NA) +
    scale_fill_manual(values=c("Control"="#4CAF50","Periodontitis"="#F44336")) +
    labs(title="Immune Cell Fraction (CIBERSORT)",
         x="", y="Estimated Fraction") +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=40, hjust=1))
  print(p)
  dev.off()

}, error=function(e) {
  cat(sprintf("  immunedeconv 오류: %s\n", e$message))
  cat("  설치: install.packages('remotes')\n")
  cat("  remotes::install_github('omnideconv/immunedeconv')\n")
})


# ══════════════════════════════════════════════════════════════
# STEP 6: 결과 요약
# ══════════════════════════════════════════════════════════════
cat("\n[6/6] 결과 요약 저장...\n")

summary_list <- list(
  n_genes_input    = nrow(counts),
  n_genes_filtered = nrow(counts_filt),
  n_samples        = ncol(counts),
  n_DEG_total      = nrow(deg_sig),
  n_DEG_up         = sum(deg_sig$direction=="Up"),
  n_DEG_down       = sum(deg_sig$direction=="Down"),
  n_WGCNA_modules  = length(unique(module_colors))-1,
  sig_WGCNA_modules= paste(sig_modules, collapse=","),
  n_hub_genes      = nrow(hub_genes_all)
)

write.csv(as.data.frame(summary_list),
          file.path(OUTDIR, "mrna_analysis_summary.csv"), row.names=FALSE)

cat("\n", "=" , rep("=",58), "\n", sep="")
cat("✅ mRNA Analysis 완료!\n")
cat(sprintf("  DEG: %d개, WGCNA Hub: %d개\n",
            nrow(deg_sig), nrow(hub_genes_all)))
cat("=" , rep("=",58), "\n\n", sep="")
