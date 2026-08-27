# ============================================================
# TCGA-COAD: Differential Expression Analysis
# Tumour vs Normal
# ============================================================

# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)  # "assay" function is stored in this library
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

# ------------------------------------------------------------
# 2. Query & download TCGA-COAD RNA-seq data (GDC)
# ------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query)
data <- GDCprepare(query)

# Confirm available count assays and select the raw unstranded counts
assayNames(data)
counts <- assay(data, "unstranded")

# ------------------------------------------------------------
# 3. Build & clean sample metadata
# ------------------------------------------------------------
metadata <- as.data.frame(colData(data))

metadata$condition <- ifelse(
  metadata$sample_type == "Primary Tumor",
  "Tumor",
  ifelse(
    metadata$sample_type == "Solid Tissue Normal",
    "Normal",
    NA
  )
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Normal", "Tumor")
)

table(metadata$condition, useNA = "ifany")

# Drop any samples with unresolved condition labels
keep_samples <- !is.na(metadata$condition)
metadata <- metadata[keep_samples, ]
counts <- counts[, keep_samples]

# Align metadata row order to count matrix column order, then verify
metadata <- metadata[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(metadata)))

table(metadata$condition)

# ------------------------------------------------------------
# 4. Filter low-expression genes
# ------------------------------------------------------------
# Keep genes with count >= 10 in at least 3 samples
keep_genes <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep_genes, ]
dim(counts_filtered)

# ------------------------------------------------------------
# 5. Differential expression analysis (DESeq2)
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

res_df <- as.data.frame(res)

# Remove genes with missing adjusted p-values
res_df <- res_df[!is.na(res_df$padj), ]

# ------------------------------------------------------------
# 6. Annotate with gene symbols
# ------------------------------------------------------------
res_df$ensembl_id <- rownames(res_df)
res_df$ensembl_id <- sub("\\..*", "", res_df$ensembl_id)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
res_df$gene_symbol <- gene_symbols

# ------------------------------------------------------------
# 7. Classify significant DEGs
# ------------------------------------------------------------
sig_genes <- res_df[
  res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1,
]
sig_genes <- sig_genes[order(sig_genes$padj), ]

upregulated <- res_df[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 1,
]
upregulated <- upregulated[order(upregulated$padj), ]

downregulated <- res_df[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < -1,
]
downregulated <- downregulated[order(downregulated$padj), ]

cat("Total significant DEGs:", nrow(sig_genes), "\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# Top 20 up/downregulated genes for a quick look
cat("\nTop 20 Upregulated Genes:\n")
print(
  upregulated[
    1:min(20, nrow(upregulated)),
    c("gene_symbol", "ensembl_id", "baseMean", "log2FoldChange", "pvalue", "padj")
  ]
)

cat("\nTop 20 Downregulated Genes:\n")
print(
  downregulated[
    1:min(20, nrow(downregulated)),
    c("gene_symbol", "ensembl_id", "baseMean", "log2FoldChange", "pvalue", "padj")
  ]
)

# ------------------------------------------------------------
# 8. Visualisation - PCA
# ------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")
# ggsave("figures/pca_plot.png", width = 8, height = 6)

# ------------------------------------------------------------
# 9. Visualisation - Volcano plot
# ------------------------------------------------------------
res_df$category <- "Not significant"
res_df$category[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 1
] <- "Upregulated"
res_df$category[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < -1
] <- "Downregulated"

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = category
  )
) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "TCGA-COAD: Tumour vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Category"
  )

print(volcano_plot)
# ggsave("figures/volcano_plot.png", plot = volcano_plot, width = 10, height = 7)

# ------------------------------------------------------------
# 10. Visualisation - Heatmap of top 50 DEGs
# ------------------------------------------------------------
top50 <- head(sig_genes[order(sig_genes$padj), ], 50)
top50_ids <- top50$ensembl_id

vsd_mat <- assay(vsd)
vsd_ids <- sub("\\..*", "", rownames(vsd_mat))

gene_index <- match(top50_ids, vsd_ids)
gene_index <- gene_index[!is.na(gene_index)]
cat("Number of top DEGs found:", length(gene_index), "\n")

heatmap_mat <- vsd_mat[gene_index, , drop = FALSE]

gene_symbols_hm <- top50$gene_symbol[match(top50_ids, vsd_ids)]
gene_symbols_hm <- gene_symbols_hm[!is.na(match(top50_ids, vsd_ids))]
gene_symbols_hm[is.na(gene_symbols_hm) | gene_symbols_hm == ""] <-
  top50_ids[is.na(gene_symbols_hm) | gene_symbols_hm == ""]
gene_symbols_hm <- make.unique(gene_symbols_hm)

rownames(heatmap_mat) <- gene_symbols_hm

annotation_col <- data.frame(Condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(vsd)

pheatmap(
  heatmap_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  main = "Top 50 DEGs - TCGA-COAD Tumour vs Normal"
  # filename = "figures/heatmap_top50.png"
)

# ------------------------------------------------------------
# 11. Functional enrichment - GO Biological Process
# ------------------------------------------------------------
up_ids <- unique(upregulated$ensembl_id[!is.na(upregulated$ensembl_id)])
down_ids <- unique(downregulated$ensembl_id[!is.na(downregulated$ensembl_id)])
background_ids <- unique(res_df$ensembl_id[!is.na(res_df$ensembl_id)])

up_entrez <- bitr(up_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_entrez <- bitr(down_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
background_entrez <- bitr(background_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

up_entrez <- unique(up_entrez$ENTREZID)
down_entrez <- unique(down_entrez$ENTREZID)
background_entrez <- unique(background_entrez$ENTREZID)

cat("Upregulated genes (Entrez):", length(up_entrez), "\n")
cat("Downregulated genes (Entrez):", length(down_entrez), "\n")
cat("Background genes (Entrez):", length(background_entrez), "\n")

GO_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

GO_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat("\n===== UPREGULATED GO TERMS =====\n")
head(as.data.frame(GO_up), 20)

cat("\n===== DOWNREGULATED GO TERMS =====\n")
head(as.data.frame(GO_down), 20)

dotplot(GO_up, showCategory = 20) +
  ggtitle("GO Biological Process - Upregulated Genes")
# ggsave("figures/go_upregulated_dotplot.png", width = 10, height = 8)

dotplot(GO_down, showCategory = 20) +
  ggtitle("GO Biological Process - Downregulated Genes")
# ggsave("figures/go_downregulated_dotplot.png", width = 10, height = 8)

# ------------------------------------------------------------
# 12. Export results
# ------------------------------------------------------------
write.csv(sig_genes, "results/TCGA_COAD_significant_DEGs.csv", row.names = FALSE)
write.csv(upregulated, "results/TCGA_COAD_upregulated_genes.csv", row.names = FALSE)
write.csv(downregulated, "results/TCGA_COAD_downregulated_genes.csv", row.names = FALSE)
write.csv(res_df, "results/TCGA_COAD_all_DESeq2_results.csv", row.names = FALSE)

write.csv(as.data.frame(GO_up), "results/TCGA_COAD_GO_Upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(GO_down), "results/TCGA_COAD_GO_Downregulated.csv", row.names = FALSE)

GO_up_top20 <- head(as.data.frame(GO_up), 20)
GO_down_top20 <- head(as.data.frame(GO_down), 20)
write.csv(GO_up_top20, "results/TCGA_COAD_GO_Upregulated_Top20.csv", row.names = FALSE)
write.csv(GO_down_top20, "results/TCGA_COAD_GO_Downregulated_Top20.csv", row.names = FALSE)

cat("\nAnalysis complete!\n")
