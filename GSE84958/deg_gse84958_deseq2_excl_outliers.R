library(DESeq2)
library(ggplot2)

OUT_DIR <- "output"
EXCLUDE <- c("MO31", "MO37")  # flagged as technical outliers via PCA (see GSE84958_PCA_by_diagnosis.png)

counts_file <- file.path(OUT_DIR, "GSE84958_subject_level_counts_filtered.csv")
meta_file   <- file.path(OUT_DIR, "GSE84958_subject_metadata_final.csv")
libsize_file <- file.path(OUT_DIR, "GSE84958_library_sizes.csv")

counts <- read.csv(counts_file, check.names = FALSE, row.names = 1)
meta   <- read.csv(meta_file, check.names = FALSE)
libsizes <- read.csv(libsize_file, check.names = FALSE)
colnames(libsizes)[1] <- "Subject_ID"

keep_samples <- setdiff(colnames(counts), EXCLUDE)
counts <- counts[, keep_samples]

meta <- meta[match(colnames(counts), meta$Subject_ID), ]
if (!all(meta$Subject_ID == colnames(counts))) {
  stop("Metadata and count matrix columns are not aligned.")
}
meta <- merge(meta, libsizes, by = "Subject_ID", sort = FALSE)
meta <- meta[match(colnames(counts), meta$Subject_ID), ]

cat("Excluded samples:", paste(EXCLUDE, collapse = ", "), "\n")
cat("Remaining samples:", ncol(counts), "\n")
cat("\nGroup counts after exclusion:\n")
print(table(meta$Diagnosis))

counts_mat <- as.matrix(counts)
mode(counts_mat) <- "integer"
meta$Diagnosis <- factor(meta$Diagnosis, levels = c("Normal", "PCOS"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta,
  design    = ~ Diagnosis
)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes retained after low-count filtering:", nrow(dds), "\n")

# --- Rerun PCA on the reduced 13-sample set ---
vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)))
pct_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- data.frame(
  Subject_ID = colnames(vsd),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Diagnosis = meta$Diagnosis,
  Library_size_millions = meta$Total_counts / 1e6
)

cat("\n=== PCA table (13 samples, MO31/MO37 excluded) ===\n")
print(pca_df[order(pca_df$Diagnosis, pca_df$PC1), ], row.names = FALSE)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis, label = Subject_ID)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title = "GSE84958 PCA (13 samples, outliers MO31/MO37 excluded)",
    x = paste0("PC1 (", pct_var[1], "% variance)"),
    y = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_minimal()

out_plot <- file.path(OUT_DIR, "GSE84958_PCA_excl_outliers.png")
ggsave(out_plot, p1, width = 7, height = 5.5, dpi = 150)
cat("\nSaved plot:", out_plot, "\n")

# --- Rerun DEG analysis on the reduced set ---
dds <- DESeq(dds)
resultsNames(dds)
coef_name <- "Diagnosis_PCOS_vs_Normal"

res_raw <- results(dds, contrast = c("Diagnosis", "PCOS", "Normal"), alpha = 0.05)
cat("\nRaw (unshrunk) log2FC range:",
    round(min(res_raw$log2FoldChange, na.rm = TRUE), 2), "to",
    round(max(res_raw$log2FoldChange, na.rm = TRUE), 2), "\n")

res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
cat("Shrunk log2FC range:",
    round(min(res$log2FoldChange, na.rm = TRUE), 2), "to",
    round(max(res$log2FoldChange, na.rm = TRUE), 2), "\n")

res_df <- as.data.frame(res)
res_df$Gene_symbol <- rownames(res_df)
res_df$stat <- res_raw$stat[match(res_df$Gene_symbol, rownames(res_raw))]
res_df <- res_df[order(res_df$pvalue), ]
res_df <- res_df[, c("Gene_symbol", setdiff(colnames(res_df), "Gene_symbol"))]

out_deg <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_deseq2_DEG_excl_outliers.csv")
write.csv(res_df, out_deg, row.names = FALSE)

sig_fdr05 <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
out_sig <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_significant_FDR05_excl_outliers.csv")
write.csv(sig_fdr05, out_sig, row.names = FALSE)

ranked <- res_df[!is.na(res_df$stat), c("Gene_symbol", "stat")]
ranked <- ranked[order(ranked$stat, decreasing = TRUE), ]
out_rnk <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_ranked_genes_stat_excl_outliers.rnk")
write.table(ranked, out_rnk, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("\nSaved files:\n")
cat(out_deg, "\n")
cat(out_sig, "\n")
cat(out_rnk, "\n")

cat("\nFinal DEG summary (outliers excluded):\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant genes FDR < 0.05:", nrow(sig_fdr05), "\n")

cat("\nTop 10 genes:\n")
print(head(res_df, 10))