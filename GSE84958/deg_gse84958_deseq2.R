library(DESeq2)

OUT_DIR <- "output"

counts_file <- file.path(OUT_DIR, "GSE84958_subject_level_counts_filtered.csv")
meta_file   <- file.path(OUT_DIR, "GSE84958_subject_metadata_final.csv")

counts <- read.csv(counts_file, check.names = FALSE, row.names = 1)
meta   <- read.csv(meta_file, check.names = FALSE)

cat("Counts matrix shape:", dim(counts), "\n")
cat("Metadata shape:", dim(meta), "\n")

# Align metadata to count matrix column order
meta <- meta[match(colnames(counts), meta$Subject_ID), ]
if (!all(meta$Subject_ID == colnames(counts))) {
  stop("Metadata and count matrix columns are not aligned.")
}

# Counts must be non-negative integers for DESeq2
counts_mat <- as.matrix(counts)
mode(counts_mat) <- "integer"

meta$Diagnosis <- factor(meta$Diagnosis, levels = c("Normal", "PCOS"))

cat("\nGroup counts:\n")
print(table(meta$Diagnosis))

# No usable numeric covariates available (Age is a bucketed range,
# Gender is constant = Female for all subjects), so the design is
# reduced to diagnosis only. This is a limitation vs. the AD model,
# which adjusted for Age_numeric + Gender.
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta,
  design    = ~ Diagnosis
)

# Pre-filter very low count genes (keep genes with at least 10 reads
# total across samples) to reduce multiple-testing burden and improve
# dispersion estimation stability given the small sample size
keep <- rowSums(counts_mat) >= 10
dds <- dds[keep, ]
cat("\nGenes retained after low-count filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

# Raw MLE log2FC estimates are unstable for low-count genes at this sample
# size (n=7 vs 8) and can blow up to biologically implausible values (e.g.
# |log2FC| > 20) when a gene is near-zero in one group. We use apeglm
# shrinkage to pull noise-driven estimates toward zero while preserving
# well-supported effects. This does not change p-values, only effect sizes.
resultsNames(dds)
coef_name <- "Diagnosis_PCOS_vs_Normal"

if (!requireNamespace("apeglm", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("apeglm", update = FALSE, ask = FALSE)
}

res_raw <- results(dds, contrast = c("Diagnosis", "PCOS", "Normal"), alpha = 0.05)
cat("\nRaw (unshrunk) log2FC range:",
    round(min(res_raw$log2FoldChange, na.rm = TRUE), 2), "to",
    round(max(res_raw$log2FoldChange, na.rm = TRUE), 2), "\n")

res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
cat("Shrunk log2FC range:",
    round(min(res$log2FoldChange, na.rm = TRUE), 2), "to",
    round(max(res$log2FoldChange, na.rm = TRUE), 2), "\n")

# Keep raw p-value/padj from the unshrunk results (lfcShrink with apeglm
# does not recompute p-values by default in a way we want to rely on here);
# we combine shrunk effect sizes with the original Wald test statistics.
res_df <- as.data.frame(res)
res_df$Gene_symbol <- rownames(res_df)
res_df$stat <- res_raw$stat[match(res_df$Gene_symbol, rownames(res_raw))]
res_df <- res_df[order(res_df$pvalue), ]
res_df <- res_df[, c("Gene_symbol", setdiff(colnames(res_df), "Gene_symbol"))]

# Additional filter: require a minimum baseMean, since very low count genes
# remain unreliable even after shrinkage. This is reported as a sensitivity
# filter, not silently applied to the "full" DEG table.
low_expr_flagged <- sum(res_df$baseMean < 10, na.rm = TRUE)
cat("\nGenes with baseMean < 10 (flagged as low-confidence, not removed):",
    low_expr_flagged, "\n")

out_deg <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_deseq2_DEG.csv")
write.csv(res_df, out_deg, row.names = FALSE)

sig_fdr05 <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
out_sig <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_significant_FDR05.csv")
write.csv(sig_fdr05, out_sig, row.names = FALSE)

sig_fdr05_lfc <- sig_fdr05[abs(sig_fdr05$log2FoldChange) >= 0.25, ]
out_sig_lfc <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_significant_FDR05_logFC025.csv")
write.csv(sig_fdr05_lfc, out_sig_lfc, row.names = FALSE)

# Ranked gene list (by test statistic) for downstream GSEA / overlap analysis,
# matching the format used for the AD .rnk file
ranked <- res_df[!is.na(res_df$stat), c("Gene_symbol", "stat")]
ranked <- ranked[order(ranked$stat, decreasing = TRUE), ]
out_rnk <- file.path(OUT_DIR, "GSE84958_PCOS_vs_Control_ranked_genes_stat.rnk")
write.table(
  ranked,
  out_rnk,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("\nSaved files:\n")
cat(out_deg, "\n")
cat(out_sig, "\n")
cat(out_sig_lfc, "\n")
cat(out_rnk, "\n")

cat("\nFinal DEG summary:\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant genes FDR < 0.05:", nrow(sig_fdr05), "\n")
cat("Significant genes FDR < 0.05 and abs(log2FC) >= 0.25:", nrow(sig_fdr05_lfc), "\n")

cat("\nTop 10 genes:\n")
print(head(res_df, 10))