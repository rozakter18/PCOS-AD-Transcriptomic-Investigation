library(limma)

OUT_DIR <- "output"

expr_file <- file.path(OUT_DIR, "GSE33000_AD_Control_log2ratio_gene_level_highestvar.csv")
meta_file <- file.path(OUT_DIR, "GSE33000_limma_design_metadata.csv")

expr <- read.csv(expr_file, check.names = FALSE)
meta <- read.csv(meta_file, check.names = FALSE)

cat("Expression file shape:", dim(expr), "\n")
cat("Metadata file shape:", dim(meta), "\n")

required_meta_cols <- c("Sample_ID", "Group", "Age_numeric", "Gender")
missing_meta <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta) > 0) {
  stop(paste("Missing metadata columns:", paste(missing_meta, collapse = ", ")))
}

sample_ids <- meta$Sample_ID
missing_samples <- setdiff(sample_ids, colnames(expr))

if (length(missing_samples) > 0) {
  stop(paste("Samples missing in expression matrix:", paste(missing_samples, collapse = ", ")))
}

expr <- expr[!is.na(expr$Gene_symbol) & expr$Gene_symbol != "", ]
expr <- expr[!duplicated(expr$Gene_symbol), ]

expr_mat <- as.matrix(expr[, sample_ids])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr$Gene_symbol

meta <- meta[match(colnames(expr_mat), meta$Sample_ID), ]

if (!all(meta$Sample_ID == colnames(expr_mat))) {
  stop("Metadata and expression columns are not aligned.")
}

meta$Group <- factor(meta$Group, levels = c("Control", "AD"))
meta$Gender <- factor(meta$Gender)

cat("\nGroup counts:\n")
print(table(meta$Group))

cat("\nGender counts:\n")
print(table(meta$Group, meta$Gender))

cat("\nAge summary:\n")
print(tapply(meta$Age_numeric, meta$Group, summary))

design <- model.matrix(~ Group + Age_numeric + Gender, data = meta)

cat("\nDesign matrix columns:\n")
print(colnames(design))

fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

res <- topTable(
  fit,
  coef = "GroupAD",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

res_out <- data.frame(
  Gene_symbol = rownames(res),
  res,
  row.names = NULL
)

out_deg <- file.path(OUT_DIR, "GSE33000_AD_vs_Control_limma_age_gender_adjusted_DEG.csv")
write.csv(res_out, out_deg, row.names = FALSE)

sig_fdr05 <- res_out[res_out$adj.P.Val < 0.05, ]
out_sig <- file.path(OUT_DIR, "GSE33000_AD_vs_Control_significant_FDR05.csv")
write.csv(sig_fdr05, out_sig, row.names = FALSE)

sig_fdr05_logfc <- res_out[res_out$adj.P.Val < 0.05 & abs(res_out$logFC) >= 0.25, ]
out_sig_logfc <- file.path(OUT_DIR, "GSE33000_AD_vs_Control_significant_FDR05_logFC025.csv")
write.csv(sig_fdr05_logfc, out_sig_logfc, row.names = FALSE)

ranked <- res_out[!is.na(res_out$t), c("Gene_symbol", "t")]
ranked <- ranked[order(ranked$t, decreasing = TRUE), ]

out_rnk <- file.path(OUT_DIR, "GSE33000_AD_vs_Control_ranked_genes_tstat.rnk")
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
cat(out_sig_logfc, "\n")
cat(out_rnk, "\n")

cat("\nFinal DEG summary:\n")
cat("Total genes tested:", nrow(res_out), "\n")
cat("Significant genes FDR < 0.05:", nrow(sig_fdr05), "\n")
cat("Significant genes FDR < 0.05 and abs(logFC) >= 0.25:", nrow(sig_fdr05_logfc), "\n")

cat("\nTop 10 genes:\n")
print(head(res_out, 10))