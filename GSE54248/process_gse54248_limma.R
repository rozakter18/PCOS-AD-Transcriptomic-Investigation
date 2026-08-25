library(limma)

OUT_DIR <- "output"
dir.create(OUT_DIR, showWarnings = FALSE)

profile_file <- "GSE54248_Raw_Data_Sample_Probe_Profile.txt/GSE54248_Raw_Data_Sample_Probe_Profile.txt"

# Sample code -> group mapping, based on series matrix order:
# Healthy Individual 1-4 = GSM1310933-936, PCOS Patient 1-4 = GSM1310937-940.
# GenomeStudio sample codes in this file are P1-P4 and N1-N4.
# We assume N = Normal/Healthy, P = PCOS (matches counts: 4 healthy, 4 PCOS).
# VERIFY this assumption before trusting downstream results.
sample_codes <- c("N1", "N2", "N3", "N4", "P1", "P2", "P3", "P4")
sample_group <- c("Control", "Control", "Control", "Control",
                   "PCOS", "PCOS", "PCOS", "PCOS")
sample_gsm   <- c("GSM1310933", "GSM1310934", "GSM1310935", "GSM1310936",
                   "GSM1310937", "GSM1310938", "GSM1310939", "GSM1310940")

cat("Assumed sample code -> group -> GSM mapping:\n")
print(data.frame(Code = sample_codes, Group = sample_group, GSM = sample_gsm))
cat("\n>>> Double-check this against the series matrix Sample_title/geo_accession",
    "\n>>> rows before trusting results downstream.\n\n")

# read.ilmn expects column prefixes like "AVG_Signal" and "Detection" with
# per-sample suffixes matching sample_codes (i.e., AVG_Signal-N1, Detection-N1, ...)
x <- read.ilmn(
  files = profile_file,
  probeid = "ProbeID",
  annotation = "TargetID",
  expr = "AVG_Signal",
  other.columns = "Detection"
)

cat("Raw data dimensions (probes x samples):", dim(x), "\n")
cat("Sample columns as read:\n")
print(colnames(x$E))

# Rename columns to Group_Code (e.g., Control_N1) for clarity downstream
match_idx <- match(sample_codes, sub(".*-", "", colnames(x$E)))
if (any(is.na(match_idx))) {
  cat("\nWARNING: could not automatically match all sample codes to columns.\n")
  cat("Columns found:", paste(colnames(x$E), collapse = ", "), "\n")
  stop("Fix sample code matching before proceeding.")
}
colnames(x$E) <- paste0(sample_group, "_", sample_codes)
colnames(x$other$Detection) <- colnames(x$E)

group <- factor(sample_group, levels = c("Control", "PCOS"))

# NOTE: In this file's export, the "Detection" column is a detection
# CONFIDENCE score (values near 1 = confidently detected), not a raw
# p-value (which would be near 0 for detected probes). We verified this
# empirically: well-expressed probes show Detection ~0.9-0.98, while
# low-signal/noise probes show lower, more scattered values. So we filter
# for HIGH confidence, not low p-value.
cat("\nDetection column summary (sanity check - should NOT look like p-values,\n")
cat("i.e. should NOT be concentrated near 0):\n")
print(summary(as.vector(x$other$Detection)))

detected <- rowSums(x$other$Detection > 0.95) >= 4  # confidently detected in at least half the samples
cat("\nProbes detected (confidence > 0.95) in >= 4/8 samples:", sum(detected), "of", nrow(x), "\n")

x_filtered <- x[detected, ]

# Background correction (normexp) + quantile normalization + log2 transform,
# the standard combined step for Illumina BeadArray data (neqc)
y <- neqc(x_filtered)

cat("\nDimensions after neqc normalization:", dim(y), "\n")

# Collapse to gene symbol level (multiple probes can map to the same gene;
# take the probe with highest average expression per gene, standard practice)
gene_symbol <- y$genes$TargetID
expr <- y$E
avg_expr <- rowMeans(expr)
keep_idx <- !duplicated(gene_symbol[order(-avg_expr)])
ord <- order(-avg_expr)
expr_ordered <- expr[ord, ]
symbol_ordered <- gene_symbol[ord]
dedup_keep <- !duplicated(symbol_ordered)
expr_gene <- expr_ordered[dedup_keep, ]
rownames(expr_gene) <- symbol_ordered[dedup_keep]

cat("Genes after collapsing probes to highest-expressed per symbol:", nrow(expr_gene), "\n")

# ---- limma differential expression: PCOS vs Control ----
design <- model.matrix(~ group)
fit <- lmFit(expr_gene, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "groupPCOS", number = Inf, adjust.method = "BH", sort.by = "P")
res_out <- data.frame(Gene_symbol = rownames(res), res, row.names = NULL)

out_deg <- file.path(OUT_DIR, "GSE54248_PCOS_vs_Control_limma_DEG.csv")
write.csv(res_out, out_deg, row.names = FALSE)

sig_fdr05 <- res_out[res_out$adj.P.Val < 0.05, ]
out_sig <- file.path(OUT_DIR, "GSE54248_PCOS_vs_Control_significant_FDR05.csv")
write.csv(sig_fdr05, out_sig, row.names = FALSE)

ranked <- res_out[!is.na(res_out$t), c("Gene_symbol", "t")]
ranked <- ranked[order(ranked$t, decreasing = TRUE), ]
out_rnk <- file.path(OUT_DIR, "GSE54248_PCOS_vs_Control_ranked_genes_tstat.rnk")
write.table(ranked, out_rnk, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("\nSaved files:\n")
cat(out_deg, "\n")
cat(out_sig, "\n")
cat(out_rnk, "\n")

cat("\nFinal DEG summary (GSE54248, blood, n=4 vs 4):\n")
cat("Total genes tested:", nrow(res_out), "\n")
cat("Significant genes FDR < 0.05:", nrow(sig_fdr05), "\n")
cat("\nNOTE: n=4 vs 4 is very small. Treat this strictly as a sensitivity/\n")
cat("validation check against the GSE84958 (adipose) findings, not a\n")
cat("standalone discovery analysis.\n")

cat("\nTop 10 genes:\n")
print(head(res_out, 10))