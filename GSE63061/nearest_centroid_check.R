library(sva)
library(limma)

set.seed(42)

# ---- 1. Reload and re-harmonize (same as the ComBat classifier script) ----
ad_expr <- readRDS("output/GSE63061_AD_CTL_expr_genelevel.rds")
ad_meta <- read.csv("output/GSE63061_AD_CTL_metadata.csv", check.names = FALSE)
ad_meta$Status <- factor(ad_meta$Status, levels = c("CTL", "AD"))
ad_expr <- as.matrix(ad_expr)

profile_file <- "../GSE54248/GSE54248_Raw_Data_Sample_Probe_Profile.txt/GSE54248_Raw_Data_Sample_Probe_Profile.txt"
sample_codes <- c("N1", "N2", "N3", "N4", "P1", "P2", "P3", "P4")
sample_group <- c("Control", "Control", "Control", "Control", "PCOS", "PCOS", "PCOS", "PCOS")

x <- read.ilmn(files = profile_file, probeid = "ProbeID", annotation = "TargetID",
                expr = "AVG_Signal", other.columns = "Detection")
colnames(x$E) <- paste0(sample_group, "_", sample_codes)
colnames(x$other$Detection) <- colnames(x$E)
detected <- rowSums(x$other$Detection > 0.95) >= 4
x_filtered <- x[detected, ]
y <- neqc(x_filtered)

pcos_gene_symbol <- y$genes$TargetID
avg_expr <- rowMeans(y$E)
ord <- order(-avg_expr)
pcos_expr_ordered <- y$E[ord, ]
symbol_ordered <- pcos_gene_symbol[ord]
dedup_keep <- !duplicated(symbol_ordered)
pcos_expr <- pcos_expr_ordered[dedup_keep, ]
rownames(pcos_expr) <- symbol_ordered[dedup_keep]
pcos_diagnosis <- factor(sample_group, levels = c("Control", "PCOS"))

common_genes <- intersect(rownames(ad_expr), rownames(pcos_expr))
ad_common <- ad_expr[common_genes, ]
pcos_common <- pcos_expr[common_genes, ]
complete_genes <- intersect(rownames(ad_common)[complete.cases(ad_common)],
                              rownames(pcos_common)[complete.cases(pcos_common)])
ad_common <- ad_common[complete_genes, ]
pcos_common <- pcos_common[complete_genes, ]

combined <- cbind(ad_common, pcos_common)
batch <- factor(c(rep("AD_blood", ncol(ad_common)), rep("PCOS_blood", ncol(pcos_common))))
mod <- model.matrix(~1, data = data.frame(batch))
combined_harmonized <- ComBat(dat = combined, batch = batch, mod = mod)

ad_harmonized   <- combined_harmonized[, colnames(ad_common)]
pcos_harmonized <- combined_harmonized[, colnames(pcos_common)]

labels <- ad_meta$Status[match(colnames(ad_harmonized), ad_meta$Chip_ID)]

cat("Data reloaded and re-harmonized. Genes:", nrow(ad_harmonized),
    "| AD-blood samples:", sum(labels == "AD"), "AD,", sum(labels == "CTL"), "CTL\n")

# ---- 2. Compute AD and CTL centroids from AD-blood data (mean expression profile per group) ----
ad_centroid <- rowMeans(ad_harmonized[, labels == "AD"])
ctl_centroid <- rowMeans(ad_harmonized[, labels == "CTL"])

cat("\nCentroids computed. Genes per centroid:", length(ad_centroid), "\n")

# ---- 3. For each PCOS-blood sample, compute correlation distance to each centroid ----
# Correlation distance = 1 - Pearson correlation; robust to overall scale shifts
# and does not require any model fitting or training.
compute_distances <- function(sample_vec, ad_cent, ctl_cent) {
  dist_ad  <- 1 - cor(sample_vec, ad_cent, method = "pearson")
  dist_ctl <- 1 - cor(sample_vec, ctl_cent, method = "pearson")
  c(dist_to_AD = dist_ad, dist_to_CTL = dist_ctl)
}

results_list <- lapply(seq_len(ncol(pcos_harmonized)), function(i) {
  compute_distances(pcos_harmonized[, i], ad_centroid, ctl_centroid)
})
dist_df <- as.data.frame(do.call(rbind, results_list))
dist_df$Sample <- colnames(pcos_harmonized)
dist_df$Diagnosis <- pcos_diagnosis
dist_df$Closer_to <- ifelse(dist_df$dist_to_AD < dist_df$dist_to_CTL, "AD", "CTL")

# A continuous "AD-likeness" analog: relative distance, so it's comparable in
# spirit to the classifier's probability score (higher = relatively closer to AD)
dist_df$AD_likeness_by_distance <- round(
  dist_df$dist_to_CTL / (dist_df$dist_to_AD + dist_df$dist_to_CTL), 3
)

dist_df <- dist_df[, c("Sample", "Diagnosis", "dist_to_AD", "dist_to_CTL", "Closer_to", "AD_likeness_by_distance")]

cat("\n=== Nearest-centroid results (independent method, no training/CV) ===\n")
print(dist_df)

cat("\n=== Summary: which centroid is each group closer to, on average? ===\n")
print(table(Diagnosis = dist_df$Diagnosis, Closer_to = dist_df$Closer_to))

cat("\nPCOS mean AD-likeness (by distance):",
    round(mean(dist_df$AD_likeness_by_distance[dist_df$Diagnosis == "PCOS"]), 3), "\n")
cat("Control mean AD-likeness (by distance):",
    round(mean(dist_df$AD_likeness_by_distance[dist_df$Diagnosis == "Control"]), 3), "\n")

cat("\n=== Cross-method comparison ===\n")
cat("Compare this table's 'Closer_to' / AD_likeness_by_distance against the\n")
cat("logistic regression classifier's AD_likeness_prob from the previous run.\n")
cat("If both methods agree in direction (which group scores higher), that is\n")
cat("a real, method-independent robustness finding.\n")

dir.create("../results/ml_blood_matched", recursive = TRUE, showWarnings = FALSE)
write.csv(dist_df, "../results/ml_blood_matched/PCOS_blood_nearest_centroid.csv", row.names = FALSE)
cat("\nSaved: results/ml_blood_matched/PCOS_blood_nearest_centroid.csv\n")
