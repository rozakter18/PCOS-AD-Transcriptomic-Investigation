library(pROC)

set.seed(42)

# ---- 1. Load AD-blood (GSE63061) gene-level expression + metadata ----
ad_expr <- readRDS("output/GSE63061_AD_CTL_expr_genelevel.rds")
ad_meta <- read.csv("output/GSE63061_AD_CTL_metadata.csv", check.names = FALSE)
ad_meta$Status <- factor(ad_meta$Status, levels = c("CTL", "AD"))

cat("AD-blood matrix:", dim(ad_expr), "\n")
cat("AD-blood value range:", round(range(ad_expr, na.rm = TRUE), 2), "\n")

# ---- 2. Load PCOS-blood (GSE54248) gene-level expression + metadata ----
# Reuse the limma DEG script's intermediate processing logic - re-derive the
# normalized expression matrix the same way (neqc output) rather than re-run
# limma; we need the *expression matrix*, not just the DEG table this time.
# For speed, re-derive from the raw probe profile using the same pipeline.
library(limma)

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
pcos_expr_all <- y$E
avg_expr <- rowMeans(pcos_expr_all)
ord <- order(-avg_expr)
pcos_expr_ordered <- pcos_expr_all[ord, ]
symbol_ordered <- pcos_gene_symbol[ord]
dedup_keep <- !duplicated(symbol_ordered)
pcos_expr <- pcos_expr_ordered[dedup_keep, ]
rownames(pcos_expr) <- symbol_ordered[dedup_keep]

pcos_diagnosis <- factor(sample_group, levels = c("Control", "PCOS"))
cat("\nPCOS-blood matrix:", dim(pcos_expr), "\n")
cat("PCOS-blood value range:", round(range(pcos_expr, na.rm = TRUE), 2), "\n")

# ---- 3. Check whether the two datasets are on a comparable scale already ----
cat("\n=== Platform-match check: are AD-blood and PCOS-blood value ranges comparable? ===\n")
cat("If yes (similar ranges, similar medians), ComBat may be unnecessary since\n")
cat("both are already log2-scale, quantile-normalized Illumina HT-12 V4.0 data.\n")
cat("AD-blood median:", round(median(as.matrix(ad_expr), na.rm = TRUE), 2), "\n")
cat("PCOS-blood median:", round(median(pcos_expr, na.rm = TRUE), 2), "\n")

# ---- 4. Restrict to common genes ----
common_genes <- intersect(rownames(ad_expr), rownames(pcos_expr))
cat("\nCommon genes across both blood datasets:", length(common_genes), "\n")

ad_expr_common <- as.matrix(ad_expr[common_genes, ])
pcos_expr_common <- pcos_expr[common_genes, ]

# ---- 5. PCA + logistic regression on AD-blood, with proper 5-fold stratified CV ----
n_pcs <- 20
labels <- ad_meta$Status[match(colnames(ad_expr_common), ad_meta$Chip_ID)]

n_folds <- 5
fold_id <- integer(length(labels))
for (lvl in levels(labels)) {
  idx <- which(labels == lvl)
  idx <- sample(idx)
  fold_id[idx] <- rep(1:n_folds, length.out = length(idx))
}

cv_pred_prob <- rep(NA, length(labels))
for (k in 1:n_folds) {
  train_idx <- which(fold_id != k)
  test_idx  <- which(fold_id == k)

  pca_k <- prcomp(t(ad_expr_common[, train_idx]), scale. = TRUE, center = TRUE)
  train_pcs_k <- pca_k$x[, 1:n_pcs]

  test_centered_k <- scale(t(ad_expr_common[, test_idx]), center = pca_k$center, scale = pca_k$scale)
  test_pcs_k <- test_centered_k %*% pca_k$rotation[, 1:n_pcs]

  train_df_k <- data.frame(Status = labels[train_idx], train_pcs_k)
  model_k <- glm(Status ~ ., data = train_df_k, family = binomial)

  test_df_k <- as.data.frame(test_pcs_k)
  cv_pred_prob[test_idx] <- predict(model_k, newdata = test_df_k, type = "response")
}

cv_auc <- roc(response = labels, predictor = cv_pred_prob, levels = c("CTL", "AD"), quiet = TRUE)
cv_pred_class <- factor(ifelse(cv_pred_prob >= 0.5, "AD", "CTL"), levels = c("CTL", "AD"))
cv_accuracy <- mean(cv_pred_class == labels)

cat("\n=== AD-BLOOD classifier: honest cross-validated performance ===\n")
cat("CV AUC:", round(as.numeric(cv_auc$auc), 3), "\n")
cat("CV Accuracy:", round(cv_accuracy, 3), "\n")
print(table(True = labels, Predicted = cv_pred_class))

# ---- 6. Train final model on ALL AD-blood data, project PCOS-blood samples ----
pca_full <- prcomp(t(ad_expr_common), scale. = TRUE, center = TRUE)
ad_pcs_full <- pca_full$x[, 1:n_pcs]
train_df_full <- data.frame(Status = labels, ad_pcs_full)
model_full <- glm(Status ~ ., data = train_df_full, family = binomial)

pcos_centered <- scale(t(pcos_expr_common), center = pca_full$center, scale = pca_full$scale)
pcos_pcs <- pcos_centered %*% pca_full$rotation[, 1:n_pcs]
pcos_pred_prob <- predict(model_full, newdata = as.data.frame(pcos_pcs), type = "response")

cat("\n=== PCOS-blood AD-likeness scores (TISSUE-MATCHED, no ComBat) ===\n")
pcos_results <- data.frame(Sample = colnames(pcos_expr_common),
                             Diagnosis = pcos_diagnosis,
                             AD_likeness_prob = round(pcos_pred_prob, 3))
print(pcos_results)

cat("\nPCOS mean:", round(mean(pcos_pred_prob[pcos_diagnosis == "PCOS"]), 3), "\n")
cat("Control mean:", round(mean(pcos_pred_prob[pcos_diagnosis == "Control"]), 3), "\n")
cat("SD across all PCOS-blood samples:", round(sd(pcos_pred_prob), 3), "\n")

dir.create("../results/ml_blood_matched", recursive = TRUE, showWarnings = FALSE)
write.csv(pcos_results, "../results/ml_blood_matched/PCOS_blood_AD_likeness_tissue_matched.csv", row.names = FALSE)
cat("\nSaved: results/ml_blood_matched/PCOS_blood_AD_likeness_tissue_matched.csv\n")
