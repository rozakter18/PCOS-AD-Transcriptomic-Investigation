library(sva)
library(dplyr)

# ---- Run from PCOS_AD root ----

# ---- 1. Load AD data (microarray, log2-ratio gene-level matrix) ----
ad_expr_file <- "GSE33000/output/GSE33000_AD_Control_log2ratio_gene_level_highestvar.csv"
ad_meta_file <- "GSE33000/output/GSE33000_limma_design_metadata.csv"

ad_expr <- read.csv(ad_expr_file, check.names = FALSE)
ad_meta <- read.csv(ad_meta_file, check.names = FALSE)

ad_expr <- ad_expr[!is.na(ad_expr$Gene_symbol) & ad_expr$Gene_symbol != "", ]
ad_expr <- ad_expr[!duplicated(ad_expr$Gene_symbol), ]
ad_mat <- as.matrix(ad_expr[, ad_meta$Sample_ID])
mode(ad_mat) <- "numeric"
rownames(ad_mat) <- ad_expr$Gene_symbol

cat("AD matrix:", dim(ad_mat), "\n")
cat("AD matrix missing values:", sum(is.na(ad_mat)), "of", length(ad_mat),
    "(", round(100 * sum(is.na(ad_mat)) / length(ad_mat), 2), "%)\n")
na_per_gene <- rowSums(is.na(ad_mat))
cat("Genes with at least 1 missing value:", sum(na_per_gene > 0), "of", nrow(ad_mat), "\n")

# ---- 2. Load PCOS data (RNA-seq raw counts), exclude the 2 flagged outliers ----
pcos_counts_file <- "GSE84958/output/GSE84958_subject_level_counts_filtered.csv"
pcos_meta_file   <- "GSE84958/output/GSE84958_subject_metadata_final.csv"
EXCLUDE <- c("MO31", "MO37")

pcos_counts <- read.csv(pcos_counts_file, check.names = FALSE, row.names = 1)
pcos_meta   <- read.csv(pcos_meta_file, check.names = FALSE)
pcos_counts <- pcos_counts[, setdiff(colnames(pcos_counts), EXCLUDE)]
pcos_meta   <- pcos_meta[pcos_meta$Subject_ID %in% colnames(pcos_counts), ]
pcos_meta   <- pcos_meta[match(colnames(pcos_counts), pcos_meta$Subject_ID), ]

# Convert raw counts -> log2 CPM, comparable scale to microarray log2 intensities
lib_sizes <- colSums(pcos_counts)
cpm <- sweep(pcos_counts, 2, lib_sizes, "/") * 1e6
pcos_mat <- log2(cpm + 1)
pcos_mat <- as.matrix(pcos_mat)

cat("PCOS matrix (after excluding outliers):", dim(pcos_mat), "\n")

# ---- 3. Restrict to common genes across both platforms ----
common_genes <- intersect(rownames(ad_mat), rownames(pcos_mat))
cat("Common genes across both platforms:", length(common_genes), "\n")

if (length(common_genes) < 500) {
  stop("Too few common genes for reliable harmonization/classification. Stopping.")
}

ad_mat_common   <- ad_mat[common_genes, ]
pcos_mat_common <- pcos_mat[common_genes, ]

# The AD log2-ratio matrix has missing values for some genes/samples (probe
# filtering artifacts). ComBat propagates these into fully-NaN rows, which
# then breaks PCA. Remove any gene with ANY missing value in either matrix
# before proceeding, rather than letting NAs silently corrupt downstream steps.
complete_ad   <- rownames(ad_mat_common)[complete.cases(ad_mat_common)]
complete_pcos <- rownames(pcos_mat_common)[complete.cases(pcos_mat_common)]
complete_genes <- intersect(complete_ad, complete_pcos)

cat("Genes with complete data (no NA) in both matrices:", length(complete_genes),
    "of", length(common_genes), "common genes\n")

if (length(complete_genes) < 500) {
  stop("Too few complete-case genes for reliable harmonization/classification. Stopping.")
}

ad_mat_common   <- ad_mat_common[complete_genes, ]
pcos_mat_common <- pcos_mat_common[complete_genes, ]

# ---- 4. Combine and run ComBat with dataset as batch ----
combined <- cbind(ad_mat_common, pcos_mat_common)
batch <- factor(c(rep("AD_microarray", ncol(ad_mat_common)),
                   rep("PCOS_rnaseq", ncol(pcos_mat_common))))

cat("\nCombined matrix before harmonization:", dim(combined), "\n")
cat("Batch counts:\n")
print(table(batch))

# No biological covariate is shared across the two batches (AD status only
# exists for the AD samples, PCOS status only for the PCOS samples), so we
# use an intercept-only model -- ComBat will remove the platform/dataset
# effect without trying to preserve a cross-dataset biological signal that
# doesn't actually exist. This is a known limitation of this approach and
# is stated explicitly rather than hidden.
mod <- model.matrix(~1, data = data.frame(batch))
combined_harmonized <- ComBat(dat = combined, batch = batch, mod = mod)

ad_harmonized   <- combined_harmonized[, colnames(ad_mat_common)]
pcos_harmonized <- combined_harmonized[, colnames(pcos_mat_common)]

cat("\nHarmonization complete.\n")
cat("AD harmonized matrix range:", round(range(ad_harmonized), 2), "\n")
cat("PCOS harmonized matrix range:", round(range(pcos_harmonized), 2), "\n")

# ---- 5. PCA on harmonized AD data, then project PCOS onto same PCs ----
ad_meta$Group <- factor(ad_meta$Group, levels = c("Control", "AD"))

pca <- prcomp(t(ad_harmonized), scale. = TRUE, center = TRUE)
n_pcs <- 20
ad_pcs <- pca$x[, 1:n_pcs]

# Project PCOS samples into the same PC space using AD's rotation/center/scale
pcos_centered <- scale(t(pcos_harmonized), center = pca$center, scale = pca$scale)
pcos_pcs <- pcos_centered %*% pca$rotation[, 1:n_pcs]

cat("\nVariance explained by first", n_pcs, "PCs:",
    round(sum((pca$sdev[1:n_pcs])^2) / sum(pca$sdev^2) * 100, 1), "%\n")

# ---- 6. Honest generalization estimate via stratified 5-fold cross-validation ----
# The self-prediction check above uses training data - this overstates real
# performance. We now hold out folds properly to get an unbiased estimate of
# how well this classifier actually generalizes to unseen AD/Control samples,
# BEFORE trusting anything about the PCOS projection.
library(pROC)

set.seed(42)
n_folds <- 5
labels <- ad_meta$Group  # factor: Control, AD

# Stratified fold assignment (preserve AD:Control ratio in each fold)
fold_id <- integer(length(labels))
for (lvl in levels(labels)) {
  idx <- which(labels == lvl)
  idx <- sample(idx)  # shuffle
  fold_id[idx] <- rep(1:n_folds, length.out = length(idx))
}

cv_pred_prob <- rep(NA, length(labels))

for (k in 1:n_folds) {
  train_idx <- which(fold_id != k)
  test_idx  <- which(fold_id == k)

  # Refit PCA on the training fold ONLY (avoids leakage - test fold must
  # never influence the PC space or the model)
  pca_k <- prcomp(t(ad_harmonized[, train_idx]), scale. = TRUE, center = TRUE)
  train_pcs_k <- pca_k$x[, 1:n_pcs]

  test_centered_k <- scale(t(ad_harmonized[, test_idx]),
                             center = pca_k$center, scale = pca_k$scale)
  test_pcs_k <- test_centered_k %*% pca_k$rotation[, 1:n_pcs]

  train_df_k <- data.frame(Group = labels[train_idx], train_pcs_k)
  model_k <- glm(Group ~ ., data = train_df_k, family = binomial)

  test_df_k <- as.data.frame(test_pcs_k)
  cv_pred_prob[test_idx] <- predict(model_k, newdata = test_df_k, type = "response")
}

cv_auc <- roc(response = labels, predictor = cv_pred_prob, levels = c("Control", "AD"), quiet = TRUE)
cv_pred_class <- factor(ifelse(cv_pred_prob >= 0.5, "AD", "Control"), levels = c("Control", "AD"))
cv_accuracy <- mean(cv_pred_class == labels)
cv_confusion <- table(True = labels, Predicted = cv_pred_class)

cat("\n=== HONEST cross-validated generalization performance (5-fold, held-out) ===\n")
cat("This is the trustworthy performance estimate -- unlike the self-prediction\n")
cat("table above, no test sample's own data was used to fit its prediction.\n\n")
cat("CV AUC:", round(as.numeric(cv_auc$auc), 3), "\n")
cat("CV Accuracy (threshold 0.5):", round(cv_accuracy, 3), "\n")
cat("\nConfusion matrix (rows=true, cols=predicted):\n")
print(cv_confusion)

cv_out <- data.frame(Sample = ad_meta$Sample_ID, True_Group = labels,
                       CV_predicted_prob_AD = round(cv_pred_prob, 3),
                       CV_predicted_class = cv_pred_class)
write.csv(cv_out, file.path("results/ml_combat_attempt", "AD_classifier_CV_performance.csv"), row.names = FALSE)
cat("\nSaved:", "results/ml_combat_attempt/AD_classifier_CV_performance.csv\n")

# ---- 7. Train logistic regression on AD (Control vs AD), predict on PCOS ----
train_df <- data.frame(Group = ad_meta$Group, ad_pcs)
model <- glm(Group ~ ., data = train_df, family = binomial)

ad_pred_prob <- predict(model, type = "response")
pcos_pred_df <- as.data.frame(pcos_pcs)
pcos_pred_prob <- predict(model, newdata = pcos_pred_df, type = "response")

cat("\n=== AD model TRAINING-DATA prediction (NOT a generalization estimate --\n")
cat("this model has already seen every one of these samples during fitting,\n")
cat("so this table is optimistically biased. See the cross-validated result\n")
cat("below for an honest performance estimate.) ===\n")
print(data.frame(Sample = ad_meta$Sample_ID, True_Group = ad_meta$Group,
                  Predicted_prob_AD = round(ad_pred_prob, 3)))

cat("\n=== PCOS sample AD-likeness scores (after ComBat harmonization) ===\n")
pcos_results <- data.frame(Subject_ID = pcos_meta$Subject_ID,
                            Diagnosis = pcos_meta$Diagnosis,
                            AD_likeness_prob = round(pcos_pred_prob, 3))
print(pcos_results)

cat("\n=== Summary statistics of PCOS AD-likeness scores ===\n")
cat("Min:", round(min(pcos_pred_prob), 3), "\n")
cat("Max:", round(max(pcos_pred_prob), 3), "\n")
cat("SD:", round(sd(pcos_pred_prob), 3), "\n")
cat("\nIf SD is near 0 and/or all scores cluster at 0 or 1, the projection is\n")
cat("still degenerate even after harmonization -- indicating the core issue\n")
cat("is sample size / model capacity mismatch, not platform scale alone.\n")
cat("If scores show real spread and relate plausibly to Diagnosis, this may\n")
cat("be a usable, reportable result.\n")

out_dir <- "results/ml_combat_attempt"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(pcos_results, file.path(out_dir, "PCOS_AD_likeness_after_ComBat.csv"), row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "PCOS_AD_likeness_after_ComBat.csv"), "\n")