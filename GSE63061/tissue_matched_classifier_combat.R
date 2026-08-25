library(sva)
library(pROC)
library(limma)

set.seed(42)

# ---- 1. Load AD-blood (GSE63061) ----
ad_expr <- readRDS("output/GSE63061_AD_CTL_expr_genelevel.rds")
ad_meta <- read.csv("output/GSE63061_AD_CTL_metadata.csv", check.names = FALSE)
ad_meta$Status <- factor(ad_meta$Status, levels = c("CTL", "AD"))
ad_expr <- as.matrix(ad_expr)

cat("AD-blood matrix:", dim(ad_expr), "\n")
cat("AD-blood missing values:", sum(is.na(ad_expr)), "\n")

# ---- 2. Re-derive PCOS-blood (GSE54248) expression matrix ----
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

cat("PCOS-blood matrix:", dim(pcos_expr), "\n")
cat("PCOS-blood missing values:", sum(is.na(pcos_expr)), "\n")

# ---- 3. Common genes, remove any with missing values in either matrix ----
common_genes <- intersect(rownames(ad_expr), rownames(pcos_expr))
ad_common <- ad_expr[common_genes, ]
pcos_common <- pcos_expr[common_genes, ]

complete_ad   <- rownames(ad_common)[complete.cases(ad_common)]
complete_pcos <- rownames(pcos_common)[complete.cases(pcos_common)]
complete_genes <- intersect(complete_ad, complete_pcos)
cat("\nCommon genes:", length(common_genes), "| complete-case genes:", length(complete_genes), "\n")

ad_common <- ad_common[complete_genes, ]
pcos_common <- pcos_common[complete_genes, ]

# ---- 4. ComBat harmonization (dataset as batch) ----
combined <- cbind(ad_common, pcos_common)
batch <- factor(c(rep("AD_blood", ncol(ad_common)), rep("PCOS_blood", ncol(pcos_common))))
cat("\nBatch counts:\n")
print(table(batch))

mod <- model.matrix(~1, data = data.frame(batch))
combined_harmonized <- ComBat(dat = combined, batch = batch, mod = mod)

ad_harmonized   <- combined_harmonized[, colnames(ad_common)]
pcos_harmonized <- combined_harmonized[, colnames(pcos_common)]

cat("\nHarmonization complete.\n")
cat("AD-blood harmonized range:", round(range(ad_harmonized), 2), "\n")
cat("PCOS-blood harmonized range:", round(range(pcos_harmonized), 2), "\n")
cat("AD-blood harmonized median:", round(median(ad_harmonized), 2), "\n")
cat("PCOS-blood harmonized median:", round(median(pcos_harmonized), 2), "\n")

# ---- 5. 5-fold stratified CV on harmonized AD-blood data ----
n_pcs <- 20
labels <- ad_meta$Status[match(colnames(ad_harmonized), ad_meta$Chip_ID)]

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

  pca_k <- prcomp(t(ad_harmonized[, train_idx]), scale. = TRUE, center = TRUE)
  train_pcs_k <- pca_k$x[, 1:n_pcs]
  test_centered_k <- scale(t(ad_harmonized[, test_idx]), center = pca_k$center, scale = pca_k$scale)
  test_pcs_k <- test_centered_k %*% pca_k$rotation[, 1:n_pcs]

  train_df_k <- data.frame(Status = labels[train_idx], train_pcs_k)
  model_k <- glm(Status ~ ., data = train_df_k, family = binomial)
  test_df_k <- as.data.frame(test_pcs_k)
  cv_pred_prob[test_idx] <- predict(model_k, newdata = test_df_k, type = "response")
}

cv_auc <- roc(response = labels, predictor = cv_pred_prob, levels = c("CTL", "AD"), quiet = TRUE)
cv_pred_class <- factor(ifelse(cv_pred_prob >= 0.5, "AD", "CTL"), levels = c("CTL", "AD"))
cv_accuracy <- mean(cv_pred_class == labels)

cat("\n=== AD-BLOOD classifier (post-ComBat): honest CV performance ===\n")
cat("CV AUC:", round(as.numeric(cv_auc$auc), 3), "\n")
cat("CV Accuracy:", round(cv_accuracy, 3), "\n")
print(table(True = labels, Predicted = cv_pred_class))

# ---- 6. Train final model on all AD-blood (harmonized), project PCOS-blood ----
pca_full <- prcomp(t(ad_harmonized), scale. = TRUE, center = TRUE)
ad_pcs_full <- pca_full$x[, 1:n_pcs]
train_df_full <- data.frame(Status = labels, ad_pcs_full)
model_full <- glm(Status ~ ., data = train_df_full, family = binomial)

pcos_centered <- scale(t(pcos_harmonized), center = pca_full$center, scale = pca_full$scale)
pcos_pcs <- pcos_centered %*% pca_full$rotation[, 1:n_pcs]
pcos_pred_prob <- predict(model_full, newdata = as.data.frame(pcos_pcs), type = "response")

cat("\n=== PCOS-blood AD-likeness scores (TISSUE-MATCHED + ComBat-harmonized) ===\n")
pcos_results <- data.frame(Sample = colnames(pcos_harmonized),
                             Diagnosis = pcos_diagnosis,
                             AD_likeness_prob = round(pcos_pred_prob, 3))
print(pcos_results)

cat("\nPCOS mean:", round(mean(pcos_pred_prob[pcos_diagnosis == "PCOS"]), 3), "\n")
cat("Control mean:", round(mean(pcos_pred_prob[pcos_diagnosis == "Control"]), 3), "\n")
cat("SD across all PCOS-blood samples:", round(sd(pcos_pred_prob), 3), "\n")

dir.create("../results/ml_blood_matched", recursive = TRUE, showWarnings = FALSE)
write.csv(pcos_results, "../results/ml_blood_matched/PCOS_blood_AD_likeness_combat.csv", row.names = FALSE)
cat("\nSaved: results/ml_blood_matched/PCOS_blood_AD_likeness_combat.csv\n")
