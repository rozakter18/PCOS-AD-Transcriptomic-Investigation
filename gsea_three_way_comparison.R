library(fgsea)
library(msigdbr)
library(dplyr)

# ---- Paths (run this script from the PCOS_AD root folder) ----
AD_RNK          <- "GSE33000/output/GSE33000_AD_vs_Control_ranked_genes_tstat.rnk"
PCOS_ADIPOSE_RNK <- "GSE84958/output/GSE84958_PCOS_vs_Control_ranked_genes_stat_excl_outliers.rnk"
PCOS_BLOOD_RNK   <- "GSE54248/output/GSE54248_PCOS_vs_Control_ranked_genes_tstat.rnk"

RESULTS_DIR <- "results/pathway"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

load_rnk <- function(path) {
  df <- read.table(path, header = FALSE, sep = "\t",
                    col.names = c("gene", "stat"), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat) & !duplicated(df$gene), ]
  ranks <- df$stat
  names(ranks) <- df$gene
  sort(ranks, decreasing = TRUE)
}

ad_ranks           <- load_rnk(AD_RNK)
pcos_adipose_ranks <- load_rnk(PCOS_ADIPOSE_RNK)
pcos_blood_ranks   <- load_rnk(PCOS_BLOOD_RNK)

cat("AD (brain) ranked genes:", length(ad_ranks), "\n")
cat("PCOS (adipose) ranked genes:", length(pcos_adipose_ranks), "\n")
cat("PCOS (blood) ranked genes:", length(pcos_blood_ranks), "\n")

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")
if (nrow(kegg) == 0) {
  kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
}
pathway_db <- bind_rows(hallmark, kegg)
pathways_list <- split(pathway_db$gene_symbol, pathway_db$gs_name)
cat("Total pathways loaded (Hallmark + KEGG):", length(pathways_list), "\n")

set.seed(42)
run_gsea <- function(ranks, label) {
  res <- fgsea(pathways = pathways_list, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
  res <- res[order(res$padj), ]
  res$leadingEdge <- sapply(res$leadingEdge, function(x) paste(x, collapse = ";"))
  out_file <- file.path(RESULTS_DIR, paste0("GSEA_", label, "_hallmark_kegg.csv"))
  write.csv(as.data.frame(res), out_file, row.names = FALSE)
  cat("\n", label, "GSEA — significant pathways (FDR < 0.05):",
      sum(res$padj < 0.05, na.rm = TRUE), "of", nrow(res), "tested\n")
  res
}

ad_gsea           <- run_gsea(ad_ranks, "AD_brain")
pcos_adipose_gsea <- run_gsea(pcos_adipose_ranks, "PCOS_adipose")
pcos_blood_gsea    <- run_gsea(pcos_blood_ranks, "PCOS_blood")

cat("\n=== Top 10 PCOS-blood pathways ===\n")
print(head(pcos_blood_gsea[, c("pathway", "NES", "padj", "size")], 10))

# ---- Three-way merge: relaxed threshold (padj < 0.25) for comparison ----
sub_ad      <- ad_gsea[ad_gsea$padj < 0.25 & !is.na(ad_gsea$padj), c("pathway", "NES", "padj")]
sub_adipose <- pcos_adipose_gsea[pcos_adipose_gsea$padj < 0.25 & !is.na(pcos_adipose_gsea$padj), c("pathway", "NES", "padj")]
sub_blood   <- pcos_blood_gsea[pcos_blood_gsea$padj < 0.25 & !is.na(pcos_blood_gsea$padj), c("pathway", "NES", "padj")]

colnames(sub_ad)      <- c("pathway", "NES_AD_brain", "padj_AD_brain")
colnames(sub_adipose) <- c("pathway", "NES_PCOS_adipose", "padj_PCOS_adipose")
colnames(sub_blood)   <- c("pathway", "NES_PCOS_blood", "padj_PCOS_blood")

# AD vs PCOS-adipose (original, tissue-mismatched comparison)
ad_vs_adipose <- merge(sub_ad, sub_adipose, by = "pathway")
ad_vs_adipose$Concordant <- sign(ad_vs_adipose$NES_AD_brain) == sign(ad_vs_adipose$NES_PCOS_adipose)

# AD vs PCOS-blood (new — closer to tissue-matched, since blood is systemic/circulating)
ad_vs_blood <- merge(sub_ad, sub_blood, by = "pathway")
ad_vs_blood$Concordant <- sign(ad_vs_blood$NES_AD_brain) == sign(ad_vs_blood$NES_PCOS_blood)

# PCOS-adipose vs PCOS-blood (does the PCOS signal itself replicate across tissues?)
adipose_vs_blood <- merge(sub_adipose, sub_blood, by = "pathway")
adipose_vs_blood$Concordant <- sign(adipose_vs_blood$NES_PCOS_adipose) == sign(adipose_vs_blood$NES_PCOS_blood)

# Full three-way overlap
three_way <- Reduce(function(x, y) merge(x, y, by = "pathway"), list(sub_ad, sub_adipose, sub_blood))
three_way$AD_vs_Adipose_concordant <- sign(three_way$NES_AD_brain) == sign(three_way$NES_PCOS_adipose)
three_way$AD_vs_Blood_concordant   <- sign(three_way$NES_AD_brain) == sign(three_way$NES_PCOS_blood)
three_way$Adipose_vs_Blood_concordant <- sign(three_way$NES_PCOS_adipose) == sign(three_way$NES_PCOS_blood)
three_way <- three_way[order(three_way$padj_AD_brain), ]

write.csv(ad_vs_adipose, file.path(RESULTS_DIR, "GSEA_AD_vs_PCOS_adipose.csv"), row.names = FALSE)
write.csv(ad_vs_blood, file.path(RESULTS_DIR, "GSEA_AD_vs_PCOS_blood.csv"), row.names = FALSE)
write.csv(adipose_vs_blood, file.path(RESULTS_DIR, "GSEA_PCOS_adipose_vs_blood.csv"), row.names = FALSE)
write.csv(three_way, file.path(RESULTS_DIR, "GSEA_three_way_comparison.csv"), row.names = FALSE)

cat("\n=== AD (brain) vs PCOS (adipose) — pathways in both (padj<0.25 each) ===\n")
cat("Total shared:", nrow(ad_vs_adipose), "| Concordant direction:", sum(ad_vs_adipose$Concordant), "\n")

cat("\n=== AD (brain) vs PCOS (blood) — pathways in both (padj<0.25 each) ===\n")
cat("Total shared:", nrow(ad_vs_blood), "| Concordant direction:", sum(ad_vs_blood$Concordant), "\n")
print(ad_vs_blood[, c("pathway", "NES_AD_brain", "NES_PCOS_blood", "Concordant")])

cat("\n=== PCOS (adipose) vs PCOS (blood) — do the two PCOS tissues agree with each other? ===\n")
cat("Total shared:", nrow(adipose_vs_blood), "| Concordant direction:", sum(adipose_vs_blood$Concordant), "\n")
print(adipose_vs_blood[, c("pathway", "NES_PCOS_adipose", "NES_PCOS_blood", "Concordant")])

cat("\n=== Three-way overlap (all three datasets, padj<0.25 each) ===\n")
print(three_way[, c("pathway", "NES_AD_brain", "NES_PCOS_adipose", "NES_PCOS_blood",
                     "AD_vs_Adipose_concordant", "AD_vs_Blood_concordant", "Adipose_vs_Blood_concordant")])

cat("\nSaved comparison files to:", RESULTS_DIR, "\n")
cat("\nInterpretation guide:\n")
cat("- A pathway concordant across AD_vs_Blood but NOT AD_vs_Adipose (or vice versa)\n")
cat("  suggests the earlier tissue mismatch (brain vs adipose) may have masked or\n")
cat("  distorted a real signal that blood (a more systemic, comparable tissue) reveals.\n")
cat("- A pathway concordant in Adipose_vs_Blood (PCOS agreeing with itself across\n")
cat("  tissues) is good evidence the PCOS-side finding is a real, replicable signal,\n")
cat("  independent of whether it matches AD.\n")