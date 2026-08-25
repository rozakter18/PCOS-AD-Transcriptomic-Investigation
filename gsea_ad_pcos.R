library(fgsea)
library(msigdbr)
library(dplyr)

# ---- Paths (run this script from the PCOS_AD root folder) ----
AD_RNK   <- "GSE33000/output/GSE33000_AD_vs_Control_ranked_genes_tstat.rnk"
PCOS_RNK <- "GSE84958/output/GSE84958_PCOS_vs_Control_ranked_genes_stat_excl_outliers.rnk"

RESULTS_DIR <- "results/pathway"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load ranked gene lists ----
load_rnk <- function(path) {
  df <- read.table(path, header = FALSE, sep = "\t",
                    col.names = c("gene", "stat"), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat) & !duplicated(df$gene), ]
  ranks <- df$stat
  names(ranks) <- df$gene
  sort(ranks, decreasing = TRUE)
}

ad_ranks   <- load_rnk(AD_RNK)
pcos_ranks <- load_rnk(PCOS_RNK)

cat("AD ranked genes:", length(ad_ranks), "\n")
cat("PCOS ranked genes:", length(pcos_ranks), "\n")

# ---- Load gene sets: Hallmark (broad processes) + KEGG (mechanistic pathways) ----
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")

# Fallback in case the KEGG subcategory name differs by msigdbr version
if (nrow(kegg) == 0) {
  kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
}

pathway_db <- bind_rows(hallmark, kegg)
pathways_list <- split(pathway_db$gene_symbol, pathway_db$gs_name)

cat("Total pathways loaded (Hallmark + KEGG):", length(pathways_list), "\n")

# ---- Run fgsea for each disease ----
set.seed(42)
run_gsea <- function(ranks, label) {
  res <- fgsea(pathways = pathways_list, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
  res <- res[order(res$padj), ]
  res$leadingEdge <- sapply(res$leadingEdge, function(x) paste(x, collapse = ";"))
  out_file <- file.path(RESULTS_DIR, paste0("GSEA_", label, "_hallmark_kegg.csv"))
  write.csv(as.data.frame(res), out_file, row.names = FALSE)
  cat("\n", label, "GSEA — significant pathways (FDR < 0.05):",
      sum(res$padj < 0.05, na.rm = TRUE), "of", nrow(res), "tested\n")
  cat("Saved:", out_file, "\n")
  res
}

ad_gsea   <- run_gsea(ad_ranks, "AD")
pcos_gsea <- run_gsea(pcos_ranks, "PCOS")

cat("\n=== Top 10 AD pathways ===\n")
print(head(ad_gsea[, c("pathway", "NES", "padj", "size")], 10))

cat("\n=== Top 10 PCOS pathways ===\n")
print(head(pcos_gsea[, c("pathway", "NES", "padj", "size")], 10))

# ---- Compare: pathways enriched in BOTH diseases ----
# Use a relaxed threshold (padj < 0.25, fgsea/GSEA convention for exploratory
# comparison — the standard threshold recommended by the Broad GSEA guidelines
# for hypothesis-generating work, distinct from the strict 0.05 used for
# single-pathway confirmation) since the PCOS side is likely underpowered.
ad_sub   <- ad_gsea[ad_gsea$padj < 0.25 & !is.na(ad_gsea$padj), c("pathway", "NES", "padj")]
pcos_sub <- pcos_gsea[pcos_gsea$padj < 0.25 & !is.na(pcos_gsea$padj), c("pathway", "NES", "padj")]

colnames(ad_sub)   <- c("pathway", "NES_AD", "padj_AD")
colnames(pcos_sub) <- c("pathway", "NES_PCOS", "padj_PCOS")

shared <- merge(ad_sub, pcos_sub, by = "pathway")
shared$Concordant_direction <- sign(shared$NES_AD) == sign(shared$NES_PCOS)
shared <- shared[order(shared$padj_AD), ]

out_shared <- file.path(RESULTS_DIR, "GSEA_shared_pathways_AD_PCOS.csv")
write.csv(shared, out_shared, row.names = FALSE)

cat("\n=== Pathways enriched in BOTH AD and PCOS (padj < 0.25 each) ===\n")
print(shared)
cat("\nSaved:", out_shared, "\n")

cat("\nConcordant direction = TRUE means the pathway trends the same way\n")
cat("(both up or both down) in AD and PCOS -- the strongest candidates\n")
cat("for a genuinely shared mechanism. Concordant = FALSE means the\n")
cat("pathway is enriched in both diseases but trending in opposite\n")
cat("directions, which is a weaker or more complex signal.\n")