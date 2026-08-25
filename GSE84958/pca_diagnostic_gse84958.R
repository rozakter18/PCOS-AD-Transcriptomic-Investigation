library(DESeq2)
library(ggplot2)

OUT_DIR <- "output"

counts_file <- file.path(OUT_DIR, "GSE84958_subject_level_counts_filtered.csv")
meta_file   <- file.path(OUT_DIR, "GSE84958_subject_metadata_final.csv")
libsize_file <- file.path(OUT_DIR, "GSE84958_library_sizes.csv")

counts <- read.csv(counts_file, check.names = FALSE, row.names = 1)
meta   <- read.csv(meta_file, check.names = FALSE)
libsizes <- read.csv(libsize_file, check.names = FALSE)
colnames(libsizes)[1] <- "Subject_ID"

meta <- meta[match(colnames(counts), meta$Subject_ID), ]
if (!all(meta$Subject_ID == colnames(counts))) {
  stop("Metadata and count matrix columns are not aligned.")
}
meta <- merge(meta, libsizes, by = "Subject_ID", sort = FALSE)
meta <- meta[match(colnames(counts), meta$Subject_ID), ]

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

# variance-stabilizing transform for PCA (better behaved than raw/log counts
# for small-n, wide dynamic range data like this)
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

cat("\n=== PCA table (PC1/PC2 by subject, diagnosis, library size) ===\n")
print(pca_df[order(pca_df$Library_size_millions), ], row.names = FALSE)

cat("\n=== Correlation between PC1 and library size ===\n")
print(cor.test(pca_df$PC1, pca_df$Library_size_millions))

cat("\n=== Correlation between PC2 and library size ===\n")
print(cor.test(pca_df$PC2, pca_df$Library_size_millions))

# Plot 1: PCA colored by Diagnosis
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis, label = Subject_ID)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title = "GSE84958 PCA — colored by Diagnosis",
    x = paste0("PC1 (", pct_var[1], "% variance)"),
    y = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_minimal()

# Plot 2: PCA colored by library size (continuous)
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Library_size_millions, label = Subject_ID)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_viridis_c(name = "Library size\n(millions)") +
  labs(
    title = "GSE84958 PCA — colored by sequencing depth",
    x = paste0("PC1 (", pct_var[1], "% variance)"),
    y = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_minimal()

out_plot1 <- file.path(OUT_DIR, "GSE84958_PCA_by_diagnosis.png")
out_plot2 <- file.path(OUT_DIR, "GSE84958_PCA_by_librarysize.png")
ggsave(out_plot1, p1, width = 7, height = 5.5, dpi = 150)
ggsave(out_plot2, p2, width = 7, height = 5.5, dpi = 150)

cat("\nSaved plots:\n")
cat(out_plot1, "\n")
cat(out_plot2, "\n")

cat("\n=== Interpretation guide ===\n")
cat("If PC1 (the axis explaining the most variance) correlates strongly\n")
cat("with library size rather than separating cleanly by Diagnosis color\n")
cat("in the first plot, sequencing depth -- not PCOS biology -- is likely\n")
cat("the dominant source of variation, and the two low-depth PCOS samples\n")
cat("(MO10, MO11) should be reviewed as candidates for exclusion.\n")