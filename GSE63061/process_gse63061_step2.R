# ---- 1. Use the locally available GPL10558 manifest file (from GSE54248 folder) ----
# instead of GEOquery, which hit a broken URL for this platform accession.
manifest_path <- "../GSE54248/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"

if (!file.exists(manifest_path)) {
  # try R1 as fallback, or current-folder copy
  alt_paths <- c(
    "../GSE54248/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz",
    "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz",
    "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"
  )
  found <- alt_paths[file.exists(alt_paths)]
  if (length(found) == 0) stop("Could not find GPL10558 manifest file. Check path.")
  manifest_path <- found[1]
}
cat("Using manifest file:", manifest_path, "\n")

# Illumina manifest files have several header lines before the actual table;
# find the line starting with "Species" or containing "Probe_Id" to locate it
con <- gzfile(manifest_path, "r")
all_lines_peek <- readLines(con, n = 20)
close(con)
header_line_idx <- grep("Probe_Id", all_lines_peek)[1]
cat("Manifest header found at line:", header_line_idx, "\n")
cat("Header preview:", substr(all_lines_peek[header_line_idx], 1, 200), "\n")

manifest <- read.delim(gzfile(manifest_path), skip = header_line_idx - 1,
                        header = TRUE, check.names = FALSE, fill = TRUE,
                        stringsAsFactors = FALSE)
cat("Manifest table dimensions:", dim(manifest), "\n")
cat("Manifest columns:", paste(head(colnames(manifest), 20), collapse = ", "), "\n")

id_col <- if ("Probe_Id" %in% colnames(manifest)) "Probe_Id" else colnames(manifest)[grep("Probe.*Id", colnames(manifest), ignore.case = TRUE)[1]]
symbol_col <- if ("Symbol" %in% colnames(manifest)) "Symbol" else colnames(manifest)[grep("Symbol", colnames(manifest), ignore.case = TRUE)[1]]
cat("Using probe ID column:", id_col, "| gene symbol column:", symbol_col, "\n")

if (is.na(id_col) || is.na(symbol_col)) stop("Could not identify ID/Symbol columns in manifest.")

probe_to_gene <- setNames(as.character(manifest[[symbol_col]]), as.character(manifest[[id_col]]))

# ---- 2. Load the AD/CTL expression matrix + metadata from step 1 ----
expr <- readRDS("output/GSE63061_AD_CTL_expr_probelevel.rds")
meta <- read.csv("output/GSE63061_AD_CTL_metadata.csv", check.names = FALSE)
cat("\nLoaded expression matrix:", dim(expr), "\n")

# Map probes to gene symbols
gene_symbols <- probe_to_gene[rownames(expr)]
mapped <- !is.na(gene_symbols) & gene_symbols != ""
cat("Probes with a valid gene symbol mapping:", sum(mapped), "of", length(gene_symbols), "\n")

expr_mapped <- expr[mapped, ]
gene_symbols_mapped <- gene_symbols[mapped]

# Collapse multiple probes per gene: keep the probe with highest mean expression
row_means <- rowMeans(expr_mapped, na.rm = TRUE)
ord <- order(gene_symbols_mapped, -row_means)
expr_ordered <- expr_mapped[ord, ]
symbols_ordered <- gene_symbols_mapped[ord]
keep <- !duplicated(symbols_ordered)

expr_gene <- expr_ordered[keep, ]
rownames(expr_gene) <- symbols_ordered[keep]

cat("Final gene-level expression matrix:", dim(expr_gene), "\n")

# ---- 3. Save final outputs ----
saveRDS(expr_gene, "output/GSE63061_AD_CTL_expr_genelevel.rds")
write.csv(data.frame(Gene_symbol = rownames(expr_gene), expr_gene, check.names = FALSE),
          "output/GSE63061_AD_CTL_expr_genelevel.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("output/GSE63061_AD_CTL_expr_genelevel.rds\n")
cat("output/GSE63061_AD_CTL_expr_genelevel.csv\n")

cat("\nSample metadata (Status column) is in output/GSE63061_AD_CTL_metadata.csv\n")
cat("Ready for classifier + CV pipeline.\n")