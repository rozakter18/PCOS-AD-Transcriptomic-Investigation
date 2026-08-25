library(dplyr)

# ---- 1. Parse series matrix for sample title, GSM, and status ----
sm_file <- "GSE63061_series_matrix.txt"
lines <- readLines(sm_file, n = 45)

parse_row <- function(line) {
  # Split on tabs, strip quotes
  parts <- strsplit(line, "\t")[[1]]
  parts <- gsub('^"|"$', '', parts)
  parts[-1]  # drop the row label (first element)
}

title_line <- grep("^!Sample_title", lines, value = TRUE)
gsm_line   <- grep("^!Sample_geo_accession", lines, value = TRUE)
char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)

titles <- parse_row(title_line)
gsms   <- parse_row(gsm_line)

# Find which characteristics row holds "status:" (should be the first one, per what we saw)
status_row <- NULL
for (cl in char_lines) {
  vals <- parse_row(cl)
  if (any(grepl("^status:", vals))) {
    status_row <- vals
    break
  }
}
if (is.null(status_row)) stop("Could not find status row in series matrix characteristics.")

status <- sub("^status: ", "", status_row)

cat("Total samples in series matrix:", length(titles), "\n")
cat("Status value counts:\n")
print(table(status))

meta <- data.frame(Chip_ID = titles, GSM = gsms, Status = status, stringsAsFactors = FALSE)

# Keep only clean AD / CTL (drop MCI, borderline MCI, OTHER, CTL to AD, MCI to CTL)
meta_clean <- meta[meta$Status %in% c("AD", "CTL"), ]
cat("\nSamples retained (AD + CTL only):", nrow(meta_clean), "\n")
print(table(meta_clean$Status))

# ---- 2. Load normalized expression matrix, match to clean samples ----
cat("\nLoading normalized expression matrix (this may take a minute, ~148MB)...\n")
expr_raw <- read.delim("GSE63061_normalized.txt", check.names = FALSE)
cat("Raw expression table dimensions:", dim(expr_raw), "\n")

id_col <- colnames(expr_raw)[1]
cat("ID column name:", id_col, "\n")

dup_count <- sum(duplicated(expr_raw[[id_col]]))
cat("Duplicate probe IDs found:", dup_count, "\n")

if (dup_count > 0) {
  # Keep the row with highest mean expression per duplicated probe ID
  numeric_cols <- setdiff(colnames(expr_raw), id_col)
  row_means <- rowMeans(expr_raw[, numeric_cols], na.rm = TRUE)
  ord <- order(expr_raw[[id_col]], -row_means)
  expr_raw <- expr_raw[ord, ]
  expr_raw <- expr_raw[!duplicated(expr_raw[[id_col]]), ]
  cat("Dimensions after deduplication (kept highest-mean row per ID):", dim(expr_raw), "\n")
}

expr <- expr_raw[, setdiff(colnames(expr_raw), id_col)]
rownames(expr) <- expr_raw[[id_col]]
cat("Expression matrix dimensions (probes x samples):", dim(expr), "\n")

# Match on Chip_ID (matrix column names) - verify all clean samples are present
missing <- setdiff(meta_clean$Chip_ID, colnames(expr))
cat("Clean-sample Chip_IDs missing from expression matrix:", length(missing), "\n")
meta_clean <- meta_clean[meta_clean$Chip_ID %in% colnames(expr), ]

expr_clean <- expr[, meta_clean$Chip_ID]
cat("Final AD/CTL expression matrix dimensions:", dim(expr_clean), "\n")

# ---- 3. Save intermediate outputs ----
dir.create("output", showWarnings = FALSE)
write.csv(meta_clean, "output/GSE63061_AD_CTL_metadata.csv", row.names = FALSE)
saveRDS(expr_clean, "output/GSE63061_AD_CTL_expr_probelevel.rds")

cat("\nSaved:\n")
cat("output/GSE63061_AD_CTL_metadata.csv\n")
cat("output/GSE63061_AD_CTL_expr_probelevel.rds\n")
cat("\nNext step: map ILMN probe IDs to gene symbols using GPL10558 annotation.\n")