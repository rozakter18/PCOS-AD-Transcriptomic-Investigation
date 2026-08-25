import pandas as pd
from pathlib import Path

COUNT_FILE = Path("output/GSE84958_subject_level_counts.csv")
META_FILE = Path("output/GSE84958_subject_metadata_final.csv")
OUT_DIR = Path("output")
OUT_DIR.mkdir(exist_ok=True)

counts = pd.read_csv(COUNT_FILE)
meta = pd.read_csv(META_FILE)

print("Count matrix shape:", counts.shape)
print("Metadata shape:", meta.shape)

genes = counts["Gene"]
sample_cols = [c for c in counts.columns if c != "Gene"]

print("\nSamples in count matrix:", len(sample_cols))
print("Samples in metadata:", meta["Subject_ID"].nunique())

missing_in_meta = set(sample_cols) - set(meta["Subject_ID"])
missing_in_counts = set(meta["Subject_ID"]) - set(sample_cols)

print("\nMissing in metadata:", missing_in_meta)
print("Missing in counts:", missing_in_counts)

# Total counts/library size per sample
lib_sizes = counts[sample_cols].sum(axis=0).sort_values(ascending=False)
print("\nLibrary sizes:")
print(lib_sizes)

lib_sizes.to_csv(OUT_DIR / "GSE84958_library_sizes.csv", header=["Total_counts"])

# Gene filtering check
count_values = counts[sample_cols]

nonzero_genes = (count_values.sum(axis=1) > 0).sum()
print("\nGenes with total count > 0:", nonzero_genes)

# Simple low-expression filter:
# keep genes with count >= 10 in at least 3 subjects
keep = (count_values >= 10).sum(axis=1) >= 3
filtered = counts.loc[keep].copy()

print("Genes kept after filter count>=10 in at least 3 subjects:", filtered.shape[0])

filtered.to_csv(OUT_DIR / "GSE84958_subject_level_counts_filtered.csv", index=False)

# Save metadata ordered according to count matrix columns
meta_ordered = meta.set_index("Subject_ID").loc[sample_cols].reset_index()
meta_ordered.to_csv(OUT_DIR / "GSE84958_subject_metadata_ordered.csv", index=False)

print("\nSaved:")
print("output/GSE84958_library_sizes.csv")
print("output/GSE84958_subject_level_counts_filtered.csv")
print("output/GSE84958_subject_metadata_ordered.csv")