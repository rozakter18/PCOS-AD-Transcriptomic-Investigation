import pandas as pd
from pathlib import Path
from collections import Counter

BASE = Path(".")
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(exist_ok=True)

MANIFEST = OUT_DIR / "GSE33000_sample_manifest_final_draft.csv"

# Find expression file automatically
raw_folder = BASE / "GSE33000_raw_data.txt"
possible_files = list(raw_folder.glob("*.txt")) + list(raw_folder.glob("*"))

if not possible_files:
    raise FileNotFoundError("No expression file found inside GSE33000_raw_data.txt folder")

EXPR_FILE = possible_files[0]
print("Using expression file:", EXPR_FILE)

# Read manifest
meta = pd.read_csv(MANIFEST)

# Keep only AD and Control
included = meta[meta["Include_or_Exclude"] == "Include"].copy()

print("Total manifest samples:", meta.shape[0])
print("Included samples:", included.shape[0])
print("Included group counts:")
print(included["Group"].value_counts())

# PFC columns are actual sample columns
sample_cols = included["Expression_column"].tolist()

# Read only annotation + included PFC columns
usecols = ["reporterID", "Gene", "Transcript"] + sample_cols

print("\nReading selected expression columns...")
expr = pd.read_csv(EXPR_FILE, sep="\t", usecols=usecols)

print("Expression matrix shape:", expr.shape)

# Check missing columns
missing_cols = set(sample_cols) - set(expr.columns)
print("Missing expression columns:", missing_cols)

# Save ordered metadata
meta_ordered = included.set_index("Expression_column").loc[sample_cols].reset_index()
meta_ordered.to_csv(OUT_DIR / "GSE33000_AD_Control_metadata_ordered.csv", index=False)

# Save probe-level expression matrix
expr.to_csv(OUT_DIR / "GSE33000_AD_Control_PFC_expression_probe_level.csv", index=False)

print("\nSaved:")
print("output/GSE33000_AD_Control_metadata_ordered.csv")
print("output/GSE33000_AD_Control_PFC_expression_probe_level.csv")

print("\nFinal summary:")
print("Rows/probes:", expr.shape[0])
print("Sample columns:", len(sample_cols))
print("AD:", (included["Group"] == "AD").sum())
print("Control:", (included["Group"] == "Control").sum())