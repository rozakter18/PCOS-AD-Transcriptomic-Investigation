import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path(".")
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(exist_ok=True)

MANIFEST = OUT_DIR / "GSE33000_sample_manifest_final_draft.csv"

raw_folder = BASE / "GSE33000_raw_data.txt"
possible_files = list(raw_folder.glob("*.txt")) + list(raw_folder.glob("*"))

if not possible_files:
    raise FileNotFoundError("No expression file found")

EXPR_FILE = possible_files[0]
print("Using expression file:", EXPR_FILE)

meta = pd.read_csv(MANIFEST)

included = meta[meta["Include_or_Exclude"] == "Include"].copy()

pfc_cols = included["Expression_column"].tolist()
ref_cols = included["Source_name_ch1_reference"].tolist()

print("Included samples:", included.shape[0])
print("Group counts:")
print(included["Group"].value_counts())

usecols = ["reporterID", "Gene", "Transcript"] + ref_cols + pfc_cols

print("\nReading paired reference + sample columns...")
df = pd.read_csv(EXPR_FILE, sep="\t", usecols=usecols)

print("Raw selected shape:", df.shape)

missing_ref = set(ref_cols) - set(df.columns)
missing_pfc = set(pfc_cols) - set(df.columns)

print("Missing reference columns:", missing_ref)
print("Missing PFC columns:", missing_pfc)

if missing_ref or missing_pfc:
    raise ValueError("Some required columns are missing.")

eps = 1e-6

ref_values = df[ref_cols].astype(float).to_numpy()
pfc_values = df[pfc_cols].astype(float).to_numpy()

logratio_values = np.log2((pfc_values + eps) / (ref_values + eps))

ratio_df = pd.DataFrame(logratio_values, columns=pfc_cols)
ratio_df.insert(0, "Transcript", df["Transcript"])
ratio_df.insert(0, "Gene", df["Gene"])
ratio_df.insert(0, "reporterID", df["reporterID"])

out_probe = OUT_DIR / "GSE33000_AD_Control_log2ratio_probe_level.csv"
ratio_df.to_csv(out_probe, index=False)

meta_ordered = included.set_index("Expression_column").loc[pfc_cols].reset_index()
meta_ordered.to_csv(OUT_DIR / "GSE33000_AD_Control_log2ratio_metadata_ordered.csv", index=False)

print("\nSaved probe-level log2-ratio matrix:")
print(out_probe)
print("Saved metadata:")
print("output/GSE33000_AD_Control_log2ratio_metadata_ordered.csv")

print("\nLog2-ratio matrix shape:", ratio_df.shape)

values = ratio_df[pfc_cols]

print("\nLog2-ratio value summary:")
print(values.stack().describe())

print("\nTotal NA values:", values.isna().sum().sum())

# Filter probes:
# keep probes with <=20% missing and non-zero variance
missing_fraction = values.isna().mean(axis=1)
probe_var = values.var(axis=1, skipna=True)

keep = (missing_fraction <= 0.20) & (probe_var > 0)

clean_df = ratio_df.loc[keep].copy()

out_clean = OUT_DIR / "GSE33000_AD_Control_log2ratio_probe_level_clean.csv"
clean_df.to_csv(out_clean, index=False)

print("\nClean probe-level shape:", clean_df.shape)
print("Saved clean probe-level matrix:")
print(out_clean)

# Gene-level collapse: keep highest-variance probe per primary gene symbol
def primary_gene(x):
    if pd.isna(x):
        return ""

    s = str(x).replace('"', "").strip()

    if not s:
        return ""

    # remove unannotated Rosetta-style IDs from gene-level collapse
    if s.startswith("RSE_"):
        return ""

    return s.split(",")[0].strip()

gene_df = clean_df.copy()
gene_df["Gene_primary"] = gene_df["Gene"].apply(primary_gene)
gene_df["ProbeVariance"] = gene_df[pfc_cols].var(axis=1, skipna=True)

gene_df = gene_df[gene_df["Gene_primary"] != ""].copy()

gene_df = gene_df.sort_values(
    by=["Gene_primary", "ProbeVariance"],
    ascending=[True, False]
)

gene_df = gene_df.drop_duplicates(subset=["Gene_primary"], keep="first")

# Put Gene_primary as first column
gene_level = gene_df[["Gene_primary", "reporterID", "Gene", "Transcript"] + pfc_cols].copy()
gene_level = gene_level.rename(columns={"Gene_primary": "Gene_symbol"})

out_gene = OUT_DIR / "GSE33000_AD_Control_log2ratio_gene_level_highestvar.csv"
gene_level.to_csv(out_gene, index=False)

print("\nGene-level collapsed shape:", gene_level.shape)
print("Saved gene-level matrix:")
print(out_gene)

print("\nFinal summary:")
print("AD:", (included["Group"] == "AD").sum())
print("Control:", (included["Group"] == "Control").sum())
print("Total included samples:", included.shape[0])