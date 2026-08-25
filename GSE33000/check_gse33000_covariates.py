import pandas as pd
from pathlib import Path

META = Path("output/GSE33000_AD_Control_log2ratio_metadata_ordered.csv")

meta = pd.read_csv(META)

print("Shape:", meta.shape)

print("\nGroup counts:")
print(meta["Group"].value_counts())

print("\nGender counts by group:")
print(pd.crosstab(meta["Group"], meta["Gender"]))

meta["Age_numeric"] = pd.to_numeric(meta["Age"], errors="coerce")

print("\nAge summary by group:")
print(meta.groupby("Group")["Age_numeric"].describe())

print("\nMissing age:", meta["Age_numeric"].isna().sum())
print("Missing gender:", meta["Gender"].isna().sum())

# Save simplified design metadata for limma
design = meta[["Expression_column", "GSM_ID", "Group", "Age_numeric", "Gender"]].copy()
design = design.rename(columns={"Expression_column": "Sample_ID"})

design.to_csv("output/GSE33000_limma_design_metadata.csv", index=False)

print("\nSaved:")
print("output/GSE33000_limma_design_metadata.csv")