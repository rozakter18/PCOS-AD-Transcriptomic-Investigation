import pandas as pd
from pathlib import Path

OUT_DIR = Path("output")

EXPR_FILE = OUT_DIR / "GSE33000_AD_Control_PFC_expression_probe_level.csv"
META_FILE = OUT_DIR / "GSE33000_AD_Control_metadata_ordered.csv"

expr = pd.read_csv(EXPR_FILE)
meta = pd.read_csv(META_FILE)

sample_cols = meta["Expression_column"].tolist()

print("Expression shape:", expr.shape)
print("Metadata shape:", meta.shape)

print("\nGroup counts:")
print(meta["Group"].value_counts())

missing_in_expr = set(sample_cols) - set(expr.columns)
missing_in_meta = set(expr.columns[3:]) - set(sample_cols)

print("\nMissing in expression:", missing_in_expr)
print("Missing in metadata:", missing_in_meta)

# Missing values
na_count = expr[sample_cols].isna().sum().sum()
print("\nTotal NA values in expression:", na_count)

# Basic expression value distribution
values = expr[sample_cols]

print("\nExpression value summary:")
print(values.stack().describe())

# Gene annotation quality
print("\nGene annotation summary:")
print("Total rows/probes:", expr.shape[0])
print("Gene missing:", expr["Gene"].isna().sum())
print("Unique Gene entries:", expr["Gene"].nunique())

# Remove probes with all-zero or all-NA values
all_na = values.isna().all(axis=1).sum()
zero_var = (values.var(axis=1) == 0).sum()

print("\nAll-NA probes:", all_na)
print("Zero-variance probes:", zero_var)

# Save clean probe-level matrix excluding all-NA and zero-variance probes
keep = (~values.isna().all(axis=1)) & (values.var(axis=1) > 0)
expr_clean = expr.loc[keep].copy()

expr_clean.to_csv(OUT_DIR / "GSE33000_AD_Control_PFC_expression_probe_level_clean.csv", index=False)

print("\nSaved clean matrix:")
print("output/GSE33000_AD_Control_PFC_expression_probe_level_clean.csv")
print("Clean shape:", expr_clean.shape)