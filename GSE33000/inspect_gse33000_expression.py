from pathlib import Path
import pandas as pd

BASE = Path(".")
RAW_FOLDER = BASE / "GSE33000_raw_data.txt"

files = []
for pattern in ["*.txt", "*.csv", "*.tsv", "*.gz"]:
    files.extend(list(RAW_FOLDER.glob(pattern)))

if not files:
    files = list(RAW_FOLDER.glob("*"))

print("Files found:")
for f in files:
    print(f.name, "| size:", round(f.stat().st_size / (1024 * 1024), 2), "MB")

if not files:
    raise FileNotFoundError("No expression/raw file found inside GSE33000_raw_data.txt folder")

fp = files[0]
print("\nInspecting first file:", fp)

print("\nFirst 20 lines:")
with open(fp, "r", encoding="utf-8", errors="ignore") as f:
    for i in range(20):
        line = f.readline()
        if not line:
            break
        print(line.rstrip()[:500])

print("\nTrying to read as table...")

try:
    df = pd.read_csv(fp, sep="\t", nrows=5)
    print("Tab-separated shape preview:", df.shape)
    print(df.head())
except Exception as e:
    print("Tab-separated read failed:", e)

try:
    df2 = pd.read_csv(fp, sep=",", nrows=5)
    print("Comma-separated shape preview:", df2.shape)
    print(df2.head())
except Exception as e:
    print("Comma-separated read failed:", e)