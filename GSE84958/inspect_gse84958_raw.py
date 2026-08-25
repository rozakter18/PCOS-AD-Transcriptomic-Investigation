import gzip
from pathlib import Path

RAW_DIR = Path("GSE84958_RAW")

files = sorted(RAW_DIR.glob("*.gz"))

print("Total raw files:", len(files))

for fp in files[:3]:
    print("\n==============================")
    print("File:", fp.name)
    print("==============================")

    with gzip.open(fp, "rt", encoding="utf-8", errors="ignore") as f:
        for i in range(10):
            line = f.readline()
            if not line:
                break
            print(line.rstrip())