import csv
from collections import Counter
from pathlib import Path

MANIFEST = Path("output/GSE84958_sample_manifest_draft.csv")

rows = []
with open(MANIFEST, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

print("Total rows:", len(rows))

print("\nDisease counts:")
print(Counter(r["Disease_status_guess"] for r in rows))

print("\nUnique Source_name values:")
for value, count in Counter(r["Source_name"] for r in rows).items():
    print(f"{count} | {value}")

print("\nUnique Characteristics values:")
all_chars = []
for r in rows:
    chars = r["Characteristics"].split(" | ")
    all_chars.extend(chars)

for value, count in Counter(all_chars).items():
    print(f"{count} | {value}")

print("\nFirst 10 sample rows:")
for r in rows[:10]:
    print("-----")
    print("GSM:", r["GSM_ID"])
    print("Title:", r["Sample_title"])
    print("Source:", r["Source_name"])
    print("Characteristics:", r["Characteristics"])
    print("Label:", r["Disease_status_guess"])
    print("Raw:", r["Raw_file"])