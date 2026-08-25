import csv
from pathlib import Path
from collections import Counter

BASE = Path(".")
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(exist_ok=True)

# Find series matrix txt file automatically
series_files = list(BASE.glob("GSE33000_series_matrix.txt/**/*.txt")) + list(BASE.glob("GSE33000_series_matrix.txt/*.txt"))

if not series_files:
    raise FileNotFoundError("Series matrix .txt file not found inside GSE33000_series_matrix.txt folder")

SERIES_FILE = series_files[0]
print("Using series matrix:", SERIES_FILE)

sample_data = {}
sample_ids = []

with open(SERIES_FILE, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith("!Sample_geo_accession"):
            parts = line.split("\t")[1:]
            sample_ids = [p.strip('"') for p in parts]
            for sid in sample_ids:
                sample_data[sid] = {"GSM_ID": sid, "Characteristics": []}

        elif line.startswith("!Sample_title"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Sample_title"] = value.strip('"')

        elif line.startswith("!Sample_source_name_ch1"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Source_name"] = value.strip('"')

        elif line.startswith("!Sample_characteristics_ch1"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Characteristics"].append(value.strip('"'))

        elif line.startswith("!Sample_description"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Description"] = value.strip('"')


def guess_group(text):
    t = text.lower()

    if "huntington" in t or "hd" in t:
        return "HD"
    if "alzheimer" in t or "ad" in t:
        return "AD"
    if "control" in t or "normal" in t or "non-demented" in t or "nondemented" in t:
        return "Control"

    return "Unknown"


rows = []

for gsm, info in sample_data.items():
    title = info.get("Sample_title", "")
    source = info.get("Source_name", "")
    desc = info.get("Description", "")
    chars = " | ".join(info.get("Characteristics", []))

    combined = " ".join([title, source, desc, chars])

    rows.append({
        "GSM_ID": gsm,
        "Sample_title": title,
        "Source_name": source,
        "Description": desc,
        "Characteristics": chars,
        "Group_guess": guess_group(combined)
    })

out_file = OUT_DIR / "GSE33000_sample_manifest_draft.csv"

with open(out_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("\nSaved draft manifest:", out_file)
print("Total samples:", len(rows))

print("\nGroup guess counts:")
print(Counter(r["Group_guess"] for r in rows))

print("\nUnique Source_name values:")
for value, count in Counter(r["Source_name"] for r in rows).items():
    print(f"{count} | {value}")

print("\nUnique Characteristics values:")
all_chars = []
for r in rows:
    if r["Characteristics"]:
        all_chars.extend(r["Characteristics"].split(" | "))

for value, count in Counter(all_chars).most_common(100):
    print(f"{count} | {value}")

print("\nFirst 10 rows:")
for r in rows[:10]:
    print("-----")
    print("GSM:", r["GSM_ID"])
    print("Title:", r["Sample_title"])
    print("Source:", r["Source_name"])
    print("Description:", r["Description"])
    print("Characteristics:", r["Characteristics"])
    print("Group guess:", r["Group_guess"])