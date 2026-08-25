import csv
import re
from pathlib import Path
from collections import Counter

BASE = Path(".")
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(exist_ok=True)

series_files = list(BASE.glob("GSE33000_series_matrix.txt/**/*.txt")) + list(BASE.glob("GSE33000_series_matrix.txt/*.txt"))

if not series_files:
    raise FileNotFoundError("Series matrix .txt file not found")

SERIES_FILE = series_files[0]
print("Using series matrix:", SERIES_FILE)

sample_data = {}
sample_ids = []

def clean(x):
    return x.strip().strip('"')

with open(SERIES_FILE, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith("!Sample_geo_accession"):
            parts = line.split("\t")[1:]
            sample_ids = [clean(p) for p in parts]
            for sid in sample_ids:
                sample_data[sid] = {
                    "GSM_ID": sid,
                    "Characteristics_ch1": [],
                    "Characteristics_ch2": []
                }

        elif line.startswith("!Sample_title"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Sample_title"] = clean(value)

        elif line.startswith("!Sample_source_name_ch1"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Source_name_ch1"] = clean(value)

        elif line.startswith("!Sample_source_name_ch2"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Source_name_ch2"] = clean(value)

        elif line.startswith("!Sample_characteristics_ch1"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Characteristics_ch1"].append(clean(value))

        elif line.startswith("!Sample_characteristics_ch2"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Characteristics_ch2"].append(clean(value))

        elif line.startswith("!Sample_description"):
            parts = line.split("\t")[1:]
            for sid, value in zip(sample_ids, parts):
                sample_data[sid]["Description"] = clean(value)


def extract_field(characteristics_list, key):
    for item in characteristics_list:
        if item.lower().startswith(key.lower() + ":"):
            return item.split(":", 1)[1].strip()
    return ""


def normalize_group(disease_status):
    t = disease_status.lower()

    if "alzheimer" in t:
        return "AD"
    if "huntington" in t:
        return "HD"
    if "control" in t or "normal" in t or "non-demented" in t or "nondemented" in t:
        return "Control"
    return "Unknown"


def clean_age(age_text):
    m = re.search(r"\d+", age_text)
    return m.group(0) if m else ""


rows = []

for gsm, info in sample_data.items():
    ch1 = info.get("Characteristics_ch1", [])
    ch2 = info.get("Characteristics_ch2", [])

    disease_status = extract_field(ch2, "disease status")
    age = extract_field(ch2, "age")
    gender = extract_field(ch2, "gender")

    group = normalize_group(disease_status)

    include = ""
    reason = ""

    if group == "AD":
        include = "Include"
    elif group == "Control":
        include = "Include"
    elif group == "HD":
        include = "Exclude"
        reason = "Huntington's disease sample; excluded from AD-vs-control analysis"
    else:
        include = "Exclude"
        reason = "Unknown disease status"

    expression_col = info.get("Source_name_ch2", "")

    rows.append({
        "GSM_ID": gsm,
        "Sample_title": info.get("Sample_title", ""),
        "Expression_column": expression_col,
        "Source_name_ch1_reference": info.get("Source_name_ch1", ""),
        "Source_name_ch2_sample": info.get("Source_name_ch2", ""),
        "Description": info.get("Description", ""),
        "Disease_status": disease_status,
        "Group": group,
        "Age": clean_age(age),
        "Age_raw": age,
        "Gender": gender,
        "Characteristics_ch1": " | ".join(ch1),
        "Characteristics_ch2": " | ".join(ch2),
        "Include_or_Exclude": include,
        "Exclusion_reason": reason
    })

out_file = OUT_DIR / "GSE33000_sample_manifest_final_draft.csv"

with open(out_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("\nSaved:", out_file)
print("Total samples:", len(rows))

print("\nDisease status counts:")
print(Counter(r["Disease_status"] for r in rows))

print("\nGroup counts:")
print(Counter(r["Group"] for r in rows))

print("\nInclude/exclude counts:")
print(Counter(r["Include_or_Exclude"] for r in rows))

print("\nIncluded group counts:")
included = [r for r in rows if r["Include_or_Exclude"] == "Include"]
print(Counter(r["Group"] for r in included))

print("\nFirst 10 rows:")
for r in rows[:10]:
    print("-----")
    print("GSM:", r["GSM_ID"])
    print("Expression column:", r["Expression_column"])
    print("Disease:", r["Disease_status"])
    print("Group:", r["Group"])
    print("Age:", r["Age_raw"])
    print("Gender:", r["Gender"])
    print("Include:", r["Include_or_Exclude"])