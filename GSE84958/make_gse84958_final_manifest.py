import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

INFILE = Path("output/GSE84958_sample_manifest_draft.csv")
OUTFILE = Path("output/GSE84958_sample_manifest_final.csv")

def extract_field(characteristics, key):
    parts = [p.strip() for p in characteristics.split(" | ")]
    for p in parts:
        if p.lower().startswith(key.lower() + ":"):
            return p.split(":", 1)[1].strip()
    return ""

def parse_raw_filename(raw_file):
    # Example: GSM2254708_MO31-1_CPre.txt.gz
    # Returns sample_code=MO31-1, subject_code=MO31, condition_code=CPre
    if not raw_file:
        return "", "", ""

    name = raw_file.replace(".txt.gz", "")
    parts = name.split("_")

    if len(parts) < 3:
        return "", "", ""

    sample_code = parts[1]
    condition_code = parts[2]

    subject_code = sample_code.split("-")[0]

    return sample_code, subject_code, condition_code

rows = []

with open(INFILE, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)

    for r in reader:
        characteristics = r["Characteristics"]

        diagnosis = extract_field(characteristics, "diagnosis")
        treatment = extract_field(characteristics, "treatment")
        tissue = extract_field(characteristics, "tissue")
        gender = extract_field(characteristics, "gender")
        age = extract_field(characteristics, "age")

        sample_code, subject_code, condition_code = parse_raw_filename(r["Raw_file"])

        include = "Include"
        reason = ""

        if treatment.lower() != "none":
            include = "Exclude"
            reason = "DHEA-treated sample; excluded from primary untreated PCOS-vs-control analysis"

        if diagnosis.lower() not in ["normal", "pcos"]:
            include = "Exclude"
            reason = "Unclear diagnosis"

        r_out = {
            "GSM_ID": r["GSM_ID"],
            "Raw_file": r["Raw_file"],
            "Sample_code": sample_code,
            "Subject_code_guess": subject_code,
            "Condition_code": condition_code,
            "Diagnosis": diagnosis,
            "Treatment": treatment,
            "Tissue": tissue,
            "Gender": gender,
            "Age": age,
            "Primary_label": "Control" if diagnosis.lower() == "normal" else "PCOS",
            "Include_or_Exclude": include,
            "Exclusion_reason": reason,
            "Characteristics": characteristics
        }

        rows.append(r_out)

with open(OUTFILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("Saved:", OUTFILE)
print("Total samples:", len(rows))

print("\nDiagnosis counts:")
print(Counter(r["Diagnosis"] for r in rows))

print("\nTreatment counts:")
print(Counter(r["Treatment"] for r in rows))

print("\nDiagnosis x Treatment counts:")
combo = Counter((r["Diagnosis"], r["Treatment"]) for r in rows)
for k, v in combo.items():
    print(k, v)

print("\nPrimary include/exclude counts:")
print(Counter(r["Include_or_Exclude"] for r in rows))

print("\nIncluded label counts:")
included = [r for r in rows if r["Include_or_Exclude"] == "Include"]
print(Counter(r["Primary_label"] for r in included))

print("\nSubject_code duplicate check:")
subject_counts = Counter(r["Subject_code_guess"] for r in rows)
duplicates = {k: v for k, v in subject_counts.items() if v > 1}
print("Total subject codes:", len(subject_counts))
print("Duplicated subject codes:", len(duplicates))
for k, v in list(duplicates.items())[:20]:
    print(k, v)