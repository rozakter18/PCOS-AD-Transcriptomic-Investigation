import csv
from pathlib import Path
from collections import defaultdict, Counter

INFILE = Path("output/GSE84958_sample_manifest_final.csv")
OUTFILE = Path("output/GSE84958_subject_group_summary.csv")

rows = []
with open(INFILE, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

groups = defaultdict(list)

for r in rows:
    subject = r["Subject_code_guess"]
    groups[subject].append(r)

summary_rows = []

for subject, items in groups.items():
    diagnoses = sorted(set(i["Diagnosis"] for i in items))
    treatments = sorted(set(i["Treatment"] for i in items))
    labels = sorted(set(i["Primary_label"] for i in items))
    include_status = sorted(set(i["Include_or_Exclude"] for i in items))
    raw_files = "; ".join(i["Raw_file"] for i in items)
    condition_codes = "; ".join(i["Condition_code"] for i in items)
    gsm_ids = "; ".join(i["GSM_ID"] for i in items)

    included_items = [i for i in items if i["Include_or_Exclude"] == "Include"]

    summary_rows.append({
        "Subject_code_guess": subject,
        "Total_files": len(items),
        "Included_files": len(included_items),
        "Diagnoses": "; ".join(diagnoses),
        "Treatments": "; ".join(treatments),
        "Labels": "; ".join(labels),
        "Include_status": "; ".join(include_status),
        "Condition_codes": condition_codes,
        "GSM_IDs": gsm_ids,
        "Raw_files": raw_files
    })

with open(OUTFILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
    writer.writeheader()
    writer.writerows(summary_rows)

print("Saved:", OUTFILE)
print("Total subject codes:", len(summary_rows))

print("\nSubject-level diagnosis/treatment summary:")
combo_counts = Counter()
included_subject_label_counts = Counter()

for s in summary_rows:
    combo = (s["Diagnoses"], s["Treatments"], s["Included_files"])
    combo_counts[combo] += 1

    if int(s["Included_files"]) > 0:
        included_subject_label_counts[s["Labels"]] += 1

for k, v in combo_counts.items():
    print(k, "=>", v)

print("\nIncluded subject label counts:")
print(included_subject_label_counts)

print("\nSubjects with included files:")
for s in summary_rows:
    if int(s["Included_files"]) > 0:
        print(
            s["Subject_code_guess"],
            "| included files:", s["Included_files"],
            "| labels:", s["Labels"],
            "| treatments:", s["Treatments"],
            "| raw:", s["Raw_files"]
        )