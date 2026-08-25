import csv
import gzip
from pathlib import Path
from collections import defaultdict, Counter

BASE = Path(".")
RAW_DIR = BASE / "GSE84958_RAW"
MANIFEST = BASE / "output" / "GSE84958_sample_manifest_final.csv"

OUT_COUNTS = BASE / "output" / "GSE84958_subject_level_counts.csv"
OUT_META = BASE / "output" / "GSE84958_subject_metadata_final.csv"

def read_count_file(path):
    counts = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            gene = parts[0]
            try:
                count = int(float(parts[1]))
            except ValueError:
                continue
            counts[gene] = count
    return counts

# Read manifest
rows = []
with open(MANIFEST, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

included = [r for r in rows if r["Include_or_Exclude"] == "Include"]

# Group files by subject
subject_files = defaultdict(list)
subject_info = {}

for r in included:
    subject = r["Subject_code_guess"]
    subject_files[subject].append(r["Raw_file"])

    subject_info[subject] = {
        "Subject_ID": subject,
        "Label": r["Primary_label"],
        "Diagnosis": r["Diagnosis"],
        "Treatment": r["Treatment"],
        "Tissue": r["Tissue"],
        "Gender": r["Gender"],
        "Age": r["Age"]
    }

# Aggregate counts by summing replicate files per subject
subject_counts = {}
all_genes = set()

for subject, files in subject_files.items():
    summed = defaultdict(int)

    for fname in files:
        fp = RAW_DIR / fname
        counts = read_count_file(fp)

        for gene, count in counts.items():
            summed[gene] += count
            all_genes.add(gene)

    subject_counts[subject] = dict(summed)

# Sort genes and subjects
genes = sorted(all_genes)
subjects = sorted(subject_counts.keys())

# Save subject-level count matrix
with open(OUT_COUNTS, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Gene"] + subjects)

    for gene in genes:
        row = [gene]
        for subject in subjects:
            row.append(subject_counts[subject].get(gene, 0))
        writer.writerow(row)

# Save subject metadata
with open(OUT_META, "w", newline="", encoding="utf-8") as f:
    fieldnames = [
        "Subject_ID",
        "Label",
        "Diagnosis",
        "Treatment",
        "Tissue",
        "Gender",
        "Age",
        "Included_File_Count",
        "Raw_files"
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    for subject in subjects:
        info = subject_info[subject]
        writer.writerow({
            **info,
            "Included_File_Count": len(subject_files[subject]),
            "Raw_files": "; ".join(subject_files[subject])
        })

print("Done.")
print("Subject-level count matrix saved:", OUT_COUNTS)
print("Subject metadata saved:", OUT_META)
print("Genes:", len(genes))
print("Subjects:", len(subjects))
print("Label counts:", Counter(subject_info[s]["Label"] for s in subjects))
print("Included file counts per subject:")
for s in subjects:
    print(s, subject_info[s]["Label"], "files:", len(subject_files[s]))