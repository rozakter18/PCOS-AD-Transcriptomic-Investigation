import os
import csv
import gzip
from pathlib import Path

BASE = Path(".")
SERIES_FILE = BASE / "GSE84958_series_matrix.txt" / "GSE84958_series_matrix.txt"
RAW_DIR = BASE / "GSE84958_RAW"
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(exist_ok=True)

def parse_series_matrix(path):
    sample_data = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    sample_ids = []

    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("!Sample_geo_accession"):
            parts = line.split("\t")[1:]
            sample_ids = [p.strip('"') for p in parts]
            for sid in sample_ids:
                sample_data[sid] = {"GSM_ID": sid}

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
            key_values = [p.strip('"') for p in parts]

            for sid, value in zip(sample_ids, key_values):
                if "Characteristics" not in sample_data[sid]:
                    sample_data[sid]["Characteristics"] = []
                sample_data[sid]["Characteristics"].append(value)

    return sample_data


def list_raw_files(raw_dir):
    raw_files = list(raw_dir.glob("*.gz"))
    raw_map = {}

    for fp in raw_files:
        name = fp.name
        gsm = name.split("_")[0]
        raw_map[gsm] = name

    return raw_map


def guess_label(text):
    t = text.lower()

    if "pcos" in t or "polycystic" in t:
        return "PCOS"
    if "control" in t or "healthy" in t or "normal" in t:
        return "Control"
    return "Unknown"


def main():
    sample_data = parse_series_matrix(SERIES_FILE)
    raw_map = list_raw_files(RAW_DIR)

    rows = []

    for gsm, info in sample_data.items():
        title = info.get("Sample_title", "")
        source = info.get("Source_name", "")
        characteristics = " | ".join(info.get("Characteristics", []))
        combined_text = " ".join([title, source, characteristics])

        disease_status = guess_label(combined_text)

        raw_file = raw_map.get(gsm, "")

        rows.append({
            "GSM_ID": gsm,
            "Sample_title": title,
            "Source_name": source,
            "Characteristics": characteristics,
            "Disease_status_guess": disease_status,
            "Raw_file": raw_file,
            "Has_raw_file": "Yes" if raw_file else "No",
            "Include_or_Exclude": "",
            "Exclusion_reason": ""
        })

    out_file = OUT_DIR / "GSE84958_sample_manifest_draft.csv"

    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print("Done.")
    print(f"Total metadata samples: {len(sample_data)}")
    print(f"Total raw files: {len(raw_map)}")
    print(f"Manifest saved to: {out_file}")

    label_counts = {}
    for r in rows:
        label = r["Disease_status_guess"]
        label_counts[label] = label_counts.get(label, 0) + 1

    print("Disease status guess counts:")
    for k, v in label_counts.items():
        print(f"  {k}: {v}")

    missing_raw = [r["GSM_ID"] for r in rows if r["Has_raw_file"] == "No"]
    if missing_raw:
        print("Samples missing raw file:")
        print(missing_raw)


if __name__ == "__main__":
    main()