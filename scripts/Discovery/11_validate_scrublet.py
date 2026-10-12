#!/usr/bin/env python3
import csv
import gzip
import math
import os
import sys
import zipfile
from pathlib import Path


SAMPLES = ["D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"]
EXPECTED_COLUMNS = ["sample_id", "barcode", "doublet_score", "predicted_doublet"]
SUMMARY_COLUMNS = [
    "sample_id",
    "cells",
    "features",
    "expected_doublet_rate",
    "scrublet_threshold",
    "predicted_doublets",
    "predicted_doublet_rate",
    "median_doublet_score",
]

analysis_root = Path(os.environ.get("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna"))
audit_dir = analysis_root / "results" / "00_audit"
scrublet_dir = analysis_root / "results" / "01_qc" / "scrublet"
validation_csv = scrublet_dir / "scrublet_validation.csv"
validation_md = scrublet_dir / "scrublet_validation.md"


def read_manifest():
    path = audit_dir / "canonical_mex_manifest.csv"
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    filtered = {row["sample_id"]: row for row in rows if row["role"] == "filtered_rsec_mex"}
    if list(filtered) != SAMPLES:
        raise ValueError(f"Canonical filtered MEX order mismatch: {list(filtered)}")
    return filtered


def read_barcodes(manifest_row):
    mex_path = Path(manifest_row["path"])
    member = manifest_row["barcode_member"]
    with zipfile.ZipFile(mex_path, "r") as archive:
        with archive.open(member, "r") as raw:
            stream = gzip.GzipFile(fileobj=raw) if member.lower().endswith(".gz") else raw
            with stream:
                return [line.decode("utf-8").rstrip("\r\n").split("\t")[0] for line in stream]


def close_float(observed, expected, tolerance=1e-12):
    return math.isclose(float(observed), float(expected), rel_tol=tolerance, abs_tol=tolerance)


def main():
    checks = []

    def add(name, passed, detail):
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    try:
        manifest = read_manifest()
        summary_path = scrublet_dir / "scrublet_summary.csv"
        with summary_path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            summary_columns = reader.fieldnames
            summary_rows = list(reader)
        add("summary_columns", summary_columns == SUMMARY_COLUMNS, str(summary_columns))
        add("summary_sample_order", [row.get("sample_id") for row in summary_rows] == SAMPLES, f"rows={len(summary_rows)}")
        summary_by_sample = {row.get("sample_id"): row for row in summary_rows}

        for sample in SAMPLES:
            row = summary_by_sample.get(sample)
            if row is None:
                add(f"{sample}_summary_present", False, "missing summary row")
                continue

            canonical_barcodes = read_barcodes(manifest[sample])
            output_path = scrublet_dir / f"{sample}_scrublet.csv.gz"
            observed_barcodes = []
            predicted_doublets = 0
            scores = []
            sample_values_ok = True
            booleans_ok = True
            scores_ok = True
            duplicate_free = True

            try:
                with gzip.open(output_path, "rt", newline="", encoding="utf-8-sig") as handle:
                    reader = csv.DictReader(handle)
                    columns_ok = reader.fieldnames == EXPECTED_COLUMNS
                    for output_row in reader:
                        observed_barcodes.append(output_row.get("barcode", ""))
                        sample_values_ok &= output_row.get("sample_id") == sample
                        call = output_row.get("predicted_doublet", "").lower()
                        booleans_ok &= call in {"true", "false"}
                        predicted_doublets += int(call == "true")
                        try:
                            score = float(output_row.get("doublet_score", "nan"))
                            scores_ok &= math.isfinite(score) and score >= 0
                            scores.append(score)
                        except ValueError:
                            scores_ok = False
                duplicate_free = len(observed_barcodes) == len(set(observed_barcodes))
            except (OSError, EOFError, csv.Error) as error:
                columns_ok = False
                sample_values_ok = False
                booleans_ok = False
                scores_ok = False
                add(f"{sample}_gzip_readable", False, repr(error))

            add(f"{sample}_columns", columns_ok, str(EXPECTED_COLUMNS))
            add(f"{sample}_sample_ids", sample_values_ok, f"rows={len(observed_barcodes)}")
            add(f"{sample}_barcodes", observed_barcodes == canonical_barcodes, f"observed={len(observed_barcodes)} canonical={len(canonical_barcodes)}")
            add(f"{sample}_barcode_unique", duplicate_free, f"unique={len(set(observed_barcodes))}")
            add(f"{sample}_boolean_calls", booleans_ok, f"predicted={predicted_doublets}")
            add(f"{sample}_finite_scores", scores_ok and len(scores) == len(canonical_barcodes), f"scores={len(scores)}")

            expected_cells = len(canonical_barcodes)
            expected_features = int(manifest[sample]["features"])
            summary_counts_ok = int(row["cells"]) == expected_cells and int(row["features"]) == expected_features
            add(f"{sample}_summary_dimensions", summary_counts_ok, f"cells={row['cells']} features={row['features']}")
            add(f"{sample}_summary_doublets", int(row["predicted_doublets"]) == predicted_doublets, f"summary={row['predicted_doublets']} observed={predicted_doublets}")
            expected_rate = predicted_doublets / expected_cells
            add(f"{sample}_summary_rate", close_float(row["predicted_doublet_rate"], expected_rate), f"summary={row['predicted_doublet_rate']} observed={expected_rate}")
            numeric_summary_ok = all(math.isfinite(float(row[column])) for column in SUMMARY_COLUMNS[3:])
            add(f"{sample}_summary_numeric", numeric_summary_ok, "finite rate, threshold, count, and median")
    except Exception as error:
        add("validator_exception", False, repr(error))

    scrublet_dir.mkdir(parents=True, exist_ok=True)
    with validation_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"])
        writer.writeheader()
        writer.writerows(checks)
    validation_md.write_text(
        "# Scrublet validation\n\n" + "\n".join(f"- {row['status']}: {row['check']} ({row['detail']})" for row in checks) + "\n",
        encoding="utf-8",
    )

    failures = [row for row in checks if row["status"] != "PASS"]
    if failures:
        for row in failures:
            print(row, file=sys.stderr)
        raise SystemExit(1)
    print(f"Scrublet validation PASS: {len(SAMPLES)} samples")


if __name__ == "__main__":
    main()
