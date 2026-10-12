#!/usr/bin/env python3
import os
import sys
from pathlib import Path

import pandas as pd


SAMPLES = ["D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"]
analysis_root = Path(os.environ.get("BD_ANALYSIS_ROOT", "/analysis"))
audit_dir = analysis_root / "results" / "00_audit"

manifest = pd.read_csv(audit_dir / "input_manifest.csv")
mex = pd.read_csv(audit_dir / "mex_h5ad_audit.csv")
rds = pd.read_csv(audit_dir / "seurat_rds_audit.csv")
metrics = pd.read_csv(audit_dir / "pipeline_metrics_by_sample.csv")
canonical = pd.read_csv(audit_dir / "canonical_mex_manifest.csv")
equivalence = pd.read_csv(audit_dir / "mex_count_equivalence.csv")
report_metrics = pd.read_csv(audit_dir / "pipeline_report_metrics.csv")
report_crosscheck = pd.read_csv(audit_dir / "pipeline_report_crosscheck.csv")
metadata = pd.read_csv(audit_dir / "10samples_metadata_locked.csv")
issues = pd.read_csv(audit_dir / "input_audit_issues.csv")

checks = []


def add(check, passed, detail):
    checks.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})


add("manifest_60_files", len(manifest) == 60, f"rows={len(manifest)}")
add("manifest_sha256_complete", manifest["sha256"].str.fullmatch(r"[0-9a-f]{64}").all(), "SHA-256 populated")
add("metadata_10_unique_samples", metadata["sample_id"].tolist() == SAMPLES, "fixed sample order")
add("metadata_5v5_severity", metadata.groupby("severity_en").size().to_dict() == {"mild": 5, "severe": 5}, "five mild/five severe")
add("mex_h5ad_10_samples", mex["sample_id"].tolist() == SAMPLES, f"rows={len(mex)}")
add("mex_h5ad_all_checks", bool(mex["all_cross_format_checks_pass"].all()), "MEX/H5AD/metrics dimensions")
add("seurat_rds_10_samples", rds["sample_id"].tolist() == SAMPLES, f"rows={len(rds)}")
add("seurat_rds_all_checks", bool(rds[["cells_unique", "features_unique", "cells_equal_filtered_mex", "features_equal_filtered_mex"]].all().all()), "RDS/MEX dimensions")
add("pipeline_metrics_10_samples", metrics["sample_id"].tolist() == SAMPLES, f"rows={len(metrics)}")
add("canonical_mex_20_archives", len(canonical) == 20 and canonical.groupby("sample_id").size().eq(2).all(), f"rows={len(canonical)}")
add("canonical_filtered_mex_primary", canonical.query("role == 'filtered_rsec_mex'")["primary_counts_source"].astype(bool).all(), "filtered RSEC MEX locked as primary counts")
add("mex_count_equivalence_10_samples", equivalence["sample_id"].tolist() == SAMPLES, f"rows={len(equivalence)}")
add("mex_count_equivalence_all_pass", equivalence["matrix_values_equal"].astype(bool).all(), "filtered MEX equals corresponding unfiltered sparse submatrix")
add("pipeline_report_metrics_10_samples", report_metrics["sample_id"].tolist() == SAMPLES, f"rows={len(report_metrics)}")
add("pipeline_report_crosscheck_all_pass", report_crosscheck["status"].eq("PASS").all(), f"checks={len(report_crosscheck)}")
add("no_input_audit_issues", len(issues) == 0, f"issues={len(issues)}")

validation = pd.DataFrame(checks)
validation.to_csv(audit_dir / "00_audit_validation.csv", index=False)

total_cells = int(mex["filtered_cells"].sum())
total_unfiltered = int(mex["unfiltered_barcodes"].sum())
report = [
    "# BD Rhapsody 10-sample input audit",
    "",
    f"- Samples: {len(SAMPLES)} (5 mild, 5 severe)",
    f"- Filtered cells across samples: {total_cells:,}",
    f"- Unfiltered barcodes across samples: {total_unfiltered:,}",
    "- Count source: filtered/unfiltered RSEC MEX archives.",
    "- Full sparse-count check: filtered MEX equals its corresponding unfiltered MEX submatrix.",
    "- Cross-checks: Pipeline Report HTML, Metrics Summary, H5AD in read-only backed mode, and Seurat RDS.",
    "- Source pipeline outputs were mounted read-only.",
    "- Locked analysis metadata contain coded sample identifiers only.",
    "- This cohort remains independent of the existing 23-sample cohort.",
    "",
    "## Validation",
    "",
]
for row in checks:
    report.append(f"- {row['status']}: {row['check']} ({row['detail']})")
(audit_dir / "input_audit_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

if validation["status"].eq("FAIL").any():
    print(validation.to_string(index=False), file=sys.stderr)
    sys.exit(1)
print(f"Stage 00 audit PASS: {total_cells} filtered cells", flush=True)
