#!/usr/bin/env python3
import csv
import gzip
import os
import sys
from collections import Counter
from pathlib import Path


SAMPLES = ["D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"]
root = Path(os.environ.get("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna"))
audit = root / "results" / "00_audit"
result = root / "results" / "01_qc"
state = root / "state"


def read_csv(path, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


checks = []


def add(name, passed, detail):
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


required = [
    "10samples_qc_summary_pre_filter.csv",
    "10samples_qc_summary.csv",
    "qc_thresholds_proposed.csv",
    "qc_thresholds_approved.csv",
    "qc_thresholds.csv",
    "qc_guardrails.csv",
    "qc_approval.csv",
    "emptydrops_summary.csv",
    "ambient_top_genes.csv",
    "cell_qc_pre_filter.csv.gz",
    "cell_qc.csv.gz",
    "cells_keep.csv.gz",
    "qc_distributions.pdf",
    "qc_distributions.png",
    "qc_count_feature_scatter.pdf",
    "qc_retention.pdf",
    "emptydrops_concordance.pdf",
    "emptydrops_concordance.png",
    "qc_review_packet.md",
    "qc_report.md",
    "r_sessionInfo.txt",
]
summary_pre = read_csv(result / "10samples_qc_summary_pre_filter.csv")
summary = read_csv(result / "10samples_qc_summary.csv")
proposed = read_csv(result / "qc_thresholds_proposed.csv")
approved = read_csv(result / "qc_thresholds_approved.csv")
guardrails = read_csv(result / "qc_guardrails.csv")
approval = read_csv(result / "qc_approval.csv")
emptydrops = read_csv(result / "emptydrops_summary.csv")
audit_rows = read_csv(audit / "mex_h5ad_audit.csv")
cell_rows = read_csv(result / "cell_qc.csv.gz", gz=True)
prefilter_rows = read_csv(result / "cell_qc_pre_filter.csv.gz", gz=True)
keep_rows = read_csv(result / "cells_keep.csv.gz", gz=True)

explicit_review = (
    len(approval) == 1
    and approval[0].get("approval_mode") == "EXPLICIT_EXCEPTION_REVIEW"
    and approval[0].get("status") == "APPROVED"
)
if explicit_review:
    required.extend(["qc_review_decision.csv", "qc_sample_flags.csv", "qc_sensitivity_analysis_plan.csv", "qc_review_r_sessionInfo.txt"])
missing = [name for name in required if not (result / name).is_file() or (result / name).stat().st_size == 0]
add("required_outputs", not missing, f"missing={missing}")
add(
    "pending_review_resolved_by_approval",
    not (state / "01_qc_review.PENDING").exists() or explicit_review,
    "explicit reviewed approval may be validated before pending marker archival",
)

add("summary_fixed_sample_order", [row["sample_id"] for row in summary] == SAMPLES, f"rows={len(summary)}")
add("prefilter_summary_fixed_sample_order", [row["sample_id"] for row in summary_pre] == SAMPLES, f"rows={len(summary_pre)}")
add("thresholds_fixed_sample_order", [row["sample_id"] for row in proposed] == SAMPLES and [row["sample_id"] for row in approved] == SAMPLES, "proposed and approved")
add("severity_5v5", Counter(row["severity"] for row in summary) == {"mild": 5, "severe": 5}, str(Counter(row["severity"] for row in summary)))
auto_approval = (
    len(approval) == 1
    and approval[0].get("approval_mode") == "AUTO_GUARDRAILS"
    and approval[0].get("status") == "APPROVED"
    and int(approval[0].get("failed_guardrails", "-1")) == 0
)
add("approval_contract", auto_approval or explicit_review, str(approval))
failed_guardrails = [row for row in guardrails if row["status"] == "FAIL"]
add("guardrails_resolved", (not failed_guardrails and auto_approval) or explicit_review, f"failed={len(failed_guardrails)}")

if explicit_review:
    decisions = read_csv(result / "qc_review_decision.csv")
    decision_keys = sorted(f"{row['sample_id']}|{row['guardrail']}" for row in decisions)
    failure_keys = sorted(f"{row['sample_id']}|{row['guardrail']}" for row in failed_guardrails)
    decisions_ok = (
        decision_keys == failure_keys
        and all(row["decision"] == "APPROVED_EXCEPTION" for row in decisions)
        and int(approval[0]["failed_guardrails"]) == len(failed_guardrails)
    )
    add("explicit_exception_decisions", decisions_ok, f"decisions={len(decisions)} failures={len(failed_guardrails)}")
    sensitivity = read_csv(result / "qc_sensitivity_analysis_plan.csv")
    scenarios = {row["scenario_id"]: row["excluded_samples"] for row in sensitivity}
    sensitivity_ok = scenarios == {
        "primary_all_10": "",
        "sensitivity_exclude_NS4": "D04",
        "sensitivity_exclude_all_qc_outliers": "D03;D04;D07",
    }
    add("prespecified_sensitivity_plan", sensitivity_ok, str(scenarios))

cell_counts = Counter(row["sample_id"] for row in cell_rows)
prefilter_counts = Counter(row["sample_id"] for row in prefilter_rows)
keep_counts = Counter(row["sample_id"] for row in keep_rows)
audit_counts = {row["sample_id"]: int(row["filtered_cells"]) for row in audit_rows}
unfiltered_counts = {row["sample_id"]: int(row["unfiltered_barcodes"]) for row in audit_rows}
add("cell_rows_match_canonical_filtered_mex", cell_counts == Counter(audit_counts), f"cells={sum(cell_counts.values())}")
add("prefilter_rows_match_canonical_filtered_mex", prefilter_counts == Counter(audit_counts), f"cells={sum(prefilter_counts.values())}")
cell_ids = [row["cell_id"] for row in cell_rows]
add("cell_ids_unique", len(cell_ids) == len(set(cell_ids)), f"unique={len(set(cell_ids))}")

summary_counts_ok = True
retention_ok = True
for row in summary:
    sample = row["sample_id"]
    summary_counts_ok &= int(float(row["cells_input"])) == cell_counts[sample]
    summary_counts_ok &= int(float(row["cells_qc_keep"])) == keep_counts[sample]
    retention = float(row["retention_rate"])
    retention_ok &= 0.50 <= retention <= 1.0 and keep_counts[sample] >= 500
add("summary_counts_match_cells", summary_counts_ok, f"retained={sum(keep_counts.values())}")
add("retention_guardrails", retention_ok, "each sample retains >=50% and >=500 cells")

doublet_ok = all(row["predicted_doublet"].lower() in {"true", "false"} and float(row["doublet_score"]) >= 0 for row in cell_rows)
add("scrublet_scores_complete", doublet_ok, f"rows={len(cell_rows)}")
emptydrops_cell_ok = all(row["emptydrops_tested"].lower() in {"true", "false"} and row["emptydrops_significant"].lower() in {"true", "false"} for row in cell_rows)
add("emptydrops_cell_diagnostics_complete", emptydrops_cell_ok, f"rows={len(cell_rows)}")
ambient_ok = all(float(row["ambient_top20_percent"]) >= 0 for row in cell_rows)
add("ambient_diagnostic_complete", ambient_ok, "top-20 soup fraction available per cell")
rank_ok = all(float(row["barcode_knee_umi"]) > 0 and float(row["barcode_inflection_umi"]) > 0 for row in proposed)
add("barcode_rank_diagnostic_complete", rank_ok, "knee and inflection available for all samples")

emptydrops_order_ok = [row["sample_id"] for row in emptydrops] == SAMPLES
emptydrops_counts_ok = emptydrops_order_ok and all(
    int(row["unfiltered_barcodes"]) == unfiltered_counts[row["sample_id"]]
    and int(row["bd_filtered_cells"]) == audit_counts[row["sample_id"]]
    for row in emptydrops
)
add("emptydrops_summary_matches_mex", emptydrops_counts_ok, f"rows={len(emptydrops)}")
emptydrops_files = [result / "emptydrops" / f"{sample}_emptydrops.csv.gz" for sample in SAMPLES]
add("emptydrops_per_barcode_outputs", all(path.is_file() and path.stat().st_size > 0 for path in emptydrops_files), "10 gzip CSV files")

with (result / "01_qc_validation.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"])
    writer.writeheader()
    writer.writerows(checks)

lines = ["# Stage 01 validation", ""] + [f"- {row['status']}: {row['check']} ({row['detail']})" for row in checks]
(result / "01_qc_validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
if any(row["status"] != "PASS" for row in checks):
    for row in checks:
        print(row, file=sys.stderr)
    raise SystemExit(1)
print(f"Stage 01 QC PASS: {len(cell_rows)} audited; {len(keep_rows)} retained")
