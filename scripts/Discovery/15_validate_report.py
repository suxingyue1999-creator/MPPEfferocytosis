#!/usr/bin/env python3
import csv
import os
import sys
from pathlib import Path


root = Path(os.environ.get("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna"))
result = root / "results" / "05_report"
state = root / "state"
checks = []


def add(check, passed, detail):
    checks.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})


stages = ["00_audit", "01_qc", "02_integrate", "03_annotate", "04_compare"]
missing_stages = [stage for stage in stages if not (state / f"{stage}.PASS").is_file()]
add("all_prior_stages_pass", not missing_stages, f"missing={missing_stages}")
with (result / "source_integrity_recheck.csv").open(newline="", encoding="utf-8-sig") as handle:
    source = list(csv.DictReader(handle))
add("source_integrity_60_pass", len(source) == 60 and all(row["status"] == "PASS" for row in source), f"rows={len(source)}")
required = [
    "final_analysis_report.md", "final_analysis_report.html", "downstream_output_manifest.csv",
    "10samples_reproducibility_bundle.tar.gz", "local_transfer_notes.md",
    "reproducibility/environment_lock.txt", "reproducibility/commands_to_rerun.md",
    "reproducibility/10samples_metadata_locked.csv",
]
missing = [name for name in required if not (result / name).is_file() or (result / name).stat().st_size == 0]
add("required_final_outputs", not missing, f"missing={missing}")
report = (result / "final_analysis_report.md").read_text(encoding="utf-8")
add("patient_replicate_statement", "patient sample is the statistical replicate" in report, "explicit pseudoreplication guardrail")
add("independent_cohort_statement", "independent BD10_2026 cohort only" in report and "23-sample cohort" in report, "cohort separation explicit")
add("local_only_statement", "performed locally in Windows/WSL" in report, "no external computation")
with (result / "downstream_output_manifest.csv").open(newline="", encoding="utf-8-sig") as handle:
    manifest = list(csv.DictReader(handle))
add("result_manifest_complete", len(manifest) >= 30 and all(len(row["sha256"]) == 64 for row in manifest), f"rows={len(manifest)}")

with (result / "05_report_validation.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"])
    writer.writeheader(); writer.writerows(checks)
(result / "05_report_validation.md").write_text("# Stage 05 validation\n\n" + "\n".join(f"- {row['status']}: {row['check']} ({row['detail']})" for row in checks) + "\n", encoding="utf-8")
if any(row["status"] != "PASS" for row in checks):
    for row in checks:
        print(row, file=sys.stderr)
    raise SystemExit(1)
print("Stage 05 final report PASS")
