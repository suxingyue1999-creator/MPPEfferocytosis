#!/usr/bin/env python3
import csv
import concurrent.futures
import gzip
import hashlib
import html
import os
import shutil
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(os.environ.get("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna"))
INPUT_ROOT = Path(os.environ.get("BD_INPUT_ROOT", "inputs/discovery/bd_rhapsody_10/work/outputs"))
SCRIPT_ROOT = Path(os.environ.get("BD_SCRIPT_ROOT", "inputs/workspace/BD_Rhapsody_setup/downstream"))
ENV_ROOT = Path(os.environ.get("BD_ENV_ROOT", str(ROOT / "env")))
HASH_WORKERS = max(1, min(int(os.environ.get("BD_HASH_WORKERS", "8")), 8))
RESULT = ROOT / "results" / "05_report"
REPRO = RESULT / "reproducibility"
RESULT.mkdir(parents=True, exist_ok=True)
REPRO.mkdir(parents=True, exist_ok=True)


def read_csv(path, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def markdown_table(rows, columns, limit=None):
    selected = rows if limit is None else rows[:limit]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in selected:
        lines.append("| " + " | ".join(str(row.get(column, "")).replace("|", "\\|") for column in columns) + " |")
    return lines


def verify_source(row):
    relative = Path(row["path"].removeprefix("/inputs/"))
    path = INPUT_ROOT / relative
    actual_size = path.stat().st_size if path.is_file() else -1
    actual_hash = sha256(path) if path.is_file() and actual_size == int(row["bytes"]) else "MISSING_OR_SIZE_CHANGED"
    return {
        "sample_id": row["sample_id"], "role": row["role"], "path": str(path),
        "expected_bytes": row["bytes"], "actual_bytes": actual_size,
        "expected_sha256": row["sha256"], "actual_sha256": actual_hash,
        "status": "PASS" if actual_size == int(row["bytes"]) and actual_hash == row["sha256"] else "FAIL",
    }


audit = read_csv(ROOT / "results" / "00_audit" / "input_manifest.csv")
with concurrent.futures.ThreadPoolExecutor(max_workers=HASH_WORKERS) as executor:
    source_checks = list(executor.map(verify_source, audit))
with (RESULT / "source_integrity_recheck.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(source_checks[0]))
    writer.writeheader(); writer.writerows(source_checks)
if any(row["status"] != "PASS" for row in source_checks):
    raise RuntimeError("One or more audited BD source files changed")

for path in sorted(SCRIPT_ROOT.iterdir()):
    if path.is_file() and path.suffix.lower() in {".py", ".r", ".sh", ".csv", ".md", ".ps1"}:
        shutil.copy2(path, REPRO / path.name)
shutil.copy2(ENV_ROOT / "environment_lock.txt", REPRO / "environment_lock.txt")
shutil.copy2(ENV_ROOT / "environment_validation.csv", REPRO / "environment_validation.csv")
shutil.copy2(ROOT / "results" / "00_audit" / "input_manifest.csv", REPRO / "input_manifest.csv")
shutil.copy2(ROOT / "results" / "00_audit" / "10samples_metadata_locked.csv", REPRO / "10samples_metadata_locked.csv")

commands = """# Exact stage commands

All commands run locally inside `Ubuntu-22.04` WSL. Source BD outputs remain read-only by analysis convention and are checksum-verified before final PASS.

```bash
bash inputs/workspace/BD_Rhapsody_setup/downstream/10_install_local_environment.sh
bash inputs/workspace/BD_Rhapsody_setup/downstream/run_stage_00_audit.sh
bash inputs/workspace/BD_Rhapsody_setup/downstream/run_stage_01_qc.sh
bash inputs/workspace/BD_Rhapsody_setup/downstream/run_stage_02_integrate.sh
bash inputs/workspace/BD_Rhapsody_setup/downstream/run_stage_03_annotate.sh
bash inputs/workspace/BD_Rhapsody_setup/downstream/run_stage_04_compare.sh
bash inputs/workspace/BD_Rhapsody_setup/downstream/run_stage_05_report.sh
```
"""
(REPRO / "commands_to_rerun.md").write_text(commands, encoding="utf-8")

with tarfile.open(RESULT / "10samples_reproducibility_bundle.tar.gz", "w:gz") as archive:
    archive.add(REPRO, arcname="10samples_reproducibility")

qc = read_csv(ROOT / "results" / "01_qc" / "10samples_qc_summary.csv")
annotations = read_csv(ROOT / "results" / "03_annotate" / "cell_type_counts_by_sample.csv")
de_summary = read_csv(ROOT / "results" / "04_compare" / "pseudobulk_de_summary.csv")
composition = read_csv(ROOT / "results" / "04_compare" / "cell_type_composition_comparison.csv")
composition.sort(key=lambda row: float(row["fdr"]))
metadata = read_csv(ROOT / "results" / "00_audit" / "10samples_metadata_locked.csv")

input_cells = sum(int(float(row["cells_input"])) for row in qc)
retained_cells = sum(int(float(row["cells_qc_keep"])) for row in qc)
broad_counts = Counter()
for row in annotations:
    broad_counts[row["cell_type_broad"]] += int(float(row["Freq"]))

report = [
    "# BD Rhapsody ten-sample single-cell analysis report",
    "",
    f"**Generated:** {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}",
    "",
    "## Scope and immutable design",
    "",
    "This report covers the independent BD10_2026 cohort only. D01 through D05 are mild and D06 through D10 are severe. `sample_id` is the only join key. Public analysis metadata contain coded sample identifiers only.",
    "",
    "All computation was performed locally in Windows/WSL. No genomic or expression data was uploaded to an external service. The 60 audited BD source artifacts were SHA-256 rechecked at finalization and remained unchanged.",
    "",
    "## Methods",
    "",
    "- Input counts: filtered and unfiltered RSEC MEX, cross-checked against H5AD, Seurat RDS, Metrics Summary, and Pipeline Report metadata.",
    "- QC: per-cell UMI, detected genes, mitochondrial/ribosomal/hemoglobin fractions; barcode-rank cell-call diagnostics; Scrublet doublet calls; SoupX ambient profiles; per-sample robust thresholds independent of severity.",
    "- Integration: Seurat v5 SCTransform within each patient, RPCA anchors across patients, and integrated values used only for PCA, UMAP, neighbors, and clustering.",
    "- Annotation: canonical marker evidence at cluster level, followed by independent reclustering of major lineages. BD experimental labels were post hoc references only.",
    "- Composition: patient-level proportions, severe-minus-mild effects, bootstrap 95% confidence intervals, Wilcoxon tests, and BH FDR.",
    "- Differential expression: raw RNA counts summed per patient and cell type. edgeR quasi-likelihood tests provide P values/FDR; limma-voom empirical-Bayes intervals provide 95% log2-fold-change CIs.",
    "- Pathways: GO Biological Process competitive tests with limma cameraPR and local org.Hs.eg.db mappings.",
    "- The patient sample is the statistical replicate throughout all group comparisons; individual cells are not independent replicates.",
    "",
    "## QC overview",
    "",
    f"- Input putative cells: {input_cells:,}",
    f"- Cells retained: {retained_cells:,} ({retained_cells / input_cells:.1%})",
    f"- Samples flagged as robust QC outliers but retained: {sum(row['sample_outlier_flag'].lower() == 'true' for row in qc)}",
    "",
]
report += markdown_table(qc, ["sample_id", "severity", "cells_input", "cells_qc_keep", "retention_rate", "median_nCount_RNA", "median_nFeature_RNA", "median_percent_mt", "predicted_doublet_rate", "soupx_rho"])
report += ["", "## Broad cell-type yield", ""]
report += markdown_table([{"cell_type": key, "cells": value} for key, value in broad_counts.most_common()], ["cell_type", "cells"])
report += ["", "## Patient-level composition comparison", ""]
report += markdown_table(composition, ["cell_type_broad", "mean_difference_severe_minus_mild", "bootstrap_ci_low", "bootstrap_ci_high", "rank_biserial", "wilcoxon_p_value", "fdr"])
report += ["", "## Pseudobulk differential-expression coverage", ""]
report += markdown_table(de_summary, ["cell_type_broad", "genes_tested", "fdr_lt_0_05", "fdr_lt_0_05_abs_logfc_ge_0_5", "n_mild_patients", "n_severe_patients"])
report += [
    "",
    "## Interpretation guardrails",
    "",
    "- There are five independent patients per group. Estimates and confidence intervals are primary; nominal or FDR significance alone is not sufficient biological evidence.",
    "- Cell type labels are automated marker-based annotations and require review of marker plots before publication or clinical interpretation.",
    "- Ambient RNA estimates are diagnostics. Primary counts and pseudobulk results were not SoupX-adjusted.",
    "- The integrated assay was never used for differential expression.",
    "- This cohort must not be merged with or matched by patient name to the separate 23-sample cohort without a new, explicitly designed cross-cohort analysis.",
    "",
    "## Reproducibility",
    "",
    "The `10samples_reproducibility_bundle.tar.gz` archive contains scripts, locked metadata, environment validation, exact commands, and the original input manifest. `downstream_output_manifest.csv` records result paths, sizes, and SHA-256 hashes. Large Seurat objects remain in their stage directories and are referenced by the manifest rather than duplicated.",
    "",
    "## Software references",
    "",
    "- Hao Y et al. Seurat v5: multimodal reference mapping and cellular mapping. Nature Biotechnology (2024).",
    "- Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq using regularized negative binomial regression. Genome Biology (2019).",
    "- Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets. Cell Systems (2019).",
    "- Young MD, Behjati S. SoupX removes ambient RNA contamination. GigaScience (2020).",
    "- Chen Y et al. edgeR quasi-likelihood methods. Nucleic Acids Research (2016).",
    "- Ritchie ME et al. limma powers differential expression analyses. Nucleic Acids Research (2015).",
]
markdown = "\n".join(report) + "\n"
(RESULT / "final_analysis_report.md").write_text(markdown, encoding="utf-8")

style = "body{font-family:Arial,sans-serif;max-width:1100px;margin:32px auto;line-height:1.5;color:#222}pre{white-space:pre-wrap;background:#f5f5f5;padding:20px;border:1px solid #ddd}"
html_report = f"<!doctype html><html><head><meta charset='utf-8'><title>BD10 single-cell report</title><style>{style}</style></head><body><pre>{html.escape(markdown)}</pre></body></html>"
(RESULT / "final_analysis_report.html").write_text(html_report, encoding="utf-8")

def hash_result(path):
    return {
        "relative_path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


result_paths = [
    path for path in sorted((ROOT / "results").rglob("*"))
    if path.is_file() and path != RESULT / "downstream_output_manifest.csv"
]
with concurrent.futures.ThreadPoolExecutor(max_workers=HASH_WORKERS) as executor:
    manifest_rows = list(executor.map(hash_result, result_paths))
with (RESULT / "downstream_output_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"])
    writer.writeheader(); writer.writerows(manifest_rows)

transfer = [
    "# Local transfer notes",
    "",
    "The analysis is complete on this machine. For inspection on a 48 GB Mac, transfer the final report, QC/UMAP PDFs, CSV tables, reproducibility bundle, and annotated Seurat RDS. Re-running SCTransform/RPCA is not required for exploration; loading the annotated object may still need substantial memory.",
    "",
    "Primary objects:",
    "- results/02_integrate/10samples_integrated_seurat_v5.rds",
    "- results/03_annotate/10samples_annotated_seurat_v5.rds",
    "- results/05_report/downstream_output_manifest.csv",
]
(RESULT / "local_transfer_notes.md").write_text("\n".join(transfer) + "\n", encoding="utf-8")
transfer_path = RESULT / "local_transfer_notes.md"
manifest_rows.append({
    "relative_path": str(transfer_path.relative_to(ROOT)),
    "bytes": transfer_path.stat().st_size,
    "sha256": sha256(transfer_path),
})
with (RESULT / "downstream_output_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"])
    writer.writeheader(); writer.writerows(manifest_rows)
