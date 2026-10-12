#!/usr/bin/env python3
import csv
import gzip
import hashlib
import html.parser
import gc
import os
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io


SAMPLES = ["D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"]
INPUT_ROOT = Path(os.environ.get("BD_INPUT_ROOT", "/inputs"))
ANALYSIS_ROOT = Path(os.environ.get("BD_ANALYSIS_ROOT", "/analysis"))
SCRIPT_ROOT = Path(os.environ.get("BD_SCRIPT_ROOT", "/scripts"))
AUDIT_DIR = ANALYSIS_ROOT / "results" / "00_audit"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def find_member(names, candidates):
    for candidate in candidates:
        matches = [name for name in names if name.lower().endswith(candidate.lower())]
        if len(matches) == 1:
            return matches[0]
    raise ValueError(f"Expected one of {candidates}; archive contains {names}")


def open_text_member(archive, member):
    raw = archive.open(member, "r")
    if member.lower().endswith(".gz"):
        return gzip.open(raw, mode="rt", encoding="utf-8")
    return raw


def read_text_lines(archive, member):
    with open_text_member(archive, member) as handle:
        lines = []
        for line in handle:
            if isinstance(line, bytes):
                line = line.decode("utf-8")
            lines.append(line.rstrip("\r\n"))
    return lines


def audit_mex(path):
    with zipfile.ZipFile(path, "r") as archive:
        names = [name for name in archive.namelist() if not name.endswith("/")]
        matrix_member = find_member(names, ["matrix.mtx.gz", "matrix.mtx"])
        feature_member = find_member(names, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
        barcode_member = find_member(names, ["barcodes.tsv.gz", "barcodes.tsv"])

        with open_text_member(archive, matrix_member) as handle:
            dimensions = None
            for line in handle:
                if isinstance(line, bytes):
                    line = line.decode("utf-8")
                if line.startswith("%"):
                    continue
                dimensions = [int(value) for value in line.strip().split()]
                break
        if dimensions is None or len(dimensions) != 3:
            raise ValueError(f"Invalid Matrix Market dimensions in {path}")

        features = read_text_lines(archive, feature_member)
        barcodes = read_text_lines(archive, barcode_member)
        bad_member = archive.testzip()
        if bad_member:
            raise ValueError(f"CRC failure in {path}: {bad_member}")

    rows, columns, nonzero = dimensions
    if rows != len(features) or columns != len(barcodes):
        raise ValueError(
            f"Dimension mismatch in {path}: matrix={rows}x{columns}, "
            f"features={len(features)}, barcodes={len(barcodes)}"
        )
    return {
        "matrix_rows": rows,
        "matrix_columns": columns,
        "matrix_nonzero": nonzero,
        "features": features,
        "barcodes": barcodes,
        "matrix_member": matrix_member,
        "feature_member": feature_member,
        "barcode_member": barcode_member,
    }


def read_mex_sparse(path, matrix_member):
    with zipfile.ZipFile(path, "r") as archive:
        raw = archive.open(matrix_member, "r")
        stream = gzip.open(raw, mode="rb") if matrix_member.lower().endswith(".gz") else raw
        with stream:
            return scipy.io.mmread(stream).tocsr()


def normalize_label(value):
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def parse_report_number(value):
    cleaned = value.strip().replace(",", "").replace("%", "")
    if not re.fullmatch(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", cleaned):
        return None
    return float(cleaned)


class ReportHTMLParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables = []
        self.table = None
        self.row = None
        self.cell = None
        self.cell_tag = None
        self.putative_candidates = []
        self.pending_putative = False

    def _finish_cell(self):
        if self.cell is not None and self.row is not None:
            value = re.sub(r"\s+", " ", " ".join(self.cell)).strip()
            self.row.append((self.cell_tag, value))
        self.cell = None
        self.cell_tag = None

    def _finish_row(self):
        self._finish_cell()
        if self.row and self.table is not None:
            self.table.append(self.row)
        self.row = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            if self.table is not None:
                self._finish_row()
                self.tables.append(self.table)
            self.table = []
        elif tag == "tr" and self.table is not None:
            self._finish_row()
            self.row = []
        elif tag in {"th", "td"} and self.table is not None:
            self._finish_cell()
            if self.row is None:
                self.row = []
            self.cell = []
            self.cell_tag = tag

    def handle_endtag(self, tag):
        if tag in {"th", "td"}:
            self._finish_cell()
        elif tag == "tr":
            self._finish_row()
        elif tag == "table" and self.table is not None:
            self._finish_row()
            self.tables.append(self.table)
            self.table = None

    def handle_data(self, data):
        stripped = re.sub(r"\s+", " ", data).strip()
        if not stripped:
            return
        if self.cell is not None:
            self.cell.append(stripped)
        if len(stripped) > 200:
            return
        normalized = normalize_label(stripped)
        if normalized == "putativecells":
            self.pending_putative = True
            return
        if self.pending_putative:
            number = parse_report_number(stripped)
            if number is not None:
                self.putative_candidates.append(number)
                self.pending_putative = False


REPORT_FIELDS = {
    "total_reads_fastq": "Total Reads in FASTQ",
    "pct_celllabel_umi_aligned_uniquely": "Pct CellLabel UMI Aligned Uniquely",
    "pct_reads_useful": "Pct Reads Useful",
    "annotated_transcriptome_pct": "Annotated Transcriptome Pct",
    "total_rsec_molecules": "Total RSEC Molecules",
    "sequencing_saturation_pct": "Sequencing Saturation",
    "pct_reads_from_putative_cells": "Pct Reads from Putative Cells",
    "mean_reads_per_cell": "Mean Reads per Cell",
    "median_reads_per_cell": "Median Reads per Cell",
    "mean_molecules_per_cell": "Mean Molecules per Cell",
    "median_molecules_per_cell": "Median Molecules per Cell",
    "mean_genes_per_cell": "Mean Bioproducts per Cell",
    "median_genes_per_cell": "Median Bioproducts per Cell",
    "total_genes_detected": "Total Bioproducts Detected",
}


def parse_pipeline_report(path):
    parser = ReportHTMLParser()
    with path.open("rb") as handle:
        blocked_tag = None
        for raw_line in handle:
            if len(raw_line) > 128 * 1024:
                head = raw_line[:4096].lower()
                tail = raw_line[-4096:].lower()
                if blocked_tag is not None:
                    closing_bytes = f"</{blocked_tag}>".encode("ascii")
                    end = tail.find(closing_bytes)
                    if end >= 0:
                        remainder = tail[end + len(closing_bytes):].decode("utf-8", errors="replace")
                        blocked_tag = None
                        parser.feed(remainder)
                else:
                    for tag in ("script", "style"):
                        opening_bytes = f"<{tag}".encode("ascii")
                        closing_bytes = f"</{tag}>".encode("ascii")
                        if opening_bytes in head and closing_bytes not in tail:
                            blocked_tag = tag
                            break
                continue
            line = raw_line.decode("utf-8", errors="replace")
            while line:
                lowered = line.lower()
                if blocked_tag is not None:
                    closing = f"</{blocked_tag}>"
                    end = lowered.find(closing)
                    if end < 0:
                        line = ""
                        continue
                    line = line[end + len(closing):]
                    blocked_tag = None
                    continue
                starts = [(lowered.find(f"<{tag}"), tag) for tag in ("script", "style")]
                starts = [(position, tag) for position, tag in starts if position >= 0]
                if not starts:
                    parser.feed(line)
                    line = ""
                    continue
                start, tag = min(starts)
                if start:
                    parser.feed(line[:start])
                closing = f"</{tag}>"
                end = lowered.find(closing, start)
                if end >= 0:
                    line = line[end + len(closing):]
                else:
                    blocked_tag = tag
                    line = ""
    parser.close()

    candidates = {field: [] for field in REPORT_FIELDS}
    target_by_label = {normalize_label(label): field for field, label in REPORT_FIELDS.items()}
    for table in parser.tables:
        header = None
        for row in table:
            values = [value for _, value in row]
            tags = [tag for tag, _ in row]
            if values and tags and all(tag == "th" for tag in tags):
                header = values
                continue
            if header is None or len(values) != len(header):
                continue
            for label, value in zip(header, values):
                field = target_by_label.get(normalize_label(label))
                number = parse_report_number(value)
                if field is not None and number is not None:
                    candidates[field].append(number)
    candidates["putative_cell_count"] = parser.putative_candidates

    metrics = {}
    candidate_text = {}
    for field, values in candidates.items():
        unique = sorted(set(values))
        candidate_text[field] = "|".join(f"{value:.12g}" for value in unique)
        metrics[field] = unique[0] if len(unique) == 1 else None
    return metrics, candidate_text


def comparison_tolerance(field, expected):
    if field in {"putative_cell_count", "total_reads_fastq", "total_rsec_molecules", "total_genes_detected"}:
        return 0.5
    if "pct" in field or field == "sequencing_saturation_pct":
        return 0.011
    return max(0.011, abs(float(expected)) * 1e-6)


def parse_metrics(path):
    lines = [line.strip("\r\n") for line in path.read_text(encoding="utf-8-sig").splitlines()]
    sections = {}
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if re.fullmatch(r"#[^#].*#", line):
            section = line.strip("#")
            cursor = index + 1
            while cursor < len(lines) and not lines[cursor].strip():
                cursor += 1
            if cursor + 1 < len(lines):
                header = next(csv.reader([lines[cursor]]))
                values = next(csv.reader([lines[cursor + 1]]))
                if len(header) == len(values):
                    sections[section] = dict(zip(header, values))
            index = cursor + 2
        else:
            index += 1
    return sections


def as_number(value):
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except ValueError:
        return value
    return int(number) if number.is_integer() else number


metadata = pd.read_csv(SCRIPT_ROOT / "10samples_metadata.csv")
metadata = metadata.loc[:, ["sample_id", "severity", "severity_en", "cohort_id"]]
if metadata["sample_id"].tolist() != SAMPLES:
    raise SystemExit("Metadata sample order/content does not match the fixed cohort")
if metadata["sample_id"].duplicated().any():
    raise SystemExit("Metadata sample IDs are not unique")
if metadata.groupby("severity_en").size().to_dict() != {"mild": 5, "severe": 5}:
    raise SystemExit("Metadata must contain exactly five mild and five severe samples")
metadata.to_csv(AUDIT_DIR / "10samples_metadata_locked.csv", index=False, encoding="utf-8-sig")

manifest_rows = []
sample_rows = []
metric_rows = []
canonical_mex_rows = []
mex_equivalence_rows = []
report_metric_rows = []
report_crosscheck_rows = []
issue_rows = []

for sample in SAMPLES:
    print(f"Auditing {sample}", flush=True)
    sample_dir = INPUT_ROOT / sample
    paths = {
        "filtered_mex": sample_dir / f"{sample}_RSEC_MolsPerCell_MEX.zip",
        "unfiltered_mex": sample_dir / f"{sample}_RSEC_MolsPerCell_Unfiltered_MEX.zip",
        "h5ad": sample_dir / f"{sample}.h5ad",
        "seurat_rds": sample_dir / f"{sample}_Seurat.rds",
        "metrics_summary": sample_dir / f"{sample}_Metrics_Summary.csv",
        "pipeline_report": sample_dir / f"{sample}_Pipeline_Report.html",
    }
    for role, path in paths.items():
        if not path.is_file():
            issue_rows.append({"sample_id": sample, "check": f"{role}_exists", "status": "FAIL", "detail": str(path)})
            continue
        manifest_rows.append(
            {
                "sample_id": sample,
                "role": role,
                "path": str(path),
                "bytes": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
                "sha256": sha256_file(path),
            }
        )

    missing = [role for role, path in paths.items() if not path.is_file()]
    if missing:
        continue

    filtered = audit_mex(paths["filtered_mex"])
    unfiltered = audit_mex(paths["unfiltered_mex"])
    for role, details, purpose in (
        ("filtered_rsec_mex", filtered, "cell_expression_counts"),
        ("unfiltered_rsec_mex", unfiltered, "emptyDrops_barcodeRanks_SoupX"),
    ):
        canonical_mex_rows.append(
            {
                "sample_id": sample,
                "role": role,
                "path": str(paths["filtered_mex"] if role == "filtered_rsec_mex" else paths["unfiltered_mex"]),
                "matrix_member": details["matrix_member"],
                "feature_member": details["feature_member"],
                "barcode_member": details["barcode_member"],
                "features": details["matrix_rows"],
                "barcodes": details["matrix_columns"],
                "nonzero": details["matrix_nonzero"],
                "primary_counts_source": role == "filtered_rsec_mex",
                "purpose": purpose,
            }
        )
    metrics = parse_metrics(paths["metrics_summary"])
    cells_metrics = metrics.get("Cells", {})
    sequencing = metrics.get("Sequencing Quality", {})
    library = metrics.get("Library Quality", {})
    alignment = metrics.get("Alignment Categories", {})
    molecules = metrics.get("Reads and Molecules", {})

    putative_cells = as_number(cells_metrics.get("Putative_Cell_Count"))
    unfiltered_feature_index = {feature: index for index, feature in enumerate(unfiltered["features"])}
    filtered_feature_positions = [
        unfiltered_feature_index.get(feature, -1) for feature in filtered["features"]
    ]
    feature_subset = all(position >= 0 for position in filtered_feature_positions)
    feature_relative_order_preserved = feature_subset and all(
        left < right
        for left, right in zip(filtered_feature_positions, filtered_feature_positions[1:])
    )
    filtered_subset = set(filtered["barcodes"]).issubset(set(unfiltered["barcodes"]))
    filtered_feature_unique = len(set(filtered["features"])) == len(filtered["features"])
    unfiltered_feature_unique = len(set(unfiltered["features"])) == len(unfiltered["features"])
    filtered_barcode_unique = len(set(filtered["barcodes"])) == len(filtered["barcodes"])
    unfiltered_barcode_unique = len(set(unfiltered["barcodes"])) == len(unfiltered["barcodes"])

    filtered_matrix = read_mex_sparse(paths["filtered_mex"], filtered["matrix_member"])
    unfiltered_matrix = read_mex_sparse(paths["unfiltered_mex"], unfiltered["matrix_member"])
    filtered_matrix_sum = float(filtered_matrix.sum())
    corresponding_matrix_sum = None
    mismatch_nonzero = None
    matrix_values_equal = False
    if all((
        feature_subset,
        filtered_subset,
        filtered_feature_unique,
        unfiltered_feature_unique,
        filtered_barcode_unique,
        unfiltered_barcode_unique,
    )):
        feature_lookup = {value: index for index, value in enumerate(unfiltered["features"])}
        barcode_lookup = {value: index for index, value in enumerate(unfiltered["barcodes"])}
        feature_positions = np.asarray([feature_lookup[value] for value in filtered["features"]], dtype=np.int64)
        barcode_positions = np.asarray([barcode_lookup[value] for value in filtered["barcodes"]], dtype=np.int64)
        corresponding = unfiltered_matrix[feature_positions, :][:, barcode_positions].tocsr()
        corresponding_matrix_sum = float(corresponding.sum())
        difference = (filtered_matrix - corresponding).tocsr()
        difference.eliminate_zeros()
        mismatch_nonzero = int(difference.nnz)
        matrix_values_equal = mismatch_nonzero == 0
        del corresponding, difference, feature_positions, barcode_positions
    mex_equivalence_rows.append(
        {
            "sample_id": sample,
            "filtered_features_unique": filtered_feature_unique,
            "unfiltered_features_unique": unfiltered_feature_unique,
            "filtered_barcodes_unique": filtered_barcode_unique,
            "unfiltered_barcodes_unique": unfiltered_barcode_unique,
            "filtered_features_in_unfiltered": feature_subset,
            "filtered_barcodes_in_unfiltered": filtered_subset,
            "filtered_matrix_sum": filtered_matrix_sum,
            "unfiltered_corresponding_matrix_sum": corresponding_matrix_sum,
            "mismatched_nonzero_entries": mismatch_nonzero,
            "matrix_values_equal": matrix_values_equal,
        }
    )
    del filtered_matrix, unfiltered_matrix
    gc.collect()
    checks = {
        "filtered_cells_equal_metrics": filtered["matrix_columns"] == putative_cells,
        "filtered_features_subset_unfiltered": feature_subset,
        "filtered_feature_relative_order_preserved": feature_relative_order_preserved,
        "filtered_barcodes_subset_unfiltered": filtered_subset,
        "filtered_mex_equals_unfiltered_submatrix": matrix_values_equal,
    }
    for check, passed in checks.items():
        if not passed:
            issue_rows.append({"sample_id": sample, "check": check, "status": "FAIL", "detail": "cross-format mismatch"})

    sample_rows.append(
        {
            "sample_id": sample,
            "filtered_cells": filtered["matrix_columns"],
            "unfiltered_barcodes": unfiltered["matrix_columns"],
            "features": filtered["matrix_rows"],
            "unfiltered_features": unfiltered["matrix_rows"],
            "filtered_nonzero": filtered["matrix_nonzero"],
            "unfiltered_nonzero": unfiltered["matrix_nonzero"],
            "h5ad_cells": None,
            "h5ad_features": None,
            "metrics_putative_cells": putative_cells,
            "filtered_features_subset_unfiltered": feature_subset,
            "filtered_feature_relative_order_preserved": feature_relative_order_preserved,
            "filtered_barcodes_subset_unfiltered": filtered_subset,
            "all_cross_format_checks_pass": all(checks.values()),
        }
    )
    metric_row = {
            "sample_id": sample,
            "total_reads_fastq": as_number(sequencing.get("Total_Reads_in_FASTQ")),
            "pct_reads_filtered_out": as_number(sequencing.get("Pct_Reads_Filtered_Out")),
            "pct_q30_filtered_r2": as_number(library.get("Pct_Q30_Bases_in_Filtered_R2")),
            "pct_celllabel_umi_aligned_uniquely": as_number(library.get("Pct_CellLabel_UMI_Aligned_Uniquely")),
            "pct_reads_useful": as_number(library.get("Pct_Reads_Useful")),
            "annotated_transcriptome_pct": as_number(alignment.get("Annotated_Transcriptome_Pct")),
            "introns_pct": as_number(alignment.get("Introns_Pct")),
            "intergenic_pct": as_number(alignment.get("Intergenic_Regions_Pct")),
            "total_rsec_molecules": as_number(molecules.get("Total_RSEC_Molecules")),
            "sequencing_saturation_pct": as_number(molecules.get("Sequencing_Saturation")),
            "putative_cell_count": putative_cells,
            "pct_reads_from_putative_cells": as_number(cells_metrics.get("Pct_Reads_from_Putative_Cells")),
            "mean_reads_per_cell": as_number(cells_metrics.get("Mean_Reads_per_Cell")),
            "median_reads_per_cell": as_number(cells_metrics.get("Median_Reads_per_Cell")),
            "mean_molecules_per_cell": as_number(cells_metrics.get("Mean_Molecules_per_Cell")),
            "median_molecules_per_cell": as_number(cells_metrics.get("Median_Molecules_per_Cell")),
            "mean_genes_per_cell": as_number(cells_metrics.get("Mean_Bioproducts_per_Cell")),
            "median_genes_per_cell": as_number(cells_metrics.get("Median_Bioproducts_per_Cell")),
            "total_genes_detected": as_number(cells_metrics.get("Total_Bioproducts_Detected")),
        }
    metric_rows.append(metric_row)

    print(f"Parsing Pipeline Report: {sample}", flush=True)
    report_metrics, report_candidates = parse_pipeline_report(paths["pipeline_report"])
    report_metric_rows.append({"sample_id": sample, **report_metrics})
    for field in ["putative_cell_count", *REPORT_FIELDS.keys()]:
        expected = metric_row.get(field)
        observed = report_metrics.get(field)
        tolerance = comparison_tolerance(field, expected) if isinstance(expected, (int, float)) else None
        passed = (
            isinstance(expected, (int, float))
            and observed is not None
            and abs(float(expected) - float(observed)) <= tolerance
        )
        report_crosscheck_rows.append(
            {
                "sample_id": sample,
                "field": field,
                "metrics_summary_value": expected,
                "pipeline_report_value": observed,
                "absolute_difference": abs(float(expected) - float(observed)) if observed is not None and isinstance(expected, (int, float)) else None,
                "tolerance": tolerance,
                "report_candidates": report_candidates.get(field, ""),
                "status": "PASS" if passed else "FAIL",
            }
        )
        if not passed:
            issue_rows.append(
                {
                    "sample_id": sample,
                    "check": f"pipeline_report_{field}",
                    "status": "FAIL",
                    "detail": f"metrics={expected}; report_candidates={report_candidates.get(field, '')}",
                }
            )

pd.DataFrame(manifest_rows).to_csv(AUDIT_DIR / "input_manifest.csv", index=False)
pd.DataFrame(sample_rows).to_csv(AUDIT_DIR / "mex_h5ad_audit.csv", index=False)
pd.DataFrame(metric_rows).to_csv(AUDIT_DIR / "pipeline_metrics_by_sample.csv", index=False)
pd.DataFrame(canonical_mex_rows).to_csv(AUDIT_DIR / "canonical_mex_manifest.csv", index=False)
pd.DataFrame(mex_equivalence_rows).to_csv(AUDIT_DIR / "mex_count_equivalence.csv", index=False)
pd.DataFrame(report_metric_rows).to_csv(AUDIT_DIR / "pipeline_report_metrics.csv", index=False)
pd.DataFrame(report_crosscheck_rows).to_csv(AUDIT_DIR / "pipeline_report_crosscheck.csv", index=False)
pd.DataFrame(issue_rows, columns=["sample_id", "check", "status", "detail"]).to_csv(
    AUDIT_DIR / "input_audit_issues.csv", index=False
)

if len(sample_rows) != 10 or issue_rows:
    print(f"Input audit failed: samples={len(sample_rows)}, issues={len(issue_rows)}", file=sys.stderr)
    sys.exit(1)
print("Input MEX/H5AD/metrics audit completed", flush=True)
