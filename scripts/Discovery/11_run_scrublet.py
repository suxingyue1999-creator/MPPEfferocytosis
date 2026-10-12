#!/usr/bin/env python3
import csv
import concurrent.futures
import gzip
import os
import zipfile
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp
import scrublet as scr


SAMPLES = ["D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"]
INPUT_ROOT = Path(os.environ.get("BD_INPUT_ROOT", "inputs/discovery/bd_rhapsody_10/work/outputs"))
ANALYSIS_ROOT = Path(os.environ.get("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna"))
OUTPUT_DIR = ANALYSIS_ROOT / "results" / "01_qc" / "scrublet"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def member(zf, basename):
    matches = [name for name in zf.namelist() if Path(name).name == basename]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {basename}, found {matches}")
    return matches[0]


def read_mex(path):
    with zipfile.ZipFile(path) as zf:
        matrix_name = member(zf, "matrix.mtx.gz")
        barcode_name = member(zf, "barcodes.tsv.gz")
        with zf.open(matrix_name) as raw, gzip.GzipFile(fileobj=raw) as stream:
            matrix = scipy.io.mmread(stream).tocsr().transpose().tocsc()
        with zf.open(barcode_name) as raw, gzip.GzipFile(fileobj=raw) as stream:
            barcodes = [line.decode("utf-8").rstrip("\n\r").split("\t")[0] for line in stream]
    if matrix.shape[0] != len(barcodes):
        raise RuntimeError(f"Barcode count {len(barcodes)} != matrix cells {matrix.shape[0]}")
    return matrix, barcodes


def run_sample(payload):
    index, sample_id = payload
    np.random.seed(20260726 + index)
    mex = INPUT_ROOT / sample_id / f"{sample_id}_RSEC_MolsPerCell_MEX.zip"
    print(f"Scrublet: loading {sample_id}", flush=True)
    counts, barcodes = read_mex(mex)
    expected_rate = min(0.12, max(0.04, 0.008 * counts.shape[0] / 1000.0))
    model = scr.Scrublet(
        counts,
        expected_doublet_rate=expected_rate,
        random_state=20260726 + index,
    )
    scores, calls = model.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=30,
        verbose=True,
    )
    threshold = float(model.threshold_)
    output = OUTPUT_DIR / f"{sample_id}_scrublet.csv.gz"
    with gzip.open(output, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "barcode", "doublet_score", "predicted_doublet"])
        for barcode, score, call in zip(barcodes, scores, calls):
            writer.writerow([sample_id, barcode, f"{float(score):.10g}", str(bool(call)).lower()])
    summary = {
        "sample_id": sample_id,
        "cells": counts.shape[0],
        "features": counts.shape[1],
        "expected_doublet_rate": expected_rate,
        "scrublet_threshold": threshold,
        "predicted_doublets": int(np.sum(calls)),
        "predicted_doublet_rate": float(np.mean(calls)),
        "median_doublet_score": float(np.median(scores)),
    }
    print(f"Scrublet: {sample_id} PASS ({np.sum(calls)}/{len(calls)} predicted)", flush=True)
    del counts, scores, calls, model
    return index, summary


def main():
    workers = max(1, min(int(os.environ.get("BD_SCRUBLET_WORKERS", "2")), len(SAMPLES)))
    print(f"Scrublet sample workers: {workers}", flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        indexed_rows = list(executor.map(run_sample, enumerate(SAMPLES)))
    summary_rows = [row for _, row in sorted(indexed_rows)]
    with (OUTPUT_DIR / "scrublet_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
