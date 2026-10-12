#!/usr/bin/env python3
"""Export a patient-level audit of efferocytosis-related proteins.

The delivered XLSX files contain an incorrect worksheet dimension (A1), so
LibreOffice is used only as a read-only converter before CSV parsing.
"""

from __future__ import annotations

import os

import csv
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(os.environ.get("MPP_PROTEOMICS_ROOT", "inputs/proteomics"))
OUT = Path("inputs/workspace/工作文件/蛋白组胞葬相关蛋白审计_2026-08-19.csv")

COMPARISONS = {
    "plastic_GroupS_vs_GroupN": ROOT
    / "支原体肺炎BALF蛋白组塑形VS非塑形/WOSP25189_20250521-塑形非塑预后/附件/4-Differential_protein/diffprotein_table/unadjDiff_GroupS_GroupN_all.xlsx",
    "plastic_GroupSE_vs_GroupNE": ROOT
    / "支原体肺炎BALF蛋白组塑形VS非塑形/WOSP25189_20250521-塑形非塑预后/附件/4-Differential_protein/diffprotein_table/unadjDiff_GroupSE_GroupNE_all.xlsx",
    "plastic_GroupSE_vs_GroupSF": ROOT
    / "支原体肺炎BALF蛋白组塑形VS非塑形/WOSP25189_20250521-塑形非塑预后/附件/4-Differential_protein/diffprotein_table/unadjDiff_GroupSE_GroupSF_all.xlsx",
    "plastic_GroupSF_vs_GroupNF": ROOT
    / "支原体肺炎BALF蛋白组塑形VS非塑形/WOSP25189_20250521-塑形非塑预后/附件/4-Differential_protein/diffprotein_table/unadjDiff_GroupSF_GroupNF_all.xlsx",
    "plastic_GroupNE_vs_GroupNF": ROOT
    / "支原体肺炎BALF蛋白组塑形VS非塑形/WOSP25189_20250521-塑形非塑预后/附件/4-Differential_protein/diffprotein_table/unadjDiff_GroupNE_GroupNF_all.xlsx",
    "severity_Group1_vs_Group2": ROOT
    / "支原体肺炎BALF蛋白组轻症VS重症/WOSP25189_20250821-轻重症预后/附件/4-Differential_protein/diffprotein_table/unadjDiff_Group1_Group2_all.xlsx",
}

TARGETS = {
    "CD36",
    "MERTK",
    "AXL",
    "GAS6",
    "MFGE8",
    "THBS1",
    "C1QA",
    "C1QB",
    "C1QC",
    "TREM2",
    "STAB1",
    "STAB2",
    "SIRPA",
    "CD47",
    "LAMP1",
    "CTSD",
    "TGFB1",
    "HIF1A",
    "APOE",
    "LGALS3",
    "ITGAV",
    "ITGB3",
}


def convert(path: Path, outdir: Path) -> Path:
    subprocess.run(
        ["soffice", "--headless", "--convert-to", "csv", "--outdir", str(outdir), str(path)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return outdir / f"{path.stem}.csv"


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="protein_efferocytosis_") as temp:
        tempdir = Path(temp)
        rows = []
        for comparison, xlsx in COMPARISONS.items():
            csv_path = convert(xlsx, tempdir)
            with csv_path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for record in reader:
                    gene = (record.get("Gene") or "").strip().upper()
                    if gene not in TARGETS:
                        continue
                    rows.append(
                        {
                            "comparison": comparison,
                            "gene": gene,
                            "accession": record.get("Accession", ""),
                            "log2FC_group1_vs_group2": record.get("log2(foldchange)", ""),
                            "p_value": record.get("P_value", ""),
                            "FDR": record.get("P_value_adjust", ""),
                            "nominal_class": record.get("Class", ""),
                        }
                    )
        rows.sort(key=lambda x: (x["comparison"], x["gene"]))
        with OUT.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
