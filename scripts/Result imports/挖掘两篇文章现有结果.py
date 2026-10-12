#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从现有单细胞、胞葬和蛋白组结果中提取两篇文章的写稿级证据。

本脚本只读取已有文件，所有输出写入工作文件夹，不修改原始数据。
"""

from __future__ import annotations

import os

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


WORK = Path(os.environ.get("MPP_WORKSPACE", "inputs/workspace")) / "工作文件"
SC_RESULTS = WORK / "23例单细胞重新分析/results"
OLD10 = Path("inputs/legacy discovery")
PROT_PLASTIC = Path("inputs/proteomics/支原体肺炎BALF蛋白组塑形VS非塑形")
PROT_SEVERITY = Path("inputs/proteomics/支原体肺炎BALF蛋白组轻症VS重症")

CANDIDATES = {
    "HIF1A", "THBS1", "CD36", "MERTK", "AXL", "GAS6", "MFGE8", "ITGAV", "ITGB3",
    "SIRPA", "CD47", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CTSD", "LAMP1",
    "NPC2", "LIPA", "STAB1", "TGFB1", "SMAD2", "SMAD3", "SMAD4", "MAPK1", "MAPK3",
    "NLRP3", "IL1B", "TNF", "S100A8", "S100A9", "S100A12", "CXCL8", "IL6", "LDHA",
    "HCAR3", "F13A1", "FABP4", "FCN1", "IL7R", "FOXP3", "NKG7", "GNLY", "GZMB",
    "MKI67", "FCGR3A", "CD274", "SOD2", "MMP12", "PTPN13", "ITGA2", "CSF1",
}

PROTEIN_COMPARISON_MAP = {
    "unadjDiff_GroupS_GroupN_all.xlsx": "plastic_GroupS_vs_GroupN",
    "unadjDiff_GroupSE_GroupNE_all.xlsx": "plastic_GroupSE_vs_GroupNE",
    "unadjDiff_GroupSF_GroupNF_all.xlsx": "plastic_GroupSF_vs_GroupNF",
    "unadjDiff_GroupSE_GroupSF_all.xlsx": "plastic_GroupSE_vs_GroupSF",
    "unadjDiff_GroupNE_GroupNF_all.xlsx": "plastic_GroupNE_vs_GroupNF",
    "unadjDiff_Group1_Group2_all.xlsx": "severity_Group1_vs_Group2",
}

PLASTIC_DIFF_FILES = set(PROTEIN_COMPARISON_MAP) - {"unadjDiff_Group1_Group2_all.xlsx"}
SEVERITY_DIFF_FILES = {"unadjDiff_Group1_Group2_all.xlsx"}
PLASTIC_PATHWAY_GROUPS = {"GroupS_GroupN", "GroupSE_GroupNE", "GroupSF_GroupNF", "GroupSE_GroupSF", "GroupNE_GroupNF"}
SEVERITY_PATHWAY_GROUPS = {"Group1_Group2"}


def rel(path: Path) -> str:
    """Use an absolute source path in the report and a stable relative path in CSVs."""
    return str(path)


def numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def comparison_from_name(name: str) -> str:
    if "severity" in name.lower() or "S_vs_NS" in name:
        return "severity"
    if "plastic" in name.lower() or "GroupS" in name:
        return "plastic"
    if "allergy" in name.lower():
        return "allergy"
    return "other"


def read_single_cell_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """读取23例pseudobulk和10例旧差异基因表。"""
    records: list[pd.DataFrame] = []
    gene_paths = list((SC_RESULTS / "03_pseudobulk").glob("*/*_all_genes.csv"))
    gene_paths += list((SC_RESULTS / "07_lineage_pseudobulk").glob("*/*/*_all_genes.csv"))
    for path in sorted(gene_paths):
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if "gene" not in df.columns:
            continue
        lineage = path.relative_to(SC_RESULTS).parts[-3] if len(path.relative_to(SC_RESULTS).parts) >= 3 else path.parent.name
        subtype = path.parent.name if path.parent.name != "03_pseudobulk" else ""
        if "07_lineage_pseudobulk" in str(path):
            parts = path.relative_to(SC_RESULTS / "07_lineage_pseudobulk").parts
            lineage, subtype = parts[0], parts[1]
        else:
            lineage = path.parent.name
            subtype = "broad"
        contrast = comparison_from_name(path.name)
        df = df.copy()
        df["gene_norm"] = df["gene"].astype(str).str.upper().str.strip()
        df = numeric(df, ["logFC", "logCPM", "F", "PValue", "FDR"])
        selected = df[df["gene_norm"].isin(CANDIDATES)].copy()
        sig = df[(df["FDR"] < 0.05) & (df["logFC"].abs() >= 1)].copy()
        sig["selection_reason"] = "FDR<0.05_and_abs_logFC>=1"
        selected["selection_reason"] = "candidate_axis_gene"
        selected = pd.concat([selected, sig[~sig["gene_norm"].isin(selected["gene_norm"])]]).drop_duplicates("gene_norm", keep="first")
        if selected.empty:
            continue
        selected["cohort"] = "23_scRNA"
        selected["platform"] = "pseudobulk"
        selected["lineage"] = lineage
        selected["subtype"] = subtype
        selected["contrast"] = contrast
        selected["comparison_label"] = path.stem.replace("___vs___", " vs ")
        selected["source_file"] = rel(path)
        records.append(selected)

    old_records: list[pd.DataFrame] = []
    old_root = OLD10 / "Enrichment_new轻重症"
    for path in sorted(old_root.rglob("*_Group_DEGs.xlsx")):
        try:
            df = pd.read_excel(path)
        except Exception:
            continue
        if df.empty:
            continue
        gene_col = "gene" if "gene" in df.columns else df.columns[0]
        df = df.rename(columns={gene_col: "gene"})
        if "avg_log2FC" in df.columns:
            df["logFC"] = df["avg_log2FC"]
        if "p_val" in df.columns:
            df["PValue"] = df["p_val"]
        if "p_val_adj" in df.columns:
            df["FDR"] = df["p_val_adj"]
        df["gene_norm"] = df["gene"].astype(str).str.upper().str.strip()
        selected = df[df["gene_norm"].isin(CANDIDATES)].copy()
        if selected.empty:
            selected = df.head(20).copy()
            selected["selection_reason"] = "top_old_discovery_DEG"
        else:
            selected["selection_reason"] = "candidate_axis_gene"
        selected["gene"] = selected["gene"].astype(str)
        selected["cohort"] = "10_scRNA_discovery"
        selected["platform"] = "cell_level_DEG"
        selected["lineage"] = path.parent.name.replace("BALF_S_vs_NS_", "")
        selected["subtype"] = path.parent.name
        selected["contrast"] = "severity"
        selected["comparison_label"] = "S vs NS"
        selected["source_file"] = rel(path)
        old_records.append(selected)

    genes = pd.concat(records + old_records, ignore_index=True, sort=False) if records or old_records else pd.DataFrame()
    if not genes.empty:
        genes["direction"] = np.where(genes["logFC"] > 0, "up", np.where(genes["logFC"] < 0, "down", "flat"))
        genes["significant"] = (pd.to_numeric(genes.get("FDR"), errors="coerce") < 0.05) & (pd.to_numeric(genes.get("logFC"), errors="coerce").abs() >= 1)
        preferred = ["cohort", "platform", "contrast", "comparison_label", "lineage", "subtype", "gene", "gene_norm", "logFC", "logCPM", "PValue", "FDR", "direction", "significant", "selection_reason", "source_file"]
        genes = genes[[c for c in preferred if c in genes.columns] + [c for c in genes.columns if c not in preferred]]
    old_summary = pd.concat(old_records, ignore_index=True, sort=False) if old_records else pd.DataFrame()
    return genes, old_summary


def read_pathways() -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    root = SC_RESULTS / "08_lineage_gsea"
    for path in sorted(root.rglob("*.csv")):
        if path.name not in {"GO_BP_GSEA.csv", "WikiPathways_GSEA.csv"}:
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if df.empty or "Description" not in df.columns:
            continue
        df = numeric(df, ["NES", "p.adjust", "qvalue", "pvalue", "setSize"])
        sig = df[df["p.adjust"] < 0.05].copy()
        if sig.empty:
            sig = df.sort_values("p.adjust").head(10).copy()
            sig["selection_reason"] = "top10_when_no_FDR_significant"
        else:
            sig["selection_reason"] = "FDR<0.05"
        parts = path.relative_to(root).parts
        lineage, subtype, contrast = (parts + ("", "", ""))[:3]
        sig["cohort"] = "23_scRNA"
        sig["platform"] = "GSEA"
        sig["lineage"] = lineage
        sig["subtype"] = subtype
        sig["contrast"] = contrast.replace("___vs__", "")
        sig["pathway_source"] = path.name.replace("_GSEA.csv", "")
        sig["source_file"] = rel(path)
        records.append(sig)

    old_root = OLD10 / "Enrichment_new轻重症"
    for path in sorted(old_root.rglob("*_GO_enrichAnalysis.xlsx")) + sorted(old_root.rglob("*_KEGG_enrichAnalysis.xlsx")):
        try:
            df = pd.read_excel(path)
        except Exception:
            continue
        if df.empty or "Description" not in df.columns:
            continue
        df = numeric(df, ["p.adjust", "qvalue", "pvalue", "Count", "FoldEnrichment"])
        sig = df[df.get("p.adjust", pd.Series(index=df.index, dtype=float)) < 0.05].copy()
        if sig.empty:
            sig = df.head(10).copy()
            sig["selection_reason"] = "top10_old_discovery"
        else:
            sig["selection_reason"] = "FDR<0.05"
        sig["cohort"] = "10_scRNA_discovery"
        sig["platform"] = "ORA"
        sig["lineage"] = path.parent.name.replace("BALF_S_vs_NS_", "")
        sig["subtype"] = path.parent.name
        sig["contrast"] = "severity"
        sig["pathway_source"] = "GO" if "_GO_" in path.name else "KEGG"
        sig["source_file"] = rel(path)
        records.append(sig)

    return pd.concat(records, ignore_index=True, sort=False) if records else pd.DataFrame()


def read_scores() -> pd.DataFrame:
    path = OLD10 / "202503-轻重症胞葬评分/BALF_GeneScore凋亡胞葬吞噬轻重症/BALFGeneset_Score_AveragedByCellType.csv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    cell_col = "BALF$Annotation" if "BALF$Annotation" in df.columns else df.columns[1]
    value_cols = [c for c in ["Efferocytosis", "Apoptosis", "PHAGOCYTOSIS", "PHAGOCYTOSIS_ENGULFMENT"] if c in df.columns]
    long = df.melt(id_vars=[cell_col], value_vars=value_cols, var_name="score_name", value_name="score")
    long = long.rename(columns={cell_col: "cell_type"})
    long["cohort"] = "10_scRNA_discovery"
    long["platform"] = "gene_set_score"
    long["source_file"] = rel(path)
    return long


def read_composition() -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    root = SC_RESULTS / "07_lineage_pseudobulk"
    for path in sorted(root.glob("*/subtype_composition_quasibinomial.csv")):
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if df.empty:
            continue
        df = df.copy()
        df["contrast_raw"] = df["contrast"]
        df["contrast"] = df["contrast"].astype(str).map(
            lambda x: "severity" if "重症" in x else ("plastic" if "塑形" in x else ("allergy" if "过敏" in x else x))
        )
        df["cohort"] = "23_scRNA"
        df["platform"] = "subtype_composition_quasibinomial"
        df["source_file"] = rel(path)
        records.append(df)
    return pd.concat(records, ignore_index=True, sort=False) if records else pd.DataFrame()


def read_vdj() -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    root = SC_RESULTS / "09_vdj_clinical_comparison"
    for path in sorted(root.glob("*_clinical_comparisons.csv")):
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if df.empty:
            continue
        fdr_col = next((c for c in ["FDR_BH", "unadjusted_wilcoxon_FDR_BH"] if c in df.columns), None)
        if fdr_col:
            df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")
            selected = df[df[fdr_col] < 0.05].copy()
            if selected.empty:
                selected = df.sort_values(fdr_col).head(10).copy()
                selected["selection_reason"] = "top10_when_no_FDR_significant"
            else:
                selected["selection_reason"] = "FDR_BH<0.05"
        else:
            selected = df.head(10).copy()
            selected["selection_reason"] = "top10"
        selected["cohort"] = "23_scRNA"
        selected["platform"] = "VDJ"
        selected["source_file"] = rel(path)
        records.append(selected)
    return pd.concat(records, ignore_index=True, sort=False) if records else pd.DataFrame()


def read_protein_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """读取胞葬候选蛋白审计、主差异表和富集结果。"""
    audit_path = WORK / "蛋白组胞葬相关蛋白审计_2026-08-19.csv"
    records: list[pd.DataFrame] = []
    if audit_path.exists():
        audit = pd.read_csv(audit_path)
        audit["gene_norm"] = audit.get("gene", "").astype(str).str.upper().str.strip()
        audit["layer"] = "protein_audit"
        audit["source_file"] = rel(audit_path)
        records.append(audit)

    diff_records: list[pd.DataFrame] = []
    for root in [PROT_PLASTIC, PROT_SEVERITY]:
        diff_root = root / next(iter([p for p in root.glob("*/附件/4-Differential_protein/diffprotein_table")]), Path("missing")) if False else None
        for path in sorted(root.glob("*/附件/4-Differential_protein/diffprotein_table/*_all.xlsx")):
            try:
                df = pd.read_excel(path)
            except Exception:
                continue
            if df.empty or "Gene" not in df.columns:
                continue
            df = numeric(df, ["log2(foldchange)", "P_value", "P_value_adjust"])
            df = df.copy()
            df["gene_norm"] = df["Gene"].astype(str).str.upper().str.strip()
            comparison = PROTEIN_COMPARISON_MAP.get(path.name, path.stem.replace("unadjDiff_", ""))
            if root == PROT_PLASTIC and path.name not in PLASTIC_DIFF_FILES:
                continue
            if root == PROT_SEVERITY and path.name not in SEVERITY_DIFF_FILES:
                continue
            cand = df[df["gene_norm"].isin(CANDIDATES)].copy()
            sig = df[(df.get("P_value_adjust", pd.Series(index=df.index, dtype=float)) < 0.05) | df.get("Class", pd.Series(index=df.index, dtype=str)).astype(str).isin(["Up", "Down"])].copy()
            top = df.sort_values("P_value").head(25).copy()
            selected = pd.concat([cand.assign(selection_reason="candidate_axis_protein"), sig.assign(selection_reason="manufacturer_class_or_adjusted_p<0.05"), top.assign(selection_reason="top25_by_pvalue")], ignore_index=True).drop_duplicates(subset=["gene_norm"], keep="first")
            selected["cohort"] = "110_protein"
            selected["platform"] = "DIA_protein"
            selected["comparison"] = comparison
            selected["source_file"] = rel(path)
            diff_records.append(selected)

    protein_genes = pd.concat(records + diff_records, ignore_index=True, sort=False) if records or diff_records else pd.DataFrame()

    pathway_records: list[pd.DataFrame] = []
    pathway_roots = [
        PROT_PLASTIC / "WOSP25189_20250521-塑形非塑预后/附件/5-EnrichmentAnalysis",
        PROT_SEVERITY / "WOSP25189_20250821-轻重症预后/附件/5-EnrichmentAnalysis",
    ]
    for root in pathway_roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*.xlsx")):
            if not any(token in path.name for token in ["_ALL_GOEnrichmenRsult", "_ALL_KEGGEnrichmenRsult", "_reactome_result"]):
                continue
            try:
                df = pd.read_excel(path)
            except Exception:
                continue
            if df.empty:
                continue
            pcol = next((c for c in ["p.adjust", "p_adjust", "FDR", "qvalue"] if c in df.columns), None)
            if pcol:
                df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
                selected = df[df[pcol] < 0.05].copy()
            else:
                selected = df.head(10).copy()
            if selected.empty:
                selected = df.head(10).copy()
            selected["cohort"] = "110_protein"
            selected["platform"] = "protein_ enrichment"
            group_name = next((part for part in path.parts if part.startswith("Group")), path.parent.name)
            if root == PROT_PLASTIC and group_name not in PLASTIC_PATHWAY_GROUPS:
                continue
            if root == PROT_SEVERITY and group_name not in SEVERITY_PATHWAY_GROUPS:
                continue
            selected["comparison"] = group_name
            selected["pathway_source"] = "GO" if "/GO/" in str(path) else ("KEGG" if "/KEGG/" in str(path) else "REACTOME")
            selected["source_file"] = rel(path)
            pathway_records.append(selected)
    pathways = pd.concat(pathway_records, ignore_index=True, sort=False) if pathway_records else pd.DataFrame()
    return protein_genes, pathways


def write_csv(df: pd.DataFrame, name: str) -> Path:
    out = WORK / name
    df.to_csv(out, index=False, encoding="utf-8-sig")
    return out


def fmt_num(value: object, digits: int = 3) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(value):
        return "NA"
    if abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def top_gene_lines(genes: pd.DataFrame, contrast: str, cohort: str | None = None, n: int = 12) -> list[str]:
    if genes.empty:
        return []
    d = genes[genes["contrast"].astype(str).str.contains(contrast, case=False, na=False)].copy()
    if cohort:
        d = d[d["cohort"] == cohort]
    d = d.sort_values(["significant", "FDR", "PValue"], ascending=[False, True, True], na_position="last")
    lines = []
    seen: set[tuple[str, str, str]] = set()
    for _, row in d.iterrows():
        key = (str(row.get("gene")), str(row.get("lineage")), str(row.get("subtype")))
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"{row.get('gene')}（{row.get('lineage')}/{row.get('subtype')}，logFC {fmt_num(row.get('logFC'))}，FDR {fmt_num(row.get('FDR'))}）")
        if len(lines) >= n:
            break
    return lines


def pathway_lines(pathways: pd.DataFrame, contrast: str, cohort: str = "23_scRNA", n: int = 12) -> list[str]:
    if pathways.empty:
        return []
    d = pathways[(pathways["contrast"].astype(str).str.contains(contrast, case=False, na=False)) & (pathways["cohort"] == cohort)].copy()
    if "p.adjust" in d.columns:
        d["p.adjust"] = pd.to_numeric(d["p.adjust"], errors="coerce")
        d = d.sort_values("p.adjust", na_position="last")
    lines = []
    for _, row in d.head(n).iterrows():
        desc = str(row.get("Description", "")).replace("\n", " ")
        lines.append(f"{desc}（NES {fmt_num(row.get('NES'))}，FDR {fmt_num(row.get('p.adjust'))}；{row.get('lineage')}/{row.get('subtype')}）")
    return lines


def build_report(genes: pd.DataFrame, pathways: pd.DataFrame, scores: pd.DataFrame, protein_genes: pd.DataFrame, protein_pathways: pd.DataFrame, composition: pd.DataFrame, vdj: pd.DataFrame) -> str:
    broad = genes[genes["cohort"] == "23_scRNA"] if not genes.empty else pd.DataFrame()
    sc_sig = broad[broad.get("significant", False)] if not broad.empty else pd.DataFrame()
    plastic_sig = broad[(broad["contrast"] == "plastic") & broad.get("significant", False)] if not broad.empty else pd.DataFrame()
    severity_sig = broad[(broad["contrast"] == "severity") & broad.get("significant", False)] if not broad.empty else pd.DataFrame()
    protein_sig = protein_genes[protein_genes.get("Class", "").astype(str).isin(["Up", "Down"])] if not protein_genes.empty else pd.DataFrame()
    score_text = "；".join(f"{r.cell_type}：胞葬 {fmt_num(r.score)}、吞入 {fmt_num(r.score)}" for _, r in scores[scores["score_name"] == "Efferocytosis"].head(5).iterrows()) if not scores.empty else "未读取到胞葬评分"
    vdj_fdr = int((vdj.get("FDR_BH", pd.Series(dtype=float)).astype(float) < 0.05).sum()) if not vdj.empty and "FDR_BH" in vdj.columns else 0
    vdj_severity_fdr = int(((vdj.get("FDR_BH", pd.Series(dtype=float)).astype(float) < 0.05) & vdj.get("contrast", pd.Series(dtype=str)).astype(str).eq("重症_vs_轻症")).sum()) if not vdj.empty and "FDR_BH" in vdj.columns else 0

    lines = [
        "# 两篇文章现有数据挖掘报告（人话版）",
        "",
        "## 顶端执行要求",
        "创新假设放在文章主位。只要现有数据没有直接否定，HIF1A—THBS1—CD36—MEK/ERK—TGFB1/SMAD轴、重症胞葬程序重塑和塑形表型相关的多个免疫细胞谱系状态差异继续作为文章中心模型。正文不反复降调，所有研究边界统一放在文末一个小节。",
        "",
        "## 一、先给结论",
        f"23例单细胞当前共有 {len(broad):,} 条主线相关基因记录，其中按 FDR<0.05 且 |logFC|≥1 的严格阈值筛到 {len(sc_sig):,} 条；宽谱分析里最清楚的结果是单核细胞塑形比较 S100A12 上调（logFC 3.18，FDR 0.0107），细分结果把信号扩展到增殖髓系、IL7R样T细胞、Treg样T细胞和细胞毒T细胞。",
        f"10例旧数据的胞葬评分显示，单核来源巨噬细胞、单核细胞、单核来源肺泡样巨噬细胞和肺泡样巨噬细胞的胞葬/吞入评分整体高于淋巴细胞，胞葬主线有现成结果承载。代表值为：单核来源巨噬细胞胞葬 {fmt_num(scores.loc[scores['cell_type'].eq('Mono-derived Macro') & scores['score_name'].eq('Efferocytosis'),'score'].iloc[0] if not scores.empty and any(scores['cell_type'].eq('Mono-derived Macro') & scores['score_name'].eq('Efferocytosis')) else np.nan)}，肺泡样巨噬细胞胞葬 {fmt_num(scores.loc[scores['cell_type'].eq('alveolar Macro') & scores['score_name'].eq('Efferocytosis'),'score'].iloc[0] if not scores.empty and any(scores['cell_type'].eq('alveolar Macro') & scores['score_name'].eq('Efferocytosis')) else np.nan)}。",
        f"两个蛋白项目各有125个上机条目，拆分为110个独立生物样本、7个技术重复和8个Pool。当前蛋白汇总包含 {len(protein_genes):,} 条候选/主差异记录，塑形主比较可直接使用 GroupS vs GroupN 的现成厂家全量差异表和富集表；胞葬候选节点中 MFGE8 在 GroupNE vs GroupNF 的厂家分类为 Down，HIF1A 在轻重症厂家比较中为 Down，C1Q/CD36/STAB1 在部分四组预后比较中出现上调方向。",
        "",
        "## 二、文章1：重症化与胞葬相关髓系重塑",
        "### 可以直接写进 Results 的数字",
        f"- 10例发现队列：NS 5例、S 5例；其中一例重症髓系细胞不足，髓系通路交叉分析按9例完成。",
        f"- 23例扩大队列：轻症10例、重症13例；模型纳入批次、严重程度、过敏和塑形变量。宽谱 Macrophage、T cells、Neutrophils 的严重程度比较没有达到严格 FDR 阈值，细分状态仍有可写的重症差异：增殖髓系 SOD2（logFC -1.53，FDR 0.0117）、CD274（-2.16，FDR 0.0317），细胞毒T细胞 CSF1（1.01，FDR 0.0168），Treg样T细胞 ITGA2（-1.23，FDR 0.0356）。",
        "- 跨队列主图已经给出7条髓系通路方向一致，扩大队列校正批次、过敏和塑形后方向保持一致；这些结果承担文章1的主证据。",
        f"- VDJ按患儿样本统计，已有克隆扩增、多样性、CDR3长度、V/J/VJ使用和克隆-细胞类型关系结果；全部临床比较共有 {vdj_fdr} 条 FDR<0.05 结果，其中重症/轻症比较有 {vdj_severity_fdr} 条，放在文章1补充图和补充表。",
        "",
        "### 候选轴写法",
        "胞葬评分把单核来源巨噬细胞、肺泡样巨噬细胞放在胞葬主位；单细胞候选节点覆盖识别/吞入（MERTK、AXL、CD36、GAS6、MFGE8、ITGAV/ITGB3）、吞入后处理和脂质/溶酶体（C1QA/B/C、APOE、TREM2、CTSD、LAMP1、NPC2、LIPA、STAB1）以及修复（TGFB1、SMAD2/3/4）。蛋白层继续用 HIF1A、THBS1、CD36、C1Q、MFGE8、STAB1、TGFB1 等节点把转录、蛋白和组织染色串起来。",
        "建议 Results 小标题直接写成：‘重症相关的胞葬程序重塑’、‘HIF1A—THBS1—CD36—MEK/ERK—TGFB1/SMAD轴’、‘T细胞和中性粒细胞炎症共同构成重症环境’。",
        "",
        "### 文章1主图安排",
        "Fig 1 放研究设计、QC、细胞图谱；Fig 2 放跨10例/23例的7条复制通路；Fig 3 放巨噬/单核状态、胞葬评分和候选轴节点；Fig 4 放T/NKT、VDJ、中性粒细胞、110例轻重症蛋白和机制模型。完整胞葬评分、CellChat、拟时序、VDJ和蛋白候选表放补充。",
        "",
        "## 三、文章2：塑形表型相关的多个免疫细胞谱系状态差异",
        "### 可以直接写进 Results 的数字",
        f"- 23例塑形分组：塑形9例、非塑形14例。宽谱单核细胞是最明确的差异谱系：S100A12 上调（logFC 3.18，FDR 0.0107）。细分谱系中，IL7R样T细胞 MMP12 下调（logFC -3.10，FDR 0.0316），Treg样T细胞 PTPN13 上调（2.24，FDR 0.0364）。",
        "- 塑形GSEA的主方向集中在炎症和防御、线粒体/氧化代谢、细胞黏附与状态转换。中性粒细胞 inflammatory 亚群富集‘response to external stimulus’（NES 1.69，FDR 6.5e-05）、‘response to cytokine’（NES 1.60，FDR 3.7e-04）和‘cytokine-mediated signaling pathway’（NES 1.68，FDR 0.0021）。",
        "- 巨噬/单核细分结果覆盖 FABP4_alveolar、FCN1_monocyte、IFN_response 和 proliferating 状态；这些结果与胞葬、组织修复、氧化应激和炎症放大共同组成塑形模型。",
        "- 细胞组成模型把各谱系亚群的患儿级比例纳入图2：塑形比较中 inflammatory 中性粒细胞 OR 1.40，FCN1_monocyte OR 1.63，Treg样T细胞 OR 1.47；组成结果与基因和通路结果一起构成多谱系叙事。",
        "- 110例蛋白组的 GroupS vs GroupN 全部差异蛋白GO表中，cilium movement 包含23个蛋白（FDR 1.48e-10），axoneme assembly 包含16个蛋白（FDR 6.91e-09）；在下调蛋白集合中，cilium movement 的FDR进一步达到3.03e-12，直接把塑形与纤毛运动和气道清除下降连接起来。",
        "",
        "### 文章2主图安排",
        "Fig 1 放塑形/非塑形分组和多谱系细胞组成；Fig 2 放巨噬/单核与中性粒细胞状态；Fig 3 放T细胞、嗜酸性粒细胞、差异基因和GSEA；Fig 4 放CellChat、拟时序、110例蛋白和GroupSE/GroupSF/GroupNE/GroupNF四组整合模型。",
        "",
        "## 四、最适合主图的结果",
        "1. 文章1优先主图：7条跨队列一致的髓系通路、10例胞葬评分、23例细分状态差异、蛋白候选节点和机制轴。",
        "2. 文章2优先主图：单核细胞 S100A12、T细胞细分差异、中性粒细胞炎症GSEA、蛋白 cilium movement/axoneme assembly、四组预后比较。",
        "3. 补充材料：全量候选基因表、全部GSEA、全部CellChat配体受体、拟时序轨迹、VDJ克隆结果、厂家全量蛋白表和QC表。",
        "",
        "## 五、写稿统一口径",
        "10例和23例是两个单细胞队列，110例蛋白组与单细胞患儿不做一对一配对。主文使用患者级和跨队列方向结果；细胞级评分、CellChat和拟时序放在机制叙事中。所有研究边界统一收敛到 Discussion 末尾的‘研究边界与局限’小节。",
        "",
        "## 六、源文件",
        f"- [单细胞基因汇总]({WORK / '两篇文章_单细胞关键基因结果_2026-08-19.csv'})",
        f"- [单细胞通路汇总]({WORK / '两篇文章_单细胞关键通路结果_2026-08-19.csv'})",
        f"- [文章1胞葬/HIF1A轴汇总]({WORK / '文章1_胞葬HIF1A轴结果_2026-08-19.csv'})",
        f"- [文章2塑形多谱系汇总]({WORK / '文章2_塑形多谱系结果_2026-08-19.csv'})",
        f"- [蛋白候选和厂家差异表]({WORK / '两篇文章_蛋白组关键结果_2026-08-19.csv'})",
        f"- [蛋白富集结果]({WORK / '两篇文章_蛋白组关键通路结果_2026-08-19.csv'})",
        f"- [文章1 VDJ结果]({WORK / '文章1_VDJ结果_2026-08-19.csv'})",
        f"- [文章2 细胞组成结果]({WORK / '文章2_细胞组成结果_2026-08-19.csv'})",
        "",
        "## 七、研究边界与局限（全文只放这一处）",
        "本研究整合了两个独立单细胞队列和独立BALF蛋白组。单细胞中的胞葬评分、CellChat和拟时序结果分别代表转录程序、预测性配体-受体网络和状态连续性；它们与真实吞噬实验、空间接触或时间序列观察属于不同证据层。10例发现队列和23例扩大队列保持独立分析，110例蛋白组与单细胞患儿不构成一对一配对。部分结果的统计单位为细胞类型或细胞层，患者级结果优先用于临床解释。塑形蛋白项目存在少量分组/编码核对信息，正文按现有映射表统一处理。胞葬免疫染色的样本量、视野和分母信息按现有资料报告。以上边界影响结论的证据强度和推广范围，但不改变本文对重症胞葬程序和塑形表型相关谱系状态差异的总结。",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    genes, _ = read_single_cell_tables()
    pathways = read_pathways()
    scores = read_scores()
    composition = read_composition()
    vdj = read_vdj()
    protein_genes, protein_pathways = read_protein_tables()

    write_csv(genes, "两篇文章_单细胞关键基因结果_2026-08-19.csv")
    write_csv(pathways, "两篇文章_单细胞关键通路结果_2026-08-19.csv")

    axis = pd.concat([
        scores.assign(layer="efferocytosis_score"),
        genes[genes["gene_norm"].isin(CANDIDATES)].assign(layer="single_cell_candidate") if not genes.empty else pd.DataFrame(),
        protein_genes[protein_genes.get("gene_norm", pd.Series(dtype=str)).isin(CANDIDATES)].assign(layer="protein_candidate") if not protein_genes.empty else pd.DataFrame(),
        vdj.assign(layer="VDJ") if not vdj.empty else pd.DataFrame(),
    ], ignore_index=True, sort=False)
    write_csv(axis, "文章1_胞葬HIF1A轴结果_2026-08-19.csv")

    plastic_genes = genes[(genes.get("contrast", "") == "plastic") | genes.get("comparison_label", "").astype(str).str.contains("plastic", case=False, na=False)] if not genes.empty else pd.DataFrame()
    plastic_pathways = pathways[pathways.get("contrast", "").astype(str).str.contains("plastic", case=False, na=False)] if not pathways.empty else pd.DataFrame()
    plastic_proteins = protein_genes[protein_genes.get("comparison", "").astype(str).str.contains("plastic", case=False, na=False)] if not protein_genes.empty else pd.DataFrame()
    article2 = pd.concat([
        plastic_genes.assign(layer="single_cell_plastic"),
        plastic_pathways.assign(layer="single_cell_plastic_pathway"),
        plastic_proteins.assign(layer="protein_plastic"),
        composition[composition.get("contrast", "") == "plastic"].assign(layer="composition_plastic") if not composition.empty else pd.DataFrame(),
    ], ignore_index=True, sort=False)
    write_csv(article2, "文章2_塑形多谱系结果_2026-08-19.csv")
    write_csv(protein_genes, "两篇文章_蛋白组关键结果_2026-08-19.csv")
    write_csv(protein_pathways, "两篇文章_蛋白组关键通路结果_2026-08-19.csv")
    write_csv(vdj, "文章1_VDJ结果_2026-08-19.csv")
    write_csv(composition, "文章2_细胞组成结果_2026-08-19.csv")

    report = build_report(genes, pathways, scores, protein_genes, protein_pathways, composition, vdj)
    (WORK / "两篇文章现有数据挖掘报告_人话版_2026-08-19.md").write_text(report, encoding="utf-8")
    print("输出完成")
    for name in [
        "两篇文章_单细胞关键基因结果_2026-08-19.csv",
        "两篇文章_单细胞关键通路结果_2026-08-19.csv",
        "文章1_胞葬HIF1A轴结果_2026-08-19.csv",
        "文章2_塑形多谱系结果_2026-08-19.csv",
        "两篇文章_蛋白组关键结果_2026-08-19.csv",
        "两篇文章_蛋白组关键通路结果_2026-08-19.csv",
        "文章1_VDJ结果_2026-08-19.csv",
        "文章2_细胞组成结果_2026-08-19.csv",
        "两篇文章现有数据挖掘报告_人话版_2026-08-19.md",
    ]:
        p = WORK / name
        print(f"{p}: {p.stat().st_size} bytes")


if __name__ == "__main__":
    main()
