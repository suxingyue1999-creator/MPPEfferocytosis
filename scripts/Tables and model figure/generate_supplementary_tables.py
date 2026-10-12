from __future__ import annotations

import os

import csv
from pathlib import Path


ROOT = Path(os.environ.get("MPP_WORKSPACE", "inputs/workspace"))
WORK = ROOT / '工作文件'
OUT = WORK / '文章1' / '最终投稿包_2026-08-27' / '补充表'


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(name: str, rows: list[dict[str, object]], fields: list[str]) -> None:
    path = OUT / name
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, '') for field in fields})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # S1: fixed cohort identity. No patient-level clinical fields are inferred here.
    write_csv(
        'S1_队列身份.csv',
        [
            {
                'cohort_id': 'BD10_2026',
                'data_layer': 'single-cell RNA discovery',
                'role_in_article1': 'discovery cohort',
                'patients_total': 10,
                'mild_patients': 5,
                'severe_patients': 5,
                'cells_after_qc': 115594,
                'comparison_direction': 'severe vs mild',
                'pairing_status': 'not paired with the 23-patient scRNA cohort or the 110-patient protein cohort',
                'trajectory_boundary': 'Any ordering analysis is a myeloid-state arrangement, not a disease-time trajectory.',
                'source_file': str(ROOT / 'BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream/10samples_scrna_mex_v4/results/01_qc/10samples_qc_summary.csv'),
                'evidence_tier': 'primary discovery dataset',
            },
            {
                'cohort_id': '23_scRNA',
                'data_layer': 'single-cell RNA independent expansion',
                'role_in_article1': 'independent expansion cohort',
                'patients_total': 23,
                'mild_patients': 10,
                'severe_patients': 13,
                'cells_after_qc': 209299,
                'comparison_direction': 'severe vs mild',
                'pairing_status': 'not paired with the 10-patient scRNA cohort or the 110-patient protein cohort',
                'trajectory_boundary': 'No disease-time trajectory is assigned.',
                'source_file': str(WORK / '23例单细胞重新分析/results/01_qc/qc_summary_by_sample.csv'),
                'evidence_tier': 'independent expansion dataset',
            },
            {
                'cohort_id': 'BALF_proteomics_110',
                'data_layer': 'BALF proteinomics',
                'role_in_article1': 'independent protein-layer support',
                'patients_total': 110,
                'mild_patients': 61,
                'severe_patients': 49,
                'cells_after_qc': '',
                'comparison_direction': 'Group1 mild vs Group2 severe; effect is Group1 minus Group2',
                'pairing_status': 'independent cohort; not one-to-one paired with either scRNA cohort',
                'trajectory_boundary': 'No cell-state or disease-time trajectory is assigned.',
                'source_file': str(WORK / '蛋白组胞葬相关蛋白审计_2026-08-19.csv'),
                'evidence_tier': 'independent protein-layer dataset',
            },
        ],
        [
            'cohort_id', 'data_layer', 'role_in_article1', 'patients_total',
            'mild_patients', 'severe_patients', 'cells_after_qc',
            'comparison_direction', 'pairing_status', 'trajectory_boundary',
            'evidence_tier',
        ],
    )

    # S2: patient-level QC for the 10-patient discovery cohort.
    qc_path = ROOT / 'BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream/10samples_scrna_mex_v4/results/01_qc/10samples_qc_summary.csv'
    flags_path = ROOT / 'BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream/10samples_scrna_mex_v4/results/01_qc/qc_sample_flags.csv'
    metadata_path = ROOT / 'BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream/10samples_scrna_mex_v4/results/00_audit/10samples_metadata_locked.csv'
    qc_rows = {row['sample_id']: row for row in read_csv(qc_path)}
    flag_rows = {row['sample_id']: row for row in read_csv(flags_path)}
    metadata_rows = {row['sample_id']: row for row in read_csv(metadata_path)}
    s2 = []
    for sample_id in sorted(qc_rows, key=lambda x: int(x[1:])):
        q = qc_rows[sample_id]
        f = flag_rows[sample_id]
        m = metadata_rows[sample_id]
        s2.append({
            'cohort_id': m['cohort_id'],
            'sample_id': sample_id,
            'severity': m['severity'],
            'severity_en': m['severity_en'],
            'cells_before_qc': q['cells_input'],
            'cells_after_qc': q['cells_qc_keep'],
            'retention_rate': q['retention_rate'],
            'median_nCount_RNA': q['median_nCount_RNA'],
            'median_nFeature_RNA': q['median_nFeature_RNA'],
            'median_percent_mt': q['median_percent_mt'],
            'predicted_doublets': q['predicted_doublets'],
            'predicted_doublet_rate': q['predicted_doublet_rate'],
            'sample_outlier_flag': q['sample_outlier_flag'],
            'qc_exception_flag': q['qc_exception_flag'],
            'qc_outlier_flag': f['qc_outlier_flag'],
            'high_doublet_flag': f['high_doublet_flag'],
            'primary_analysis': f['primary_analysis'],
            'sensitivity_action': f['sensitivity_action'],
            'source_file_qc': str(qc_path),
            'source_file_flags': str(flags_path),
            'source_file_metadata': str(metadata_path),
            'evidence_tier': 'patient-level QC audit',
        })
    write_csv('S2_10例逐患儿QC.csv', s2, [
        'cohort_id', 'sample_id', 'severity', 'severity_en',
        'cells_before_qc', 'cells_after_qc', 'retention_rate',
        'median_nCount_RNA', 'median_nFeature_RNA', 'median_percent_mt',
        'predicted_doublets', 'predicted_doublet_rate', 'sample_outlier_flag',
        'qc_exception_flag', 'qc_outlier_flag', 'high_doublet_flag',
        'primary_analysis', 'sensitivity_action', 'evidence_tier',
    ])

    # S3: seven pathways with discovery and independent-expansion effects.
    robustness_path = WORK / '两篇文章_跨队列通路稳健性_2026-08-19.csv'
    crosscheck_path = WORK / '10例与23例重症通路交叉核对_2026-08-20.csv'
    meta_path = ROOT / 'BD_Rhapsody_setup/figure_work_package/cross_cohort_results/04_replicated_mechanism/replicated_pathway_score_meta_analysis.csv'
    camera_path = ROOT / 'BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream/10_23_cross_cohort_v1/results/03_meta_pathway/shared_GO_BP_cameraPR.csv.gz'
    robustness = read_csv(robustness_path)
    crosscheck = {row['go_id']: row for row in read_csv(crosscheck_path)}
    pathway_meta = {row['go_id']: row for row in read_csv(meta_path)}
    import gzip
    with gzip.open(camera_path, 'rt', encoding='utf-8', newline='') as handle:
        camera_rows = list(csv.DictReader(handle))
    camera = {
        (row['cohort'], row['harmonized_lineage'], row['go_id']): row
        for row in camera_rows
        if row['harmonized_lineage'] == 'Macrophage_Monocyte'
    }
    go_ids = [
        'GO:0002250', 'GO:0009615', 'GO:0031295', 'GO:0033209',
        'GO:0050862', 'GO:0050871', 'GO:0071222',
    ]
    s3 = []
    for go_id in go_ids:
        by_model = {(r['cohort'], r['model_key']): r for r in robustness if r['go_id'] == go_id}
        c = crosscheck[go_id]
        d = by_model[('discovery_10', 'discovery_primary')]
        v = by_model[('validation_23', 'validation_primary')]
        va = by_model[('validation_23', 'validation_adjusted')]
        pm = pathway_meta[go_id]
        d_camera = camera[('discovery_10', 'Macrophage_Monocyte', go_id)]
        e_camera = camera[('validation_23', 'Macrophage_Monocyte', go_id)]
        s3.append({
            'go_id': go_id,
            'term_en': d['term'],
            'term_cn': c['term_cn'],
            'comparison_direction': 'severe vs mild',
            'discovery10_effect': c['discovery10_formal_effect'],
            'discovery10_loo_min': d['min_loo_effect'],
            'discovery10_loo_max': d['max_loo_effect'],
            'discovery10_all_loo_negative': d['all_negative'],
            'discovery10_cameraPR_p': d_camera['PValue'],
            'discovery10_cameraPR_fdr': d_camera['FDR'],
            'expanded23_cameraPR_p': e_camera['PValue'],
            'expanded23_cameraPR_fdr': e_camera['FDR'],
            'expanded23_primary_effect': c['validation23_primary_effect'],
            'expanded23_primary_loo_min': v['min_loo_effect'],
            'expanded23_primary_loo_max': v['max_loo_effect'],
            'expanded23_adjusted_effect': c['validation23_adjusted_effect'],
            'expanded23_adjusted_loo_min': va['min_loo_effect'],
            'expanded23_adjusted_loo_max': va['max_loo_effect'],
            'expanded23_all_loo_negative': v['all_negative'],
            'go_bp_tests_per_cohort': 1188,
            'common_gene_universe': 15029,
            'discovery_pathway_model_n': 9,
            'discovery_pathway_model_mild_n': 5,
            'discovery_pathway_model_severe_n': 4,
            'expanded23_pathway_model_n': 23,
            'expanded23_pathway_model_mild_n': 10,
            'expanded23_pathway_model_severe_n': 13,
            'discovery_effect_se': pm['standard_error_discovery'],
            'expanded23_effect_se': pm['standard_error_validation'],
            'fixed_effect': pm['fixed_effect'],
            'fixed_se': pm['fixed_se'],
            'fixed_p_value': pm['fixed_p_value'],
            'fixed_ci95_low': str(float(pm['fixed_effect']) - 1.96 * float(pm['fixed_se'])),
            'fixed_ci95_high': str(float(pm['fixed_effect']) + 1.96 * float(pm['fixed_se'])),
            'random_effect': pm['random_effect'],
            'random_se': pm['random_se'],
            'random_p_value': pm['random_p_value'],
            'random_ci95_low': str(float(pm['random_effect']) - 1.96 * float(pm['random_se'])),
            'random_ci95_high': str(float(pm['random_effect']) + 1.96 * float(pm['random_se'])),
            'tau2': pm['tau2'],
            'i2_percent': pm['i2_percent'],
            'fixed_fdr': pm['fixed_fdr'],
            'random_fdr': pm['random_fdr'],
            'strict_cross_cohort_rule': 'Both cohorts must show the same direction and cameraPR BH-FDR < 0.05 within the respective 1,188-term GO Biological Process family.',
            'direction_interpretation': 'Negative effect means the pathway score or model effect is lower in severe than mild samples.',
            'source_file_robustness': str(robustness_path),
            'source_file_crosscheck': str(crosscheck_path),
            'evidence_tier': 'cross-cohort directionally replicated pathway result',
        })
    write_csv('S3_7条跨队列通路.csv', s3, [
        'go_id', 'term_en', 'term_cn', 'comparison_direction',
        'discovery10_effect', 'discovery10_loo_min', 'discovery10_loo_max',
        'discovery10_all_loo_negative', 'discovery10_cameraPR_p',
        'discovery10_cameraPR_fdr', 'expanded23_cameraPR_p',
        'expanded23_cameraPR_fdr', 'expanded23_primary_effect',
        'expanded23_primary_loo_min', 'expanded23_primary_loo_max',
        'expanded23_adjusted_effect', 'expanded23_adjusted_loo_min',
        'expanded23_adjusted_loo_max', 'expanded23_all_loo_negative',
            'go_bp_tests_per_cohort', 'common_gene_universe',
        'discovery_pathway_model_n', 'discovery_pathway_model_mild_n',
        'discovery_pathway_model_severe_n', 'expanded23_pathway_model_n',
        'expanded23_pathway_model_mild_n', 'expanded23_pathway_model_severe_n',
        'discovery_effect_se', 'expanded23_effect_se',
        'fixed_effect', 'fixed_se', 'fixed_p_value', 'fixed_ci95_low',
        'fixed_ci95_high', 'random_effect', 'random_se', 'random_p_value',
        'random_ci95_low', 'random_ci95_high', 'tau2', 'i2_percent',
        'fixed_fdr', 'random_fdr',
        'strict_cross_cohort_rule', 'direction_interpretation', 'evidence_tier',
    ])

    # S4: cell-type mean program scores, used only to localize the source of the program.
    score_path = WORK / '文章1_胞葬HIF1A轴结果_2026-08-19.csv'
    scores = [row for row in read_csv(score_path) if 'BALF_GeneScore' in row['source_file']]
    macrophage_like = {'Mono-derived Macro', 'Mono-derived alveolar Macro', 'alveolar Macro', 'Monocytes'}
    s4 = []
    for row in scores:
        s4.append({
            'cell_type': row['cell_type'],
            'score_name': row['score_name'],
            'score': row['score'],
            'cohort': row['cohort'],
            'platform': row['platform'],
            'layer': row['layer'],
            'comparison_direction': 'cell-type averaged program score; no severity contrast in this source table',
            'cell_type_group': 'macrophage/monocyte' if row['cell_type'] in macrophage_like else 'other immune cell type',
            'source_file': row['source_file'],
            'evidence_tier': 'program-source localization',
            'interpretation_boundary': 'This is a cell-type mean transcriptional program score. It is not a direct efferocytosis or phagocytic flux measurement and contains no patient-level score.',
        })
    write_csv('S4_胞葬评分来源.csv', s4, [
        'cell_type', 'score_name', 'score', 'cohort', 'platform', 'layer',
        'comparison_direction', 'cell_type_group', 'evidence_tier',
        'interpretation_boundary',
    ])

    # S5/S6: split the 10-patient formal candidate-axis table into main and supplementary nodes.
    node_path = WORK / '10例单细胞重症候选轴_2026-08-20.csv'
    nodes = read_csv(node_path)
    main_nodes = {'HIF1A', 'THBS1', 'CD36', 'MERTK', 'TGFB1', 'GAS6', 'SPP1'}
    modules = {
        'HIF1A': 'oxygen-response candidate upstream node',
        'THBS1': 'target-cell contact candidate node',
        'CD36': 'recognition/contact and lipid-uptake candidate node',
        'MERTK': 'recognition and engulfment candidate node',
        'TGFB1': 'repair and inflammatory-resolution candidate node',
        'GAS6': 'recognition and engulfment ligand candidate node',
        'SPP1': 'myeloid tissue-response and interaction candidate node',
        'CTSD': 'post-engulfment lysosomal-processing candidate node',
        'LAMP1': 'post-engulfment lysosomal candidate node',
        'APOE': 'lipid-processing candidate node',
        'C1QA': 'complement and efferocytosis-associated candidate node',
        'TREM2': 'lipid-processing and myeloid-state candidate node',
        'STAB1': 'clearance and repair candidate node',
    }
    discovery_gene_tests = {
        'Macrophage_Monocyte': 15280,
        'T_NK': 15892,
        'Neutrophil': 10715,
    }
    discovery_n = {
        'Macrophage_Monocyte': ('5', '4'),
        'T_NK': ('5', '5'),
        'Neutrophil': ('5', '4'),
    }
    def enrich_node(row: dict[str, str], figure_role: str) -> dict[str, str]:
        out = dict(row)
        lineage = row.get('cell_type', '')
        mild_n, severe_n = discovery_n.get(lineage, ('NA', 'NA'))
        out.update({
            'functional_module': modules.get(row['gene'], 'supplementary immune-environment candidate node'),
            'figure5_role': figure_role,
            'evidence_tier': '10-patient discovery candidate-axis result',
            'fdr_statement': 'FDR is retained exactly as reported; a candidate direction is not equivalent to FDR significance.',
            'source_file': str(node_path),
            'statistical_unit': 'patient-level pseudobulk; one aggregated count profile per patient and cell group',
            'model': 'limma-voom empirical Bayes model used in the 10-patient discovery pseudobulk result',
            'mild_n': mild_n,
            'severe_n': severe_n,
            'genes_tested_fdr_family': discovery_gene_tests.get(lineage, 'NA'),
            'fdr_family': 'all genes tested for the same cell group and severe-versus-mild contrast',
            'zero_missing_rule': 'Genes without estimable statistics are absent from the source result; no zero was recoded as missing.',
        })
        return out
    s5 = [enrich_node(row, 'main candidate model') for row in nodes if row['gene'] in main_nodes]
    s6 = [enrich_node(row, 'supplementary node') for row in nodes if row['gene'] not in main_nodes]
    node_fields = [
        'cohort', 'cell_type', 'contrast', 'gene', 'log2FC', 'PValue', 'FDR',
        'ci95_low', 'ci95_high', 'interpretation', 'article_use',
        'functional_module', 'figure5_role', 'evidence_tier', 'fdr_statement',
        'statistical_unit', 'model', 'mild_n', 'severe_n',
        'genes_tested_fdr_family', 'fdr_family', 'zero_missing_rule',
    ]
    write_csv('S5_HIF1A候选节点.csv', s5, node_fields)
    write_csv('S6_胞葬处理与旁路候选节点.csv', s6, node_fields)

    # S7: full 27-row article-1 VDJ result table, with the original source fields untouched.
    vdj_path = WORK / '文章1_VDJ结果_2026-08-19.csv'
    vdj_rows = read_csv(vdj_path)
    s7 = []
    cohort_group_sizes = {
        '重症_vs_轻症': ('10', '13'),
        '过敏_vs_不过敏': ('15', '8'),
        '塑形_vs_非塑形': ('14', '9'),
    }
    for row in vdj_rows:
        out = dict(row)
        if not out.get('reference_n') or not out.get('case_n'):
            group_sizes = cohort_group_sizes.get(out.get('contrast', ''))
            if group_sizes and out.get('samples_used') == '23':
                out['reference_n'], out['case_n'] = group_sizes
                out['comparison_status'] = 'group sizes derived from the locked 23-patient cohort labels'
                out['missing_reason'] = ''
            else:
                out['reference_n'] = out.get('reference_n') or 'NA'
                out['case_n'] = out.get('case_n') or 'NA'
                out['comparison_status'] = 'not derivable from the source table'
                out['missing_reason'] = 'The source table reports 15 available samples but does not provide reference-versus-case allocation.'
        else:
            out['comparison_status'] = 'reported in the source table'
            out['missing_reason'] = ''
        out['cohort_definition'] = 'independent expansion single-cell cohort: 10 mild patients and 13 severe patients'
        out['statistical_unit'] = 'patient-sample level comparison in the 23-patient independent expansion scRNA VDJ analysis'
        out['comparison_direction_note'] = 'For severity_vs_mild rows, positive log2 odds ratio indicates higher usage in severe samples; other clinical contrasts retain their original labels.'
        out['minimum_receptor_cells_rule'] = 'At least 20 receptor-bearing cells per patient for the corresponding receptor analysis.'
        out['minimum_available_patients_rule'] = 'At least 8 patients with available data to establish a clinical comparison model.'
        family_tests = {
            'V_gene': '202 features per clinical contrast after filtering',
            'J_gene': '88 features per clinical contrast after filtering',
            'VJ_pair': '1,140 features per clinical contrast after filtering',
            'mean_cdr3_aa_length': '21 receptor-locus metrics per clinical contrast',
            'expanded_fraction_ge5': '33 cell-type and receptor combinations per clinical contrast',
            'max_clone_size': '36 repertoire metrics per clinical contrast',
            'd50_clone_count': '36 repertoire metrics per clinical contrast',
            'inv_simpson': '36 repertoire metrics per clinical contrast',
            'simpson_diversity': '36 repertoire metrics per clinical contrast',
            'shannon': '36 repertoire metrics per clinical contrast',
            'top10_clone_fraction': '36 repertoire metrics per clinical contrast',
            'max_clone_fraction': '36 repertoire metrics per clinical contrast',
            'pielou_evenness': '36 repertoire metrics per clinical contrast',
            'clonality': '36 repertoire metrics per clinical contrast',
            'd50_clone_fraction': '36 repertoire metrics per clinical contrast',
        }
        if out.get('feature_type') in ('V_gene', 'J_gene', 'VJ_pair'):
            family_key = out['feature_type']
        elif out.get('metric'):
            family_key = out['metric']
        elif out.get('outcome'):
            family_key = out['outcome']
        else:
            family_key = ''
        out['fdr_family_tests'] = family_tests.get(family_key, 'not specified in source report')
        out['zero_missing_rule'] = 'Zero usage is retained as a valid fraction; unavailable features or summaries remain missing and are not recoded as zero.'
        out['fdr_family'] = 'BH correction was applied separately within each clinical contrast and result family.'
        out['evidence_tier'] = 'VDJ clinical comparison result'
        out['interpretation_boundary'] = 'Public CDR3 and receptor-signature resources are descriptive and do not establish Mycoplasma-specific antigen recognition.'
        s7.append(out)
    vdj_fields = [
        'feature_type', 'receptor', 'locus', 'feature', 'total_count',
        'samples_present', 'contrast', 'reference', 'case', 'samples_used',
        'reference_mean_fraction', 'case_mean_fraction',
        'raw_mean_fraction_difference', 'log2_odds_ratio',
        'adjusted_standard_error_log_odds', 'adjusted_p_value', 'FDR_BH',
        'selection_reason', 'cohort', 'platform', 'metric', 'reference_n',
        'case_n', 'reference_median', 'case_median', 'raw_median_difference',
        'transformation', 'adjusted_beta', 'adjusted_standard_error',
        'unadjusted_wilcoxon_p_value', 'celltype_label', 'outcome',
        'min_receptor_cells_per_sample', 'reference_median_fraction',
        'case_median_fraction', 'unadjusted_wilcoxon_FDR_BH',
        'comparison_status', 'missing_reason',
        'cohort_definition', 'statistical_unit', 'comparison_direction_note',
        'minimum_receptor_cells_rule', 'minimum_available_patients_rule',
        'fdr_family_tests', 'zero_missing_rule', 'fdr_family',
        'evidence_tier', 'interpretation_boundary',
    ]
    write_csv('S7_VDJ完整结果.csv', s7, vdj_fields)

    # S8: all 21 severity-comparison proteins and the direction convention.
    protein_path = WORK / '蛋白组胞葬相关蛋白审计_2026-08-19.csv'
    proteins = [row for row in read_csv(protein_path) if row['comparison'] == 'severity_Group1_vs_Group2']
    s8 = []
    for row in proteins:
        effect = float(row['log2FC_group1_vs_group2'])
        if effect < 0:
            clinical_direction = 'higher in severe BALF overall'
        elif effect > 0:
            clinical_direction = 'higher in mild BALF overall'
        else:
            clinical_direction = 'no directional difference'
        out = dict(row)
        out.update({
            'group_definition': 'Group1 mild, n=61; Group2 severe, n=49',
            'comparison_direction': 'Group1 minus Group2',
            'clinical_direction': clinical_direction,
            'statistical_unit': 'independent BALF proteinomics participant; not paired to either scRNA cohort',
            'fdr_method': 'BH FDR as reported in the source audit table',
            'evidence_tier': 'independent protein-layer directional support',
            'interpretation_boundary': 'Negative Group1−Group2 values indicate higher overall protein abundance in severe BALF. FDR values are not rewritten as significant.',
            'source_file': str(protein_path),
            'effective_n_group1_mild': 61,
            'effective_n_group2_severe': 49,
            'vendor_entries_total': 11005,
            'records_with_accession': 10983,
            'records_with_complete_gene_and_statistics': 10980,
            'fdr_test_universe': 10980,
            'statistical_model': 'Unpaired Welch t-test on processed log2 expression values, as reported by the supplier.',
            'normalization_and_batch_handling': 'Supplier-processed log2 expression; technical repeats and pooled samples used for QC and batch monitoring, not group testing; no batch correction was reported.',
            'missing_zero_rule': 'Supplier pipeline imputed missing values as zero; candidate extraction retained the supplied statistics without recoding additional values.',
        })
        s8.append(out)
    write_csv('S8_110例蛋白方向.csv', s8, [
        'comparison', 'gene', 'accession', 'log2FC_group1_vs_group2',
        'p_value', 'FDR', 'nominal_class', 'group_definition',
        'comparison_direction', 'clinical_direction', 'statistical_unit',
        'fdr_method', 'evidence_tier', 'interpretation_boundary',
        'effective_n_group1_mild', 'effective_n_group2_severe',
        'vendor_entries_total', 'records_with_accession',
        'records_with_complete_gene_and_statistics', 'fdr_test_universe',
        'statistical_model', 'normalization_and_batch_handling',
        'missing_zero_rule',
    ])

    # S9: patient-level CellChat object analysis for the 23-patient expansion cohort.
    # The object table is joined to the locked sample-level severity labels so the
    # supplementary file reports the actual statistical unit without exposing
    # project-only workflow fields.
    cellchat_dir = WORK / '23例单细胞重新分析/results/06_cellchat_by_sample'
    object_path = cellchat_dir / 'cellchat_object_validation.csv'
    sample_info_path = cellchat_dir / 'sample_info.csv'
    sample_info = {row['sample_id']: row for row in read_csv(sample_info_path)}
    s9 = []
    for row in read_csv(object_path):
        info = sample_info.get(row['sample_id'], {})
        out = {
            'sample_id': row['sample_id'],
            'severity_group': info.get('severity_group', ''),
            'class': row['class'],
            'cells': row['cells'],
            'cell_groups': row['cell_groups'],
            'interactions': row['interactions'],
            'network_rows': row['net_count_rows'],
            'network_columns': row['net_count_cols'],
            'pathways': row['pathways'],
            'has_probability_matrix': row['has_prob'],
            'statistical_unit': 'one CellChat object per patient sample',
            'evidence_tier': 'patient-level object analysis',
        }
        s9.append(out)
    write_csv('S9_CellChat对象级分析.csv', s9, [
        'sample_id', 'severity_group', 'class', 'cells', 'cell_groups',
        'interactions', 'network_rows', 'network_columns', 'pathways',
        'has_probability_matrix', 'statistical_unit', 'evidence_tier',
    ])

    def add_cellchat_comparison(row: dict[str, str], feature: str) -> dict[str, object]:
        # In the source analysis Severe=group1 and Mild=group2. The source
        # log2_ratio is therefore Mild/Severe and is retained without changing
        # its numeric value.
        source_missing = not row.get('n_group2')
        out = {
            feature: row[feature],
            'group1': 'Severe',
            'group2': 'Mild',
            'n_group1_severe': row['n_group1'],
            'n_group2_mild': row['n_group2'] or 'NA',
            'mean_group1_severe': row['mean_group1'],
            'mean_group2_mild': row['mean_group2'] or 'NA',
            'log2_ratio_mild_over_severe': row['log2_ratio'] or 'NA',
            'p_value_nominal': row['p_value'] or 'NA',
            'fdr_bh': row['fdr'] or 'NA',
            'comparison_status': 'not estimable from the source table' if source_missing else 'reported in the source table',
            'missing_reason': 'The source table contains no mild-group count or mild-group summary for this pathway; no group comparison was estimable.' if source_missing else '',
            'ratio_interpretation': 'Positive values indicate higher inferred communication in Mild than Severe; negative values indicate higher inferred communication in Severe.',
            'statistical_unit': 'patient-level CellChat object; Wilcoxon comparison across patient objects',
            'evidence_tier': 'patient-level candidate interaction network',
        }
        return out

    pathway_rows = read_csv(cellchat_dir / 'pathway_severity_comparison.csv')
    s10 = [add_cellchat_comparison(row, 'pathway_name') for row in pathway_rows]
    write_csv('S10_CellChat通路重症比较.csv', s10, [
        'pathway_name', 'group1', 'group2', 'n_group1_severe', 'n_group2_mild',
        'mean_group1_severe', 'mean_group2_mild',
        'log2_ratio_mild_over_severe', 'p_value_nominal', 'fdr_bh',
        'comparison_status', 'missing_reason',
        'ratio_interpretation', 'statistical_unit', 'evidence_tier',
    ])

    source_target_rows = read_csv(cellchat_dir / 'source_target_severity_comparison.csv')
    s11 = []
    for row in source_target_rows:
        out = add_cellchat_comparison(row, 'source')
        out.pop('source')
        out = {
            'source_cell_group': row['source'],
            'target_cell_group': row['target'],
            **out,
        }
        s11.append(out)
    write_csv('S11_CellChat来源靶细胞重症比较.csv', s11, [
        'source_cell_group', 'target_cell_group', 'group1', 'group2',
        'n_group1_severe', 'n_group2_mild', 'mean_group1_severe',
        'mean_group2_mild', 'log2_ratio_mild_over_severe', 'p_value_nominal',
        'fdr_bh', 'ratio_interpretation', 'statistical_unit', 'evidence_tier',
    ])

    readme = '''# 文章1补充表索引与统计口径

本目录仅保存文章1补充表，不修改正文、图文件或文章2结果。补充表保留可用于复核的结果字段，内部追溯路径集中记录在本索引的来源文件部分。

| 表 | 文件 | 内容 | 行数 | 统计单位 | 比较方向 | 校正方法 | 证据等级 |
|---|---|---|---:|---|---|---|---|
| S1 | S1_队列身份.csv | 10例发现、23例扩大、110例蛋白组的队列角色与身份 | 3 | 队列 | 轻重症身份固定 | 不适用 | 队列设计锁定信息 |
| S2 | S2_10例逐患儿QC.csv | 10例逐样本输入细胞、QC后细胞、保留率和质量标记 | 10 | 患儿样本 | 仅描述QC流程 | 预设QC标记；无多重检验 | 患儿级QC审计 |
| S3 | S3_7条跨队列通路.csv | 两个单细胞队列方向一致的7条通路 | 7 | 患儿级通路效应或通路分数 | severe vs mild；负值表示重症侧更低 | 每队列1,188条GO BP检验，15,029个共同基因；两队列cameraPR方向一致且BH-FDR<0.05 | 跨队列方向重复 |
| S4 | S4_胞葬评分来源.csv | 10例细胞类型平均Efferocytosis、Apoptosis、吞噬和吞入程序评分 | 56 | 细胞类型平均 | 来源表无轻重症比较 | 基因集评分；无患者级显著性校正 | 程序来源定位 |
| S5 | S5_HIF1A候选节点.csv | 10例髓系候选网络主节点 | 7 | Macrophage_Monocyte 细胞群的候选节点结果 | severe vs mild | BH FDR，原值保留 | 发现队列候选机制 |
| S6 | S6_胞葬处理与旁路候选节点.csv | 吞入后处理、脂质处理、修复及旁路免疫候选节点 | 13 | Macrophage_Monocyte、T_NK或Neutrophil细胞群的候选节点结果 | severe vs mild | BH FDR，原值保留 | 发现队列补充候选结果 |
| S7 | S7_VDJ完整结果.csv | 文章1使用的23例VDJ完整汇总结果 | 27 | 患儿样本级临床比较 | 轻症10例、重症13例；原始 contrast 保留 | 每例至少20个相应受体细胞，至少8例有可用数据；各比较家族分别BH校正 | VDJ临床比较 |
| S8 | S8_110例蛋白方向.csv | 110例BALF蛋白组21个胞葬相关候选蛋白 | 21 | 独立BALF样本级蛋白比较 | Group1轻症 − Group2重症；负值代表重症BALF整体更高 | BH FDR，原值保留 | 独立蛋白层方向支持 |
| S9 | S9_CellChat对象级分析.csv | 23例扩大单细胞队列的患儿级CellChat对象分析 | 23 | 患儿级CellChat对象 | 逐患儿对象，无组间比较 | 对象纳入与网络结构记录 | 患儿级对象分析 |
| S10 | S10_CellChat通路重症比较.csv | CellChat通路的重症与轻症比较 | 116 | 患儿级CellChat对象 | group1=Severe，group2=Mild；log2比值为Mild/Severe | BH FDR，原值保留 | 候选互作网络 |
| S11 | S11_CellChat来源靶细胞重症比较.csv | CellChat来源至靶细胞的重症与轻症比较 | 36 | 患儿级CellChat对象 | group1=Severe，group2=Mild；log2比值为Mild/Severe | BH FDR，原值保留 | 候选互作网络 |
| S12 | S12_两单细胞队列逐例身份去重审计.csv | 两个单细胞队列逐例身份、日期、去标识化病案号哈希和分析覆盖 | 33 | 患儿样本 | 队列身份固定 | 去重审计；不适用 | 队列独立性审计 |
| S13 | S13_逐患儿QC与分析覆盖.csv | 两个单细胞队列QC前后细胞数及各分析实际患儿分母 | 7 | 患儿或队列分析覆盖 | 仅描述覆盖 | 不适用 | 分母审计 |
| S14 | S14_统计模型缺失规则与FDR全集审计.csv | 各分析家族统计单位、模型、协变量、有效n、缺失和零值规则及BH测试全集 | 137 | 分析家族 | 依源表定义 | 各家族预设BH校正 | 统计口径审计 |

## 来源文件

- S1：`inputs/workspace/工作文件/BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream/10samples_scrna_mex_v4/results/01_qc/10samples_qc_summary.csv`；`inputs/workspace/工作文件/23例单细胞重新分析/results/01_qc/qc_summary_by_sample.csv`；`inputs/workspace/工作文件/蛋白组胞葬相关蛋白审计_2026-08-19.csv`
- S2：10例 `10samples_qc_summary.csv`、`qc_sample_flags.csv` 和 `10samples_metadata_locked.csv`
- S3：`inputs/workspace/工作文件/两篇文章_跨队列通路稳健性_2026-08-19.csv`；`inputs/workspace/工作文件/10例与23例重症通路交叉核对_2026-08-20.csv`
- S4：`inputs/workspace/工作文件/文章1_胞葬HIF1A轴结果_2026-08-19.csv` 中 `source_file` 含 `BALF_GeneScore` 的记录
- S5、S6：`inputs/workspace/工作文件/10例单细胞重症候选轴_2026-08-20.csv`
- S7：`inputs/workspace/工作文件/文章1_VDJ结果_2026-08-19.csv`
- S8：`inputs/workspace/工作文件/蛋白组胞葬相关蛋白审计_2026-08-19.csv` 中 `comparison=severity_Group1_vs_Group2` 的记录
- S9：`inputs/workspace/工作文件/23例单细胞重新分析/results/06_cellchat_by_sample/cellchat_object_validation.csv` 与同目录 `sample_info.csv`
- S10：`inputs/workspace/工作文件/23例单细胞重新分析/results/06_cellchat_by_sample/pathway_severity_comparison.csv`
- S11：`inputs/workspace/工作文件/23例单细胞重新分析/results/06_cellchat_by_sample/source_target_severity_comparison.csv`
- S12：`audit_article1_cohorts_statistics.py`读取10例、23例临床表、QC表、伪批量元数据、CellChat对象和V(D)J连接表后生成
- S13：`audit_article1_cohorts_statistics.py`汇总两个队列的QC表、通路评分、伪批量元数据、CellChat对象和V(D)J连接表后生成
- S14：`audit_article1_cohorts_statistics.py`汇总各分析源表的统计模型、有效分母、缺失规则、零值规则和BH测试全集后生成

## 解释边界

- 10例、23例和110例保持为三个独立数据身份，不能合并成33例，也不能把蛋白组写成单细胞的一对一配对验证。
- 10例拟时序或状态排序仅表示髓系状态排列，不能解释为真实病程时间线。
- S4 的胞葬评分是转录程序评分，用于定位程序来源，不能解释为真实吞噬通量。
- S5 和 S6 的候选节点保留原始 P 值和 BH FDR。FDR 未达到阈值时，表中保留方向和候选机制含义，不将其改写为显著结果。
- MAPK1 和 MAPK3 若在机制图中出现，只作为 MEK和ERK候选中继的转录节点。本次 S5 和 S6 仅收录正式候选网络源表的真实记录。
- VDJ 的公共 CDR3 和受体特征资源只能作描述性结果，不能解释为支原体特异性受体。
- CellChat 严格 BH-FDR<0.05 未检出通过通路时，只能作为候选互作网络，不能写成已证实的细胞通讯。
- S8 的 Group1−Group2 负值表示重症 BALF 整体蛋白相对更高；这属于跨层方向支持，不能改写为蛋白组 FDR 显著结果。
- S9 的每一行对应一个23例扩大队列患儿的CellChat对象。S10和S11中group1为重症、group2为轻症，`log2_ratio_mild_over_severe`为轻症相对重症的比值；正值表示轻症推断通信较高，负值表示重症推断通信较高。严格BH-FDR<0.05没有通路通过检验，因此两表只支持患儿级候选互作网络。
- S12至S14为投稿前审计表，补充队列独立性、逐患儿QC和分析覆盖、统计模型及FDR测试全集；审计表不新增样本或分析结果。
'''
    (OUT / 'README_补充表索引.md').write_text(readme, encoding='utf-8')


if __name__ == '__main__':
    main()
