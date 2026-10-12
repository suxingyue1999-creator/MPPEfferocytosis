import os
from pathlib import Path
import csv


BASE = Path(os.environ.get("MPP_PUBLIC_ROOT", "inputs/public"))
OUT = Path('inputs/workspace/工作文件/文章1/最终投稿包_2026-08-27/中文_20260829/supplementary_tables')


def read_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    gse274 = read_csv(BASE / 'results/GSE274269_candidate_nodes.csv')
    gse179 = read_csv(BASE / 'results/GSE179051_candidate_nodes.csv')
    fields = [
        'dataset', 'comparison', 'sample_count', 'group_1', 'group_2', 'gene',
        'detected_n', 'detected_rate', 'mean_expression_group_1',
        'mean_expression_group_2', 'logFC_group_2_minus_group_1', 'P_value',
        'whole_dataset_FDR', 'candidate_family_tested_n', 'candidate_family_FDR',
    ]
    rows = []
    for row in gse274:
        rows.append({
            'dataset': row['dataset'],
            'comparison': row['comparison'],
            'sample_count': '174',
            'group_1': 'P1-1 (n=108)',
            'group_2': 'P1-2 (n=66)',
            'gene': row['gene'],
            'detected_n': row['detected_n'],
            'detected_rate': row['detected_rate'],
            'mean_expression_group_1': row['mean_logCPM_P1_1'],
            'mean_expression_group_2': row['mean_logCPM_P1_2'],
            'logFC_group_2_minus_group_1': row['logFC_P1_2_minus_P1_1'],
            'P_value': row['P_value'],
            'whole_dataset_FDR': row['FDR'],
            'candidate_family_tested_n': row['candidate_family_tested_n'],
            'candidate_family_FDR': row['FDR_candidate_tested'],
        })
    for row in gse179:
        rows.append({
            'dataset': row['dataset'],
            'comparison': row['comparison'],
            'sample_count': '6',
            'group_1': 'PBS (n=2)',
            'group_2': 'M. pneumoniae (n=2)',
            'gene': row['gene'],
            'detected_n': row['detected_n'],
            'detected_rate': row['detected_rate'],
            'mean_expression_group_1': row['mean_logCPM_PBS'],
            'mean_expression_group_2': row['mean_logCPM_M_pneumoniae'],
            'logFC_group_2_minus_group_1': row['logFC_M_pneumoniae_minus_PBS'],
            'P_value': row['P_value'],
            'whole_dataset_FDR': row['FDR'],
            'candidate_family_tested_n': row['candidate_family_tested_n'],
            'candidate_family_FDR': row['FDR_candidate_tested'],
        })
    write_csv(OUT / 'Table_S15_GEO_candidate_nodes.csv', fields, rows)

    prjna = read_csv(BASE / 'results/PRJNA792730/analysis/PRJNA792730_candidate_nodes.csv')
    prjna_fields = [
        'dataset', 'comparison', 'comparison_role', 'timepoint_group_1',
        'timepoint_group_2', 'gene', 'detected_n', 'detected_rate',
        'mean_logCPM_group_1', 'mean_logCPM_group_2', 'logFC_group_1_minus_group_2',
        'P_value', 'whole_transcriptome_FDR', 'candidate_family_tested_n',
        'candidate_family_FDR',
    ]
    prjna_rows = []
    labels = {
        '12h_vs_0h': ('12 h (n=3)', '0 h (n=3)'),
        '24h_vs_0h': ('24 h (n=3)', '0 h (n=3)'),
        '24h_vs_12h': ('24 h (n=3)', '12 h (n=3)'),
    }
    for row in prjna:
        g1, g2 = labels[row['comparison']]
        time_1, time_2 = {
            '12h_vs_0h': ('12h', '0h'),
            '24h_vs_0h': ('24h', '0h'),
            '24h_vs_12h': ('24h', '12h'),
        }[row['comparison']]
        prjna_rows.append({
            'dataset': row['dataset'],
            'comparison': row['comparison'],
            'comparison_role': row['comparison_role'],
            'timepoint_group_1': g1,
            'timepoint_group_2': g2,
            'gene': row['gene'],
            'detected_n': row['detected_n'],
            'detected_rate': row['detected_rate'],
            'mean_logCPM_group_1': row.get(f'mean_logCPM_{time_1}', 'NA'),
            'mean_logCPM_group_2': row.get(f'mean_logCPM_{time_2}', 'NA'),
            'logFC_group_1_minus_group_2': row['logFC'],
            'P_value': row['P_value'],
            'whole_transcriptome_FDR': row['FDR'],
            'candidate_family_tested_n': row['candidate_family_tested_n'],
            'candidate_family_FDR': row['FDR_candidate_tested'],
        })
    write_csv(OUT / 'Table_S16_PRJNA792730_candidate_nodes.csv', prjna_fields, prjna_rows)


if __name__ == '__main__':
    main()
