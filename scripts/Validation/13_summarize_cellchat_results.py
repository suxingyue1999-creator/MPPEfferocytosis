import os
#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path('inputs/workspace/工作文件/23例单细胞重新分析/results/06_cellchat_by_sample')
OUT = ROOT / 'cellchat_nominal_top_candidates.csv'
REPORT = Path('inputs/workspace/工作文件/23例CellChat结果解读_人话版_2026-08-20.md')

comparisons = {
    'severity': ('pathway_severity_comparison.csv', 'Severe', 'Mild'),
    'plasticity': ('pathway_plasticity_comparison.csv', 'NonPlastic', 'Plastic'),
    'allergy': ('pathway_allergy_comparison.csv', 'Allergy', 'NonAllergy'),
}

rows = []
for comparison, (filename, group1, group2) in comparisons.items():
    df = pd.read_csv(ROOT / filename)
    df['comparison'] = comparison
    df['group1'] = group1
    df['group2'] = group2
    df['ratio_interpretation'] = f'{group2}/{group1}'
    df = df.rename(columns={'log2_ratio': 'log2_group2_over_group1'})
    rows.append(df.sort_values(['p_value', 'pathway_name']).head(20))

out = pd.concat(rows, ignore_index=True)
out = out[['comparison', 'group1', 'group2', 'ratio_interpretation', 'pathway_name',
           'n_group1', 'n_group2', 'mean_group1', 'mean_group2',
           'log2_group2_over_group1', 'p_value', 'fdr']]
out.to_csv(OUT, index=False, encoding='utf-8-sig')

audit = pd.read_csv(ROOT / 'run_audit.csv')
audit_map = dict(zip(audit['item'], audit['value']))

def top_lines(comparison):
    d = out[out['comparison'] == comparison].head(8)
    return [
        f"- `{r.pathway_name}`：{r.group2}相对{r.group1} log2比值 {r.log2_group2_over_group1:.2f}，名义P={r.p_value:.4g}，FDR={r.fdr:.4g}。"
        for r in d.itertuples()
    ]

report = [
    '# 23例正式样本级CellChat结果解读（人话版）', '',
    '## 先说结论', '',
    f"23例当前RNA对象已按患儿逐例重算CellChat：{int(audit_map['samples'])}例、{int(audit_map['cells']):,}个细胞、{int(audit_map['sample_objects'])}个CellChat对象，汇总得到{int(audit_map['communication_records']):,}条配体-受体通信记录。每个患儿单独计算，再把患儿作为统计单位比较临床组。", '',
    '这批结果可以用来画“网络候选机制图”和支持文章叙事；它不适合把某一条通路直接写成已经被临床组统计学确认的结论，因为样本级比较在当前23例中没有通路达到BH-FDR<0.05。', '',
    '## 三种比较的候选通路', '',
    '### 重症相对轻症', '',
    '重症组为group1，轻症组为group2；表中的比值是“轻症/重症”。',
    *top_lines('severity'), '',
    '其中PLAU和SPP1是当前重症比较中最靠前的候选通路，方向为轻症高于重症；CD80、CXCL和COLLAGEN也有同方向信号。', '',
    '### 塑形相对非塑形', '',
    '非塑形为group1，塑形为group2；表中的比值是“塑形/非塑形”。',
    *top_lines('plasticity'), '',
    '当前塑形比较最靠前的是CD6、VEGF和SPP1；巨噬细胞到B细胞的总通信量也呈塑形组升高方向，适合放入文章2的网络候选图。', '',
    '### 过敏相对非过敏', '',
    '过敏为group1，非过敏为group2；表中的比值是“非过敏/过敏”。',
    *top_lines('allergy'), '',
    'PTPRM和NEGR是名义P值靠前的候选通路。过敏比较作为背景分层使用，不承担两篇文章的主结论。', '',
    '## 文章里怎么用', '',
    '- 文章1：把PLAU、SPP1、CCL、CXCL和巨噬细胞-淋巴细胞/粒细胞互作放到重症胞葬环境的候选网络中；主证据仍由跨队列髓系通路、胞葬评分、候选轴基因和蛋白节点承担。',
    '- 文章2：把CD6、VEGF、SPP1、MHC-II、CCL以及巨噬细胞-单核细胞-中性粒细胞的通信变化放到塑形表型相关的候选网络图；蛋白层并列呈现纤毛/轴丝相关集合与ROS/过氧化氢相关集合的富集方向。',
    '- CellChat图注统一写“sample-level inferred communication”或“患儿级推断通信网络”，避免把预测通信写成直接实验验证。', '',
    '## 运行验收', '',
    '- 23/23个患儿对象成功保存；每例包含6类细胞身份、6×6网络矩阵和通路概率。',
    '- 个别样本的上皮细胞或中性粒细胞低于`min.cells=20`，CellChat自动移除这些细胞群相关通信；这些警告已写入`run.log`，不影响其余细胞群的对象保存。',
    '- 汇总输出：`communication_by_sample.csv`、`pathway_by_sample.csv`、`source_target_by_sample.csv`以及三种临床比较表。', '',
    f'- [CellChat候选通路表]({OUT})',
    f'- [CellChat对象级验收表]({ROOT / "cellchat_object_validation.csv"})',
    f'- [CellChat结果目录]({ROOT})',
]
REPORT.write_text('\n'.join(report) + '\n', encoding='utf-8')
print(OUT)
print(REPORT)
