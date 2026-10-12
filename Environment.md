# Recorded analysis environments

The versions below come from existing stage-specific session records, not from a
new analysis run. They should not be replaced with one global environment claim.

| Stage | Recorded environment | Main dependencies |
|---|---|---|
| Historical discovery/cross-cohort stages | R 4.1.2, Ubuntu 22.04.5; cross-cohort GO stage records limma 3.50.0, org.Hs.eg.db 3.14.0, AnnotationDbi 1.56.2, ggplot2 3.3.5 | Matrix, Seurat, edgeR, limma, SoupX, DropletUtils, digest, AnnotationDbi, org.Hs.eg.db; Python numpy, pandas, scipy, scrublet |
| Validation RNA | R 4.5.2, Seurat 5.4.0, SeuratObject 5.3.0, Matrix 1.7-4, ggplot2 4.0.1, future 1.69.0, patchwork 1.3.2 | The listed packages; lineage gene-set analysis also uses clusterProfiler and org.Hs.eg.db |
| Sample-specific interaction analysis | Validation R environment | CellChat, Seurat, data.table, future; CellChatDB.human supplied by the original CellChat package |
| GEO count-matrix analyses | R 4.5.2, edgeR 4.8.2, limma 3.66.0, org.Hs.eg.db 3.22.0, AnnotationDbi 1.72.0, data.table 1.18.2.1 | The listed packages |
| PRJNA792730 | Salmon 1.10.3, GENCODE v44 | R edgeR, data.table, jsonlite; figure code also uses ggplot2, patchwork and scales; Perl and gzip |
| Protein/result imports | Python | pandas, numpy, openpyxl; protein spreadsheet conversion uses LibreOffice `soffice` |
| Historical figures | R | data.table, Matrix, Seurat, ggplot2, patchwork, ragg, svglite, grid |

The PRJNA quantification shell script retains the original macOS `stat -f` and
`md5 -q` checks. On another operating system, equivalent file-size/checksum tools
are needed before running it; the scientific quantification options are unchanged.
R package and annotation-database versions matter for exact gene filtering and GO
gene-set membership. No dependency installation is performed by this release.

Source records:

- Discovery `results/01_qc`, `02_integrate`, `03_annotate` and `04_compare` records.
- Cross-cohort `results/03_meta_pathway/sessionInfo.txt`.
- Validation `results/02_rna_reanalysis/sessionInfo.txt` and its other stage records.
- Public-data `results/sessionInfo.txt` and PRJNA `analysis/sessionInfo.txt`.

These records are located in the authors' historical workspace; data-bearing
delivery archives and analysis objects were intentionally not copied into the
code repository. The preparation checks used locally available parsers only and
do not establish that all recorded dependency combinations are installed now.
