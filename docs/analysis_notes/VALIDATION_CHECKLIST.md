# Modular-workflow validation checklist

Do not delete historical files until this checklist is complete.

## GSE255075

- [ ] Phase 1 reads the intended six samples.
- [ ] DESeq2 retained-gene count matches the historical pipeline.
- [ ] Number of genes with adjusted P < 0.05 matches the historical analysis.
- [ ] Direction and approximate magnitude of manuscript-critical fold changes match.
- [ ] PCA sample placement is consistent with the legacy output.
- [ ] Volcano plot labels the expected RAS/GLUT genes.
- [ ] Phase 2 RAS/GLUT heatmap contains the expected genes.
- [ ] Correlation table reproduces historical Pearson r values.
- [ ] Phase 3 GSEA/KEGG results reproduce manuscript-used pathways.

## GSE249311

- [ ] Exactly 30 samples are loaded.
- [ ] Group counts are Control 11 / GDMA1 5 / GDMA2 9 / T2DM 5.
- [ ] No sample is dropped during metadata joins.
- [ ] SLC2 boxplots contain all 30 samples.
- [ ] DESeq2 contrast tables use `padj` for multiplicity-controlled inference.
- [ ] ACE2–SLC2A9 Pearson r and nominal P reproduce the historical result.
- [ ] Correlation heatmap uses nominal symbols and CSV exports BH-adjusted P values.

## Institutional protein cohort

- [ ] Exactly 18 placentas are loaded.
- [ ] Group counts are Control 3 / GDM-diet 7 / GDM-insulin 8.
- [ ] One row corresponds to one placenta/pregnancy.
- [ ] ACE2–GLUT1 full-cohort r reproduces the manuscript value.
- [ ] ACE2–GLUT3 full-cohort r reproduces the manuscript value.
- [ ] Sensitivity results reproduce the supplementary table.
- [ ] Figure 4 visually matches the approved manuscript figure.

## After validation

- [ ] Copy manuscript-approved figures into `figures/main/` and `figures/supplementary/`.
- [ ] Archive or delete superseded loose PNGs.
- [ ] Archive old output directories.
- [ ] Commit only source code, documentation, and intentionally tracked final figures.
- [ ] Verify `git status` contains no patient-level spreadsheet.
