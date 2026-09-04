# Placental ACE2–GLUT analysis in gestational diabetes mellitus

This repository contains reproducible analysis code associated with the study:

**Placental ACE2 and glucose transporter expression in gestational diabetes mellitus: an exploratory observational pilot study**

## Scope

The project combines:

- placental RNA-seq analysis of **GSE255075**
- placental RNA-seq analysis of **GSE249311**
- targeted renin–angiotensin system (RAS) and SLC2/GLUT analyses
- pathway enrichment and pathway visualization
- cross-dataset RAS–GLUT concordance analyses
- exploratory ACE2–GLUT protein correlations in an institutional placental cohort

The institutional cohort is observational and small. Clinical subgroup analyses are therefore treated as exploratory and are not intended to establish causal effects of insulin or aspirin.

## Repository structure

```text
ACE2-GLUT/
├── ACE2-GLUT.Rproj
├── README.md
├── .gitignore
├── R/
│   ├── 00_project_setup.R
│   ├── 01_GSE255075_differential_expression.R
│   ├── 02_GSE255075_RAS_GLUT_analysis.R
│   ├── 03_GSE255075_pathway_analysis.R
│   ├── 04_GSE249311_differential_expression.R
│   ├── 05_GSE249311_RAS_GLUT_analysis.R
│   ├── 06_cross_dataset_concordance.R
│   ├── 07_clinical_ACE2_GLUT_correlations.R
│   └── archive/
├── data/
│   ├── raw/public/
│   ├── private/clinical/
│   └── processed/
├── outputs/
├── figures/
└── docs/
```

## Run order

Open `ACE2-GLUT.Rproj` in RStudio and run scripts from the project root:

1. `R/01_GSE255075_differential_expression.R`
2. `R/02_GSE255075_RAS_GLUT_analysis.R`
3. `R/03_GSE255075_pathway_analysis.R`
4. `R/04_GSE249311_differential_expression.R`
5. `R/05_GSE249311_RAS_GLUT_analysis.R`
6. `R/06_cross_dataset_concordance.R`
7. `R/07_clinical_ACE2_GLUT_correlations.R`

Each analysis script sources `R/00_project_setup.R`.

## Required input files

Place public inputs at:

```text
data/raw/public/GSE255075/GSE255075_2.csv
data/raw/public/GSE249311/GSE249311_count_matrix.xlsx
```

Place the institutional protein spreadsheet locally at:

```text
data/private/clinical/GDM_aspirin_protein_data.xlsx
```

`data/private/` is excluded from Git and must never be pushed to the public repository.

## Public datasets

- GSE255075
- GSE249311

Public RNA-seq data should be obtained from NCBI Gene Expression Omnibus or from the source files used in the validated historical analysis.

## Clinical data

Participant-level clinical and placental protein data are not included in this public repository because they are controlled human-subject data. Access is subject to institutional and ethical requirements.

## Statistical conventions

### RNA-seq differential expression

- DESeq2
- design: `~ condition`
- low-expression filter: total count >= 10 across samples
- Benjamini–Hochberg multiple-testing correction
- VST values are used for visualization/correlation, not for DE inference

### GSE249311 labels

The source-dataset categories are represented as:

- `Control`
- `GDMA1` — lifestyle-managed GDM
- `GDMA2` — medication-requiring GDM
- `T2DM` — pregestational type 2 diabetes

`GDMA2` is **not** treated as synonymous with insulin treatment.

### Correlations

RAS–SLC2 correlation heatmaps retain nominal Pearson P-value symbols for continuity with the manuscript analysis, while BH-adjusted P values are also exported in the underlying result tables.

### Institutional placental cohort

The biological unit is **one placenta per pregnancy**. Technical replicates are not treated as independent biological observations.

## Package requirements

Core packages include:

- here
- readxl
- DESeq2
- apeglm
- EnhancedVolcano
- gage
- pathview
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- pheatmap
- ComplexHeatmap
- circlize
- ggplot2
- dplyr
- tidyr
- tibble
- stringr
- openxlsx

The scripts stop with an informative message if required packages are missing.

## Reproducibility notes

Generated intermediate objects are written to `data/processed/`, and numerical results are written to `outputs/`. These directories are ignored by Git because they are reproducible from the source data and scripts.

Final manuscript figures may be copied intentionally into:

```text
figures/main/
figures/supplementary/
```

Historical outputs should be retained locally in `figures/archive/` or `outputs/archive/` until the modular workflow has been validated against the legacy pipeline.

## Recommended migration rule

**Reorganize first → reproduce historical results → validate → only then delete legacy duplicates.**
