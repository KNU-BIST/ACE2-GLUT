# Placental ACE2–GLUT analysis in gestational diabetes mellitus

This repository contains the validated reproducible analysis workflow associated with the study:

**Placental ACE2 and glucose transporter expression in gestational diabetes mellitus: an exploratory observational pilot study**

The project combines placental protein measurements from a small institutional cohort with reanalysis of two public placental RNA-sequencing datasets.

## Study scope

The repository contains code for:

- differential-expression analysis of **GSE255075**
- targeted RAS–SLC2/GLUT analysis of **GSE255075**
- exploratory pathway analysis of **GSE255075**
- differential-expression analysis of **GSE249311**
- targeted RAS–SLC2/GLUT analysis of **GSE249311**
- descriptive cross-dataset comparison of RAS–SLC2 correlations
- exploratory ACE2–GLUT protein correlations in an institutional placental cohort

The institutional cohort contains 18 independent placentas:

- normoglycemic control: n = 3
- GDM_diet: n = 7
- GDM_insulin: n = 8

The clinical cohort is observational and small. Treatment category and aspirin exposure were not randomized. Clinical subgroup analyses are therefore exploratory and are not intended to establish causal effects of insulin or aspirin.

## Analysis run order

Open `ACE2-GLUT.Rproj` in RStudio and run the scripts from the project root in numerical order:

```text
R/01_GSE255075_differential_expression.R
R/02_GSE255075_RAS_GLUT_analysis.R
R/03_GSE255075_pathway_analysis.R
R/04_GSE249311_differential_expression.R
R/05_GSE249311_RAS_GLUT_analysis.R
R/06_cross_dataset_concordance.R
R/07_clinical_ACE2_GLUT_correlations.R
```

Each analysis script sources `R/00_project_setup.R`.

The validated workflow was reproduced using **R 4.5.3**.

## Required input files

Public inputs are expected locally at:

```text
data/raw/public/GSE255075/GSE255075_2.csv
data/raw/public/GSE249311/GSE249311_count_matrix.xlsx
```

The institutional protein dataset is expected locally at:

```text
data/private/clinical/GDM_aspirin_protein_data.xlsx
```

Participant-level institutional data are controlled human-subject data and are not included in this public repository.

The complete `data/private/` directory is excluded from Git and must never be committed.

Public raw input files are also kept locally rather than version-controlled. See `data/README.md` for dataset provenance and analysis-specific notes.

## Public datasets

### GSE255075

The validated analysis uses 3 Control/Normal and 3 GDM placentas.

The low-expression filter retains genes with total count >= 10 across the six samples.

Differential expression is analyzed with DESeq2 using `design = ~ condition`.

Log2 fold changes used for the manuscript-facing selected-gene results are apeglm-shrunken estimates.

VST values used for downstream visualization and correlation are generated with `blind = FALSE`.

### GSE249311

The full published GSE249311 bulk RNA-sequencing study contains more samples than the historical count matrix used in this project.

The validated reproducibility workflow deliberately retains the same historical 30-sample matrix:

```text
Control:                       n = 11
GDM_diet / source GDMA1:       n = 5
GDM_insulin / source GDMA2:    n = 9
T2DM:                          n = 5
Total:                         n = 30
```

The original public-dataset categories are:

```text
Control
GDMA1 — lifestyle-managed GDM
GDMA2 — medication-requiring GDM
T2DM  — pregestational type 2 diabetes
```

For consistency with the manuscript and historical analysis workflow:

```text
GDMA1 -> GDM_diet
GDMA2 -> GDM_insulin
```

`GDM_insulin` is retained as a manuscript/historical analysis label for the source medication-requiring GDM category. It must **not** be interpreted as evidence that every GDMA2 participant received insulin.

The historical GSE249311 low-expression filter is total count > 10 across samples. Duplicate source feature identifiers are retained as separate rows using deterministic `make.unique()` identifiers to reproduce the historical workflow. VST values are generated with `blind = TRUE`. The historical DESeq2 independent-filtering setting `alpha = 0.1` is retained for reproducibility, while significance remains defined using Benjamini–Hochberg adjusted P < 0.05.

## Statistical conventions

Formal transcriptome-level differential-expression inference uses DESeq2 with Benjamini–Hochberg multiple-testing correction. VST-normalized values are used for visualization, exploratory summaries, heatmaps, and correlations, not as a substitute for count-based DE inference.

RAS–SLC2 correlation analyses use Pearson correlations on VST-normalized expression. Heatmap significance symbols retain nominal P values for continuity with the historical analysis, while BH-adjusted P values are exported separately.

The validated GSE249311 ACE2–SLC2A9 association is:

```text
Pearson r = 0.405
two-sided P = 0.026
95% CI = 0.053 to 0.668
BH-adjusted P = 0.140
```

The association is therefore nominal and does not survive correction across the tested RAS–SLC2 correlation family.

For the institutional cohort, the biological unit is one placenta per pregnancy. Technical replicates are summarized at the placenta level and are not treated as independent biological observations. Protein subgroup comparisons and ACE2–GLUT correlations are exploratory and observational and do not establish causality.

## Manuscript-facing validated results

The exact validated numerical results used in the manuscript are summarized in:

```text
docs/MANUSCRIPT_RESULTS.md
```

## Figures

Validated historical correlation-plot filenames currently include:

```text
figures/main/Figure4A_ACE2_GLUT1.png
figures/main/Figure4B_ACE2_GLUT3.png
```

These correspond to the ACE2–GLUT1 and ACE2–GLUT3 panels used in the manuscript correlation figure. The historical filenames are intentionally retained to avoid modifying the already validated Script 07 solely for manuscript renumbering.

## Generated outputs

Generated intermediate R objects are written to `data/processed/`. Generated numerical outputs are written to `outputs/`. Generated exploratory figures are written under `figures/exploratory/`.

These directories are excluded from Git because their contents can be regenerated from the required inputs and analysis scripts.

Selected validated manuscript-facing figures may be intentionally retained in `figures/main/` and `figures/supplementary/`.

## Data availability

The public RNA-sequencing datasets are available through the NCBI Gene Expression Omnibus under accessions GSE255075 and GSE249311.

The institutional participant-level dataset is not publicly distributed because of participant privacy and Biobank governance requirements. De-identified institutional data may be available subject to corresponding-author request and institutional approval, as described in the manuscript.

## Reproducibility boundary

The repository reproduces the analysis used for the manuscript.

The public transcriptomic datasets are used as contextual evidence and are not interpreted as mechanistic validation of the institutional ACE2–GLUT protein associations.

Likewise, subgroup differences associated with GDM treatment category or aspirin exposure are observational and should not be interpreted as independent pharmacological effects.

## Core R packages

Core packages include `here`, `readxl`, `DESeq2`, `apeglm`, `EnhancedVolcano`, `gage`, `pathview`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `pheatmap`, `ComplexHeatmap`, `circlize`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `stringr`, `AnnotationDbi`, `openxlsx`, and `SummarizedExperiment`.

## Citation

If using this analysis workflow, please cite the associated manuscript once publication details become available.

Repository: https://github.com/KNU-BIST/ACE2-GLUT
