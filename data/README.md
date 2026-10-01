# Data directory

This directory contains the expected locations for public RNA-sequencing inputs, controlled institutional data, and generated processed objects used by the validated ACE2–GLUT analysis workflow.

## `raw/public/`

Publicly available or historically retained RNA-sequencing inputs are stored locally under:

```text
raw/public/GSE255075/
raw/public/GSE249311/
```

Expected files:

```text
raw/public/GSE255075/GSE255075_2.csv
raw/public/GSE249311/GSE249311_count_matrix.xlsx
```

Public raw data files are intentionally not version-controlled in this repository.

### GSE255075

The validated workflow analyzes:

```text
Normal / Control: n = 3
GDM:              n = 3
```

The input file is `GSE255075_2.csv`.

Validated low-expression criterion:

```text
total count >= 10 across all six samples
```

After filtering, the analysis retains:

```text
16,839 genes
```

### GSE249311

The published GSE249311 bulk RNA-sequencing study contains more samples than the historical matrix used in this project.

The validated reproducibility workflow deliberately uses the same historical 30-sample matrix:

```text
Control:                       n = 11
GDM_diet / source GDMA1:       n = 5
GDM_insulin / source GDMA2:    n = 9
T2DM:                          n = 5
Total:                         n = 30
```

Source terminology:

```text
GDMA1 = lifestyle-managed GDM
GDMA2 = medication-requiring GDM
```

Manuscript-facing terminology:

```text
GDMA1 -> GDM_diet
GDMA2 -> GDM_insulin
```

`GDM_insulin` is retained as a manuscript/historical analysis label. The source GDMA2 category should not be interpreted as proof that every participant received insulin.

The historical input matrix is `GSE249311_count_matrix.xlsx`.

Validated low-expression criterion:

```text
total count > 10 across all 30 samples
```

After filtering, the historical workflow retains:

```text
27,536 feature rows
```

Duplicate source feature labels are intentionally retained as separate rows using `make.unique()` identifiers because this reproduces the historical analysis.

## `private/clinical/`

Controlled participant-level institutional data are expected locally at:

```text
private/clinical/GDM_aspirin_protein_data.xlsx
```

The institutional cohort contains:

```text
Control:       n = 3
GDM_diet:      n = 7
GDM_insulin:   n = 8
Total:         n = 18 placentas
```

The biological unit is one placenta per pregnancy.

Technical assay replicates are summarized at the placenta level and are not treated as independent biological observations.

Participant-level clinical and placental protein data are controlled human-subject data and are not distributed through this public repository.

The complete directory `data/private/` is excluded by `.gitignore` and must never be committed.

## `processed/`

Generated DESeq2, VST, metadata, and related `.rds` objects are written under `processed/`.

These files connect sequential analysis phases but can be regenerated from the required inputs and scripts, so the directory is excluded from Git.

## Privacy and reproducibility

No participant-level institutional clinical dataset should be copied into the public repository.

Public dataset analyses can be reproduced after placing the required source files in the paths described above and running the numbered R scripts from the repository root.
