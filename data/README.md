# Data directory

## `raw/public/`

Publicly available RNA-seq inputs.

Expected files:

```text
raw/public/GSE255075/GSE255075_2.csv
raw/public/GSE249311/GSE249311_count_matrix.xlsx
```

## `private/clinical/`

Controlled participant-level institutional data.

Expected local file:

```text
private/clinical/GDM_aspirin_protein_data.xlsx
```

This directory is excluded by `.gitignore` and must not be committed.

## `processed/`

Generated `.rds` objects and processed metadata used to connect analysis phases.

This directory is reproducible and excluded from Git.
