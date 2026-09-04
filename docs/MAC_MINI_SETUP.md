# Mac mini migration and local-folder cleanup

## 1. Freeze the transferred Windows folder

Keep the transferred project as a read-only legacy backup, for example:

```text
~/Archive/GDM_legacy_windows_2026-09-04/
```

Do not reorganize or delete this copy until the modular workflow has been validated.

## 2. Clone the GitHub repository fresh

```bash
mkdir -p ~/Rprojects
cd ~/Rprojects
git clone https://github.com/KNU-BIST/ACE2-GLUT.git
cd ACE2-GLUT
```

Copy the clean project files from this refactoring package into the cloned repository.

## 3. Copy source data into the new structure

### GSE255075

Copy the validated source file to:

```text
data/raw/public/GSE255075/GSE255075_2.csv
```

### GSE249311

Copy:

```text
GSE249311_count_matrix.xlsx
```

to:

```text
data/raw/public/GSE249311/GSE249311_count_matrix.xlsx
```

### Institutional clinical data

Copy:

```text
GDM_aspirin_protein_data.xlsx
```

to:

```text
data/private/clinical/GDM_aspirin_protein_data.xlsx
```

Do not commit this file.

## 4. Move historical generated files

Move old result folders such as:

```text
GSE249311_outputs/
GSE249311_results_nominal/
GSE255075_outputs/
```

to a local archive outside Git, for example:

```text
outputs/archive_legacy/
```

Move old loose PNG figures to:

```text
figures/archive/
```

Do not delete them until the new workflow reproduces the manuscript-critical results.

## 5. Disable R workspace restoration

In RStudio:

**Tools → Global Options → General**

Set:

- Restore `.RData` into workspace at startup: **OFF**
- Save workspace to `.RData` on exit: **Never**

Then remove from the clean project:

```text
.RData
.Rhistory
.RDataTmp*
```

## 6. Open the project

Open:

```text
ACE2-GLUT.Rproj
```

Do not use `setwd()` in analysis scripts.

All paths are built with `here::here()` and should work on macOS, Windows, and Linux.

## 7. Validate before deleting legacy files

Recommended validation sequence:

1. Run Script 01.
2. Compare DESeq2 DEG counts and key genes with legacy outputs.
3. Run Script 02.
4. Compare RAS/GLUT heatmaps and correlation tables.
5. Run Script 03.
6. Compare pathway results used in the manuscript.
7. Run Scripts 04–05.
8. Confirm GSE249311 boxplots, DESeq2 contrasts, and ACE2–SLC2A9 correlation.
9. Run Script 06 and compare cross-dataset concordance.
10. Run Script 07 and verify Figure 4 correlations.

Only after these checks should duplicate historical outputs be deleted.

## 8. First Git commit after validation

```bash
git status
git add README.md .gitignore ACE2-GLUT.Rproj R/ data/README.md docs/
git status
git commit -m "Reorganize project into reproducible modular analysis workflow"
git pull --rebase origin main
git push origin main
```
