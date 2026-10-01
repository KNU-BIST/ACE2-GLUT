# Locked manuscript-facing analysis results

This document records the validated numerical results used in the manuscript:

**Placental ACE2 and glucose transporter expression in gestational diabetes mellitus: an exploratory observational pilot study**

Its purpose is to provide a transparent audit trail between manuscript statements and the validated R analysis workflow.

No new statistical analyses are defined in this file.

---

## 1. GSE255075 differential expression

Analysis source: `R/01_GSE255075_differential_expression.R`

Dataset analyzed:

```text
Normal / Control: n = 3
GDM:              n = 3
Total:            n = 6
```

Low-expression filter:

```text
total count >= 10 across samples
```

Genes retained:

```text
16,839
```

Differential-expression summary:

```text
BH-adjusted P < 0.05: 3,124 genes
BH-adjusted P < 0.05 and |shrunken log2FC| > 1: 981 genes
```

Selected manuscript-facing results:

| Gene | Shrunken log2FC, GDM vs Normal | BH-adjusted P |
|---|---:|---:|
| ACE2 | +2.059 | 6.37e-05 |
| SLC2A4 | -1.411 | 6.77e-04 |
| AGTR1 | -0.622 | 0.00292 |
| REN | +1.502 | 0.0110 |
| SLC2A9 | +0.570 | 0.1057 |

ACE2, SLC2A4, AGTR1, and REN meet the manuscript differential-expression criterion of BH-adjusted P < 0.05. SLC2A9 does not.

The public dataset is used as transcriptomic context and does not reproduce the institutional clinical treatment strata.

---

## 2. GSE249311 reproducibility matrix

Primary analysis source: `R/04_GSE249311_differential_expression.R`

Targeted analysis source: `R/05_GSE249311_RAS_GLUT_analysis.R`

Validated historical matrix:

```text
Control:                       n = 11
GDM_diet / source GDMA1:       n = 5
GDM_insulin / source GDMA2:    n = 9
T2DM:                          n = 5
Total:                         n = 30
```

Source-category meaning:

```text
GDMA1 = lifestyle-managed GDM
GDMA2 = medication-requiring GDM
```

Manuscript-facing mapping:

```text
GDMA1 -> GDM_diet
GDMA2 -> GDM_insulin
```

`GDM_insulin` is retained as a manuscript/historical analysis label for the source medication-requiring GDM category. It does not establish verified insulin exposure for every public-dataset participant.

Historical low-expression filter:

```text
total count > 10
```

Feature rows retained:

```text
27,536
```

Selected manuscript-facing DESeq2 results:

| Contrast | Gene | log2FC | BH-adjusted P |
|---|---|---:|---:|
| GDM_diet vs Control | SLC2A1 | -0.629 | 0.000508 |
| GDM_diet vs Control | ACE2 | +1.126 | 0.0373 |
| GDM_insulin vs Control | SLC2A1 | -0.592 | 0.000249 |
| GDM_insulin vs Control | ACE2 | +0.514 | 0.348 |

SLC2A1 is lower in both manuscript-labeled GDM_diet and GDM_insulin groups relative to Control.

ACE2 is higher in GDM_diet relative to Control but does not meet the FDR significance criterion in GDM_insulin relative to Control.

These public data cannot establish an insulin-specific treatment effect.

---

## 3. GSE249311 ACE2–SLC2A9 correlation

Analysis source: `R/05_GSE249311_RAS_GLUT_analysis.R`

Validated two-sided Pearson correlation across all 30 samples:

```text
n = 30
Pearson r = 0.4054653
two-sided P = 0.0262212
95% CI = 0.05292761 to 0.66813554
BH-adjusted P across the tested RAS–SLC2 correlation family = 0.1398464
```

Manuscript rounding:

```text
r = 0.41
P = 0.026
95% CI = 0.053 to 0.668
BH-adjusted P = 0.140
```

The ACE2–SLC2A9 association is nominally significant using the unadjusted two-sided Pearson test but does not survive multiple-testing correction across the tested RAS–SLC2 correlation family.

The previously reported historical values of approximately r = 0.45 and P = 0.013 are not used in the revised manuscript because they were not reproducible from the validated workflow.

---

## 4. Institutional ACE2–GLUT protein correlations

Analysis source: `R/07_clinical_ACE2_GLUT_correlations.R`

Cohort:

```text
Control:       n = 3
GDM_diet:      n = 7
GDM_insulin:   n = 8
Total:         n = 18
```

All observations represent independent placentas from separate pregnancies.

Full-cohort correlations:

| Comparison | n | Pearson r | P |
|---|---:|---:|---:|
| ACE2 vs GLUT1 | 18 | 0.8108 | 4.48e-05 |
| ACE2 vs GLUT3 | 18 | 0.7048 | 0.001091 |
| ACE2 vs GLUT4 | 18 | -0.1903 | 0.4495 |
| ACE2 vs GLUT9 | 18 | 0.0359 | 0.8875 |

Manuscript rounding:

```text
ACE2–GLUT1: r = 0.81, P < 0.001
ACE2–GLUT3: r = 0.70, P = 0.001
ACE2–GLUT4: r = -0.19, P = 0.449
ACE2–GLUT9: r = 0.04, P = 0.887
```

The strongest institutional finding is the positive association of ACE2 protein abundance with GLUT1 and GLUT3 measured by ELISA.

These observational correlations do not establish regulatory direction, causality, or altered transplacental glucose transport.

---

## 5. Institutional sensitivity analyses

Analysis source: `R/07_clinical_ACE2_GLUT_correlations.R`

### GDM_diet-excluded sensitivity analysis

Control + GDM_insulin, n = 11:

| Comparison | n | Pearson r | P |
|---|---:|---:|---:|
| ACE2 vs GLUT1 | 11 | 0.6058 | 0.0482 |
| ACE2 vs GLUT3 | 11 | 0.4503 | 0.1646 |

The ACE2–GLUT1 correlation is attenuated but remains nominally significant after excluding GDM_diet.

The ACE2–GLUT3 correlation is attenuated and is not statistically significant in this sensitivity subset.

### Within-group exploratory correlations

GDM_diet:

```text
ACE2 vs GLUT1: n = 7, r = 0.9332, P = 0.00213
ACE2 vs GLUT3: n = 7, r = 0.8337, P = 0.0198
```

GDM_insulin:

```text
ACE2 vs GLUT1: n = 8, r = 0.6184, P = 0.1022
ACE2 vs GLUT3: n = 8, r = 0.4662, P = 0.2443
```

Control:

```text
n = 3
```

Quantitative within-Control correlations are not reported because the sample size is too small for a stable estimate.

All within-group analyses are exploratory.

---

## 6. Protein subgroup comparisons

The selected subgroup comparisons shown in the manuscript protein figures are two-sided unpaired Student t tests.

The complete set of 40 comparisons is provided in the manuscript supplementary table.

These analyses are unadjusted for multiplicity, include very small and unbalanced subgroup sizes, do not include covariate adjustment, and do not establish insulin or aspirin treatment effects.

Reported P values should therefore be described as nominal exploratory findings.

Aspirin exposure was not randomized and is susceptible to confounding by indication.

GDM treatment category is also associated with underlying glycemic severity and should not be interpreted as an independent pharmacological exposure.

---

## 7. Public-dataset interpretation boundary

The public transcriptomic datasets provide contextual evidence that ACE2, selected SLC2 genes, and RAS-related genes differ across diabetic pregnancy phenotypes.

They do not provide direct replication of the institutional ACE2–GLUT1 or ACE2–GLUT3 protein correlations.

Differences between GSE255075 and GSE249311 should be interpreted as evidence of biological and clinical heterogeneity rather than forced into a single common mechanistic model.

---

## 8. Manuscript figure-file mapping

The validated Script 07 outputs:

```text
figures/main/Figure4A_ACE2_GLUT1.png
figures/main/Figure4B_ACE2_GLUT3.png
```

These historical filenames correspond to the two ACE2–GLUT correlation panels used in the manuscript correlation figure.

The filenames are intentionally retained to avoid altering the validated script solely for manuscript figure renumbering.

---

## 9. Validated software environment

Primary validation environment:

```text
R version 4.5.3
```

Package and session details generated by individual scripts should be retained with the corresponding local analysis logs.

---

## 10. Reproducibility status

Scripts 01–07 were reproduced and validated before manuscript revision.

The manuscript-facing numerical results above represent the locked values to be used for submission unless a formally documented reanalysis is performed.
