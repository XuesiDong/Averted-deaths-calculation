# Six-cancer screening: averted deaths and probabilistic sensitivity analysis

This repository provides a concise R implementation for projecting **five-year index-cancer deaths potentially averted** under Chinese and multi-country screening stage distributions for six cancer sites: oesophageal, stomach, colorectal, liver, lung, and female breast cancer.

The script reproduces the **main-analysis point-estimate calculations** (including the original component-level rounding) and performs a **probabilistic sensitivity analysis (PSA)** with results for each cancer and **Overall**.

## Files

- `screening_psa.R` — analysis script with comments and model parameters.
- `stage.data.xlsx` — **required input; not included**. Supply your own study-specific stage-distribution workbook before running.
- `screening_psa_results.xlsx` — generated output containing `Results` and `Parameters` sheets.

## Requirements and run

R packages: `readxl`, `dplyr`, `tidyr`, `openxlsx`.

```r
install.packages(c("readxl", "dplyr", "tidyr", "openxlsx"))
source("screening_psa.R")
```

Run from a folder containing the script and Excel input, or edit `input_file` in the script. The input worksheet must have these columns:

| Column       | Required values                                                                                    |
| ------------ | -------------------------------------------------------------------------------------------------- |
| `Cancer`     | `Esophageal`/`Oesophageal`, `Stomach`, `Colorectal`/`Colorectum`, `Liver`, `Lung`, `Female breast` |
| `Group`      | `Non-screening`, `Chinese post-screening`, `Multi-country post-screening`                          |
| `Stage`      | `I`, `II`, `III`, `IV`                                                                             |
| `Percentage` | Nonnegative stage values (proportions or percentages; each cancer/group is normalised to sum to 1) |

Provide one row per cancer × group × stage, with no duplicates. The input workbook and underlying cohort/meta-analysis data are not supplied here.

## Method in brief

1. Read stage distributions for routine clinical diagnosis and two post-screening scenarios; normalise distributions for PSA while retaining the input values for the original main-analysis point estimates.
2. Calculate expected five-year index-cancer deaths as **incident cancer cases × the sum of stage proportions weighted by stage-specific death risks**.
3. Compute deaths potentially averted under Chinese and multi-country screening distributions, with and without observed adherence. Nonadherent individuals are assigned the routine-diagnosis stage distribution.
4. For PSA, sample each cancer's four stage-specific death risks on the **logit scale**, using their point estimates and standard errors. Keep cases, stage distributions and adherence fixed. Use the same sampled risks across comparison scenarios within a simulation.
5. Repeat **1,000 simulations** (`set.seed(123)`). Report the 2.5th and 97.5th percentiles as **95% uncertainty intervals (UI)**. For Overall, sum cancer outcomes **within each simulation**, then take percentiles of those totals.

The mortality risks, standard errors, incidence counts and adherence values are defined in `params` in the script and should be checked against the underlying study before reuse.

## Output columns

`Cancer`, `Scenario`, `Cases`, `Baseline_Deaths`, `Averted_Deaths`, `Death_Reduction_Pct`, `PSA_Mean`, `UI_Lower`, `UI_Upper`. The four scenarios are `China`, `Multi`, `China_Adherence`, and `Multi_Adherence`.

## Data and reproducibility

The analysis is deterministic given identical inputs, installed packages, simulation count and random seed. Publicly releasing the script **does not imply permission to distribute** the original input dataset; share only data you are authorised to release.
