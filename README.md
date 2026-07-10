# LTBI_2023

Code and selected analysis-ready data for the study *Global, Regional and National Trends in Latent Tuberculosis Infection Prevalence From 1990 to 2023, With Projections to 2050*.

## Contents

- `script/`: R and Python scripts for data preparation, trend analysis, joinpoint regression, national analyses, forecasting, decomposition, diagnostics, and appendix construction.
- `Data/CleanGBD/`: project-specific, cleaned GBD 2023 LTBI prevalence estimates used by downstream analyses.
- `Data/iso_code.csv` and `Data/who_region.csv`: country and WHO-region lookup tables used by the national analyses.
- `outcome/`: tabular analytic results, including Table 1, supplementary tables, national AAPCs, and negative-binomial forecast outputs.

## Requirements

- R 4.4 or later
- Python 3.10 or later, with `pandas` and `numpy`
- Joinpoint Regression Program callable through the R package `nih.joinpoint`
- R packages: `Cairo`, `cowplot`, `dplyr`, `ggplot2`, `loo`, `nih.joinpoint`, `paletteer`, `patchwork`, `purrr`, `readr`, `rstan`, `segmented`, `sf`, `tidyr`, and `tidyverse`

The forecasting scripts also require a working Stan toolchain compatible with `rstan`.

## Reproduce the analyses

Run commands from the repository root. The following scripts can run directly from the included cleaned data, assuming the listed R packages and Joinpoint software are available:

```powershell
Rscript script/2_a_global_trend.R
Rscript script/2_b_age_trend.R
Rscript script/2_b_SDI_trend.R
Rscript script/2_a_table.R
Rscript script/3_a_national_aapc.R
Rscript script/3_c_national_joinpoint_sensitivity.R
```

For the full pipeline, first obtain the external data described in `DATA_AVAILABILITY.md`, then run:

```powershell
Rscript script/4_a_forecast.R
Rscript script/4_b_decomposition.R
Rscript script/4_c_forecast_model_diagnostics.R
Rscript script/4_d_plot_negbin_forecast.R
Rscript script/2_c_fine_age_sensitivity.R
Rscript script/3_b_national_figure.R
Rscript script/99_build_appendix.R
```

`script/0_gbd_export_metadata.R` and `script/1_data_prepare.py` are included for provenance and reprocessing from an independently downloaded GBD export. The raw GBD export is intentionally not included in this repository.