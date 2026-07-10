#####################################
## @Description: Build appendix from generated figures and tables
## @Date: 2026-06-02
#####################################

set_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
    setwd(dirname(dirname(script_path)))
  }
}

set_project_root()

library(dplyr)
library(readr)

source("./script/config.R")
source("./script/function.R")

appendix_path <- "../appendix.md"
lines_out <- character()

add <- function(...) {
  lines_out <<- c(lines_out, ...)
}

add_blank <- function() {
  add("")
}

add_page_break <- function() {
  add('<div style="page-break-after: always;"></div>')
  add_blank()
}

remove_trailing_page_break <- function() {
  while (length(lines_out) > 0 && tail(lines_out, 1) == "") {
    lines_out <<- head(lines_out, -1)
  }

  page_break <- '<div style="page-break-after: always;"></div>'
  if (length(lines_out) > 0 && tail(lines_out, 1) == page_break) {
    lines_out <<- head(lines_out, -1)
  }
}

add_table <- function(caption, table_file) {
  if (!file.exists(table_file)) {
    stop("Missing generated table: ", table_file)
  }
  add(paste0("**", caption, "**"))
  add_blank()
  add(readLines(table_file, warn = FALSE, encoding = "UTF-8"))
  add_blank()
  add_page_break()
}

add_figure <- function(image_file, image_md, alt, title, legend) {
  if (!file.exists(image_file)) {
    stop("Missing generated figure: ", image_file)
  }
  add(paste0("![", alt, "](", image_md, ")"))
  add_blank()
  add(paste0("**", title, "** ", legend))
  add_blank()
  add_page_break()
}

fig_file <- function(name) file.path(appendix_dir, name)
fig_md <- function(name) file.path("outcome", "appendix", name)
tab <- function(name) file.path(appendix_dir, name)

write_age_group_aapc_summary <- function() {
  in_file <- file.path(appendix_dir, "table_s_age_group_aapc.csv")
  if (!file.exists(in_file)) {
    stop("Missing age-group AAPC source table: ", in_file)
  }

  age_levels <- c("<20 years", "20-54 years", "55+ years")
  source_data <- readr::read_csv(in_file, show_col_types = FALSE)

  format_p <- function(p) {
    dplyr::case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "<0.001",
      TRUE ~ formatC(p, format = "f", digits = 3)
    )
  }

  rate_rows <- source_data |>
    dplyr::filter(aapc_index == "Full Range", Year == "1990~2023", Measure == "Rate") |>
    dplyr::mutate(index = factor(index, levels = age_levels)) |>
    dplyr::arrange(index) |>
    dplyr::transmute(
      `Age group` = as.character(index),
      `AAPC for prevalence rate, % (95% CI)` = gsub("~", " to ", Value, fixed = TRUE),
      `P value for rate AAPC` = format_p(p_value)
    )

  number_rows <- source_data |>
    dplyr::filter(aapc_index == "Full Range", Year == "1990~2023", Measure == "Number") |>
    dplyr::mutate(index = factor(index, levels = age_levels)) |>
    dplyr::arrange(index) |>
    dplyr::transmute(
      `Age group` = as.character(index),
      `AAPC for prevalent cases, % (95% CI)` = gsub("~", " to ", Value, fixed = TRUE),
      `P value for case-count AAPC` = format_p(p_value)
    )

  out <- dplyr::left_join(rate_rows, number_rows, by = "Age group")

  readr::write_csv(out, file.path(appendix_dir, "table_s4_age_group_aapc_summary.csv"))
  write_markdown_table(out, file.path(appendix_dir, "table_s4_age_group_aapc_summary.md"))
}

write_age_group_aapc_summary()

add_figure(
  fig_file("fig_s_sex_trend_prevalence.png"),
  fig_md("fig_s_sex_trend_prevalence.png"),
  "Supplementary Figure S1. Sex-specific global trends in LTBI prevalence, 1990-2023",
  "Supplementary Figure S1. Sex-specific trends in global latent tuberculosis infection prevalence, 1990-2023.",
  "Age-standardized prevalence rates per 100,000 population and prevalent cases in millions are presented for males (A, B) and females (C, D), respectively. Points and shaded bands denote GBD 2023 estimates with 95% uncertainty intervals; fitted lines denote joinpoint regression estimates, with segment labels indicating annual percentage change (APC) and 95% confidence intervals."
)

add_figure(
  fig_file("fig_s_age_trend_prevalence.png"),
  fig_md("fig_s_age_trend_prevalence.png"),
  "Supplementary Figure S2. Age-specific global trends in LTBI prevalence, 1990-2023",
  "Supplementary Figure S2. Age-specific trends in global latent tuberculosis infection prevalence, 1990-2023.",
  "Age-specific prevalence rates per 100,000 population and prevalent cases in millions are presented for people younger than 20 years (A, B), adults aged 20-54 years (C, D), and adults aged 55 years or older (E, F). Points and shaded bands denote GBD 2023 estimates with 95% uncertainty intervals; fitted lines denote joinpoint regression estimates, with segment labels indicating APC and 95% confidence intervals."
)

add_figure(
  fig_file("fig_s_sdi_rate_prevalence.png"),
  fig_md("fig_s_sdi_rate_prevalence.png"),
  "Supplementary Figure S3. Age-standardized LTBI prevalence rates by SDI quintile, 1990-2023",
  "Supplementary Figure S3. Age-standardized latent tuberculosis infection prevalence rates by socio-demographic index (SDI) quintile, 1990-2023.",
  "Rates are expressed per 100,000 population for (A) high-middle SDI, (B) high SDI, (C) low-middle SDI, (D) low SDI, and (E) middle SDI groups. Points and shaded bands denote GBD 2023 estimates with 95% uncertainty intervals; fitted lines denote joinpoint regression estimates, with segment labels indicating APC and 95% confidence intervals."
)

add_figure(
  fig_file("fig_s_sdi_number_prevalence.png"),
  fig_md("fig_s_sdi_number_prevalence.png"),
  "Supplementary Figure S4. LTBI prevalent cases by SDI quintile, 1990-2023",
  "Supplementary Figure S4. Number of prevalent latent tuberculosis infection cases by socio-demographic index (SDI) quintile, 1990-2023.",
  "Prevalent cases are expressed in millions for (A) high-middle SDI, (B) high SDI, (C) low-middle SDI, (D) low SDI, and (E) middle SDI groups. Points and shaded bands denote GBD 2023 estimates with 95% uncertainty intervals; fitted lines denote joinpoint regression estimates, with segment labels indicating APC and 95% confidence intervals."
)

add_figure(
  fig_file("fig_s_region_rate_prevalence.png"),
  fig_md("fig_s_region_rate_prevalence.png"),
  "Supplementary Figure S5. Age-standardized LTBI prevalence rates by WHO region, 1990-2023",
  "Supplementary Figure S5. Age-standardized latent tuberculosis infection prevalence rates by WHO region, 1990-2023.",
  "Rates are expressed per 100,000 population for (A) African Region, (B) Region of the Americas, (C) Eastern Mediterranean Region, (D) European Region, (E) South-East Asia Region, and (F) Western Pacific Region. Points and shaded bands denote GBD 2023 estimates with 95% uncertainty intervals; fitted lines denote joinpoint regression estimates, with segment labels indicating APC and 95% confidence intervals."
)

add_figure(
  fig_file("fig_s_region_number_prevalence.png"),
  fig_md("fig_s_region_number_prevalence.png"),
  "Supplementary Figure S6. LTBI prevalent cases by WHO region, 1990-2023",
  "Supplementary Figure S6. Number of prevalent latent tuberculosis infection cases by WHO region, 1990-2023.",
  "Prevalent cases are expressed in millions for (A) African Region, (B) Region of the Americas, (C) Eastern Mediterranean Region, (D) European Region, (E) South-East Asia Region, and (F) Western Pacific Region. Points and shaded bands denote GBD 2023 estimates with 95% uncertainty intervals; fitted lines denote joinpoint regression estimates, with segment labels indicating APC and 95% confidence intervals."
)

add_figure(
  fig_file("fig_s_fine_age_sensitivity.png"),
  fig_md("fig_s_fine_age_sensitivity.png"),
  "Supplementary Figure S7. Fine-age sensitivity analysis of global LTBI prevalent cases",
  "Supplementary Figure S7. Fine-age sensitivity analysis of global latent tuberculosis infection prevalent cases.",
  "The neutral within-band allocation combines GBD 2023 broad-age prevalence rates with UN population denominators regrouped as younger than 20, 20-39, 40-54, 55-69, and 70 years or older. This sensitivity analysis preserves the observed or projected prevalence rate within each broad GBD age stratum while illustrating the contribution of finer age groups to global prevalent cases."
)

add("## Supplementary Tables")
add_blank()

add_table(
  "Supplementary Table S1. Joinpoint-derived AAPC (1990-2023) and segment-specific APC estimates for global age-standardized LTBI prevalence rates and prevalent case counts, 1990-2023.",
  tab("table_s1_global_trend.md")
)

add_table(
  "Supplementary Table S2. Sensitivity of national age-standardized LTBI prevalence AAPC estimates to alternative maximum joinpoint settings, 2010-2023.",
  tab("table_s_national_joinpoint_sensitivity_summary.md")
)

add_table(
  "Supplementary Table S3. Decomposition of changes in global prevalent LTBI cases into population growth, population ageing, and epidemiologic change.",
  tab("table_s_decomposition_components.md")
)

add_table(
  "Supplementary Table S4. Age-specific AAPC estimates for global LTBI prevalence rates and prevalent case counts, 1990-2023.",
  tab("table_s4_age_group_aapc_summary.md")
)

add_table(
  "Supplementary Table S5. Fine-age sensitivity profile for global LTBI prevalent cases under equal within-band prevalence allocation.",
  tab("table_s_fine_age_sensitivity_profiles.md")
)

add_table(
  "Supplementary Table S6. Older-age contrast scenarios in the fine-age sensitivity analysis while preserving the aggregate prevalence rate among adults aged 55 years or older.",
  tab("table_s_fine_age_sensitivity_oldest_scenarios.md")
)

add_table(
  "Supplementary Table S7. Convergence diagnostics for Bayesian APC forecast models.",
  tab("table_s_forecast_convergence_diagnostics.md")
)

add_table(
  "Supplementary Table S8. PSIS-LOO comparison of Poisson and negative-binomial likelihoods in Bayesian APC forecast models.",
  tab("table_s_forecast_loo_model_comparison.md")
)

add_table(
  "Supplementary Table S9. Posterior predictive checks for Poisson and negative-binomial Bayesian APC forecast models.",
  tab("table_s_forecast_posterior_predictive_checks.md")
)

add_table(
  "Supplementary Table S10. Forecast summary from the negative-binomial Bayesian APC model for 2024 and 2050.",
  tab("table_s_forecast_negbin_summary.md")
)

remove_trailing_page_break()
writeLines(lines_out, appendix_path, useBytes = TRUE)
message("Wrote ", appendix_path)
