#####################################
## @Description: National joinpoint sensitivity analysis
## @Author: Li Kangguo
## @Date: 2026-03-06
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

library(nih.joinpoint)
library(tidyverse)

source("./script/config.R")
source("./script/function.R")
source("./script/joinpoint_setting.R")

dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

measure_target <- "Prevalence"
sensitivity_max_joinpoints <- c(0, 1, 2, 3)
outlier_locations <- c("United States of America", "Sri Lanka")

df_map_iso <- readr::read_csv(iso_code_file, show_col_types = FALSE)

national_dir <- file.path(database_dir, "national_by_year")
csvs <- list.files(national_dir, pattern = "^national_.*\\.csv$", full.names = TRUE)
df_raw_rate <- purrr::map_dfr(csvs, ~ readr::read_csv(.x, show_col_types = FALSE))

df_all_rate <- df_raw_rate |>
  filter(
    age_name == target_age_global,
    sex_name == "Both",
    location %in% df_map_iso$location_id,
    measure_name == measure_target,
    year %in% 2010:2023
  ) |>
  select(location, location_name, measure_name, year, val, lower, upper) |>
  arrange(location_name, year)

year_start <- min(df_all_rate$year, na.rm = TRUE)
year_end <- max(df_all_rate$year, na.rm = TRUE)

fit_by_country_rate <- function(d, max_joinpoints) {
  joinpoint(
    d,
    year,
    val,
    by = location_name,
    run_opt = build_run_options(max_joinpoints, "age-adjusted rate"),
    export_opt = build_export_opt(year_start, year_end)
  )
}

run_sensitivity <- function(max_joinpoints) {
  model_rate <- fit_by_country_rate(df_all_rate, max_joinpoints)

  get_aapc(model_rate) |>
    mutate(
      max_joinpoints_setting = max_joinpoints,
      selected_joinpoints = as.integer(joinpoint_model),
      aapc = as.numeric(aapc),
      aapc_c_i_low = as.numeric(aapc_c_i_low),
      aapc_c_i_high = as.numeric(aapc_c_i_high)
    ) |>
    filter(aapc_index == "Full Range")
}

df_sensitivity <- purrr::map_dfr(sensitivity_max_joinpoints, run_sensitivity)

distribution_summary <- df_sensitivity |>
  count(max_joinpoints_setting, selected_joinpoints, name = "countries") |>
  group_by(max_joinpoints_setting) |>
  mutate(proportion = countries / sum(countries)) |>
  ungroup() |>
  arrange(max_joinpoints_setting, selected_joinpoints)

distribution_compact <- distribution_summary |>
  mutate(item = sprintf("%d: %d (%.1f%%)", selected_joinpoints, countries, 100 * proportion)) |>
  group_by(max_joinpoints_setting) |>
  summarise(selected_joinpoints_distribution = paste(item, collapse = "; "), .groups = "drop")

baseline_max_joinpoints <- max(sensitivity_max_joinpoints)
baseline_aapc <- df_sensitivity |>
  filter(max_joinpoints_setting == baseline_max_joinpoints) |>
  select(location_name, aapc_baseline = aapc)

comparison_summary <- df_sensitivity |>
  left_join(baseline_aapc, by = "location_name") |>
  mutate(delta_vs_baseline = aapc - aapc_baseline) |>
  group_by(max_joinpoints_setting) |>
  summarise(
    countries = n(),
    median_aapc = median(aapc, na.rm = TRUE),
    iqr_low = quantile(aapc, probs = 0.25, na.rm = TRUE),
    iqr_high = quantile(aapc, probs = 0.75, na.rm = TRUE),
    median_abs_delta_vs_baseline = median(abs(delta_vs_baseline), na.rm = TRUE),
    max_abs_delta_vs_baseline = max(abs(delta_vs_baseline), na.rm = TRUE),
    .groups = "drop"
  )

outlier_summary <- df_sensitivity |>
  filter(location_name %in% outlier_locations) |>
  mutate(
    aapc_95_ci = sprintf("%.2f (%.2f to %.2f)", aapc, aapc_c_i_low, aapc_c_i_high)
  ) |>
  select(
    location_name,
    max_joinpoints_setting,
    selected_joinpoints,
    aapc_95_ci,
    p_value,
    p_value_label
  ) |>
  arrange(location_name, max_joinpoints_setting)

summary_markdown <- comparison_summary |>
  left_join(distribution_compact, by = "max_joinpoints_setting") |>
  transmute(
    `Max. allowed joinpoints` = max_joinpoints_setting,
    `Selected joinpoints, n (%)` = selected_joinpoints_distribution,
    `Countries, n` = countries,
    `Median AAPC (%)` = sprintf("%.2f", median_aapc),
    `IQR of AAPC (%)` = sprintf("%.2f to %.2f", iqr_low, iqr_high),
    `Median absolute change vs primary model` = sprintf("%.3f", median_abs_delta_vs_baseline),
    `Maximum absolute change vs primary model` = sprintf("%.3f", max_abs_delta_vs_baseline)
  )

outlier_markdown <- outlier_summary |>
  transmute(
    Country = location_name,
    `Max. allowed joinpoints` = max_joinpoints_setting,
    `Selected joinpoints` = selected_joinpoints,
    `AAPC (95% CI)` = aapc_95_ci,
    `P value` = as.character(p_value)
  )

readr::write_csv(
  df_sensitivity,
  file.path(appendix_dir, "table_s_national_joinpoint_sensitivity_all.csv")
)
readr::write_csv(
  distribution_summary,
  file.path(appendix_dir, "table_s_national_joinpoint_sensitivity_distribution.csv")
)
readr::write_csv(
  comparison_summary,
  file.path(appendix_dir, "table_s_national_joinpoint_sensitivity_summary.csv")
)
readr::write_csv(
  outlier_summary,
  file.path(appendix_dir, "table_s_national_joinpoint_sensitivity_outliers.csv")
)

write_markdown_table(
  summary_markdown,
  file.path(appendix_dir, "table_s_national_joinpoint_sensitivity_summary.md")
)

write_markdown_table(
  outlier_markdown,
  file.path(appendix_dir, "table_s_national_joinpoint_sensitivity_outliers.md")
)

message("Done: national joinpoint sensitivity outputs written to ", appendix_dir)
