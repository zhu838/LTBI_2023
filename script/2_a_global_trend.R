# loading packages --------------------------------------------------------

set_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
    setwd(dirname(dirname(script_path)))
  }
}

set_project_root()

# devtools::install_github("DanChaltiel/nih.joinpoint")
library(nih.joinpoint)
library(segmented)
library(tidyverse)
library(patchwork)
library(paletteer)
library(Cairo)

source("./script/config.R")
source("./script/function.R")
source("./script/joinpoint_setting.R")

dir.create(outcome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

# data --------------------------------------------------------------------

df_number <- readr::read_csv(file.path(database_dir, "global_regional_number.csv"), show_col_types = FALSE)
df_rate <- readr::read_csv(file.path(database_dir, "global_regional_rate.csv"), show_col_types = FALSE)

measures_available <- intersect(target_measures, sort(unique(c(df_number$measure_name, df_rate$measure_name))))
if (length(measures_available) == 0) {
  stop("No target measures found in prepared data.")
}

message("Measures used in global trend: ", paste(measures_available, collapse = ", "))

df_global_rate <- df_rate |>
  filter(
    location_name == target_location_global,
    age_name == target_age_global,
    measure_name %in% measures_available
  ) |>
  arrange(measure_name, sex_name, year)

# Number data has a GBD-provided `All ages` total (there is no
# age-standardized count). Use it directly: summing all rows would double-count
# because the export also contains its component age-specific rows.
df_global_number <- df_number |>
  filter(
    location_name == target_location_global,
    age_name == "All ages",
    measure_name %in% measures_available
  ) |>
  arrange(measure_name, sex_name, year) |>
  mutate(across(c(val, lower, upper), ~round(., 0)))

if (nrow(df_global_rate) == 0 || nrow(df_global_number) == 0) {
  stop("No rows after global filters. Check config.R target_age/sex/location.")
}

year_min <- min(c(df_global_rate$year, df_global_number$year), na.rm = TRUE)
year_max <- max(c(df_global_rate$year, df_global_number$year), na.rm = TRUE)
export_opt_new <- build_export_opt(year_min, year_max)

make_ylabel <- function(measure, metric) {
  if (metric == "Rate") {
    if (measure == "Prevalence") return("Prevalence rate")
    return(paste0(measure, " rate"))
  }
  if (metric == "Number") {
    if (measure == "Prevalence") return("Number of prevalent cases (million)")
    return(paste0("Number of ", measure, " (million)"))
  }
  paste(measure, metric)
}

fit_joinpoint <- function(d, metric_name) {
  run_opt <- if (metric_name == "Rate") run_opt_rate else run_opt_number
  joinpoint(
    d,
    year,
    val,
    run_opt = run_opt,
    export_opt = export_opt_new
  )
}

# model + plots -----------------------------------------------------------

aapc_rows <- list()
apc_rows <- list()

sex_levels_target <- c("Both", "Male", "Female")
sex_levels_data <- unique(c(df_global_rate$sex_name, df_global_number$sex_name))
sex_levels_available <- sex_levels_target[sex_levels_target %in% sex_levels_data]
if (length(sex_levels_available) == 0) {
  stop("No sex rows found after global filters.")
}
message("Sex levels plotted in Figure 1: ", paste(sex_levels_available, collapse = ", "))

plot_list_measure <- list()
plot_list_measure_sex <- list()

for (m in measures_available) {
  # AAPC tables: keep consistent with previous behavior (Both sex only)
  d_rate_both <- df_global_rate |>
    filter(measure_name == m, sex_name == target_sex_global) |>
    select(year, val, lower, upper)
  d_number_both <- df_global_number |>
    filter(measure_name == m, sex_name == target_sex_global) |>
    select(year, val, lower, upper)

  model_rate_both <- fit_joinpoint(d_rate_both, "Rate")
  model_number_both <- fit_joinpoint(d_number_both, "Number")

  aapc_rows[[paste0(m, "_rate")]] <- get_aapc(model_rate_both) |>
    mutate(var = m, Measure = "Rate")
  aapc_rows[[paste0(m, "_number")]] <- get_aapc(model_number_both) |>
    mutate(var = m, Measure = "Number")
  apc_rows[[paste0(m, "_rate")]] <- get_apc(model_rate_both) |>
    mutate(var = m, Measure = "Rate")
  apc_rows[[paste0(m, "_number")]] <- get_apc(model_number_both) |>
    mutate(var = m, Measure = "Number")

  # Figure 1 (main): overall (Both) only
  plot_rate_both <- plot_apc(model_rate_both, d_rate_both, use_scientific_10 = TRUE) +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
    labs(
      title = "A",
      x = "Year",
      y = make_ylabel(m, "Rate"),
      fill = "APC (95% CI)"
    )

  plot_number_both <- plot_apc(model_number_both, d_number_both, use_scientific_10 = FALSE, y_divisor = 1e6) +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
    labs(
      title = "B",
      x = "Year",
      y = make_ylabel(m, "Number"),
      fill = "APC (95% CI)"
    )

  plot_list_measure[[m]] <- patchwork::wrap_plots(
    list(plot_rate_both, plot_number_both),
    ncol = 2
  ) + patchwork::plot_layout(widths = c(1, 1))

  # Supplementary figure (appendix): rate + number, stratified by sex
  # Exclude overall 'Both' because it's shown in the main figure
  sex_for_appendix <- setdiff(sex_levels_available, target_sex_global)
  if (length(sex_for_appendix) == 0) {
    message("Skipping appendix sex-stratified plot for ", m, " (only 'Both' available)")
  } else {
    plot_sex_rate <- list()
    plot_sex_number <- list()

    for (i in seq_along(sex_for_appendix)) {
      sx <- sex_for_appendix[[i]]
      letter_rate <- LETTERS[(i - 1) * 2 + 1]
      letter_number <- LETTERS[(i - 1) * 2 + 2]

      d_rate <- df_global_rate |>
        filter(measure_name == m, sex_name == sx) |>
        select(year, val, lower, upper)
      d_number <- df_global_number |>
        filter(measure_name == m, sex_name == sx) |>
        select(year, val, lower, upper)

      model_rate <- fit_joinpoint(d_rate, "Rate")
      model_number <- fit_joinpoint(d_number, "Number")

      message(
        sprintf(
          "JP segments | %s | %s | Rate=%d | Number=%d",
          m,
          sx,
          nrow(model_rate$apc),
          nrow(model_number$apc)
        )
      )

      plot_sex_rate[[sx]] <- plot_apc(model_rate, d_rate, use_scientific_10 = TRUE) +
        guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
        labs(
          title = paste0(letter_rate, ". ", sx),
          x = "Year",
          y = make_ylabel(m, "Rate"),
          fill = "APC (95% CI)"
        )

      plot_sex_number[[sx]] <- plot_apc(model_number, d_number, use_scientific_10 = FALSE, y_divisor = 1e6) +
        guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
        labs(
          title = paste0(letter_number, ". ", sx),
          x = "Year",
          y = make_ylabel(m, "Number"),
          fill = "APC (95% CI)"
        )
    }

    # Arrange as rows (sex) × 2 columns (Rate/Number)
    fig_m <- patchwork::wrap_plots(plot_sex_rate, ncol = 1) |
      patchwork::wrap_plots(plot_sex_number, ncol = 1)
    fig_m <- fig_m + patchwork::plot_layout(widths = c(1, 1))

    plot_list_measure_sex[[m]] <- fig_m
  }
}

# Figure 1 (main): overall trend (Both) only
fig_1 <- patchwork::wrap_plots(plot_list_measure, ncol = length(plot_list_measure))

ggsave(
  file.path(outcome_dir, "fig_1_global_trend.png"),
  plot = fig_1,
  width = 14,
  height = 4
)

# Supplementary: sex-stratified global trend (combined)
if (length(plot_list_measure_sex) > 0) {
  n_panels <- length(plot_list_measure_sex)
  ncol_s <- min(2, n_panels)
  nrow_s <- ceiling(n_panels / ncol_s)

  fig_s_sex <- patchwork::wrap_plots(plot_list_measure_sex, ncol = ncol_s)

  ggsave(
    file.path(appendix_dir, "fig_s_sex_trend_prevalence.png"),
    plot = fig_s_sex,
    width = 14,
    height = 8 * nrow_s
  )
}

# AAPC outputs ------------------------------------------------------------

df_aapc <- bind_rows(aapc_rows)
df_apc <- bind_rows(apc_rows)

readr::write_csv(df_aapc, file.path(appendix_dir, "table_s1_global_trend.csv"))

markdown_table <- bind_rows(
  df_aapc |>
    filter(aapc_index == "Full Range") |>
    transmute(
      Outcome = var,
      Measure,
      Statistic = "AAPC",
      Interval = paste(start_obs, end_obs, sep = " to "),
      `Estimate (95% CI)` = gsub("~", " to ", Value, fixed = TRUE),
      `P value` = case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ formatC(p_value, format = "f", digits = 3)
      )
    ),
  df_apc |>
    transmute(
      Outcome = var,
      Measure,
      Statistic = "APC",
      Interval,
      `Estimate (95% CI)`,
      `P value` = p_value_label
    )
) |>
  mutate(Measure = factor(Measure, levels = c("Rate", "Number"))) |>
  arrange(Measure, Statistic, Interval) |>
  knitr::kable(format = "markdown")

write(markdown_table, file.path(appendix_dir, "table_s1_global_trend.md"))

message("Done: global trend figure(s) + AAPC tables written to ", outcome_dir)

