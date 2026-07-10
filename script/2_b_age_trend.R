#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-03-04 19:17:33
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-04 20:52:47
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
library(segmented)
library(tidyverse)
library(patchwork)

source("./script/config.R")
source("./script/function.R")
source("./script/joinpoint_setting.R")

dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

df_number <- readr::read_csv(file.path(database_dir, "global_regional_number.csv"), show_col_types = FALSE)
df_rate <- readr::read_csv(file.path(database_dir, "global_regional_rate.csv"), show_col_types = FALSE)

measures_available <- intersect(target_measures, sort(unique(c(df_number$measure_name, df_rate$measure_name))))
if (length(measures_available) == 0) {
  stop("No target measures found in prepared data.")
}

age_groups <- intersect(age_groups_for_age_trend, sort(unique(c(df_rate$age_name, df_number$age_name))))
if (length(age_groups) == 0) {
  stop("No configured age groups exist in data. Check config.R age_groups_for_age_trend")
}

df_global_rate <- df_rate |>
  filter(
    location_name == target_location_global,
    sex_name == target_sex_global,
    age_name %in% age_groups,
    measure_name %in% measures_available
  ) |>
  arrange(measure_name, age_name, year)

df_global_number <- df_number |>
  filter(
    location_name == target_location_global,
    sex_name == target_sex_global,
    age_name %in% age_groups,
    measure_name %in% measures_available
  ) |>
  arrange(measure_name, age_name, year)

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

year_min <- min(c(df_global_rate$year, df_global_number$year), na.rm = TRUE)
year_max <- max(c(df_global_rate$year, df_global_number$year), na.rm = TRUE)
export_opt_new <- build_export_opt(year_min, year_max)

fit_joinpoint <- function(d, metric_name) {
  run_opt <- if (metric_name == "Rate") run_opt_age_specific_rate else run_opt_number
  joinpoint(
    d,
    year,
    val,
    run_opt = run_opt,
    export_opt = export_opt_new
  )
}

prep_for_jp <- function(d) {
  d |>
    dplyr::filter(!is.na(year), !is.na(val), is.finite(val), val > 0) |>
    dplyr::arrange(year) |>
    dplyr::select(year, val, lower, upper)
}

for (m in measures_available) {
  panels <- list()
  for (i in seq_along(age_groups)) {
    ag <- age_groups[[i]]
    letter_rate <- LETTERS[(i - 1) * 2 + 1]
    letter_number <- LETTERS[(i - 1) * 2 + 2]

    d_rate <- df_global_rate |>
      dplyr::filter(measure_name == m, age_name == ag) |>
      prep_for_jp()

    if (nrow(d_rate) < 3) {
      warning("Too few rows for joinpoint (Rate): ", m, " | ", ag)
      p_rate <- plot_val(
        data = df_global_rate,
        measure = m,
        filter_col = "age_name",
        filter_val = ag,
        ylab = make_ylabel(m, "Rate")
      )
      panels[[length(panels) + 1]] <- p_rate +
        ggplot2::labs(
          title = paste0(letter_rate, ". ", ag),
          x = "Year",
          y = paste0(make_ylabel(m, "Rate"), "\n", ag)
        )
    } else {
      model_rate <- fit_joinpoint(d_rate, "Rate")
      panels[[length(panels) + 1]] <- plot_apc(model_rate, d_rate, use_scientific_10 = TRUE) +
        guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
        labs(
          title = paste0(letter_rate, ". ", ag),
          x = "Year",
          y = paste0(make_ylabel(m, "Rate"), "\n", ag),
          fill = "APC (95% CI)"
        )
    }

    d_number <- df_global_number |>
      dplyr::filter(measure_name == m, age_name == ag) |>
      prep_for_jp()

    if (nrow(d_number) < 3) {
      warning("Too few rows for joinpoint (Number): ", m, " | ", ag)
      p_number <- plot_val(
        data = df_global_number,
        measure = m,
        filter_col = "age_name",
        filter_val = ag,
        ylab = make_ylabel(m, "Number")
      )
      panels[[length(panels) + 1]] <- p_number +
        ggplot2::labs(
          title = paste0(letter_number, ". ", ag),
          x = "Year",
          y = paste0(make_ylabel(m, "Number"), "\n", ag)
        )
    } else {
      model_number <- fit_joinpoint(d_number, "Number")
      panels[[length(panels) + 1]] <- plot_apc(model_number, d_number, use_scientific_10 = FALSE, y_divisor = 1e6) +
        guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
        labs(
          title = paste0(letter_number, ". ", ag),
          x = "Year",
          y = paste0(make_ylabel(m, "Number"), "\n", ag),
          fill = "APC (95% CI)"
        )
    }
  }

  fig_age <- patchwork::wrap_plots(panels, ncol = 2)
  ggsave(
    file.path(appendix_dir, paste0("fig_s_age_trend_", tolower(gsub(" ", "_", m)), ".png")),
    plot = fig_age,
    width = 14,
    height = 10
  )
}

# -----------------------------------------------------------------
# AAPC by age group (Global, Both sex, Prevalence)
# Saved for reuse in Table 1 and other summaries
# -----------------------------------------------------------------

measure_target <- "Prevalence"

df_age_rate_all <- df_global_rate |>
  dplyr::filter(measure_name == measure_target) |>
  dplyr::rename(Index = age_name)

df_age_number_all <- df_global_number |>
  dplyr::filter(measure_name == measure_target) |>
  dplyr::rename(Index = age_name)

if (nrow(df_age_rate_all) > 0 && nrow(df_age_number_all) > 0) {
  model_rate_age <- joinpoint(
    df_age_rate_all,
    year,
    val,
    by = Index,
    run_opt = run_opt_age_specific_rate,
    export_opt = export_opt_new
  )

  model_number_age <- joinpoint(
    df_age_number_all,
    year,
    val,
    by = Index,
    run_opt = run_opt_number,
    export_opt = export_opt_new
  )

  df_aapc_age <- dplyr::bind_rows(
    get_aapc(model_number_age) |>
      dplyr::mutate(Label = measure_target, Measure = "Number"),
    get_aapc(model_rate_age) |>
      dplyr::mutate(Label = measure_target, Measure = "Rate")
  )

  readr::write_csv(
    df_aapc_age,
    file.path(appendix_dir, "table_s_age_group_aapc.csv")
  )
}

message("Done: age trend figures and age-group AAPC written to ", appendix_dir)
