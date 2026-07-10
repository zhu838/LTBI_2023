# SDI + WHO region trends (adapted from adult_pertussis/script/2_SDI_trend.R)

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
library(paletteer)
library(Cairo)

source("./script/config.R")
source("./script/function.R")
source("./script/joinpoint_setting.R")

dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

df_number_raw <- readr::read_csv(file.path(database_dir, "global_regional_number.csv"), show_col_types = FALSE)
df_rate <- readr::read_csv(file.path(database_dir, "global_regional_rate.csv"), show_col_types = FALSE)

measures_available <- intersect(target_measures, sort(unique(c(df_number_raw$measure_name, df_rate$measure_name))))
if (length(measures_available) == 0) stop("No target measures found in prepared data.")

# location groups derived from current exports
locations_sdi <- df_rate |>
  distinct(location_name) |>
  filter(str_detect(location_name, regex("SDI", ignore_case = TRUE))) |>
  pull(location_name) |>
  sort()

# The GBD export uses "Region of the Americas", which does not end in
# "Region". Keep an explicit WHO-region order to avoid dropping that group.
who_regions <- c(
  "African Region",
  "Region of the Americas",
  "Eastern Mediterranean Region",
  "European Region",
  "South-East Asia Region",
  "Western Pacific Region"
)
locations_region <- who_regions[who_regions %in% unique(df_rate$location_name)]

locations_all <- c(locations_sdi, locations_region)
if (length(locations_all) == 0) stop("No SDI/Region locations found in data.")

# Use the GBD-provided All ages total. Summing all age rows would double-count
# because the export contains both the All ages record and its age-specific rows.
df_number <- df_number_raw |>
  filter(age_name == "All ages")

df_rate_main <- df_rate |>
  filter(sex_name == "Both", age_name == target_age_global, location_name %in% locations_all, measure_name %in% measures_available) |>
  arrange(measure_name, location_name, year)

df_number_main <- df_number |>
  filter(sex_name == "Both", location_name %in% locations_all, measure_name %in% measures_available) |>
  arrange(measure_name, location_name, year)

year_min <- min(c(df_rate_main$year, df_number_main$year), na.rm = TRUE)
year_max <- max(c(df_rate_main$year, df_number_main$year), na.rm = TRUE)
export_opt_new <- build_export_opt(year_min, year_max)

run_by_location <- function(d, metric_type) {
  run_opt <- if (metric_type == "Rate") run_opt_rate else run_opt_number
  joinpoint(d, year, val, by = location_name, run_opt = run_opt, export_opt = export_opt_new)
}

prep_for_jp <- function(d) {
  d |>
    dplyr::filter(!is.na(year), !is.na(val), is.finite(val), val > 0) |>
    dplyr::arrange(year) |>
    dplyr::select(year, val, lower, upper)
}

fit_joinpoint_single <- function(d, metric_type) {
  run_opt <- if (metric_type == "Rate") run_opt_rate else run_opt_number
  joinpoint(d, year, val, run_opt = run_opt, export_opt = export_opt_new)
}

for (m in measures_available) {
  model_rate <- run_by_location(filter(df_rate_main, measure_name == m), "Rate")
  model_number <- run_by_location(filter(df_number_main, measure_name == m), "Number")

  df_aapc <- bind_rows(
    get_aapc(model_number) |> mutate(Label = m, Measure = "Number"),
    get_aapc(model_rate) |> mutate(Label = m, Measure = "Rate")
  )

  readr::write_csv(df_aapc, file.path(appendix_dir, paste0("table_s_sdi_region_aapc_", tolower(gsub(" ", "_", m)), ".csv")))

  sdi_name <- locations_sdi
  reg_name <- locations_region
  if (length(sdi_name) > 0) {
    # Rate figure (SDI) with APC
    fig_sdi_rate <- list()
    for (i in seq_along(sdi_name)) {
      loc <- sdi_name[[i]]
      letter <- LETTERS[i]
      d_rate <- df_rate_main |>
        dplyr::filter(measure_name == m, location_name == loc) |>
        prep_for_jp()

      if (nrow(d_rate) < 3) {
        warning("Too few rows for joinpoint (Rate): ", m, " | ", loc)
        p_rate <- plot_val(
          data = df_rate_main,
          measure = m,
          filter_col = "location_name",
          filter_val = loc,
          ylab = paste0(m, " rate")
        )
        fig_sdi_rate[[loc]] <- p_rate +
          ggplot2::labs(
            title = paste0(letter, ". ", loc),
            x = "Year"
          )
      } else {
        model_rate_loc <- fit_joinpoint_single(d_rate, "Rate")
        fig_sdi_rate[[loc]] <- plot_apc(model_rate_loc, d_rate, use_scientific_10 = TRUE) +
          guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
          labs(
            title = paste0(letter, ". ", loc),
            x = "Year",
            y = paste0(m, " rate"),
            fill = "APC (95% CI)"
          )
      }
    }
    fig_sdi_rate <- patchwork::wrap_plots(fig_sdi_rate, ncol = 2)
    ggsave(
      file.path(appendix_dir, paste0("fig_s_sdi_rate_", tolower(gsub(" ", "_", m)), ".png")),
      plot = fig_sdi_rate,
      width = 14,
      height = 12
    )

    # Number figure (SDI) with APC
    fig_sdi_num <- list()
    for (i in seq_along(sdi_name)) {
      loc <- sdi_name[[i]]
      letter <- LETTERS[i]
      d_num <- df_number_main |>
        dplyr::filter(measure_name == m, location_name == loc) |>
        prep_for_jp()

      if (nrow(d_num) < 3) {
        warning("Too few rows for joinpoint (Number): ", m, " | ", loc)
        p_num <- plot_val(
          data = df_number_main,
          measure = m,
          filter_col = "location_name",
          filter_val = loc,
          ylab = paste0("Number of ", m)
        )
        fig_sdi_num[[loc]] <- p_num +
          ggplot2::labs(
            title = paste0(letter, ". ", loc),
            x = "Year"
          )
      } else {
        model_num_loc <- fit_joinpoint_single(d_num, "Number")
        fig_sdi_num[[loc]] <- plot_apc(model_num_loc, d_num, use_scientific_10 = FALSE, y_divisor = 1e6) +
          guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
          labs(
            title = paste0(letter, ". ", loc),
            x = "Year",
            y = paste0("Number of ", m, " (million)"),
            fill = "APC (95% CI)"
          )
      }
    }
    fig_sdi_num <- patchwork::wrap_plots(fig_sdi_num, ncol = 2)
    ggsave(
      file.path(appendix_dir, paste0("fig_s_sdi_number_", tolower(gsub(" ", "_", m)), ".png")),
      plot = fig_sdi_num,
      width = 14,
      height = 12
    )
  }
  if (length(reg_name) > 0) {
    # Rate figure (Region) with APC
    fig_reg_rate <- list()
    for (i in seq_along(reg_name)) {
      loc <- reg_name[[i]]
      letter <- LETTERS[i]
      d_rate <- df_rate_main |>
        dplyr::filter(measure_name == m, location_name == loc) |>
        prep_for_jp()

      if (nrow(d_rate) < 3) {
        warning("Too few rows for joinpoint (Rate): ", m, " | ", loc)
        p_rate <- plot_val(
          data = df_rate_main,
          measure = m,
          filter_col = "location_name",
          filter_val = loc,
          ylab = paste0(m, " rate")
        )
        fig_reg_rate[[loc]] <- p_rate +
          ggplot2::labs(
            title = paste0(letter, ". ", loc),
            x = "Year"
          )
      } else {
        model_rate_loc <- fit_joinpoint_single(d_rate, "Rate")
        fig_reg_rate[[loc]] <- plot_apc(model_rate_loc, d_rate, use_scientific_10 = TRUE) +
          guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
          labs(
            title = paste0(letter, ". ", loc),
            x = "Year",
            y = paste0(m, " rate"),
            fill = "APC (95% CI)"
          )
      }
    }
    fig_reg_rate <- patchwork::wrap_plots(fig_reg_rate, ncol = 2)
    ggsave(
      file.path(appendix_dir, paste0("fig_s_region_rate_", tolower(gsub(" ", "_", m)), ".png")),
      plot = fig_reg_rate,
      width = 14,
      height = 12
    )

    # Number figure (Region) with APC
    fig_reg_num <- list()
    for (i in seq_along(reg_name)) {
      loc <- reg_name[[i]]
      letter <- LETTERS[i]
      d_num <- df_number_main |>
        dplyr::filter(measure_name == m, location_name == loc) |>
        prep_for_jp()

      if (nrow(d_num) < 3) {
        warning("Too few rows for joinpoint (Number): ", m, " | ", loc)
        p_num <- plot_val(
          data = df_number_main,
          measure = m,
          filter_col = "location_name",
          filter_val = loc,
          ylab = paste0("Number of ", m)
        )
        fig_reg_num[[loc]] <- p_num +
          ggplot2::labs(
            title = paste0(letter, ". ", loc),
            x = "Year"
          )
      } else {
        model_num_loc <- fit_joinpoint_single(d_num, "Number")
        fig_reg_num[[loc]] <- plot_apc(model_num_loc, d_num, use_scientific_10 = FALSE, y_divisor = 1e6) +
          guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
          labs(
            title = paste0(letter, ". ", loc),
            x = "Year",
            y = paste0("Number of ", m, " (million)"),
            fill = "APC (95% CI)"
          )
      }
    }
    fig_reg_num <- patchwork::wrap_plots(fig_reg_num, ncol = 2)
    ggsave(
      file.path(appendix_dir, paste0("fig_s_region_number_", tolower(gsub(" ", "_", m)), ".png")),
      plot = fig_reg_num,
      width = 14,
      height = 12
    )
  }
}

message("Done: SDI/Region joinpoint + plots -> ", appendix_dir)
