#####################################
## @Description: Fine-age sensitivity analysis within GBD LTBI age strata
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
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)

source("./script/config.R")
source("./script/function.R")

dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

fine_age_order <- c("<20", "20-39", "40-54", "55-69", "70+")
broad_age_order <- c("<20", "20-54", "55+")

scenarios <- tibble(
  scenario = c(
    "older_lower_0.8",
    "equal_within_band_1.0",
    "older_higher_1.2",
    "older_much_higher_1.5"
  ),
  older_to_55_69_rate_ratio = c(0.8, 1.0, 1.2, 1.5)
)

fine_age_group <- function(age_start) {
  case_when(
    age_start < 20 ~ "<20",
    age_start >= 20 & age_start < 40 ~ "20-39",
    age_start >= 40 & age_start < 55 ~ "40-54",
    age_start >= 55 & age_start < 70 ~ "55-69",
    age_start >= 70 ~ "70+",
    TRUE ~ NA_character_
  )
}

broad_age_group <- function(age_start) {
  case_when(
    age_start < 20 ~ "<20",
    age_start >= 20 & age_start < 55 ~ "20-54",
    age_start >= 55 ~ "55+",
    TRUE ~ NA_character_
  )
}

load_fine_population <- function(years) {
  pop_file <- "./Data/unpopulation_dataportal_20260305182132.csv"
  iso <- readr::read_csv(iso_code_file, show_col_types = FALSE)
  iso3_set <- toupper(as.character(iso$ISO3))

  readr::read_csv(
    pop_file,
    col_select = c(Iso3, Time, Variant, Sex, AgeStart, Value),
    show_col_types = FALSE
  ) |>
    mutate(
      Iso3 = toupper(as.character(Iso3)),
      year = as.integer(Time),
      age_start = as.integer(AgeStart),
      population = as.numeric(Value)
    ) |>
    filter(
      Iso3 %in% iso3_set,
      year %in% years,
      Variant == "Median",
      Sex == "Both sexes"
    ) |>
    mutate(
      fine_age_group = fine_age_group(age_start),
      broad_age_group = broad_age_group(age_start)
    ) |>
    filter(!is.na(fine_age_group), !is.na(broad_age_group)) |>
    group_by(year, broad_age_group, fine_age_group) |>
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop")
}

load_observed_broad_rates <- function(years) {
  readr::read_csv(file.path(database_dir, "global_regional_rate.csv"), show_col_types = FALSE) |>
    filter(
      measure_name == "Prevalence",
      metric_name == "Rate",
      location_name == target_location_global,
      sex_name == target_sex_global,
      age_name %in% c("<20 years", "20-54 years", "55+ years"),
      year %in% years
    ) |>
    transmute(
      year = as.integer(year),
      broad_age_group = recode(
        age_name,
        "<20 years" = "<20",
        "20-54 years" = "20-54",
        "55+ years" = "55+"
      ),
      prevalence_per_100k = as.numeric(val),
      source = "GBD observed"
    )
}

load_projected_broad_rates <- function(years) {
  candidates <- c(
    file.path(outcome_dir, "forecast", "predictions_rstan_negbin_Global_2024_2050.csv"),
    file.path(outcome_dir, "forecast", "predictions_rstan_Global_2024_2050_corrected_pop.csv"),
    file.path(outcome_dir, "forecast", "predictions_rstan_Global_2024_2050.csv")
  )
  forecast_file <- candidates[file.exists(candidates)][[1]]
  if (is.na(forecast_file)) {
    stop("No forecast prediction CSV found.")
  }

  readr::read_csv(forecast_file, show_col_types = FALSE) |>
    filter(year %in% years) |>
    transmute(
      year = as.integer(year),
      broad_age_group = as.character(age_group),
      prevalence_per_100k = as.numeric(pred_prev_mean),
      source = "APC projected"
    )
}

split_rates_for_year <- function(pop_year, rate_year, older_ratio) {
  purrr::map_dfr(broad_age_order, function(broad) {
    pop_b <- pop_year |> filter(broad_age_group == broad)
    if (nrow(pop_b) == 0) {
      return(tibble())
    }

    broad_rate <- rate_year |>
      filter(broad_age_group == broad) |>
      pull(prevalence_per_100k)

    if (length(broad_rate) != 1) {
      stop("Missing or duplicated broad-age rate for ", broad, ".")
    }

    if (broad == "55+") {
      pop_55_69 <- pop_b |> filter(fine_age_group == "55-69") |> pull(population) |> sum(na.rm = TRUE)
      pop_70 <- pop_b |> filter(fine_age_group == "70+") |> pull(population) |> sum(na.rm = TRUE)
      total <- pop_55_69 + pop_70
      rate_55_69 <- broad_rate / ((pop_55_69 / total) + older_ratio * (pop_70 / total))
      pop_b <- pop_b |>
        mutate(
          prevalence_per_100k = if_else(
            fine_age_group == "70+",
            older_ratio * rate_55_69,
            rate_55_69
          )
        )
    } else {
      pop_b <- pop_b |>
        mutate(prevalence_per_100k = broad_rate)
    }

    pop_b |>
      mutate(cases = prevalence_per_100k / 1e5 * population)
  })
}

build_sensitivity_profiles <- function() {
  years_observed <- c(1990, 2023)
  years_projected <- c(2024, 2050)
  all_years <- c(years_observed, years_projected)

  pop <- load_fine_population(all_years)
  rates <- bind_rows(
    load_observed_broad_rates(years_observed),
    load_projected_broad_rates(years_projected)
  )

  profiles <- purrr::pmap_dfr(
    tidyr::expand_grid(scenarios, year = all_years),
    function(scenario, older_to_55_69_rate_ratio, year) {
      split_rates_for_year(
        pop |> filter(year == !!year),
        rates |> filter(year == !!year),
        older_to_55_69_rate_ratio
      ) |>
        mutate(
          scenario = scenario,
          older_to_55_69_rate_ratio = older_to_55_69_rate_ratio,
          source = rates |> filter(year == !!year) |> pull(source) |> unique() |> first()
        )
    }
  )

  profiles |>
    group_by(scenario, year) |>
    mutate(
      total_cases = sum(cases, na.rm = TRUE),
      case_share = cases / total_cases
    ) |>
    ungroup() |>
    mutate(
      fine_age_group = factor(fine_age_group, levels = fine_age_order),
      broad_age_group = factor(broad_age_group, levels = broad_age_order),
      scenario = factor(scenario, levels = scenarios$scenario)
    ) |>
    arrange(scenario, year, fine_age_group)
}

format_neutral_table <- function(profiles) {
  profiles |>
    filter(scenario == "equal_within_band_1.0") |>
    mutate(
      `Population (million)` = sprintf("%.2f", population / 1e6),
      `Prevalence per 100,000` = sprintf("%.2f", prevalence_per_100k),
      `Cases (million)` = sprintf("%.2f", cases / 1e6),
      `Case share (%)` = sprintf("%.1f", 100 * case_share)
    ) |>
    transmute(
      Source = source,
      Year = year,
      `Fine age group` = as.character(fine_age_group),
      `Population (million)`,
      `Prevalence per 100,000`,
      `Cases (million)`,
      `Case share (%)`
    )
}

format_oldest_table <- function(profiles) {
  profiles |>
    filter(as.character(fine_age_group) == "70+", year %in% c(2023, 2050)) |>
    mutate(
      `70+ to 55-69 rate ratio` = sprintf("%.1f", older_to_55_69_rate_ratio),
      `Prevalence per 100,000` = sprintf("%.2f", prevalence_per_100k),
      `Cases (million)` = sprintf("%.2f", cases / 1e6),
      `Case share (%)` = sprintf("%.1f", 100 * case_share)
    ) |>
    arrange(year, scenario) |>
    transmute(
      Scenario = as.character(scenario),
      `70+ to 55-69 rate ratio`,
      Source = source,
      Year = year,
      `Fine age group` = as.character(fine_age_group),
      `Prevalence per 100,000`,
      `Cases (million)`,
      `Case share (%)`
    )
}

plot_fine_age_sensitivity <- function(profiles, out_path) {
  neutral <- profiles |>
    filter(scenario == "equal_within_band_1.0") |>
    mutate(
      year = factor(year, levels = c(1990, 2023, 2024, 2050)),
      cases_million = cases / 1e6,
      fine_age_group = factor(as.character(fine_age_group), levels = fine_age_order)
    )

  p <- ggplot(neutral, aes(x = year, y = cases_million, fill = fine_age_group)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = c(
        "<20" = "#4C78A8",
        "20-39" = "#72B7B2",
        "40-54" = "#54A24B",
        "55-69" = "#F58518",
        "70+" = "#B279A2"
      ),
      name = "Fine age group"
    ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(x = NULL, y = "Prevalent cases (million)") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title.position = "top",
      plot.margin = margin(8, 8, 8, 8)
    )

  ggsave(out_path, p, width = 8.5, height = 5.2, dpi = 300)
}

profiles <- build_sensitivity_profiles()

profiles_for_csv <- profiles |>
  mutate(
    fine_age_group = as.character(fine_age_group),
    broad_age_group = as.character(broad_age_group),
    scenario = as.character(scenario)
  )

readr::write_csv(
  profiles_for_csv,
  file.path(appendix_dir, "table_s_fine_age_sensitivity_profiles.csv")
)

oldest_for_csv <- profiles_for_csv |>
  filter(fine_age_group == "70+", year %in% c(2023, 2050)) |>
  arrange(year, factor(scenario, levels = scenarios$scenario))

readr::write_csv(
  oldest_for_csv,
  file.path(appendix_dir, "table_s_fine_age_sensitivity_oldest_scenarios.csv")
)

write_markdown_table(
  format_neutral_table(profiles),
  file.path(appendix_dir, "table_s_fine_age_sensitivity_profiles.md")
)

write_markdown_table(
  format_oldest_table(profiles),
  file.path(appendix_dir, "table_s_fine_age_sensitivity_oldest_scenarios.md")
)

plot_fine_age_sensitivity(
  profiles,
  file.path(appendix_dir, "fig_s_fine_age_sensitivity.png")
)

summary_table <- profiles_for_csv |>
  filter(scenario == "equal_within_band_1.0", fine_age_group == "70+") |>
  transmute(year, neutral_70_plus_case_share_pct = 100 * case_share)

writeLines(
  c(
    "{",
    paste(
      apply(summary_table, 1, function(row) {
        paste0('  "', row[["year"]], '": ', sprintf("%.6f", as.numeric(row[["neutral_70_plus_case_share_pct"]])))
      }),
      collapse = ",\n"
    ),
    "}"
  ),
  file.path(appendix_dir, "fine_age_sensitivity_summary.json"),
  useBytes = TRUE
)

message("Wrote fine-age sensitivity tables and Figure S7.")
