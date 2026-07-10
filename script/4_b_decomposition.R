#####################################
## @Description: Decomposition of LTBI case change
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

library(dplyr)
library(purrr)
library(readr)
library(tidyr)

source("./script/config.R")
source("./script/function.R")

dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

rate_file <- file.path(database_dir, "global_regional_rate.csv")
number_file <- file.path(database_dir, "global_regional_number.csv")
forecast_file_candidates <- c(
  # The manuscript and Figure 3 report the negative-binomial APC model.
  # Use the same output for decomposition so that case-count changes reconcile.
  file.path(outcome_dir, "forecast", "predictions_rstan_negbin_Global_2024_2050.csv")
)
forecast_file <- forecast_file_candidates[file.exists(forecast_file_candidates)][1]

if (is.na(forecast_file) || !file.exists(forecast_file)) {
  stop("Missing negative-binomial forecast file required for projected decomposition.")
}

target_age_bins <- c("<20", "20-54", "55+")

compose_cases <- function(total_population, age_share, prevalence) {
  total_population * sum(age_share * prevalence)
}

decompose_three_factor <- function(base_profile, target_profile) {
  factors_base <- list(
    N = base_profile$total_population[[1]],
    s = base_profile$age_share,
    r = base_profile$prevalence
  )
  factors_target <- list(
    N = target_profile$total_population[[1]],
    s = target_profile$age_share,
    r = target_profile$prevalence
  )

  permutations <- list(
    c("N", "s", "r"),
    c("N", "r", "s"),
    c("s", "N", "r"),
    c("s", "r", "N"),
    c("r", "N", "s"),
    c("r", "s", "N")
  )

  contrib <- c(N = 0, s = 0, r = 0)

  for (perm in permutations) {
    current <- factors_base
    current_value <- compose_cases(current$N, current$s, current$r)

    for (factor_name in perm) {
      updated <- current
      updated[[factor_name]] <- factors_target[[factor_name]]
      updated_value <- compose_cases(updated$N, updated$s, updated$r)
      contrib[[factor_name]] <- contrib[[factor_name]] + (updated_value - current_value)
      current <- updated
      current_value <- updated_value
    }
  }

  contrib <- contrib / length(permutations)

  tribble(
    ~component, ~contribution_cases,
    "Population growth", contrib[["N"]],
    "Population ageing", contrib[["s"]],
    "Epidemiologic change", contrib[["r"]]
  )
}

df_rate_raw <- readr::read_csv(rate_file, show_col_types = FALSE)
df_number_raw <- readr::read_csv(number_file, show_col_types = FALSE)

observed_rates <- df_rate_raw |>
  filter(
    measure_name == "Prevalence",
    metric_name == "Rate",
    location_name == target_location_global,
    sex_name == target_sex_global,
    age_name %in% c("<20 years", "20-54 years", "55+ years")
  ) |>
  transmute(
    year = as.integer(year),
    age_group = case_when(
      age_name == "<20 years" ~ "<20",
      age_name == "20-54 years" ~ "20-54",
      age_name == "55+ years" ~ "55+",
      TRUE ~ age_name
    ),
    prevalence_per_100k = as.numeric(val),
    prevalence = as.numeric(val) / 1e5
  )

observed_numbers <- df_number_raw |>
  filter(
    measure_name == "Prevalence",
    metric_name == "Number",
    location_name == target_location_global,
    sex_name == target_sex_global,
    age_name %in% c("<20 years", "20-54 years", "55+ years")
  ) |>
  transmute(
    year = as.integer(year),
    age_group = case_when(
      age_name == "<20 years" ~ "<20",
      age_name == "20-54 years" ~ "20-54",
      age_name == "55+ years" ~ "55+",
      TRUE ~ age_name
    ),
    cases = as.numeric(val)
  )

observed_profiles <- observed_rates |>
  inner_join(observed_numbers, by = c("year", "age_group")) |>
  mutate(population = cases / prevalence)

build_profile <- function(data, year_value, prevalence_col) {
  year_data <- data |>
    filter(year == year_value, age_group %in% target_age_bins) |>
    arrange(factor(age_group, levels = target_age_bins))

  if (!"cases" %in% names(year_data)) {
    year_data <- year_data |>
      mutate(cases = .data[[prevalence_col]] * population)
  }

  total_population <- sum(year_data$population, na.rm = TRUE)

  year_data |>
    mutate(
      total_population = total_population,
      age_share = population / total_population,
      prevalence = .data[[prevalence_col]],
      cases_million = cases / 1e6
    )
}

run_decomposition <- function(start_profile, end_profile, label, start_year, end_year) {
  contributions <- decompose_three_factor(start_profile, end_profile)
  total_change_cases <- sum(end_profile$cases) - sum(start_profile$cases)

  contributions |>
    mutate(
      comparison = label,
      start_year = start_year,
      end_year = end_year,
      total_change_cases = total_change_cases,
      contribution_million = contribution_cases / 1e6,
      contribution_pct = contribution_cases / total_change_cases
    ) |>
    select(comparison, start_year, end_year, component, contribution_cases,
           contribution_million, contribution_pct, total_change_cases)
}

profile_1990 <- build_profile(observed_profiles, 1990, "prevalence")
profile_2023 <- build_profile(observed_profiles, 2023, "prevalence")

decomposition_results <- list(
  run_decomposition(profile_1990, profile_2023, "Observed 1990 to 2023", 1990, 2023)
)

profile_exports <- list(
  profile_1990 |> mutate(comparison = "Observed 1990 to 2023", profile_year = 1990),
  profile_2023 |> mutate(comparison = "Observed 1990 to 2023", profile_year = 2023)
)

if (!is.na(forecast_file) && file.exists(forecast_file)) {
  forecast_profiles <- readr::read_csv(forecast_file, show_col_types = FALSE) |>
    transmute(
      year = as.integer(year),
      age_group = age_group,
      population = as.numeric(population),
      prevalence_per_100k = as.numeric(pred_prev_mean),
      prevalence = as.numeric(pred_prev_mean) / 1e5
    )

  profile_2024 <- build_profile(forecast_profiles, 2024, "prevalence")
  profile_2050 <- build_profile(forecast_profiles, 2050, "prevalence")

  decomposition_results[[length(decomposition_results) + 1]] <-
    run_decomposition(profile_2024, profile_2050, "Projected 2024 to 2050", 2024, 2050)

  profile_exports[[length(profile_exports) + 1]] <-
    profile_2024 |> mutate(comparison = "Projected 2024 to 2050", profile_year = 2024)
  profile_exports[[length(profile_exports) + 1]] <-
    profile_2050 |> mutate(comparison = "Projected 2024 to 2050", profile_year = 2050)
}

df_decomposition <- bind_rows(decomposition_results)
df_profiles <- bind_rows(profile_exports) |>
  select(comparison, profile_year, age_group, total_population, population, age_share,
         prevalence_per_100k, prevalence, cases, cases_million)

decomposition_markdown <- df_decomposition |>
  transmute(
    `Comparison period` = comparison,
    Component = component,
    `Contribution (million cases)` = sprintf("%.2f", contribution_million),
    `Contribution (% of total change)` = sprintf("%.1f%%", 100 * contribution_pct),
    `Total change (million cases)` = sprintf("%.2f", total_change_cases / 1e6)
  )

profiles_markdown <- df_profiles |>
  transmute(
    `Comparison period` = comparison,
    Year = profile_year,
    `Age group` = age_group,
    `Total population` = sprintf("%.0f", total_population),
    `Population in age group` = sprintf("%.0f", population),
    `Age share` = sprintf("%.4f", age_share),
    `Prevalence per 100,000` = sprintf("%.2f", prevalence_per_100k),
    `Cases (million)` = sprintf("%.2f", cases_million)
  )

readr::write_csv(
  df_profiles,
  file.path(appendix_dir, "table_s_decomposition_profiles.csv")
)
readr::write_csv(
  df_decomposition,
  file.path(appendix_dir, "table_s_decomposition_components.csv")
)

write_markdown_table(
  decomposition_markdown,
  file.path(appendix_dir, "table_s_decomposition_components.md")
)

write_markdown_table(
  profiles_markdown,
  file.path(appendix_dir, "table_s_decomposition_profiles.md")
)

message("Done: decomposition outputs written to ", appendix_dir)
