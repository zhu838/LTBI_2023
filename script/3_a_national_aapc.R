#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-03-05 17:20:58
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-05 17:21:05
#####################################
# National-level AAPC for age-standardized prevalence rate (2010–2023)

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

source("./script/config.R")
source("./script/function.R")
source("./script/joinpoint_setting.R")

# -------------------------------------------------------------------------
# Prepared data: age-standardized national prevalence rate
# -------------------------------------------------------------------------

df_map_iso <- readr::read_csv(iso_code_file, show_col_types = FALSE)

# Prefer per-year national files if available (written by 1_data_prepare.py)
national_dir <- file.path(database_dir, "national_by_year")
csvs <- list.files(national_dir, pattern = "^national_.*\\.csv$", full.names = TRUE)
df_raw_rate <- purrr::map_dfr(csvs, ~ readr::read_csv(.x, show_col_types = FALSE))

measures_available <- intersect(target_measures, sort(unique(df_raw_rate$measure_name)))
if (length(measures_available) == 0) {
  stop("No target measures found in prepared data.")
}

# Age-standardized prevalence rate for all countries with ISO mapping
df_all_rate <- df_raw_rate |>
  filter(age_name == target_age_global,
         sex_name == 'Both',
         location %in% df_map_iso$location_id,
         year %in% 2010:2023,
         measure_name %in% measures_available) |>
  select(location, location_name, measure_name, year, val, lower, upper)

year_min_data <- min(df_all_rate$year, na.rm = TRUE)
year_max_data <- max(df_all_rate$year, na.rm = TRUE)

# AAPC window: 2010–last available year (expected: 2023)
year_start <- max(2010, year_min_data)
year_end <- year_max_data

if (year_end < year_start + 2) {
  stop("Not enough years after 2010 to compute AAPC.")
}

export_opt_new <- build_export_opt(year_start, year_end)

# -------------------------------------------------------------------------
# Joinpoint by country (rate only) and extract AAPC 2010–year_end
# -------------------------------------------------------------------------

fit_by_country_rate <- function(d) {
  joinpoint(d, year, val, by = location_name, run_opt = run_opt_rate, export_opt = export_opt_new)
}

aapc_rows <- list()
for (m in measures_available) {
  model_rate <- fit_by_country_rate(filter(df_all_rate, measure_name == m))
  aapc_rows[[m]] <- get_aapc(model_rate) |>
    mutate(Label = m, Measure = "Rate")
}

if (length(aapc_rows) == 0) {
  stop("No AAPC models were fitted.")
}

df_aapc_all <- bind_rows(aapc_rows) |>
  mutate(aapc = as.numeric(aapc),
         aapc_c_i_low = as.numeric(aapc_c_i_low),
         aapc_c_i_high = as.numeric(aapc_c_i_high)) |> 
  filter(aapc_index == 'Full Range')

out_file <- file.path(outcome_dir, sprintf("national_aapc_%d_%d.csv", year_start, year_end))
readr::write_csv(df_aapc_all, out_file)

message("Saved: ", out_file)
