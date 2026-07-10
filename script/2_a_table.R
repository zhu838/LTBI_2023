
# -----------------------------------------------------------------
# Table 1: Prevalence summary across key groups
# Columns: Group, 1990/2023 prevalence rate & number, 1990-2023 AAPC
# -----------------------------------------------------------------

rm(list = ls())

source("./script/config.R")
source("./script/function.R")
source("./script/joinpoint_setting.R")

library(nih.joinpoint)
library(segmented)
library(tidyverse)

# -----------------------------------------------------------------
# Load base data used across groups
# -----------------------------------------------------------------

df_raw_number <- readr::read_csv(
  file.path(database_dir, "global_regional_number.csv"),
  show_col_types = FALSE
)

df_raw_rate <- readr::read_csv(
  file.path(database_dir, "global_regional_rate.csv"),
  show_col_types = FALSE
)

df_number_total <- df_raw_number |>
  # Avoid double-counting: `global_regional_number.csv` contains both age-specific
  # rows and an "All ages" total. For total prevalent cases we use the "All ages"
  # rows only, while retaining summation across any duplicated exports.
  filter(age_name == "All ages") |>
  group_by(measure_name, location_name, sex_name, year, metric_name, age_name) |>
  summarise(
    val = sum(val, na.rm = TRUE),
    lower = sum(lower, na.rm = TRUE),
    upper = sum(upper, na.rm = TRUE),
    .groups = "drop"
  )

df_rate_main <- df_raw_rate |>
  filter(age_name == target_age_global)

locations_sdi <- df_rate_main |>
  distinct(location_name) |>
  filter(stringr::str_detect(location_name, regex("SDI", ignore_case = TRUE))) |>
  pull(location_name) |>
  sort()

preferred_sdi_order <- c(
  "Low SDI",
  "Low-middle SDI",
  "Middle SDI",
  "High-middle SDI",
  "High SDI"
)

sdi_levels <- preferred_sdi_order[preferred_sdi_order %in% locations_sdi]

# The source label for the Americas is "Region of the Americas" rather than
# a name ending in "Region"; use an explicit WHO-region order.
who_regions <- c(
  "African Region",
  "Region of the Americas",
  "Eastern Mediterranean Region",
  "European Region",
  "South-East Asia Region",
  "Western Pacific Region"
)
locations_region <- who_regions[who_regions %in% unique(df_rate_main$location_name)]

age_groups <- intersect(age_groups_for_age_trend, sort(unique(df_raw_number$age_name)))

fit_by_index <- function(d, metric_type) {
  run_opt <- if (metric_type == "Rate") run_opt_rate else run_opt_number
  joinpoint(d, year, val, by = Index, run_opt = run_opt, export_opt = build_export_opt(1990, 2023))
}

fmt_num <- function(x, digits = 2, scale = 1) {
  format(round(x / scale, digits), big.mark = ",", trim = TRUE, nsmall = digits)
}

# -----------------------------------------------------------------
# Read all appendix CSVs once and reuse
# -----------------------------------------------------------------

csv_files <- list.files(appendix_dir, pattern = "\\.csv$", full.names = TRUE)
csv_list <- purrr::map(csv_files, ~ readr::read_csv(.x, show_col_types = FALSE))
names(csv_list) <- basename(csv_files)

measure_target <- "Prevalence"
years_target <- c(1990, 2023)

# ---- AAPC: Global (Both sexes, age-standardized) ----

df_aapc_global <- csv_list[["table_s1_global_trend.csv"]] |>
  filter(var == measure_target,
         start_obs == years_target[1],
         end_obs == years_target[2]) |>
  mutate(Group = "Global") |>
  select(Group, Measure, aapc, aapc_c_i_low, aapc_c_i_high, p_value_label)

# ---- AAPC: Sex (Global) ----

df_aapc_sex_prev <- csv_list[["table_s_sex_group_aapc.csv"]] |>
  filter(Label == measure_target,
         start_obs == years_target[1],
         end_obs == years_target[2]) |>
  select(Group = index, Measure, aapc, aapc_c_i_low, aapc_c_i_high, p_value_label)

# ---- AAPC: Age groups (Global, Both sexes) ----

df_aapc_age_prev <- csv_list[["table_s_age_group_aapc.csv"]] |>
  filter(Label == measure_target,
         start_obs == years_target[1],
         end_obs == years_target[2]) |>
  select(Group = index, Measure, aapc, aapc_c_i_low, aapc_c_i_high, p_value_label)

# ---- AAPC: SDI & Region (Both sexes, age-standardized) ----

df_aapc_sdi_region <- csv_list[["table_s_sdi_region_aapc_prevalence.csv"]] |>
  filter(Label == measure_target,
         start_obs == years_target[1],
         end_obs == years_target[2]) |>
  select(Group = location_name, Measure, aapc, aapc_c_i_low, aapc_c_i_high, p_value_label)

# Combine all AAPC into wide form (Rate / Number columns)

df_aapc_all <- bind_rows(df_aapc_global,
                         df_aapc_sex_prev,
                         df_aapc_age_prev,
                         df_aapc_sdi_region) |>
  unique() |>
  mutate(label = paste0(aapc, " (", aapc_c_i_low, " to ", aapc_c_i_high, ")", dplyr::coalesce(p_value_label, "")),
         Measure = if_else(Measure == "Rate", "rate_aapc_1990_2023", "number_aapc_1990_2023")) |>
  select(Group, Measure, label) |>
  tidyr::pivot_wider(names_from = Measure, values_from = label)

# ---- 1990 & 2023 values for rate and number (with 95% UI) ----

# Global (Both sexes, age-standardized)
df_values_global <- tibble::tibble(
  Group = "Global",
  rate_1990 = df_rate_main |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[1]
    ) |>
    pull(val) |>
    first(),
  rate_1990_lower = df_rate_main |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[1]
    ) |>
    pull(lower) |>
    first(),
  rate_1990_upper = df_rate_main |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[1]
    ) |>
    pull(upper) |>
    first(),
  rate_2023 = df_rate_main |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[2]
    ) |>
    pull(val) |>
    first(),
  rate_2023_lower = df_rate_main |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[2]
    ) |>
    pull(lower) |>
    first(),
  rate_2023_upper = df_rate_main |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[2]
    ) |>
    pull(upper) |>
    first(),
  number_1990 = df_number_total |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[1]
    ) |>
    pull(val) |>
    first(),
  number_1990_lower = df_number_total |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[1]
    ) |>
    pull(lower) |>
    first(),
  number_1990_upper = df_number_total |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[1]
    ) |>
    pull(upper) |>
    first(),
  number_2023 = df_number_total |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[2]
    ) |>
    pull(val) |>
    first(),
  number_2023_lower = df_number_total |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[2]
    ) |>
    pull(lower) |>
    first(),
  number_2023_upper = df_number_total |>
    filter(
      location_name == target_location_global,
      sex_name == "Both",
      measure_name == measure_target,
      year == years_target[2]
    ) |>
    pull(upper) |>
    first()
)

# Sex (Global)
df_rate_sex_vals <- df_rate_main |>
  filter(
    location_name == target_location_global,
    sex_name %in% c("Male", "Female"),
    measure_name == measure_target,
    year %in% years_target
  ) |>
  transmute(
    Group = sex_name,
    year,
    rate = val,
    rate_lower = lower,
    rate_upper = upper
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(rate, rate_lower, rate_upper),
    names_sep = "_"
  ) |>
  transmute(
    Group,
    rate_1990 = rate_1990,
    rate_1990_lower = rate_lower_1990,
    rate_1990_upper = rate_upper_1990,
    rate_2023 = rate_2023,
    rate_2023_lower = rate_lower_2023,
    rate_2023_upper = rate_upper_2023
  )

df_number_sex_vals <- df_number_total |>
  filter(
    location_name == target_location_global,
    sex_name %in% c("Male", "Female"),
    measure_name == measure_target,
    year %in% years_target
  ) |>
  transmute(
    Group = sex_name,
    year,
    number = val,
    number_lower = lower,
    number_upper = upper
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(number, number_lower, number_upper),
    names_sep = "_"
  ) |>
  transmute(
    Group,
    number_1990 = number_1990,
    number_1990_lower = number_lower_1990,
    number_1990_upper = number_upper_1990,
    number_2023 = number_2023,
    number_2023_lower = number_lower_2023,
    number_2023_upper = number_upper_2023
  )

df_values_sex <- left_join(df_rate_sex_vals, df_number_sex_vals, by = "Group")

# Age groups (Global, Both sexes)
df_rate_age_vals <- df_raw_rate |>
  filter(
    location_name == target_location_global,
    sex_name == "Both",
    measure_name == measure_target,
    age_name %in% age_groups,
    year %in% years_target
  ) |>
  transmute(
    Group = age_name,
    year,
    rate = val,
    rate_lower = lower,
    rate_upper = upper
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(rate, rate_lower, rate_upper),
    names_sep = "_"
  ) |>
  transmute(
    Group,
    rate_1990 = rate_1990,
    rate_1990_lower = rate_lower_1990,
    rate_1990_upper = rate_upper_1990,
    rate_2023 = rate_2023,
    rate_2023_lower = rate_lower_2023,
    rate_2023_upper = rate_upper_2023
  )

df_number_age_vals <- df_raw_number |>
  filter(
    location_name == target_location_global,
    sex_name == "Both",
    measure_name == measure_target,
    age_name %in% age_groups,
    year %in% years_target
  ) |>
  transmute(
    Group = age_name,
    year,
    number = val,
    number_lower = lower,
    number_upper = upper
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(number, number_lower, number_upper),
    names_sep = "_"
  ) |>
  transmute(
    Group,
    number_1990 = number_1990,
    number_1990_lower = number_lower_1990,
    number_1990_upper = number_upper_1990,
    number_2023 = number_2023,
    number_2023_lower = number_lower_2023,
    number_2023_upper = number_upper_2023
  )

df_values_age <- left_join(df_rate_age_vals, df_number_age_vals, by = "Group")

# SDI & Region (Both sexes, age-standardized)
locs_all <- c(locations_sdi, locations_region)

df_rate_sdi_region_vals <- df_rate_main |>
  filter(
    sex_name == "Both",
    measure_name == measure_target,
    location_name %in% locs_all,
    year %in% years_target
  ) |>
  transmute(
    Group = location_name,
    year,
    rate = val,
    rate_lower = lower,
    rate_upper = upper
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(rate, rate_lower, rate_upper),
    names_sep = "_"
  ) |>
  transmute(
    Group,
    rate_1990 = rate_1990,
    rate_1990_lower = rate_lower_1990,
    rate_1990_upper = rate_upper_1990,
    rate_2023 = rate_2023,
    rate_2023_lower = rate_lower_2023,
    rate_2023_upper = rate_upper_2023
  )

df_number_sdi_region_vals <- df_number_total |>
  filter(
    sex_name == "Both",
    measure_name == measure_target,
    location_name %in% locs_all,
    year %in% years_target
  ) |>
  transmute(
    Group = location_name,
    year,
    number = val,
    number_lower = lower,
    number_upper = upper
  ) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(number, number_lower, number_upper),
    names_sep = "_"
  ) |>
  transmute(
    Group,
    number_1990 = number_1990,
    number_1990_lower = number_lower_1990,
    number_1990_upper = number_upper_1990,
    number_2023 = number_2023,
    number_2023_lower = number_lower_2023,
    number_2023_upper = number_upper_2023
  )

df_values_sdi_region <- left_join(df_rate_sdi_region_vals, df_number_sdi_region_vals, by = "Group")

# Combine all groups and attach AAPC
df_values_all <- bind_rows(
  df_values_global,
  df_values_sex,
  df_values_age,
  df_values_sdi_region
)

df_table_1 <- df_values_all |>
  left_join(df_aapc_all, by = "Group") |>
  mutate(
    `1990_prevalence_rate` = sprintf(
      "%s (%s to %s)",
      fmt_num(rate_1990, digits = 2, scale = 1),
      fmt_num(rate_1990_lower, digits = 2, scale = 1),
      fmt_num(rate_1990_upper, digits = 2, scale = 1)
    ),
    `2023_prevalence_rate` = sprintf(
      "%s (%s to %s)",
      fmt_num(rate_2023, digits = 2, scale = 1),
      fmt_num(rate_2023_lower, digits = 2, scale = 1),
      fmt_num(rate_2023_upper, digits = 2, scale = 1)
    ),
    `1990_prevalence_number` = sprintf(
      "%s (%s to %s)",
      fmt_num(number_1990, digits = 2, scale = 1e6),
      fmt_num(number_1990_lower, digits = 2, scale = 1e6),
      fmt_num(number_1990_upper, digits = 2, scale = 1e6)
    ),
    `2023_prevalence_number` = sprintf(
      "%s (%s to %s)",
      fmt_num(number_2023, digits = 2, scale = 1e6),
      fmt_num(number_2023_lower, digits = 2, scale = 1e6),
      fmt_num(number_2023_upper, digits = 2, scale = 1e6)
    ),
    `1990_2023_rate_AAPC` = rate_aapc_1990_2023,
    `1990_2023_number_AAPC` = number_aapc_1990_2023
  ) |>
  select(
    Group,
    `1990_prevalence_rate`,
    `2023_prevalence_rate`,
    `1990_2023_rate_AAPC`,
    `1990_prevalence_number`,
    `2023_prevalence_number`,
    `1990_2023_number_AAPC`
  )

# Order rows: Global -> Sex -> Age -> SDI -> Region
group_levels <- c(
  "Global",
  "Female", "Male",
  sort(age_groups),
  sdi_levels,
  locations_region
)

df_table_1 <- df_table_1 |>
  mutate(Group = factor(Group, levels = group_levels)) |>
  arrange(Group)

readr::write_csv(df_table_1, file.path(appendix_dir, "table_1_prevalence_summary.csv"))

# -----------------------------------------------------------------
# Markdown version with nicer column names and group separators
# -----------------------------------------------------------------

df_global <- df_table_1 |> dplyr::filter(Group == "Global")
df_sex    <- df_table_1 |> dplyr::filter(Group %in% c("Female", "Male"))
df_age    <- df_table_1 |> dplyr::filter(Group %in% age_groups)
df_sdi    <- df_table_1 |> dplyr::filter(Group %in% sdi_levels)
df_region <- df_table_1 |> dplyr::filter(Group %in% locations_region)

make_block <- function(label, df_section) {
  if (nrow(df_section) == 0) return(df_section)
  cols <- names(df_section)
  header <- tibble::as_tibble(setNames(as.list(rep("", length(cols))), cols))
  header$Group[1] <- label
  dplyr::bind_rows(header, df_section)
}

df_table_1_md <- dplyr::bind_rows(
  df_global,
  make_block("Sex group", df_sex),
  make_block("Age group", df_age),
  make_block("SDI group", df_sdi),
  make_block("Region group", df_region)
)

markdown_table <- knitr::kable(
  df_table_1_md,
  format = "markdown",
  col.names = c(
    "Group",
    "1990 Prevalence rate (per 100,000, 95% UI)",
    "2023 Prevalence rate (per 100,000, 95% UI)",
    "1990–2023 rate AAPC (95% CI)",
    "1990 Prevalent cases (million, 95% UI)",
    "2023 Prevalent cases (million, 95% UI)",
    "1990–2023 number AAPC (95% CI)"
  ),
  align = "lcccccc",
  escape = FALSE
)

write(markdown_table, file.path(outcome_dir, "table_1_prevalence_summary.md"))

message("Done: regional summary outputs and Table 1 -> ", outcome_dir)
