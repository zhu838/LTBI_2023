#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-03-04 19:17:33
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-04 20:35:32
#####################################
library(nih.joinpoint)

# If Joinpoint is identifying too many segments (overfitting), reduce these.
# Start conservative: 1 joinpoint (2 segments). Increase to 2–3 only if needed.
max_joinpoints_rate <- 3
max_joinpoints_number <- 3

build_run_options <- function(max_joinpoints, dependent_variable_type) {
  run_options(
    model = "ln",
    max_joinpoints = max_joinpoints,
    model_selection_method = "permutation test",
    permutation_signif_lvl = 0.05,
    n_permutations = 4499,
    min_obs_end = 2,
    min_obs_between = 2,
    het_error = "constant variance",
    ci_method = "parametric",
    dependent_variable_type = dependent_variable_type,
    n_cores = parallel::detectCores()
  )
}

run_opt_rate <- build_run_options(
  max_joinpoints = max_joinpoints_rate,
  # Input rates are GBD age-standardized prevalence rates, not crude rates.
  dependent_variable_type = "age-adjusted rate"
)

# Age-specific rates are not age-adjusted; they retain the crude-rate setting.
run_opt_age_specific_rate <- build_run_options(
  max_joinpoints = max_joinpoints_rate,
  dependent_variable_type = "crude rate"
)

run_opt_number <- build_run_options(
  max_joinpoints = max_joinpoints_number,
  dependent_variable_type = "count"
)

build_export_opt <- function(year_min, year_max) {
  # Choose interpretable ranges. If data covers 2019, also report pre/post-2019.
  year_min <- as.integer(year_min)
  year_max <- as.integer(year_max)

  if (year_max <= year_min) {
    stop("Invalid year range: year_max must be > year_min")
  }

  pivot_2019 <- if (year_min < 2019 && year_max > 2019) 2019 else NA_integer_

  pivot_1 <- min(year_min + 9L, year_max - 1L)
  pivot_2 <- if (!is.na(pivot_2019)) {
    min(pivot_2019, year_max - 1L)
  } else {
    min(year_min + 19L, year_max - 1L)
  }

  # Ensure increasing pivots
  if (pivot_2 <= pivot_1) {
    pivot_2 <- year_max - 1L
  }

  ranges <- list(
    c(year_min, pivot_1),
    c(pivot_1, pivot_2),
    c(pivot_2, year_max),
    c(year_min, pivot_2),
    c(year_min, year_max)
  )

  # Keep only valid, unique ranges (start < end)
  ranges <- ranges[vapply(ranges, function(r) r[2] > r[1], logical(1))]
  ranges_key <- vapply(ranges, function(r) paste0(r[1], "-", r[2]), character(1))
  ranges <- ranges[!duplicated(ranges_key)]

  base <- export_options(aapc_full_range = TRUE, export_aapc = TRUE)
  extra <- unlist(lapply(seq_along(ranges), function(i) {
    c(
      paste0("AAPC Start Range", i, "=", ranges[[i]][1]),
      paste0("AAPC End Range", i, "=", ranges[[i]][2])
    )
  }))

  paste(c(base, extra), collapse = "\n")
}
