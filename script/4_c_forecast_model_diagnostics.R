#####################################
## @Description: Forecast APC model diagnostics and validation
## @Date: 2026-05-29
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
library(rstan)
library(loo)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))

source("./script/config.R")
source("./script/function.R")

out_dir <- file.path(outcome_dir, "forecast", "diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

target_location <- target_location_global
target_age_bins <- c("<20", "20-54", "55+")
observed_years <- 2010:2023
projection_years <- 2024:2050
all_model_years <- 2010:2050
pop_file <- "./Data/unpopulation_dataportal_20260305182132.csv"
gr_rate_file <- file.path(database_dir, "global_regional_rate.csv")

build_model_data <- function() {
  pop_raw <- readr::read_csv(pop_file, show_col_types = FALSE)
  df_iso_pop <- readr::read_csv(iso_code_file, show_col_types = FALSE)
  country_iso3 <- toupper(as.character(df_iso_pop$ISO3))

  pop_age3 <- pop_raw |>
    mutate(Iso3 = toupper(as.character(Iso3))) |>
    filter(Iso3 %in% country_iso3) |>
    transmute(
      year = as.integer(Time),
      age_start = as.integer(AgeStart),
      population = as.numeric(Value),
      age_group = case_when(
        age_start < 20 ~ "<20",
        age_start >= 20 & age_start < 55 ~ "20-54",
        age_start >= 55 ~ "55+",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(age_group), year %in% all_model_years) |>
    group_by(year, age_group) |>
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop") |>
    mutate(location = target_location) |>
    select(location, year, age_group, population)

  df_raw_rate <- readr::read_csv(gr_rate_file, show_col_types = FALSE)
  df_rate <- df_raw_rate |>
    filter(
      measure_name %in% target_measures,
      metric_name == "Rate",
      location_name == target_location_global,
      sex_name == target_sex_global
    ) |>
    transmute(
      location = location_name,
      year = as.integer(year),
      age_group = case_when(
        age_name == "<20 years" ~ "<20",
        age_name == "20-54 years" ~ "20-54",
        age_name == "55+ years" ~ "55+",
        TRUE ~ age_name
      ),
      prevalence_rate = as.numeric(val)
    ) |>
    filter(year %in% observed_years, age_group %in% target_age_bins)

  if (max(df_rate$prevalence_rate, na.rm = TRUE) > 10) {
    df_rate <- df_rate |> mutate(prevalence_rate = prevalence_rate / 1e5)
  }

  hist_data <- df_rate |>
    left_join(pop_age3 |> filter(year %in% observed_years), by = c("location", "year", "age_group")) |>
    mutate(cases = round(prevalence_rate * population))

  future_grid <- tidyr::expand_grid(
    location = target_location,
    year = projection_years,
    age_group = target_age_bins
  ) |>
    left_join(pop_age3 |> filter(year %in% projection_years), by = c("location", "year", "age_group")) |>
    mutate(prevalence_rate = NA_real_, cases = NA_integer_)

  all_data <- bind_rows(hist_data, future_grid) |>
    arrange(year, factor(age_group, levels = target_age_bins)) |>
    mutate(
      age_index = as.integer(factor(age_group, levels = target_age_bins)),
      period_index = as.integer(factor(year, levels = sort(unique(year))))
    )

  all_data <- all_data |>
    mutate(
      cohort_raw = period_index - age_index + length(target_age_bins),
      cohort_index = as.integer(factor(cohort_raw))
    )

  if (any(is.na(all_data$population)) || any(all_data$population <= 0)) {
    stop("Population denominators are missing or non-positive.")
  }

  y_vec <- as.integer(all_data$cases)
  y_vec[is.na(y_vec)] <- 0L
  is_obs_vec <- as.integer(!is.na(all_data$cases))

  stan_data <- list(
    N = nrow(all_data),
    y = y_vec,
    is_obs = is_obs_vec,
    log_E = log(all_data$population),
    A = length(target_age_bins),
    T = length(unique(all_data$period_index)),
    K = length(unique(all_data$cohort_index)),
    age_index = as.integer(all_data$age_index),
    period_index = as.integer(all_data$period_index),
    cohort_index = as.integer(all_data$cohort_index)
  )

  list(all_data = all_data, stan_data = stan_data, obs_idx = which(is_obs_vec == 1L))
}

stan_code <- "
data {
  int<lower=1> N;
  int<lower=0> y[N];
  int<lower=0,upper=1> is_obs[N];
  int<lower=0,upper=1> use_negbin;
  vector[N] log_E;
  int<lower=1> A;
  int<lower=1> T;
  int<lower=1> K;
  int<lower=1,upper=A> age_index[N];
  int<lower=1,upper=T> period_index[N];
  int<lower=1,upper=K> cohort_index[N];
}
parameters {
  real alpha;
  vector[A] age_raw;
  vector[T-1] d1_period;
  vector[K-1] d1_cohort;
  real<lower=0> sigma_age;
  real<lower=0> sigma_period;
  real<lower=0> sigma_cohort;
  real log_omega;
}
transformed parameters {
  vector[A] age_effect;
  vector[T] period_effect;
  vector[K] cohort_effect;
  vector[N] eta;
  real omega;
  real phi;

  omega = exp(log_omega);
  phi = inv(omega);
  age_effect = (age_raw - mean(age_raw)) * sigma_age;

  period_effect[1] = 0;
  period_effect[2] = d1_period[1] * sigma_period;
  for (t in 3:T)
    period_effect[t] = 2 * period_effect[t-1] - period_effect[t-2]
                       + (d1_period[t-1] - d1_period[t-2]) * sigma_period;

  cohort_effect[1] = 0;
  cohort_effect[2] = d1_cohort[1] * sigma_cohort;
  for (k in 3:K)
    cohort_effect[k] = 2 * cohort_effect[k-1] - cohort_effect[k-2]
                       + (d1_cohort[k-1] - d1_cohort[k-2]) * sigma_cohort;

  for (n in 1:N) {
    eta[n] = alpha
             + age_effect[age_index[n]]
             + period_effect[period_index[n]]
             + cohort_effect[cohort_index[n]]
             + log_E[n];
  }
}
model {
  alpha ~ normal(0, 5);
  age_raw ~ normal(0, 1);
  d1_period ~ normal(0, 1);
  d1_cohort ~ normal(0, 1);
  sigma_age ~ exponential(1);
  sigma_period ~ exponential(1);
  sigma_cohort ~ exponential(1);
  log_omega ~ normal(log(1e-4), 1);

  for (n in 1:N) {
    if (is_obs[n] == 1) {
      if (use_negbin == 1)
        y[n] ~ neg_binomial_2_log(eta[n], phi);
      else
        y[n] ~ poisson_log(eta[n]);
    }
  }
}
generated quantities {
  vector[N] mu;
  vector[N] log_lik;
  for (n in 1:N) {
    mu[n] = exp(eta[n]);
    if (is_obs[n] == 1) {
      if (use_negbin == 1)
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta[n], phi);
      else
        log_lik[n] = poisson_log_lpmf(y[n] | eta[n]);
    } else {
      log_lik[n] = 0;
    }
  }
}
"

fit_model <- function(sm, stan_data, use_negbin, seed) {
  stan_data$use_negbin <- as.integer(use_negbin)
  model_label <- ifelse(use_negbin == 1, "negative binomial", "Poisson")
  adapt_delta <- ifelse(use_negbin == 1, 0.999, 0.95)
  max_treedepth <- ifelse(use_negbin == 1, 15, 12)
  message("Fitting ", model_label, " APC model")
  fit <- rstan::sampling(
    object = sm,
    data = stan_data,
    chains = 4,
    iter = 3000,
    warmup = 1500,
    thin = 1,
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = 250
  )
  fit
}

compile_model <- function() {
  sm <- rstan::stan_model(model_code = stan_code)
  sm
}

summarise_convergence <- function(fit, model_name, configured_max_treedepth) {
  base_pars <- c(
    "alpha", "age_raw", "d1_period", "d1_cohort",
    "sigma_age", "sigma_period", "sigma_cohort"
  )
  pars <- if (model_name == "Negative binomial") c(base_pars, "omega", "phi") else base_pars
  summ <- rstan::summary(fit, pars = pars)$summary
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
  max_treedepth_observed <- max(vapply(sampler_params, function(x) max(x[, "treedepth__"]), numeric(1)))
  max_treedepth_hits <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= configured_max_treedepth), numeric(1)))

  tibble(
    model = model_name,
    diagnostic = c("Maximum R-hat", "Minimum effective sample size", "Divergent transitions", "Maximum treedepth observed", "Configured maximum-treedepth transitions"),
    value = c(
      max(summ[, "Rhat"], na.rm = TRUE),
      min(summ[, "n_eff"], na.rm = TRUE),
      divergences,
      max_treedepth_observed,
      max_treedepth_hits
    )
  )
}

summarise_key_parameters <- function(fit, model_name) {
  pars <- c("alpha", "sigma_age", "sigma_period", "sigma_cohort")
  if (model_name == "Negative binomial") pars <- c(pars, "omega", "phi")
  summ <- as.data.frame(rstan::summary(fit, pars = pars)$summary)
  summ$parameter <- rownames(summ)
  summ |>
    transmute(
      model = model_name,
      parameter,
      mean = mean,
      sd = sd,
      `2.5%` = `2.5%`,
      `97.5%` = `97.5%`,
      n_eff = n_eff,
      Rhat = Rhat
    )
}

run_loo <- function(fit, obs_idx) {
  extracted <- rstan::extract(fit, pars = "log_lik")
  ll <- extracted$log_lik[, obs_idx, drop = FALSE]
  loo::loo(ll)
}

loo_summary_table <- function(loo_list) {
  out <- bind_rows(
    as_tibble(loo_list$poisson$estimates, rownames = "estimate") |>
      filter(estimate %in% c("elpd_loo", "p_loo", "looic")) |>
      mutate(model = "Poisson likelihood"),
    as_tibble(loo_list$negbin$estimates, rownames = "estimate") |>
      filter(estimate %in% c("elpd_loo", "p_loo", "looic")) |>
      mutate(model = "Negative binomial likelihood")
  ) |>
    select(model, estimate, Estimate, SE)

  pareto <- tibble(
    model = c("Poisson likelihood", "Negative binomial likelihood"),
    estimate = "max_pareto_k",
    Estimate = c(max(loo_list$poisson$diagnostics$pareto_k), max(loo_list$negbin$diagnostics$pareto_k)),
    SE = NA_real_
  )
  bind_rows(out, pareto)
}

ppc_summary <- function(fit, obs_idx, y_obs, model_name, n_draws_keep = 2000) {
  extracted <- rstan::extract(fit, pars = c("mu", "phi"))
  mu <- extracted$mu[, obs_idx, drop = FALSE]
  if (nrow(mu) > n_draws_keep) {
    set.seed(20260529)
    keep <- sort(sample(seq_len(nrow(mu)), n_draws_keep))
    mu <- mu[keep, , drop = FALSE]
    phi <- extracted$phi[keep]
  } else {
    phi <- extracted$phi
  }

  y_obs_mat <- matrix(y_obs, nrow = nrow(mu), ncol = length(y_obs), byrow = TRUE)

  if (model_name == "Poisson") {
    y_rep <- matrix(
      stats::rpois(length(mu), lambda = pmax(as.vector(mu), 1e-9)),
      nrow = nrow(mu),
      ncol = ncol(mu)
    )
  } else {
    y_rep <- matrix(NA_real_, nrow = nrow(mu), ncol = ncol(mu))
    for (i in seq_len(nrow(mu))) {
      y_rep[i, ] <- stats::rnbinom(ncol(mu), size = phi[i], mu = pmax(mu[i, ], 1e-9))
    }
  }

  pearson_obs <- rowSums((y_obs_mat - mu)^2 / pmax(mu, 1))
  pearson_rep <- rowSums((y_rep - mu)^2 / pmax(mu, 1))
  max_abs_obs <- apply(abs((y_obs_mat - mu) / sqrt(pmax(mu, 1))), 1, max)
  max_abs_rep <- apply(abs((y_rep - mu) / sqrt(pmax(mu, 1))), 1, max)
  total_obs <- sum(y_obs)
  total_rep <- rowSums(y_rep)

  row_lwr <- apply(y_rep, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  row_upr <- apply(y_rep, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  coverage <- mean(y_obs >= row_lwr & y_obs <= row_upr)
  mape <- mean(abs(colMeans(mu) - y_obs) / pmax(y_obs, 1)) * 100

  tibble(
    model = model_name,
    check = c("Pearson discrepancy", "Maximum absolute standardized residual", "Total observed count", "Row-level 95% predictive coverage", "Mean absolute percentage error"),
    observed = c(mean(pearson_obs), mean(max_abs_obs), total_obs, coverage, mape),
    posterior_predictive_mean = c(mean(pearson_rep), mean(max_abs_rep), mean(total_rep), NA_real_, NA_real_),
    posterior_predictive_lwr = c(
      quantile(pearson_rep, 0.025),
      quantile(max_abs_rep, 0.025),
      quantile(total_rep, 0.025),
      NA_real_,
      NA_real_
    ),
    posterior_predictive_upr = c(
      quantile(pearson_rep, 0.975),
      quantile(max_abs_rep, 0.975),
      quantile(total_rep, 0.975),
      NA_real_,
      NA_real_
    ),
    posterior_predictive_p = c(
      mean(pearson_rep >= pearson_obs),
      mean(max_abs_rep >= max_abs_obs),
      mean(total_rep >= total_obs),
      NA_real_,
      NA_real_
    )
  )
}

format_numeric_table <- function(df) {
  df |>
    mutate(across(where(is.numeric), ~ ifelse(is.na(.x), NA_character_, sprintf("%.3f", .x))))
}

write_negbin_predictions <- function(fit, all_data) {
  post_mu <- rstan::extract(fit, pars = "mu")$mu
  pred <- all_data |>
    mutate(
      fitted_mean = colMeans(post_mu),
      fitted_lwr = apply(post_mu, 2, stats::quantile, probs = 0.025),
      fitted_upr = apply(post_mu, 2, stats::quantile, probs = 0.975),
      pred_prev_mean = fitted_mean / population * 1e5,
      pred_prev_lwr = fitted_lwr / population * 1e5,
      pred_prev_upr = fitted_upr / population * 1e5
    )

  pred_file <- file.path(outcome_dir, "forecast", "predictions_rstan_negbin_Global_2024_2050.csv")
  readr::write_csv(
    pred |>
      filter(year %in% projection_years) |>
      select(location, year, age_group, population,
             pred_prev_mean, pred_prev_lwr, pred_prev_upr,
             fitted_mean, fitted_lwr, fitted_upr),
    pred_file
  )

  all_pred_file <- file.path(out_dir, "predictions_rstan_negbin_Global_2010_2050_all.csv")
  readr::write_csv(pred, all_pred_file)

  summary_rows <- pred |>
    filter(year %in% c(2024, 2050)) |>
    group_by(year) |>
    summarise(
      total_population = sum(population),
      cases_mean = sum(fitted_mean),
      cases_lwr = sum(fitted_lwr),
      cases_upr = sum(fitted_upr),
      rate_mean = cases_mean / total_population * 1e5,
      rate_lwr = cases_lwr / total_population * 1e5,
      rate_upr = cases_upr / total_population * 1e5,
      .groups = "drop"
    )

  summary_file <- file.path(out_dir, "forecast_summary_negbin_2024_2050.csv")
  readr::write_csv(summary_rows, summary_file)

  summary_md <- summary_rows |>
    transmute(
      Year = year,
      `Population (billion)` = sprintf("%.3f", total_population / 1e9),
      `Prevalent cases (billion), mean` = sprintf("%.2f", cases_mean / 1e9),
      `Prevalent cases (billion), 95% CrI` = sprintf("%.2f to %.2f", cases_lwr / 1e9, cases_upr / 1e9),
      `Crude prevalence per 100,000, mean` = sprintf("%.2f", rate_mean),
      `Crude prevalence per 100,000, 95% CrI` = sprintf("%.2f to %.2f", rate_lwr, rate_upr)
    )
  write_markdown_table(summary_md, file.path(appendix_dir, "table_s_forecast_negbin_summary.md"))
  readr::write_csv(summary_md, file.path(appendix_dir, "table_s_forecast_negbin_summary.csv"))

  list(pred = pred, summary = summary_rows)
}

input <- build_model_data()
sm <- compile_model()
poisson_fit_path <- file.path(out_dir, "forecast_apc_poisson_fit.rds")
negbin_fit_path <- file.path(out_dir, "forecast_apc_negbin_fit.rds")

if (Sys.getenv("CACHE_POISSON_FIT", "0") == "1" && file.exists(poisson_fit_path)) {
  message("Loading cached Poisson APC fit from ", poisson_fit_path)
  fit_poisson <- readRDS(poisson_fit_path)
} else {
  fit_poisson <- fit_model(sm, input$stan_data, use_negbin = 0, seed = 20260529)
}
fit_negbin <- fit_model(sm, input$stan_data, use_negbin = 1, seed = 20260530)
saveRDS(fit_poisson, poisson_fit_path)
saveRDS(fit_negbin, negbin_fit_path)

conv <- bind_rows(
  summarise_convergence(fit_poisson, "Poisson", configured_max_treedepth = 12),
  summarise_convergence(fit_negbin, "Negative binomial", configured_max_treedepth = 15)
)
key_pars <- bind_rows(
  summarise_key_parameters(fit_poisson, "Poisson"),
  summarise_key_parameters(fit_negbin, "Negative binomial")
)
loo_list <- list(
  poisson = run_loo(fit_poisson, input$obs_idx),
  negbin = run_loo(fit_negbin, input$obs_idx)
)
loo_table <- loo_summary_table(loo_list)

y_obs <- input$all_data$cases[input$obs_idx]
ppc <- bind_rows(
  ppc_summary(fit_poisson, input$obs_idx, y_obs, "Poisson"),
  ppc_summary(fit_negbin, input$obs_idx, y_obs, "Negative binomial")
)
negbin_predictions <- write_negbin_predictions(fit_negbin, input$all_data)

loo_compare_tbl <- as.data.frame(loo::loo_compare(loo_list$poisson, loo_list$negbin))
loo_compare_tbl$model <- rownames(loo_compare_tbl)
loo_compare_tbl <- loo_compare_tbl |>
  select(model, everything())

readr::write_csv(conv, file.path(out_dir, "forecast_convergence_diagnostics.csv"))
readr::write_csv(key_pars, file.path(out_dir, "forecast_key_parameter_summary.csv"))
readr::write_csv(loo_table, file.path(out_dir, "forecast_loo_model_comparison.csv"))
readr::write_csv(loo_compare_tbl, file.path(out_dir, "forecast_loo_compare.csv"))
readr::write_csv(ppc, file.path(out_dir, "forecast_posterior_predictive_checks.csv"))

readr::write_csv(conv, file.path(appendix_dir, "table_s_forecast_convergence_diagnostics.csv"))
readr::write_csv(loo_table, file.path(appendix_dir, "table_s_forecast_loo_model_comparison.csv"))
readr::write_csv(ppc, file.path(appendix_dir, "table_s_forecast_posterior_predictive_checks.csv"))

write_markdown_table(format_numeric_table(conv), file.path(appendix_dir, "table_s_forecast_convergence_diagnostics.md"))
write_markdown_table(format_numeric_table(loo_table), file.path(appendix_dir, "table_s_forecast_loo_model_comparison.md"))
write_markdown_table(format_numeric_table(ppc), file.path(appendix_dir, "table_s_forecast_posterior_predictive_checks.md"))

summary_lines <- c(
  "# Forecast model diagnostics",
  "",
  "## Convergence",
  readLines(file.path(appendix_dir, "table_s_forecast_convergence_diagnostics.md")),
  "",
  "## PSIS-LOO and overdispersion comparison",
  readLines(file.path(appendix_dir, "table_s_forecast_loo_model_comparison.md")),
  "",
  "## Posterior predictive checks",
  readLines(file.path(appendix_dir, "table_s_forecast_posterior_predictive_checks.md")),
  "",
  "## Negative binomial forecast summary",
  readLines(file.path(appendix_dir, "table_s_forecast_negbin_summary.md"))
)
writeLines(summary_lines, file.path(out_dir, "forecast_model_diagnostics.md"))

message("Done: forecast diagnostics written to ", out_dir)
