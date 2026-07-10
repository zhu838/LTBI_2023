
# setup -------------------------------------------------------------------

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
library(patchwork)
library(paletteer)
library(rstan)

# RStan recommended options
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("./script/config.R")
source("./script/function.R")

## Data files (adjust if needed)
pop_file <- "./Data/unpopulation_dataportal_20260305182132.csv"
gr_rate_file <- file.path(database_dir, "global_regional_rate.csv")

cat("Population file:", pop_file, "\n")
cat("Global/regional rate file:", gr_rate_file, "\n")

# load data ---------------------------------------------------------------

# Read population projections (5-year age groups, 2010-2050)
pop_raw <- readr::read_csv(pop_file, show_col_types = FALSE)

## GBD prevalence (Global, Both sexes, Prevalence rate) -----------------

df_raw_rate <- readr::read_csv(gr_rate_file, show_col_types = FALSE)

df_rate <- df_raw_rate |>
  filter(
    measure_name %in% target_measures,
    metric_name == "Rate",
    location_name == target_location_global,
    sex_name == target_sex_global
  ) |>
  rename(
    age_group = age_name,
    prevalence_rate = val
  ) |>
  mutate(
    year = as.integer(year),
    age_group = case_when(
      age_group == "<20 years" ~ "<20",
      age_group == "20-54 years" ~ "20-54",
      age_group == "55+ years" ~ "55+",
      TRUE ~ age_group
    ),
    location = location_name
  ) |>
  select(location, year, age_group, prevalence_rate)

## 4. Choose target location and age groups ------------------

target_location <- target_location_global  # e.g. "Global" from config.R

# Forecast age groups aligned to GBD: <20, 20-54, 55+
target_age_bins <- c("<20", "20-54", "55+")

age_groups_prev <- sort(unique(df_rate$age_group))
cat("Age groups (prevalence, standardized):", paste(age_groups_prev, collapse = ", "), "\n")

## 5. Aggregate population 5-year groups to GBD age bands ----

# Load country-level LocationId list (same file used by 3_a_national_aapc.R)
# This filters pop_raw to country rows only, excluding regional/world aggregates
df_iso_pop <- readr::read_csv(iso_code_file, show_col_types = FALSE)
country_iso3 <- toupper(as.character(df_iso_pop$ISO3))

# 5-year groups -> three GBD bands
# Value is in persons (UN WPP Data Portal). Filter to country rows only,
# then sum across countries to get world totals per age band.
pop_age3 <- pop_raw |>
  mutate(Iso3 = toupper(as.character(Iso3))) |>
  filter(Iso3 %in% country_iso3) |>
  transmute(year      = as.integer(Time),
            age_start = as.integer(AgeStart),
            # Value is in persons
            population = as.numeric(Value),
            age_group = case_when(age_start < 20              ~ "<20",
                                  age_start >= 20 & age_start < 55 ~ "20-54",
                                  age_start >= 55              ~ "55+",
                                  TRUE ~ NA_character_)) |>
  filter(!is.na(age_group)) |>
  group_by(year, age_group) |>
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop") |>
  mutate(location = target_location) |>
  select(location, year, age_group, population)

# Harmonized population table used below
pop_harmonized <- pop_age3
cat(sprintf("Population (persons) after unit fix - Global 2023 sample:\n"))
pop_harmonized |> filter(year == 2023) |> print()

## 6. Prepare prevalence data for target location and target years
df_rate_sub <- df_rate |>
	filter(location == target_location,
	       year >= 2010, year <= 2023,
	       age_group %in% target_age_bins) |>
	# try to normalize prevalence_rate: if per 100k convert to proportion
	mutate(prevalence_rate = as.numeric(prevalence_rate))

# If prevalence values look large (>1), assume per 100k
if (max(df_rate_sub$prevalence_rate, na.rm = TRUE) > 10){
	df_rate_sub <- df_rate_sub |> mutate(prevalence_rate = prevalence_rate / 1e5)
}

if (!all(target_age_bins %in% unique(df_rate_sub$age_group))) stop("GBD prevalence age groups do not cover all target bands")

## 7. Combine prevalence and population, compute cases -------------
# Merge historical population (2010-2023) with prevalence
pop_hist <- pop_harmonized |> filter(year >= 2010 & year <= 2023, location == target_location)

hist_data <- df_rate_sub |>
	left_join(pop_hist, by = c("location", "year", "age_group")) |>
	mutate(cases = prevalence_rate * population,
				 cases = round(cases))

## 8. Create future grid (2024-2050) and combine ---------------------
years_proj <- 2024:2050
pop_proj <- pop_harmonized |> filter(year %in% years_proj, location == target_location)

future_grid <- expand_grid(location = target_location, year = years_proj, age_group = target_age_bins) |>
	left_join(pop_proj, by = c("location","year","age_group")) |>
	mutate(prevalence_rate = NA_real_, cases = NA_integer_)

all_data <- bind_rows(hist_data, future_grid) |> arrange(year, age_group)

## 9. APC indexing ----------------------------------------------
all_data <- all_data |>
	mutate(age_index = as.integer(factor(age_group, levels = target_age_bins)),
				 period_index = as.integer(factor(year, levels = sort(unique(year)))))

# cohort raw index (period - age + A)
A <- length(target_age_bins)
all_data <- all_data |>
	mutate(cohort_raw = period_index - age_index + A,
				 cohort_index = as.integer(factor(cohort_raw)))

## 10. Build APC model in Stan (Poisson with offset) -----------------

# Check population availability (no NA allowed in offset)
if (any(is.na(all_data$population) & all_data$year >= 2024)) {
  stop("Missing population for some projection years - cannot run predictions. Check population WPP file and harmonization.")
}
if (any(is.na(all_data$population))) {
  stop("Missing population for some historical rows - cannot run APC model. Check population WPP file and harmonization.")
}

# Prepare outcome and observation indicator without NA
y_vec <- as.integer(all_data$cases)
y_vec[is.na(y_vec)] <- 0L
is_obs_vec <- ifelse(is.na(all_data$cases), 0L, 1L)

# Stan data list
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

## Diagnostics and safety checks before compiling Stan model
T_val <- stan_data$T
K_val <- stan_data$K
A_val <- stan_data$A
n_obs <- sum(as.integer(stan_data$is_obs))
min_pop <- min(all_data$population, na.rm = TRUE)
max_pop <- max(all_data$population, na.rm = TRUE)

cat(sprintf("Stan data diagnostics: N=%d, A=%d, T=%d, K=%d, observed_rows=%d\n", stan_data$N, A_val, T_val, K_val, n_obs))
cat(sprintf("Population range: min=%g, max=%g\n", min_pop, max_pop))

if (T_val < 2) stop(sprintf("Insufficient period levels for RW1: T=%d (need >=2)", T_val))
if (K_val < 2) stop(sprintf("Insufficient cohort levels for RW1: K=%d (need >=2)", K_val))
if (!all(is.finite(stan_data$log_E))) stop("Non-finite values found in log_E (check population > 0)")
if (any(all_data$population <= 0, na.rm = TRUE)) stop("Non-positive population values found (must be >0)")

## Ensure is_obs and y are integer vectors of correct length
stan_data$is_obs <- as.integer(stan_data$is_obs)
stan_data$y <- as.integer(stan_data$y)

# Stan model code: APC with RW1 priors for period & cohort, IID age
stan_code <- "
data {
  int<lower=1> N;
  int<lower=0> y[N];
  int<lower=0,upper=1> is_obs[N];
  vector[N] log_E;
  int<lower=1> A;
  int<lower=1> T;
  int<lower=1> K;
  int<lower=1,upper=A> age_index[N];
  int<lower=1,upper=T> period_index[N];
  int<lower=1,upper=K> cohort_index[N];
}
parameters {
  real alpha;                  // intercept
  vector[A] age_raw;           // age effects (centered later)
  // RW2 for period: second-order differences -> preserves trend momentum into projection
  vector[T-1] d1_period;       // first-order increments
  // RW2 for cohort
  vector[K-1] d1_cohort;       // first-order increments
  real<lower=0> sigma_age;
  real<lower=0> sigma_period;  // scale of RW2 second differences
  real<lower=0> sigma_cohort;
}
transformed parameters {
  vector[A] age_effect;
  vector[T] period_effect;
  vector[K] cohort_effect;
  vector[N] eta;

  // Center age effects to improve identifiability
  age_effect = (age_raw - mean(age_raw)) * sigma_age;

  // RW2 for period: period_effect[t] = 2*period_effect[t-1] - period_effect[t-2] + eps
  // Equivalently, first accumulate d1 (RW1 on increments), then accumulate into levels
  period_effect[1] = 0;
  period_effect[2] = d1_period[1] * sigma_period;
  for (t in 3:T)
    period_effect[t] = 2 * period_effect[t-1] - period_effect[t-2]
                       + (d1_period[t-1] - d1_period[t-2]) * sigma_period;

  // RW2 for cohort
  cohort_effect[1] = 0;
  cohort_effect[2] = d1_cohort[1] * sigma_cohort;
  for (k in 3:K)
    cohort_effect[k] = 2 * cohort_effect[k-1] - cohort_effect[k-2]
                       + (d1_cohort[k-1] - d1_cohort[k-2]) * sigma_cohort;

  // Linear predictor including offset log_E
  for (n in 1:N) {
    eta[n] = alpha
             + age_effect[age_index[n]]
             + period_effect[period_index[n]]
             + cohort_effect[cohort_index[n]]
             + log_E[n];
  }
}
model {
  // Priors
  alpha ~ normal(0, 5);
  age_raw ~ normal(0, 1);
  // RW2: second differences are i.i.d. Normal
  d1_period ~ normal(0, 1);
  d1_cohort ~ normal(0, 1);
  sigma_age ~ exponential(1);
  sigma_period ~ exponential(1);
  sigma_cohort ~ exponential(1);

  // Likelihood only for observed counts
  for (n in 1:N) {
    if (is_obs[n] == 1)
      y[n] ~ poisson_log(eta[n]);
  }
}
generated quantities {
  vector[N] mu;      // posterior mean counts per row
  for (n in 1:N) {
    mu[n] = exp(eta[n]);
  }
}
"

## 11. Compile and fit Stan model ----------------------------------

## Compile Stan model (capture compilation errors)
sm <- NULL
sm <- tryCatch({
  rstan::stan_model(model_code = stan_code)
}, error = function(e){
  cat("Stan model compilation failed:\n")
  message(e$message)
  stop(e)
})

## Run sampling with clearer error reporting; if chains fail, print diagnostics
fit <- tryCatch({
  rstan::sampling(
    object = sm,
    data = stan_data,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    thin = 2,
    seed = 1234,
    control = list(adapt_delta = 0.9, max_treedepth = 12)
  )
}, error = function(e){
  cat("Stan sampling failed. Error message:\n")
  message(e$message)
  cat("You can try running with 'chains=1' and smaller iter/warmup to debug interactively.\n")
  stop(e)
})

print(fit, pars = c("alpha", "sigma_age", "sigma_period", "sigma_cohort"))

 ## 11b. Convergence diagnostics -------------------------------------------

 ## Extract R-hat (potential scale reduction factor)
 rhat_vals <- rstan::summary(fit)$summary[, "Rhat"]
 max_rhat <- max(rhat_vals, na.rm = TRUE)
 cat(sprintf("\n=== Stan Convergence Diagnostics ===\n"))
 cat(sprintf("Max R-hat across all parameters: %.4f\n", max_rhat))
 if (max_rhat > 1.05) {
   warning(sprintf("Max R-hat = %.4f exceeds 1.05 threshold; chains may not have converged.", max_rhat))
 }

 ## Extract effective sample size
 ess_vals <- rstan::summary(fit)$summary[, "n_eff"]
 min_ess <- min(ess_vals, na.rm = TRUE)
 cat(sprintf("Min effective sample size across all parameters: %.0f\n", min_ess))
 if (min_ess < 400) {
   warning(sprintf("Min effective sample size = %.0f is below 400; consider increasing iterations.", min_ess))
 }

 ## Check for divergent transitions
 sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
 divergent_per_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
 total_divergent <- sum(divergent_per_chain)
 cat(sprintf("Total divergent transitions across all chains: %d\n", total_divergent))
 if (total_divergent > 0) {
   warning(sprintf("%d divergent transitions detected; consider increasing adapt_delta.", total_divergent))
 }

 ## Print diagnostics summary for key parameters
 cat("\nR-hat and effective sample size for key parameters:\n")
 key_pars <- c("alpha", "sigma_age", "sigma_period", "sigma_cohort")
 summ <- rstan::summary(fit, pars = key_pars)$summary
 print(summ[, c("mean", "sd", "n_eff", "Rhat")], digits = 4)
 cat("=====================================\n\n")

## 12. Posterior summaries and prevalence/cases ---------------------

post_mu <- rstan::extract(fit, "mu")$mu   # iterations x N

mu_mean <- colMeans(post_mu)
mu_lwr  <- apply(post_mu, 2, quantile, probs = 0.025)
mu_upr  <- apply(post_mu, 2, quantile, probs = 0.975)

all_data <- all_data |>
	mutate(fitted_mean = mu_mean,
			 fitted_lwr  = mu_lwr,
			 fitted_upr  = mu_upr,
			 # divide posterior case counts by population -> proportion -> *1e5 -> per 100k
			 pred_prev_mean = (fitted_mean / population) * 1e5,
			 pred_prev_lwr  = (fitted_lwr  / population) * 1e5,
			 pred_prev_upr  = (fitted_upr  / population) * 1e5)

results_proj <- all_data |> filter(year >= 2024)

# Save CSV outputs (posterior mean & 95% CrI)
out_dir <- file.path(outcome_dir, "forecast")
if (!dir.exists(out_dir)) dir.create(out_dir)

write_csv(
	results_proj |> select(location, year, age_group, population,
										 pred_prev_mean, pred_prev_lwr, pred_prev_upr,
										 fitted_mean, fitted_lwr, fitted_upr),
	file.path(out_dir, sprintf("predictions_rstan_%s_2024_2050.csv", target_location))
)

## 13. Figure 3: APC forecast – 4-panel layout --------------------------

age_pal <- as.character(paletteer::paletteer_d("MetBrewer::Hiroshige",
                                                    n = 3, direction = -1))
age_colors  <- stats::setNames(age_pal, c("<20", "20-54", "55+"))
overall_col <- "#00798CFF"

common_x <- scale_x_continuous(breaks = seq(2010, 2050, 10),
                                expand = expansion(add = c(0.5, 1)))

common_theme <- function(show_legend = FALSE) {
  theme_bw() + theme(
    plot.title.position  = "plot",
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    legend.position      = if (show_legend) "bottom" else "none",
    legend.title         = element_text(size = 9),
    legend.text          = element_text(size = 8)
  )
}

df_overall <- all_data |>
  group_by(year) |>
  summarise(
    total_pop      = sum(population,  na.rm = TRUE),
    total_fit_mean = sum(fitted_mean, na.rm = TRUE),
    total_fit_lwr  = sum(fitted_lwr,  na.rm = TRUE),
    total_fit_upr  = sum(fitted_upr,  na.rm = TRUE),
    total_obs      = sum(cases,       na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    # rate per 100k
    rate_mean    = (total_fit_mean / total_pop) * 1e5,
    rate_lwr     = (total_fit_lwr  / total_pop) * 1e5,
    rate_upr     = (total_fit_upr  / total_pop) * 1e5,
    rate_obs     = if_else(year <= 2023,
                                  (total_obs / total_pop) * 1e5, NA_real_),
    # cases in millions
    cases_mean_m = total_fit_mean / 1e6,
    cases_lwr_m  = total_fit_lwr  / 1e6,
    cases_upr_m  = total_fit_upr  / 1e6,
    cases_obs_m  = if_else(year <= 2023, total_obs / 1e6, NA_real_)
  )

plot_df <- all_data |>
  mutate(
    obs_prev      = if_else(year <= 2023, prevalence_rate * 1e5, NA_real_),
    obs_cases_m   = if_else(year <= 2023, cases / 1e6,           NA_real_),
    fitted_mean_m = fitted_mean / 1e6,
    fitted_lwr_m  = fitted_lwr  / 1e6,
    fitted_upr_m  = fitted_upr  / 1e6
  )

y_A <- scales::pretty_breaks(n = 6)(
  range(c(df_overall$rate_obs, df_overall$rate_mean,
          df_overall$rate_lwr, df_overall$rate_upr), na.rm = TRUE))

fig3_A <- ggplot() +
  geom_vline(xintercept = 2023.5,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = df_overall,
              aes(x = year, ymin = rate_lwr, ymax = rate_upr),
              fill = overall_col, alpha = 0.25) +
  geom_point(data = filter(df_overall, year <= 2023),
             aes(x = year, y = rate_obs),
             colour = overall_col, size = 1.4, shape = 16) +
  geom_line(data = filter(df_overall, year <= 2023),
            aes(x = year, y = rate_mean),
            colour = overall_col, linewidth = 0.7) +
  geom_line(data = filter(df_overall, year >= 2023),
            aes(x = year, y = rate_mean),
            colour = overall_col, linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 2016.5, y = max(y_A),
           label = "Historical", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  annotate("text", x = 2037, y = max(y_A),
           label = "Forecasted", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  common_x +
  scale_y_continuous(limits = range(y_A), breaks = y_A,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "A", x = "Year", y = "Prevalence rate (per 100,000)") +
  common_theme()

y_B <- scales::pretty_breaks(n = 6)(
  range(c(plot_df$obs_prev, plot_df$pred_prev_mean,
          plot_df$pred_prev_lwr, plot_df$pred_prev_upr), na.rm = TRUE))

fig3_B <- ggplot() +
  geom_vline(xintercept = 2023.5,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = plot_df,
              aes(x = year, ymin = pred_prev_lwr, ymax = pred_prev_upr,
                  fill = age_group),
              alpha = 0.20) +
  geom_point(data = filter(plot_df, year <= 2023),
             aes(x = year, y = obs_prev, colour = age_group),
             size = 1.4, shape = 16) +
  geom_line(data = filter(plot_df, year <= 2023),
            aes(x = year, y = pred_prev_mean, colour = age_group),
            linewidth = 0.7) +
  geom_line(data = filter(plot_df, year >= 2023),
            aes(x = year, y = pred_prev_mean, colour = age_group),
            linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 2016.5, y = max(y_B),
           label = "Historical", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  annotate("text", x = 2037, y = max(y_B),
           label = "Forecasted", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  scale_colour_manual(values = age_colors, name = "Age group") +
  scale_fill_manual(values = age_colors,   name = "Age group") +
  common_x +
  scale_y_continuous(limits = range(y_B), breaks = y_B,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "B", x = "Year", y = "Prevalence rate (per 100,000)") +
  common_theme(show_legend = TRUE)

y_C <- scales::pretty_breaks(n = 6)(
  range(c(df_overall$cases_obs_m, df_overall$cases_mean_m,
          df_overall$cases_lwr_m, df_overall$cases_upr_m), na.rm = TRUE))

fig3_C <- ggplot() +
  geom_vline(xintercept = 2023.5,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = df_overall,
              aes(x = year, ymin = cases_lwr_m, ymax = cases_upr_m),
              fill = overall_col, alpha = 0.25) +
  geom_point(data = filter(df_overall, year <= 2023),
             aes(x = year, y = cases_obs_m),
             colour = overall_col, size = 1.4, shape = 16) +
  geom_line(data = filter(df_overall, year <= 2023),
            aes(x = year, y = cases_mean_m),
            colour = overall_col, linewidth = 0.7) +
  geom_line(data = filter(df_overall, year >= 2023),
            aes(x = year, y = cases_mean_m),
            colour = overall_col, linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 2016.5, y = max(y_C),
           label = "Historical", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  annotate("text", x = 2037, y = max(y_C),
           label = "Forecasted", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  common_x +
  scale_y_continuous(limits = range(y_C), breaks = y_C,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "C", x = "Year", y = "Number of prevalent cases (million)") +
  common_theme()

y_D <- scales::pretty_breaks(n = 6)(
  range(c(plot_df$obs_cases_m, plot_df$fitted_mean_m,
          plot_df$fitted_lwr_m, plot_df$fitted_upr_m), na.rm = TRUE))

fig3_D <- ggplot() +
  geom_vline(xintercept = 2023.5,
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = plot_df,
              aes(x = year, ymin = fitted_lwr_m, ymax = fitted_upr_m,
                  fill = age_group),
              alpha = 0.20) +
  geom_point(data = filter(plot_df, year <= 2023),
             aes(x = year, y = obs_cases_m, colour = age_group),
             size = 1.4, shape = 16) +
  geom_line(data = filter(plot_df, year <= 2023),
            aes(x = year, y = fitted_mean_m, colour = age_group),
            linewidth = 0.7) +
  geom_line(data = filter(plot_df, year >= 2023),
            aes(x = year, y = fitted_mean_m, colour = age_group),
            linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 2016.5, y = max(y_D),
           label = "Historical", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  annotate("text", x = 2037, y = max(y_D),
           label = "Forecasted", colour = "grey45", size = 3,
           vjust = 1.6, hjust = 0.5) +
  scale_colour_manual(values = age_colors, name = "Age group") +
  scale_fill_manual(values = age_colors,   name = "Age group") +
  common_x +
  scale_y_continuous(limits = range(y_D), breaks = y_D,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "D", x = "Year", y = "Number of prevalent cases (million)") +
  common_theme(show_legend = TRUE)

fig_3 <- patchwork::wrap_plots(list(fig3_A, fig3_B, fig3_C, fig3_D),
                               ncol = 2, nrow = 2) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(outcome_dir, "fig_3_forecast_Global.png"),
  plot     = fig_3,
  width    = 14,
  height   = 10
)

cat("Done. Figure 3 and predictions saved to:", out_dir, "\n")

