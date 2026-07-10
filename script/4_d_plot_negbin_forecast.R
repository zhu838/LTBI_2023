#####################################
## @Description: Plot negative-binomial APC forecast
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
library(readr)
library(ggplot2)
library(patchwork)
library(paletteer)

source("./script/config.R")

pred_file <- file.path(outcome_dir, "forecast", "diagnostics", "predictions_rstan_negbin_Global_2010_2050_all.csv")
if (!file.exists(pred_file)) {
  stop("Missing negative-binomial forecast predictions: ", pred_file)
}

all_data <- readr::read_csv(pred_file, show_col_types = FALSE) |>
  mutate(age_group = factor(age_group, levels = c("<20", "20-54", "55+")))

age_colors <- c("<20" = "#4C78A8", "20-54" = "#F58518", "55+" = "#E45756")
overall_col <- "#00798CFF"

common_x <- scale_x_continuous(breaks = seq(2010, 2050, 10), expand = expansion(add = c(0.5, 1)))
common_theme <- function(show_legend = FALSE) {
  theme_bw() +
    theme(
      plot.title.position = "plot",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

df_overall <- all_data |>
  group_by(year) |>
  summarise(
    total_pop = sum(population, na.rm = TRUE),
    total_fit_mean = sum(fitted_mean, na.rm = TRUE),
    total_fit_lwr = sum(fitted_lwr, na.rm = TRUE),
    total_fit_upr = sum(fitted_upr, na.rm = TRUE),
    total_obs = sum(cases, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    rate_mean = total_fit_mean / total_pop * 1e5,
    rate_lwr = total_fit_lwr / total_pop * 1e5,
    rate_upr = total_fit_upr / total_pop * 1e5,
    rate_obs = if_else(year <= 2023, total_obs / total_pop * 1e5, NA_real_),
    cases_mean_m = total_fit_mean / 1e6,
    cases_lwr_m = total_fit_lwr / 1e6,
    cases_upr_m = total_fit_upr / 1e6,
    cases_obs_m = if_else(year <= 2023, total_obs / 1e6, NA_real_)
  )

plot_df <- all_data |>
  mutate(
    obs_prev = if_else(year <= 2023, prevalence_rate * 1e5, NA_real_),
    obs_cases_m = if_else(year <= 2023, cases / 1e6, NA_real_),
    fitted_mean_m = fitted_mean / 1e6,
    fitted_lwr_m = fitted_lwr / 1e6,
    fitted_upr_m = fitted_upr / 1e6
  )

make_breaks <- function(...) {
  scales::pretty_breaks(n = 6)(range(c(...), na.rm = TRUE))
}

y_A <- make_breaks(df_overall$rate_obs, df_overall$rate_mean, df_overall$rate_lwr, df_overall$rate_upr)
fig3_A <- ggplot() +
  geom_vline(xintercept = 2023.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = df_overall, aes(x = year, ymin = rate_lwr, ymax = rate_upr), fill = overall_col, alpha = 0.25) +
  geom_point(data = filter(df_overall, year <= 2023), aes(x = year, y = rate_obs), colour = overall_col, size = 1.4) +
  geom_line(data = filter(df_overall, year <= 2023), aes(x = year, y = rate_mean), colour = overall_col, linewidth = 0.7) +
  geom_line(data = filter(df_overall, year >= 2023), aes(x = year, y = rate_mean), colour = overall_col, linetype = "dashed", linewidth = 0.7) +
  common_x +
  scale_y_continuous(limits = range(y_A), breaks = y_A, expand = expansion(mult = c(0, 0))) +
  labs(title = "A", x = "Year", y = "Prevalence rate (per 100,000)") +
  common_theme()

y_B <- make_breaks(plot_df$obs_prev, plot_df$pred_prev_mean, plot_df$pred_prev_lwr, plot_df$pred_prev_upr)
fig3_B <- ggplot() +
  geom_vline(xintercept = 2023.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = plot_df, aes(x = year, ymin = pred_prev_lwr, ymax = pred_prev_upr, fill = age_group), alpha = 0.20) +
  geom_point(data = filter(plot_df, year <= 2023), aes(x = year, y = obs_prev, colour = age_group), size = 1.4) +
  geom_line(data = filter(plot_df, year <= 2023), aes(x = year, y = pred_prev_mean, colour = age_group), linewidth = 0.7) +
  geom_line(data = filter(plot_df, year >= 2023), aes(x = year, y = pred_prev_mean, colour = age_group), linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = age_colors, name = "Age group") +
  scale_fill_manual(values = age_colors, name = "Age group") +
  common_x +
  scale_y_continuous(limits = range(y_B), breaks = y_B, expand = expansion(mult = c(0, 0))) +
  labs(title = "B", x = "Year", y = "Prevalence rate (per 100,000)") +
  common_theme(show_legend = TRUE)

y_C <- make_breaks(df_overall$cases_obs_m, df_overall$cases_mean_m, df_overall$cases_lwr_m, df_overall$cases_upr_m)
fig3_C <- ggplot() +
  geom_vline(xintercept = 2023.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = df_overall, aes(x = year, ymin = cases_lwr_m, ymax = cases_upr_m), fill = overall_col, alpha = 0.25) +
  geom_point(data = filter(df_overall, year <= 2023), aes(x = year, y = cases_obs_m), colour = overall_col, size = 1.4) +
  geom_line(data = filter(df_overall, year <= 2023), aes(x = year, y = cases_mean_m), colour = overall_col, linewidth = 0.7) +
  geom_line(data = filter(df_overall, year >= 2023), aes(x = year, y = cases_mean_m), colour = overall_col, linetype = "dashed", linewidth = 0.7) +
  common_x +
  scale_y_continuous(limits = range(y_C), breaks = y_C, expand = expansion(mult = c(0, 0))) +
  labs(title = "C", x = "Year", y = "Number of prevalent cases (million)") +
  common_theme()

y_D <- make_breaks(plot_df$obs_cases_m, plot_df$fitted_mean_m, plot_df$fitted_lwr_m, plot_df$fitted_upr_m)
fig3_D <- ggplot() +
  geom_vline(xintercept = 2023.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = plot_df, aes(x = year, ymin = fitted_lwr_m, ymax = fitted_upr_m, fill = age_group), alpha = 0.20) +
  geom_point(data = filter(plot_df, year <= 2023), aes(x = year, y = obs_cases_m, colour = age_group), size = 1.4) +
  geom_line(data = filter(plot_df, year <= 2023), aes(x = year, y = fitted_mean_m, colour = age_group), linewidth = 0.7) +
  geom_line(data = filter(plot_df, year >= 2023), aes(x = year, y = fitted_mean_m, colour = age_group), linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = age_colors, name = "Age group") +
  scale_fill_manual(values = age_colors, name = "Age group") +
  common_x +
  scale_y_continuous(limits = range(y_D), breaks = y_D, expand = expansion(mult = c(0, 0))) +
  labs(title = "D", x = "Year", y = "Number of prevalent cases (million)") +
  common_theme(show_legend = TRUE)

fig_3 <- patchwork::wrap_plots(list(fig3_A, fig3_B, fig3_C, fig3_D), ncol = 2, nrow = 2) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(outcome_dir, "fig_3_forecast_Global.png"), fig_3, width = 14, height = 10)
ggsave(file.path(outcome_dir, "fig_3_forecast_Global_negbin.png"), fig_3, width = 14, height = 10)

message("Done: negative-binomial Figure 3 written to ", outcome_dir)
