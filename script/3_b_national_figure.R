# National-level figures: WHO region panels and global maps
# Uses AAPC (2010–last year) computed by 3_a_national_aapc.R

set_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
    setwd(dirname(dirname(script_path)))
  }
}

set_project_root()

library(tidyverse)
library(patchwork)
library(paletteer)
library(Cairo)
library(sf)
library(cowplot)

source("./script/config.R")
source("./script/function.R")

# Load assets and data ---------------------------------------------

# ISO codes and maps
df_map_iso <- readr::read_csv(iso_code_file, show_col_types = FALSE)
df_map <- sf::st_read(map_global_shp, quiet = TRUE)
df_map_border <- sf::st_read(map_china_border_shp, quiet = TRUE)

# WHO region mapping: assumes columns `World regions according to WHO` and `Code`
df_region_map <- readr::read_csv(who_region_file, show_col_types = FALSE) |>
  rename(Region = `World regions according to WHO`, ISO3 = Code)

# Age-standardized rate data: prefer per-year national files
national_dir <- file.path(database_dir, "national_by_year")
csvs <- list.files(national_dir, pattern = "^national_.*\\.csv$", full.names = TRUE)
df_national_all <- purrr::map_dfr(csvs, ~ readr::read_csv(.x, show_col_types = FALSE))

measures_available <- intersect(target_measures, sort(unique(df_national_all$measure_name)))
if (length(measures_available) == 0) stop("No target measures found in prepared data.")
measure_use <- measures_available[[1]]

df_all_rate <- df_national_all |>
  filter(age_name == target_age_global,
         sex_name == 'Both',
         location %in% df_map_iso$location_id,
         year %in% 2010:2023,
         measure_name %in% measures_available) |>
  select(location, location_name, measure_name, year, val, lower, upper)

# Determine AAPC window (should match 3_a_national_aapc.R)
year_min_data <- min(df_all_rate$year, na.rm = TRUE)
year_max_data <- max(df_all_rate$year, na.rm = TRUE)
year_start <- max(2010, year_min_data)
year_end <- year_max_data

# AAPC results from 3_a_national_aapc.R
aapc_file <- file.path(outcome_dir, sprintf("national_aapc_%d_%d.csv", year_start, year_end))
if (!file.exists(aapc_file)) {
  stop("AAPC file not found. Please run 3_a_national_aapc.R first.")
}

df_aapc <- readr::read_csv(aapc_file, show_col_types = FALSE)

# Extract 2010 & last-year rates
df_rate_2010 <- df_all_rate |>
  filter(year == year_start) |>
  select(location_name, val)

df_rate_last <- df_all_rate |>
  filter(year == year_end) |>
  select(location_name, val)

df_aapc_rate <- df_aapc |>
  select(location_name, val = aapc)

# Helper: jitter plot by WHO region (one metric) ------------------------------------------

plot_region_panel <- function(data, value_col, title, legend_title, pal_name = "MetBrewer::Hiroshige", accuracy = NULL) {
  data2 <- data |>
    left_join(df_map_iso, by = c("location_name" = "location_name_1")) |>
    left_join(df_region_map, by = "ISO3")
  # remove literal "(WHO)" from Region labels if present and set factor order
  region_levels <- dplyr::pull(df_region_map, Region) |> unique() |> stringr::str_remove_all("\\s*\\(WHO\\)\\s*")
  data2 <- data2 |>
    dplyr::mutate(Region = stringr::str_remove_all(Region, "\\s*\\(WHO\\)\\s*")) |>
    dplyr::mutate(Region = factor(Region, levels = region_levels))
  if (!"Economy" %in% names(data2)) {
    data2$Economy <- data2$location_name
  } else {
    data2$Economy[is.na(data2$Economy)] <- data2$location_name[is.na(data2$Economy)]
  }
  
  data2 <- data2 |>
    filter(!is.na(Region), !is.na(.data[[value_col]]))
  
  if (nrow(data2) == 0) {
    stop("No data available for ", title)
  }
  
  breaks <- pretty(data2[[value_col]], n = 5)
  
  ggplot(data2) +
    geom_jitter(
      aes(y = Region, x = .data[[value_col]], color = .data[[value_col]]),
      height = 0.2, width = 0
    ) +
    scale_color_gradientn(
      colors = paletteer::paletteer_d(pal_name, direction = -1),
      limits = range(breaks), breaks = breaks,
      labels = if (is.null(accuracy)) scales::label_number(big.mark = ",") else scales::label_number(big.mark = ",", accuracy = accuracy)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0)),
      limits = range(breaks), breaks = breaks,
      labels = if (is.null(accuracy)) scales::label_number(big.mark = ",") else scales::label_number(big.mark = ",", accuracy = accuracy)
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title.position = "top"
    ) +
    labs(
      y = NULL,
      x = legend_title,
      color = legend_title,
      title = title
    ) +
    guides(color = guide_colorbar(barwidth = 12))
}

plot_world_map <- function(data, value_col, title, legend_title, pal_name = "MetBrewer::Hiroshige", accuracy = NULL) {
  d_iso <- data |>
    left_join(df_map_iso, by = c("location_name" = "location_name_1")) |>
    select(ISO3, val = .data[[value_col]])
  
  check_missing <- d_iso$ISO3[!d_iso$ISO3 %in% df_map$SOC]
  if (length(check_missing) > 0) {
    message("Locations missing in map shapefile: ", paste(unique(check_missing), collapse = ", "))
  }
  
  data_map <- df_map |>
    left_join(d_iso, by = c("SOC" = "ISO3"))
  
  if (all(is.na(data_map$val))) {
    stop("No mapped values for ", title)
  }
  
  legend_breaks <- pretty(data_map$val, n = 7)
  
  ggplot(data_map) +
    geom_sf(data = df_map_border, color = "grey", fill = NA) +
    geom_sf(aes(fill = val), show.legend = F) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title.position = "plot"
    ) +
    scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-60, 75)) +
    scale_fill_gradientn(
      colors = paletteer::paletteer_d(pal_name, direction = -1),
      breaks = legend_breaks,
      na.value = "white"
    ) +
    labs(title = title, fill = legend_title)
}

# point panel -------------------------------------------------------------

fig_region_2010 <- plot_region_panel(
  data = df_rate_2010,
  value_col = "val",
  title = sprintf("A: %d", year_start),
  legend_title = sprintf("%s rate", 'Age-standardized prevalence'),
  pal_name = "Redmonder::dPBIRdBu"
)

fig_region_last <- plot_region_panel(
  data = df_rate_last,
  value_col = "val",
  title = sprintf("B: %d", year_end),
  legend_title = sprintf("%s rate", 'Age-standardized prevalence'),
  pal_name = "Redmonder::dPBIRdBu"
)

fig_region_aapc <- plot_region_panel(
  data = df_aapc_rate,
  value_col = "val",
  title = sprintf("C: AAPC (%d–%d)", year_start, year_end),
  legend_title = "AAPC (%)",
  pal_name = "Redmonder::dPBIPuGn",
  accuracy = 0.01
)

# map panel ---------------------------------------------------------------

fig_map_2010 <- plot_world_map(
  data = df_rate_2010,
  value_col = "val",
  title = sprintf("D: %d", year_start),
  legend_title = sprintf("%s rate", measure_use),
  pal_name = "Redmonder::dPBIRdBu"
)

fig_map_last <- plot_world_map(
  data = df_rate_last,
  value_col = "val",
  title = sprintf("E: %d", year_end),
  legend_title = sprintf("%s rate", measure_use),
  pal_name = "Redmonder::dPBIRdBu"
)

fig_map_aapc <- plot_world_map(
  data = df_aapc_rate,
  value_col = "val",
  title = sprintf("F: AAPC (%d–%d)", year_start, year_end),
  legend_title = "AAPC (%)",
  pal_name = "Redmonder::dPBIPuGn",
  accuracy = 0.01
)


fig_region_panels <- fig_region_2010 + fig_region_last + fig_region_aapc +
  patchwork::plot_layout(ncol = 3, guides = 'collect', axes = 'collect_y')&
  theme(legend.position = "bottom",
        legend.margin = margin(t = 5, b = 5, r = 15, l = 15),
        plot.margin = margin(r = 15, t = 5, l = 15, b = 5))

fig <- cowplot::plot_grid(
  fig_region_panels,
  fig_map_2010,
  fig_map_last,
  fig_map_aapc,
  ncol = 2
)

ggsave(
  file.path(outcome_dir, "fig_2_national_maps.png"),
  plot = fig,
  width = 20,
  height = 9
)

message("Done: national WHO-region panels and maps saved under ", outcome_dir)
