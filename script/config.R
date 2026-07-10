raw_data_dir <- "./Data/RawGBD"
database_dir <- "./Data/CleanGBD"
outcome_dir <- "./outcome"
appendix_dir <- file.path(outcome_dir, "appendix")

# The cleaned project data and lookup tables are included in this public release.
# Raw GBD exports, the UN population input, and map shapefiles are intentionally
# excluded; see ../DATA_AVAILABILITY.md before adding them locally.
shared_assets_dir <- "./Data"
iso_code_file <- "./Data/iso_code.csv"
# Map inputs are needed only by script/3_b_national_figure.R and must be obtained
# separately under terms suitable for the intended use.
map_global_shp <- file.path(shared_assets_dir, "Map GS(2021)648 - geojson", "globalmap.shp")
map_china_border_shp <- file.path(shared_assets_dir, "Map GS(2021)648 - geojson", "china_border.shp")

# Optional WHO region mapping for national summaries
# Expected columns: at minimum ISO3, and either Region or WHO_region
who_region_file <- "./Data/who_region.csv"

# Target indicators for this project.
# If the raw exports in `raw_data_dir` don't contain some of them,
# scripts will skip the missing measures.
target_measures <- c("Prevalence")

# Global trend settings
target_location_global <- "Global"
target_sex_global <- "Both"
target_age_global <- "Age-standardized"

# Age trend settings (used by `2_b_age_trend.R`)
age_groups_for_age_trend <- c("<20 years", "20-54 years", "55+ years")

# Forecast settings (used by `5_forecast_trend.R`)
forecast_h <- 2

