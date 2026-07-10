#####################################
## @Description: Summarize local GBD 2023 export metadata
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
library(readr)
library(purrr)

source("./script/config.R")
source("./script/function.R")

dir.create(appendix_dir, recursive = TRUE, showWarnings = FALSE)

zip_export_times <- function(zip_path) {
  if (!file.exists(zip_path)) {
    return(as.POSIXct(character()))
  }

  info <- utils::unzip(zip_path, list = TRUE)
  keep <- grepl("\\.csv$|citation\\.txt$", tolower(info$Name))
  info$Date[keep]
}

summarize_csv_export <- function(csv_path) {
  df <- readr::read_csv(csv_path, show_col_types = FALSE)
  zip_path <- sub("\\.csv$", ".zip", csv_path)
  zip_times <- zip_export_times(zip_path)
  years <- sort(unique(as.integer(df$year)))

  tibble(
    file = basename(csv_path),
    n_rows = nrow(df),
    year_range = if (length(years) > 0) paste0(min(years), "-", max(years)) else "",
    n_locations = n_distinct(df$location_name),
    sex = paste(sort(unique(df$sex_name)), collapse = "; "),
    metrics = paste(sort(unique(df$metric_name)), collapse = "; "),
    age_strata = paste(sort(unique(df$age_name)), collapse = "; "),
    zip_timestamp_min = if (length(zip_times) > 0) format(min(zip_times), "%Y-%m-%d %H:%M:%S") else "",
    zip_timestamp_max = if (length(zip_times) > 0) format(max(zip_times), "%Y-%m-%d %H:%M:%S") else ""
  )
}

csv_files <- sort(list.files(
  raw_data_dir,
  pattern = "^IHME-GBD_2023_DATA-.*\\.csv$",
  full.names = TRUE
))

if (length(csv_files) == 0) {
  stop("No GBD CSV exports found in ", raw_data_dir)
}

export_table <- purrr::map_dfr(csv_files, summarize_csv_export)

readr::write_csv(export_table, file.path(appendix_dir, "table_s_gbd_export_metadata.csv"))

write_markdown_table(
  export_table,
  file.path(appendix_dir, "table_s_gbd_export_metadata.md")
)

all_times <- do.call(
  c,
  purrr::map(csv_files, ~ zip_export_times(sub("\\.csv$", ".zip", .x)))
)
all_dates <- sort(unique(format(all_times, "%Y-%m-%d")))
available_age_strata <- export_table$age_strata |>
  strsplit("; ", fixed = TRUE) |>
  unlist(use.names = FALSE) |>
  unique() |>
  sort()

summary_lines <- c(
  "{",
  paste0('  "raw_dir": "', normalizePath(raw_data_dir, winslash = "/", mustWork = FALSE), '",'),
  paste0('  "n_exports": ', nrow(export_table), ","),
  paste0('  "export_dates": [', paste(sprintf('"%s"', all_dates), collapse = ", "), "],"),
  paste0(
    '  "final_export_timestamp": "',
    if (length(all_times) > 0) format(max(all_times), "%Y-%m-%d %H:%M:%S") else "",
    '",'
  ),
  paste0(
    '  "available_age_strata": [',
    paste(sprintf('"%s"', available_age_strata), collapse = ", "),
    "]"
  ),
  "}"
)

writeLines(summary_lines, file.path(appendix_dir, "gbd_export_metadata_summary.json"), useBytes = TRUE)

message("Wrote GBD export metadata table and summary.")
message("GBD export dates: ", paste(all_dates, collapse = ", "))
message(
  "Final export timestamp: ",
  if (length(all_times) > 0) format(max(all_times), "%Y-%m-%d %H:%M:%S") else ""
)
