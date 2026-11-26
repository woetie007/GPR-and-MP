###############################################################################
# 01_data.R
# Builds `final_data`, the single monthly Euro Area dataset used by every stage.
#
# DATA_SOURCE = "frozen"   (default) reads data_thesis_r0838656.csv, the exact
#                          `final_data` written by shock_analysis_and_purging.R
#                          right after Merge.R. Results are reproducible and no
#                          internet access is required.
# DATA_SOURCE = "download" re-runs Merge.R (Eurostat, ECB SDW, FRED and the
#                          hfdshocks package). Needs internet access; the series
#                          may have been revised since the thesis was written.
###############################################################################

if (DATA_SOURCE == "download") {
  message("Downloading and merging raw data via Merge.R ...")
  dl_env <- new.env()
  old_wd <- setwd(PROJECT_DIR)
  sys.source(file.path(PROJECT_DIR, "Merge.R"), envir = dl_env)
  setwd(old_wd)
  final_data <- dl_env$final_data
  rm(dl_env)
  write_csv(final_data, file.path(DIRS$data, "final_data_downloaded.csv"))
} else {
  final_data <- read_csv(file.path(PROJECT_DIR, "data_thesis_r0838656.csv"),
                         show_col_types = FALSE, guess_max = 10000)
}

final_data <- final_data %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= as.Date("1999-01-01"), Date <= as.Date("2025-12-01")) %>%
  arrange(Country, Date)

message(sprintf("  final_data: %d rows x %d columns (%s to %s)",
                nrow(final_data), ncol(final_data),
                min(final_data$Date), max(final_data$Date)))

# Euro-Area-only view used by the time-series stages
ea_data <- final_data %>% filter(Country == "EA")
