#!/usr/bin/env Rscript
###############################################################################
# run_pipeline.R
# Master's Dissertation, Wout Baetens - GPR and monetary policy transmission
#
# Runs the whole project end to end and writes every table (LaTeX) and figure
# (PNG + PDF) to results/, grouped by subject:
#
#   results/descriptive_statistics  Tables 1-4, A6 and Figures B1-B6
#   results/unit_root               Tables A2-A3 (+ series / break plots)
#   results/shock_diagnostics       Tables A4-A5 and Figure A1
#   results/VAR                     Figures 1, 2, C1-C21 (+ IRF data)
#   results/LP                      Figure 3 (+ coefficient table)
#   results/logs                    console output of every stage
#
# Usage (from the project folder):
#   Rscript run_pipeline.R                       # everything
#   Rscript run_pipeline.R descriptive unit_root # only some stages
#   Rscript run_pipeline.R --quick               # fast BVAR settings (testing)
#   Rscript run_pipeline.R --download            # rebuild data from the web
#   Rscript run_pipeline.R --cores=2             # cores for the BVARs
#   Rscript run_pipeline.R VAR --models=Figure1_STBVAR_Baseline
#   Rscript run_pipeline.R shocks --exo-lags=2     # 2-lag exogeneity test
#
# Stages: descriptive, unit_root, shocks, VAR, LP
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
get_flag  <- function(flag) any(args == flag)
get_value <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", key, "="), "", hit[1]) else default
}

#==========================================
# Configuration
#==========================================
# Project folder: the folder of this script (Rscript or source()), otherwise
# the working directory (e.g. the RStudio project folder)
PROJECT_DIR <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(dirname(sub("^--file=", "", f[1]))))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(dirname(of)))
  }
  getwd()
})
if (!file.exists(file.path(PROJECT_DIR, "pipeline", "00_setup.R"))) {
  stop("Cannot find pipeline/00_setup.R in ", PROJECT_DIR,
       ". Run the script from the project folder (open the .Rproj or use setwd()).")
}

DATA_SOURCE <- if (get_flag("--download")) "download" else "frozen"
FIG_FORMATS <- strsplit(get_value("--formats", "png,pdf"), ",")[[1]]
SEED        <- as.integer(get_value("--seed", "20250101"))
# Default: all cores but one, at most 8 (each BVAR worker needs a few hundred MB)
N_CORES     <- as.integer(get_value("--cores", max(1, min(8, parallel::detectCores() - 1))))
EXO_LAGS    <- as.integer(get_value("--exo-lags", "1"))  # 1 = thesis Table A4
VAR_MODELS  <- { m <- get_value("--models", ""); if (nzchar(m)) strsplit(m, ",")[[1]] else NULL }

# Gibbs sampler: thesis setting is 4000 burn-in and 4000 saved draws (every 4th)
if (get_flag("--quick")) {
  BVAR_BURNIN <- 200
  BVAR_NDRAWS <- 200
} else {
  BVAR_BURNIN <- as.integer(get_value("--burnin", "4000"))
  BVAR_NDRAWS <- as.integer(get_value("--ndraws", "4000"))
}

ALL_STAGES <- c("descriptive", "unit_root", "shocks", "VAR", "LP")
stages <- intersect(ALL_STAGES, args)
if (length(stages) == 0) stages <- ALL_STAGES

STAGE_FILES <- c(
  descriptive = "02_descriptive_statistics.R",
  unit_root   = "03_unit_root.R",
  shocks      = "04_shock_diagnostics.R",
  VAR         = "05_var.R",
  LP          = "06_lp.R"
)

#==========================================
# Run
#==========================================
setwd(PROJECT_DIR)

# Install missing CRAN packages before anything else
REQUIRED_PACKAGES <- c("dplyr", "tidyr", "readr", "lubridate", "ggplot2", "patchwork",
                       "sandwich", "urca", "strucchange", "dynlm", "FinTS", "moments",
                       "xtable", "knitr", "kableExtra", "gt", "stargazer", "corrplot")
missing_pkgs <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace,
                                          logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  repos <- getOption("repos")
  if (is.null(repos) || any(repos == "@CRAN@")) repos <- "https://cloud.r-project.org"
  install.packages(missing_pkgs, repos = repos)
  still_missing <- missing_pkgs[!vapply(missing_pkgs, requireNamespace, logical(1),
                                        quietly = TRUE)]
  if (length(still_missing)) {
    stop("Could not install: ", paste(still_missing, collapse = ", "),
         ". Install them manually and rerun.")
  }
}

source(file.path(PROJECT_DIR, "pipeline", "00_setup.R"))

message("==============================================================")
message(" GPR and MP pipeline")
message(" data: ", DATA_SOURCE, " | stages: ", paste(stages, collapse = ", "))
message(" BVAR: burn-in ", BVAR_BURNIN, ", draws ", BVAR_NDRAWS,
        " | cores: ", N_CORES, " | seed: ", SEED)
message("==============================================================")

timings <- data.frame(Stage = character(), Minutes = numeric())
t_start <- Sys.time()

source(file.path(PROJECT_DIR, "pipeline", "01_data.R"))

for (st in stages) {
  message("\n>>> Stage: ", st)
  t0 <- Sys.time()
  set.seed(SEED)
  log_file <- file.path(DIRS$logs, paste0(st, ".log"))
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  tryCatch(
    source(file.path(PROJECT_DIR, "pipeline", STAGE_FILES[[st]]), echo = FALSE),
    finally = { sink(); close(log_con) }
  )
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  timings <- rbind(timings, data.frame(Stage = st, Minutes = round(mins, 2)))
  message(sprintf("<<< %s finished in %.1f min", st, mins))
}

write.csv(timings, file.path(DIRS$logs, "timings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(DIRS$logs, "sessionInfo.txt"))
message(sprintf("\nPipeline finished in %.1f min. Results in %s",
                as.numeric(difftime(Sys.time(), t_start, units = "mins")), RESULTS_DIR))
