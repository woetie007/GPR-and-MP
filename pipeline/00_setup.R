###############################################################################
# 00_setup.R
# Packages, paths and small output helpers shared by every pipeline stage.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(sandwich)
  library(urca)
  library(strucchange)
  library(dynlm)
  library(FinTS)
  library(moments)
  library(xtable)
  library(knitr)
  library(kableExtra)
  library(gt)
  library(stargazer)
  library(corrplot)
})

# dplyr must win over stats for lag()/filter(): the original scripts relied on
# conflicted::conflicts_prefer(dplyr::lag / dplyr::filter / dplyr::select).
lag    <- dplyr::lag
filter <- dplyr::filter
select <- dplyr::select
first  <- dplyr::first

# ---------------------------------------------------------------------------
# Results folder layout (one sub-folder per subject)
# ---------------------------------------------------------------------------
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DIRS <- list(
  data        = file.path(RESULTS_DIR, "data"),
  descriptive = file.path(RESULTS_DIR, "descriptive_statistics"),
  unit_root   = file.path(RESULTS_DIR, "unit_root"),
  shocks      = file.path(RESULTS_DIR, "shock_diagnostics"),
  var         = file.path(RESULTS_DIR, "VAR"),
  lp          = file.path(RESULTS_DIR, "LP"),
  logs        = file.path(RESULTS_DIR, "logs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

# Write a LaTeX string (or character vector) to <dir>/<name>.tex
save_tex <- function(x, dir, name) {
  path <- file.path(dir, paste0(name, ".tex"))
  writeLines(as.character(x), path, useBytes = TRUE)
  message("  wrote ", sub(paste0(PROJECT_DIR, "/"), "", path))
  invisible(path)
}

# Print an xtable object straight into a .tex file
save_xtable <- function(xt, dir, name, ...) {
  path <- file.path(dir, paste0(name, ".tex"))
  print(xt, file = path, ...)
  message("  wrote ", sub(paste0(PROJECT_DIR, "/"), "", path))
  invisible(path)
}

# Save a ggplot / patchwork object in every format listed in FIG_FORMATS
save_plot <- function(plot, dir, name, width, height, dpi = 300) {
  for (fmt in FIG_FORMATS) {
    path <- file.path(dir, paste0(name, ".", fmt))
    if (fmt == "pdf") {
      ggsave(path, plot, width = width, height = height, device = cairo_pdf)
    } else {
      ggsave(path, plot, width = width, height = height, dpi = dpi)
    }
  }
  message("  wrote ", name, " [", paste(FIG_FORMATS, collapse = ", "), "]")
  invisible(NULL)
}

# Save a base-graphics plot (plot.ts, corrplot, strucchange, urca ...)
save_base_plot <- function(expr, dir, name, width = 8, height = 5, res = 300) {
  expr <- substitute(expr)
  env  <- parent.frame()
  for (fmt in FIG_FORMATS) {
    path <- file.path(dir, paste0(name, ".", fmt))
    if (fmt == "pdf") {
      cairo_pdf(path, width = width, height = height)
    } else {
      png(path, width = width, height = height, units = "in", res = res)
    }
    tryCatch(eval(expr, env), finally = dev.off())
  }
  invisible(NULL)
}

# Run code while diverting its console output to a log file
with_log <- function(path, expr) {
  con <- file(path, open = "wt")
  sink(con)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  })
  force(expr)
}

# Format a vector of ts-time break dates as "YYYY:Mmm" (as in thesis Table A2)
format_ym <- function(t) {
  yr <- floor(t + 1e-8)
  mo <- round((t - yr) * 12) + 1
  sprintf("%d:M%02d", yr, mo)
}
