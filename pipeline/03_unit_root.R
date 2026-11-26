###############################################################################
# 03_unit_root.R
# Port of Unit_root_testing_with_notes.R + unit_root_testing_table.R
# Thesis: Table A2 (Bai-Perron), Table A3 (ADF / Lee-Strazicich)
#
# Every series goes through the same steps as in the original scripts:
#   1. BIC lag selection for an AR(p) (ARMA_IC_selection_CSS, Pmax = 15)
#   2. AR(p) by CSS + Ljung-Box on the residuals
#   3. ADF sequence (Enders 2010): trend -> drift -> none, ADF lags = p - 1
#   4. Bai-Perron multiple breaks (h = 0.15)
#   5. Zivot-Andrews (model "both", lag 2) where the original ran it
#   6. Two-break Lee-Strazicich (GTOS, pn = 0.1)
# The specification of each series (AR order, ADF lags, last ADF type, ...)
# is copied from the original scripts.
###############################################################################
out      <- DIRS$unit_root
fig_dir  <- file.path(out, "figures")
log_dir  <- file.path(out, "test_output")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(log_dir, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "TS_functions.R"))
source(file.path(PROJECT_DIR, "additional_functions.R"))
source(file.path(PROJECT_DIR, "LeeStrazicichUnitRootTest.R"))

#=====================
# Prepare data
#=====================
# Common sample of Unit_root_testing_with_notes.R (2000:M2 onwards)
unit_root_data <- ea_data %>%
  select(Date, df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, Log_GPRA, Log_GPRT,
         Log_HICP, Log_Core_HICP, Log_Oil, Spliced_EONIA_3M, Log_VSTOXX, Log_EPU,
         df_Core_HICP, df_True_Policy_Rate, Spliced_EONIA_1Y, Unemployment_Rate,
         ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX,
         Target_Shock, FG_Shock, QE_Shock) %>%
  drop_na()

to_ts <- function(df, col) {
  d0 <- min(df$Date)
  ts(df[[col]], start = c(year(d0), month(d0)), frequency = 12)
}
# Series with a shorter history (PCCI from 2001, inflation swaps from 2008)
own_sample_ts <- function(col) {
  df <- ea_data %>% select(Date, all_of(col)) %>% drop_na()
  to_ts(df, col)
}

#=====================
# Test specification
#=====================
# label, series, AR order(s) for the Ljung-Box check, ADF lags, ADF type
# sequence (the last one is reported), Bai-Perron with trend?, run ZA?, LS model
ur_spec <- list(
  list(id = "HICP",        label = "HICP",                y = to_ts(unit_root_data, "Log_HICP"),         ar = 6,        adf_lags = 5,  adf_types = c("trend", "drift"),         bp_trend = TRUE,  za = TRUE,  ls = "break"),
  list(id = "core_HICP",   label = "Core HICP",           y = to_ts(unit_root_data, "Log_Core_HICP"),    ar = 7,        adf_lags = 6,  adf_types = c("trend", "drift", "none"), bp_trend = TRUE,  za = TRUE,  ls = "break"),
  list(id = "PCCI",        label = "PCCI",                y = own_sample_ts("PCCI"),                     ar = 11,       adf_lags = 10, adf_types = "trend",                     bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "IP",          label = "Ind. Production",     y = to_ts(unit_root_data, "Log_IP"),           ar = 3,        adf_lags = 2,  adf_types = "trend",                     bp_trend = TRUE,  za = TRUE,  ls = "break"),
  list(id = "Unemployment_Rate", label = "Unemployment",  y = to_ts(unit_root_data, "Unemployment_Rate"),ar = 8,        adf_lags = 7,  adf_types = c("trend", "drift", "none"), bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "ESI",         label = "ESI",                 y = to_ts(unit_root_data, "ESI"),              ar = 2,        adf_lags = 1,  adf_types = "trend",                     bp_trend = FALSE, za = FALSE, ls = "crash"),
  list(id = "Oil_prices",  label = "Oil Prices",          y = to_ts(unit_root_data, "Log_Oil"),          ar = 2,        adf_lags = 1,  adf_types = c("trend", "drift", "none"), bp_trend = FALSE, za = FALSE, ls = "crash"),
  list(id = "EONIA_1Y",    label = "EONIA / €STR 1Y (%)", y = to_ts(unit_root_data, "Spliced_EONIA_1Y"), ar = 3,        adf_lags = 2,  adf_types = c("trend", "drift", "none"), bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "True_Policy_Rate", label = "True Policy Rate", y = to_ts(unit_root_data, "True_Policy_Rate"), ar = c(4, 9), adf_lags = 8, adf_types = c("trend", "drift", "none"), bp_trend = FALSE, za = TRUE, ls = "crash"),
  list(id = "GPR",         label = "GPR",                 y = to_ts(unit_root_data, "Log_GPR"),          ar = 1,        adf_lags = 0,  adf_types = "trend",                     bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "VIX",         label = "VIX",                 y = to_ts(unit_root_data, "Log_VIX"),          ar = 1,        adf_lags = 0,  adf_types = "trend",                     bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "VSTOXX",      label = "VSTOXX",              y = to_ts(unit_root_data, "Log_VSTOXX"),       ar = 2,        adf_lags = 1,  adf_types = "trend",                     bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "EU_ILS_1Y",   label = "EU ILS 1Y",           y = own_sample_ts("EU_ILS_1Y"),                ar = 3,        adf_lags = 2,  adf_types = c("trend", "drift", "none"), bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "EU_ILS_5Y5Y", label = "EU ILS 5Y5Y",         y = own_sample_ts("EU_ILS_5Y5Y"),              ar = 2,        adf_lags = 1,  adf_types = c("trend", "drift", "none"), bp_trend = FALSE, za = TRUE,  ls = "crash"),
  list(id = "EPU",         label = "EPU",                 y = to_ts(unit_root_data, "Log_EPU"),          ar = 4,        adf_lags = 3,  adf_types = "trend",                     bp_trend = FALSE, za = TRUE,  ls = "crash")
)

#=====================
# Series plots
#=====================
for (s in ur_spec) {
  save_base_plot(plot.ts(s$y, main = s$label, ylab = s$label, xlab = "Time"),
                 fig_dir, paste0("series_", s$id), width = 8, height = 4.5)
}

#=====================
# Run the tests
#=====================
run_ur_tests <- function(s) {
  y <- s$y
  res <- list(spec = s)

  with_log(file.path(log_dir, paste0(s$id, ".txt")), {
    cat("=====================================================\n")
    cat(" Unit root testing:", s$label, "\n")
    cat(" Sample:", format_ym(time(y)[1]), "-", format_ym(tail(time(y), 1)),
        sprintf("(T = %d)\n", length(y)))
    cat("=====================================================\n\n")

    # 1. BIC lag selection
    cat("--- AR lag selection (CSS, Pmax = 15), sorted by BIC ---\n")
    ic <- ARMA_IC_selection_CSS(15, 0, y)
    print(head(ic, 5))
    res$bic_p <- ic$p[1]

    # 2. AR(p) + Ljung-Box on the residuals
    for (p in s$ar) {
      fit <- arima(y, order = c(p, 0, 0), method = "CSS",
                   include.mean = TRUE, transform.pars = FALSE)
      lb <- LjungBox(residuals(fit), 20, p)
      cat(sprintf("\n--- Ljung-Box on AR(%d) residuals ---\n", p))
      stargazer(lb, type = "text", summary = FALSE)
    }
    res$ar_used <- tail(s$ar, 1)
    res$lb_p20  <- as.character(lb$`p-value`[20])

    # 3. ADF sequence; the last specification is the one reported
    for (type in s$adf_types) {
      adf <- ur.df(y, type = type, lags = s$adf_lags)
      cat(sprintf("\n--- ADF(%d), type = %s ---\n", s$adf_lags, type))
      print(summary(adf))
    }
    res$adf      <- adf
    res$adf_type <- tail(s$adf_types, 1)

    # 4. Bai-Perron
    if (s$bp_trend) {
      time_trend <- time(y)
      bp <- breakpoints(y ~ time_trend, h = 0.15)
    } else {
      bp <- breakpoints(y ~ 1, h = 0.15)
    }
    cat("\n--- Bai-Perron (h = 0.15) ---\n")
    print(summary(bp))
    res$bp <- bp

    # 5. Zivot-Andrews
    if (s$za) {
      za <- ur.za(y, model = "both", lag = 2)
      cat("\n--- Zivot-Andrews (model = both, lag = 2) ---\n")
      print(summary(za))
      res$za <- za
    }

    # 6. Lee-Strazicich, two breaks
    cat(sprintf("\n--- Lee-Strazicich, 2 breaks, model = %s, GTOS ---\n", s$ls))
    res$ls <- ur.ls(y = y, model = s$ls, breaks = 2, method = "GTOS",
                    pn = 0.1, print.results = "print")
  })

  # Figures: Bai-Perron fit and Zivot-Andrews statistic path
  save_base_plot({
    plot(y, main = paste("Bai-Perron Structural Breaks:", s$label), ylab = s$label)
    lines(bp)
  }, fig_dir, paste0("bai_perron_", s$id), width = 8, height = 4.5)
  if (s$za) {
    save_base_plot(plot(res$za), fig_dir, paste0("zivot_andrews_", s$id),
                   width = 8, height = 4.5)
  }
  res
}

# The tests are deterministic, so running the series in parallel gives the
# same results; the two-break Lee-Strazicich grid search dominates run time.
message(sprintf("  running unit root tests for %d series on %d cores ...",
                length(ur_spec), N_CORES))
ur_results <- par_map(ur_spec, run_ur_tests)
names(ur_results) <- sapply(ur_spec, `[[`, "id")

#==========================================
# Tables (original xtable layout)
#==========================================
results_table <- bind_rows(lapply(ur_results, function(r) {
  bd <- breakdates(r$bp)
  n_breaks <- if (is.na(bd[1])) 0 else length(bd)
  ls_stat <- r$ls[[1]]$`t-stat`
  data.frame(
    Variable  = r$spec$label,
    ADF_Spec  = tools::toTitleCase(r$adf_type),
    ADF_Lags  = r$spec$adf_lags,
    ADF_Stat  = round(r$adf@teststat[1], 2),
    ADF_CV5   = r$adf@cval[1, "5pct"],
    LS_Stat   = ifelse(!is.null(ls_stat), round(ls_stat, 2), NA),
    LS_Break1 = format_ym(time(r$spec$y)[r$ls[[1]]$`First break`]),
    LS_Break2 = format_ym(time(r$spec$y)[r$ls[[1]]$`Second break`]),
    BP_Breaks = n_breaks,
    BP_Dates  = if (n_breaks == 0) "None" else paste(format_ym(bd), collapse = ", "),
    BP_Dates_Decimal = if (n_breaks == 0) "None" else paste(round(bd, 2), collapse = ", "),
    BIC_p     = r$bic_p,
    AR_p      = r$ar_used,
    LB_p20    = r$lb_p20,
    stringsAsFactors = FALSE
  )
}))
write_csv(results_table, file.path(out, "unit_root_results.csv"))

# Table A3: unit root tests
ur_table <- results_table %>%
  transmute(Variable = sprintf("%s (%s)", Variable, ADF_Spec),
            `ADF Statistic` = ADF_Stat,
            `LS Statistic`  = LS_Stat)
save_xtable(
  xtable(ur_table, caption = "Unit Root Test Results", label = "tab:unit_roots",
         align = "llcc"),
  out, "TableA3_Unit_Root_Tests", include.rownames = FALSE
)

# Detailed version with critical value, lags and LS break dates
ur_table_full <- results_table %>%
  transmute(Variable,
            `ADF spec.` = ADF_Spec,
            `ADF lags`  = ADF_Lags,
            `ADF stat.` = ADF_Stat,
            `ADF 5\\% c.v.` = ADF_CV5,
            `LS stat.`  = LS_Stat,
            `LS break 1` = LS_Break1,
            `LS break 2` = LS_Break2)
save_xtable(
  xtable(ur_table_full, caption = "Unit Root Test Results (detailed)",
         label = "tab:unit_roots_detailed", align = "llcrrrrcc", digits = 2),
  out, "TableA3_Unit_Root_Tests_Detailed", include.rownames = FALSE,
  sanitize.colnames.function = identity
)

# Table A2: Bai-Perron
bp_table <- results_table %>%
  transmute(Variable,
            `Opt. Breaks (BIC)` = BP_Breaks,
            `Identified Break Dates` = BP_Dates)
save_xtable(
  xtable(bp_table, caption = "Bai-Perron Multiple Structural Break Analysis",
         label = "tab:bai_perron", align = "llcl"),
  out, "TableA2_Bai_Perron", include.rownames = FALSE
)

# Lag selection overview
lag_table <- results_table %>%
  transmute(Variable,
            `BIC AR order` = BIC_p,
            `AR order used` = AR_p,
            `LB p-value (lag 20)` = LB_p20)
save_xtable(
  xtable(lag_table, caption = "AR Lag Selection and Ljung-Box Residual Check",
         label = "tab:ur_lag_selection", align = "llccc"),
  out, "Unit_Root_Lag_Selection", include.rownames = FALSE
)
