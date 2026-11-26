###############################################################################
# 04_shock_diagnostics.R
# Port of shock_analysis_and_purging.R
# Thesis: Table A4 (ADF / ARCH / exogeneity), Table A5 (Ljung-Box), Figure A1
###############################################################################
out     <- DIRS$shocks
log_dir <- file.path(out, "test_output")
dir.create(log_dir, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "TS_functions.R"))
source(file.path(PROJECT_DIR, "additional_functions.R"))

shock_data <- ea_data %>%
  select(Date, Log_IP, Log_HICP, Log_Oil, Spliced_EONIA_1Y, OIS_3M, OIS_1Y,
         STOXX50, Target_Shock, FG_Shock, QE_Shock) %>%
  drop_na()

mk <- function(col) ts(shock_data[[col]], start = c(2000, 2), frequency = 12)
OIS_3M       <- mk("OIS_3M")
OIS_1Y       <- mk("OIS_1Y")
STOXX50      <- mk("STOXX50")
Target_Shock <- mk("Target_Shock")
FG_Shock     <- mk("FG_Shock")
QE_Shock     <- mk("QE_Shock")
HICP         <- mk("Log_HICP")
IP           <- mk("Log_IP")
Oil_prices   <- mk("Log_Oil")
EONIA_1Y     <- mk("Spliced_EONIA_1Y")

ts_data <- ts.union(OIS_3M, OIS_1Y, STOXX50, Target_Shock, FG_Shock, QE_Shock,
                    HICP, IP, Oil_prices, EONIA_1Y)

# AR order (BIC, adjusted for serial correlation) from the original script;
# the ADF uses p - 1 lags (with a floor of 0)
shock_spec <- list(
  OIS_3M       = list(label = "OIS 3M",       ar = 2, adf_lags = 1),
  OIS_1Y       = list(label = "OIS 1Y",       ar = 2, adf_lags = 1),
  STOXX50      = list(label = "STOXX50",      ar = 0, adf_lags = 0),
  Target_Shock = list(label = "Target Shock", ar = 3, adf_lags = 2),
  FG_Shock     = list(label = "FG Shock",     ar = 4, adf_lags = 3),
  QE_Shock     = list(label = "QE Shock",     ar = 0, adf_lags = 0)
)

exo_formula_rhs <- "L(IP, 1) + L(IP, 2) + L(HICP, 1) + L(HICP, 2) + L(Oil_prices, 1) + L(Oil_prices, 2) + L(EONIA_1Y, 1) + L(EONIA_1Y, 2)"

shock_res <- list()
for (nm in names(shock_spec)) {
  s <- shock_spec[[nm]]
  y <- get(nm)
  message("  shock diagnostics: ", s$label)
  r <- list()
  with_log(file.path(log_dir, paste0(nm, ".txt")), {
    cat("=====================================================\n")
    cat(" Shock diagnostics:", s$label, "\n")
    cat("=====================================================\n\n")
    cat("--- AR lag selection (CSS, Pmax = 15), sorted by BIC ---\n")
    print(head(ARMA_IC_selection_CSS(15, 0, y), 5))

    fit <- arima(y, order = c(s$ar, 0, 0), method = "CSS",
                 include.mean = TRUE, transform.pars = FALSE)
    r$ar_resid <- residuals(fit)
    cat(sprintf("\n--- Ljung-Box on AR(%d) residuals ---\n", s$ar))
    stargazer(LjungBox(r$ar_resid, 20, s$ar), type = "text", summary = FALSE)

    r$adf <- ur.df(y, type = "trend", lags = s$adf_lags)
    cat(sprintf("\n--- ADF(%d), type = trend ---\n", s$adf_lags))
    print(summary(r$adf))

    r$lb <- LjungBox(y, 12, 0)
    cat("\n--- Ljung-Box on the raw series (lags 1-12) ---\n")
    stargazer(r$lb, type = "text", summary = FALSE)

    r$arch <- ArchTest(y, lags = 6, demean = TRUE)
    cat("\n--- ARCH-LM (6 lags) ---\n")
    print(r$arch)

    exo <- dynlm(as.formula(paste(nm, "~", exo_formula_rhs)), data = ts_data)
    fs  <- summary(exo)$fstatistic
    r$exo_p <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
    cat("\n--- Exogeneity regression on lagged IP, HICP, Oil, EONIA 1Y ---\n")
    print(summary(exo))
    cat(s$label, "exogeneity p-value:", r$exo_p, "\n")
  })
  shock_res[[nm]] <- r
}

#==================================
# Table A4: main diagnostic tests
#==================================
summary_table <- data.frame(
  Variable = sapply(shock_spec, `[[`, "label"),
  ADF = sapply(shock_res, function(r) r$adf@teststat[1]),
  ARCH = round(sapply(shock_res, function(r) r$arch$p.value), 4),
  Exogeneity = round(sapply(shock_res, function(r) r$exo_p), 4),
  row.names = NULL
)
colnames(summary_table) <- c("Variable", "ADF (H0: unit root)", "ARCH (p-value)",
                             "Exogeneity (p-value)")
write_csv(summary_table, file.path(out, "shock_diagnostics.csv"))
save_xtable(
  xtable(summary_table,
         caption = "Diagnostic Tests for High-Frequency Monetary Policy Surprises",
         label = "tab:shock_diagnostics", align = "llccc", digits = c(0, 0, 3, 4, 4)),
  out, "TableA4_Shock_Diagnostics", include.rownames = FALSE
)

#==================================
# Table A5: Ljung-Box, lags 1-12
#==================================
lb_panel <- function(ids, panel_title) {
  wide <- data.frame(Lag = 1:12)
  for (nm in ids) {
    lb <- shock_res[[nm]]$lb
    wide[[paste0(shock_spec[[nm]]$label, " Stat")]]    <- sprintf("%.3f", lb$`Ljung-Box`)
    wide[[paste0(shock_spec[[nm]]$label, " p-value")]] <- ifelse(
      lb$`p-value` == "<0.0001", "$<$0.0001", sprintf("%.4f", as.numeric(lb$`p-value`)))
  }
  hdr <- c(" " = 1, setNames(rep(2, length(ids)), sapply(shock_spec[ids], `[[`, "label")))
  kbl(wide, format = "latex", booktabs = TRUE, escape = FALSE, align = "c",
      col.names = c("Lag", rep(c("Stat", "p-value"), length(ids))),
      caption = paste("Ljung-Box Test for Serial Correlation (Lags 1--12):", panel_title)) %>%
    add_header_above(hdr) %>%
    kable_styling(latex_options = "hold_position", font_size = 9)
}
save_tex(c(lb_panel(c("OIS_3M", "OIS_1Y", "STOXX50"), "Panel A: OIS Rates and Equity Proxy"),
           "",
           lb_panel(c("Target_Shock", "FG_Shock", "QE_Shock"), "Panel B: Structural Monetary Policy Shocks")),
         out, "TableA5_Ljung_Box")

#==================================
# Figure A1: correlation of shocks
#==================================
shocks_matrix <- cbind(Target_Shock, FG_Shock, QE_Shock, OIS_3M, OIS_1Y, STOXX50)
colnames(shocks_matrix) <- c("Target", "FG", "QE", "OIS 3M", "OIS 1Y", "STOXX50")
ortho_cor_matrix <- cor(shocks_matrix, use = "complete.obs")

save_base_plot(
  corrplot(ortho_cor_matrix, method = "number", type = "upper",
           tl.col = "black", tl.srt = 60,
           title = "Orthogonality of Monetary Policy Surprises", mar = c(0, 0, 1, 0)),
  out, "FigureA1_Shock_Correlations", width = 7, height = 7
)
save_tex(
  kbl(ortho_cor_matrix, format = "latex", booktabs = TRUE, digits = 3,
      caption = "Correlations of Monetary Policy Surprises") %>%
    kable_styling(latex_options = "hold_position"),
  out, "FigureA1_Shock_Correlations_table"
)
