###############################################################################
# 06_lp.R
# Port of ST_LP_main.R
# Thesis: Figure 3 (State-dependent LP with the extracted pure MP shock)
#
# Uses the median structural pure MP shock of the linear JK (2020) BVAR from
# stage 05. If that stage was skipped, the shock saved by an earlier run in
# results/VAR/baseline_bvar_structural_mp_shock.csv is used instead.
###############################################################################
out     <- DIRS$lp
log_dir <- file.path(out, "test_output")
dir.create(log_dir, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "TS_functions.R"))
source(file.path(PROJECT_DIR, "additional_functions.R"))
source(file.path(PROJECT_DIR, "ST_LP_functions.R"))

if (!exists("var_data_eff")) {
  shock_file <- file.path(DIRS$var, "baseline_bvar_structural_mp_shock.csv")
  if (!file.exists(shock_file)) {
    stop("No structural MP shock available: run the VAR stage first.")
  }
  message("  using structural MP shock from ", basename(shock_file))
  var_data_eff <- read_csv(shock_file, show_col_types = FALSE)
}

#----------------------------------------------
# Purge the MP shock of serial correlation
#----------------------------------------------
mp_shock <- ts(var_data_eff$`Pure MP Shock`, start = c(2001, 2), frequency = 12)
with_log(file.path(log_dir, "mp_shock_purging.txt"), {
  cat("--- Ljung-Box on the extracted pure MP shock (lags 1-12) ---\n")
  stargazer(LjungBox(mp_shock, 12, 0), type = "text", summary = FALSE)
  cat("\n--- AR lag selection (CSS, Pmax = 15), sorted by BIC ---\n")
  print(head(ARMA_IC_selection_CSS(15, 0, mp_shock), 5))
  AR_MP <- arima(mp_shock, order = c(1, 0, 0), method = "CSS",
                 include.mean = TRUE, transform.pars = FALSE)
  print(AR_MP)
  cat("\n--- Ljung-Box on AR(1) residuals ---\n")
  stargazer(LjungBox(residuals(AR_MP), 20, 1), type = "text", summary = FALSE)
})
var_data_eff$`Purged MP Shock` <- as.numeric(residuals(AR_MP))
save_base_plot(acf(mp_shock, main = "ACF of the extracted pure MP shock"),
               out, "MP_Shock_ACF", width = 7, height = 4.5)

#==================
# Settings
#==================
ynames     <- c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP")
mnames     <- "Purged MP Shock"
threshname <- "Log_GPR"
data_names <- c(mnames, ynames, threshname)
p          <- 12

Y_matrix <- as.matrix(var_data_eff[, data_names])

lp_results_list <- lapply(ynames, function(v) {
  run_st_lp(data = Y_matrix, target_var = v, shock_var = mnames,
            gpr_var = threshname, p_lags = p, horizon = 24)
})
lp_plot_data <- do.call(rbind, lp_results_list)

lp_long <- bind_rows(
  lp_plot_data %>%
    select(Horizon, Variable, Median = Beta_Low, Lower68 = Lower68_Low, Upper68 = Upper68_Low) %>%
    mutate(Regime = "Low GPR"),
  lp_plot_data %>%
    select(Horizon, Variable, Median = Beta_High, Lower68 = Lower68_High, Upper68 = Upper68_High) %>%
    mutate(Regime = "High GPR")
) %>%
  mutate(Variable = factor(Variable, levels = ynames))
write_csv(lp_long, file.path(out, "Figure3_STLP_irf_data.csv"))

T_eff <- nrow(Y_matrix) - p

p_lp <- ggplot(lp_long, aes(x = Horizon, color = Regime, fill = Regime)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  geom_line(aes(y = Median), linewidth = 1) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("High GPR" = "firebrick", "Low GPR" = "steelblue")) +
  scale_fill_manual(values = c("High GPR" = "firebrick", "Low GPR" = "steelblue")) +
  labs(
    title = "State-Dependent Local Projection: Response to a Pure Monetary Policy Shock",
    subtitle = "High vs. Low Geopolitical Risk Regimes",
    caption = sprintf("Effective Observations (T): %d | Lags (p): %d", T_eff, p),
    x = "Months after shock",
    y = "Impulse Response"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )
save_plot(p_lp, out, "Figure3_STLP_Pure_MP_Shock", width = 10, height = 7)

# LaTeX table of the LP coefficients at selected horizons
lp_tab <- lp_long %>%
  filter(Horizon %in% c(0, 3, 6, 12, 18, 24)) %>%
  mutate(Cell = sprintf("%.3f [%.3f, %.3f]", Median, Lower68, Upper68)) %>%
  select(Variable, Regime, Horizon, Cell) %>%
  pivot_wider(names_from = Horizon, values_from = Cell, names_prefix = "h=") %>%
  arrange(Variable, desc(Regime))
save_tex(
  kbl(lp_tab, format = "latex", booktabs = TRUE, align = "llcccccc",
      caption = "State-Dependent LP responses to a pure MP shock (point estimate and 68\\% Newey-West band)",
      escape = TRUE) %>%
    kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 8) %>%
    collapse_rows(columns = 1, latex_hline = "major"),
  out, "Figure3_STLP_Coefficients"
)
