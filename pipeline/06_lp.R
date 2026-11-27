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

#==============================================================================
# Regression tables
#   For every response variable the horizon-h regression is
#     y_{t+h} - y_t = F_{t-1}[a_H + b_H e_t + sum_l g_H,l y_{t-l}]
#                   + (1 - F_{t-1})[a_L + b_L e_t + sum_l g_L,l y_{t-l}] + u_{t+h}
#   (level at h = 0). Standard errors are Newey-West with h lags.
#==============================================================================
reg_dir <- file.path(out, "regression_tables")
dir.create(reg_dir, showWarnings = FALSE)
names(lp_results_list) <- ynames

stars <- function(p) ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$",
                     ifelse(p < 0.10, "$^{*}$", "")))
tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  gsub("([%&_#$])", "\\\\\\1", x)
}

# Coefficients, NW standard errors and the High - Low difference per horizon
lp_coef_table <- function(var) {
  fits <- attr(lp_results_list[[var]], "fits")
  bind_rows(lapply(fits, function(f) {
    nm <- names(f$coef)
    get <- function(k) f$coef[[paste0("X_matrix", k)]]
    se  <- function(k) sqrt(f$vcov_nw[paste0("X_matrix", k), paste0("X_matrix", k)])
    iH <- "X_matrixshock_high"; iL <- "X_matrixshock_low"
    d    <- get("shock_high") - get("shock_low")
    d_se <- sqrt(f$vcov_nw[iH, iH] + f$vcov_nw[iL, iL] - 2 * f$vcov_nw[iH, iL])
    tibble(
      Variable = var, Horizon = f$horizon,
      b_high = get("shock_high"), se_high = se("shock_high"),
      b_low  = get("shock_low"),  se_low  = se("shock_low"),
      diff = d, se_diff = d_se,
      c_high = get("F_high"), se_c_high = se("F_high"),
      c_low  = get("F_low"),  se_c_low  = se("F_low"),
      N = f$nobs, R2 = f$r2, adj_R2 = f$adj_r2, NW_lag = f$nw_lag,
      n_controls = f$n_controls
    )
  })) %>%
    mutate(p_high = 2 * pnorm(-abs(b_high / se_high)),
           p_low  = 2 * pnorm(-abs(b_low / se_low)),
           p_diff = 2 * pnorm(-abs(diff / se_diff)),
           p_c_high = 2 * pnorm(-abs(c_high / se_c_high)),
           p_c_low  = 2 * pnorm(-abs(c_low / se_c_low)))
}
lp_coefs <- bind_rows(lapply(ynames, lp_coef_table))
write_csv(lp_coefs, file.path(reg_dir, "STLP_regression_coefficients.csv"))

fmt_coef <- function(b, p) sprintf("%.3f%s", b, stars(p))
fmt_se   <- function(s) sprintf("(%.3f)", s)

# (1) Stargazer-style table per variable: selected horizons as columns
sel_h <- c(0, 3, 6, 12, 18, 24)
lp_var_table <- function(var, label) {
  d <- lp_coefs %>% filter(Variable == var, Horizon %in% sel_h) %>% arrange(Horizon)
  rows <- list(
    c("Shock $\\times F_{t-1}$ (High GPR)", fmt_coef(d$b_high, d$p_high)),
    c("",                                    fmt_se(d$se_high)),
    c("Shock $\\times (1-F_{t-1})$ (Low GPR)", fmt_coef(d$b_low, d$p_low)),
    c("",                                    fmt_se(d$se_low)),
    c("Difference High $-$ Low",             fmt_coef(d$diff, d$p_diff)),
    c("",                                    fmt_se(d$se_diff)),
    c("$F_{t-1}$ (High-GPR constant)",       fmt_coef(d$c_high, d$p_c_high)),
    c("",                                    fmt_se(d$se_c_high)),
    c("$1-F_{t-1}$ (Low-GPR constant)",      fmt_coef(d$c_low, d$p_c_low)),
    c("",                                    fmt_se(d$se_c_low)),
    c(sprintf("Lags of %s $\\times$ regime", tex_escape(var)),
      rep(sprintf("%d $\\times$ 2", p), length(sel_h))),
    c("Observations",                        as.character(d$N)),
    c("$R^2$",                               sprintf("%.3f", d$R2)),
    c("Adjusted $R^2$",                      sprintf("%.3f", d$adj_R2)),
    c("Newey-West lags",                     as.character(d$NW_lag))
  )
  tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(tab) <- c("", paste0("$h=", sel_h, "$"))
  kbl(tab, format = "latex", booktabs = TRUE, escape = FALSE,
      align = c("l", rep("c", length(sel_h))), linesep = "",
      caption = sprintf("State-dependent local projection: response of %s to a pure MP shock",
                        tex_escape(var)),
      label = paste0("stlp_", label)) %>%
    kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 9) %>%
    add_header_above(c(" " = 1, "Horizon (months after the shock)" = length(sel_h))) %>%
    row_spec(10, hline_after = TRUE) %>%
    footnote(general = paste(
      "Dependent variable: $y_{t+h}-y_t$ for $h>0$ and $y_t$ for $h=0$.",
      "Shock: AR(1)-purged pure monetary policy shock extracted from the linear JK (2020) BVAR.",
      "$F_{t-1}$ is the logistic transition weight of log GPR ($c=10$, threshold at the median).",
      "Newey-West standard errors in parentheses;",
      "$^{***}$, $^{**}$, $^{*}$: significant at 1\\\\%, 5\\\\%, 10\\\\%."),
      threeparttable = TRUE, escape = FALSE, general_title = "Notes:")
}
lp_labels <- c("eonia_1y", "ind_production", "oil_prices", "hicp")
for (i in seq_along(ynames)) {
  save_tex(lp_var_table(ynames[i], lp_labels[i]), reg_dir,
           paste0("STLP_Regression_", c("EONIA_1Y", "Ind_Production", "Oil_Prices", "HICP")[i]))
}

# (2) All horizons in one long table
all_h <- lp_coefs %>%
  transmute(Variable = tex_escape(Variable), `$h$` = Horizon,
            `$\\beta_H$` = fmt_coef(b_high, p_high), `s.e.$_H$` = sprintf("%.3f", se_high),
            `$\\beta_L$` = fmt_coef(b_low, p_low),  `s.e.$_L$` = sprintf("%.3f", se_low),
            `$\\beta_H-\\beta_L$` = fmt_coef(diff, p_diff), `s.e.$_{H-L}$` = sprintf("%.3f", se_diff),
            N = N, `$R^2$` = sprintf("%.3f", R2))
save_tex(
  kbl(all_h, format = "latex", booktabs = TRUE, escape = FALSE, longtable = TRUE,
      align = "lrcccccccc", linesep = "",
      caption = "State-dependent local projections: shock coefficients at all horizons",
      label = "stlp_all_horizons") %>%
    kable_styling(latex_options = c("repeat_header"), font_size = 8) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "top") %>%
    footnote(general = paste(
      "$\\\\beta_H$ and $\\\\beta_L$: response in the high- and low-GPR regime;",
      "Newey-West standard errors with $h$ lags;",
      "$^{***}$, $^{**}$, $^{*}$: significant at 1\\\\%, 5\\\\%, 10\\\\%."),
      escape = FALSE, general_title = "Notes:"),
  reg_dir, "STLP_Regression_All_Horizons"
)
