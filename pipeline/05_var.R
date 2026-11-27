###############################################################################
# 05_var.R
# Ports of jk_replication_main.R, jk_SD_BVAR_main.R and chol_ident_SD_BVAR_main.R
# Thesis: Figure 1 (baseline ST-BVAR), Figure 2 (ST-BVAR, HFI + Cholesky),
#         Figure C1 (linear JK-2020 BVAR) and Figures C2-C21 (robustness)
#
# Every model gets its own fixed RNG seed, so results are reproducible and do
# not depend on the order in which models are run or on the number of cores.
###############################################################################
out      <- DIRS$var
irf_dir  <- file.path(out, "irf_data")
log_dir  <- file.path(out, "logs")
dir.create(irf_dir, showWarnings = FALSE)
dir.create(log_dir, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "jk_replication_function.R"))  # VAR_withiid, resirfssign, ...
source(file.path(PROJECT_DIR, "jk_SD_BVAR_functions.R"))     # VAR_withiid_weighted

#======================================
# Variable labels used in the VARs
#======================================
var_labels <- c(
  `Ind. Production`              = "Log_IP",
  `Diff. Oil Prices`             = "df_Oil",
  `True Policy Rate`             = "True_Policy_Rate",
  `HICP`                         = "Log_HICP",
  `Core HICP`                    = "Log_Core_HICP",
  `Diff. HICP`                   = "df_HICP",
  `Oil Prices`                   = "Log_Oil",
  `EONIA/€STR 3M`                = "Spliced_EONIA_3M",
  `EONIA/€STR 1Y`                = "Spliced_EONIA_1Y",
  `EONIA/€STR 2Y`                = "Spliced_EONIA_2Y",
  `Diff. EONIA/€STR 1Y`          = "df_EONIA_1Y",
  `ESI`                          = "ESI",
  `PCCI`                         = "PCCI",
  `Diff. Unemployment Rate`      = "df_Unemployment_Rate",
  `Diff. Ind. Prod.`             = "df_IP",
  `1Y Infl. Linked Swap`         = "EU_ILS_1Y",
  `5Y5Y Infl. Linked Swap`       = "EU_ILS_5Y5Y",
  `Diff. 1Y Infl. Linked Swap`   = "df_ILS_1Y",
  `Diff. 5Y5Y Infl. Linked Swap` = "df_ILS_5Y5Y",
  `OIS 1Y surprise`              = "OIS_1Y",
  `EPU`                          = "Log_EPU",
  OIS_3M = "OIS_3M", STOXX50 = "STOXX50",
  Log_GPR = "Log_GPR", Log_GPRA = "Log_GPRA", Log_GPRT = "Log_GPRT",
  Log_VIX = "Log_VIX", Log_VSTOXX = "Log_VSTOXX", Log_EPU = "Log_EPU"
)

# Build the estimation sample: rename columns to their labels and drop rows
# with a missing value in any of the `sample_vars` (as the drop_na() in the
# original scripts did on their select() lists).
build_var_data <- function(sample_vars) {
  sample_vars <- unique(sample_vars)
  cols <- var_labels[sample_vars]
  ea_data %>%
    select(Date, all_of(setNames(cols, sample_vars))) %>%
    drop_na()
}

# Variables of the original select() lists (they define the common sample)
base_vars_linear <- c("Ind. Production", "Diff. Oil Prices", "True Policy Rate",
                      "Log_GPR", "Log_GPRA", "Log_GPRT", "HICP", "Core HICP",
                      "Oil Prices", "EONIA/€STR 3M", "ESI", "OIS_3M",
                      "OIS 1Y surprise", "STOXX50", "Log_VIX", "EONIA/€STR 1Y",
                      "EPU", "Diff. Unemployment Rate", "Diff. HICP")
base_vars_chol   <- setdiff(base_vars_linear, "Diff. HICP")
base_vars_st     <- c("Ind. Production", "Diff. Oil Prices", "True Policy Rate",
                      "Log_GPR", "Log_GPRA", "Log_GPRT", "HICP", "Core HICP",
                      "Diff. HICP", "Oil Prices", "EONIA/€STR 3M", "ESI", "OIS_3M",
                      "OIS 1Y surprise", "STOXX50", "Log_VIX", "EONIA/€STR 1Y",
                      "Log_VSTOXX", "Diff. EONIA/€STR 1Y", "Log_EPU",
                      "Diff. Unemployment Rate", "Diff. Ind. Prod.")

#======================================
# Common settings (thesis Section 3.3)
#======================================
mnames <- c("OIS_3M", "STOXX50")
prior_base <- list(
  lags = 12,
  minnesota_tightness = 0.2,
  minnesota_decay = 1,
  Nm = length(mnames),
  minnesota_mvector = c(0, 0, 0.9, 0.9, 0.9, 0.9, 0.9)  # 0 for surprises, 0.9 for levels
)
gssettings <- list(
  ndraws = BVAR_NDRAWS,
  burnin = BVAR_BURNIN,
  saveevery = 4,
  computemarglik = 0
)
MAlags  <- 36
max_try <- 1000
imonpol <- 1
inews   <- 2
dims    <- c(1, 2)
# JK (2020) "sgnm2" sign restrictions on impact
test_restr <- function(irfs) {
  (irfs[1, imonpol, 1] > 0) && (irfs[2, imonpol, 1] < 0) &&
    (irfs[1, inews, 1] > 0)   && (irfs[2, inews, 1] > 0)
}

# Smooth-transition weights F(z_{t-d}) on the effective sample
st_weights <- function(thresh, P, d = 1, c = 10, gamma) {
  z_t <- thresh[(P + 1 - d):(length(thresh) - d)]
  F_high <- 1 / (1 + exp(-c * (z_t - gamma) / sd(z_t)))
  list(F_high = F_high, F_low = 1 - F_high)
}

regime_cols <- c("Low GPR" = "steelblue", "High GPR" = "firebrick")

plot_st_irf <- function(plot_data, ynames, title, subtitle, T_eff, P) {
  ggplot(plot_data, aes(x = Horizon, y = Median, color = Regime, fill = Regime)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = regime_cols) +
    scale_fill_manual(values = regime_cols) +
    labs(title = title, subtitle = subtitle,
         caption = sprintf("Effective Observations (T): %d | Lags (p): %d", T_eff, P),
         x = "Months after Shock", y = "Impulse Response") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      plot.caption = element_text(face = "italic", size = 10, color = "gray30")
    )
}

regime_irfs <- function(irf_low, irf_high, data_names, ynames) {
  df_low  <- extract_matlab_irfs(irf_low,  data_names, shock_idx = 1)
  df_high <- extract_matlab_irfs(irf_high, data_names, shock_idx = 1)
  df_low$Regime  <- "Low GPR"
  df_high$Regime <- "High GPR"
  rbind(df_low, df_high) %>%
    filter(Variable %in% ynames) %>%
    mutate(Variable = factor(Variable, levels = ynames),
           Horizon = as.numeric(Horizon))
}

#==============================================================================
# Posterior summaries for the regression tables (computed inside each worker,
# uses no random numbers, so it does not affect the estimates)
#==============================================================================
posterior_summary <- function(res, irfs, data_names, weights = NULL, T_eff) {
  N <- length(data_names)
  P <- res$lags
  # rows of beta: lag 1 (all variables), ..., lag P, constant
  regressors <- c(paste0(rep(data_names, P), " (t-", rep(1:P, each = N), ")"), "Constant")
  q <- function(a, pr) apply(a, c(1, 2), quantile, probs = pr, na.rm = TRUE)
  bd <- res$beta_draws
  coef <- list(median = apply(bd, c(1, 2), median), sd = apply(bd, c(1, 2), sd),
               q05 = q(bd, 0.05), q16 = q(bd, 0.16), q84 = q(bd, 0.84), q95 = q(bd, 0.95))
  coef <- lapply(coef, function(m) { dimnames(m) <- list(regressors, data_names); m })
  sigma_med <- apply(res$sigma_draws, c(1, 2), median)
  dimnames(sigma_med) <- list(data_names, data_names)
  imp <- irfs[, 1:2, 1, , drop = FALSE]          # impact responses to shocks 1 and 2
  impact <- lapply(c(q16 = 0.16, median = 0.5, q84 = 0.84), function(pr) {
    m <- apply(imp, c(1, 2), quantile, probs = pr, na.rm = TRUE)
    dimnames(m) <- list(data_names, c("Shock 1", "Shock 2")); m
  })
  list(coef = coef, sigma = sigma_med, impact = impact, T_eff = T_eff,
       eff_obs = if (is.null(weights)) T_eff else sum(weights),
       n_draws = sum(!is.na(irfs[1, 1, 1, ])))
}

#==============================================================================
# 1) Linear JK (2020) BVAR, sample 2000-2025 (Figure C1)
#    Also delivers the structural pure MP shock used in the local projections.
#==============================================================================
run_linear_bvar <- function(seed) {
  set.seed(seed)
  ynames <- c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP")
  data_names <- c(mnames, ynames)
  N <- length(data_names)
  var_data <- build_var_data(c(base_vars_linear, data_names))
  Y_matrix <- as.matrix(var_data[, data_names])
  prior <- prior_base

  res <- VAR_withiid(Y_matrix, prior, gssettings)
  irfs_draws <- resirfssign(res, MAlags, dims, test_restr, rep(1, N), max_try)

  plot_data <- extract_matlab_irfs(irfs_draws, data_names, shock_idx = 1) %>%
    filter(Variable %in% ynames) %>%
    mutate(Variable = factor(Variable, levels = ynames))
  T_eff <- nrow(Y_matrix) - prior$lags

  # Structural shocks: epsilon_t = A0^{-1} u_t, median over posterior draws
  P <- prior$lags
  T_full <- nrow(Y_matrix)
  Y_t <- Y_matrix[(P + 1):T_full, , drop = FALSE]
  X_t <- matrix(1, nrow = T_eff, ncol = 1)
  for (l in 1:P) X_t <- cbind(X_t, Y_matrix[(P + 1 - l):(T_full - l), , drop = FALSE])
  X_t_reordered <- cbind(X_t[, -1, drop = FALSE], X_t[, 1, drop = FALSE])
  nd <- dim(res$beta_draws)[3]
  shock_draws <- matrix(NA, nrow = T_eff, ncol = nd)
  for (i in 1:nd) {
    A0_Q <- irfs_draws[, , 1, i]
    if (any(is.na(A0_Q))) next
    U_i <- Y_t - X_t_reordered %*% res$beta_draws[, , i]
    shock_draws[, i] <- t(solve(A0_Q, t(U_i)))[, 1]
  }
  var_data_eff <- var_data[(P + 1):nrow(var_data), ]
  var_data_eff$`Pure MP Shock` <- apply(shock_draws, 1, median, na.rm = TRUE)

  list(plot_data = plot_data, T_eff = T_eff, P = P, var_data_eff = var_data_eff,
       posterior = list(Linear = posterior_summary(res, irfs_draws, data_names, T_eff = T_eff)),
       shock_names = c("Pure MP", "CB information"))
}

#==============================================================================
# 2) Smooth-transition BVARs, HFI + sign restrictions (Figure 1, C2-C21)
#==============================================================================
yn <- list(
  base      = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP"),
  core      = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "Core HICP"),
  pcci      = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "PCCI"),
  esi       = c("EONIA/€STR 1Y", "ESI", "Oil Prices", "HICP"),
  unemp     = c("EONIA/€STR 1Y", "Diff. Unemployment Rate", "Oil Prices", "HICP"),
  tpr       = c("True Policy Rate", "Ind. Production", "Oil Prices", "HICP"),
  eonia3m   = c("EONIA/€STR 3M", "Ind. Production", "Oil Prices", "HICP"),
  eonia2y   = c("EONIA/€STR 2Y", "Ind. Production", "Oil Prices", "HICP"),
  ils5y5y   = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP", "5Y5Y Infl. Linked Swap"),
  ils1y     = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP", "1Y Infl. Linked Swap"),
  dils5y5y  = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP", "Diff. 5Y5Y Infl. Linked Swap"),
  dils1y    = c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP", "Diff. 1Y Infl. Linked Swap"),
  diffs     = c("Diff. EONIA/€STR 1Y", "Diff. Ind. Prod.", "Diff. Oil Prices", "Diff. HICP")
)

# name, figure file, y-variables, threshold, gamma quantile (NULL = median), c, title
st_models <- list(
  list(file = "Figure1_STBVAR_Baseline",            y = yn$base,     thresh = "Log_GPR",    c = 10, gq = NULL, desc = "baseline"),
  list(file = "FigureC2_STBVAR_ESI",                y = yn$esi,      thresh = "Log_GPR",    c = 10, gq = NULL, desc = "Economic Sentiment Index"),
  list(file = "FigureC3_STBVAR_Unemployment",       y = yn$unemp,    thresh = "Log_GPR",    c = 10, gq = NULL, desc = "unemployment rate"),
  list(file = "FigureC4_STBVAR_Core_HICP",          y = yn$core,     thresh = "Log_GPR",    c = 10, gq = NULL, desc = "Core HICP"),
  list(file = "FigureC5_STBVAR_PCCI",               y = yn$pcci,     thresh = "Log_GPR",    c = 10, gq = NULL, desc = "PCCI"),
  list(file = "FigureC6_STBVAR_True_Policy_Rate",   y = yn$tpr,      thresh = "Log_GPR",    c = 10, gq = NULL, desc = "True Policy Rate"),
  list(file = "FigureC7_STBVAR_EONIA_3M",           y = yn$eonia3m,  thresh = "Log_GPR",    c = 10, gq = NULL, desc = "3-month EONIA/€STR"),
  list(file = "FigureC8_STBVAR_EONIA_2Y",           y = yn$eonia2y,  thresh = "Log_GPR",    c = 10, gq = NULL, desc = "2-year EONIA/€STR"),
  list(file = "FigureC9_STBVAR_ILS_1Y",             y = yn$ils1y,    thresh = "Log_GPR",    c = 10, gq = NULL, desc = "1Y inflation linked swap"),
  list(file = "FigureC10_STBVAR_ILS_5Y5Y",          y = yn$ils5y5y,  thresh = "Log_GPR",    c = 10, gq = NULL, desc = "5Y5Y inflation linked swap"),
  list(file = "FigureC11_STBVAR_First_Differences", y = yn$diffs,    thresh = "Log_GPR",    c = 10, gq = NULL, desc = "first differences"),
  list(file = "FigureC12_STBVAR_Diff_ILS_1Y",       y = yn$dils1y,   thresh = "Log_GPR",    c = 10, gq = NULL, desc = "diff. 1Y inflation swap"),
  list(file = "FigureC13_STBVAR_Diff_ILS_5Y5Y",     y = yn$dils5y5y, thresh = "Log_GPR",    c = 10, gq = NULL, desc = "diff. 5Y5Y inflation swap"),
  list(file = "FigureC14_STBVAR_Threshold_GPRT",    y = yn$base,     thresh = "Log_GPRT",   c = 10, gq = NULL, desc = "GPRT threshold"),
  list(file = "FigureC15_STBVAR_Threshold_GPRA",    y = yn$base,     thresh = "Log_GPRA",   c = 10, gq = NULL, desc = "GPRA threshold"),
  list(file = "FigureC16_STBVAR_Threshold_VIX",     y = yn$base,     thresh = "Log_VIX",    c = 10, gq = NULL, desc = "VIX threshold"),
  list(file = "FigureC17_STBVAR_Threshold_VSTOXX",  y = yn$base,     thresh = "Log_VSTOXX", c = 10, gq = NULL, desc = "VSTOXX threshold"),
  list(file = "FigureC18_STBVAR_Threshold_EPU",     y = yn$base,     thresh = "Log_EPU",    c = 10, gq = NULL, desc = "EPU threshold"),
  list(file = "FigureC19_STBVAR_c5",                y = yn$base,     thresh = "Log_GPR",    c = 5,  gq = NULL, desc = "transition speed c = 5"),
  list(file = "FigureC20_STBVAR_c20",               y = yn$base,     thresh = "Log_GPR",    c = 20, gq = NULL, desc = "transition speed c = 20"),
  list(file = "FigureC21_STBVAR_Gamma70",           y = yn$base,     thresh = "Log_GPR",    c = 10, gq = 0.7,  desc = "GPR threshold at 70%")
)
st_models <- Filter(function(m) is.null(VAR_MODELS) || m$file %in% VAR_MODELS, st_models)

run_st_bvar <- function(m, seed) {
  set.seed(seed)
  data_names <- c(mnames, m$y)
  N <- length(data_names)
  var_data <- build_var_data(c(base_vars_st, data_names, m$thresh))
  Y_matrix <- as.matrix(var_data[, data_names])
  prior <- prior_base
  P <- prior$lags

  thresh_var <- var_data[[m$thresh]]
  gamma <- if (is.null(m$gq)) median(thresh_var, na.rm = TRUE) else
    quantile(thresh_var, m$gq, na.rm = TRUE)
  w <- st_weights(thresh_var, P, d = 1, c = m$c, gamma = gamma)

  cat("\nEstimating Low GPR Regime...\n")
  res_low  <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = w$F_low)
  cat("\nEstimating High GPR Regime...\n")
  res_high <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = w$F_high)

  cat("\nApplying Sign Restrictions - Low GPR Regime...\n")
  irf_low  <- resirfssign(res_low,  MAlags, dims, test_restr, rep(1, N), max_try)
  cat("\nApplying Sign Restrictions - High GPR Regime...\n")
  irf_high <- resirfssign(res_high, MAlags, dims, test_restr, rep(1, N), max_try)

  T_eff <- nrow(Y_matrix) - P
  list(plot_data = regime_irfs(irf_low, irf_high, data_names, m$y),
       T_eff = T_eff, P = P,
       sample = range(var_data$Date), gamma = unname(gamma),
       posterior = list(
         `Low GPR`  = posterior_summary(res_low,  irf_low,  data_names, w$F_low,  T_eff),
         `High GPR` = posterior_summary(res_high, irf_high, data_names, w$F_high, T_eff)),
       shock_names = c("Pure MP", "CB information"))
}

#==============================================================================
# 3) Smooth-transition BVAR with standard HFI + Cholesky (Figure 2)
#==============================================================================
run_chol_st_bvar <- function(seed) {
  set.seed(seed)
  ynames <- c("EONIA/€STR 1Y", "Ind. Production", "Oil Prices", "HICP")
  data_names <- c(mnames, ynames)
  N <- length(data_names)
  var_data <- build_var_data(c(base_vars_chol, data_names))
  Y_matrix <- as.matrix(var_data[, data_names])
  prior <- prior_base
  P <- prior$lags

  w <- st_weights(var_data$Log_GPR, P, d = 1, c = 10,
                  gamma = median(var_data$Log_GPR, na.rm = TRUE))

  cat("\nEstimating Low GPR Regime...\n")
  res_low  <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = w$F_low)
  cat("\nEstimating High GPR Regime...\n")
  res_high <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = w$F_high)

  chol_irfs <- function(res) {
    nd <- dim(res$beta_draws)[3]
    irfs <- array(NA, dim = c(N, N, MAlags, nd))
    for (i in 1:nd) {
      B_array <- array(t(res$beta_draws[1:(N * P), , i]), dim = c(N, N, P))
      irfs[, , , i] <- impulsdtrf(B_array, chol(res$sigma_draws[, , i]), MAlags)
    }
    irfs
  }
  cat("Computing State-Dependent Cholesky IRFs...\n")
  irf_low  <- chol_irfs(res_low)
  irf_high <- chol_irfs(res_high)
  T_eff <- nrow(Y_matrix) - P
  list(plot_data = regime_irfs(irf_low, irf_high, data_names, ynames),
       T_eff = T_eff, P = P,
       posterior = list(
         `Low GPR`  = posterior_summary(res_low,  irf_low,  data_names, w$F_low,  T_eff),
         `High GPR` = posterior_summary(res_high, irf_high, data_names, w$F_high, T_eff)),
       shock_names = c("OIS 3M (Cholesky, ordered first)", "STOXX50 (Cholesky, ordered second)"))
}

#==============================================================================
# Run everything (in parallel where possible)
#==============================================================================
jobs <- c(
  list(list(type = "linear", file = "FigureC1_Linear_BVAR_Baseline", seed = SEED)),
  lapply(seq_along(st_models), function(i)
    c(st_models[[i]], type = "st", seed = SEED + i)),
  list(list(type = "chol", file = "Figure2_STBVAR_HFI_Cholesky", seed = SEED + 100))
)
jobs <- Filter(function(j) is.null(VAR_MODELS) || j$file %in% VAR_MODELS ||
                 j$type == "linear", jobs)

run_job <- function(j) {
  t0 <- Sys.time()
  res <- with_log(file.path(log_dir, paste0(j$file, ".txt")), {
    switch(j$type,
           linear = run_linear_bvar(j$seed),
           st     = run_st_bvar(j, j$seed),
           chol   = run_chol_st_bvar(j$seed))
  })
  res$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  res
}

message(sprintf("  estimating %d BVAR models on %d cores (burn-in %d, %d draws)...",
                length(jobs), N_CORES, BVAR_BURNIN, BVAR_NDRAWS))
var_results <- par_map(jobs, run_job)
names(var_results) <- sapply(jobs, `[[`, "file")

#==============================================================================
# Figures and IRF data
#==============================================================================
model_info <- list()
for (j in jobs) {
  r <- var_results[[j$file]]
  write_csv(r$plot_data, file.path(irf_dir, paste0(j$file, ".csv")))

  if (j$type == "linear") {
    p <- ggplot(r$plot_data, aes(x = Horizon, y = Median)) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
      geom_ribbon(aes(ymin = Lower68, ymax = Upper68), fill = "#2c3e50", alpha = 0.2) +
      geom_line(color = "#2c3e50", linewidth = 1.2) +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      theme_minimal(base_size = 14) +
      labs(title = "Baseline BVAR(sample:2000-2025)",
           subtitle = "Impulse Responses to a 1-Standard-Deviation Pure MP Shock",
           caption = sprintf("Effective Observations (T): %d | Lags (p): %d", r$T_eff, r$P),
           x = "Months after Shock", y = "Response") +
      theme(plot.title = element_text(face = "bold", size = 16),
            strip.text = element_text(face = "bold", size = 12),
            plot.caption = element_text(face = "italic", size = 10, color = "gray30"),
            panel.grid.minor = element_blank())
    write_csv(r$var_data_eff, file.path(out, "baseline_bvar_structural_mp_shock.csv"))
    # exact binary copy, used when the LP stage is run on its own
    saveRDS(r$var_data_eff, file.path(out, "baseline_bvar_structural_mp_shock.rds"))
    desc <- "Linear JK (2020) BVAR, sign restrictions"
  } else if (j$type == "st") {
    p <- plot_st_irf(r$plot_data, j$y,
                     title = sprintf("State-Dependent BVAR: Pure MP Shock (%s)", j$desc),
                     subtitle = "HF & Sign Restriction Identification: High vs. Low Geopolitical Risk Regimes",
                     T_eff = r$T_eff, P = r$P)
    desc <- sprintf("ST-BVAR, %s (threshold %s, c = %g, gamma = %s)", j$desc,
                    sub("Log_", "", j$thresh), j$c,
                    if (is.null(j$gq)) "median" else paste0(100 * j$gq, "%"))
  } else {
    p <- plot_st_irf(r$plot_data, unique(as.character(r$plot_data$Variable)),
                     title = "State-Dependent Responses to a Monetary Policy Shock",
                     subtitle = "HFI Choleski Identification: High vs. Low Geopolitical Risk Regimes",
                     T_eff = r$T_eff, P = r$P)
    desc <- "ST-BVAR, standard HFI with Cholesky identification"
  }
  save_plot(p, out, j$file, width = 10, height = 7)
  model_info[[j$file]] <- data.frame(Figure = j$file, Model = desc,
                                     T = r$T_eff, Lags = r$P, Seed = j$seed,
                                     Minutes = round(r$minutes, 1))
}

model_overview <- bind_rows(model_info)
write_csv(model_overview, file.path(out, "model_overview.csv"))
save_xtable(
  xtable(model_overview %>% select(-Minutes) %>%
           mutate(across(c(Figure, Model), ~ gsub("([%&_#$])", "\\\\\\1", .x))),
         caption = sprintf("Estimated BVAR specifications (Gibbs sampler: %d burn-in, %d saved draws, every 4th draw kept; 68\\%% bands)",
                           BVAR_BURNIN, BVAR_NDRAWS),
         label = "tab:var_models", align = c("l", "l", "p{6cm}", "r", "r", "r"),
         digits = 0),
  out, "VAR_Model_Overview", include.rownames = FALSE,
  sanitize.text.function = identity, size = "\\scriptsize"
)

# Hand the structural shock to the LP stage
var_data_eff <- var_results[["FigureC1_Linear_BVAR_Baseline"]]$var_data_eff

#==============================================================================
# Regression tables: posterior reduced-form coefficients and impact responses
#==============================================================================
reg_dir <- file.path(out, "regression_tables")
dir.create(reg_dir, showWarnings = FALSE)

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  gsub("([%&_#$])", "\\\\\\1", x)
}
cred_mark <- function(lo68, hi68, lo90, hi90) {
  ifelse(lo90 > 0 | hi90 < 0, "$^{**}$", ifelse(lo68 > 0 | hi68 < 0, "$^{*}$", ""))
}

# Long table of posterior medians (posterior s.d.) for the macro equations
coef_table_tex <- function(ps, regime, title, label, Nm) {
  cf <- ps$coef
  eqs <- colnames(cf$median)[(Nm + 1):ncol(cf$median)]   # surprise equations are block-exogenous
  rows <- list()
  for (r in rownames(cf$median)) {
    est <- sprintf("%.3f%s", cf$median[r, eqs],
                   cred_mark(cf$q16[r, eqs], cf$q84[r, eqs], cf$q05[r, eqs], cf$q95[r, eqs]))
    rows[[length(rows) + 1]] <- c(tex_escape(r), est)
    rows[[length(rows) + 1]] <- c("", sprintf("(%.3f)", cf$sd[r, eqs]))
  }
  foot <- list(
    c("Observations ($T$)", rep(as.character(ps$T_eff), length(eqs))),
    c("Effective observations ($\\sum$ weights)", rep(sprintf("%.1f", ps$eff_obs), length(eqs))),
    c("Residual s.d. (posterior median)", sprintf("%.3f", sqrt(diag(ps$sigma))[eqs])),
    c("Saved posterior draws", rep(as.character(gssettings$ndraws), length(eqs)))
  )
  tab <- as.data.frame(do.call(rbind, c(rows, foot)), stringsAsFactors = FALSE)
  colnames(tab) <- c("Regressor", tex_escape(eqs))
  kbl(tab, format = "latex", booktabs = TRUE, escape = FALSE, longtable = TRUE,
      align = c("l", rep("c", length(eqs))), linesep = "",
      caption = sprintf("%s: posterior reduced-form coefficients (%s)", tex_escape(title), regime),
      label = label) %>%
    kable_styling(latex_options = "repeat_header", font_size = 7,
                  repeat_header_continued = "\\textit{(continued on next page)}") %>%
    add_header_above(c(" " = 1, "Equation (dependent variable)" = length(eqs))) %>%
    row_spec(2 * nrow(cf$median), hline_after = TRUE) %>%
    footnote(general = c(
      "Posterior medians with posterior standard deviations in parentheses.",
      "$^{**}$ / $^{*}$: the 90\\\\% / 68\\\\% posterior credible interval excludes zero.",
      "The surprise equations (OIS 3M, STOXX50) are block-exogenous:",
      "their lag coefficients and constants are restricted to zero.",
      "Minnesota prior with tightness 0.2 and lag decay 1."),
      escape = FALSE, general_title = "Notes:")
}

# Structural impact responses (median and 68% band) for each regime
impact_table_tex <- function(r, title, label) {
  vars <- rownames(r$posterior[[1]]$impact$median)
  cols <- list(); hdr <- c(" " = 1)
  for (reg in names(r$posterior)) {
    im <- r$posterior[[reg]]$impact
    for (k in 1:2) {
      cols[[paste(reg, k)]] <- sprintf("%.3f [%.3f, %.3f]", im$median[, k], im$q16[, k], im$q84[, k])
    }
    hdr <- c(hdr, setNames(2, reg))
  }
  tab <- data.frame(Variable = tex_escape(vars), cols, check.names = FALSE)
  colnames(tab) <- c("Variable", rep(r$shock_names, length(r$posterior)))
  if (length(r$posterior) == 1) hdr <- NULL
  k <- kbl(tab, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
           align = c("l", rep("c", ncol(tab) - 1)),
           caption = sprintf("%s: structural impact responses", tex_escape(title)),
           label = label) %>%
    kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 8)
  if (!is.null(hdr)) k <- add_header_above(k, hdr)
  footnote(k, general = paste(
    "Posterior median impact response ($h=0$) to a one-standard-deviation shock,",
    "68\\\\% credible band in brackets."),
    escape = FALSE, general_title = "Notes:", threeparttable = TRUE)
}

for (j in jobs) {
  r <- var_results[[j$file]]
  title <- switch(j$type,
                  linear = "Linear JK (2020) BVAR",
                  st     = sprintf("ST-BVAR (%s)", j$desc),
                  chol   = "ST-BVAR with standard HFI (Cholesky)")
  parts <- c(sprintf("%% Regression tables for %s (%s)", j$file, title), "")
  for (reg in names(r$posterior)) {
    parts <- c(parts,
               coef_table_tex(r$posterior[[reg]], reg, title,
                              paste0("var_coef_", tolower(gsub("[^A-Za-z0-9]", "", paste(j$file, reg)))),
                              Nm = length(mnames)),
               "")
  }
  parts <- c(parts, impact_table_tex(r, title,
                                     paste0("var_impact_", tolower(gsub("[^A-Za-z0-9]", "", j$file)))))
  save_tex(parts, reg_dir, paste0(j$file, "_Regression"))

  # Machine-readable version of the coefficient summaries
  coef_long <- bind_rows(lapply(names(r$posterior), function(reg) {
    cf <- r$posterior[[reg]]$coef
    bind_rows(lapply(names(cf), function(stat) {
      as.data.frame(as.table(cf[[stat]]), stringsAsFactors = FALSE) %>%
        setNames(c("Regressor", "Equation", "Value")) %>% mutate(Statistic = stat)
    })) %>% mutate(Regime = reg)
  })) %>%
    pivot_wider(names_from = Statistic, values_from = Value)
  write_csv(coef_long, file.path(reg_dir, paste0(j$file, "_coefficients.csv")))
}
