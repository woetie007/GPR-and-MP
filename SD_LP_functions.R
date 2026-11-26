# ==============================================================================
# STATE-DEPENDENT LOCAL PROJECTIONS (SD-LP) WITH JK (2020) IDENTIFICATION
# ==============================================================================

#' Estimate State-Dependent Local Projections across GPR Regimes
#'
#' @param Y Matrix/Dataframe of endogenous variables (T x N).
#'          ORDER MUST BE: [1] MP_shocks, [2] stocks_shocks, [3] i_ECB, [4] IP, [5] pi, [6] Oil
#' @param GPR Vector of length T: Exogenous Geopolitical Risk index
#' @param p Integer: Number of control lags in the LP regressions (default: 4)
#' @param d Integer: Delay lag for exogenous GPR transition variable (default: 1)
#' @param gamma Numeric: Threshold location parameter for GPR
#' @param c Numeric: Smoothness parameter (c > 0; higher = sharper transition)
#' @param horizon Integer: LP forecast horizon H (default: 24 months)
#' @param shock_type Character: "Pure MP" or "CBI"
#' @return A data frame containing median IRFs and Newey-West standard errors
estimate_sd_lp_jk2020 <- function(Y, GPR, p = 4, d = 1,
                                  gamma = median(GPR, na.rm = TRUE), c = 3.5,
                                  horizon = 24, shock_type = "Pure MP") {
  
  Y_mat   <- as.matrix(Y)
  T_total <- nrow(Y_mat)
  N       <- ncol(Y_mat)
  var_names <- colnames(Y_mat)
  if (is.null(var_names)) {
    var_names <- c("MP_shocks", "stocks_shocks", "i_ECB", "IP", "pi", "Oil")
  }
  
  # ----------------------------------------------------------------------------
  # Extract JK (2020) Structural Shock via Contemporaneous Rotation
  # ----------------------------------------------------------------------------
  # Fit a reduced-form VAR(1) or covariance on the surprise variables to 
  # identify the orthogonalized structural shock series E_t
  
  Sigma_u <- cov(Y_mat[, 1:2]) # Covariance of [MP_shocks, stocks_shocks]
  P_chol  <- t(chol(Sigma_u))
  
  # Draw orthogonal Haar rotation matrices Q until JK (2020) signs hold
  valid_Q <- FALSE
  while (!valid_Q) {
    X_rand    <- matrix(rnorm(4), nrow = 2, ncol = 2)
    qr_decomp <- qr(X_rand)
    Q         <- qr.Q(qr_decomp) %*% diag(sign(diag(qr.R(qr_decomp))))
    A0_sub    <- P_chol %*% Q
    
    # Check JK (2020) signs:
    # Col 1 (Pure MP): A0[1,1] > 0, A0[2,1] < 0
    # Col 2 (CBI):     A0[1,2] > 0, A0[2,2] > 0
    is_pure_mp <- (A0_sub[1, 1] > 0) && (A0_sub[2, 1] < 0)
    is_cbi     <- (A0_sub[1, 2] > 0) && (A0_sub[2, 2] > 0)
    
    if (!is_pure_mp && (A0_sub[1, 1] < 0) && (A0_sub[2, 1] > 0)) {
      Q[, 1] <- -Q[, 1]; is_pure_mp <- TRUE
    }
    if (!is_cbi && (A0_sub[1, 2] < 0) && (A0_sub[2, 2] < 0)) {
      Q[, 2] <- -Q[, 2]; is_cbi <- TRUE
    }
    if (is_pure_mp && is_cbi) valid_Q <- TRUE
  }
  
  # Invert to compute structural shock series: epsilon_t = Q' P^{-1} u_t
  u_surprises <- t(Y_mat[, 1:2])
  struct_shocks <- t(t(Q) %*% solve(P_chol) %*% u_surprises)
  
  # Select the target structural shock series (standardized to unit variance)
  shock_idx <- ifelse(shock_type == "Pure MP", 1, 2)
  shock_vec <- as.numeric(struct_shocks[, shock_idx])
  shock_vec <- shock_vec / sd(shock_vec, na.rm = TRUE)
  
  # ----------------------------------------------------------------------------
  # Construct Exogenous Logistic Transition Probabilities F(z_t)
  # ----------------------------------------------------------------------------
  z_t <- lag(GPR, d)
  sigma_z <- sd(z_t, na.rm = TRUE)
  F_high  <- 1 / (1 + exp(-c * (z_t - gamma) / sigma_z))
  F_low   <- 1 - F_high
  
  # ----------------------------------------------------------------------------
  # Run Local Projections for each Horizon h and Variable n
  # ----------------------------------------------------------------------------
  results_list <- list()
  
  for (n_idx in 1:N) { # Loop over macro variables: i_ECB, IP, pi, Oil
    target_name <- var_names[n_idx]
    
    for (h in 0:horizon) {
      
      # Dependent variable: Forward shift Y_{i, t+h}
      # (For h=0, it is the contemporaneous response Y_{i, t})
      y_lh <- dplyr::lead(Y_mat[, n_idx], h)
      
      # Base estimation dataframe
      df_reg <- data.frame(
        y        = y_lh,
        shock    = shock_vec,
        F_high   = F_high,
        F_low    = F_low
      )
      
      # Add p lagged controls of ALL endogenous variables + lagged GPR
      for (l in 1:p) {
        for (v in 1:N) {
          df_reg[[sprintf("lag_%s_%d", var_names[v], l)]] <- dplyr::lag(Y_mat[, v], l)
        }
        df_reg[[sprintf("lag_GPR_%d", l)]] <- dplyr::lag(GPR, l)
      }
      
      # Remove rows with NAs created by leading/lagging
      df_clean <- na.omit(df_reg)
      
      # State-Dependent LP Regression Equation:
      # y_{t+h} = F_low * alpha_L + F_high * alpha_H 
      #           + F_low * beta_L * shock_t + F_high * beta_H * shock_t
      #           + Gamma * Control_Lags + e_{t+h}
      
      control_cols <- grep("^lag_", colnames(df_clean), value = TRUE)
      formula_str  <- sprintf("y ~ 0 + F_low + F_high + I(F_low * shock) + I(F_high * shock) + %s",
                              paste(control_cols, collapse = " + "))
      
      fit <- lm(as.formula(formula_str), data = df_clean)
      
      # Newey-West HAC Standard Errors (lag truncation = h + 1)
      nw_cov <- sandwich::NeweyWest(fit, lag = h + 1, prewhite = FALSE)
      coefs  <- coef(fit)
      se_nw  <- sqrt(diag(nw_cov))
      
      # Extract Low-GPR coefficient and HAC standard error
      beta_low <- coefs["I(F_low * shock)"]
      se_low   <- se_nw["I(F_low * shock)"]
      
      # Extract High-GPR coefficient and HAC standard error
      beta_high <- coefs["I(F_high * shock)"]
      se_high   <- se_nw["I(F_high * shock)"]
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Horizon  = h,
        Variable = target_name,
        Shock    = shock_type,
        Regime   = "Low GPR Regime",
        Beta     = as.numeric(beta_low),
        SE       = as.numeric(se_low),
        Lower68  = as.numeric(beta_low - 1.0 * se_low),
        Upper68  = as.numeric(beta_low + 1.0 * se_low),
        Lower90  = as.numeric(beta_low - 1.645 * se_low),
        Upper90  = as.numeric(beta_low + 1.645 * se_low)
      )
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Horizon  = h,
        Variable = target_name,
        Shock    = shock_type,
        Regime   = "High GPR Regime",
        Beta     = as.numeric(beta_high),
        SE       = as.numeric(se_high),
        Lower68  = as.numeric(beta_high - 1.0 * se_high),
        Upper68  = as.numeric(beta_high + 1.0 * se_high),
        Lower90  = as.numeric(beta_high - 1.645 * se_high),
        Upper90  = as.numeric(beta_high + 1.645 * se_high)
      )
    }
  }
  
  return(do.call(rbind, results_list))
}


