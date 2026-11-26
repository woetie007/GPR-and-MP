run_st_lp <- function(data, target_var, shock_var, gpr_var, p_lags = 4, horizon = 24, c_speed = 10, d_delay = 1) {
  
  # Ensure data is treated as a data frame and subset using base R (works for both matrices and data frames)
  cols_needed <- c(target_var, shock_var, gpr_var)
  df <- as.data.frame(data)[, cols_needed]
  df <- na.omit(df)
  
  T_total <- nrow(df)
  
  # 1. Calculate Transition Weights F(z_{t-d})
  z_t <- df[[gpr_var]][1:(T_total - d_delay)]
  gamma_val <- median(z_t, na.rm = TRUE)
  sigma_z <- sd(z_t, na.rm = TRUE)
  
  # Align weights with the usable sample
  F_high_full <- 1 / (1 + exp(-c_speed * (z_t - gamma_val) / sigma_z))
  F_low_full  <- 1 - F_high_full
  
  # Trim initial periods for lags
  start_idx <- p_lags + 1 + d_delay
  end_idx <- T_total - horizon
  
  results_list <- list()
  
  # 2. Horizon-by-Horizon Estimation Loop
  for (h in 0:horizon) {
    
    # Dependent variable: cumulative change (or level if specified)
    if (h == 0) {
      y_dep <- df[[target_var]][start_idx:end_idx]
    } else {
      y_dep <- df[[target_var]][(start_idx + h):(end_idx + h)] - df[[target_var]][start_idx:end_idx]
    }
    
    # Shock at time t
    shock_t <- df[[shock_var]][start_idx:end_idx]
    
    # Weights at time t-1 (or t-d)
    F_high <- F_high_full[(start_idx - d_delay):(end_idx - d_delay)]
    F_low  <- F_low_full[(start_idx - d_delay):(end_idx - d_delay)]
    
    # Build controls (lags of target variable only - omitting shock lags for degrees of freedom)
    controls_high <- matrix(0, nrow = length(y_dep), ncol = 0)
    controls_low  <- matrix(0, nrow = length(y_dep), ncol = 0)
    
    for (l in 1:p_lags) {
      lag_target <- lag(df[[target_var]], l)[start_idx:end_idx]
      
      controls_high <- cbind(controls_high, lag_target * F_high)
      controls_low  <- cbind(controls_low,  lag_target * F_low)
    }
    
    # Interaction Regressors
    shock_high <- shock_t * F_high
    shock_low  <- shock_t * F_low
    
    # Regressors Matrix
    X_matrix <- cbind(F_high, F_low, shock_high, shock_low, controls_high, controls_low)
    
    # OLS Regression for Horizon h
    fit <- lm(y_dep ~ 0 + X_matrix)
    
    # Extract coefficients and Newey-West standard errors (correcting for MA(h) serial correlation)
    nw_vcov <- NeweyWest(fit, lag = h, prewhite = FALSE)
    nw_se   <- sqrt(diag(nw_vcov))
    
    coefs <- coef(fit)
    
    # Safely extract coefficients
    beta_high <- if("X_matrixshock_high" %in% names(coefs)) coefs["X_matrixshock_high"] else NA
    se_high   <- if("X_matrixshock_high" %in% names(coefs)) nw_se["X_matrixshock_high"] else NA
    
    beta_low  <- if("X_matrixshock_low" %in% names(coefs)) coefs["X_matrixshock_low"] else NA
    se_low    <- if("X_matrixshock_low" %in% names(coefs)) nw_se["X_matrixshock_low"] else NA
    
    results_list[[h + 1]] <- data.frame(
      Horizon = h,
      Variable = target_var,
      Beta_Low = beta_low,
      Lower68_Low = beta_low - 1 * se_low,
      Upper68_Low = beta_low + 1 * se_low,
      Beta_High = beta_high,
      Lower68_High = beta_high - 1 * se_high,
      Upper68_High = beta_high + 1 * se_high
    )
  }
  
  return(do.call(rbind, results_list))
}