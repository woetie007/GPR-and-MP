# ==============================================================================
# HELPER: NORMALIZE IRF ARRAY TO A SPECIFIC IMPACT TARGET
# ==============================================================================

normalize_irf_array <- function(irf_array, target_var_idx, target_size, shock_idx = 1) {
  # irf_array dimensions: [n_draws, horizon + 1, n_vars, n_shocks]
  n_draws <- dim(irf_array)[1]
  
  for (s in 1:n_draws) {
    # Get the contemporaneous impact (Month 0, which is index 1) 
    # of the target variable to the selected shock
    impact_val <- irf_array[s, 1, target_var_idx, shock_idx]
    
    # Scale the entire response path across all variables and horizons
    # (The >= 0.02 floor from the BVAR ensures we never divide by zero)
    if (abs(impact_val) > 1e-8) {
      scale_factor <- target_size / impact_val
      irf_array[s, , , shock_idx] <- irf_array[s, , , shock_idx] * scale_factor
    }
  }
  
  return(irf_array)
}

# ==============================================================================
# UNIFIED HYPERPARAMETER OPTIMIZATION FOR SMOOTH TRANSITION BVAR
# Jointly selects (p, d, gamma, c, lambda1, lambda3) via Marginal Data Density
# ==============================================================================

library(dplyr)

#' Jointly Optimize All ST-BVAR Hyperparameters via Pooled Log-MDD
#'
#' @param Y Matrix of endogenous variables (T_total x N)
#' @param GPR Exogenous Geopolitical Risk vector of length T_total
#' @param p_grid Vector of candidate VAR lag lengths (e.g., c(4, 5, 6))
#' @param d_grid Vector of candidate GPR delay lags (e.g., c(1, 2))
#' @param gamma_pcts Vector of candidate GPR percentile thresholds (e.g., c(0.50, 0.60, 0.70))
#' @param c_grid Vector of candidate transition speeds (e.g., c(3.0, 5.0, 10.0))
#' @param lambda1_grid Vector of candidate overall tightness values
#' @param lambda3_grid Vector of candidate lag decay rates
#' @param prior_means Vector of length N for own 1st-lag prior means
#' @return A list containing the top ranked models table and the winning hyperparameters
optimize_all_st_bvar_hyperparameters <- function(
    Y, GPR,
    p_grid       = c(4, 6),
    d_grid       = c(1, 2),
    gamma_pcts   = c(0.50, 0.65, 0.75),
    c_grid       = c(3.5, 5.0, 10.0),
    lambda1_grid = c(0.08, 0.10, 0.12, 0.15, 0.20),
    lambda3_grid = c(1.5, 2.0, 3.0),
    prior_means  = c(0, 0, 0.9, 0.9, 0.9, 0.9)
) {
  
  T_total <- nrow(Y)
  N       <- ncol(Y)
  
  # CRITICAL: Lock common evaluation window to max(p_grid) so MDDs are comparable across lags
  p_max   <- max(p_grid)
  T_eval  <- T_total - p_max
  
  # Helper: Multivariate log-Gamma function ln(Gamma_N(x))
  lmvgamma <- function(x, dim_n) {
    0.25 * dim_n * (dim_n - 1) * log(pi) + sum(lgamma(x + (1 - 1:dim_n) / 2))
  }
  
  # Helper: Residual variance s_sq from AR(1) for Minnesota scaling
  s_sq <- numeric(N)
  for (i in 1:N) {
    y_i <- Y[-1, i]
    x_i <- cbind(1, Y[-T_total, i])
    valid <- complete.cases(x_i, y_i)
    b_i <- solve(t(x_i[valid, ]) %*% x_i[valid, ]) %*% t(x_i[valid, ]) %*% y_i[valid]
    e_i <- y_i[valid] - x_i[valid, ] %*% b_i
    s_sq[i] <- sum(e_i^2) / (sum(valid) - 2)
  }
  
  # Helper: Compute Closed-Form Analytical Log-MDD for a single weighted regime
  eval_regime_mdd <- function(Y_mat, p_val, weights, lam1, lam3) {
    T_eff  <- sum(weights)
    Y_sub  <- Y_mat[(p_max + 1):T_total, , drop = FALSE]
    
    # Build lag regressor matrix X_sub for current lag order p_val
    X_sub <- matrix(1, nrow = nrow(Y_sub), ncol = 1)
    for (l in 1:p_val) {
      l_start <- (p_max - l + 1)
      l_end   <- (T_total - l)
      X_sub   <- cbind(X_sub, Y_mat[l_start:l_end, , drop = FALSE])
    }
    
    W_sqrt <- diag(sqrt(weights))
    Y_w    <- W_sqrt %*% Y_sub
    X_w    <- W_sqrt %*% X_sub
    
    K   <- 1 + N * p_val
    nu0 <- N + 2
    S0  <- diag(s_sq) * (nu0 - N - 1)
    
    B0 <- matrix(0, nrow = K, ncol = N)
    for (i in 1:N) { B0[1 + i, i] <- prior_means[i] }
    
    v0_diag <- numeric(K)
    v0_diag[1] <- 1e4
    idx <- 2
    for (l in 1:p_val) {
      for (j in 1:N) {
        v0_diag[idx] <- (lam1 / ((l^lam3) * sqrt(s_sq[j])))^2
        idx <- idx + 1
      }
    }
    V0     <- diag(v0_diag)
    V0_inv <- diag(1 / v0_diag)
    
    XtX_w  <- t(X_w) %*% X_w
    Vn     <- solve(V0_inv + XtX_w)
    Vn_inv <- V0_inv + XtX_w
    
    Bn   <- Vn %*% (V0_inv %*% B0 + t(X_w) %*% Y_w)
    Sn   <- S0 + t(Y_w) %*% Y_w + t(B0) %*% V0_inv %*% B0 - t(Bn) %*% Vn_inv %*% Bn
    nu_n <- max(round(T_eff) + nu0 - K, N + 2)
    
    log_det_Vn <- as.numeric(determinant(Vn, logarithm = TRUE)$modulus)
    log_det_V0 <- as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    log_det_S0 <- as.numeric(determinant(S0, logarithm = TRUE)$modulus)
    log_det_Sn <- as.numeric(determinant(Sn, logarithm = TRUE)$modulus)
    
    log_mdd <- -0.5 * T_eff * N * log(pi) +
      0.5 * N * (log_det_Vn - log_det_V0) +
      0.5 * nu0 * log_det_S0 -
      0.5 * nu_n * log_det_Sn +
      lmvgamma(nu_n / 2, N) - lmvgamma(nu0 / 2, N)
    return(log_mdd)
  }
  
  # Calculate total grid combinations for progress tracking
  total_models <- length(p_grid) * length(d_grid) * length(gamma_pcts) * 
    length(c_grid) * length(lambda1_grid) * length(lambda3_grid)
  cat(sprintf("Starting joint optimization across %d candidate ST-BVAR specifications...\n", total_models))
  
  results_list <- list()
  model_id     <- 0
  
  for (p_val in p_grid) {
    for (d_val in d_grid) {
      for (g_pct in gamma_pcts) {
        g_val <- as.numeric(quantile(GPR, probs = g_pct, na.rm = TRUE))
        
        for (c_val in c_grid) {
          # Compute logistic transition weights over the common evaluation window
          z_t     <- GPR[(p_max + 1 - d_val):(T_total - d_val)]
          sigma_z <- sd(z_t, na.rm = TRUE)
          F_high  <- 1 / (1 + exp(-c_val * (z_t - g_val) / sigma_z))
          F_low   <- 1 - F_high
          
          for (lam1 in lambda1_grid) {
            for (lam3 in lambda3_grid) {
              model_id <- model_id + 1
              
              mdd_low  <- eval_regime_mdd(Y, p_val, F_low,  lam1, lam3)
              mdd_high <- eval_regime_mdd(Y, p_val, F_high, lam1, lam3)
              pooled   <- mdd_low + mdd_high
              
              results_list[[length(results_list) + 1]] <- data.frame(
                Model_ID       = model_id,
                p              = p_val,
                d              = d_val,
                Gamma_Pct      = g_pct,
                Gamma_Val      = g_val,
                c              = c_val,
                Lambda1        = lam1,
                Lambda3        = lam3,
                Pooled_Log_MDD = pooled
              )
            }
          }
        }
      }
    }
  }
  
  df_res <- do.call(rbind, results_list)
  
  # Calculate Posterior Model Probabilities (PMP)
  max_mdd <- max(df_res$Pooled_Log_MDD, na.rm = TRUE)
  df_res$PMP <- exp(df_res$Pooled_Log_MDD - max_mdd) / sum(exp(df_res$Pooled_Log_MDD - max_mdd))
  
  df_sorted <- df_res %>% arrange(desc(Pooled_Log_MDD))
  best_mod  <- df_sorted[1, ]
  
  cat("\n==============================================================================\n")
  cat("OPTIMAL JOINT ST-BVAR SPECIFICATION SELECTED BY MARGINAL DATA DENSITY\n")
  cat("==============================================================================\n")
  cat(sprintf(" -> Lags (p*):                  %d months\n", best_mod$p))
  cat(sprintf(" -> Delay (d*):                 %d month(s)\n", best_mod$d))
  cat(sprintf(" -> Threshold (gamma*):         %.2f (%.0fth percentile)\n", best_mod$Gamma_Val, best_mod$Gamma_Pct * 100))
  cat(sprintf(" -> Transition Speed (c*):      %.1f\n", best_mod$c))
  cat(sprintf(" -> Overall Tightness (lam1*):  %.2f\n", best_mod$Lambda1))
  cat(sprintf(" -> Lag Decay Rate (lam3*):     %.1f\n", best_mod$Lambda3))
  cat(sprintf(" -> Posterior Model Prob (PMP): %.1f%%\n", best_mod$PMP * 100))
  cat("==============================================================================\n")
  
  return(list(
    all_results  = df_sorted,
    top_10       = head(df_sorted, 10),
    optimal      = best_mod
  ))
}

# ==============================================================================
# EMPIRICAL CHECK FOR PRIOR MEANS (AR(1) PERSISTENCE)
# ==============================================================================

check_empirical_ar1 <- function(Y_matrix) {
  var_names <- colnames(Y_matrix)
  ar1_values <- numeric(ncol(Y_matrix))
  
  for (i in 1:ncol(Y_matrix)) {
    y_t <- Y_matrix[-1, i]
    y_l1 <- Y_matrix[-nrow(Y_matrix), i]
    # Simple AR(1) slope: corr(y_t, y_t-1)
    ar1_values[i] <- cor(y_t, y_l1, use = "complete.obs")
  }
  
  df_ar1 <- data.frame(
    Variable       = var_names,
    Empirical_AR1  = round(ar1_values, 3),
    Recommended_B0 = ifelse(abs(ar1_values) < 0.35, 0.0, 0.9)
  )
  return(df_ar1)
}

# ==============================================================================
# CHECK VAR STABILITY (REJECT EXPLOSIVE DRAW IF MAX EIGENVALUE >= 1.0)
# ==============================================================================

is_stable_var <- function(B, p, N) {
  # B is (K x N) where K = 1 + N*p. Remove intercept:
  A_lags <- t(B[-1, , drop = FALSE])
  
  # Construct the (N*p) x (N*p) Companion Matrix
  companion <- matrix(0, nrow = N * p, ncol = N * p)
  companion[1:N, ] <- A_lags
  
  if (p > 1) {
    companion[(N + 1):(N * p), 1:(N * (p - 1))] <- diag(N * (p - 1))
  }
  
  # Compute eigenvalues and check if max modulus is strictly inside unit circle
  max_eig <- max(abs(eigen(companion, only.values = TRUE)$values))
  return(max_eig < 0.998) # Use 0.998 to allow high persistence without explosion
}

# ==============================================================================
# SMOOTH TRANSITION BAYESIAN VAR WITH JK (2020) SIGN RESTRICTIONS
# ==============================================================================

library(MASS)

#' Estimate Smooth Transition JK (2020) Sign-Restricted BVAR
#'
#' @param Y Matrix of endogenous variables (T x N). 
#' @param GPR Exogenous Geopolitical Risk vector of length T
#' @param p Integer: Number of VAR lags
#' @param d Integer: Lag delay for GPR transition variable (default d = 1)
#' @param gamma Numeric: Threshold location parameter for GPR (e.g., sample mean or median)
#' @param c Numeric: Smoothness parameter (c > 0; higher = sharper transition)
#' @param n_draws Integer: Number of successful posterior structural draws to keep
#' @param horizon Integer: IRF forecast horizon (months)
#' @return A list containing arrays of IRFs for Low-GPR and High-GPR regimes
estimate_st_jk2020_bvar <- function(Y, GPR, p = 6, d = 1, 
                                    gamma = median(GPR), c = 10.0, 
                                    lambda1     = 0.15,
                                    lambda3     = 2.0,
                                    prior_means = c(0, 0, 0.9, 0.9, 0.9, 0.9), 
                                    n_draws = 1000, horizon = 24) {
  
  T_total <- nrow(Y)
  N <- ncol(Y)
  
  # ----------------------------------------------------------------------------
  # Helper: Construct lagged regressor matrix
  # ----------------------------------------------------------------------------
  build_var_matrices <- function(data, lags) {
    T_obs <- nrow(data)
    Y_trim <- data[(lags + 1):T_obs, , drop = FALSE]
    X_mat <- matrix(1, nrow = T_obs - lags, ncol = 1) # Intercept column
    
    for (l in 1:lags) {
      X_mat <- cbind(X_mat, data[(lags - l + 1):(T_obs - l), , drop = FALSE])
    }
    return(list(Y = Y_trim, X = X_mat))
  }
  
  var_data <- build_var_matrices(Y, p)
  Y_t <- var_data$Y
  X_t <- var_data$X
  T_eff <- nrow(Y_t)
  
  # ----------------------------------------------------------------------------
  # Logistic transition function F(z_t) IN [0, 1]
  # ----------------------------------------------------------------------------
  # Align lagged exogenous GPR transition variable z_t = GPR_{t-d}
  z_t <- GPR[(p + 1 - d):(T_total - d)]
  
  # Standardize z_t by sample standard deviation for numerical scale invariance
  sigma_z <- sd(z_t)
  F_high <- 1 / (1 + exp(-c * (z_t - gamma) / sigma_z))
  F_low  <- 1 - F_high
  
  # ----------------------------------------------------------------------------
  # Helper: Weighted Normal-Inverse-Wishart (WLS) Posterior Draw
  # ----------------------------------------------------------------------------
  draw_weighted_niw_posterior <- function(Y_sub, X_sub, weights, 
                                          lambda1 = 0.15, lambda3 = 2.0,
                                          prior_means = c(0, 0, 0.9, 0.9, 0.9, 0.9)) {
    T_sub <- nrow(Y_sub)
    K     <- ncol(X_sub)
    eff_T <- sum(weights)
    
    # Calculate AR(1) residual variances (s_sq) for Minnesota scale standardization
    s_sq <- numeric(N)
    for (i in 1:N) {
      y_i <- Y_sub[, i]
      x_i <- cbind(1, dplyr::lag(y_i, 1))
      valid <- complete.cases(x_i, y_i)
      b_i <- solve(t(x_i[valid, ]) %*% x_i[valid, ]) %*% t(x_i[valid, ]) %*% y_i[valid]
      e_i <- y_i[valid] - x_i[valid, ] %*% b_i
      s_sq[i] <- sum(e_i^2) / (sum(valid) - 2)
    }
    
    # Build Minnesota Prior Mean Matrix B0 (K x N)
    B0 <- matrix(0, nrow = K, ncol = N)
    for (i in 1:N) {
      B0[1 + i, i] <- prior_means[i] # First own-lag prior expectation
    }
    
    # Build Minnesota Prior Covariance V0 (K x K diagonal matrix)
    v0_diag <- numeric(K)
    v0_diag[1] <- 1e4 # Uninformative intercept
    
    idx <- 2
    for (l in 1:p) {
      for (j in 1:N) {
        # Standard Minnesota decay formula using lambda1 and lambda3
        v0_diag[idx] <- (lambda1 / ((l^lambda3) * sqrt(s_sq[j])))^2
        idx <- idx + 1
      }
    }
    V0_inv <- diag(1 / v0_diag)
    
    # Apply smooth WLS weights to data
    W_sqrt <- diag(sqrt(weights))
    Y_w    <- W_sqrt %*% Y_sub
    X_w    <- W_sqrt %*% X_sub
    XtX_w  <- t(X_w) %*% X_w
    
    # Compute Conjugate Posterior Hyperparameters (Vn, Bn, Sn, nu_n)
    Vn     <- solve(V0_inv + XtX_w)
    Vn_inv <- V0_inv + XtX_w
    
    Bn   <- Vn %*% (V0_inv %*% B0 + t(X_w) %*% Y_w)
    
    # Degrees of freedom and prior scale S0
    nu0  <- N + 2
    S0   <- diag(s_sq) * (nu0 - N - 1)
    Sn   <- S0 + crossprod(Y_w) + t(B0) %*% V0_inv %*% B0 - t(Bn) %*% Vn_inv %*% Bn
    nu_n <- nu_n <- max(round(eff_T) + nu0, N + 2)
    
    # Draw Sigma ~ Inverse-Wishart(Sn, nu_n)
    C_chol <- chol(Sn)
    Z_w    <- matrix(rnorm(nu_n * N), nrow = nu_n, ncol = N)
    Sigma_draw <- t(C_chol) %*% solve(t(Z_w) %*% Z_w) %*% C_chol
    
    # Draw B ~ Matrix-Normal(Bn, Sigma_draw, Vn)
    Z_b    <- matrix(rnorm(K * N), nrow = K, ncol = N)
    B_draw <- Bn + t(chol(Vn)) %*% Z_b %*% chol(Sigma_draw)
    
    return(list(B = B_draw, Sigma = Sigma_draw))
  }
  
  # ----------------------------------------------------------------------------
  # Helper: Compute IRFs
  # ----------------------------------------------------------------------------
  compute_irf <- function(B, A0, lags, horiz) {
    A_lags <- t(B[-1, , drop = FALSE]) # Remove intercept, transpose to companion shape
    irf_array <- array(0, dim = c(horiz + 1, N, N))
    irf_array[1, , ] <- A0
    
    if (horiz >= 1) {
      for (h in 1:horiz) {
        response_h <- matrix(0, nrow = N, ncol = N)
        for (l in 1:min(h, lags)) {
          col_start <- (l - 1) * N + 1
          col_end   <- l * N
          A_l <- A_lags[, col_start:col_end]
          response_h <- response_h + A_l %*% irf_array[h - l + 1, , ]
        }
        irf_array[h + 1, , ] <- response_h
      }
    }
    return(irf_array)
  }
  
  # ----------------------------------------------------------------------------
  # Main MCMC + Smooth Transition + JK (2020) Identification Loop
  # ----------------------------------------------------------------------------
  irf_draws_low  <- array(0, dim = c(n_draws, horizon + 1, N, N))
  irf_draws_high <- array(0, dim = c(n_draws, horizon + 1, N, N))
  
  cat("Starting Smooth Transition MCMC + JK (2020) Identification...\n")
  
  for (regime_name in c("Low", "High")) {
    
    draw_count <- 0
    regime_weights <- if (regime_name == "Low") F_low else F_high
    total_attempts <- 0
    
    while (draw_count < n_draws) {
      
      total_attempts <- total_attempts + 1
      
      # Print progress every 100 attempts:
      if (total_attempts %% 100 == 0) {
        cat(sprintf("-> Attempts: %d | Accepted: %d (Hit Rate: %.1f%%)\n", 
                    total_attempts, draw_count, (draw_count / total_attempts) * 100))
      }
      # Keep drawing NIW parameters until we get a non-explosive VAR system
      is_stable <- FALSE
      while (!is_stable) {
        niw <- draw_weighted_niw_posterior(
          Y_t, X_t, regime_weights, 
          lambda1     = lambda1, 
          lambda3     = lambda3,
          prior_means = prior_means
        )
        is_stable <- is_stable_var(niw$B, p, N)
      }
      P_chol <- t(chol(niw$Sigma)) # Lower-triangular Cholesky factor
      
      # Draw random Haar orthonormal matrix Q via QR decomposition
      X_rand <- matrix(rnorm(N * N), nrow = N, ncol = N)
      qr_decomp <- qr(X_rand)
      Q <- qr.Q(qr_decomp)
      R_qr <- qr.R(qr_decomp)
      Q <- Q %*% diag(sign(diag(R_qr))) # Positive diagonal normalization
      
      # Candidate structural contemporaneous impact matrix
      A0 <- P_chol %*% Q
      
      # Apply JK (2020) Sign Restrictions on Top-Left 2x2 Block
      # Col 1: Pure MP Shock -> (+) MP_shocks [1,1], (-) stocks_shocks [2,1]
      # Col 2: CBI Shock     -> (+) MP_shocks [1,2], (+) stocks_shocks [2,2]
      
      is_pure_mp <- (A0[1, 1] > 0) && (A0[2, 1] < 0) 
      is_cbi     <- (A0[1, 2] > 0) && (A0[2, 2] > 0) 
      
      # Correct column sign inversions if structural shocks are flipped by -1
      if (!is_pure_mp || !is_cbi) {
        if ((A0[1, 1] < 0) && (A0[2, 1] > 0)) {
          A0[, 1] <- -A0[, 1]
          is_pure_mp <- TRUE
        }
        if ((A0[1, 2] < 0) && (A0[2, 2] < 0)) {
          A0[, 2] <- -A0[, 2]
          is_cbi <- TRUE
        }
      }
      
      # Retain draw if JK (2020) restrictions hold simultaneously
      if (is_pure_mp && is_cbi) {
        draw_count <- draw_count + 1
        irf_mat <- compute_irf(niw$B, A0, p, horizon)
        
        if (regime_name == "Low") {
          irf_draws_low[draw_count, , , ] <- irf_mat
        } else {
          irf_draws_high[draw_count, , , ] <- irf_mat
        }
      }
    }
    cat(sprintf("-> Completed %d valid structural draws for %s GPR Regime.\n", 
                n_draws, regime_name))
  }
  
  return(list(
    IRF_Low_GPR      = irf_draws_low,
    IRF_High_GPR     = irf_draws_high,
    Transition_Weights = data.frame(F_low = F_low, F_high = F_high)
  ))
}


# ==============================================================================
# HELPER: EXTRACT IRF QUANTILES (68% CREDIBLE BANDS)
# ==============================================================================

extract_irf_quantiles <- function(irf_array) {
  # irf_array dimensions: [draws, horizon + 1, variables, shocks]
  # We extract Shock 1 (Pure Monetary Policy Shock)
  
  apply(irf_array[, , , 1], MARGIN = c(2, 3), function(x) {
    c(
      Lower68 = quantile(x, 0.16, na.rm = TRUE),
      Median  = median(x, na.rm = TRUE),
      Upper68 = quantile(x, 0.84, na.rm = TRUE)
    )
  })
}

# ==============================================================================
# HELPER: FLATTEN QUANTILES INTO A GGPLOT-FRIENDLY DATAFRAME
# ==============================================================================

build_plot_df <- function(quant_array, regime_name, var_names) {
  df_list <- list()
  
  for (i in seq_along(var_names)) {
    # quant_array dimensions are: [3 (Lower, Median, Upper), horizon + 1, N (Variables)]
    tmp_df <- data.frame(
      Horizon  = 0:(dim(quant_array)[2] - 1),
      Variable = var_names[i],
      Regime   = regime_name,
      Lower68  = quant_array[1, , i],
      Median   = quant_array[2, , i],
      Upper68  = quant_array[3, , i]
    )
    df_list[[i]] <- tmp_df
  }
  
  return(do.call(rbind, df_list))
}
# ==============================================================================
# SMOOTH TRANSITION BAYESIAN VAR: INDEPENDENT NORMAL-WISHART (GIBBS SAMPLER)
# Exact Litterman (1986) Framework with Cross-Variable Tightness
# ==============================================================================

library(MASS)

#' Estimate Smooth Transition JK (2020) BVAR via Gibbs Sampling
#'
#' @param Y Matrix of endogenous variables (T x N). 
#'          ORDER MUST BE: [1] OIS_3M, [2] STOXX50, [3] IP, [4] Oil, [5] pi, [6] i_ECB
#' @param GPR Exogenous Geopolitical Risk vector of length T
#' @param p Integer: Number of VAR lags 
#' @param d Integer: Lag delay for GPR transition variable
#' @param gamma Numeric: Threshold location parameter for GPR 
#' @param c Numeric: Smoothness parameter
#' @param lambda1 Numeric: Overall Tightness (e.g., 0.20)
#' @param lambda2 Numeric: Cross-Variable Tightness (e.g., 0.50)
#' @param lambda3 Numeric: Lag Decay (e.g., 2.0)
#' @param prior_means Vector of length N: 0 for surprises, 0.9 for levels
#' @param n_draws Integer: Number of valid structural draws to keep
#' @param burn_in Integer: Number of initial Gibbs draws to discard
#' @param horizon Integer: IRF forecast horizon (months)
estimate_st_jk2020_gibbs <- function(Y, GPR, p = 6, d = 1, 
                                     gamma = median(GPR), c = 10.0, 
                                     lambda1 = 0.20, lambda2 = 0.30, lambda3 = 2.0,
                                     prior_means = c(0, 0, 1, 1, 1, 1), 
                                     n_draws = 1000, burn_in = 500, horizon = 24) {
  
  T_total <- nrow(Y)
  N <- ncol(Y)
  
  # ----------------------------------------------------------------------------
  # 1. Construct Lagged Matrices
  # ----------------------------------------------------------------------------
  build_var_matrices <- function(data, lags) {
    T_obs <- nrow(data)
    Y_trim <- data[(lags + 1):T_obs, , drop = FALSE]
    X_mat <- matrix(1, nrow = T_obs - lags, ncol = 1) 
    for (l in 1:lags) {
      X_mat <- cbind(X_mat, data[(lags - l + 1):(T_obs - l), , drop = FALSE])
    }
    return(list(Y = Y_trim, X = X_mat))
  }
  
  var_data <- build_var_matrices(Y, p)
  Y_t <- var_data$Y
  X_t <- var_data$X
  K <- ncol(X_t)
  
  # Logistic transition F(z_t)
  z_t <- GPR[(p + 1 - d):(T_total - d)]
  F_high <- 1 / (1 + exp(-c * (z_t - gamma) / sd(z_t)))
  F_low  <- 1 - F_high
  
  # ----------------------------------------------------------------------------
  # 2. Build Independent Normal-Wishart Prior (Litterman 1986)
  # ----------------------------------------------------------------------------
  s_sq <- numeric(N)
  for (i in 1:N) {
    y_i <- Y_t[, i]
    x_i <- cbind(1, dplyr::lag(y_i, 1))
    valid <- complete.cases(x_i, y_i)
    b_i <- solve(t(x_i[valid, ]) %*% x_i[valid, ]) %*% t(x_i[valid, ]) %*% y_i[valid]
    s_sq[i] <- sum((y_i[valid] - x_i[valid, ] %*% b_i)^2) / (sum(valid) - 2)
  }
  
  # Prior Mean Vector for vectorized B (KN x 1)
  B0_mat <- matrix(0, nrow = K, ncol = N)
  for (i in 1:N) { B0_mat[1 + i, i] <- prior_means[i] }
  b0 <- as.vector(B0_mat)
  
  # Prior Covariance V0 for vectorized B (KN x KN Diagonal)
  v0_diag <- numeric(K * N)
  idx <- 1
  for (i in 1:N) { # Equation i
    v0_diag[idx] <- 1e4 * s_sq[i] # Intercept
    idx <- idx + 1
    for (l in 1:p) {
      for (j in 1:N) { # Variable j
        if (i == j) {
          # Own Lag Tightness
          v0_diag[idx] <- (lambda1 / (l^lambda3))^2
        } else {
          # Cross-Variable Tightness (lambda2)
          v0_diag[idx] <- ((lambda1 * lambda2) / (l^lambda3))^2 * (s_sq[i] / s_sq[j])
        }
        idx <- idx + 1
      }
    }
  }
  V0_inv <- diag(1 / v0_diag)
  
  # Prior for Sigma (Inverse-Wishart)
  nu0 <- N + 2
  S0  <- diag(s_sq) * (nu0 - N - 1)
  
  # Helper function to compute IRFs
  compute_irf <- function(B_mat, A0, lags, horiz) {
    A_lags <- t(B_mat[-1, , drop = FALSE]) 
    irf_array <- array(0, dim = c(horiz + 1, N, N))
    irf_array[1, , ] <- A0
    if (horiz >= 1) {
      for (h in 1:horiz) {
        resp <- matrix(0, nrow = N, ncol = N)
        for (l in 1:min(h, lags)) {
          A_l <- A_lags[, ((l - 1) * N + 1):(l * N)]
          resp <- resp + A_l %*% irf_array[h - l + 1, , ]
        }
        irf_array[h + 1, , ] <- resp
      }
    }
    return(irf_array)
  }
  
  # ----------------------------------------------------------------------------
  # 3. Main MCMC Loop (Gibbs Sampler -> Sign Restriction Filter)
  # ----------------------------------------------------------------------------
  irf_draws_low  <- array(0, dim = c(n_draws, horizon + 1, N, N))
  irf_draws_high <- array(0, dim = c(n_draws, horizon + 1, N, N))
  
  cat("Starting Independent Normal-Wishart Gibbs Sampler...\n")
  
  for (regime_name in c("Low", "High")) {
    
    weights <- if (regime_name == "Low") F_low else F_high
    W_sqrt  <- diag(sqrt(weights))
    Y_w     <- W_sqrt %*% Y_t
    X_w     <- W_sqrt %*% X_t
    XtX     <- crossprod(X_w)
    eff_T   <- sum(weights)
    nu_bar  <- nu0 + eff_T
    
    # Initialize Sigma with OLS estimate
    B_ols <- solve(XtX) %*% crossprod(X_w, Y_w)
    Sigma_current <- crossprod(Y_w - X_w %*% B_ols) / (eff_T - K)
    
    valid_draws <- 0
    total_steps <- 0
    
    while (valid_draws < n_draws) {
      total_steps <- total_steps + 1
      
      # A. Draw B conditional on Sigma
      Sigma_inv <- solve(Sigma_current)
      V_bar_inv <- V0_inv + kronecker(Sigma_inv, XtX)
      V_bar     <- chol2inv(chol(V_bar_inv))
      
      # Vectorized OLS logic for the mean
      term2 <- as.vector(crossprod(X_w, Y_w) %*% Sigma_inv)
      b_bar <- V_bar %*% (V0_inv %*% b0 + term2)
      
      # Ensure non-explosive draws (Truncated Normal)
      is_stable <- FALSE
      while (!is_stable) {
        B_vec <- b_bar + t(chol(V_bar)) %*% rnorm(K * N)
        B_mat <- matrix(B_vec, nrow = K, ncol = N)
        is_stable <- is_stable_var(B_mat, p, N)
      }
      
      # B. Draw Sigma conditional on B
      E_mat <- Y_w - X_w %*% B_mat
      S_bar <- S0 + crossprod(E_mat)
      Sigma_current <- solve(rWishart(1, df = nu_bar, Sigma = chol2inv(chol(S_bar)))[,,1])
      
      # C. After burn-in, apply Haar rotation and JK (2020) Sign Restrictions
      if (total_steps > burn_in) {
        P_chol <- t(chol(Sigma_current))
        
        X_rand <- matrix(rnorm(N * N), nrow = N, ncol = N)
        qr_decomp <- qr(X_rand)
        Q <- qr.Q(qr_decomp) %*% diag(sign(diag(qr.R(qr_decomp))))
        
        A0 <- P_chol %*% Q
        
        is_pure_mp <- (A0[1, 1] > 0) && (A0[2, 1] < 0)
        is_cbi     <- (A0[1, 2] > 0) && (A0[2, 2] > 0)
        
        if (!is_pure_mp || !is_cbi) {
          if ((A0[1, 1] < 0) && (A0[2, 1] > 0)) { A0[, 1] <- -A0[, 1]; is_pure_mp <- TRUE }
          if ((A0[1, 2] < 0) && (A0[2, 2] < 0)) { A0[, 2] <- -A0[, 2]; is_cbi <- TRUE }
        }
        
        if (is_pure_mp && is_cbi) {
          valid_draws <- valid_draws + 1
          irf_mat <- compute_irf(B_mat, A0, p, horizon)
          
          if (regime_name == "Low") { irf_draws_low[valid_draws, , , ] <- irf_mat } 
          else { irf_draws_high[valid_draws, , , ] <- irf_mat }
          
          if (valid_draws %% 100 == 0) {
            cat(sprintf("-> Valid Draws for %s: %d (Total Gibbs Steps: %d)\n", 
                        regime_name, valid_draws, total_steps))
          }
        }
      }
    }
  }
  
  return(list(
    IRF_Low_GPR      = irf_draws_low,
    IRF_High_GPR     = irf_draws_high,
    Transition_Weights = data.frame(F_low = F_low, F_high = F_high)
  ))
}


# ==============================================================================
# ORIGINAL JAROCIŃSKI & KARADI (2020) LINEAR BAYESIAN VAR
# Independent Normal-Wishart Prior via Gibbs Sampling
# ==============================================================================

library(MASS)

#' Estimate Linear JK (2020) Sign-Restricted BVAR
#'
#' @param Y Matrix of endogenous variables (T x N). 
#'          ORDER MUST BE: [1] OIS_3M, [2] STOXX50, [3] Log_IP, [4] Log_Oil, [5] Log_HICP, [6] i_ECB
#' @param p Integer: Number of VAR lags 
#' @param lambda1 Numeric: Overall Tightness (e.g., 0.20)
#' @param lambda2 Numeric: Cross-Variable Tightness (e.g., 0.50)
#' @param lambda3 Numeric: Lag Decay (e.g., 2.0)
#' @param prior_means Vector of length N: 0 for surprises, 0.9 for levels
#' @param n_draws Integer: Number of valid structural draws to keep
#' @param burn_in Integer: Number of initial Gibbs draws to discard
#' @param horizon Integer: IRF forecast horizon (months)
#' @return An array of valid structural IRFs
estimate_jk2020_bvar <- function(Y, p = 6, 
                                 lambda1 = 0.20, lambda2 = 0.50, lambda3 = 2.0,
                                 prior_means = c(0, 0, 0.9, 0.9, 0.9, 0.9), 
                                 n_draws = 1000, burn_in = 500, horizon = 24) {
  
  T_total <- nrow(Y)
  N <- ncol(Y)
  
  # ----------------------------------------------------------------------------
  # 1. Construct Lagged Matrices (Linear, unweighted)
  # ----------------------------------------------------------------------------
  build_var_matrices <- function(data, lags) {
    T_obs <- nrow(data)
    Y_trim <- data[(lags + 1):T_obs, , drop = FALSE]
    X_mat <- matrix(1, nrow = T_obs - lags, ncol = 1) 
    for (l in 1:lags) {
      X_mat <- cbind(X_mat, data[(lags - l + 1):(T_obs - l), , drop = FALSE])
    }
    return(list(Y = Y_trim, X = X_mat))
  }
  
  var_data <- build_var_matrices(Y, p)
  Y_t <- var_data$Y
  X_t <- var_data$X
  
  T_eff <- nrow(Y_t)
  K <- ncol(X_t)
  
  XtX <- crossprod(X_t)
  XtY <- crossprod(X_t, Y_t)
  
  # ----------------------------------------------------------------------------
  # 2. Build Independent Normal-Wishart Prior (Litterman 1986 Framework)
  # ----------------------------------------------------------------------------
  # Estimate univariate AR(1) residual variances (s_sq)
  s_sq <- numeric(N)
  for (i in 1:N) {
    y_i <- Y_t[, i]
    x_i <- cbind(1, dplyr::lag(y_i, 1))
    valid <- complete.cases(x_i, y_i)
    b_i <- solve(t(x_i[valid, ]) %*% x_i[valid, ]) %*% t(x_i[valid, ]) %*% y_i[valid]
    s_sq[i] <- sum((y_i[valid] - x_i[valid, ] %*% b_i)^2) / (sum(valid) - 2)
  }
  
  # Prior Mean Vector for vectorized B (KN x 1)
  B0_mat <- matrix(0, nrow = K, ncol = N)
  for (i in 1:N) { B0_mat[1 + i, i] <- prior_means[i] }
  b0 <- as.vector(B0_mat)
  
  # Prior Covariance V0 for vectorized B (KN x KN Diagonal)
  v0_diag <- numeric(K * N)
  idx <- 1
  for (i in 1:N) { 
    v0_diag[idx] <- 1e4 * s_sq[i] # Intercept
    idx <- idx + 1
    for (l in 1:p) {
      for (j in 1:N) { 
        if (i == j) {
          # Own Lag Tightness
          v0_diag[idx] <- (lambda1 / (l^lambda3))^2
        } else {
          # Cross-Variable Tightness (lambda2)
          v0_diag[idx] <- ((lambda1 * lambda2) / (l^lambda3))^2 * (s_sq[i] / s_sq[j])
        }
        idx <- idx + 1
      }
    }
  }
  V0_inv <- diag(1 / v0_diag)
  
  # Prior for Sigma (Inverse-Wishart)
  nu0 <- N + 2
  S0  <- diag(s_sq) * (nu0 - N - 1)
  nu_bar <- nu0 + T_eff
  
  # ----------------------------------------------------------------------------
  # Helper: Compute IRFs
  # ----------------------------------------------------------------------------
  compute_irf <- function(B_mat, A0, lags, horiz) {
    A_lags <- t(B_mat[-1, , drop = FALSE]) 
    irf_array <- array(0, dim = c(horiz + 1, N, N))
    irf_array[1, , ] <- A0
    if (horiz >= 1) {
      for (h in 1:horiz) {
        resp <- matrix(0, nrow = N, ncol = N)
        for (l in 1:min(h, lags)) {
          A_l <- A_lags[, ((l - 1) * N + 1):(l * N)]
          resp <- resp + A_l %*% irf_array[h - l + 1, , ]
        }
        irf_array[h + 1, , ] <- resp
      }
    }
    return(irf_array)
  }
  
  # ----------------------------------------------------------------------------
  # 3. Main MCMC Loop (Gibbs Sampler -> Sign Restriction Filter)
  # ----------------------------------------------------------------------------
  irf_draws <- array(0, dim = c(n_draws, horizon + 1, N, N))
  
  cat("Starting Original JK (2020) Gibbs Sampler...\n")
  
  # Initialize Sigma with OLS estimate
  B_ols <- solve(XtX) %*% XtY
  Sigma_current <- crossprod(Y_t - X_t %*% B_ols) / (T_eff - K)
  
  valid_draws <- 0
  step <- 0
  
  while (valid_draws < n_draws) {
    step <- step + 1
    
    # A. Draw B conditional on Sigma
    Sigma_inv <- solve(Sigma_current)
    V_bar_inv <- V0_inv + kronecker(Sigma_inv, XtX)
    V_bar     <- chol2inv(chol(V_bar_inv))
    
    # Vectorized OLS logic for the mean
    term2 <- as.vector(XtY %*% Sigma_inv)
    b_bar <- V_bar %*% (V0_inv %*% b0 + term2)
    
    # Ensure non-explosive draws (Truncated Normal)
    is_stable <- FALSE
    while (!is_stable) {
      B_vec <- b_bar + t(chol(V_bar)) %*% rnorm(K * N)
      B_mat <- matrix(B_vec, nrow = K, ncol = N)
      
      # Check max eigenvalue
      A_lags <- t(B_mat[-1, , drop = FALSE])
      companion <- matrix(0, nrow = N * p, ncol = N * p)
      companion[1:N, ] <- A_lags
      if (p > 1) { companion[(N + 1):(N * p), 1:(N * (p - 1))] <- diag(N * (p - 1)) }
      is_stable <- (max(abs(eigen(companion, only.values = TRUE)$values)) < 0.998)
    }
    
    # B. Draw Sigma conditional on B
    E_mat <- Y_t - X_t %*% B_mat
    S_bar <- S0 + crossprod(E_mat)
    Sigma_current <- solve(rWishart(1, df = nu_bar, Sigma = chol2inv(chol(S_bar)))[,,1])
    
    # C. After burn-in, apply Haar rotation and JK (2020) Sign Restrictions
    if (step > burn_in) {
      P_chol <- t(chol(Sigma_current))
      
      X_rand <- matrix(rnorm(N * N), nrow = N, ncol = N)
      qr_decomp <- qr(X_rand)
      Q <- qr.Q(qr_decomp) %*% diag(sign(diag(qr.R(qr_decomp))))
      
      A0 <- P_chol %*% Q
      
      # Pure MP: OIS > 0, Stocks < 0
      is_pure_mp <- (A0[1, 1] > 0) && (A0[2, 1] < 0)
      # CBI: OIS > 0, Stocks > 0
      is_cbi     <- (A0[1, 2] > 0) && (A0[2, 2] > 0)
      
      if (!is_pure_mp || !is_cbi) {
        if ((A0[1, 1] < 0) && (A0[2, 1] > 0)) { A0[, 1] <- -A0[, 1]; is_pure_mp <- TRUE }
        if ((A0[1, 2] < 0) && (A0[2, 2] < 0)) { A0[, 2] <- -A0[, 2]; is_cbi <- TRUE }
      }
      
      if (is_pure_mp && is_cbi) {
        target_impact <- 25   # 25bp OIS surprise, adjust to your units
        A0[, 1] <- A0[, 1] * (target_impact / A0[1, 1])
        A0[, 2] <- A0[, 2] * (target_impact / A0[1, 2])
        valid_draws <- valid_draws + 1
        irf_draws[valid_draws, , , ] <- compute_irf(B_mat, A0, p, horizon)
        
        if (valid_draws %% 100 == 0) {
          cat(sprintf("-> Valid Structural Draws: %d / %d (Total Gibbs Steps: %d)\n", 
                      valid_draws, n_draws, step))
        }
      }
    }
  }
  
  cat("Estimation Complete.\n")
  return(irf_draws)
}



