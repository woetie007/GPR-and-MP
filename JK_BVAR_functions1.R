# ==============================================================================
# JOINT 2D GRID SEARCH FOR BVAR LAG ORDER (p) AND OVERALL TIGHTNESS (lambda1)
# ==============================================================================

#' Compute Closed-Form Log-MDD over a 2D Grid of Lags and Overall Tightness
#'
#' @param Y Matrix of endogenous variables (T x N).
#'          ORDER: [1] MP_shocks, [2] stocks_shocks, [3] i_ECB, [4] IP, [5] pi, [6] Oil
#' @param lags_to_test Integer vector of lag orders to evaluate (e.g., 1:12)
#' @param lambda1_grid Numeric vector of candidate overall tightness values
#' @param prior_means Vector of length N for own 1st-lag prior means.
#'                    Default sets 0 for high-frequency surprises, 0.9 for macro variables.
#' @param lambda3 Numeric: Lag decay rate (standard default: 2.0 for harmonic decay)
#' @return A list containing the full grid results, optimal pair, and a 2D heatmap ggplot
compute_bvar_2d_grid_mdd <- function(
    Y, 
    lags_to_test = 1:12, 
    lambda1_grid = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50, 1.00),
    prior_means  = c(0, 0, 0.9, 0.9, 0.9, 0.9),
    lambda3      = 2.0
) {
  
  T_total <- nrow(Y)
  N <- ncol(Y)
  p_max <- max(lags_to_test)
  
  # ----------------------------------------------------------------------------
  # SAMPLE ALIGNMENT: Ensure identical effective sample across ALL models
  # ----------------------------------------------------------------------------
  T_eff  <- T_total - p_max
  Y_eval <- Y[(p_max + 1):T_total, , drop = FALSE]
  
  # Helper: Multivariate log-Gamma function ln(Gamma_N(x))
  lmvgamma <- function(x, d) {
    0.25 * d * (d - 1) * log(pi) + sum(lgamma(x + (1 - 1:d) / 2))
  }
  
  # Helper: Residual variance s_i^2 from AR(1) for Minnesota scaling
  s_sq <- numeric(N)
  for (i in 1:N) {
    y_i <- Y[-1, i]
    x_i <- cbind(1, Y[-T_total, i])
    b_i <- solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_i
    e_i <- y_i - x_i %*% b_i
    s_sq[i] <- sum(e_i^2) / (length(y_i) - 2)
  }
  
  # ----------------------------------------------------------------------------
  # 2D Evaluation Loop
  # ----------------------------------------------------------------------------
  results_list <- list()
  
  for (p in lags_to_test) {
    K <- 1 + N * p  # Intercept + N*p lags
    
    # Build regressor matrix X_eval aligned to Y_eval
    X_eval <- matrix(1, nrow = T_eff, ncol = 1)
    for (l in 1:p) {
      l_start <- (p_max - l + 1)
      l_end   <- (T_total - l)
      X_eval  <- cbind(X_eval, Y[l_start:l_end, , drop = FALSE])
    }
    
    XtX <- t(X_eval) %*% X_eval
    
    for (lam1 in lambda1_grid) {
      
      # Construct Conjugate Minnesota Prior (B0, V0, S0, nu0)
      nu0 <- N + 2
      S0  <- diag(s_sq) * (nu0 - N - 1)
      
      # Prior mean B0 (K x N)
      B0 <- matrix(0, nrow = K, ncol = N)
      for (i in 1:N) {
        B0[1 + i, i] <- prior_means[i]
      }
      
      # Diagonal prior covariance V0 (K x K)
      v0_diag <- numeric(K)
      v0_diag[1] <- 1e4  # Uninformative intercept
      
      idx <- 2
      for (l in 1:p) {
        for (j in 1:N) {
          # Minnesota decay formula using candidate lam1
          v0_diag[idx] <- (lam1 / ((l^lambda3) * sqrt(s_sq[j])))^2
          idx <- idx + 1
        }
      }
      
      V0     <- diag(v0_diag)
      V0_inv <- diag(1 / v0_diag)
      
      # Compute posterior hyperparameters
      Vn     <- solve(V0_inv + XtX)
      Vn_inv <- V0_inv + XtX
      
      Bn   <- Vn %*% (V0_inv %*% B0 + t(X_eval) %*% Y_eval)
      Sn   <- S0 + t(Y_eval) %*% Y_eval + t(B0) %*% V0_inv %*% B0 - t(Bn) %*% Vn_inv %*% Bn
      nu_n <- nu0 + T_eff
      
      # Exact Anallytical Log Marginal Data Density
      log_det_Vn <- as.numeric(determinant(Vn, logarithm = TRUE)$modulus)
      log_det_V0 <- as.numeric(determinant(V0, logarithm = TRUE)$modulus)
      log_det_S0 <- as.numeric(determinant(S0, logarithm = TRUE)$modulus)
      log_det_Sn <- as.numeric(determinant(Sn, logarithm = TRUE)$modulus)
      
      log_mdd <- -0.5 * T_eff * N * log(pi) +
        0.5 * N * (log_det_Vn - log_det_V0) +
        0.5 * nu0 * log_det_S0 -
        0.5 * nu_n * log_det_Sn +
        lmvgamma(nu_n / 2, N) - 
        lmvgamma(nu0 / 2, N)
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Lags    = p,
        Lambda1 = lam1,
        Log_MDD = log_mdd
      )
    }
  }
  
  df_res <- do.call(rbind, results_list)
  
  # ----------------------------------------------------------------------------
  # Compute Posterior Model Probabilities (PMP) over the 2D Grid
  # ----------------------------------------------------------------------------
  max_mdd <- max(df_res$Log_MDD)
  df_res$PMP <- exp(df_res$Log_MDD - max_mdd) / sum(exp(df_res$Log_MDD - max_mdd))
  
  # Identify the optimal pair (p*, lambda1*)
  best_row <- df_res[which.max(df_res$Log_MDD), ]
  opt_p    <- best_row$Lags
  opt_lam  <- best_row$Lambda1
  
  # ----------------------------------------------------------------------------
  # Generate 2D Heatmap
  # ----------------------------------------------------------------------------
  p_heatmap <- ggplot(df_res, aes(x = factor(Lags), y = factor(Lambda1), fill = PMP)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # Outline the winning model cell in distinct yellow/black
    geom_tile(
      data = best_row, 
      aes(x = factor(Lags), y = factor(Lambda1)), 
      fill = NA, color = "#ff7f00", linewidth = 1.8
    ) +
    geom_text(
      aes(label = sprintf("%.1f%%", PMP * 100)), 
      color = ifelse(df_res$PMP > 0.15, "white", "black"), 
      size  = 3.2, 
      fontface = "bold"
    ) +
    scale_fill_gradient(
      low = "#f7fcf0", high = "#084081", 
      labels = scales::percent_format(accuracy = 1),
      name = "Posterior\nProbability"
    ) +
    labs(
      title    = "2D Grid Search: Joint Selection of Lag Order & Overall Tightness",
      subtitle = sprintf("Optimal Pair Selected by MDD: p* = %d lags, lambda1* = %.2f (PMP = %.1f%%)", 
                         opt_p, opt_lam, best_row$PMP * 100),
      x        = "Candidate Lag Order (p)",
      y        = "Overall Tightness Hyperparameter (lambda1)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold"),
      axis.title   = element_text(face = "bold"),
      panel.grid   = element_blank()
    )
  
  return(list(
    summary_table = df_res,
    optimal_pair  = list(p = opt_p, lambda1 = opt_lam, pmp = best_row$PMP),
    plot          = p_heatmap
  ))
}
# ==============================================================================
# OPTIMAL GRID SEARCH FOR TRANSITION PARAMETERS (gamma, c, d) VIA POOLED MDD
# ==============================================================================
#' Compute Closed-Form Log-MDD for a Weighted Regime (WLS-MDD)
#'
#' @param Y Matrix of endogenous variables (T x N)
#' @param p Integer: VAR lag order
#' @param weights Vector of length T-p: Logistic transition weights F(z_t) or 1-F(z_t)
#' @param prior_means Vector of length N for own 1st-lag prior means
#' @param lambda1 Numeric: Overall tightness parameter
#' @param lambda3 Numeric: Lag decay rate
#' @return Numeric: Log Marginal Data Density for the weighted subsample
compute_weighted_log_mdd <- function(Y, p = 6, weights,
                                     prior_means = c(0, 0, 0.9, 0.9, 0.9, 0.9),
                                     lambda1 = 0.2, lambda3 = 2.0) {
  
  T_total <- nrow(Y)
  N       <- ncol(Y)
  T_eff   <- T_total - p
  
  # Trim sample to matching evaluation window
  Y_eval <- Y[(p + 1):T_total, , drop = FALSE]
  
  # Helper: Multivariate log-Gamma function ln(Gamma_N(x))
  lmvgamma <- function(x, d) {
    0.25 * d * (d - 1) * log(pi) + sum(lgamma(x + (1 - 1:d) / 2))
  }
  
  # Helper: Residual variance s_i^2 from AR(1) for Minnesota scaling
  s_sq <- numeric(N)
  for (i in 1:N) {
    y_i <- Y[-1, i]
    x_i <- cbind(1, Y[-T_total, i])
    b_i <- solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_i
    e_i <- y_i - x_i %*% b_i
    s_sq[i] <- sum(e_i^2) / (length(y_i) - 2)
  }
  
  # Build regressor matrix X_eval
  X_eval <- matrix(1, nrow = T_eff, ncol = 1)
  for (l in 1:p) {
    l_start <- (p - l + 1)
    l_end   <- (T_total - l)
    X_eval  <- cbind(X_eval, Y[l_start:l_end, , drop = FALSE])
  }
  
  # Effective sample size and WLS weighting matrix W
  eff_T  <- sum(weights)
  W_sqrt <- diag(sqrt(weights))
  
  Y_w <- W_sqrt %*% Y_eval
  X_w <- W_sqrt %*% X_eval
  
  # Construct Conjugate Minnesota Prior (B0, V0, S0, nu0)
  nu0 <- N + 2
  S0  <- diag(s_sq) * (nu0 - N - 1)
  
  K  <- 1 + N * p
  B0 <- matrix(0, nrow = K, ncol = N)
  for (i in 1:N) {
    B0[1 + i, i] <- prior_means[i]
  }
  
  v0_diag <- numeric(K)
  v0_diag[1] <- 1e4  # Uninformative intercept
  
  idx <- 2
  for (l in 1:p) {
    for (j in 1:N) {
      v0_diag[idx] <- (lambda1 / ((l^lambda3) * sqrt(s_sq[j])))^2
      idx <- idx + 1
    }
  }
  V0     <- diag(v0_diag)
  V0_inv <- diag(1 / v0_diag)
  
  # Compute posterior hyperparameters under WLS
  XtX_w  <- t(X_w) %*% X_w
  Vn     <- solve(V0_inv + XtX_w)
  Vn_inv <- V0_inv + XtX_w
  
  Bn   <- Vn %*% (V0_inv %*% B0 + t(X_w) %*% Y_w)
  Sn   <- S0 + t(Y_w) %*% Y_w + t(B0) %*% V0_inv %*% B0 - t(Bn) %*% Vn_inv %*% Bn
  nu_n <- nu0 + eff_T
  
  # Exact Anallytical Log Marginal Data Density
  log_det_Vn <- as.numeric(determinant(Vn, logarithm = TRUE)$modulus)
  log_det_V0 <- as.numeric(determinant(V0, logarithm = TRUE)$modulus)
  log_det_S0 <- as.numeric(determinant(S0, logarithm = TRUE)$modulus)
  log_det_Sn <- as.numeric(determinant(Sn, logarithm = TRUE)$modulus)
  
  log_mdd <- -0.5 * eff_T * N * log(pi) +
    0.5 * N * (log_det_Vn - log_det_V0) +
    0.5 * nu0 * log_det_S0 -
    0.5 * nu_n * log_det_Sn +
    lmvgamma(nu_n / 2, N) - 
    lmvgamma(nu0 / 2, N)
  
  return(log_mdd)
}


select_optimal_st_parameters <- function(Y, GPR, p = 6, 
                                         d_grid     = c(1, 2),
                                         gamma_pcts = c(0.40, 0.50, 0.60, 0.70, 0.80),
                                         c_grid     = c(1.5, 3.0, 5.0, 7.0)) {
  
  T_total <- nrow(Y)
  N       <- ncol(Y)
  results <- list()
  
  for (d_cand in d_grid) {
    # Aligned exogenous transition variable
    z_t     <- GPR[(p + 1 - d_cand):(T_total - d_cand)]
    sigma_z <- sd(z_t)
    
    for (pct in gamma_pcts) {
      gamma_cand <- quantile(GPR, probs = pct, na.rm = TRUE)
      
      for (c_cand in c_grid) {
        
        # Compute logistic transition weights F(z_t)
        F_high <- 1 / (1 + exp(-c_cand * (z_t - gamma_cand) / sigma_z))
        F_low  <- 1 - F_high
        
        # Compute closed-form Weighted Log-MDD for Low and High regimes
        # (Using the analytical conjugate WLS-MDD formula from our previous step)
        log_mdd_low  <- compute_weighted_log_mdd(Y, p, weights = F_low)
        log_mdd_high <- compute_weighted_log_mdd(Y, p, weights = F_high)
        
        # Total pooled model evidence
        pooled_log_mdd <- log_mdd_low + log_mdd_high
        
        results[[length(results) + 1]] <- data.frame(
          d              = d_cand,
          Gamma_Pct      = pct,
          Gamma_Val      = as.numeric(gamma_cand),
          c              = c_cand,
          Pooled_Log_MDD = pooled_log_mdd
        )
      }
    }
  }
  
  df_res <- do.call(rbind, results)
  
  # Return the optimal triplet with the highest evidence score
  best_param <- df_res[which.max(df_res$Pooled_Log_MDD), ]
  
  return(list(
    grid_results = df_res,
    optimal      = best_param
  ))
}

# ==============================================================================
# OPTIMIZE MINNESOTA PRIORS (lambda1, lambda3) VIA GIANNONE, LENZA & PRIMICERI (2015)
# ==============================================================================

library(ggplot2)
library(dplyr)

#' Optimize Minnesota Prior Hyperparameters for a Smooth Transition BVAR
#'
#' @param Y Matrix of endogenous variables (T x N)
#' @param GPR Exogenous Geopolitical Risk vector of length T
#' @param p Integer: Selected VAR lag order (e.g., 6)
#' @param d Integer: Delay parameter for GPR
#' @param gamma Numeric: Threshold location parameter
#' @param c Numeric: Transition speed parameter
#' @param lambda1_grid Vector of candidate overall tightness values
#' @param lambda3_grid Vector of candidate lag decay rate values
#' @param prior_means Vector of length N for own 1st-lag prior means
#' @return List containing summary table, optimal parameters, and a heatmap ggplot
optimize_minnesota_priors <- function(
    Y, GPR, p = 6, d = 1, 
    gamma = median(GPR, na.rm = TRUE), c = 10.0,
    lambda1_grid = c(0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.50),
    lambda3_grid = c(1.0, 1.5, 2.0, 3.0),
    prior_means  = c(0, 0, 0.9, 0.9, 0.9, 0.9)
) {
  
  T_total <- nrow(Y)
  N       <- ncol(Y)
  
  # Compute logistic transition weights F_low and F_high
  z_t     <- GPR[(p + 1 - d):(T_total - d)]
  sigma_z <- sd(z_t, na.rm = TRUE)
  F_high  <- 1 / (1 + exp(-c * (z_t - gamma) / sigma_z))
  F_low   <- 1 - F_high
  
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
  
  # Helper: Analytical Closed-Form Log-MDD for a single weighted regime
  eval_regime_mdd <- function(weights, lam1, lam3) {
    T_eff <- sum(weights)
    Y_eval <- Y[(p + 1):T_total, , drop = FALSE]
    
    X_eval <- matrix(1, nrow = nrow(Y_eval), ncol = 1)
    for (l in 1:p) {
      l_start <- (p - l + 1)
      l_end   <- (T_total - l)
      X_eval  <- cbind(X_eval, Y[l_start:l_end, , drop = FALSE])
    }
    
    W_sqrt <- diag(sqrt(weights))
    Y_w    <- W_sqrt %*% Y_eval
    X_w    <- W_sqrt %*% X_eval
    
    K  <- 1 + N * p
    nu0 <- N + 2
    S0  <- diag(s_sq) * (nu0 - N - 1)
    
    B0 <- matrix(0, nrow = K, ncol = N)
    for (i in 1:N) { B0[1 + i, i] <- prior_means[i] }
    
    v0_diag <- numeric(K)
    v0_diag[1] <- 1e4
    idx <- 2
    for (l in 1:p) {
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
    nu_n <- max(round(nu0 + eff_T), N + 2)
    
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
  
  # Grid Search Loop over lambda1 and lambda3
  results_list <- list()
  for (lam1 in lambda1_grid) {
    for (lam3 in lambda3_grid) {
      mdd_low  <- eval_regime_mdd(F_low,  lam1, lam3)
      mdd_high <- eval_regime_mdd(F_high, lam1, lam3)
      pooled   <- mdd_low + mdd_high
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Lambda1        = lam1,
        Lambda3        = lam3,
        Log_MDD_Low    = mdd_low,
        Log_MDD_High   = mdd_high,
        Pooled_Log_MDD = pooled
      )
    }
  }
  
  df_res <- do.call(rbind, results_list)
  
  # Compute Posterior Model Probabilities (PMP) via Softmax
  max_mdd <- max(df_res$Pooled_Log_MDD)
  df_res$PMP <- exp(df_res$Pooled_Log_MDD - max_mdd) / sum(exp(df_res$Pooled_Log_MDD - max_mdd))
  
  best_row <- df_res[which.max(df_res$Pooled_Log_MDD), ]
  
  # Heatmap Plot
  p_chart <- ggplot(df_res, aes(x = factor(Lambda1), y = factor(Lambda3), fill = PMP)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_tile(data = best_row, aes(x = factor(Lambda1), y = factor(Lambda3)), 
              fill = NA, color = "#ff7f00", linewidth = 1.8) +
    geom_text(aes(label = sprintf("%.1f%%", PMP * 100)), 
              color = ifelse(df_res$PMP > 0.15, "white", "black"), 
              size  = 3.3, fontface = "bold") +
    scale_fill_gradient(low = "#f7fcf0", high = "#084081", 
                        labels = scales::percent_format(accuracy = 1),
                        name = "Posterior\nProbability") +
    labs(
      title    = "Optimal Minnesota Prior Selection (Giannone, Lenza & Primiceri, 2015)",
      subtitle = sprintf("Optimal Priors Chosen: lambda1* = %.2f, lambda3* = %.1f (PMP = %.1f%%)",
                         best_row$Lambda1, best_row$Lambda3, best_row$PMP * 100),
      x        = "Overall Tightness (lambda1)",
      y        = "Lag Decay Rate (lambda3)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())
  
  return(list(
    summary_table = df_res,
    optimal       = best_row,
    plot          = p_chart
  ))
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
# CHECK VAR STABILITY (REJECT EXPLOSIVE DRAW IF MAX EIGENVALUE >= 0.998)
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
  return(max_eig < 0.998) 
}

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
# SMOOTH TRANSITION BAYESIAN VAR WITH JK (2020) SIGN RESTRICTIONS
# ==============================================================================

library(MASS)

#' Estimate Smooth Transition JK (2020) Sign-Restricted BVAR
#'
#' @param Y Matrix of endogenous variables (T x N). 
#'          ORDER MUST BE: [1] MP_shocks, [2] stocks_shocks, [3] i_ECB, [4] IP, [5] pi, [6] Oil
#' @param GPR Exogenous Geopolitical Risk vector of length T
#' @param p Integer: Number of VAR lags
#' @param d Integer: Lag delay for GPR transition variable (default d = 1)
#' @param gamma Numeric: Threshold location parameter for GPR
#' @param c Numeric: Smoothness parameter (c > 0; higher = sharper transition)
#' @param n_draws Integer: Number of successful posterior structural draws to keep
#' @param horizon Integer: IRF forecast horizon (months)
#' @return A list containing arrays of IRFs for Low-GPR and High-GPR regimes
estimate_st_jk2020_bvar <- function(Y, GPR, p = 6, d = 1, 
                                    gamma = median(GPR), c = 10.0, 
                                    lambda1     = 0.12,
                                    lambda3     = 2.0,
                                    lambda4     = 1.0,  
                                    lambda5     = 1.0,  
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
  
  # ----------------------------------------------------------------------------
  # Logistic transition function F(z_t) IN [0, 1]
  # ----------------------------------------------------------------------------
  z_t <- GPR[(p + 1 - d):(T_total - d)]
  sigma_z <- sd(z_t)
  F_high <- 1 / (1 + exp(-c * (z_t - gamma) / sigma_z))
  F_low  <- 1 - F_high
  
  # ----------------------------------------------------------------------------
  # Helper: Weighted Normal-Inverse-Wishart (WLS) Posterior Draw
  # ----------------------------------------------------------------------------
  draw_weighted_niw_posterior <- function(Y_sub, X_sub, weights, 
                                          lambda1 = 0.15, lambda3 = 2.0,
                                          lambda4 = 1.0, lambda5 = 1.0, 
                                          prior_means = c(0, 0, 0.9, 0.9, 0.9, 0.9)) {
    T_sub <- nrow(Y_sub)
    K     <- ncol(X_sub)
    eff_T <- sum(weights)
    
    # 1. Calculate AR(1) residual variances (s_sq) for Minnesota scale standardization
    s_sq <- numeric(N)
    for (i in 1:N) {
      y_i <- Y_sub[, i]
      x_i <- cbind(1, dplyr::lag(y_i, 1))
      valid <- complete.cases(x_i, y_i)
      b_i <- solve(t(x_i[valid, ]) %*% x_i[valid, ]) %*% t(x_i[valid, ]) %*% y_i[valid]
      e_i <- y_i[valid] - x_i[valid, ] %*% b_i
      s_sq[i] <- sum(e_i^2) / (sum(valid) - 2)
    }
    
    # 2. Build Minnesota Prior Mean Matrix B0 (K x N)
    B0 <- matrix(0, nrow = K, ncol = N)
    for (i in 1:N) { B0[1 + i, i] <- prior_means[i] }
    
    # 3. Build Minnesota Prior Covariance V0 (K x K diagonal matrix)
    v0_diag <- numeric(K)
    v0_diag[1] <- 1e4 
    idx <- 2
    for (l in 1:p) {
      for (j in 1:N) {
        v0_diag[idx] <- (lambda1 / ((l^lambda3) * sqrt(s_sq[j])))^2
        idx <- idx + 1
      }
    }
    V0_inv <- diag(1 / v0_diag)
    
    # 4. Apply smooth WLS weights to data
    W_sqrt <- diag(sqrt(weights))
    Y_w    <- W_sqrt %*% Y_sub
    X_w    <- W_sqrt %*% X_sub
    
    # ========================================================================
    # 5. SIMS-ZHA DUMMY OBSERVATION PRIORS (SOC & DIO)
    # ========================================================================
    y_bar <- colMeans(Y_sub)
    y_bar[1:2] <- 0  # CRITICAL: No unit root priors on surprises
    
    Y_soc <- diag(y_bar) / lambda4
    X_soc <- matrix(0, nrow = N, ncol = K)
    for (l in 1:p) {
      col_start <- 2 + (l - 1) * N
      col_end   <- 1 + l * N
      X_soc[, col_start:col_end] <- diag(y_bar) / lambda4
    }
    
    Y_dio <- matrix(y_bar / lambda5, nrow = 1, ncol = N)
    X_dio <- matrix(0, nrow = 1, ncol = K)
    X_dio[1, 1] <- 1 / lambda5 
    for (l in 1:p) {
      col_start <- 2 + (l - 1) * N
      col_end   <- 1 + l * N
      X_dio[1, col_start:col_end] <- y_bar / lambda5
    }
    
    Y_w <- rbind(Y_w, Y_soc, Y_dio)
    X_w <- rbind(X_w, X_soc, X_dio)
    
    # 6. Compute Conjugate Posterior Hyperparameters
    XtX_w  <- t(X_w) %*% X_w
    Vn     <- solve(V0_inv + XtX_w)
    Vn_inv <- V0_inv + XtX_w
    Bn   <- Vn %*% (V0_inv %*% B0 + t(X_w) %*% Y_w)
    
    nu0  <- N + 2
    S0   <- diag(s_sq) * (nu0 - N - 1)
    Sn   <- S0 + crossprod(Y_w) + t(B0) %*% V0_inv %*% B0 - t(Bn) %*% Vn_inv %*% Bn
    nu_n <- round(eff_T) + nu0
    
    C_chol <- chol(Sn)
    Z_w    <- matrix(rnorm(nu_n * N), nrow = nu_n, ncol = N)
    Sigma_draw <- t(C_chol) %*% solve(t(Z_w) %*% Z_w) %*% C_chol
    
    Z_b    <- matrix(rnorm(K * N), nrow = K, ncol = N)
    B_draw <- Bn + t(chol(Vn)) %*% Z_b %*% chol(Sigma_draw)
    
    return(list(B = B_draw, Sigma = Sigma_draw))
  }
  
  # ----------------------------------------------------------------------------
  # Helper: Compute IRFs
  # ----------------------------------------------------------------------------
  compute_irf <- function(B, A0, lags, horiz) {
    A_lags <- t(B[-1, , drop = FALSE]) 
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
      
      if (total_attempts %% 100 == 0) {
        cat(sprintf("-> Attempts: %d | Accepted: %d (Hit Rate: %.1f%%)\n", 
                    total_attempts, draw_count, (draw_count / total_attempts) * 100))
      }
      
      is_stable <- FALSE
      while (!is_stable) {
        niw <- draw_weighted_niw_posterior(
          Y_t, X_t, regime_weights, 
          lambda1     = lambda1, 
          lambda3     = lambda3,
          lambda4     = lambda4,
          lambda5     = lambda5,
          prior_means = prior_means
        )
        is_stable <- is_stable_var(niw$B, p, N)
      }
      P_chol <- t(chol(niw$Sigma)) 
      
      X_rand <- matrix(rnorm(N * N), nrow = N, ncol = N)
      qr_decomp <- qr(X_rand)
      Q <- qr.Q(qr_decomp)
      R_qr <- qr.R(qr_decomp)
      Q <- Q %*% diag(sign(diag(R_qr))) 
      
      A0 <- P_chol %*% Q
      
      is_pure_mp <- (A0[1, 1] > 0) && (A0[2, 1] < 0) && (A0[3, 1] >= 0.02)
      is_cbi     <- (A0[1, 2] > 0) && (A0[2, 2] > 0) && (A0[3, 2] >= 0.02)
      
      if (!is_pure_mp || !is_cbi) {
        if ((A0[1, 1] < 0) && (A0[2, 1] > 0) && (A0[3, 1] <= -0.02)) {
          A0[, 1] <- -A0[, 1]
          is_pure_mp <- TRUE
        }
        if ((A0[1, 2] < 0) && (A0[2, 2] < 0) && (A0[3, 2] <= -0.02)) {
          A0[, 2] <- -A0[, 2]
          is_cbi <- TRUE
        }
      }
      
      if (is_pure_mp && is_cbi) {
        irf_test <- compute_irf(niw$B, A0, p, horiz = 2)
        
        is_persistent <- (irf_test[2, 3, 1] >= -0.01)
        
        if (is_persistent) {
          draw_count <- draw_count + 1
          irf_mat <- compute_irf(niw$B, A0, p, horizon)
          
          if (regime_name == "Low") {
            irf_draws_low[draw_count, , , ] <- irf_mat
          } else {
            irf_draws_high[draw_count, , , ] <- irf_mat
          }
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