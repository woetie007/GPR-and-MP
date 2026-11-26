
#==================================
# Christopher Sims' IRF Calculator
#==================================
impulsdtrf <- function(B, smat, nstep) {
  # B: N x N x nlags array of reduced-form coefficients
  # smat: Initial impact matrix (typically transposed Cholesky factor)
  # Returns: N x N x nstep array of impulse responses
  
  dims <- dim(B)
  neq <- dims[1]
  nvar <- dims[2]
  nlag <- dims[3]
  
  response <- array(0, dim = c(neq, nvar, nstep))
  response[, , 1] <- t(smat)
  
  if (nstep > 1) {
    for (it in 2:nstep) {
      for (ilag in 1:min(nlag, it - 1)) {
        response[, , it] <- response[, , it] + B[, , ilag] %*% response[, , it - ilag]
      }
    }
  }
  return(response)
}

#=======================
# The VAR Gibbs Sampler
#=======================
VAR_withiid <- function(Y, prior, gssettings) {
  T_full <- nrow(Y)
  N <- ncol(Y)
  Nm <- prior$Nm
  Ny <- N - Nm
  P <- prior$lags
  T_eff <- T_full - P
  
  K <- P * N + 1 # Number of regressors (lags + intercept)
  
  # Prepare data matrices
  Y_t <- Y[(P + 1):T_full, , drop = FALSE]
  X_t <- matrix(1, nrow = T_eff, ncol = 1)
  for (l in 1:P) {
    X_t <- cbind(X_t, Y[(P + 1 - l):(T_full - l), , drop = FALSE])
  }
  
  # Calculate AR(1) residual standard errors for the prior scale (sigma)
  prior_sigma <- numeric(N)
  for (n in 1:N) {
    y_n <- Y_t[, n]
    x_n <- cbind(1, Y[(P):(T_full - 1), n])
    valid <- complete.cases(x_n, y_n)
    b_n <- solve(t(x_n[valid, ]) %*% x_n[valid, ]) %*% t(x_n[valid, ]) %*% y_n[valid]
    prior_sigma[n] <- sd(y_n[valid] - x_n[valid, ] %*% b_n)
  }
  
  # Build prior covariance matrix (Q0) for Minnesota prior
  temp1 <- rep((1:P)^(-prior$minnesota_decay), each = N)
  temp2 <- rep(prior_sigma, times = P)
  
  Q0_list <- list()
  for(i in 1:N) {
    temp3 <- 1 / prior_sigma[i]
    Q0_col <- prior$minnesota_tightness * temp1 * temp2 * temp3
    Q0_list[[i]] <- c(1e5, Q0_col) # 1e5 for the intercept
  }
  Q0_mat <- do.call(cbind, Q0_list)
  Q0_mat <- Q0_mat^2
  
  # Drop equations for the surprises (m)
  Q0_mat <- Q0_mat[, (Nm + 1):N, drop = FALSE] 
  
  Q_inv <- diag(1 / as.vector(Q0_mat))
  
  # Build prior mean matrix (B)
  prior_B <- matrix(0, nrow = K, ncol = N)
  mvec <- prior$minnesota_mvector
  for (n in 1:N) {
    prior_B[1 + n, n] <- mvec[n] 
  }
  prior_B <- prior_B[, (Nm + 1):N, drop = FALSE]
  prior_QinvB <- Q_inv %*% as.vector(prior_B)
  
  # Prior for Sigma (Inverse-Wishart)
  v_prior <- N + 2
  S_prior <- diag(prior_sigma^2 * (v_prior - N - 1))
  v_post <- T_eff + v_prior
  
  # Initialize Gibbs Sampler
  ndraws <- gssettings$ndraws
  nalldraws <- gssettings$burnin + ndraws * gssettings$saveevery
  
  beta_draws <- array(NA, dim = c(K, N, ndraws))
  sigma_draws <- array(NA, dim = c(N, N, ndraws))
  
  # Starting value for B
  BB <- matrix(0, nrow = K, ncol = N)
  BB[, (Nm + 1):N] <- prior_B
  
  it_save <- 0
  
  cat("Starting Gibbs Sampler...\n")
  # X'X does not change across Gibbs iterations: compute it once
  XtX <- crossprod(X_t)

  for (draw in 1:nalldraws) {
    
    # Draw Sigma
    U <- Y_t - X_t %*% BB
    S_post <- crossprod(U) + S_prior
    S_post_chol <- t(chol(S_post))
    Z <- matrix(rnorm(v_post * N), nrow = v_post, ncol = N)
    Sigma <- S_post_chol %*% solve(t(Z) %*% Z) %*% t(S_post_chol)
    
    # Draw B
    C_sig <- t(chol(Sigma))
    Sig_yy_inv <- t(solve(C_sig[(Nm + 1):N, (Nm + 1):N])) %*% solve(C_sig[(Nm + 1):N, (Nm + 1):N])
    
    A_var <- Q_inv + kronecker(Sig_yy_inv, XtX)
    
    y_st <- Y_t[, (Nm + 1):N] - Y_t[, 1:Nm] %*% solve(t(C_sig[1:Nm, 1:Nm])) %*% t(C_sig[(Nm + 1):N, 1:Nm])
    a_vec <- prior_QinvB + as.vector(crossprod(X_t, y_st) %*% Sig_yy_inv)
    
    C_var <- chol(A_var)
    B_vec <- solve(C_var, solve(t(C_var), a_vec) + rnorm(K * Ny))
    BB[, (Nm + 1):N] <- matrix(B_vec, nrow = K, ncol = Ny)
    
    # Save draws
    if (draw > gssettings$burnin && (draw %% gssettings$saveevery == 0)) {
      it_save <- it_save + 1
      
      BB_reordered <- rbind(BB[-1, ], BB[1, ])
      
      beta_draws[, , it_save] <- BB_reordered
      sigma_draws[, , it_save] <- Sigma
    }
  }
  
  return(list(beta_draws = beta_draws, sigma_draws = sigma_draws, lags = P))
}

#================================
# Sign Restrictions and Rotation
#================================
resirfssign <- function(res, nstep, dims, test_restr, b_normalize, max_attempts) {
  dims_var <- dim(res$beta_draws)
  K <- dims_var[1]
  N <- dims_var[2]
  ndraws <- dims_var[3]
  P <- res$lags
  
  irfs_draws <- array(NA, dim = c(N, N, nstep, ndraws))
  failures <- 0
  
  cat("Computing sign restrictions...\n")
  for (i in 1:ndraws) {
    betadraw <- res$beta_draws[1:(N * P), , i]
    sigmadraw <- res$sigma_draws[, , i]
    
    # Cholesky IRF
    smat <- chol(sigmadraw)
    
    # Convert B into the 3D array
    B_array <- array(t(betadraw), dim = c(N, N, P))
    irfchol <- impulsdtrf(B_array, smat, nstep)
    
    success <- FALSE
    for (a in 1:max_attempts) {
      Q <- diag(N)
      
      # Random orthogonal matrix via QR decomposition
      rand_mat <- matrix(rnorm(length(dims)^2), length(dims), length(dims))
      qr_dec <- qr(rand_mat)
      QQ <- qr.Q(qr_dec)
      Q[dims, dims] <- QQ
      
      # Normalize shocks to have the desired sign on impact
      impact <- irfchol[, , 1] %*% Q
      toflip <- which(diag(impact) * b_normalize < 0)
      if (length(toflip) > 0) {
        Q[, toflip] <- -Q[, toflip]
      }
      
      # The sign restrictions only involve impact responses: check them on
      # the impact matrix first and rotate the whole IRF array only for an
      # accepted candidate (same draws and results, far fewer matrix products)
      impactcand <- array(irfchol[, , 1] %*% Q, dim = c(N, N, 1))
      
      if (test_restr(impactcand)) {
        respcand <- array(0, dim = c(N, N, nstep))
        for(step in 1:nstep) {
          respcand[, , step] <- irfchol[, , step] %*% Q
        }
        irfs_draws[, , , i] <- respcand
        success <- TRUE
        break
      }
    }
    if (!success) {
      failures <- failures + 1
    }
  }
  
  cat(sprintf("Total failures finding valid rotation: %d\n", failures))
  return(irfs_draws)
}


#==============================================================================
# Helper: Extract Quantiles From MATLAB-Style 4D Array
#==============================================================================
extract_matlab_irfs <- function(irfs_array, var_names, shock_idx = 1) {
  
  # Subset to the requested shock -> Shape becomes [Variables, Horizon, Draws]
  shock_data <- irfs_array[, shock_idx, , ]
  
  # Calculate quantiles across the 'Draws' dimension
  quantiles <- apply(shock_data, MARGIN = c(1, 2), FUN = function(x) {
    c(
      Lower68 = quantile(x, 0.16, na.rm = TRUE),
      Median  = median(x, na.rm = TRUE),
      Upper68 = quantile(x, 0.84, na.rm = TRUE)
    )
  })
  # quantiles shape is now [3 (stats), Variables, Horizon]
  
  # Flatten into a ggplot-friendly data frame
  N <- dim(quantiles)[2]
  H <- dim(quantiles)[3]
  
  df_list <- list()
  for(i in 1:N) {
    df_list[[i]] <- data.frame(
      Horizon  = 0:(H - 1),
      Variable = var_names[i],
      Lower68  = quantiles[1, i, ],
      Median   = quantiles[2, i, ],
      Upper68  = quantiles[3, i, ]
    )
  }
  
  return(do.call(rbind, df_list))
}