

ARMA_IC_selection_CSSML <- function(Pmax, Qmax, series) {
  # Validate inputs
  if (!is.numeric(Pmax) || !is.numeric(Qmax) || Pmax < 0 || Qmax < 0) {
    stop("Pmax and Qmax must be non-negative numeric values.")
  }
  
  total_models <- (Pmax + 1) * (Qmax + 1)  # Total number of models
  aic_bic_results <- matrix(NA, ncol = 7, nrow = total_models)
  colnames(aic_bic_results) <- c("p", "q", "AIC", "BIC", "RSS", "SampleSize", "Method")
  
  k <- 1
  for (p in 0:Pmax) {
    for (q in 0:Qmax) {
      # Fit the model with error handling
      fit <- tryCatch(
        arima(series, order = c(p, 0, q), method = "CSS-ML", include.mean = TRUE, transform.pars = FALSE),
        error = function(e) NULL # Return NULL if fitting fails
      )
      
      if (!is.null(fit)) {  # Only process if the model fit succeeds
        # Calculate RSS (Residual Sum of Squares)
        rss <- sum(fit$residuals^2, na.rm = TRUE)
        # Get the sample size
        sample_size <- length(fit$residuals)
        # Store results
        aic_bic_results[k, ] <- c(p, q, AIC(fit, k = 2), AIC(fit, k = log(sample_size)), rss, sample_size, "CSS-ML")
      } else {
        # Assign NA if model fitting failed
        aic_bic_results[k, ] <- c(p, q, NA, NA, NA, NA, "CSS-ML")
      }
      k <- k + 1
    }
  }
  
  # Convert results to a data frame
  aic_bic_results <- as.data.frame(aic_bic_results)
  # Convert numeric columns back to numeric type (R matrices coerce everything to character)
  numeric_cols <- c("p", "q", "AIC", "BIC", "RSS", "SampleSize")
  aic_bic_results[numeric_cols] <- lapply(aic_bic_results[numeric_cols], as.numeric)
  # Order according to BIC and reset the index after sorting
  aic_bic_results <- aic_bic_results[order(aic_bic_results$BIC), ]
  row.names(aic_bic_results) <- NULL
  # Return the results table
  return(aic_bic_results)
}

ARMA_IC_selection_CSS <- function(Pmax, Qmax, series) {
  # Validate inputs
  if (!is.numeric(Pmax) || !is.numeric(Qmax) || Pmax < 0 || Qmax < 0) {
    stop("Pmax and Qmax must be non-negative numeric values.")
  }
  if (!is.numeric(series) && !is.ts(series)) {
    stop("series must be a numeric vector or a time series object.")
  }
  
  total_models <- (Pmax + 1) * (Qmax + 1)  # Total number of models
  aic_bic_results <- matrix(NA, ncol = 7, nrow = total_models)
  colnames(aic_bic_results) <- c("p", "q", "AIC", "BIC", "RSS", "SampleSize", "Method")
  
  k <- 1
  for (p in 0:Pmax) {
    for (q in 0:Qmax) {
      # Fit the model with error handling
      fit <- tryCatch(
        arima(series, order = c(p, 0, q), method = "CSS", n.cond = Pmax, include.mean = TRUE, transform.pars = FALSE),
        error = function(e) NULL # Return NULL if fitting fails
      )
      
      if (!is.null(fit)) {  # Only process if the model fit succeeds
        # Calculate RSS (Residual Sum of Squares)
        rss <- sum(residuals(fit)^2, na.rm = TRUE)
        # Get the sample size
        sample_size <- length(fit$residuals) - Pmax
        # Calculate AIC and BIC manually (CSS does not calculate them directly)
        aic <- (fit$nobs - Pmax) * log(rss) + 2 * length(fit$coef)
        bic <- (fit$nobs - Pmax) * log(rss) + length(fit$coef) * log(fit$nobs - Pmax)
        aic_bic_results[k, ] <- c(p, q, aic, bic, rss, sample_size, "CSS")
      } else {
        # Assign NA if model fitting failed
        aic_bic_results[k, ] <- c(p, q, NA, NA, NA, NA, "CSS")
      }
      k <- k + 1
    }
  }
  
  # Convert results to a data frame
  aic_bic_results <- as.data.frame(aic_bic_results)
  # Convert numeric columns back to numeric type (R matrices coerce everything to character)
  numeric_cols <- c("p", "q", "AIC", "BIC", "RSS", "SampleSize")
  aic_bic_results[numeric_cols] <- lapply(aic_bic_results[numeric_cols], as.numeric)
  
  # Order according to BIC and reset the index after sorting
  aic_bic_results <- aic_bic_results[order(aic_bic_results$BIC), ]
  row.names(aic_bic_results) <- NULL
  
  # Return the results table
  return(aic_bic_results)
}




find_optimal_lags_for_ADF <- function(timeseries, max_lags) {
  # Ensure the `timeseries` is a time series object
  if (!is.ts(timeseries)) {
    stop("Input 'timeseries' must be a time series object.")
  }
  
  # Initialize a results data frame
  results <- data.frame(
    Lags = integer(),
    AIC = numeric(),
    BIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over the number of lags
  for (lags in 1:max_lags) {
    # Fit the AR model with the current number of lags
    model <- tryCatch({
      arima(timeseries, order = c(lags, 0, 0), method = "CSS-ML", include.mean = TRUE) # Force ML for better handling of non-stationarity
    }, error = function(e) {
      return(NULL) # Handle cases where model fitting fails
    })
    
    if (!is.null(model)) {
      # Extract AIC and BIC
      aic <- AIC(model)
      bic <- BIC(model)
      
      # Append the results to the data frame
      results <- rbind(results, data.frame(Lags = lags, AIC = aic, BIC = bic))
    }
  }
  
  if (nrow(results) == 0) {
    stop("No valid models were fitted. Check your time series.")
  }
  
  # Sort the results first by BIC, then by AIC
  results <- results[order(results$BIC, results$AIC), ]
  
  # Get the optimal number of lags for AIC and BIC
  optimal_aic_lags <- results$Lags[which.min(results$AIC)]
  optimal_bic_lags <- results$Lags[which.min(results$BIC)]
  
  # Display the sorted results table
  print("Sorted Results:")
  print(results)
  
  # Return the optimal lags
  return(list(
    Optimal_AIC_Lags = optimal_aic_lags,
    Optimal_BIC_Lags = optimal_bic_lags
  ))
}
