
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("TS_functions.R")
source("joint grid search.R")
source("JK_BVAR_functions.R")

var_data <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, 
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M,
         PCCI, ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, VIX)%>%
  drop_na()
var_data1 <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, 
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M,
         PCCI, ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, VIX,
         EU_ILS_5Y5Y, EU_ILS_1Y)%>%
  drop_na()

Log_GPR <- var_data$Log_GPR
VIX <- var_data$VIX

endog_vars1 <- c("OIS_3M", "STOXX50", "True_Policy_Rate","Log_IP", "Log_HICP", "Log_Oil" )
endog_vars2 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "PCCI", "True_Policy_Rate" )
endog_vars3 <- c("OIS_3M", "STOXX50","df_Unemployment_Rate","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars4 <- c("OIS_3M", "STOXX50","ESI","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars5 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "Log_HICP", "EA_2Y_Yield" )
endog_vars6 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "Log_HICP", "Spliced_EONIA_3M" )
endog_vars7 <- c("OIS_3M", "STOXX50","Log_IP","EU_ILS_5Y5Y", "Log_HICP", "True_Policy_Rate")

Y_matrix <- as.matrix(var_data[, endog_vars1])

#---------------------------------------------
# Parameter Estimation with grid Search
#---------------------------------------------

joint_opt <- optimize_all_st_bvar_hyperparameters(
  Y            = Y_matrix,
  GPR          = Log_GPR,
  p_grid       = c( 5, 6, 7, 8),                                    
  d_grid       = c(1, 2),                                    # Try 1 and 2 month delay
  gamma_pcts   = c(0.3, 0.4, 0.50),                        # Percentile thresholds
  c_grid       = c(7.0 , 10.0),                          # Transition speeds
  lambda1_grid = c(0.05, 0.08, 0.10),            # Shrinkage tightness
  lambda3_grid = c(1.5, 2.0, 3.0),                           # Lag decay rates
  prior_means  = c(0, 0, 0.9, 0.9, 0.9, 0.9)
)

# View the Top 5 specifications ranked by pooled log marginal likelihood
print(joint_opt$top_10[1:5, c("p", "d", "Gamma_Pct", "c", "Lambda1", "Lambda3", "Pooled_Log_MDD", "PMP")])
opt_p <-joint_opt$optimal$p
opt_c <- joint_opt$optimal$c
opt_d <- joint_opt$optimal$d
opt_gamma <- joint_opt$optimal$Gamma_Val
opt_lambda1 <- joint_opt$optimal$Lambda1
opt_lambda3 <-joint_opt$optimal$Lambda3



# ----------------------------------------------------------------------------------------
# Estimate the Smooth Transititon JK (2020) BVAR 
# ----------------------------------------------------------------------------------------

results <- estimate_st_jk2020_bvar(
  Y           = Y_matrix,
  GPR         = Log_GPR,
  p           = opt_p,              
  d           = opt_d,              
  gamma       = opt_gamma,     
  c           = opt_c,            
  n_draws     = 1000,           # Keep this at 1000 for final smooth bands
  horizon     = 24,            
  lambda1     = 0.25,           # Override opt_lambda1 (Standard is 0.15 - 0.25)
  lambda3     = 1.5,            # Override opt_lambda3 
  lambda4     = 10.0,           # Relax Sum-of-Coefficients (SoC)
  lambda5     = 10.0,           # Relax Dummy Initial Observation (DIO)
  
  prior_means = c(0, 0, 0.9, 0.90, 0.90, 0.90)
)

results$IRF_Low_GPR <- normalize_irf_array(
  irf_array      = results$IRF_Low_GPR, 
  target_var_idx = 3, 
  target_size    = 0.25
)

results$IRF_High_GPR <- normalize_irf_array(
  irf_array      = results$IRF_High_GPR, 
  target_var_idx = 3, 
  target_size    = 0.25
)

extract_irf_quantiles <- function(irf_array) {
  # irf_array dimensions: [draws, horizon, variable, shock]
  # We extract Shock 1 (Pure Monetary Policy Shock)
  
  apply(irf_array[, , , 1], MARGIN = c(2, 3), function(x) {
    c(
      Lower68 = quantile(x, 0.16, na.rm = TRUE),
      Median  = median(x, na.rm = TRUE),
      Upper68 = quantile(x, 0.84, na.rm = TRUE)
    )
  })
}

quantiles_low  <- extract_irf_quantiles(results$IRF_Low_GPR)
quantiles_high <- extract_irf_quantiles(results$IRF_High_GPR)
var_names <- c("MP_shocks", "stocks_shocks", "i_ECB", "IP", "pi", "Oil")

# 2. Helper function to flatten the array
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

# 3. Build and combine the data frames
df_low  <- build_plot_df(quantiles_low, "Low GPR (Peacetime)", var_names)
df_high <- build_plot_df(quantiles_high, "High GPR (Crisis)", var_names)

plot_data <- rbind(df_low, df_high)

plot_data <- plot_data %>%
  # Keep only the macroeconomic variables
  filter(Variable %in% c("i_ECB", "IP", "pi", "Oil")) %>%
  # Lock the factor levels in the exact order you want them displayed
  mutate(Variable = factor(Variable, levels = c("i_ECB", "IP", "pi", "Oil")))


# ==============================================================================
# RENDER PLOT
# ==============================================================================
p_irf <- ggplot(plot_data, aes(x = Horizon, y = Median, color = Regime, fill = Regime)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  
  # Set ncol = 2 or ncol = 4 depending on how you want the 4 graphs arranged
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  
  scale_color_manual(values = c("Low GPR (Peacetime)" = "#1f78b4", "High GPR (Crisis)" = "#e31a1c")) +
  scale_fill_manual( values = c("Low GPR (Peacetime)" = "#1f78b4", "High GPR (Crisis)" = "#e31a1c")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Impulse Responses to a +25bp Pure Monetary Policy Shock",
    subtitle = "Comparing Low vs. High Geopolitical Risk Regimes (68% Credible Intervals)",
    x = "Months after Shock",
    y = "Response"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_irf)