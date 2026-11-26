
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("TS_functions.R")
source("TVAR_functions.R")
source("joint grid search.R")
source("JK_BVAR_functions.R")

var_data <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, Log_GPRA, Log_GPRT,
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M
         , ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX,
         jk_OIS_3M_monthly, jk_STOXX50_monthly, Spliced_EONIA_1Y,Spliced_EONIA_2Y )%>%
  drop_na()
var_data1 <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, 
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M,
         PCCI, ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, VIX,
         EU_ILS_5Y5Y, EU_ILS_1Y)%>%
  drop_na()

Threshold_var <- var_data$Log_GPR

endog_vars1 <- c("jk_OIS_3M_monthly", "jk_STOXX50_monthly", "Log_IP","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars2 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "PCCI", "True_Policy_Rate" )
endog_vars3 <- c("OIS_3M", "STOXX50","df_Unemployment_Rate","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars4 <- c("OIS_3M", "STOXX50","ESI","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars5 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "Log_HICP", "EA_2Y_Yield" )
endog_vars6 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "Log_HICP", "Spliced_EONIA_1Y" )
endog_vars7 <- c("OIS_3M", "STOXX50","Log_IP","EU_ILS_5Y5Y", "Log_HICP", "True_Policy_Rate")

Y_matrix <- as.matrix(var_data[, endog_vars1])

#---------------------------------------------
# Parameter Estimation with grid Search
#---------------------------------------------

joint_opt <- optimize_all_st_bvar_hyperparameters(
  Y            = Y_matrix,
  GPR          = Threshold_var,
  p_grid       = c(10,11,12),                                    # Try 4 and 6 lags
  d_grid       = c(1, 2),                                    # Try 1 and 2 month delay
  gamma_pcts   = c(0.4, 0.5),                        # Percentile thresholds
  c_grid       = c(7.0 , 10.0),                          # Transition speeds
  lambda1_grid = c(0.25),            # Shrinkage tightness
  lambda3_grid = c(1.5),                           # Lag decay rates
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
# 1.IR Shock, 2.Stocks Shock, 3.Economic activity,4.Oil prices, 5.Inflation, 6.Policy Rate
# ----------------------------------------------------------------------------------------

results <- estimate_st_jk2020_bvar(
  Y       = Y_matrix,
  GPR     = Threshold_var,
  p       = opt_p,              # VAR lag order
  d       = opt_d,              # GPR transition lag delay
  gamma   = opt_gamma,     # Threshold location
  c       = opt_c,            # Transition speed
  n_draws = 100,           # Number of valid posterior structural draws
  horizon = 24,            # IRF horizon (24 months)
  lambda1 = opt_lambda1,
  lambda3 = opt_lambda3,
  prior_means = c(0, 0, 0.9, 0.90, 0.90, 0.90)
)

sum(results$Transition_Weights$F_high)
sum(results$Transition_Weights$F_low)
results$IRF_Low_GPR <- normalize_irf_array(
  irf_array      = results$IRF_Low_GPR, 
  target_var_idx = 1,     
  target_size    = 0.25
)

results$IRF_High_GPR <- normalize_irf_array(
  irf_array      = results$IRF_High_GPR, 
  target_var_idx = 1,     
  target_size    = 0.25
)


# -------------------------------
# IRFs and Confidence Bands
# -------------------------------
# Helper function to summarize an IRF array into a clean ggplot dataframe
summarize_irf_draws <- function(irf_array, regime_name, var_names) {
  # irf_array dimensions: [draws, horizon+1, response_var, shock_var]
  H <- dim(irf_array)[2] - 1
  N <- dim(irf_array)[3]
  
  df_list <- list()
  idx <- 1
  
  for (shock_idx in 1:2) { # 1 = Pure MP Shock, 2 = CBI Shock
    for (resp_idx in 1:N) {
      for (h in 0:H) {
        draws_h <- irf_array[, h + 1, resp_idx, shock_idx]
        
        df_list[[idx]] <- data.frame(
          Horizon  = h,
          Regime   = regime_name,
          Shock    = ifelse(shock_idx == 1, "Pure MP Shock", "CBI Shock"),
          Response = var_names[resp_idx],
          Median   = median(draws_h),
          Lower68  = quantile(draws_h, 0.16),
          Upper68  = quantile(draws_h, 0.84),
          Lower90  = quantile(draws_h, 0.05),
          Upper90  = quantile(draws_h, 0.95)
        )
        idx <- idx + 1
      }
    }
  }
  return(do.call(rbind, df_list))
}

var_labels <- c("MP Surprises", "Stock Surprises", 
                "Economic Activity", "Oil Prices", "Inflation", "Policy Rate")

df_low  <- summarize_irf_draws(results$IRF_Low_GPR,  "Low GPR Regime",  var_labels)
df_high <- summarize_irf_draws(results$IRF_High_GPR, "High GPR Regime", var_labels)

irf_all <- rbind(df_low, df_high)

# ------------------------------------------------------------------------------
# Plot & Compare Regimes for Macro Variables
# ------------------------------------------------------------------------------
# Filter to inspect Industrial Production & Inflation responses to a Pure MP Shock
plot_data_pure <- irf_all %>%
  filter(Shock == "Pure MP Shock",
         Response %in% c("Economic Activity", "Oil Prices", "Inflation", "Policy Rate"))

p_comparison_pure <- ggplot(plot_data_pure, aes(x = Horizon, y = Median, color = Regime, fill = Regime)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  # Add 68% credible bands
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Response, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Low GPR Regime" = "#1f78b4", "High GPR Regime" = "#e31a1c")) +
  scale_fill_manual(values  = c("Low GPR Regime" = "#1f78b4", "High GPR Regime" = "#e31a1c")) +
  labs(
    title    = "Macroeconomic Transmission of a Pure ECB Monetary Policy Shock(Base Model)",
    subtitle = "Jarocinski & Karadi (2020) Sign Restrictions across Smooth Geopolitical Risk Regimes",
    x        = "Horizon (Months)",
    y        = "Response"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    strip.text      = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_comparison_pure)




#############################################
# Complete state dependent JK(2020) framework
#############################################
results_st_bvar <- estimate_st_jk2020_gibbs(
  Y           = Y_matrix,
  GPR         = Threshold_var,
  p           = opt_p,              # Keep your optimal lags from the grid search
  d           = opt_d,              # Keep your optimal delay
  gamma       = opt_gamma,          # Keep your optimal threshold
  c           = opt_c,              # Keep your optimal transition speed
  n_draws     = 500,               # 1000 valid structural draws for smooth bands
  burn_in     = 500,                # Let the Markov Chain converge first
  horizon     = 24,            
  
  # The Litterman (1986) / JK (2020) Hyperparameters:
  lambda1     = 0.20,               # Overall tightness
  lambda2     = 0.30,               # Cross-variable tightness (Only possible with Gibbs!)
  lambda3     = 2.0,                # Lag decay
  
  prior_means = c(0, 0, 1, 1, 1, 1) 
)

quantiles_low  <- extract_irf_quantiles(results_st_bvar$IRF_Low_GPR)
quantiles_high <- extract_irf_quantiles(results_st_bvar$IRF_High_GPR)

var_names <- c("OIS_3M", "STOXX50", "Log_IP", "Log_Oil", "Log_HICP", "True_Policy_Rate")

df_low  <- build_plot_df(quantiles_low, "Low GPR (Peacetime)", var_names)
df_high <- build_plot_df(quantiles_high, "High GPR (Crisis)", var_names)

plot_data <- rbind(df_low, df_high)

# ==============================================================================
# 3. FILTER VARIABLES & PLOT
# ==============================================================================
plot_data <- plot_data %>%
  filter(Variable %in% c("Log_IP", "Log_Oil", "Log_HICP", "True_Policy_Rate")) %>%
  mutate(Variable = factor(Variable, levels = c("True_Policy_Rate", "Log_IP", "Log_HICP", "Log_Oil")))

p_1sd <- ggplot(plot_data, aes(x = Horizon, y = Median, color = Regime, fill = Regime)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  
  # free_y is critical here because a 1-SD shock will have completely different 
  # numeric scales for interest rates vs. log variables
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  
  scale_color_manual(values = c("Low GPR (Peacetime)" = "#1f78b4", "High GPR (Crisis)" = "#e31a1c")) +
  scale_fill_manual( values = c("Low GPR (Peacetime)" = "#1f78b4", "High GPR (Crisis)" = "#e31a1c")) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Impulse Responses to a 1-Standard-Deviation MP Shock",
    subtitle = "Comparing Low vs. High GPR Regimes (68% Credible Intervals)",
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

print(p_1sd)

#==================================
# JK 2020 without state dependence
#==================================
var_data <- final_data %>% filter(Country == "EA")%>% 
  filter(Date >= as.Date("2000-01-01")) %>%
  filter(Date <= as.Date("2016-12-01")) %>%
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, Log_GPRA, Log_GPRT,
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M
         , ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX,
         jk_OIS_3M_monthly, jk_STOXX50_monthly, Spliced_EONIA_1Y,Spliced_EONIA_2Y )%>%
  drop_na()
var_data1 <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, 
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M,
         PCCI, ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, VIX,
         EU_ILS_5Y5Y, EU_ILS_1Y)%>%
  drop_na()

Threshold_var <- var_data$Log_GPR

endog_vars1 <- c("jk_OIS_3M_monthly", "jk_STOXX50_monthly", "Log_IP","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars2 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "PCCI", "True_Policy_Rate" )
endog_vars3 <- c("OIS_3M", "STOXX50","df_Unemployment_Rate","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars4 <- c("OIS_3M", "STOXX50","ESI","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars5 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "Log_HICP", "EA_2Y_Yield" )
endog_vars6 <- c("jk_OIS_3M_monthly", "jk_STOXX50_monthly","Log_IP","Log_Oil", "Log_HICP", "Spliced_EONIA_1Y" )
endog_vars7 <- c("OIS_3M", "STOXX50","Log_IP","EU_ILS_5Y5Y", "Log_HICP", "True_Policy_Rate")

Y_matrix <- as.matrix(var_data[, endog_vars6])
# -------------------------------
# 1. ESTIMATE THE BASELINE MODEL
# -------------------------------
# Note: Since this is linear, we drop the GPR arguments entirely.
results_linear <- estimate_jk2020_bvar(
  Y           = Y_matrix,
  p           = 12,               
  n_draws     = 500,               
  burn_in     = 500,                
  horizon     = 36,            
  
  lambda1     = 0.20,               
  lambda2     = 0.30,               
  lambda3     = 2.0,                
  
  prior_means = c(0, 0, 1, 1, 1, 1) 
)

# ==============================================================================
# 2. EXTRACT & FORMAT
# ==============================================================================
# We skip the normalizer so we map out the pure 1-standard-deviation structural shocks.
quantiles_linear <- extract_irf_quantiles(results_linear)
var_names <- c("OIS_3M", "STOXX50", "Log_IP", "Log_Oil", "Log_HICP", "True_Policy_Rate")

# We can reuse your existing build_plot_df by passing "Baseline" as the regime name
df_linear <- build_plot_df(quantiles_linear, "Linear Baseline", var_names)

plot_data_linear <- df_linear %>%
  filter(Variable %in% c("Log_IP", "Log_Oil", "Log_HICP", "True_Policy_Rate")) %>%
  mutate(Variable = factor(Variable, levels = c("True_Policy_Rate", "Log_IP", "Log_HICP", "Log_Oil")))

# ==============================================================================
# 3. PLOT
# ==============================================================================
p_linear <- ggplot(plot_data_linear, aes(x = Horizon, y = Median)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  
  # Standard 68% credible bands
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), fill = "#2c3e50", alpha = 0.2) +
  geom_line(color = "#2c3e50", linewidth = 1.2) +
  
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Original JK (2020) Baseline BVAR",
    subtitle = "Impulse Responses to a 1-Standard-Deviation Pure Monetary Policy Shock",
    x = "Months after Shock",
    y = "Response"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    strip.text       = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

print(p_linear)


