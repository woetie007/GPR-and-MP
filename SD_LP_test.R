
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("TS_functions.R")
source("joint grid search.R")
source("SD_LP_functions.R")

var_data <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, Log_GPRA, Log_GPRT,
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M
         , ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX)%>%
  drop_na()
var_data1 <- final_data %>% filter(Country == "EA")%>% 
  select(df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, 
         Log_HICP,Log_Core_HICP, Log_Oil, EA_2Y_Yield, Spliced_EONIA_3M,
         PCCI, ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX,
         EU_ILS_5Y5Y, EU_ILS_1Y)%>%
  drop_na()

Threshold_var <- var_data$Log_GPR

endog_vars1 <- c("OIS_3M", "STOXX50", "True_Policy_Rate","Log_IP", "Log_HICP", "Log_Oil" )
endog_vars2 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "PCCI", "True_Policy_Rate" )
endog_vars3 <- c("OIS_3M", "STOXX50","df_Unemployment_Rate","Log_Oil", "Log_HICP", "True_Policy_Rate" )
endog_vars4 <- c("OIS_3M", "STOXX50","Spliced_EONIA_3M", "ESI", "Log_HICP","Log_Oil" )
endog_vars5 <- c("OIS_3M", "STOXX50","Log_IP","Log_Oil", "Log_HICP", "EA_2Y_Yield" )
endog_vars6 <- c("OIS_3M", "STOXX50","Spliced_EONIA_3M","Log_IP", "Log_HICP","Log_Oil")
endog_vars7 <- c("OIS_3M", "STOXX50","Log_IP","EU_ILS_5Y5Y", "Log_HICP", "True_Policy_Rate")

Y_matrix <- as.matrix(var_data[, endog_vars1])

#===========================================
# Run the State-Dependent Local-Projection
#===========================================
cat("Estimating SD-LP Models...\n")

lp_results <- estimate_sd_lp_jk2020(
  Y          = Y_matrix,
  GPR        = Threshold_var,
  p          = 4,               
  d          = opt_d,           # Use optimal delay from your grid search
  gamma      = opt_gamma,       # Use optimal threshold from your grid search
  c          = opt_c,           # Use optimal transition speed
  horizon    = 24,
  shock_type = "Pure MP"
)

#--------------------------------------
# Scale to +25bps to compare with BVAR
#--------------------------------------
# Find the exact Month 0 (impact) response of True_Policy_Rate in both regimes
impact_low <- lp_results %>% 
  filter(Variable == "OIS_3M", Horizon == 0, Regime == "Low GPR Regime") %>% 
  pull(Beta)

impact_high <- lp_results %>% 
  filter(Variable == "OIS_3M", Horizon == 0, Regime == "High GPR Regime") %>% 
  pull(Beta)

# Print them out to verify they are healthy, non-zero numbers!
cat(sprintf("OIS_3M Impact - Low Regime: %.4f\n", impact_low))
cat(sprintf("OIS_3M Impact - High Regime: %.4f\n", impact_high))

# Multiply all estimates and bands by (0.25 / OIS_impact)
lp_scaled <- lp_results %>%
  mutate(
    scale_factor = ifelse(Regime == "Low GPR Regime", 0.25 / impact_low, 0.25 / impact_high),
    Beta    = Beta * scale_factor,
    Lower68 = Lower68 * scale_factor,
    Upper68 = Upper68 * scale_factor,
    Lower90 = Lower90 * scale_factor,
    Upper90 = Upper90 * scale_factor
  )

#---------------------------
# Filter Variables & Plot
#---------------------------
# Keep only True_Policy_Rate, Log_IP, and Log_HICP
lp_plot_data <- lp_scaled %>%
  filter(Variable %in% c("True_Policy_Rate", "Log_IP", "Log_HICP", "Log_Oil")) %>%
  mutate(Variable = factor(Variable, levels = c("True_Policy_Rate", "Log_IP", "Log_HICP", "Log_Oil")))

p_lp <- ggplot(lp_plot_data, aes(x = Horizon, y = Beta, color = Regime, fill = Regime)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Low GPR Regime" = "#1f78b4", "High GPR Regime" = "#e31a1c")) +
  scale_fill_manual( values = c("Low GPR Regime" = "#1f78b4", "High GPR Regime" = "#e31a1c")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "State-Dependent LP: Impulse Responses to a +25bp MP Shock",
    subtitle = "Comparing Low vs. High GPR Regimes (Scaled to OIS_3M, 68% HAC Bands)",
    x = "Months after Shock",
    y = "Response"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    strip.text       = element_text(face = "bold", size = 14),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_lp)

