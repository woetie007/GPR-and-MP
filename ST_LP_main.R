
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("ST_LP_functions.R") # Functions for the Smooth Transition Local Projection
source("jk_replication_main.R") # You need to run this script to get the data and the structural MP shocks from the BVAR


mp_shock <- ts(var_data_eff$`Pure MP Shock`, start = c(2001, 2), frequency = 12)
LB_test <- LjungBox(mp_shock, 12, 0)
stargazer(LB_test, type = "text", summary = FALSE)
acf(mp_shock)
ARMA_IC_selection_CSS(15,0,mp_shock)
AR_MP <- arima(mp_shock, order = c(1, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_residuals_MP <- residuals(AR_MP)
Qstat_AR_residuals <- LjungBox(AR_residuals_MP, 20, 1)
stargazer(Qstat_AR_residuals, type = "text", summary = FALSE)
var_data_eff$`Purged MP Shock` <- AR_residuals_MP


#==================
# Settings & Setup
#==================
ynames1 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP")

#===========
# Variables
#===========
mnames1 <- "Purged MP Shock"
mnames <- mnames1
ynames <- ynames1
threshname <- "Log_GPR"
data_names <- c(mnames, ynames, threshname)
N <- length(data_names)





#----------------
# Prepare matrix
#----------------
Y_matrix <- as.matrix(var_data_eff[, data_names])
p=12
variables_to_run <- ynames
lp_results_list <- lapply(variables_to_run, function(v) {
  run_st_lp(
    data = Y_matrix,
    target_var = v,
    shock_var = mnames,
    gpr_var = threshname,
    p_lags =p,
    horizon = 24
  )
})

lp_plot_data <- do.call(rbind, lp_results_list)

lp_long <- bind_rows(
  lp_plot_data %>% 
    select(Horizon, Variable, Median = Beta_Low, Lower68 = Lower68_Low, Upper68 = Upper68_Low) %>% 
    mutate(Regime = "Low GPR"),
  lp_plot_data %>% 
    select(Horizon, Variable, Median = Beta_High, Lower68 = Lower68_High, Upper68 = Upper68_High) %>% 
    mutate(Regime = "High GPR")
)

T_eff <- nrow(Y_matrix) - p


p_lp <- ggplot(lp_long, aes(x = Horizon, color = Regime, fill = Regime)) +
  
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
  
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  
  geom_line(aes(y = Median), linewidth = 1) +
  
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  
  scale_x_continuous(expand = c(0, 0)) +
  
  scale_color_manual(values = c("High GPR" = "firebrick", 
                                "Low GPR" = "steelblue")) +
  scale_fill_manual(values = c("High GPR" = "firebrick", 
                               "Low GPR" = "steelblue")) +
  
  labs(
    title = "State-Dependent Local Projection: Response to a Pure Monetary Policy Shock",
    subtitle = "High vs. Low Geopolitical Risk Regimes",
    caption = sprintf("Effective Observations (T): %d | Lags (p): %d", T_eff, p),
    x = "Months after shock",
    y = "Impulse Response"
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_lp)