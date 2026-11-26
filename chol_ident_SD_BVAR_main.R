
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("jk_replication_function.R")

#======
# Data
#======
var_data <- final_data %>% filter(Country == "EA")%>% 
  select(Date, `Ind. Production`=Log_IP, `Diff. Oil Prices`=df_Oil,`True Policy Rate` =True_Policy_Rate, 
         Log_GPR,Log_GPRA, Log_GPRT,`HICP`=Log_HICP,`Core HICP`=Log_Core_HICP, 
         `Oil Prices` = Log_Oil, `EONIA/€STR 3M`=Spliced_EONIA_3M, , ESI, OIS_3M, 
         `OIS 1Y surprise`=OIS_1Y, STOXX50,  Log_VIX,`EONIA/€STR 1Y`=Spliced_EONIA_1Y,
         `EPU`= Log_EPU, `Diff. Unemployment Rate`= df_Unemployment_Rate)%>% 
  drop_na()


#==================
# Settings & Setup
#==================
ynames1 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP")

#===========
# Variables
#===========
mnames <- c("OIS_3M", "STOXX50")
ynames <- ynames1
data_names <- c(mnames, ynames)
N <- length(data_names)

#================
# Prior Settings
#================
prior <- list(
  lags = 12,
  minnesota_tightness = 0.2,
  minnesota_decay = 1,
  Nm = length(mnames)
)

#========================
# Gibbs Sampler Settings
#========================
gssettings <- list(
  ndraws = 4000,
  burnin = 4000,
  saveevery = 4,
  computemarglik = 0
)

#----------------
# Prepare matrix
#----------------
Y_matrix <- as.matrix(var_data[, data_names])
prior$minnesota_mvector <- c(0, 0, 0.9, 0.9, 0.9, 0.9, 0.9)

#--------------------------------------------------
# Calculate Smooth Transition Weights from Log_GPR
#--------------------------------------------------
d_delay <- 1
c_speed <- 10.0
gamma_threshold <- median(var_data$Log_GPR, na.rm = TRUE)

z_t <- var_data$Log_GPR[(prior$lags + 1 - d_delay):(nrow(var_data) - d_delay)]

F_high <- 1 / (1 + exp(-c_speed * (z_t - gamma_threshold) / sd(z_t)))
F_low  <- 1 - F_high

#==============================================================================
# Estimate the State-Dependent VARs
#==============================================================================
cat("\nEstimating Low GPR Regime...\n")
res_low  <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = F_low)

cat("\nEstimating High GPR Regime...\n")
res_high <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = F_high)


# ==============================================================================
# State-Dependent Cholesky Identification
# ==============================================================================

MAlags <- 36
ndraws <- dim(res_high$beta_draws)[3]

irfs_high_draws <- array(NA, dim = c(N, N, MAlags, ndraws))
irfs_low_draws  <- array(NA, dim = c(N, N, MAlags, ndraws))

cat("Computing State-Dependent Cholesky IRFs...\n")
for (i in 1:ndraws) {
  
  betadraw_h <- res_high$beta_draws[1:(N * P), , i]
  sigmadraw_h <- res_high$sigma_draws[, , i]
  
  smat_h <- chol(sigmadraw_h)
  B_array_h <- array(t(betadraw_h), dim = c(N, N, P))
  
  irfs_high_draws[, , , i] <- impulsdtrf(B_array_h, smat_h, MAlags)
  
  betadraw_l <- res_low$beta_draws[1:(N * P), , i]
  sigmadraw_l <- res_low$sigma_draws[, , i]
  
  smat_l <- chol(sigmadraw_l)
  B_array_l <- array(t(betadraw_l), dim = c(N, N, P))
  
  irfs_low_draws[, , , i] <- impulsdtrf(B_array_l, smat_l, MAlags)
}
cat("Done.\n")

# ==============================================================================
# Extract Quantiles for Plotting
# ==============================================================================
df_low  <- extract_matlab_irfs(irfs_low_draws, data_names, shock_idx = 1)
df_high <- extract_matlab_irfs(irfs_high_draws, data_names, shock_idx = 1)

df_low$Regime  <- "Low GPR"
df_high$Regime <- "High GPR"
plot_data <- rbind(df_low, df_high)

plot_data <- plot_data %>%
  filter(Variable %in% ynames) %>%
  mutate(Variable = factor(Variable, levels = ynames))

plot_data$Horizon <- as.numeric(plot_data$Horizon)

T_eff <- nrow(Y_matrix) - prior$lags

p <- ggplot(plot_data, aes(x = Horizon, color = Regime, fill = Regime)) +
  
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
  
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  
  geom_line(aes(y = Median), linewidth = 1) +
  
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  
  scale_x_continuous(expand = c(0, 0)) +
  
  scale_color_manual(values = c("High GPR" = "firebrick", "Low GPR" = "steelblue")) +
  scale_fill_manual(values = c("High GPR" = "firebrick", "Low GPR" = "steelblue")) +
  
  labs(
    title = "State-Dependent Responses to a Monetary Policy Shock",
    subtitle = "HFI Choleski Identification: High vs. Low Geopolitical Risk Regimes",
    caption = sprintf("Effective Observations (T): %d | Lags (p): %d", T_eff, prior$lags),
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

print(p)
