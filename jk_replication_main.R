#############################################################
# REPLICATION OF JAROCINSKI & KARADI (2020) - EURO AREA 
##############################################################
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
         `EPU`= Log_EPU, `Diff. Unemployment Rate`= df_Unemployment_Rate, `Diff. HICP`=df_HICP
  )%>% 
  # Adding PCCI , `EONIA/€STR 5Y` = Spliced_EONIA_5Y and other rates gets rid of obs
  drop_na()
#==================
# Settings & Setup
#==================
idscheme <- "sgnm2" 

#===========
# Variables
#===========
mnames <- c("OIS_3M", "STOXX50")
ynames <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP")
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
prior$minnesota_mvector <- c(0, 0, 0.9, 0.9, 0.9, 0.9, 0.9) # 0 for surprises, 0.9 for levels

#==================
# Estimate the VAR
#==================
res <- VAR_withiid(Y_matrix, prior, gssettings)

#--------------------------------------
# Identification with Sign Restriction
#--------------------------------------
MAlags <- 36

if (idscheme == "sgnm2") {
  
  dims <- c(1, 2)
  imonpol <- 1
  inews <- 2
  
  test_restr <- function(irfs) {
    (irfs[1, imonpol, 1] > 0) && (irfs[2, imonpol, 1] < 0) && 
      (irfs[1, inews, 1] > 0)   && (irfs[2, inews, 1] > 0)
  }
  
  b_normalize <- rep(1, N)
  max_try <- 1000
  
  irfs_draws <- resirfssign(res, MAlags, dims, test_restr, b_normalize, max_try)
  
} else if (idscheme == "sgnm2strg") {
  
  dims <- c(1, 2)
  iyld <- 3  
  imonpol <- 1
  inews <- 2
  
  test_restr <- function(irfs) {
    (irfs[1, imonpol, 1] > 0) && (irfs[2, imonpol, 1] < 0) && (irfs[iyld, imonpol, 1] > 0.01) && 
      (irfs[1, inews, 1] > 0)   && (irfs[2, inews, 1] > 0)
  }
  
  b_normalize <- rep(1, N)
  max_try <- 1000
  
   irfs_draws <- resirfssign(res, MAlags, dims, test_restr, b_normalize, max_try)
}


plot_data <- extract_matlab_irfs(irfs_draws, data_names , shock_idx = 1)

plot_data <- plot_data %>%
  filter(Variable %in% ynames) %>%
  mutate(Variable = factor(Variable, levels = ynames))

T_eff <- nrow(Y_matrix) - prior$lags

p_baseline <- ggplot(plot_data, aes(x = Horizon, y = Median)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), fill = "#2c3e50", alpha = 0.2) +
  geom_line(color = "#2c3e50", linewidth = 1.2) +
  
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Baseline BVAR(sample:2000-2025)",
    subtitle = "Impulse Responses to a 1-Standard-Deviation Pure MP Shock",
    caption = sprintf("Effective Observations (T): %d | Lags (p): %d", T_eff, prior$lags),
    
    x = "Months after Shock",
    y = "Response"
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    strip.text       = element_text(face = "bold", size = 12),
    plot.caption     = element_text(face = "italic", size = 10, color = "gray30"), 
    panel.grid.minor = element_blank()
  )

print(p_baseline)

# =============================
# Extraxt structural shocks
# =============================
extract_structural_shocks <- function(res, irfs_draws, Y, prior, shock_idx = 1) {
  P <- prior$lags
  T_full <- nrow(Y)
  T_eff <- T_full - P
  
  Y_t <- Y[(P + 1):T_full, , drop = FALSE]
  X_t <- matrix(1, nrow = T_eff, ncol = 1)
  for (l in 1:P) {
    X_t <- cbind(X_t, Y[(P + 1 - l):(T_full - l), , drop = FALSE])
  }
  X_t_reordered <- cbind(X_t[, -1, drop = FALSE], X_t[, 1, drop = FALSE])
  
  ndraws <- dim(res$beta_draws)[3]
  N <- ncol(Y)
  
  shock_history_draws <- matrix(NA, nrow = T_eff, ncol = ndraws)
  
  for (i in 1:ndraws) {
    beta_i <- res$beta_draws[, , i]
    U_i <- Y_t - X_t_reordered %*% beta_i
    
    A0_Q <- irfs_draws[, , 1, i]
    
    if (any(is.na(A0_Q))) next
    
    epsilon_i <- t(solve(A0_Q, t(U_i)))
    
    shock_history_draws[, i] <- epsilon_i[, shock_idx]
  }
  
  median_shock_path <- apply(shock_history_draws, 1, median, na.rm = TRUE)
  
  return(median_shock_path)
}

median_mp_shock <- extract_structural_shocks(res, irfs_draws, Y_matrix, prior, shock_idx = 1)

effective_dates_idx <- (prior$lags + 1):nrow(var_data)
var_data_eff <- var_data[effective_dates_idx, ]
var_data_eff$`Pure MP Shock` <- median_mp_shock


