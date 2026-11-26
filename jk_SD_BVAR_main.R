####################################################################
# State Dependent version of Jarocinski & Karadi (2020) - Euro Area
####################################################################
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("jk_SD_BVAR_functions.R")

#======
# Data
#======
var_data <- final_data %>% filter(Country == "EA")%>% 
  select(Date, `Ind. Production`=Log_IP, `Diff. Oil Prices`=df_Oil,`True Policy Rate` =True_Policy_Rate, 
         Log_GPR,Log_GPRA, Log_GPRT,`HICP`=Log_HICP,`Core HICP`=Log_Core_HICP, `Diff. HICP`= df_HICP, `Diff. 1Y Infl. Linked Swap` = df_ILS_1Y,
         `Oil Prices` = Log_Oil, `EONIA/€STR 3M`=Spliced_EONIA_3M, , ESI, OIS_3M, Log_GPRT, Log_GPRA, `Diff. 5Y5Y Infl. Linked Swap`= df_ILS_5Y5Y,
         `OIS 1Y surprise`=OIS_1Y, STOXX50,  Log_VIX,`EONIA/€STR 1Y`=Spliced_EONIA_1Y, Log_VSTOXX, `Diff. EONIA/€STR 1Y`= df_EONIA_1Y,
          Log_EPU, `Diff. Unemployment Rate`= df_Unemployment_Rate, `Diff. Ind. Prod.`=df_IP)%>% # Adding PCCI or inflation linked swaps gets rid of  observations
  drop_na()

#==================
# Settings & Setup
#==================
ynames1 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP") #done
ynames2 <- c("EONIA/€STR 1Y","Ind. Production","Oil Prices" ,"Core HICP")#done
ynames3 <- c("EONIA/€STR 1Y","Ind. Production","Oil Prices" ,"PCCI")#done
ynames4 <- c("EONIA/€STR 1Y", "ESI","Oil Prices" ,"HICP")#done
ynames5 <- c("EONIA/€STR 1Y", "Diff. Unemployment Rate","Oil Prices", "HICP")#done
ynames6 <- c("True Policy Rate", "Ind. Production", "Oil Prices", "HICP")#done
ynames7 <- c("EONIA/€STR 3M", "Ind. Production","Oil Prices", "HICP") #done
ynames8 <- c("EONIA/€STR 2Y", "Ind. Production","Oil Prices", "HICP") #done
ynames9 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP", "5Y5Y Infl. Linked Swap") #done
ynames10 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP", "1Y Infl. Linked Swap") #done
ynames11 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP", "Diff. 5Y5Y Infl. Linked Swap") #done
ynames12 <- c("EONIA/€STR 1Y", "Ind. Production","Oil Prices", "HICP", "Diff. 1Y Infl. Linked Swap") #done
ynames13 <- c("Diff. EONIA/€STR 1Y", "Diff. Ind. Prod.","Diff. Oil Prices", "Diff. HICP")#done
#===========
# Variables
#===========
mnames <- c("OIS_3M", "STOXX50")
ynames <- ynames12
thresh_var <- var_data$Log_GPR
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
  ndraws = 2000,
  burnin = 2000,
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
c_speed <- 10
gamma_threshold <- quantile(thresh_var,0.7, na.rm = TRUE)

z_t <- thresh_var[(prior$lags + 1 - d_delay):(nrow(var_data) - d_delay)]

F_high <- 1 / (1 + exp(-c_speed * (z_t - gamma_threshold) / sd(z_t)))
F_low  <- 1 - F_high

#==============================================================================
# Estimate the State-Dependent VARs
#==============================================================================
cat("\nEstimating Low GPR Regime...\n")
res_low  <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = F_low)

cat("\nEstimating High GPR Regime...\n")
res_high <- VAR_withiid_weighted(Y_matrix, prior, gssettings, weights = F_high)

#--------------------------------------
# Identification with Sign Restriction
#--------------------------------------
MAlags <- 36

dims <- c(1, 2)
imonpol <- 1
inews <- 2

test_restr <- function(irfs) {
  (irfs[1, imonpol, 1] > 0) && (irfs[2, imonpol, 1] < 0) && 
    (irfs[1, inews, 1] > 0)   && (irfs[2, inews, 1] > 0)
}

b_normalize <- rep(1, N)
max_try <- 1000

cat("\nApplying Sign Restrictions - Low GPR Regime...\n")
irfs_draws_low  <- resirfssign(res_low, MAlags, dims, test_restr, b_normalize, max_try)

cat("\nApplying Sign Restrictions - High GPR Regime...\n")
irfs_draws_high <- resirfssign(res_high, MAlags, dims, test_restr, b_normalize, max_try)

#==============================================================================
# Extract, Merge, and Plot
#==============================================================================
df_low  <- extract_matlab_irfs(irfs_draws_low, data_names, shock_idx = 1)
df_high <- extract_matlab_irfs(irfs_draws_high, data_names, shock_idx = 1)

df_low$Regime  <- "Low GPR"
df_high$Regime <- "High GPR"
plot_data <- rbind(df_low, df_high)

plot_data <- plot_data %>%
  filter(Variable %in% ynames) %>%
  mutate(Variable = factor(Variable, levels = ynames))

plot_data$Horizon <- as.numeric(plot_data$Horizon)

T_eff <- nrow(Y_matrix) - prior$lags

p_st <- ggplot(plot_data, aes(x = Horizon, y = Median, color = Regime, fill = Regime)) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  
  geom_ribbon(aes(ymin = Lower68, ymax = Upper68), alpha = 0.2, color = NA) +
  
  geom_line(linewidth = 1) +
  
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  
  scale_color_manual(values = c("Low GPR" = "steelblue", "High GPR" = "firebrick")) +
  scale_fill_manual(values = c("Low GPR" = "steelblue", "High GPR" = "firebrick")) +
  
  labs(
    title = "State-Dependent BVAR: Pure MP Shock (diff. 1Y inflation swap model)",
    subtitle = "HF & Sign Restriction Identification: High vs. Low Geopolitical Risk Regimes",
    caption = sprintf("Effective Observations (T): %d | Lags (p): %d", T_eff, prior$lags),
    x = "Months after Shock",
    y = "Impulse Response"
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(face = "bold", size = 16),
    plot.caption = element_text(face = "italic", size = 10, color = "gray30")
  )

print(p_st)

