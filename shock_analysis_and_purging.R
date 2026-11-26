###############################################################
# UNIT ROOT TESTING
###############################################################
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("TS_functions.R")
source("additional_functions.R")
source("LeeStrazicichUnitRootTest.R")
#=====================
# Prepare & plot Data
#=====================
write_csv(final_data, "data_thesis_r0838656.csv")
write_xlsx(final_data, "data_thesis_r0838656.xlsx")
unit_root_data <- final_data %>% filter(Country == "EA")%>% 
  select(Date, Log_IP, Log_HICP, Log_Oil, Spliced_EONIA_1Y, OIS_3M, OIS_1Y, STOXX50, Target_Shock, FG_Shock, QE_Shock)%>%
  drop_na()


OIS_3M              <- ts(unit_root_data$OIS_3M, start = c(2000, 2), frequency = 12)
OIS_1Y              <- ts(unit_root_data$OIS_1Y, start = c(2000, 2), frequency = 12)
STOXX50             <- ts(unit_root_data$STOXX50, start = c(2000, 2), frequency = 12)
Target_Shock             <- ts(unit_root_data$Target_Shock, start = c(2000, 2), frequency = 12)
FG_Shock             <- ts(unit_root_data$FG_Shock, start = c(2000, 2), frequency = 12)
QE_Shock             <- ts(unit_root_data$QE_Shock, start = c(2000, 2), frequency = 12)
HICP                <- ts(unit_root_data$Log_HICP, start = c(2000, 2), frequency = 12)
IP                  <- ts(unit_root_data$Log_IP, start = c(2000, 2), frequency = 12)
Oil_prices          <- ts(unit_root_data$Log_Oil, start = c(2000, 2), frequency = 12)
EONIA_1Y            <- ts(unit_root_data$Spliced_EONIA_1Y, start = c(2000, 2), frequency = 12)


ts_data <- ts.union(OIS_3M, OIS_1Y, STOXX50, Target_Shock, FG_Shock, QE_Shock, 
                    HICP, IP, Oil_prices, EONIA_1Y)

#--------
# OIS_3M
#--------

ARMA_IC_selection_CSS(15,0,OIS_3M) #BIC suggests AR(2)
AR_OIS_3M <- arima(OIS_3M, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_residuals_OIS_3M <- residuals(AR_OIS_3M)

## Ljung-Box
Qstat_AR_residuals_OIS_3M <- LjungBox(AR_residuals_OIS_3M, 20, 2)
stargazer(Qstat_AR_residuals_OIS_3M, type = "text", summary = FALSE) 

#ADF(1) (p-1)
adf_OIS_3M <- ur.df(OIS_3M, type = "trend", lags = 1)
summary(adf_OIS_3M)
# gamma=0? no, no unit root

LB_test_OIS_3M <- LjungBox(OIS_3M, 12, 0)
stargazer(LB_test_OIS_3M, type = "text", summary = FALSE) # serial correlation

arch_test_OIS_3M <- ArchTest(OIS_3M, lags = 6, demean = TRUE)
print(arch_test_OIS_3M) # There is heteroskedasticity


#--------
# OIS_1Y
#--------

ARMA_IC_selection_CSS(15,0,OIS_1Y) #BIC suggests AR(2)
AR_OIS_1Y <- arima(OIS_1Y, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_residuals_OIS_1Y <- residuals(AR_OIS_1Y)

## Ljung-Box
Qstat_AR_residuals_OIS_1Y <- LjungBox(AR_residuals_OIS_1Y, 20, 2)
stargazer(Qstat_AR_residuals_OIS_1Y, type = "text", summary = FALSE) 

#ADF(1) (p-1)
adf_OIS_1Y <- ur.df(OIS_1Y, type = "trend", lags = 1)
summary(adf_OIS_1Y)
# gamma=0? no, no unit root

LB_test_OIS_1Y <- LjungBox(OIS_1Y, 12, 0)
stargazer(LB_test_OIS_1Y, type = "text", summary = FALSE) # serial correlation

arch_test_OIS_1Y <- ArchTest(OIS_1Y, lags = 6, demean = TRUE)
print(arch_test_OIS_1Y) # There is heteroskedasticity


#--------
# STOXX50
#--------

ARMA_IC_selection_CSS(15,0,STOXX50) #BIC suggests AR(0)
AR_STOXX50 <- arima(STOXX50, order = c(0, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_residuals_STOXX50 <- residuals(AR_STOXX50)

## Ljung-Box
Qstat_AR_residuals_STOXX50 <- LjungBox(AR_residuals_STOXX50, 20, 0)
stargazer(Qstat_AR_residuals_STOXX50, type = "text", summary = FALSE) 

#ADF(0) (p-1)
adf_STOXX50 <- ur.df(STOXX50, type = "trend", lags = 0)
summary(adf_STOXX50)
# gamma=0? no, no unit root

LB_test_STOXX50 <- LjungBox(STOXX50, 12, 0)
stargazer(LB_test_STOXX50, type = "text", summary = FALSE) # No serial correlation

arch_test_STOXX50 <- ArchTest(STOXX50, lags = 6, demean = TRUE)
print(arch_test_STOXX50) # There is no heteroskedasticity


#-------------
# Target_Shock
#-------------

ARMA_IC_selection_CSS(15,0,Target_Shock) #BIC suggests AR(3)
AR_Target_Shock <- arima(Target_Shock, order = c(3, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_residuals_Target_Shock <- residuals(AR_Target_Shock)

## Ljung-Box
Qstat_AR_residuals_Target_Shock <- LjungBox(AR_residuals_Target_Shock, 20, 3)
stargazer(Qstat_AR_residuals_Target_Shock, type = "text", summary = FALSE) 

#ADF(2) (p-1)
adf_Target_Shock <- ur.df(Target_Shock, type = "trend", lags = 2)
summary(adf_Target_Shock)
# gamma=0? no, no unit root

LB_test_Target_Shock <- LjungBox(Target_Shock, 12, 0)
stargazer(LB_test_Target_Shock, type = "text", summary = FALSE) # serial correlation

arch_test_Target_Shock <- ArchTest(Target_Shock, lags = 6, demean = TRUE)
print(arch_test_Target_Shock) # There is heteroskedasticity


#-----------
# FG_Shock
#-----------

ARMA_IC_selection_CSS(15,0,FG_Shock) #BIC suggests AR(0), next AR(2), next AR(4)
AR_FG_Shock <- arima(FG_Shock, order = c(4, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_residuals_FG_Shock <- residuals(AR_FG_Shock)

## Ljung-Box
Qstat_AR_residuals_FG_Shock <- LjungBox(AR_residuals_FG_Shock, 20, 4)
stargazer(Qstat_AR_residuals_FG_Shock, type = "text", summary = FALSE)# still serial correlation with AR(0) and at A(2)

#ADF(3) (p-1)
adf_FG_Shock <- ur.df(FG_Shock, type = "trend", lags = 3)
summary(adf_FG_Shock)
# gamma=0? no, no unit root

LB_test_FG_Shock <- LjungBox(FG_Shock, 12, 0)
stargazer(LB_test_FG_Shock, type = "text", summary = FALSE) # serial correlation

arch_test_FG_Shock <- ArchTest(FG_Shock, lags = 6, demean = TRUE)
print(arch_test_FG_Shock) # There is heteroskedasticity


#--------
# QE_Shock
#--------

ARMA_IC_selection_CSS(15,0,QE_Shock) #BIC suggests AR(0)
AR_QE_Shock <- arima(QE_Shock, order = c(0, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_residuals_QE_Shock <- residuals(AR_QE_Shock)

## Ljung-Box
Qstat_AR_residuals_QE_Shock <- LjungBox(AR_residuals_QE_Shock, 20, 0)
stargazer(Qstat_AR_residuals_QE_Shock, type = "text", summary = FALSE) 

#ADF(0) (p-1)
adf_QE_Shock <- ur.df(QE_Shock, type = "trend", lags = 0)
summary(adf_QE_Shock)
# gamma=0? no, no unit root

LB_test_QE_Shock <- LjungBox(QE_Shock, 12, 0)
stargazer(LB_test_QE_Shock, type = "text", summary = FALSE) # No serial correlation

arch_test_QE_Shock <- ArchTest(QE_Shock, lags = 6, demean = TRUE)
print(arch_test_QE_Shock) # no heteroskedasticity







#-----------------------
# Exogeneity: Target_Shock
#-----------------------
exo_Target <- dynlm(Target_Shock ~ L(IP, 1)+L(IP,2) + L(HICP, 1)+L(HICP,2) + L(Oil_prices, 1)+L(Oil_prices, 2) + L(EONIA_1Y, 1)+L(EONIA_1Y, 2), data = ts_data)
summary_exo_Target <- summary(exo_Target)
print(summary_exo_Target)

f_stat_Target <- summary_exo_Target$fstatistic
p_val_Target <- pf(f_stat_Target[1], f_stat_Target[2], f_stat_Target[3], lower.tail = FALSE)
cat("Target_Shock Exogeneity p-value:", p_val_Target, "\n") 

#-----------------------
# Exogeneity: FG_Shock
#-----------------------
exo_FG <- dynlm(FG_Shock ~ L(IP, 1)+L(IP,2) + L(HICP, 1)+L(HICP,2) + L(Oil_prices, 1)+L(Oil_prices, 2) + L(EONIA_1Y, 1)+L(EONIA_1Y, 2), data = ts_data)
summary_exo_FG <- summary(exo_FG)
print(summary_exo_FG)

f_stat_FG <- summary_exo_FG$fstatistic
p_val_FG <- pf(f_stat_FG[1], f_stat_FG[2], f_stat_FG[3], lower.tail = FALSE)
cat("FG_Shock Exogeneity p-value:", p_val_FG, "\n")

#-----------------------
# Exogeneity: QE_Shock
#-----------------------
exo_QE <- dynlm(QE_Shock ~ L(IP, 1)+L(IP,2) + L(HICP, 1)+L(HICP,2) + L(Oil_prices, 1)+L(Oil_prices, 2) + L(EONIA_1Y, 1)+L(EONIA_1Y, 2), data = ts_data)
summary_exo_QE <- summary(exo_QE)
print(summary_exo_QE)

f_stat_QE <- summary_exo_QE$fstatistic
p_val_QE <- pf(f_stat_QE[1], f_stat_QE[2], f_stat_QE[3], lower.tail = FALSE)
cat("QE_Shock Exogeneity p-value:", p_val_QE, "\n")

#-----------------------
# Exogeneity: OIS_3M, OIS_1Y, STOXX50
#-----------------------
exo_OIS_3M <- dynlm(OIS_3M ~ L(IP, 1)+L(IP,2) + L(HICP, 1)+L(HICP,2) + L(Oil_prices, 1)+L(Oil_prices, 2) + L(EONIA_1Y, 1)+L(EONIA_1Y, 2), data = ts_data)
p_val_OIS_3M <- pf(summary(exo_OIS_3M)$fstatistic[1], summary(exo_OIS_3M)$fstatistic[2], summary(exo_OIS_3M)$fstatistic[3], lower.tail = FALSE)
cat("OIS_3M Exogeneity p-value:", p_val_OIS_3M, "\n")

exo_OIS_1Y <- dynlm(OIS_1Y ~ L(IP, 1)+L(IP,2) + L(HICP, 1)+L(HICP,2) + L(Oil_prices, 1)+L(Oil_prices, 2) + L(EONIA_1Y, 1)+L(EONIA_1Y, 2), data = ts_data)
p_val_OIS_1Y <- pf(summary(exo_OIS_1Y)$fstatistic[1], summary(exo_OIS_1Y)$fstatistic[2], summary(exo_OIS_1Y)$fstatistic[3], lower.tail = FALSE)
cat("OIS_1Y Exogeneity p-value:", p_val_OIS_1Y, "\n")

exo_STOXX50 <- dynlm(STOXX50 ~ L(IP, 1)+L(IP,2) + L(HICP, 1)+L(HICP,2) + L(Oil_prices, 1)+L(Oil_prices, 2) + L(EONIA_1Y, 1)+L(EONIA_1Y, 2), data = ts_data)
p_val_STOXX50 <- pf(summary(exo_STOXX50)$fstatistic[1], summary(exo_STOXX50)$fstatistic[2], summary(exo_STOXX50)$fstatistic[3], lower.tail = FALSE)
cat("STOXX50 Exogeneity p-value:", p_val_STOXX50, "\n")
#==================================
# Correlation Matrix
#=================================
library(corrplot)
shocks_matrix <- cbind(Target_Shock, FG_Shock, QE_Shock, OIS_3M, OIS_1Y, STOXX50)
colnames(shocks_matrix) <- c("Target", "FG", "QE", "OIS 3M", "OIS 1Y", "STOXX50")

ortho_cor_matrix <- cor(shocks_matrix, use = "complete.obs")

print("Correlation Matrix of Shocks (Check for values close to 0 between structural shocks):")
print(round(ortho_cor_matrix, 3))

corrplot(ortho_cor_matrix, method = "number", type = "upper", 
         tl.col = "black", tl.srt = 60, 
         title = "Orthogonality of Monetary Policy Surprises", mar=c(0,0,1,0))





summary_table <- data.frame(
  Variable = c("OIS_3M", "OIS_1Y", "STOXX50", "Target_Shock", "FG_Shock", "QE_Shock"),
  ADF_Test_Stat = c(
    adf_OIS_3M@teststat[1], 
    adf_OIS_1Y@teststat[1], 
    adf_STOXX50@teststat[1],
    adf_Target_Shock@teststat[1], 
    adf_FG_Shock@teststat[1], 
    adf_QE_Shock@teststat[1]
  ),
  ARCH_p_value = round(c(
    arch_test_OIS_3M$p.value, 
    arch_test_OIS_1Y$p.value, 
    arch_test_STOXX50$p.value,
    arch_test_Target_Shock$p.value, 
    arch_test_FG_Shock$p.value, 
    arch_test_QE_Shock$p.value
  ), 4),
  Exogeneity_p_value = round(c(
    p_val_OIS_3M, 
    p_val_OIS_1Y, 
    p_val_STOXX50, 
    p_val_Target, 
    p_val_FG, 
    p_val_QE
  ), 4)
)

print("--- MAIN DIAGNOSTIC TESTS SUMMARY ---")
print(summary_table)


LB_test_OIS_3M <- LjungBox(OIS_3M, 12, 0)
stargazer(LB_test_OIS_3M, type = "text", summary = FALSE)
LB_test_OIS_1Y <- LjungBox(OIS_1Y, 12, 0)
stargazer(LB_test_OIS_1Y, type = "text", summary = FALSE)
LB_test_STOXX50 <- LjungBox(STOXX50, 12, 0)
stargazer(LB_test_STOXX50, type = "text", summary = FALSE)
LB_test_Target_Shock <- LjungBox(Target_Shock, 12, 0)
stargazer(LB_test_Target_Shock, type = "text", summary = FALSE)
LB_test_FG_Shock <- LjungBox(FG_Shock, 12, 0)
stargazer(LB_test_FG_Shock, type = "text", summary = FALSE)
LB_test_QE_Shock <- LjungBox(QE_Shock, 12, 0)
stargazer(LB_test_QE_Shock, type = "text", summary = FALSE)






#----------------------------------
# purging and adding to final data
#----------------------------------
unit_root_data <- unit_root_data %>%
  mutate(
    Purged_OIS_3M       = as.numeric(AR_residuals_OIS_3M),
    Purged_OIS_1Y       = as.numeric(AR_residuals_OIS_1Y),
    Purged_Target_Shock = as.numeric(AR_residuals_Target_Shock),
    Purged_FG_Shock     = as.numeric(AR_residuals_FG_Shock),
  )

purged_shocks_only <- unit_root_data %>%
  select(Date, starts_with("Purged_"))

final_data <- final_data %>%
  left_join(purged_shocks_only, by = "Date")


