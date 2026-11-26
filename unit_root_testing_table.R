#==================
# Unit Root Testing
#==================


source("LeeStrazizichUnitRootTest.R") # Ensure this is sourced in your session
unit_root_data <- final_data %>% filter(Country == "EA")%>% 
  select(Date,df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, Log_GPRA, Log_GPRT,
         Log_HICP,Log_Core_HICP, Log_Oil, Spliced_EONIA_3M, Log_VSTOXX,PCCI,
         df_Core_HICP, df_True_Policy_Rate, Spliced_EONIA_1Y, Unemployment_Rate
         , ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX)%>% # Adding PCCI gets rid of 14 observations
  drop_na()
#-----
# HICP
#-----
ARMA_IC_selection_CSS(15,0,HICP) 
AR_HICP <- arima(HICP, order = c(6, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) 
AR_HICP_residuals <- residuals(AR_HICP)
Qstat_AR_HICP_residuals <- LjungBox(AR_HICP_residuals,20,6)
stargazer(Qstat_AR_HICP_residuals,type = "text", summary = FALSE) 

adf_HICP <- ur.df(HICP, type = "trend", lags =5)
summary(adf_HICP)
adf_HICP <- ur.df(HICP, type = "drift", lags =5)
summary(adf_HICP)

time_trend <- time(HICP)
bp_HICP <- breakpoints(HICP ~ time_trend, h = 0.15)
summary(bp_HICP) 
breakdates(bp_HICP)
plot(HICP, main = "Bai-Perron Structural Breaks: HICP")
lines(bp_HICP)

za_HICP <- ur.za(HICP, model = "both", lag = 2)
summary(za_HICP) 
plot(za_HICP) 

ls_HICP <- ur.ls(y = HICP, model = "break", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" ) 

#----------
# Core HICP
#----------
ARMA_IC_selection_CSS(15,0,core_HICP) 
AR_core_HICP <- arima(core_HICP, order = c(7, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) 
AR_core_HICP_residuals <- residuals(AR_core_HICP)
Qstat_AR_core_HICP_residuals <- LjungBox(AR_core_HICP_residuals, 20, 7)
stargazer(Qstat_AR_core_HICP_residuals,type = "text", summary = FALSE) 

adf_core_HICP <- ur.df(core_HICP, type = "trend", lags =6)
summary(adf_core_HICP)
adf_core_HICP <- ur.df(core_HICP, type = "drift", lags =6)
summary(adf_core_HICP)
adf_core_HICP <- ur.df(core_HICP, type = "none", lags =6)
summary(adf_core_HICP)

time_trend <- time(core_HICP)
bp_core_HICP <- breakpoints(core_HICP ~ time_trend, h = 0.15)
summary(bp_core_HICP) 
breakdates(bp_core_HICP)
plot(core_HICP, main = "Bai-Perron Structural Breaks: Core HICP")
lines(bp_core_HICP)

za_core_HICP <- ur.za(core_HICP, model = "both", lag = 2)
summary(za_core_HICP) 
plot(za_core_HICP)

ls_core_HICP <- ur.ls(y = core_HICP, model = "break", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#------------
# PCCI
#------------
ARMA_IC_selection_CSS(15,0,PCCI) 
AR_PCCI <- arima(PCCI, order=c(11,0,0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_PCCI_residuals <- residuals(AR_PCCI)
Qstat_AR_PCCI_residuals <- LjungBox(AR_PCCI_residuals,20,11)
stargazer(Qstat_AR_PCCI_residuals,type = "text", summary = FALSE)

adf_PCCI <- ur.df(PCCI, type = "trend", lags =10)
summary(adf_PCCI)

bp_PCCI <- breakpoints(PCCI ~ 1, h = 0.15)
summary(bp_PCCI)
breakdates(bp_PCCI)
plot(PCCI, main = "Bai-Perron Structural Breaks: PCCI")
lines(bp_PCCI)

za_PCCI <- ur.za(PCCI, model = "both", lag = 2)
summary(za_PCCI)
plot(za_PCCI)

ls_PCCI <- ur.ls(y = PCCI, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#-----
# IP
#-----
ARMA_IC_selection_CSS(15,0,IP) 
AR_IP <- arima(IP, order = c(3, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_IP_residuals <- residuals(AR_IP)
Qstat_AR_IP_residuals <- LjungBox(AR_IP_residuals,20,3)
stargazer(Qstat_AR_IP_residuals,type = "text", summary = FALSE)

adf_IP <- ur.df(IP, type = "trend", lags =2)
summary(adf_IP)

time_trend <- time(IP)
bp_IP <- breakpoints(IP ~ time_trend, h = 0.15)
summary(bp_IP) 
breakdates(bp_IP)
plot(IP, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_IP)

za_IP <- ur.za(IP, model = "both", lag = 2)
summary(za_IP) 
plot(za_IP)

ls_IP <- ur.ls(y = IP, model = "break", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#------------------
# Unemployment Rate
#------------------
ARMA_IC_selection_CSS(15,0,Unemployment_Rate) 
AR_Unemployment_Rate <- arima(Unemployment_Rate, order = c(8, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_Unemployment_Rate_residuals <- residuals(AR_Unemployment_Rate)
Qstat_AR_Unemployment_Rate_residuals <- LjungBox(AR_Unemployment_Rate_residuals,20,8)
stargazer(Qstat_AR_Unemployment_Rate_residuals,type = "text", summary = FALSE)

adf_Unemployment_Rate <- ur.df(Unemployment_Rate, type = "trend", lags =7)
summary(adf_Unemployment_Rate)
adf_Unemployment_Rate <- ur.df(Unemployment_Rate, type = "drift", lags =7)
summary(adf_Unemployment_Rate)
adf_Unemployment_Rate <- ur.df(Unemployment_Rate, type = "none", lags =7)
summary(adf_Unemployment_Rate)

bp_Unemployment_Rate <- breakpoints(Unemployment_Rate ~ 1, h = 0.15)
summary(bp_Unemployment_Rate) 
breakdates(bp_Unemployment_Rate)
plot(Unemployment_Rate, main = "Bai-Perron Structural Breaks: Unemployment Rate")
lines(bp_Unemployment_Rate)

za_Unemployment_Rate <- ur.za(Unemployment_Rate, model = "both", lag = 2)
summary(za_Unemployment_Rate) 
plot(za_Unemployment_Rate)

ls_Unemployment_Rate <- ur.ls(y = Unemployment_Rate, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#-------
# ESI
#-------
ARMA_IC_selection_CSS(15,0,ESI) 
AR_ESI <- arima(ESI, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_ESI_residuals <- residuals(AR_ESI)
Qstat_AR_ESI_residuals <- LjungBox(AR_ESI_residuals,20,2)
stargazer(Qstat_AR_ESI_residuals,type = "text", summary = FALSE)

adf_ESI <- ur.df(ESI, type = "trend", lags =1)
summary(adf_ESI)

bp_ESI <- breakpoints(ESI ~ 1, h = 0.15)
summary(bp_ESI) 
breakdates(bp_ESI)
plot(ESI, main = "Bai-Perron Structural Breaks: ESI")
lines(bp_ESI)

ls_ESI <- ur.ls(y = ESI, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )


#-----------
# Oil Prices
#-----------
ARMA_IC_selection_CSS(15,0,Oil_prices)
AR_oil <- arima(Oil_prices, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) 
AR_oil_residuals <- residuals(AR_oil)
Qstat_AR_oil_residuals <- LjungBox(AR_oil_residuals,20,2)
stargazer(Qstat_AR_oil_residuals,type = "text", summary = FALSE) 

adf_Oil_prices <- ur.df(Oil_prices, type = "trend", lags =1)
summary(adf_Oil_prices)
adf_Oil_prices <- ur.df(Oil_prices, type = "drift", lags =1)
summary(adf_Oil_prices)
adf_Oil_prices <- ur.df(Oil_prices, type = "none", lags =1)
summary(adf_Oil_prices)

bp_Oil_prices <- breakpoints(Oil_prices ~ 1, h = 0.15)
summary(bp_Oil_prices) 
breakdates(bp_Oil_prices)
plot(Oil_prices, main = "Bai-Perron Structural Breaks: Oil Prices")
lines(bp_Oil_prices)

ls_Oil_prices <- ur.ls(y = Oil_prices, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )


#------------
# EONIA_1Y
#------------
ARMA_IC_selection_CSS(15,0,EONIA_1Y) 
AR_IR <- arima(EONIA_1Y, order = c(3, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) 
AR_IR_residuals <- residuals(AR_IR)
Qstat_AR_IR_residuals <- LjungBox(AR_IR_residuals,20,3)
stargazer(Qstat_AR_IR_residuals,type = "text", summary = FALSE) 

adf_EONIA_1Y <- ur.df(EONIA_1Y, type = "trend", lags =2)
summary(adf_EONIA_1Y)
adf_EONIA_1Y <- ur.df(EONIA_1Y, type = "drift", lags =2)
summary(adf_EONIA_1Y)
adf_EONIA_1Y <- ur.df(EONIA_1Y, type = "none", lags =2)
summary(adf_EONIA_1Y)

bp_EONIA_1Y <- breakpoints(EONIA_1Y ~ 1, h = 0.15)
summary(bp_EONIA_1Y) 
breakdates(bp_EONIA_1Y)
plot(EONIA_1Y, main = "Bai-Perron Structural Breaks: EONIA 1Y")
lines(bp_EONIA_1Y)

za_EONIA_1Y <- ur.za(EONIA_1Y, model = "both", lag = 2)
summary(za_EONIA_1Y)
plot(za_EONIA_1Y)

ls_EONIA_1Y <- ur.ls(y = EONIA_1Y, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )


#-----------------
# True Policy rate
#-----------------
ARMA_IC_selection_CSS(15,0,True_Policy_Rate)
AR_True_Policy_Rate <- arima(True_Policy_Rate, order = c(4, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_True_Policy_Rate_residuals <- residuals(AR_True_Policy_Rate)
Qstat_AR_True_Policy_Rate_residuals <- LjungBox(AR_True_Policy_Rate_residuals,20,4)
stargazer(Qstat_AR_True_Policy_Rate_residuals,type = "text", summary = FALSE)

AR_True_Policy_Rate <- arima(True_Policy_Rate, order = c(9, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_True_Policy_Rate_residuals <- residuals(AR_True_Policy_Rate)
Qstat_AR_True_Policy_Rate_residuals <- LjungBox(AR_True_Policy_Rate_residuals,20,9)
stargazer(Qstat_AR_True_Policy_Rate_residuals,type = "text", summary = FALSE)

adf_True_Policy_Rate <- ur.df(True_Policy_Rate, type = "trend", lags =8)
summary(adf_True_Policy_Rate)
adf_True_Policy_Rate <- ur.df(True_Policy_Rate, type = "drift", lags =8)
summary(adf_True_Policy_Rate)
adf_True_Policy_Rate <- ur.df(True_Policy_Rate, type = "none", lags =8)
summary(adf_True_Policy_Rate)

bp_True_Policy_Rate <- breakpoints(True_Policy_Rate ~ 1, h = 0.15)
summary(bp_True_Policy_Rate) 
breakdates(bp_True_Policy_Rate)
plot(True_Policy_Rate, main = "Bai-Perron Structural Breaks: True Policy Rate")
lines(bp_True_Policy_Rate)

za_True_Policy_Rate <- ur.za(True_Policy_Rate, model = "both", lag = 2)
summary(za_True_Policy_Rate)
plot(za_True_Policy_Rate)

ls_True_Policy_Rate <- ur.ls(y = True_Policy_Rate, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#------
# GPR
#------
ARMA_IC_selection_CSS(15,0,GPR) 
AR_GPR <- arima(GPR, order = c(1, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) 
AR_GPR_residuals <- residuals(AR_GPR)
Qstat_AR_GPR_residuals <- LjungBox(AR_GPR_residuals,20,1)
stargazer(Qstat_AR_GPR_residuals,type = "text", summary = FALSE) 

adf_GPR <- ur.df(GPR, type = "trend", lags =0)
summary(adf_GPR)

bp_GPR <- breakpoints(GPR ~ 1, h = 0.15)
summary(bp_GPR) 
breakdates(bp_GPR)
plot(GPR, main = "Bai-Perron Structural Breaks: GPR")
lines(bp_GPR)

za_GPR <- ur.za(GPR, model = "both", lag = 2)
summary(za_GPR) 
plot(za_GPR)

ls_GPR <- ur.ls(y = GPR, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" ) 

#------
# VIX
#------
ARMA_IC_selection_CSS(15,0,VIX) 
AR_VIX <- arima(VIX, order = c(1, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_VIX_residuals <- residuals(AR_VIX)
Qstat_AR_VIX_residuals <- LjungBox(AR_VIX_residuals,20,1)
stargazer(Qstat_AR_VIX_residuals,type = "text", summary = FALSE)

adf_VIX <- ur.df(VIX, type = "trend", lags =0)
summary(adf_VIX)

bp_VIX <- breakpoints(VIX ~ 1, h = 0.15)
summary(bp_VIX) 
breakdates(bp_VIX)
plot(VIX, main = "Bai-Perron Structural Breaks: VIX")
lines(bp_VIX)

za_VIX <- ur.za(VIX, model = "both", lag = 2)
summary(za_VIX)
plot(za_VIX)

ls_VIX <- ur.ls(y = VIX, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#--------
# VSTOXX
#--------
ARMA_IC_selection_CSS(15,0,VSTOXX) 
AR_VSTOXX <- arima(VSTOXX, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_VSTOXX_residuals <- residuals(AR_VSTOXX)
Qstat_AR_VSTOXX_residuals <- LjungBox(AR_VSTOXX_residuals,20,2)
stargazer(Qstat_AR_VSTOXX_residuals,type = "text", summary = FALSE)

adf_VSTOXX <- ur.df(VSTOXX, type = "trend", lags =1)
summary(adf_VSTOXX)

bp_VSTOXX <- breakpoints(VSTOXX ~ 1, h = 0.15)
summary(bp_VSTOXX) 
breakdates(bp_VSTOXX)
plot(VSTOXX, main = "Bai-Perron Structural Breaks: VSTOXX")
lines(bp_VSTOXX)

za_VSTOXX <- ur.za(VSTOXX, model = "both", lag = 2)
summary(za_VSTOXX)
plot(za_VSTOXX)

ls_VSTOXX <- ur.ls(y = VSTOXX, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#------------
# EU_ILS_1Y
#------------
ARMA_IC_selection_CSS(15,0,EU_ILS_1Y) 
AR_EU_ILS_1Y <- arima(EU_ILS_1Y, order = c(3, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_EU_ILS_1Y_residuals <- residuals(AR_EU_ILS_1Y)
Qstat_AR_EU_ILS_1Y_residuals <- LjungBox(AR_EU_ILS_1Y_residuals,20,3)
stargazer(Qstat_AR_EU_ILS_1Y_residuals,type = "text", summary = FALSE)

adf_EU_ILS_1Y <- ur.df(EU_ILS_1Y, type = "trend", lags = 2)
summary(adf_EU_ILS_1Y)

adf_EU_ILS_1Y <- ur.df(EU_ILS_1Y, type = "drift", lags = 2)
summary(adf_EU_ILS_1Y)

adf_EU_ILS_1Y <- ur.df(EU_ILS_1Y, type = "none", lags = 2)
summary(adf_EU_ILS_1Y)

bp_EU_ILS_1Y <- breakpoints(EU_ILS_1Y ~ 1, h = 0.15)
summary(bp_EU_ILS_1Y) 
breakdates(bp_EU_ILS_1Y)
plot(EU_ILS_1Y, main = "Bai-Perron Structural Breaks: EU ILS 1Y")
lines(bp_EU_ILS_1Y)

za_EU_ILS_1Y <- ur.za(EU_ILS_1Y, model = "both", lag = 2)
summary(za_EU_ILS_1Y)
plot(za_EU_ILS_1Y)

ls_EU_ILS_1Y <- ur.ls(y = EU_ILS_1Y, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )


#--------------
# EU_ILS_5Y5Y
#--------------
ARMA_IC_selection_CSS(15,0,EU_ILS_5Y5Y) 
AR_EU_ILS_5Y5Y <- arima(EU_ILS_5Y5Y, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_EU_ILS_5Y5Y_residuals <- residuals(AR_EU_ILS_5Y5Y)
Qstat_AR_EU_ILS_5Y5Y_residuals <- LjungBox(AR_EU_ILS_5Y5Y_residuals,20,2)
stargazer(Qstat_AR_EU_ILS_5Y5Y_residuals,type = "text", summary = FALSE)

adf_EU_ILS_5Y5Y <- ur.df(EU_ILS_5Y5Y, type = "trend", lags = 1)
summary(adf_EU_ILS_5Y5Y)
adf_EU_ILS_5Y5Y <- ur.df(EU_ILS_5Y5Y, type = "drift", lags = 1)
summary(adf_EU_ILS_5Y5Y)
adf_EU_ILS_5Y5Y <- ur.df(EU_ILS_5Y5Y, type = "none", lags = 1)
summary(adf_EU_ILS_5Y5Y)

bp_EU_ILS_5Y5Y <- breakpoints(EU_ILS_5Y5Y ~ 1, h = 0.15)
summary(bp_EU_ILS_5Y5Y) 
breakdates(bp_EU_ILS_5Y5Y)
plot(EU_ILS_5Y5Y, main = "Bai-Perron Structural Breaks: EU ILS 5Y5Y")
lines(bp_EU_ILS_5Y5Y)

za_EU_ILS_5Y5Y <- ur.za(EU_ILS_5Y5Y, model = "both", lag = 2)
summary(za_EU_ILS_5Y5Y)
plot(za_EU_ILS_5Y5Y)

ls_EU_ILS_5Y5Y <- ur.ls(y = EU_ILS_5Y5Y, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

#-----
# EPU
#-----
ARMA_IC_selection_CSS(15,0,EPU)
AR_EPU <- arima(EPU, order = c(4, 0, 0),
                method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_EPU_residuals <- residuals(AR_EPU)
Qstat_AR_EPU_residuals <- LjungBox(AR_EPU_residuals,20,4)
stargazer(Qstat_AR_EPU_residuals,type = "text", summary = FALSE)

adf_EPU <- ur.df(EPU, type = "trend", lags =3)
summary(adf_EPU)

bp_EPU <- breakpoints(EPU ~ 1, h = 0.15)
summary(bp_EPU) 
breakdates(bp_EPU)
plot(EPU, main = "Bai-Perron Structural Breaks: EPU")
lines(bp_EPU)

za_EPU <- ur.za(EPU, model = "both", lag = 2)
summary(za_EPU)
plot(za_EPU)

ls_EPU <- ur.ls(y = EPU, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, print.results = "print" )

# ==========================================
# AUTOMATED TABLE GENERATION (xtable)
# ==========================================

# Helper function to cleanly extract Bai-Perron dates
get_bp_dates <- function(bp_obj) {
  # If the object is missing or wasn't run, return "Not Run" to prevent a crash
  if (!exists(deparse(substitute(bp_obj)))) return("Not Run")
  
  dates <- breakdates(bp_obj)
  if (is.na(dates[1])) return("None")
  return(paste(round(dates, 2), collapse = ", "))
}

# Initialize empty DataFrame (Lag columns removed)
results_table <- data.frame(
  Variable = character(),
  ADF_Stat = numeric(),
  LS_Stat = numeric(),
  BP_Breaks = numeric(),
  BP_Dates = character(),
  stringsAsFactors = FALSE
)

# Helper function to dynamically add rows (adf_lag parameter removed)
add_row <- function(var_name, adf_obj, ls_obj, bp_obj) {
  bp_dates <- get_bp_dates(bp_obj)
  bp_breaks <- length(breakdates(bp_obj))
  if(is.na(breakdates(bp_obj)[1])) bp_breaks <- 0
  
  new_row <- data.frame(
    Variable = var_name,
    ADF_Stat = round(adf_obj@teststat[1], 2),
    LS_Stat = ifelse(!is.null(ls_obj[[1]]$`t-stat`), round(ls_obj[[1]]$`t-stat`, 2), NA),
    BP_Breaks = bp_breaks,
    BP_Dates = bp_dates,
    stringsAsFactors = FALSE
  )
  return(new_row)
}

# Populate table with the assigned objects (manual lag inputs removed)
results_table <- rbind(results_table, add_row("HICP", adf_HICP, ls_HICP, bp_HICP))
results_table <- rbind(results_table, add_row("Core HICP", adf_core_HICP, ls_core_HICP, bp_core_HICP))
results_table <- rbind(results_table, add_row("PCCI", adf_PCCI, ls_PCCI, bp_PCCI))
results_table <- rbind(results_table, add_row("Ind. Production", adf_IP, ls_IP, bp_IP))
results_table <- rbind(results_table, add_row("Unemployment", adf_Unemployment_Rate, ls_Unemployment_Rate, bp_Unemployment_Rate))
results_table <- rbind(results_table, add_row("ESI", adf_ESI, ls_ESI, bp_ESI))
results_table <- rbind(results_table, add_row("Oil Prices", adf_Oil_prices, ls_Oil_prices, bp_Oil_prices))
results_table <- rbind(results_table, add_row("EONIA 1Y", adf_EONIA_1Y, ls_EONIA_1Y, bp_EONIA_1Y))
results_table <- rbind(results_table, add_row("True Policy Rate", adf_True_Policy_Rate, ls_True_Policy_Rate, bp_True_Policy_Rate))
results_table <- rbind(results_table, add_row("GPR", adf_GPR, ls_GPR, bp_GPR))
results_table <- rbind(results_table, add_row("VIX", adf_VIX, ls_VIX, bp_VIX))
results_table <- rbind(results_table, add_row("VSTOXX", adf_VSTOXX, ls_VSTOXX, bp_VSTOXX))
results_table <- rbind(results_table, add_row("EU ILS 1Y", adf_EU_ILS_1Y, ls_EU_ILS_1Y, bp_EU_ILS_1Y))
results_table <- rbind(results_table, add_row("EU ILS 5Y5Y", adf_EU_ILS_5Y5Y, ls_EU_ILS_5Y5Y, bp_EU_ILS_5Y5Y))
results_table <- rbind(results_table, add_row("EPU", adf_EPU, ls_EPU, bp_EPU))

# ==========================================
# SPLIT INTO TWO SEPARATE LATEX TABLES
# ==========================================

# 1. Unit Root Table (ADF and LS)
ur_table <- results_table[, c("Variable", "ADF_Stat", "LS_Stat")]
# Clean up the column names for LaTeX
colnames(ur_table) <- c("Variable", "ADF Statistic", "LS Statistic") 

latex_ur <- xtable(
  ur_table, 
  caption = "Unit Root Test Results", 
  label = "tab:unit_roots",
  align = "llcc" # 1 hidden rowname col + 3 data cols
)

# 2. Bai-Perron Table
bp_table <- results_table[, c("Variable", "BP_Breaks", "BP_Dates")]
# Clean up the column names for LaTeX
colnames(bp_table) <- c("Variable", "Opt. Breaks (BIC)", "Identified Break Dates")

latex_bp <- xtable(
  bp_table, 
  caption = "Bai-Perron Multiple Structural Break Analysis", 
  label = "tab:bai_perron",
  align = "llcl" # 'l' for left-aligned dates
)

# ==========================================
# PRINT BOTH TO CONSOLE
# ==========================================

cat("\n% --- TABLE 1: UNIT ROOT TESTS ---\n")
print(latex_ur, include.rownames = FALSE)

cat("\n% --- TABLE 2: BAI-PERRON BREAKS ---\n")
print(latex_bp, include.rownames = FALSE)