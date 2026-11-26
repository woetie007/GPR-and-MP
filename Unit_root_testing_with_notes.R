###############################################################
# UNIT ROOT TESTING
###############################################################

#=========
#Packages
#=========
library(dplyr)
library(lubridate)
library(readxl)
library(countrycode)
library(xts)
library(tidyr)
library(plotly)
library(lpirfs)
library(AER)
library(ggplot2)
library(vars)
library(sandwich)
library(conflicted)
library(car)
library(stargazer)
library(strucchange)
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
source("Merge.R")
source("TS_functions.R")
source("TVAR_functions.R")
source("joint grid search.R")
source("additional_functions.R")
source("LeeStrazicichUnitRootTest.R")
#=====================
# Prepare & plot Data
#=====================

unit_root_data <- final_data %>% filter(Country == "EA")%>% 
  select(Date,df_HICP, Log_IP, df_Oil, True_Policy_Rate, Log_GPR, Log_GPRA, Log_GPRT,
         Log_HICP,Log_Core_HICP, Log_Oil, Spliced_EONIA_3M, Log_VSTOXX, Log_EPU,
         df_Core_HICP, df_True_Policy_Rate, Spliced_EONIA_1Y, Unemployment_Rate
         , ESI, df_Unemployment_Rate, OIS_3M, OIS_1Y, STOXX50, Log_VIX, Target_Shock, FG_Shock, QE_Shock)%>% # Adding PCCI gets rid of 14 observations
  drop_na()

# Convert each variable into a time series object starting Feb 2000
HICP                <- ts(unit_root_data$Log_HICP, start = c(2000, 2), frequency = 12)
core_HICP           <- ts(unit_root_data$Log_Core_HICP, start = c(2000, 2), frequency = 12)
IP                  <- ts(unit_root_data$Log_IP, start = c(2000, 2), frequency = 12)
Oil_prices          <- ts(unit_root_data$Log_Oil, start = c(2000, 2), frequency = 12)
EONIA_1Y            <- ts(unit_root_data$Spliced_EONIA_1Y, start = c(2000, 2), frequency = 12)
GPR                 <- ts(unit_root_data$Log_GPR, start = c(2000, 2), frequency = 12)
Unemployment_Rate   <- ts(unit_root_data$Unemployment_Rate, start = c(2000, 2), frequency = 12)
ESI                 <- ts(unit_root_data$ESI, start = c(2000, 2), frequency = 12)
True_Policy_Rate    <- ts(unit_root_data$True_Policy_Rate, start = c(2000, 2), frequency = 12)
PCCI               <- ts(unit_root_data$PCCI, start = c(2001, 4), frequency = 12)
df_HICP             <- ts(unit_root_data$df_HICP, start = c(2000, 2), frequency = 12)
df_Core_HICP        <- ts(unit_root_data$df_Core_HICP, start = c(2000, 2), frequency = 12)
df_True_Policy_Rate <- ts(unit_root_data$df_True_Policy_Rate, start = c(2000, 2), frequency = 12)
OIS_3M              <- ts(unit_root_data$OIS_3M, start = c(2000, 2), frequency = 12)
OIS_1Y              <- ts(unit_root_data$OIS_1Y, start = c(2000, 2), frequency = 12)
STOXX50             <- ts(unit_root_data$STOXX50, start = c(2000, 2), frequency = 12)
Target_Shock             <- ts(unit_root_data$Target_Shock, start = c(2000, 2), frequency = 12)
FG_Shock             <- ts(unit_root_data$FG_Shock, start = c(2000, 2), frequency = 12)
QE_Shock             <- ts(unit_root_data$QE_Shock, start = c(2000, 2), frequency = 12)
STOXX50             <- ts(unit_root_data$QE_Shock, start = c(2000, 2), frequency = 12)
VIX                 <- ts(unit_root_data$Log_VIX, start = c(2000, 2), frequency = 12)
VSTOXX              <- ts(unit_root_data$Log_VSTOXX, start = c(2000, 2), frequency = 12)
EU_ILS_1Y           <- ts(unit_root_data$EU_ILS_1Y, start = c(2008, 8), frequency = 12)
EU_ILS_5Y5Y         <- ts(unit_root_data$EU_ILS_5Y5Y, start = c(2008, 8), frequency = 12)
EPU                 <- ts(unit_root_data$Log_EPU, start = c(2000, 2), frequency = 12)

plot.ts(HICP, main = "HICP", ylab = "Log HICP", xlab = "Time")
plot.ts(core_HICP, main = "Core HICP", ylab = "Log Core HICP", xlab = "Time")
#plot.ts(PCCI, main="PCCI", ylab="PCCI", xlab="Time")
plot.ts(IP, main = "Industrial Production", ylab = "Log IP", xlab = "Time")
plot.ts(ESI, main="Economic Sentiment Indicator", ylab= "ESI", xlab="Time")
plot.ts(Unemployment_Rate, main="Unemployment Rate", ylab="Unemployment Rate", xlab="Time")
plot.ts(Oil_prices, main = "Oil Prices", ylab = "Log Oil Prices", xlab = "Time")
plot.ts(True_Policy_Rate, main="True Policy Rate", ylab="True Policy Rate", xlab="Time")
plot.ts(EONIA_1Y, main = "EONIA 1Y", ylab = "EONIA 1Y", xlab = "Time")
plot.ts(GPR, main = "Geopolitical Risk", ylab = "Log GPR", xlab = "Time")
plot.ts(EPU, main = "Economic Policy Uncertainty", ylab = "Log EPU", xlab = "Time")
plot.ts(VSTOXX, main = "VSTOXX", ylab = "VSTOXX", xlab = "Time")
plot.ts(VIX, main = "VIX", ylab = "VIX", xlab = "Time")
plot.ts(df_HICP, main = "df HICP", ylab = "df HICP", xlab = "Time")
plot.ts(df_Core_HICP, main = "df Core HICP", ylab = "df Core HICP", xlab = "Time")
plot.ts(df_True_Policy_Rate, main = "df True_Policy_Rate", ylab = "df True Policy Rate", xlab = "Time")
plot.ts(STOXX50, main = "STOXX50", ylab = "STOXX50 surprise", xlab = "Time")
plot.ts(OIS_1Y, main = "OIS 1Y surprise", ylab = "OIS 1Y surprise", xlab = "Time")
plot.ts(OIS_3M, main = "OIS_3M surprise", ylab = "OIS 3M surprise", xlab = "Time")
plot.ts(EU_ILS_1Y, main = "EU ILS 1Y", ylab = "EU ILS 1Y", xlab = "Time")
plot.ts(EU_ILS_5Y5Y, main = "EU ILS 5Y5Y", ylab = "EU ILS 5Y5Y", xlab = "Time")






#==================
# Unit Root Testing
#==================
#-----
#HICP
#-----
ARMA_IC_selection_CSS(15,0,HICP) #BIC suggests AR(5)
AR_HICP <- arima(HICP, order = c(6, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_HICP_residuals <- residuals(AR_HICP)
## Ljung-Box
Qstat_AR_HICP_residuals <- LjungBox(AR_HICP_residuals,20,6)
stargazer(Qstat_AR_HICP_residuals,type = "text", summary = FALSE) # with AR(5) still autocorrelation


#ADF(5) (p-1)
adf <- ur.df(HICP, type = "trend", lags =5)
summary(adf)
# gamma=0? yes
# gamma= a2=0? yes, kick out trend
adf <- ur.df(HICP, type = "drift", lags =5)
summary(adf)
# gamma=0? yes
# gamma= a0=0? No
# gamma =0 with normal dist? Yes, series is a random walk drift

# Bai perron test
bp_test <- breakpoints(HICP ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # BIC indicates 5 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(HICP, main = "Bai-Perron Structural Breaks: HICP")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(HICP, model = "both", lag = 2)

# View the full results
summary(za_test) # unit root with structural break at t=263

# Plot the test to visually see where the algorithm identified the break
plot(za_test) 

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = HICP, model = "break", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" ) 
# trend-stationary around breaks at t=142 and 248,p=13 but at 10%

#----------
# Core HICP
#----------
ARMA_IC_selection_CSS(15,0,core_HICP) #BIC suggest AR(7)
AR_core_HICP <- arima(core_HICP, order = c(7, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_core_HICP_residuals <- residuals(AR_core_HICP)
## Ljung-Box
Qstat_AR_core_HICP_residuals <- LjungBox(AR_core_HICP_residuals, 20, 7)
stargazer(Qstat_AR_core_HICP_residuals,type = "text", summary = FALSE) #No serial correlation

#ADF(6)(p-1)
adf <- ur.df(core_HICP, type = "trend", lags =6)
summary(adf)
# gamma=0? yes
# gamma= a2=0? yes, kick out trend
adf <- ur.df(core_HICP, type = "drift", lags =6)
summary(adf)
# gamma=0? yes
# gamma= a0=0? yes, kick out drift
adf <- ur.df(core_HICP, type = "none", lags =6)
summary(adf)
# gamma=0? no, conclude no unit root

#Bai-perron test
bp_test <- breakpoints(core_HICP ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # BIC indicates 5 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(core_HICP, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(core_HICP, model = "both", lag = 2)

# View the full results
summary(za_test) #unit root with break at t=231

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = core_HICP, model = "break", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )
# p= 12, trend stationary with breaks at t= 142 and 245
#------------
# PCCI
#------------

ARMA_IC_selection_CSS(15,0,PCCI) #BIC suggests p=11
AR_PCCI <- arima(PCCI, order=c(11,0,0),
                 method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_PCCI_residuals <- residuals(AR_PCCI)
## Ljung-Box
Qstat_AR_PCCI_residuals <- LjungBox(AR_PCCI_residuals,20,11)
stargazer(Qstat_AR_PCCI_residuals,type = "text", summary = FALSE)# no serial correlation

#ADF(10)(p-1)
adf <- ur.df(PCCI, type = "trend", lags =10)
summary(adf)
# gamma=0? No, conclude no unit root

#Bai-perron test
bp_test <- breakpoints(PCCI ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test)

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(PCCI, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(PCCI, model = "both", lag = 2)

# View the full results
summary(za_test)

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = PCCI, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )

#-----
# IP
#-----
ARMA_IC_selection_CSS(15,0,IP) #BIC and AIC suggests AR(3)
AR_IP <- arima(IP, order = c(3, 0, 0),
 method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_IP_residuals <- residuals(AR_IP)
## Ljung-Box
Qstat_AR_IP_residuals <- LjungBox(AR_IP_residuals,20,3)
stargazer(Qstat_AR_IP_residuals,type = "text", summary = FALSE)# no serial correlation

#ADF(2)(p-1)
adf <- ur.df(IP, type = "trend", lags =2)
summary(adf)
#gamma=0? no, conclude no unit root

#Bai-perron test
bp_test <- breakpoints(IP ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # BIC indicates 2 or 4 structural breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(IP, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(IP, model = "both", lag = 2)

# View the full results
summary(za_test) # no unit root, break at t=104, GFC

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = IP, model = "break", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )
# p =2, unit root with breaks at t= 103 and 122

#------------------
# Unemployment Rate
#------------------

ARMA_IC_selection_CSS(15,0,Unemployment_Rate) #BIC suggests AR(3), AIC suggests AR(10)
AR_Unemployment_Rate <- arima(Unemployment_Rate, order = c(8, 0, 0),
               method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_Unemployment_Rate_residuals <- residuals(AR_Unemployment_Rate)
## Ljung-Box
Qstat_AR_Unemployment_Rate_residuals <- LjungBox(AR_Unemployment_Rate_residuals,20,8)
stargazer(Qstat_AR_Unemployment_Rate_residuals,type = "text", summary = FALSE)# No serial correlation at p=8 


#ADF(7)(p-1)
adf <- ur.df(Unemployment_Rate, type = "trend", lags =7)
summary(adf)
# gamma=0? Yes
# gamma=a2=0? yes, drop trend
adf <- ur.df(Unemployment_Rate, type = "drift", lags =7)
summary(adf)
#gamma=0? yes
#gamma=a0=0? yes
adf <- ur.df(Unemployment_Rate, type = "none", lags =7)
summary(adf)
#gamma=0? yes, conclude unit root
#Bai-perron test
bp_test <- breakpoints(Unemployment_Rate ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 5 structural breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(Unemployment_Rate, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(Unemployment_Rate, model = "both", lag = 2)

# View the full results
summary(za_test) # unit root, break at t=102, GFC

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = Unemployment_Rate, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )
# p=8; unit root with breaks at t=109 and 187

#-------
# ESI
#-------

ARMA_IC_selection_CSS(15,0,ESI) #BIC suggests AR(2), AIC suggests AR(10)
AR_ESI <- arima(ESI, order = c(2, 0, 0),
               method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_ESI_residuals <- residuals(AR_ESI)
## Ljung-Box
Qstat_AR_ESI_residuals <- LjungBox(AR_ESI_residuals,20,2)
stargazer(Qstat_AR_ESI_residuals,type = "text", summary = FALSE)# no serial correlation

#ADF(1)(p-1)
adf <- ur.df(ESI, type = "trend", lags =1)
summary(adf)
#gamma=0? no, conclude no unit root
#Bai-perron test
bp_test <- breakpoints(ESI ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 3 structural breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(ESI, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = ESI, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )
# p=9, unit root with breaks at t=163 and t=266

#-----------
# Oil Prices
#-----------
ARMA_IC_selection_CSS(15,0,Oil_prices) #BIC suggests AR(2), AIC suggests AR(3)
AR_oil <- arima(Oil_prices, order = c(2, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_oil_residuals <- residuals(AR_oil)
## Ljung-Box
Qstat_AR_oil_residuals <- LjungBox(AR_oil_residuals,20,2)
stargazer(Qstat_AR_oil_residuals,type = "text", summary = FALSE) #No serial correlation

##ADF(1)(p-1)
adf <- ur.df(Oil_prices, type = "trend", lags =1)
summary(adf)
#gamma=0? yes
# gamma=a2=0? yes, kick out trend
adf <- ur.df(Oil_prices, type = "drift", lags =1)
summary(adf)
#gamma=0? yes
# gamma=a0=0? yes but significant at 10%, random walk with drift? 
adf <- ur.df(Oil_prices, type = "none", lags =1)
summary(adf)
#gamma=0? yes, conclude unit root
#Bai-perron test
bp_test <- breakpoints(Oil_prices ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 structural breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(Oil_prices, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = Oil_prices, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )
# unit root, structural breaks at t=178 and 190, p=1
#------------
# EONIA_1Y
#------------

ARMA_IC_selection_CSS(15,0,EONIA_1Y) #BIC suggests AR(3), AIC suggests AR(5)
AR_IR <- arima(EONIA_1Y, order = c(3, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_IR_residuals <- residuals(AR_IR)
## Ljung-Box
Qstat_AR_IR_residuals <- LjungBox(AR_IR_residuals,20,3)
stargazer(Qstat_AR_IR_residuals,type = "text", summary = FALSE) #No serial correlation

#ADF(2)(p-1)
adf <- ur.df(EONIA_1Y, type = "trend", lags =2)
summary(adf)
#gamma=0? yes
# gamma=a2=0? yes, kick out trend
adf <- ur.df(EONIA_1Y, type = "drift", lags =2)
summary(adf)
#gamma=0? yes
# gamma=a0=0? yes, kick out drift
adf <- ur.df(EONIA_1Y, type = "none", lags =2)
summary(adf)
#gamma=0? no, conclude no unit root at 10%
#Bai-perron test
bp_test <- breakpoints(EONIA_1Y ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 3 structural breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(EONIA_1Y, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(EONIA_1Y, model = "both", lag = 2)

# View the full results
summary(za_test)

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = EONIA_1Y, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )

#-----------------
# True Policy rate
#-----------------

ARMA_IC_selection_CSS(15,0,True_Policy_Rate) #BIC suggests AR(4), AIC suggests AR(9)
AR_True_Policy_Rate <- arima(True_Policy_Rate, order = c(4, 0, 0),
                method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_True_Policy_Rate_residuals <- residuals(AR_True_Policy_Rate)
## Ljung-Box
Qstat_AR_True_Policy_Rate_residuals <- LjungBox(AR_True_Policy_Rate_residuals,20,4)
stargazer(Qstat_AR_True_Policy_Rate_residuals,type = "text", summary = FALSE)# still serial correlation

# AIC suggestion
AR_True_Policy_Rate <- arima(True_Policy_Rate, order = c(9, 0, 0),
                             method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_True_Policy_Rate_residuals <- residuals(AR_True_Policy_Rate)
## Ljung-Box
Qstat_AR_True_Policy_Rate_residuals <- LjungBox(AR_True_Policy_Rate_residuals,20,9)
stargazer(Qstat_AR_True_Policy_Rate_residuals,type = "text", summary = FALSE)# No autocorrelation

#ADF(8)(p-1)
adf <- ur.df(True_Policy_Rate, type = "trend", lags =8)
summary(adf)
#gamma=0? yes
# gamma=a2=0? yes, kick out trend
adf <- ur.df(True_Policy_Rate, type = "drift", lags =8)
summary(adf)
# gamma=0? yes
#gamma=a0=0? yes, kick out drift
adf <- ur.df(True_Policy_Rate, type = "none", lags =8)
summary(adf)# gamma = 0? No, conclude no unit rot

#Bai-perron test
bp_test <- breakpoints(True_Policy_Rate ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(True_Policy_Rate, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(True_Policy_Rate, model = "both", lag = 2)

# View the full results
summary(za_test)

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = True_Policy_Rate, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )

#------
# GPR
#------
ARMA_IC_selection_CSS(15,0,GPR) #BIC suggests AR(1), AIC suggests AR(1)
AR_GPR <- arima(GPR, order = c(1, 0, 0), method = "CSS", include.mean = TRUE, transform.pars = FALSE) #BIC
AR_GPR_residuals <- residuals(AR_GPR)
## Ljung-Box
Qstat_AR_GPR_residuals <- LjungBox(AR_GPR_residuals,20,1)
stargazer(Qstat_AR_GPR_residuals,type = "text", summary = FALSE) #No serial correlation

#ADF(0)(p-1)
adf <- ur.df(GPR, type = "trend", lags =0)
summary(adf)
#gamma=0? no, conclude no unit root

#Bai-perron test
bp_test <- breakpoints(GPR ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 2 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(GPR, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Zivot Andrews test
za_test <- ur.za(GPR, model = "both", lag = 2)

# View the full results
summary(za_test) # no unit root, break at t=38

# Plot the test to visually see where the algorithm identified the break
plot(za_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = GPR, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" ) # no unit root, break at t=81 and 216, p=NA


#------
# VIX
#------
ARMA_IC_selection_CSS(15,0,VIX) #BIC suggests AR(1)
AR_VIX <- arima(VIX, order = c(1, 0, 0),
                             method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_VIX_residuals <- residuals(AR_VIX)
## Ljung-Box
Qstat_AR_VIX_residuals <- LjungBox(AR_VIX_residuals,20,1)
stargazer(Qstat_AR_VIX_residuals,type = "text", summary = FALSE)

#ADF(0)(p-1)
adf <- ur.df(VIX, type = "trend", lags =0)
summary(adf)
#gamma=0? no, conclude no unit root

#Bai-perron test
bp_test <- breakpoints(VIX ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(VIX, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = VIX, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )

#--------
# VSTOXX
#--------

ARMA_IC_selection_CSS(15,0,VSTOXX) #BIC suggests AR(2)
AR_VSTOXX <- arima(VSTOXX, order = c(2, 0, 0),
                             method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_VSTOXX_residuals <- residuals(AR_VSTOXX)
## Ljung-Box
Qstat_AR_VSTOXX_residuals <- LjungBox(AR_VSTOXX_residuals,20,2)
stargazer(Qstat_AR_VSTOXX_residuals,type = "text", summary = FALSE)

#ADF(1)(p-1)
adf <- ur.df(VSTOXX, type = "trend", lags =1)
summary(adf)
#gamma=0? No, conclude no unit rot

#Bai-perron test
bp_test <- breakpoints(VSTOXX ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(VSTOXX, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = VSTOXX, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )



#-----------
# EU_ILS_1Y
#-----------
ARMA_IC_selection_CSS(15,0,EU_ILS_1Y) #BIC suggests AR(3)
AR_EU_ILS_1Y <- arima(EU_ILS_1Y, order = c(3, 0, 0),
                   method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_EU_ILS_1Y_residuals <- residuals(AR_EU_ILS_1Y)
## Ljung-Box
Qstat_AR_EU_ILS_1Y_residuals <- LjungBox(AR_EU_ILS_1Y_residuals,20,3)
stargazer(Qstat_AR_EU_ILS_1Y_residuals,type = "text", summary = FALSE)

#ADF(2)(p-1)
adf <- ur.df(EU_ILS_1Y, type = "trend", lags =2)
summary(adf)
#gamma=0? yes
#gamma =a2=0? yes, kick out trend
adf <- ur.df(EU_ILS_1Y, type = "drift", lags =2)
summary(adf)
#gamma=0? yes
#gamma=a0=0?Yes, kick out drift
adf <- ur.df(EU_ILS_1Y, type = "none", lags =2)
summary(adf)
#gamma=0? yes, conclude unit root

#Bai-perron test
bp_test <- breakpoints(EU_ILS_1Y ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(EU_ILS_1Y, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = EU_ILS_1Y, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )


#-------------
# EU_ILS_5Y5Y
#-------------
ARMA_IC_selection_CSS(15,0,EU_ILS_5Y5Y) #BIC suggests AR(2)
AR_EU_ILS_5Y5Y <- arima(EU_ILS_5Y5Y, order = c(2, 0, 0),
                      method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_EU_ILS_5Y5Y_residuals <- residuals(AR_EU_ILS_5Y5Y)
## Ljung-Box
Qstat_AR_EU_ILS_5Y5Y_residuals <- LjungBox(AR_EU_ILS_5Y5Y_residuals,20,2)
stargazer(Qstat_AR_EU_ILS_5Y5Y_residuals,type = "text", summary = FALSE)

#ADF(1)(p-1)
adf <- ur.df(EU_ILS_5Y5Y, type = "trend", lags =1)
summary(adf)
#gamma=0? yes
# gamma=a2=0? yes, kick out trend
adf <- ur.df(EU_ILS_5Y5Y, type = "drift", lags =1)
summary(adf)
#gamma =0? yes
#gamma=a0=0? yes, kick out drift
adf <- ur.df(EU_ILS_5Y5Y, type = "none", lags =1)
summary(adf)
#gamma=0?yes, conclude unit root 

#Bai-perron test
bp_test <- breakpoints(EU_ILS_5Y5Y ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(EU_ILS_5Y5Y, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = EU_ILS_5Y5Y, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )


#-----
# EPU
#-----

ARMA_IC_selection_CSS(15,0,EPU) #BIC suggests AR(4)
AR_EPU <- arima(EPU, order = c(4, 0, 0),
                        method = "CSS", include.mean = TRUE, transform.pars = FALSE)
AR_EPU_residuals <- residuals(AR_EPU)
## Ljung-Box
Qstat_AR_EPU_residuals <- LjungBox(AR_EPU_residuals,20,4)
stargazer(Qstat_AR_EPU_residuals,type = "text", summary = FALSE)

#ADF(3)(p-1)
adf <- ur.df(EPU, type = "trend", lags =3)
summary(adf)
#gamma=0? No , conclude no unit root 

#Bai-perron test
bp_test <- breakpoints(EPU ~ 1, h = 0.15)

# View the optimal number of breaks and their specific dates (selected via BIC)
summary(bp_test) # 4 breaks

# To extract the exact break dates:
breakdates(bp_test)

# Plot the original series with the Bai-Perron break dates overlaid
plot(EPU, main = "Bai-Perron Structural Breaks: Industrial Production")
lines(bp_test)

# Two-Break Lee-Strazicich Test
ls_test <- ur.ls(y = EPU, model = "crash", breaks = 2, method = "GTOS", pn = 0.1, , print.results = "print" )
