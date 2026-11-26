#################################################
# Master's Dissertation, Wout Baetens, r0838656
#################################################
#Packages
library(eurostat)
library(dplyr)
library(tidyr)
library(conflicted)
library(ecb)
library(lubridate)
library(readxl)
library(countrycode)
library(xts)
library(plotly)
library(lpirfs)
library(AER)
library(ggplot2)
library(vars)
library(sandwich)
library(car)
library(urca)
library(strucchange)
library(xtable)
library(stargazer)
library(FinTS)
library(dynlm)
library(readr)
library(writexl)
library(MASS)
library(quantmod)
library(moments)
library(knitr)
library(kableExtra)

# gathering all the different datasets
source("GPR data.R")
source("Euro Area data.R")
source("monetary policy shocks calculation.R")
source("ECB policy rates.R")
source("FRED data.R")
source("jk_2020_shocks.R")

#Merging the datasets

final_data <- long_macro_data %>%
  left_join(clean_gpr_long, by = c("Country", "Date")) %>%
  left_join(EPU, by = c("Date")) %>%
  left_join(monthly_shocks, by = "Date") %>%
  left_join(policy_rates, by = "Date") %>%
  left_join(oil_data, by = "Date")%>%
  left_join(vix_clean, by = "Date") %>%
  left_join(df_monthly_shocks, by="Date")
final_data <- final_data %>%
  filter(Date >= as.Date("1999-01-01")) %>%
  filter(Date <= as.Date("2025-12-01")) %>%
  arrange(Country, Date)

summary(final_data$Date)

colSums(is.na(final_data))
final_data %>% count(Country)





