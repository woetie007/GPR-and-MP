
conflicts_prefer(ecb::get_data)
conflicts_prefer(dplyr::filter)
# DOWNLOAD RAW DATA 

# MRO (Daily)
mro_raw <- get_data("FM.D.U2.EUR.4F.KR.MRR_RT.LEV") %>%
  mutate(Date = as.Date(obstime)) %>%
  select(Date, MRO = obsvalue)

# DFR (Daily)
deposit_raw <- get_data("FM.D.U2.EUR.4F.KR.DFR.LEV") %>%
  mutate(Date = as.Date(obstime)) %>%
  select(Date, Deposit_Rate = obsvalue)

# EONIA (Monthly)
eonia_raw <- get_data("FM.M.U2.EUR.4F.MM.EONIA.HSTA") %>%
  mutate(Date = paste0(obstime, "-01"),
         Date = as.Date(Date)) %>%
  select(Date, EONIA = obsvalue)

# €STR (Daily)
estr_raw <- get_data("EST.B.EU000A2X2A25.WT") %>%
  mutate(Date = as.Date(obstime)) %>%
  select(Date, ESTR = obsvalue)

# AGGREGATE TO MONTHLY

# MRO & DFR
daily_policy_rates <- mro_raw %>%
  full_join(deposit_raw, by = "Date") %>%
  filter(Date >= as.Date("2000-01-01"))

monthly_policy_rates <- daily_policy_rates %>%
  mutate(year_month = floor_date(Date, "month")) %>%
  group_by(year_month) %>%
  summarise(
    MRO = mean(MRO, na.rm = TRUE),
    Deposit_Rate = mean(Deposit_Rate, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  rename(Date = year_month)

# €STR
estr_monthly <- estr_raw %>%
  mutate(Date = floor_date(Date, "month")) %>%
  group_by(Date) %>%
  summarise(ESTR = mean(ESTR, na.rm = TRUE)) %>%
  ungroup()


# Clean Shadow Rate

shadow_rates <- read_excel("shadowrate_ECB.xls")
clean_shadow <- shadow_rates %>%
  mutate(
    Date_Text = as.character(Date), 
    Date_Full = paste0(Date_Text, "01"), 
    Date = ymd(Date_Full) 
  ) %>%
  select(Date, Shadow_Rate)

# Inflation Swap rates and EONIA/€STR maturities 

ILS_EONIA <- read_excel("Inflation swap and euro rates.xlsx")
ILS_EONIA <- ILS_EONIA %>%
  mutate(Date = floor_date(Date, "month"),
         ILS_10Y = EU_ILS_10Y/100,
         ILS_5Y = EU_ILS_5Y/100,
         EU_ILS_5Y5Y= (((((1+ILS_10Y)^10)/((1+ILS_5Y)^5)))^(1/5)-1)*100,
         df_ILS_5Y5Y= EU_ILS_5Y5Y - lag(EU_ILS_5Y5Y),
         df_ILS_1Y= EU_ILS_1Y- lag(EU_ILS_1Y))

# Merge and Splice 

policy_rates <- monthly_policy_rates %>%
  full_join(eonia_raw, by = "Date") %>%
  full_join(estr_monthly, by = "Date") %>% # Merge ESTR here
  filter(Date >= as.Date("2000-01-01") & Date <= as.Date("2025-12-01")) %>%
  arrange(Date) %>%
  left_join(clean_shadow, by = "Date") %>%
  left_join(ILS_EONIA, by = "Date")

# Build the final True Policy Rate
policy_rates <- policy_rates %>% 
  mutate(
    Log_VSTOXX= log(VSTOXX_Volatility_Index)*100,
    # Create the continuous EONIA/€STR series
    Spliced_EONIA_3M = case_when(
      !is.na(EONIA) ~ EONIA,
      is.na(EONIA) & !is.na(ESTR) ~ ESTR + 0.085,
      TRUE ~ NA_real_
    ),
    Spliced_EONIA_1Y = case_when(
      !is.na(EONIA_1Y) ~ EONIA_1Y,
      is.na(EONIA_1Y) & !is.na(ESTR_1Y) ~ ESTR_1Y + 0.085,
      TRUE ~ NA_real_
    ) ,
    Spliced_EONIA_2Y = case_when(
      !is.na(EONIA_2Y) ~ EONIA_2Y,
      is.na(EONIA_2Y) & !is.na(ESTR_2Y) ~ ESTR_2Y + 0.085,
      TRUE ~ NA_real_
    ),
    Spliced_EONIA_5Y = case_when(
      !is.na(EONIA_5Y) ~ EONIA_5Y,
      is.na(EONIA_5Y) & !is.na(ESTR_5Y) ~ ESTR_5Y + 0.085,
      TRUE ~ NA_real_
    ),
    Spliced_EONIA_10Y = case_when(
      !is.na(EONIA_10Y) ~ EONIA_10Y,
      is.na(EONIA_10Y) & !is.na(ESTR_10Y) ~ ESTR_10Y + 0.085,
      TRUE ~ NA_real_
    ) ,
    # Calculate the True Policy Rate using the Spliced series and the Shadow Rate
    True_Policy_Rate = ifelse(is.na(Shadow_Rate), Spliced_EONIA_3M, pmin(Spliced_EONIA_3M, Shadow_Rate)),
    df_True_Policy_Rate = True_Policy_Rate - lag(True_Policy_Rate, 1),
    df_EONIA_1Y = Spliced_EONIA_1Y-lag(Spliced_EONIA_1Y)
  ) %>%
  select(Date, EONIA, ESTR, Spliced_EONIA_3M,Spliced_EONIA_1Y, Spliced_EONIA_2Y, df_EONIA_1Y,
         Spliced_EONIA_5Y, Spliced_EONIA_10Y, MRO, Deposit_Rate, Shadow_Rate, df_ILS_5Y5Y, df_ILS_1Y,
         True_Policy_Rate, df_True_Policy_Rate, EU_ILS_2Y =`EU_ILS_ 2Y`, EU_ILS_5Y, EU_ILS_10Y,
         EU_ILS_1Y, EU_ILS_5Y5Y, Log_VSTOXX, VSTOXX = VSTOXX_Volatility_Index)
