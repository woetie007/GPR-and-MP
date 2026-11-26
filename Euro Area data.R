conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

########################
# Industrial Production
########################

# Download the raw database for Industrial Production
raw_ip <- get_eurostat("sts_inpr_m", time_format = "date", cache = FALSE)

# Filter for the exact macroeconomic definition you need
clean_ip <- raw_ip %>%
  filter(
    indic_bt == "PRD",        # Production volume
    nace_r2 == "B-D",         # Industry (excluding construction)
    s_adj == "SCA",           # Seasonally & Calendar Adjusted
    unit == "I21",            # Base year 2021 = 100
    # Select enitities
    geo %in% c("EA21") 
  ) %>%
  select(geo, TIME_PERIOD, values) %>%
  rename(
    Country = geo, 
    Date = TIME_PERIOD, 
    IP_Index = values
  )%>%
  mutate(Country = if_else(Country == "EA21", "EA", Country))

# View the final, perfectly formatted dataset
head(clean_ip)

#######################
# Unemployment Rate
#######################

raw_une <- get_eurostat("une_rt_m", time_format = "date")
unique(raw_une$geo)   # See all available country/area codes
unique(raw_une$age)   # See all available age bracket codes
# Filter for EA20, Seasonally Adjusted, Total Age (15-74), Total Sex, Percentage Rate
clean_une <- raw_une %>%
  filter(
    geo == "EA21",
    s_adj == "SA",
    age == "TOTAL",
    sex == "T",
    unit == "PC_ACT",
    TIME_PERIOD >= as.Date("1999-01-01") & TIME_PERIOD <= as.Date("2025-12-31")
  ) %>%
  rename(Country = geo) %>%
  mutate(Country = if_else(Country == "EA21", "EA", Country)) %>%
  select(Date = TIME_PERIOD,Country, Unemployment_Rate = values)

head(clean_une)

##############################
# Economic Sentiment Indicator
##############################

# Download the monthly Economic Sentiment Indicator table
raw_esi <- get_eurostat("ei_bssi_m_r2", time_format = "date")

# Filter for Euro Area (EA20 or EA21), Seasonally Adjusted ESI Index
clean_esi <- raw_esi %>%
  filter(
    geo == "EA21",          # Replace with "EA21" if using the expanded aggregate
    indic == "BS-ESI-I",    # ESI index (long-term mean = 100)
    s_adj == "SA",           # Seasonally adjusted
    TIME_PERIOD >= as.Date("1999-01-01") & TIME_PERIOD <= as.Date("2025-12-31")
  ) %>%
  rename(Country = geo) %>%
  mutate(Country = if_else(Country == "EA21", "EA", Country)) %>%
  select(Date = TIME_PERIOD,Country, ESI = values) %>%
  arrange(Date)

head(clean_esi)




##############
# HICP
##############

countries <- c("U2") 

# Construct the ECB filter query
geo_filter <- paste(countries, collapse = "+")
ecb_key    <- sprintf("ICP.M.%s.Y.000000.3.INX", geo_filter)

# Download official Seasonally Adjusted HICP
clean_HICP <- get_data(ecb_key) %>%
  select(ref_area, obstime, obsvalue) %>%
  rename(
    Country = ref_area,
    Date = obstime,
    HICP_Index = obsvalue
  ) %>%
  mutate(Date = as.Date(paste0(Date, "-01")),
         Country = if_else(Country == "U2", "EA", Country)) %>%
  arrange(Country, Date)


############
# Core HICP
############

ecb_core_key <- sprintf("ICP.M.%s.Y.XEF000.3.INX", geo_filter)

# Download official Seasonally Adjusted Core HICP
clean_core_HICP <- get_data(ecb_core_key) %>%
  select(ref_area, obstime, obsvalue) %>%
  rename(
    Country = ref_area,
    Date = obstime,
    Core_HICP_Index = obsvalue
  ) %>%
  mutate(Date = as.Date(paste0(Date, "-01")),
         Country = if_else(Country == "U2", "EA", Country)) %>%
  arrange(Country, Date)

################
# PCCI
################

ecb_PCCI_key <- sprintf("ICP.M.U2.N.PCCI00.3.3MM", geo_filter)

# Download official Seasonally Adjusted Core HICP
clean_PCCI <- get_data(ecb_PCCI_key) %>%
  select(ref_area, obstime, obsvalue) %>%
  rename(
    Country = ref_area,
    Date = obstime,
    PCCI = obsvalue
  ) %>%
  mutate(Date = as.Date(paste0(Date, "-01")),
         Country = if_else(Country == "U2", "EA", Country)) %>%
  arrange(Country, Date)

#################
# Merging
#################

long_macro_data <- clean_ip %>%
  # full_join ensures we don't lose a row if one indicator is temporarily missing
  full_join(clean_HICP, by = c("Country", "Date")) %>%
  full_join(clean_core_HICP, by = c("Country", "Date")) %>%
  full_join(clean_PCCI, by = c("Country", "Date")) %>%
  full_join(clean_une, by = c("Country", "Date"))%>%
  full_join(clean_esi, by = c("Country", "Date"))%>%
  arrange(Country, Date)

sort(unique(long_macro_data$Country))

# transforming data to logs and first differences

long_macro_data <- long_macro_data %>%
  group_by(Country) %>%
  arrange(Date) %>%
  mutate(
    Log_IP = log(IP_Index)*100,
    df_IP = Log_IP - lag(Log_IP, 1),
    Log_HICP = log(HICP_Index)*100,
    df_HICP = Log_HICP - lag(Log_HICP, 1),
    Log_Core_HICP = log(Core_HICP_Index)*100,
    df_Core_HICP = Log_Core_HICP - lag(Log_Core_HICP, 1),
    df_Unemployment_Rate = Unemployment_Rate - lag(Unemployment_Rate, 1)
  ) %>%
  select(Country, Date, IP_Index, Log_IP, HICP_Index, Log_HICP, Core_HICP_Index, 
         Log_Core_HICP, df_IP, df_HICP, df_Core_HICP, PCCI, Unemployment_Rate, 
         df_Unemployment_Rate, ESI) %>%
  ungroup()
  
