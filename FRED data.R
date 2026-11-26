
# Download FRED data on oil prices
url_oil <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=MCOILBRENTEU"
oil_raw <- read.csv(url_oil)

# Clean and Format
oil_data <- oil_raw %>%
  mutate(
    Date = as.Date(observation_date),
    Brent_Price = as.numeric(MCOILBRENTEU)
  ) %>%
  filter(!is.na(Brent_Price)) %>%
  mutate(
    Date = floor_date(Date, "month"),
    Log_Oil = log(Brent_Price)*100,
    df_Oil = Log_Oil-lag(Log_Oil, 1)
  ) %>%
  select(Date, Log_Oil, df_Oil)


# VIX

url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=VIXCLS"
vix_raw <- read.csv(url, na.strings = ".")

# Clean and convert dates
vix_raw$Date <- as.Date(vix_raw$observation_date)
vix_clean    <- na.omit(vix_raw)

vix_clean$Date <- as.Date(format(vix_clean$Date, "%Y-%m-01"))

#Aggregate the monthly mean by the new 1st-of-the-month date
vix_clean <- aggregate(VIXCLS ~ Date, data = vix_clean, FUN = mean)

# Rename the column for clarity
names(vix_clean)[2] <- "VIX"

vix_clean <- vix_clean %>%
  mutate(Log_VIX= log(VIX)*100)

# View class and result
class(vix_clean)
head(vix_clean)




