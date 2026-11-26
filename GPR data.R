
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)


GPR <- read_excel("GPR.xls")
EPU <- read_excel("Europe_Policy_Uncertainty_Data.xlsx")
# Convert the date column to Date format
GPR$Date <- as.Date(GPR$month)
GPR <- dplyr::select(GPR, -month)
long_gpr <- GPR %>%
  # Pivot from wide to long
  pivot_longer(
    cols = starts_with("GPRC_"),   # Only pivot the country columns
    names_to = "Country",          # New column name for the countries
    names_prefix = "GPRC_",        # strips "GPRC_" from the values
    values_to = "GPRC"      # New column name for the actual values
  ) %>%
  # Convert 3-letter codes to Eurostat 2-letter codes
  mutate(
    Country = countrycode(
      sourcevar = Country, 
      origin = "iso3c",            # GPR uses standard 3-letter ISO codes
      destination = "eurostat"     # Converts to Eurostat's specific 2-letter format
    )
  ) %>%
  arrange(Country, Date)
sort(unique(long_gpr$Country))

# EA20 row
ea20_gpr <- long_gpr %>%
  distinct(Date, .keep_all = TRUE) %>%    # Grabs just one row per month
  mutate(
    Country = "EA",                     # Name this new region EA20
    GPRC = GPR                            # Copy the total GPR into the GPRC column
  ) %>%
  select(Date, Country, GPR, GPRT, GPRA)

# Stack them together into one perfect dataset
clean_gpr_long <- bind_rows(ea20_gpr) %>%
  mutate(Date = floor_date(as.Date(Date), "month"))
clean_gpr_long <- clean_gpr_long %>%
  mutate(Log_GPR = log(GPR)*100,
         df_GPR = Log_GPR - lag(Log_GPR, 1),
         Log_GPRT = log(GPRT)*100,
         df_GPRT = Log_GPRT - lag(Log_GPRT, 1),
         Log_GPRA= log(GPRA)*100,
         df_GPRA = Log_GPRA - lag(Log_GPRA, 1))



#=========
# EPU
#=========

EPU <- EPU %>%
  # Combine year, month, and day into a properly formatted Date object
  mutate(Date = make_date(Year, Month, 1),
         Log_EPU = log(EPU)*100) %>%
  # Keep only the columns you want
  select(Date, EPU, Log_EPU)

# View the result to confirm
head(EPU)



