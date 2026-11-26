##########################################################
# EA-MPD adjustment code with Martin Baumgaertner package
##########################################################

if (!requireNamespace("hfdshocks", quietly = TRUE)) {
  devtools::install_github("https://github.com/martinbaumgaertner/hfdshocks.git")
}
library(hfdshocks)
EA_MPD <- read_excel("EA-MPD.xlsx", col_types = c("date", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric"))

raw_output<-ecb_shocks(date_range=c("2000-01-01","2025-12-31"),crisis_date="2008-09-04",path=".",
           reproduce=F,extended=F,extended_release_date="2015-12-03",
           include_loadings=T,return_data = "all")

# Extract the factors and rename the overlapping QE columns
release_shocks <- raw_output$factors$release

conference_shocks <- raw_output$factors$conference %>%
  rename(QE_Shock = qe)

EA_MPD <- EA_MPD %>%
  select(date, OIS_3M, OIS_1Y, OIS_2Y, STOXX50)

# Merge them into a single daily dataframe
daily_shocks <- EA_MPD %>%
  full_join(release_shocks, by="date")%>% 
  full_join (conference_shocks, by = "date")

# Clean up any NAs (if a shock didn't occur in one window, it is a 0)
daily_shocks[is.na(daily_shocks)] <- 0

# Aggregate from Daily to Monthly
monthly_shocks <- daily_shocks %>%
  mutate(date = as.Date(date)) %>%
  mutate(Date = floor_date(date, "month")) %>%
  
  # Group by month and sum ALL the shocks
  group_by(Date) %>%
  summarise(
    Target_Shock = sum(target, na.rm = TRUE),
    Timing_Shock = sum(timing, na.rm = TRUE),
    FG_Shock     = sum(fg, na.rm = TRUE),
    QE_Shock      = sum(QE_Shock, na.rm = TRUE),
    OIS_3M       = sum(OIS_3M, na.rm =TRUE),
    OIS_1Y       = sum(OIS_1Y, na.rm =TRUE),
    OIS_2Y       = sum(OIS_2Y, na.rm =TRUE),
    STOXX50      = sum(STOXX50, na.rm=TRUE)
  ) %>%
  ungroup() %>%
  
# Fill in the missing months with ZEROS
  complete(
    Date = seq.Date(from = as.Date("2000-01-01"), 
                          to = as.Date("2025-12-01"), 
                          by = "month"),
    fill = list(
      Target_Shock = 0, 
      Timing_Shock = 0, 
      FG_Shock = 0, 
      QE_Shock = 0, 
      OIS_3M       = 0,
      OIS_1Y       = 0,
      OIS_2Y       = 0,
      STOXX50      = 0
    )
  )

# Check the final output
head(monthly_shocks)
summary(monthly_shocks)
# Create a correlation matrix of just the shock columns
cor(monthly_shocks %>% dplyr::select(-Date))


