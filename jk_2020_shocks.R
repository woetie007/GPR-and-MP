
df_release <- read_excel("release and conference data.xlsx", 
                         sheet = "df_release", col_types = c("date", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric"))

df_conference <- read_excel("release and conference data.xlsx", 
                            sheet = "df_conference", col_types = c("date", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric"))

#Merging conference and release data
df_merged <- full_join(
  df_release, 
  df_conference, 
  by = "date", 
  suffix = c("_rel", "_conf")
)

# Sum the surprises, treating NAs (missing windows) as 0
df_jk2020_shocks <- df_merged %>%
  mutate(
    # Replace NA with 0 so the math doesn't return NA for missing windows
    OIS_3M_total  = replace_na(OIS_3M_rel, 0) + replace_na(OIS_3M_conf, 0),
    OIS_1Y_total  = replace_na(OIS_1Y_rel, 0) + replace_na(OIS_1Y_conf, 0),
    STOXX50_total = replace_na(STOXX50_rel, 0) + replace_na(STOXX50_conf, 0)
  ) %>%
  # Select the final columns needed for your state-dependent regressions
  select(
    date, 
    OIS_3M_total, 
    OIS_1Y_total, 
    STOXX50_total
  ) %>%
  arrange(date)

# View the first few rows to verify the summation
head(df_jk2020_shocks)

df_monthly_shocks <- df_jk2020_shocks %>%
  mutate(
    # Floor the date to the first of the month
    month_yr = floor_date(date, "month"),
  ) %>%
  group_by(month_yr) %>%
  summarise(
    # Sum the daily surprises for the month
    jk_OIS_3M_monthly  = sum(OIS_3M_total, na.rm = TRUE),
    jk_OIS_1Y_monthly  = sum(OIS_1Y_total, na.rm = TRUE),
    jk_STOXX50_monthly = sum(STOXX50_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Force a continuous sequence of months, filling missing ones with 0
  complete(
    month_yr = seq.Date(min(month_yr), max(month_yr), by = "month"),
    fill = list(
      jk_OIS_3M_monthly = 0,
      jk_OIS_1Y_monthly = 0,
      jk_STOXX50_monthly = 0
    )
  ) %>%
  arrange(month_yr) %>%
  select(Date=month_yr, jk_OIS_3M_monthly, jk_OIS_1Y_monthly, jk_STOXX50_monthly )

# View the complete aggregated monthly data
head(df_monthly_shocks)

library(readr)
data <- read_csv("data.csv")
