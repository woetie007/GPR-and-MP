###################################################
# Data statistics and summary tables
###################################################
conflicted::conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(ecb::get_data)
conflicted::conflicts_prefer(plotly::layout)
conflicts_prefer(dplyr::first)
source("Merge.R")

#==========================
# Summary statistics
#==========================

# ---------------
# Prepare Data
# ---------------
summary_prep <- final_data %>%
  select(
    `OIS 3 Month Shock` = OIS_3M,
    `OIS 1 Year Shock` = OIS_1Y,
    `OIS 2 Year Shock` = OIS_2Y,
    `STOXX50 Shock` = STOXX50,
    `Target Shock` = Target_Shock,
    `Forward Guidance Shock`  = FG_Shock,
    `QE Shock` = QE_Shock,
    `HICP (Index: 2015=100)`  = HICP_Index,
    `Core HICP (Index: 2015=100)` = Core_HICP_Index,
    `PCCI` = PCCI,
    `Ind. Production (2021=100)` = IP_Index,
    `Economic Sentiment Index (ESI)` = ESI,
    `Unemployment Rate (%)` = Unemployment_Rate,
    `Inflation Linked Rate 1Y (%)`= EU_ILS_1Y,
    `Inflation Linked Forward Rate 5Y5Y (%)`= EU_ILS_5Y5Y,
    `Inflation Linked Rate 5Y (%)`= EU_ILS_5Y, 
    `Inflation Linked Rate 10Y (%)`= EU_ILS_10Y,
    `EONIA / €STR 3M (%)` = Spliced_EONIA_3M,
    `EONIA / €STR 1Y (%)` = Spliced_EONIA_1Y,
    `EONIA / €STR 2Y (%)` = Spliced_EONIA_2Y,
    `EONIA / €STR 5Y (%)` = Spliced_EONIA_5Y,
    `EONIA / €STR 10Y (%)` = Spliced_EONIA_10Y,
    `EONIA / €STR / Shadow Rate 3M (%)`   = True_Policy_Rate,
    `MRO (%)` = MRO,
    `DFR (%)` = Deposit_Rate,
    `Global GPR Index`        = GPR,
    `Global Geopolitical Acts Index` = GPRA,
    `Global Geopolitical Threats Index` = GPRT,
    `EPU` = EPU,
    `VIX` = VIX,
    `VSTOXX` = VSTOXX
  )

#-----------------------------
# Create correlation matrix
#-----------------------------
corr_data <- summary_prep %>%
  select(
    -`DFR (%)`,
    -`Global Geopolitical Acts Index`,
    -`Global Geopolitical Threats Index`,
    -`EONIA / €STR 10Y (%)`,
    -`EONIA / €STR 5Y (%)`,
    -`Inflation Linked Rate 10Y (%)`,
    -`Inflation Linked Rate 5Y (%)`,
    -`OIS 2 Year Shock`
  )

cor_matrix <- cor(corr_data, use = "pairwise.complete.obs")

print(cor_matrix)

latex_table <- cor_matrix %>%
  kable(
    format = "latex", 
    booktabs = TRUE,          
    digits = 2,               
    caption = "Correlation Matrix of Macroeconomic and Financial Variables",
    label = "cor_matrix"      
  ) %>%
  kable_styling(
    latex_options = c("scale_down", "hold_position"), 
    font_size = 8
  ) %>%
  landscape()                 

cat(latex_table)

# ------------------------------------------
# Compute the Summary Statistics Matrix
# ------------------------------------------
summary_stats <- summary_prep %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Variable", 
    values_to = "Value"
  ) %>%
  group_by(Variable) %>%
  summarise(
    Obs = sum(!is.na(Value)),
    Mean      = mean(Value, na.rm = TRUE),
    Median    = median(Value, na.rm=TRUE),
    `Std. Dev.` = sd(Value, na.rm = TRUE),
    Min       = min(Value, na.rm = TRUE),
    Max       = max(Value, na.rm = TRUE),
    `Skew.`      = skewness(Value, na.rm = TRUE),
    `Kurt.` = kurtosis(Value, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  slice(match(c(
    "Target Shock", 
    "Forward Guidance Shock", 
    "QE Shock", 
    "OIS 3 Month Shock", 
    "OIS 1 Year Shock", 
    "OIS 2 Year Shock", 
    "STOXX50 Shock",
    "HICP (Index: 2015=100)", 
    "Core HICP (Index: 2015=100)", 
    "PCCI", 
    "Ind. Production (2021=100)", 
    "Economic Sentiment Index (ESI)", 
    "Unemployment Rate (%)",
    "Inflation Linked Rate 1Y (%)", 
    "Inflation Linked Forward Rate 5Y5Y (%)",
    "Inflation Linked Rate 5Y (%)", 
    "Inflation Linked Rate 10Y (%)",
    "EONIA / €STR 3M (%)",
    "EONIA / €STR 1Y (%)",
    "EONIA / €STR 2Y (%)",
    "EONIA / €STR 5Y (%)",
    "EONIA / €STR 10Y (%)",
    "EONIA / €STR / Shadow Rate 3M (%)", 
    "MRO (%)", 
    "DFR (%)",
    "Global GPR Index", 
    "Global Geopolitical Acts Index", 
    "Global Geopolitical Threats Index",
    "EPU",
    "VIX",
    "VSTOXX"
  ), Variable))

# ------------------------------------
# Add Category Grouping Headers
# ------------------------------------
summary_stats_grouped <- summary_stats %>%
  mutate(
    Category = case_when(
      Variable %in% c("HICP (Index: 2015=100)", "Core HICP (Index: 2015=100)", "PCCI", 
                      "Ind. Production (2021=100)", "Economic Sentiment Index (ESI)", 
                      "Unemployment Rate (%)", "Inflation Linked Rate 1Y (%)", "Inflation Linked Forward Rate 5Y5Y (%)",
                      "Inflation Linked Rate 5Y (%)", "Inflation Linked Rate 10Y (%)" ) ~ "Macroeconomic Fundamentals",
      Variable %in% c("EONIA / €STR 3M (%)", "EONIA / €STR 1Y (%)", "EONIA / €STR 2Y (%)",
                      "EONIA / €STR 5Y (%)", "EONIA / €STR 10Y (%)", "EONIA / €STR / Shadow Rate 3M (%)",
                      "MRO (%)", "DFR (%)") ~ "Policy & Interest Rates",
      Variable %in% c("Global GPR Index", "Global Geopolitical Acts Index", 
                      "Global Geopolitical Threats Index","EPU", "VIX", "VSTOXX") ~ "Uncertainty & Risk Indices",
      Variable %in% c("Target Shock", "Forward Guidance Shock", "QE Shock") ~ "Monetary Policy Shocks (EA-MPD)",
      Variable %in% c("OIS 3 Month Shock", "OIS 1 Year Shock",
                      "OIS 2 Year Shock", "STOXX50 Shock") ~ "High-Frequency Market Reactions (EA-MPD)"
    )
  )

# ------------------------------------------------
# Split and Generate the Formatted Tables 
# ------------------------------------------------

table_normal <- summary_stats_grouped %>%
  filter(Category %in% c("Macroeconomic Fundamentals", "Policy & Interest Rates", "Uncertainty & Risk Indices")) %>%
  mutate(Category = factor(Category, levels = c(
    "Policy & Interest Rates",
    "Macroeconomic Fundamentals",
    "Uncertainty & Risk Indices"
  ))) %>%
  arrange(Category) %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Table 1A: Descriptive Statistics - Macroeconomic & Financial Variables**"),
    subtitle = "Euro Area Sample (2000:M1 – 2025:M10)"
  ) %>%
  fmt_number(
    columns = c(Mean, Median, `Std. Dev.`, Min, Max, `Skew.`, `Kurt.`),
    decimals = 2
  ) %>%
  fmt_integer(columns = Obs) %>%
  cols_align(align = "left", columns = Variable) %>%
  cols_align(align = "center", columns = c(Obs, Mean, Median, `Std. Dev.`, Min, Max,  `Skew.`, `Kurt.`)) %>%
  tab_footnote(
    footnote = "In the empirical estimations, macroeconomic (HICP, Core HICP and Ind. Production) and uncertainty indices are log-transformed."
  ) %>%
  tab_options(
    table.width = pct(90),
    heading.title.font.size = px(16),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    row_group.background.color = "#f2f2f2"
  )

table_shocks <- summary_stats_grouped %>%
  filter(Category %in% c("Monetary Policy Shocks (EA-MPD)", "High-Frequency Market Reactions (EA-MPD)")) %>%
  mutate(Category = factor(Category, levels = c(
    "High-Frequency Market Reactions (EA-MPD)",
    "Monetary Policy Shocks (EA-MPD)"
  ))) %>%
  arrange(Category) %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Table 1B: Descriptive Statistics - High-Frequency Reactions & MP Shocks**"),
    subtitle = "Euro Area Sample (2000:M1 – 2025:M10)"
  ) %>%
  fmt_number(
    columns = c(Mean, Median, `Std. Dev.`, Min, Max, `Skew.`, `Kurt.`),
    decimals = 3
  ) %>%
  fmt_integer(columns = Obs) %>%
  cols_align(align = "left", columns = Variable) %>%
  cols_align(align = "center", columns = c(Obs, Mean, Median, `Std. Dev.`, Min, Max,  `Skew.`, `Kurt.`)) %>%
  tab_footnote(
    footnote = "Monetary policy shocks and asset price reactions are measured in basis points or percentage points."
  ) %>%
  tab_options(
    table.width = pct(90),
    heading.title.font.size = px(16),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    row_group.background.color = "#f2f2f2"
  )

table_normal
table_shocks

gtsave(table_normal, filename = "Table1A_Descriptive_Stats_Normal.tex")
gtsave(table_shocks, filename = "Table1B_Descriptive_Stats_Shocks.tex")


#=========================
# Plots of all variables
#=========================

plot_df <- final_data %>%
  filter(Country == "EA") %>%
  mutate(Date = as.Date(Date))

events <- data.frame(
  Date = as.Date(c("2008-09-15", "2020-03-01")),
  Event = c("GFC", "COVID-19")
)

# ==============================================================================
# FIGURE 1A: Macroeconomic Fundamentals & Policy Rates
# ==============================================================================
data_2000_A <- plot_df %>%
  select(Date, 
         `Log HICP` = Log_HICP, 
         `Log Core HICP` = Log_Core_HICP, 
         `Log IP` = Log_IP,
         `Unemployment Rate` = Unemployment_Rate, 
         `ESI` = ESI, 
         `True Policy Rate` = True_Policy_Rate, 
         `df HICP` = df_HICP,
         `df Core HICP` = df_Core_HICP, 
         `df True Policy Rate` = df_True_Policy_Rate) %>%
  pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value") %>%
  drop_na()

p_2000_A <- ggplot(data_2000_A, aes(x = Date, y = Value)) +
  geom_line(color = "midnightblue") +
  geom_vline(data = events, aes(xintercept = Date), color = "red", linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(subtitle = "Red dashed lines indicate GFC (Sept 2008) and COVID-19 (Mar 2020)",
       x = NULL, y = NULL) +
  theme(strip.background = element_rect(fill = "#f2f2f2"),
        strip.text = element_text(face = "bold"))

print(p_2000_A)
ggsave("Figure1A_Macro_Rates.png", p_2000_A, width = 12, height = 8, dpi = 300)


# ==============================================================================
# FIGURE 1B: Uncertainty & Financial Markets (Excluding Shocks)
# ==============================================================================
data_2000_B <- plot_df %>%
  select(Date, 
         `Log Oil Prices` = Log_Oil, 
         `EONIA 1Y` = Spliced_EONIA_1Y, 
         `Log GPR` = Log_GPR,
         `Log VIX` = Log_VIX,
         `Log VSTOXX` = Log_VSTOXX, 
         `Log EPU` = Log_EPU) %>%
  pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value") %>%
  drop_na()

p_2000_B <- ggplot(data_2000_B, aes(x = Date, y = Value)) +
  geom_line(color = "darkred") +
  geom_vline(data = events, aes(xintercept = Date), color = "red", linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) + # 6 variables = 2 rows of 3
  theme_bw() +
  labs(subtitle = "Red dashed lines indicate GFC (Sept 2008) and COVID-19 (Mar 2020)",
       x = NULL, y = NULL) +
  theme(strip.background = element_rect(fill = "#f2f2f2"),
        strip.text = element_text(face = "bold"))

print(p_2000_B)
ggsave("Figure1B_Markets_Uncertainty.png", p_2000_B, width = 12, height = 6, dpi = 300)


# ==============================================================================
# FIGURE 1C: High-Frequency Reactions & Monetary Policy Shocks
# ==============================================================================
data_shocks <- plot_df %>%
  select(Date, 
         `OIS 3M Surprise` = OIS_3M, 
         `OIS 1Y Surprise` = OIS_1Y, 
         `STOXX50 Surprise` = STOXX50,
         `Target Shock` = Target_Shock,
         `Forward Guidance Shock` = FG_Shock,
         `QE Shock` = QE_Shock) %>%
  pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value") %>%
  drop_na()

p_shocks <- ggplot(data_shocks, aes(x = Date, y = Value)) +
  geom_line(color = "purple4") +
  geom_vline(data = events, aes(xintercept = Date), color = "red", linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) + # 6 variables = 2 rows of 3
  theme_bw() +
  labs(subtitle = "Red dashed lines indicate GFC (Sept 2008) and COVID-19 (Mar 2020)",
       x = NULL, y = NULL) +
  theme(strip.background = element_rect(fill = "#f2f2f2"),
        strip.text = element_text(face = "bold"))

print(p_shocks)
ggsave("Figure1C_Shocks.png", p_shocks, width = 12, height = 6, dpi = 300)

# ==============================================================================
# FIGURE 2: PCCI (Starts in 2001)
# ==============================================================================
data_pcci <- plot_df %>%
  select(Date, `PCCI`) %>%
  drop_na()

p_pcci <- ggplot(data_pcci, aes(x = Date, y = PCCI)) +
  geom_line(color = "forestgreen", linewidth = 1) +
  geom_vline(data = events, aes(xintercept = Date), color = "red", linetype = "dashed", alpha = 0.7) +
  theme_bw() +
  labs(subtitle = "Red dashed lines indicate GFC (Sept 2008) and COVID-19 (Mar 2020)",
       x = NULL, y = "Index Value")

print(p_pcci)
ggsave("Figure2_PCCI.png", p_pcci, width = 8, height = 4, dpi = 300)


# ==============================================================================
# FIGURE 3: Inflation-Linked Swaps (Starts in 2008)
# ==============================================================================
data_ils <- plot_df %>%
  select(Date, 
         `EU ILS 1Y` = EU_ILS_1Y, 
         `EU ILS 5Y5Y` = EU_ILS_5Y5Y) %>%
  pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value") %>%
  drop_na()

p_ils <- ggplot(data_ils, aes(x = Date, y = Value)) +
  geom_line(color = "darkorange3", linewidth = 1) +
  geom_vline(data = events, aes(xintercept = Date), color = "red", linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(subtitle = "Red dashed lines indicate GFC (Sept 2008) and COVID-19 (Mar 2020)",
       x = NULL, y = "Rate (%)") +
  theme(strip.background = element_rect(fill = "#f2f2f2"),
        strip.text = element_text(face = "bold"))

print(p_ils)
ggsave("Figure3_ILS_Rates.png", p_ils, width = 10, height = 4, dpi = 300)



#================================
# High vs Low GPR statistics
#================================

data_stats <- final_data %>%
  arrange(Date) %>%
  mutate(
    gpr_mean = median(GPR, na.rm = TRUE),
    regime   = ifelse(GPR > gpr_mean, "High GPR", "Low GPR"),
    
    rate_change     = MRO - lag(MRO),
    abs_rate_change = abs(rate_change),
    
    is_move = ifelse(abs_rate_change > 0.001 & !is.na(rate_change), 1, 0),
    
    spell_id = cumsum(regime != lag(regime, default = first(regime)))
  )

# -----------------------------------------------------
# Episode-Level Aggregates (Consecutive GPR Spells)
# -----------------------------------------------------
episode_stats <- data_stats %>%
  group_by(spell_id, regime) %>%
  summarise(
    spell_length_months = n(),
    num_rate_moves      = sum(is_move, na.rm = TRUE),
    net_rate_change     = sum(rate_change, na.rm = TRUE),
    cum_abs_change      = sum(abs_rate_change, na.rm = TRUE), 
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Compute Summary Table Statistics
# ------------------------------------------------------------------------------

monthly_summary <- data_stats %>%
  filter(!is.na(rate_change)) %>%
  group_by(regime) %>%
  summarise(
    avg_monthly_change = mean(abs_rate_change, na.rm = TRUE),
    
    avg_active_move_size = mean(abs_rate_change[is_move == 1], na.rm = TRUE),
    
    pct_active_months = mean(is_move, na.rm = TRUE) * 100,
    
    pct_hikes = mean(rate_change > 0.001, na.rm = TRUE) * 100,
    pct_cuts  = mean(rate_change < -0.001, na.rm = TRUE) * 100,
    
    .groups = "drop"
  )

episode_summary <- episode_stats %>%
  group_by(regime) %>%
  summarise(
    avg_spell_duration_months  = mean(spell_length_months),
    max_spell_duration_months  = max(spell_length_months),
    max_moves_in_single_spell  = max(num_rate_moves),
    avg_net_change_per_spell   = mean(net_rate_change),
    max_cum_change_in_spell    = max(cum_abs_change),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Combine into a Clean Introductory Table
# ------------------------------------------------------------------------------
intro_stats_table <- monthly_summary %>%
  left_join(episode_summary, by = "regime") %>%
  pivot_longer(-regime, names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = regime, values_from = Value)

print(intro_stats_table)

library(kableExtra)
library(tidyverse)

intro_stats_table %>%
  mutate(
    Metric = case_when(
      Metric == "avg_monthly_change" ~ "Avg. Monthly Change (Absolute)",
      Metric == "avg_active_move_size" ~ "Avg. Active Move Size",
      Metric == "pct_active_months" ~ "Frequency of Active Months (%)",
      Metric == "pct_hikes" ~ "Percentage of Hikes (%)",
      Metric == "pct_cuts" ~ "Percentage of Cuts (%)",
      Metric == "avg_spell_duration_months" ~ "Avg. Spell Duration (Months)",
      Metric == "max_spell_duration_months" ~ "Max. Spell Duration (Months)",
      Metric == "max_moves_in_single_spell" ~ "Max. Moves in Single Spell",
      Metric == "avg_net_change_per_spell" ~ "Avg. Net Change per Spell",
      Metric == "max_cum_change_in_spell" ~ "Max. Cumulative Change in Spell",
      TRUE ~ Metric
    )
  ) %>%
  arrange(match(Metric, c(
    "Avg. Monthly Change (Absolute)", 
    "Avg. Active Move Size",
    "Frequency of Active Months (%)", 
    "Percentage of Hikes (%)", 
    "Percentage of Cuts (%)",
    "Avg. Spell Duration (Months)", 
    "Max. Spell Duration (Months)",
    "Max. Moves in Single Spell", 
    "Avg. Net Change per Spell", 
    "Max. Cumulative Change in Spell"
  ))) %>%
  kbl(format = "latex", booktabs = TRUE, align = "lcc", digits = 3,
      caption = "Monetary Policy Dynamics: High vs. Low Geopolitical Risk Regimes") %>%
  kable_styling(latex_options = c("hold_position"), font_size = 10) %>%
  pack_rows("Monthly / Action Level Metrics", 1, 5, latex_gap_space = "2.2pt") %>%
  pack_rows("Episode / Consecutive Period Metrics", 6, 10, latex_gap_space = "2.2pt") %>%
  footnote(general = "High and Low GPR regimes are defined based on the median of the Global Geopolitical Risk (GPR) Index. Active moves are defined as absolute rate changes strictly greater than zero.", 
           threeparttable = TRUE) %>%
  cat()




#================================
# High vs Low GPR statistics for Macro Variables (Inflation & IP Growth)
#================================

data_stats_macro <- final_data %>%
  arrange(Date) %>%
  mutate(
    # Define High vs Low GPR Threshold
    gpr_median = median(GPR, na.rm = TRUE),
    regime     = ifelse(GPR > gpr_median, "High GPR", "Low GPR"),
    
    # Create an ID for each consecutive GPR spell/episode
    spell_id = cumsum(regime != lag(regime, default = first(regime)))
  )

# -----------------------------------------------------
# Episode-Level Aggregates (Consecutive GPR Spells)
# -----------------------------------------------------
episode_stats_macro <- data_stats_macro %>%
  group_by(spell_id, regime) %>%
  summarise(
    spell_length_months = n(),
    cum_inflation = sum(df_HICP, na.rm = TRUE),
    cum_ip_growth = sum(df_IP, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------
# Compute Summary Table Statistics
# ------------------------------------

monthly_summary_macro <- data_stats_macro %>%
  group_by(regime) %>%
  summarise(
    `Inflation: Avg Monthly Change`        = mean(df_HICP, na.rm = TRUE),
    `Inflation: Volatility (Std. Dev.)`    = sd(df_HICP, na.rm = TRUE),
    `Inflation: Minimum Monthly Change`    = min(df_HICP, na.rm = TRUE),
    `Inflation: Maximum Monthly Change`    = max(df_HICP, na.rm = TRUE),
    `Inflation: % Months Positive (Rise)`  = mean(df_HICP > 0, na.rm = TRUE) * 100,
    
    `IP Growth: Avg Monthly Change`        = mean(df_IP, na.rm = TRUE),
    `IP Growth: Volatility (Std. Dev.)`    = sd(df_IP, na.rm = TRUE),
    `IP Growth: Minimum Monthly Change`    = min(df_IP, na.rm = TRUE),
    `IP Growth: Maximum Monthly Change`    = max(df_IP, na.rm = TRUE),
    `IP Growth: % Months Positive (Rise)`  = mean(df_IP > 0, na.rm = TRUE) * 100,
    
    .groups = "drop"
  ) %>%
  pivot_longer(-regime, names_to = "Metric", values_to = "Value")

episode_summary_macro <- episode_stats_macro %>%
  group_by(regime) %>%
  summarise(
    `Inflation: Avg Cumulative Change per Spell` = mean(cum_inflation, na.rm = TRUE),
    `IP Growth: Avg Cumulative Change per Spell` = mean(cum_ip_growth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-regime, names_to = "Metric", values_to = "Value")

# ----------------------------------------
# Combine into a Clean Summary Table
# ----------------------------------------
macro_stats_table <- bind_rows(monthly_summary_macro, episode_summary_macro) %>%
  pivot_wider(names_from = regime, values_from = Value) %>%
  arrange(Metric)

print(macro_stats_table)
library(kableExtra)
library(tidyverse)
macro_stats_table %>%
  mutate(Metric = str_remove(Metric, "Inflation: |IP Growth: ")) %>%
  kbl(format = "latex", booktabs = TRUE, align = "lcc", 
      caption = "Macroeconomic Dynamics: High vs. Low Geopolitical Risk Regimes") %>%
  kable_styling(latex_options = c("hold_position"), font_size = 10) %>%
  pack_rows("Inflation (HICP)", 7, 12, latex_gap_space = "2.2pt") %>%
  pack_rows("Industrial Production Growth", 1, 6, latex_gap_space = "2.2pt") %>%
  footnote(general = "High and Low GPR regimes are defined based on the median of the Global Geopolitical Risk (GPR) Index.", 
           threeparttable = TRUE) %>%
  cat()

#=======================================
# Plotting the variables together by dates
#========================================

library(tidyverse)
library(patchwork)
library(lubridate)

# ------------------------------------------------------------------------------
# Prepare the Data & Shading Regions
# ------------------------------------------------------------------------------
plot_data <- final_data %>%
  mutate(Date = as.Date(Date)) %>% 
  arrange(Date) %>%
  mutate(
    gpr_median = median(GPR, na.rm = TRUE),
    regime = ifelse(GPR > gpr_median, "High GPR Regime", "Low GPR Regime"),
    
    spell_id = cumsum(regime != lag(regime, default = first(regime)))
  )

shading_data <- plot_data %>%
  group_by(spell_id, regime) %>%
  summarise(
    start_date = min(Date, na.rm = TRUE),
    end_date = max(Date, na.rm = TRUE) + months(1), 
    .groups = "drop"
  )

regime_colors <- c("High GPR Regime" = "red", "Low GPR Regime" = "blue")

event_lines <- as.Date(c("2008-09-15", "2020-03-01"))

# ------------------------------------------------------------------------------
# Build the Individual Plots
# ------------------------------------------------------------------------------

p_ip <- ggplot() +
  geom_rect(data = shading_data, 
            aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf, fill = regime), 
            alpha = 0.25) +
  geom_vline(xintercept = event_lines, linetype = "dashed", color = "black") +
  geom_line(data = plot_data, aes(x = Date, y = IP_Index), color = "forestgreen", linewidth = 1) +
  scale_fill_manual(values = regime_colors) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = c(0, 0)) +
  theme_bw() +
  labs(
    title = "Industrial Production (IP), Headline Inflation (HICP) Policy Rate, Euro Area",
    y = "IP (Index)", 
    x = NULL
  ) +
  theme(
    legend.position = c(0.15, 0.85), 
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_line(color = "grey70")
  )

p_hicp <- ggplot() +
  geom_rect(data = shading_data, 
            aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf, fill = regime), 
            alpha = 0.25) +
  geom_vline(xintercept = event_lines, linetype = "dashed", color = "black") +
  geom_line(data = plot_data, aes(x = Date, y = df_HICP), color = "mediumblue", linewidth = 1) +
  scale_fill_manual(values = regime_colors) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = c(0, 0)) +
  theme_bw() +
  labs(title = "HICP (Headline Inflation)", y = "HICP (%)", x = NULL) +
  theme(
    legend.position = "none", 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid.major = element_line(color = "grey70")
  )

p_rate <- ggplot() +
  geom_rect(data = shading_data, 
            aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf, fill = regime), 
            alpha = 0.25) +
  geom_vline(xintercept = event_lines, linetype = "dashed", color = "black") +
  geom_line(data = plot_data, aes(x = Date, y = MRO), color = "black", linewidth = 1) +
  scale_fill_manual(values = regime_colors) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = c(0, 0)) +
  theme_bw() +
  labs(title = "Monetary Policy Rate", y = "Policy Rate (%)", x = "Year") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey70")
  )

# ------------------------------------------------------------------------------
# Combine and Render the Plot
# ------------------------------------------------------------------------------

final_plot <- p_ip / p_hicp / p_rate
print(final_plot)

ggsave("Regime_Macro_Plot.png", final_plot, width = 12, height = 8, dpi = 300)
