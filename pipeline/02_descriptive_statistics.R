###############################################################################
# 02_descriptive_statistics.R
# Port of Desc_Stat_and_Stylized_Facts.R
# Thesis: Tables 1-4, Table A6, Figures B1-B6
###############################################################################
out <- DIRS$descriptive

#==========================
# Summary statistics
#==========================
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
# Correlation matrix (Table A6)
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
save_tex(latex_table, out, "TableA6_Correlation_Matrix")

# ------------------------------------------
# Summary statistics matrix (Tables 1 and 4)
# ------------------------------------------
var_order <- c(
  "Target Shock", "Forward Guidance Shock", "QE Shock",
  "OIS 3 Month Shock", "OIS 1 Year Shock", "OIS 2 Year Shock", "STOXX50 Shock",
  "HICP (Index: 2015=100)", "Core HICP (Index: 2015=100)", "PCCI",
  "Ind. Production (2021=100)", "Economic Sentiment Index (ESI)",
  "Unemployment Rate (%)", "Inflation Linked Rate 1Y (%)",
  "Inflation Linked Forward Rate 5Y5Y (%)", "Inflation Linked Rate 5Y (%)",
  "Inflation Linked Rate 10Y (%)", "EONIA / €STR 3M (%)", "EONIA / €STR 1Y (%)",
  "EONIA / €STR 2Y (%)", "EONIA / €STR 5Y (%)", "EONIA / €STR 10Y (%)",
  "EONIA / €STR / Shadow Rate 3M (%)", "MRO (%)", "DFR (%)",
  "Global GPR Index", "Global Geopolitical Acts Index",
  "Global Geopolitical Threats Index", "EPU", "VIX", "VSTOXX"
)

summary_stats <- summary_prep %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(
    Obs = sum(!is.na(Value)),
    Mean      = mean(Value, na.rm = TRUE),
    Median    = median(Value, na.rm = TRUE),
    `Std. Dev.` = sd(Value, na.rm = TRUE),
    Min       = min(Value, na.rm = TRUE),
    Max       = max(Value, na.rm = TRUE),
    `Skew.`   = skewness(Value, na.rm = TRUE),
    `Kurt.`   = kurtosis(Value, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  slice(match(var_order, Variable))

summary_stats_grouped <- summary_stats %>%
  mutate(
    Category = case_when(
      Variable %in% c("HICP (Index: 2015=100)", "Core HICP (Index: 2015=100)", "PCCI",
                      "Ind. Production (2021=100)", "Economic Sentiment Index (ESI)",
                      "Unemployment Rate (%)", "Inflation Linked Rate 1Y (%)",
                      "Inflation Linked Forward Rate 5Y5Y (%)",
                      "Inflation Linked Rate 5Y (%)", "Inflation Linked Rate 10Y (%)") ~ "Macroeconomic Fundamentals",
      Variable %in% c("EONIA / €STR 3M (%)", "EONIA / €STR 1Y (%)", "EONIA / €STR 2Y (%)",
                      "EONIA / €STR 5Y (%)", "EONIA / €STR 10Y (%)",
                      "EONIA / €STR / Shadow Rate 3M (%)",
                      "MRO (%)", "DFR (%)") ~ "Policy & Interest Rates",
      Variable %in% c("Global GPR Index", "Global Geopolitical Acts Index",
                      "Global Geopolitical Threats Index", "EPU", "VIX", "VSTOXX") ~ "Uncertainty & Risk Indices",
      Variable %in% c("Target Shock", "Forward Guidance Shock", "QE Shock") ~ "Monetary Policy Shocks (EA-MPD)",
      Variable %in% c("OIS 3 Month Shock", "OIS 1 Year Shock",
                      "OIS 2 Year Shock", "STOXX50 Shock") ~ "High-Frequency Market Reactions (EA-MPD)"
    )
  )

make_gt <- function(df, levels, title, decimals, footnote) {
  df %>%
    filter(Category %in% levels) %>%
    mutate(Category = factor(Category, levels = levels)) %>%
    arrange(Category) %>%
    gt(groupname_col = "Category") %>%
    tab_header(title = md(title), subtitle = "Euro Area Sample (2000:M1 – 2025:M10)") %>%
    fmt_number(columns = c(Mean, Median, `Std. Dev.`, Min, Max, `Skew.`, `Kurt.`),
               decimals = decimals) %>%
    fmt_integer(columns = Obs) %>%
    cols_align(align = "left", columns = Variable) %>%
    cols_align(align = "center",
               columns = c(Obs, Mean, Median, `Std. Dev.`, Min, Max, `Skew.`, `Kurt.`)) %>%
    tab_footnote(footnote = footnote) %>%
    tab_options(
      table.width = pct(90),
      heading.title.font.size = px(16),
      column_labels.font.weight = "bold",
      row_group.font.weight = "bold",
      row_group.background.color = "#f2f2f2"
    )
}

table_normal <- make_gt(
  summary_stats_grouped,
  c("Policy & Interest Rates", "Macroeconomic Fundamentals", "Uncertainty & Risk Indices"),
  "**Table 1A: Descriptive Statistics - Macroeconomic & Financial Variables**",
  2,
  "In the empirical estimations, macroeconomic (HICP, Core HICP and Ind. Production) and uncertainty indices are log-transformed."
)
table_shocks <- make_gt(
  summary_stats_grouped,
  c("High-Frequency Market Reactions (EA-MPD)", "Monetary Policy Shocks (EA-MPD)"),
  "**Table 1B: Descriptive Statistics - High-Frequency Reactions & MP Shocks**",
  3,
  "Monetary policy shocks and asset price reactions are measured in basis points or percentage points."
)

gtsave(table_normal, filename = file.path(out, "Table1_Descriptive_Stats_Macro_Financial.tex"))
gtsave(table_shocks, filename = file.path(out, "Table4_Descriptive_Stats_HF_Shocks.tex"))
message("  wrote Table1 / Table4 descriptive statistics")

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

facet_theme <- theme(strip.background = element_rect(fill = "#f2f2f2"),
                     strip.text = element_text(face = "bold"))
event_subtitle <- "Red dashed lines indicate GFC (Sept 2008) and COVID-19 (Mar 2020)"

facet_plot <- function(df, colour, ncol) {
  ggplot(df, aes(x = Date, y = Value)) +
    geom_line(color = colour) +
    geom_vline(data = events, aes(xintercept = Date), color = "red",
               linetype = "dashed", alpha = 0.7) +
    facet_wrap(~ Variable, scales = "free_y", ncol = ncol) +
    theme_bw() +
    labs(subtitle = event_subtitle, x = NULL, y = NULL) +
    facet_theme
}

# FIGURE B1: Macroeconomic Fundamentals & Policy Rates
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
save_plot(facet_plot(data_2000_A, "midnightblue", 3), out,
          "FigureB1_Macro_Rates", width = 12, height = 8)

# FIGURE B4: Uncertainty & Financial Markets
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
save_plot(facet_plot(data_2000_B, "darkred", 3), out,
          "FigureB4_Markets_Uncertainty", width = 12, height = 6)

# FIGURE B6: High-Frequency Reactions & Monetary Policy Shocks
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
save_plot(facet_plot(data_shocks, "purple4", 3), out,
          "FigureB6_Shocks", width = 12, height = 6)

# FIGURE B2: PCCI (starts in 2001)
data_pcci <- plot_df %>% select(Date, PCCI) %>% drop_na()
p_pcci <- ggplot(data_pcci, aes(x = Date, y = PCCI)) +
  geom_line(color = "forestgreen", linewidth = 1) +
  geom_vline(data = events, aes(xintercept = Date), color = "red",
             linetype = "dashed", alpha = 0.7) +
  theme_bw() +
  labs(subtitle = event_subtitle, x = NULL, y = "Index Value")
save_plot(p_pcci, out, "FigureB2_PCCI", width = 8, height = 4)

# FIGURE B3: Inflation-Linked Swaps (start in 2008)
data_ils <- plot_df %>%
  select(Date, `EU ILS 1Y` = EU_ILS_1Y, `EU ILS 5Y5Y` = EU_ILS_5Y5Y) %>%
  pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value") %>%
  drop_na()
p_ils <- ggplot(data_ils, aes(x = Date, y = Value)) +
  geom_line(color = "darkorange3", linewidth = 1) +
  geom_vline(data = events, aes(xintercept = Date), color = "red",
             linetype = "dashed", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(subtitle = event_subtitle, x = NULL, y = "Rate (%)") +
  facet_theme
save_plot(p_ils, out, "FigureB3_ILS_Rates", width = 10, height = 4)

#================================================
# Table 3: Monetary policy, High vs Low GPR
#================================================
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

episode_stats <- data_stats %>%
  group_by(spell_id, regime) %>%
  summarise(
    spell_length_months = n(),
    num_rate_moves      = sum(is_move, na.rm = TRUE),
    net_rate_change     = sum(rate_change, na.rm = TRUE),
    cum_abs_change      = sum(abs_rate_change, na.rm = TRUE),
    .groups = "drop"
  )

monthly_summary <- data_stats %>%
  filter(!is.na(rate_change)) %>%
  group_by(regime) %>%
  summarise(
    avg_monthly_change   = mean(abs_rate_change, na.rm = TRUE),
    avg_active_move_size = mean(abs_rate_change[is_move == 1], na.rm = TRUE),
    pct_active_months    = mean(is_move, na.rm = TRUE) * 100,
    pct_hikes = mean(rate_change > 0.001, na.rm = TRUE) * 100,
    pct_cuts  = mean(rate_change < -0.001, na.rm = TRUE) * 100,
    .groups = "drop"
  )

episode_summary <- episode_stats %>%
  group_by(regime) %>%
  summarise(
    avg_spell_duration_months = mean(spell_length_months),
    max_spell_duration_months = max(spell_length_months),
    max_moves_in_single_spell = max(num_rate_moves),
    avg_net_change_per_spell  = mean(net_rate_change),
    max_cum_change_in_spell   = max(cum_abs_change),
    .groups = "drop"
  )

intro_stats_table <- monthly_summary %>%
  left_join(episode_summary, by = "regime") %>%
  pivot_longer(-regime, names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = regime, values_from = Value)

metric_labels <- c(
  avg_monthly_change        = "Avg. Monthly Change (Absolute)",
  avg_active_move_size      = "Avg. Active Move Size",
  pct_active_months         = "Frequency of Active Months (%)",
  pct_hikes                 = "Percentage of Hikes (%)",
  pct_cuts                  = "Percentage of Cuts (%)",
  avg_spell_duration_months = "Avg. Spell Duration (Months)",
  max_spell_duration_months = "Max. Spell Duration (Months)",
  max_moves_in_single_spell = "Max. Moves in Single Spell",
  avg_net_change_per_spell  = "Avg. Net Change per Spell",
  max_cum_change_in_spell   = "Max. Cumulative Change in Spell"
)

table3 <- intro_stats_table %>%
  mutate(Metric = metric_labels[Metric]) %>%
  arrange(match(Metric, metric_labels)) %>%
  kbl(format = "latex", booktabs = TRUE, align = "lcc", digits = 3,
      caption = "Monetary Policy Dynamics: High vs. Low Geopolitical Risk Regimes") %>%
  kable_styling(latex_options = c("hold_position"), font_size = 10) %>%
  pack_rows("Monthly / Action Level Metrics", 1, 5, latex_gap_space = "2.2pt") %>%
  pack_rows("Episode / Consecutive Period Metrics", 6, 10, latex_gap_space = "2.2pt") %>%
  footnote(general = "High and Low GPR regimes are defined based on the median of the Global Geopolitical Risk (GPR) Index. Active moves are defined as absolute rate changes strictly greater than zero.",
           threeparttable = TRUE)
save_tex(table3, out, "Table3_MP_Dynamics_GPR_Regimes")

#================================================
# Table 2: Macro variables, High vs Low GPR
#================================================
data_stats_macro <- final_data %>%
  arrange(Date) %>%
  mutate(
    gpr_median = median(GPR, na.rm = TRUE),
    regime     = ifelse(GPR > gpr_median, "High GPR", "Low GPR"),
    spell_id   = cumsum(regime != lag(regime, default = first(regime)))
  )

episode_stats_macro <- data_stats_macro %>%
  group_by(spell_id, regime) %>%
  summarise(
    spell_length_months = n(),
    cum_inflation = sum(df_HICP, na.rm = TRUE),
    cum_ip_growth = sum(df_IP, na.rm = TRUE),
    .groups = "drop"
  )

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

macro_stats_table <- bind_rows(monthly_summary_macro, episode_summary_macro) %>%
  pivot_wider(names_from = regime, values_from = Value) %>%
  arrange(Metric)

table2 <- macro_stats_table %>%
  mutate(Metric = sub("Inflation: |IP Growth: ", "", Metric)) %>%
  kbl(format = "latex", booktabs = TRUE, align = "lcc", digits = 3,
      caption = "Macroeconomic Dynamics: High vs. Low Geopolitical Risk Regimes") %>%
  kable_styling(latex_options = c("hold_position"), font_size = 10) %>%
  pack_rows("Inflation (HICP)", 7, 12, latex_gap_space = "2.2pt") %>%
  pack_rows("Industrial Production Growth", 1, 6, latex_gap_space = "2.2pt") %>%
  footnote(general = "High and Low GPR regimes are defined based on the median of the Global Geopolitical Risk (GPR) Index.",
           threeparttable = TRUE)
save_tex(table2, out, "Table2_Macro_Dynamics_GPR_Regimes")

#=======================================================
# Figure B5: IP, HICP and MRO with GPR regime shading
#=======================================================
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

regime_panel <- function(yvar, colour, title, ylab, xlab = NULL, show_x = FALSE) {
  ggplot() +
    geom_rect(data = shading_data,
              aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf, fill = regime),
              alpha = 0.25) +
    geom_vline(xintercept = event_lines, linetype = "dashed", color = "black") +
    geom_line(data = plot_data, aes(x = Date, y = .data[[yvar]]), color = colour, linewidth = 1) +
    scale_fill_manual(values = regime_colors) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = c(0, 0)) +
    theme_bw() +
    labs(title = title, y = ylab, x = xlab) +
    theme(
      legend.position = "none",
      axis.text.x  = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      panel.grid.major = element_line(color = "grey70")
    )
}

p_ip <- regime_panel("IP_Index", "forestgreen",
                     "Industrial Production (IP), Headline Inflation (HICP) Policy Rate, Euro Area",
                     "IP (Index)") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"))
p_hicp <- regime_panel("df_HICP", "mediumblue", "HICP (Headline Inflation)", "HICP (%)")
p_rate <- regime_panel("MRO", "black", "Monetary Policy Rate", "Policy Rate (%)",
                       xlab = "Year", show_x = TRUE)

save_plot(p_ip / p_hicp / p_rate, out, "FigureB5_Regime_Macro_Plot", width = 12, height = 8)
