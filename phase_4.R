# ==============================================================================
# phase_4.R
# ==============================================================================

# 1. SETUP & CONSTANTS
library(tidyverse)
library(demography)
library(StMoMo) 
library(readxl)
library(scales)

colors_issurance <- list(
  primary   = "#2c3e50", accent    = "#e74c3c", highlight = "#2980b9", 
  ground    = "#7f8c8d", green_pf  = "#27ae60", silver_pf = "#bdc3c7", gold_pf   = "#f1c40f"
)

theme_issurance <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12),
      axis.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    )
}

valuation_year <- 2023
R              <- 65
amount         <- 1
proj_years     <- 100
future_years   <- valuation_year:(valuation_year + proj_years - 1)
portfolio_ages <- c(30, 40, 50, 60, 70, 80, 90)
n_sims         <- 1000

df_portfolios <- data.frame(
  Age = portfolio_ages,
  Green  = c(500, 1200, 2000, 1800, 1500, 300, 0),
  Silver = c(300, 850, 1400, 1800, 1650, 550, 50),
  Gold   = c(100, 500, 800, 1800, 1800, 800, 100)
) %>% pivot_longer(cols = c(Green, Silver, Gold), names_to = "Portfolio", values_to = "Participants")

# 2. ACTUARIAL FUNCTIONS
kannisto <- function(mhat, est.ages, proj.ages) {
  years <- rownames(mhat)
  mhat.proj <- matrix(NA, nrow = nrow(mhat), ncol = length(proj.ages), dimnames = list(rownames(mhat), proj.ages))
  for (t in years){
    mhat.est <- mhat[as.character(t), as.character(est.ages)] 
    ols.est  <- lm(log(mhat.est / (1 - mhat.est)) ~ est.ages) 
    logit.proj <- coef(ols.est)[1] + coef(ols.est)[2] * proj.ages
    mhat.proj[as.character(t),] <- exp(logit.proj) / (1 + exp(logit.proj))
  }
  cbind(mhat, mhat.proj)
}

LA_cohort <- function(mhat_all, Age, Year, amount, rfr, R = 65){
  n_years_needed <- max(as.integer(colnames(mhat_all))) - Age
  mtx <- sapply(1:n_years_needed, function(k) mhat_all[as.character(Year + k - 1), as.character(Age + k - 1)])
  Tpxt <- cumprod(exp(-mtx))
  n <- length(Tpxt)
  i_rates <- rfr$i[1:n] 
  i_rates[is.na(i_rates)] <- tail(rfr$i[!is.na(rfr$i)], 1) 
  cf_0 <- ifelse(Age >= R, amount, 0)
  return(cf_0 + sum(amount * ifelse(Age + (1:n) >= R, 1, 0) * Tpxt / ((1 + i_rates)^(1:n))))
}

expected_cashflow_cohort <- function(mhat_all, Age, Year, R = 65, amount = 1) {
  n_years_needed <- max(as.integer(colnames(mhat_all))) - Age
  mtx <- sapply(1:n_years_needed, function(k) mhat_all[as.character(Year + k - 1), as.character(Age + k - 1)])
  Tpxt <- cumprod(exp(-mtx))
  cf_0 <- ifelse(Age >= R, amount, 0)
  future_cf <- amount * ifelse(Age + (1:length(Tpxt)) >= R, 1, 0) * Tpxt
  return(c(cf_0, future_cf))
}

calc_portfolio_cashflows <- function(target_portfolio, mhat_matrix) {
  cf_matrix <- matrix(0, nrow = length(portfolio_ages), ncol = proj_years)
  for (i in seq_along(portfolio_ages)) {
    n_participants <- df_portfolios %>% filter(Portfolio == target_portfolio, Age == portfolio_ages[i]) %>% pull(Participants)
    if (n_participants > 0) {
      cf_single <- expected_cashflow_cohort(mhat_matrix, portfolio_ages[i], valuation_year)
      cf_single <- if(length(cf_single) < proj_years) c(cf_single, rep(0, proj_years - length(cf_single))) else cf_single[1:proj_years]
      cf_matrix[i, ] <- cf_single * n_participants
    }
  }
  return(colSums(cf_matrix))
}

calc_total_pv <- function(cf_vector, rfr_curve) {
  i_rates <- rfr_curve$i[1:(proj_years - 1)]
  i_rates[is.na(i_rates)] <- tail(rfr_curve$i[!is.na(rfr_curve$i)], 1)
  return(sum(cf_vector * c(1, 1 / ((1 + i_rates)^(1:(proj_years - 1))))))
}

# 3. DATA PREP & BASELINE FORECASTING
nld_data <- read.table("./HMD Data/Deaths_1x1.txt", skip=2, header=T) %>% left_join(read.table("./HMD Data/Exposures_1x1.txt", skip=2, header=T), by=c("Year", "Age")) %>% filter(Age != "110+") %>% mutate(Age=as.numeric(Age), mx=Total.x/Total.y)
Dxt <- nld_data %>% filter(Age %in% 0:90) %>% select(Year, Age, Total.x) %>% pivot_wider(names_from = Year, values_from = Total.x) %>% column_to_rownames("Age") %>% as.matrix()
Ext <- nld_data %>% filter(Age %in% 0:90) %>% select(Year, Age, Total.y) %>% pivot_wider(names_from = Year, values_from = Total.y) %>% column_to_rownames("Age") %>% as.matrix()

LCfit <- fit(lc(), Dxt = Dxt, Ext = Ext, ages = 0:90, years = as.numeric(colnames(Dxt)))
mhat_closed <- kannisto(t(fitted(LCfit, type = "rates")), 80:90, 90:120)
mhat_future <- t(forecast(LCfit, h = 100)$rates)

mhat_all_closed <- rbind(mhat_closed, kannisto(mhat_future, 80:90, 90:120))
mhat_all_long <- rbind(mhat_closed, kannisto(mhat_future * 0.80, 80:90, 90:120))
mhat_all_mort <- rbind(mhat_closed, kannisto(mhat_future * 1.15, 80:90, 90:120))

rfr_files <- list.files(path = "./EIOPA Data", pattern = "\\.xlsx$", full.names = TRUE)
if(length(rfr_files) == 0) stop("ERROR: No EIOPA Excel files found!")
rfr_curve_base <- bind_rows(lapply(rfr_files, function(file) { read_excel(file, sheet = "RFR_spot_with_VA", skip = 1) %>% select(T = 1, i = starts_with("Nether")) %>% mutate(T = as.numeric(T), i = as.numeric(i), RFR_Year = str_extract(file, "\\d{4}")) %>% drop_na() })) %>% filter(RFR_Year == "2026")

# 4. CASHFLOWS & PRESENT VALUES
cf_green_base  <- calc_portfolio_cashflows("Green", mhat_all_closed)
cf_silver_base <- calc_portfolio_cashflows("Silver", mhat_all_closed)
cf_gold_base   <- calc_portfolio_cashflows("Gold", mhat_all_closed)

pv_green_base  <- calc_total_pv(cf_green_base, rfr_curve_base)
pv_silver_base <- calc_total_pv(cf_silver_base, rfr_curve_base)
pv_gold_base   <- calc_total_pv(cf_gold_base, rfr_curve_base)

pv_green_long  <- calc_total_pv(calc_portfolio_cashflows("Green", mhat_all_long), rfr_curve_base)
pv_silver_long <- calc_total_pv(calc_portfolio_cashflows("Silver", mhat_all_long), rfr_curve_base)
pv_gold_long   <- calc_total_pv(calc_portfolio_cashflows("Gold", mhat_all_long), rfr_curve_base)

pv_green_mort  <- calc_total_pv(calc_portfolio_cashflows("Green", mhat_all_mort), rfr_curve_base)
pv_silver_mort <- calc_total_pv(calc_portfolio_cashflows("Silver", mhat_all_mort), rfr_curve_base)
pv_gold_mort   <- calc_total_pv(calc_portfolio_cashflows("Gold", mhat_all_mort), rfr_curve_base)

rfr_up   <- rfr_curve_base %>% mutate(i = i + 0.01)
rfr_down <- rfr_curve_base %>% mutate(i = i - 0.01)

pv_green_up    <- calc_total_pv(cf_green_base, rfr_up)
pv_silver_up   <- calc_total_pv(cf_silver_base, rfr_up)
pv_gold_up     <- calc_total_pv(cf_gold_base, rfr_up)

pv_green_down  <- calc_total_pv(cf_green_base, rfr_down)
pv_silver_down <- calc_total_pv(cf_silver_base, rfr_down)
pv_gold_down   <- calc_total_pv(cf_gold_base, rfr_down)

# 5. STOCHASTIC SIMULATIONS (For Uncertainty & VaR)
set.seed(2026) 
LC_sim <- simulate(LCfit, h = 100, nsim = n_sims)
mhat_sim_list <- lapply(1:n_sims, function(s) rbind(mhat_closed, kannisto(t(LC_sim$rates[,,s]), 80:90, 90:120)))

df_epv_all_sims <- expand_grid(Age = portfolio_ages, Simulation = 1:n_sims) %>%
  rowwise() %>% mutate(EPV = LA_cohort(mhat_sim_list[[Simulation]], Age, valuation_year, amount, rfr_curve_base, R)) %>% ungroup()

df_portfolio_sims <- df_epv_all_sims %>%
  inner_join(df_portfolios, by = "Age", relationship = "many-to-many") %>%
  mutate(Liability_Component = EPV * Participants) %>% group_by(Simulation, Portfolio) %>%
  summarize(Total_PV = sum(Liability_Component), .groups = 'drop')

# 6. PLOTTING
plot_sim_cf_multi <- ggplot() +
  geom_line(data = bind_rows(lapply(c(50, 70), function(age) bind_rows(lapply(1:100, function(s) {
    cf <- expected_cashflow_cohort(mhat_sim_list[[s]], age, valuation_year)
    data.frame(Year = future_years, Cashflow = if(length(cf) < proj_years) c(cf, rep(0, proj_years - length(cf))) else cf[1:proj_years], Simulation = s, Age = as.factor(age))
  })))), aes(x = Year, y = Cashflow, group = interaction(Simulation, Age), color = Age), alpha = 0.1) +
  geom_line(data = bind_rows(lapply(c(50, 70), function(age) {
    cf <- expected_cashflow_cohort(mhat_all_closed, age, valuation_year)
    data.frame(Year = future_years, Cashflow = if(length(cf) < proj_years) c(cf, rep(0, proj_years - length(cf))) else cf[1:proj_years], Age = as.factor(age))
  })), aes(x = Year, y = Cashflow, group = Age), color = "black", linewidth = 1.2) +
  scale_color_manual(values = c("50" = colors_issurance$highlight, "70" = colors_issurance$accent)) + theme_issurance() + theme(legend.position = "bottom") + labs(title = "Stochastic Cashflows: Active (Age 50) vs. Retired (Age 70)", x = "Year", y = "Expected Payout (EUR)", color = "Participant Age")

plot_portfolio_cf <- ggplot(data.frame(Year = rep(future_years, 3), Cashflow = c(cf_green_base, cf_silver_base, cf_gold_base), Portfolio = factor(rep(c("Green", "Silver", "Gold"), each = proj_years), levels = c("Green", "Silver", "Gold"))), aes(x = Year, y = Cashflow, color = Portfolio)) +
  geom_line(linewidth = 1.5) + scale_color_manual(values = c("Green" = colors_issurance$green_pf, "Silver" = colors_issurance$silver_pf, "Gold" = colors_issurance$gold_pf)) +
  scale_y_continuous(labels = comma) + theme_issurance() + theme(legend.position = "right") + labs(title = "Expected Portfolio Cashflows Over Time", x = "Year", y = "Expected Payout")

plot_portfolio_uncertainty <- ggplot(df_portfolio_sims, aes(x = Total_PV, fill = Portfolio, color = Portfolio)) +
  geom_density(alpha = 0.6) + scale_fill_manual(values = c("Green" = colors_issurance$green_pf, "Silver" = colors_issurance$silver_pf, "Gold" = colors_issurance$gold_pf)) +
  scale_color_manual(values = c("Green" = colors_issurance$green_pf, "Silver" = colors_issurance$silver_pf, "Gold" = colors_issurance$gold_pf)) +
  scale_x_continuous(labels = comma) + theme_issurance() + theme(legend.position = "bottom") + labs(title = "Uncertainty of Total Liabilities by Portfolio", x = "Total Present Value", y = "Density")

plot_stress_bio_pct <- ggplot(data.frame(Portfolio = factor(rep(c("Green", "Silver", "Gold"), each = 2), levels = c("Green", "Silver", "Gold")), Scenario = rep(c("Longevity Shock (0.8x)", "Mortality Shock (1.15x)"), times = 3), Impact_Pct = c((pv_green_long - pv_green_base)/pv_green_base, (pv_green_mort - pv_green_base)/pv_green_base, (pv_silver_long - pv_silver_base)/pv_silver_base, (pv_silver_mort - pv_silver_base)/pv_silver_base, (pv_gold_long - pv_gold_base)/pv_gold_base, (pv_gold_mort - pv_gold_base)/pv_gold_base)), aes(x = Portfolio, y = Impact_Pct, fill = Scenario)) +
  geom_hline(yintercept = 0, color = colors_issurance$primary, linewidth = 1) + geom_col(position = "dodge", color = colors_issurance$primary, alpha = 0.95) +
  geom_text(aes(label = scales::percent(Impact_Pct, accuracy = 0.1), vjust = ifelse(Impact_Pct > 0, -0.5, 1.5)), position = position_dodge(width = 0.9), fontface = "bold", color = colors_issurance$primary) +
  scale_fill_manual(values = c("Longevity Shock (0.8x)" = colors_issurance$accent, "Mortality Shock (1.15x)" = colors_issurance$green_pf)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.15, 0.15))) + theme_issurance() + theme(legend.position = "bottom") + labs(title = "Relative Sensitivity to Mortality Shocks", x = "Portfolio", y = "Impact on Total Liability (%)", fill = "Scenario")

plot_stress_ir_pct <- ggplot(data.frame(Portfolio = factor(rep(c("Green", "Silver", "Gold"), each = 2), levels = c("Green", "Silver", "Gold")), Scenario = rep(c("Rates UP (+100bps)", "Rates DOWN (-100bps)"), times = 3), Impact_Pct = c((pv_green_up - pv_green_base)/pv_green_base, (pv_green_down - pv_green_base)/pv_green_base, (pv_silver_up - pv_silver_base)/pv_silver_base, (pv_silver_down - pv_silver_base)/pv_silver_base, (pv_gold_up - pv_gold_base)/pv_gold_base, (pv_gold_down - pv_gold_base)/pv_gold_base)), aes(x = Portfolio, y = Impact_Pct, fill = Scenario)) +
  geom_hline(yintercept = 0, color = colors_issurance$primary, linewidth = 1) + geom_col(position = "dodge", color = colors_issurance$primary, alpha = 0.95) +
  geom_text(aes(label = scales::percent(Impact_Pct, accuracy = 0.1), vjust = ifelse(Impact_Pct > 0, -0.5, 1.5)), position = position_dodge(width = 0.9), fontface = "bold", color = colors_issurance$primary) +
  scale_fill_manual(values = c("Rates UP (+100bps)" = colors_issurance$highlight, "Rates DOWN (-100bps)" = colors_issurance$accent)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.15, 0.15))) + theme_issurance() + theme(legend.position = "bottom") + labs(title = "Relative Sensitivity to Interest Rate Shocks", x = "Portfolio", y = "Impact on Total Liability (%)", fill = "Scenario")

# 7. METRICS & OUTPUTS
portfolio_metrics <- df_portfolio_sims %>% group_by(Portfolio) %>% summarize(Best_Estimate = mean(Total_PV), VaR_98 = quantile(Total_PV, 0.98), .groups = 'drop')

cat("\n========================================================================\n")
cat("PHASE 4: PORTFOLIO LIABILITY METRICS (CRO REPORT)\n")
cat("========================================================================\n")
for(i in 1:nrow(portfolio_metrics)) {
  pf <- portfolio_metrics$Portfolio[i]
  cat(sprintf("%-10s Portfolio | Best Estimate: %11s EUR | 98%% VaR: %11s EUR\n", pf, format(round(portfolio_metrics$Best_Estimate[i], 0), big.mark = ","), format(round(portfolio_metrics$VaR_98[i], 0), big.mark = ",")))
}
cat("========================================================================\n\n")

print(plot_sim_cf_multi)
print(plot_portfolio_cf)
print(plot_portfolio_uncertainty)
print(plot_stress_bio_pct)
print(plot_stress_ir_pct)