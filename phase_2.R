# ==============================================================================
# phase_2.R
# ==============================================================================

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
      plot.caption = element_text(face = "italic", size = 10, hjust = 0, margin = margin(t = 15)),
      axis.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    )
}

valuation_year <- 2023
R              <- 65      
amount         <- 1       
proj_years     <- 100     
all_pf_ages    <- 30:90   

kannisto <- function(mhat, est.ages, proj.ages) {
  years <- rownames(mhat)
  mhat.proj <- matrix(NA, nrow = nrow(mhat), ncol = length(proj.ages), dimnames = list(rownames(mhat), proj.ages))
  for (t in years){
    mhat.est       <- mhat[as.character(t), as.character(est.ages)] 
    logit.mhat.est <- log(mhat.est / (1 - mhat.est)) 
    ols.est        <- lm(logit.mhat.est ~ est.ages) 
    logit.proj <- coef(ols.est)[1] + coef(ols.est)[2] * proj.ages
    mhat.proj[as.character(t),] <- exp(logit.proj) / (1 + exp(logit.proj))
  }
  cbind(mhat, mhat.proj)
}

LE_period <- function(mhat, Age, Year){
  ages  <- as.integer(colnames(mhat))
  mtx   <- mhat[as.character(Year), as.character(Age:max(ages))]
  term1 <- (1 - exp(-mtx[1])) / mtx[1]
  term2 <- 0
  for(k in 1:(length(mtx) - 1)){
    term2 <- term2 + (prod(exp(-mtx[1:k])) * ((1 - exp(-mtx[k+1])) / mtx[k+1]))
  }
  return(term1 + term2)
}

LE_cohort <- function(mhat_all, Age, Year){
  ages  <- as.integer(colnames(mhat_all))
  n_years_needed <- max(ages) - Age + 1
  target_years <- Year:(Year + n_years_needed - 1)
  mtx <- sapply(1:n_years_needed, function(k) mhat_all[as.character(target_years[k]), as.character(Age + k - 1)])
  term1 <- (1 - exp(-mtx[1])) / mtx[1]
  term2 <- 0
  for(k in 1:(length(mtx) - 1)){
    term2 <- term2 + (prod(exp(-mtx[1:k])) * ((1 - exp(-mtx[k+1])) / mtx[k+1]))
  }
  return(term1 + term2)
}

LA_period <- function(mhat, Age, Year, amount, rfr, R = 65){
  ages  <- as.integer(colnames(mhat))
  Tpxt <- cumprod(exp(-mhat[as.character(Year), as.character(Age:(max(ages)-1))]))
  n <- length(Tpxt)
  i_rates <- rfr$i[1:n] 
  i_rates[is.na(i_rates)] <- tail(rfr$i[!is.na(rfr$i)], 1) 
  v <- 1 / ((1 + i_rates)^(1:n))
  cf_0 <- ifelse(Age >= R, amount, 0)
  indicator <- ifelse(Age + (1:n) >= R, 1, 0)
  return(cf_0 + sum(amount * indicator * Tpxt * v))
}

LA_cohort <- function(mhat_all, Age, Year, amount, rfr, R = 65){
  ages  <- as.integer(colnames(mhat_all))
  n_years_needed <- max(ages) - Age
  target_years <- Year:(Year + n_years_needed - 1)
  mtx <- sapply(1:n_years_needed, function(k) mhat_all[as.character(target_years[k]), as.character(Age + k - 1)])
  Tpxt <- cumprod(exp(-mtx))
  n <- length(Tpxt)
  i_rates <- rfr$i[1:n] 
  i_rates[is.na(i_rates)] <- tail(rfr$i[!is.na(rfr$i)], 1) 
  v <- 1 / ((1 + i_rates)^(1:n))
  cf_0 <- ifelse(Age >= R, amount, 0)
  indicator <- ifelse(Age + (1:n) >= R, 1, 0)
  return(cf_0 + sum(amount * indicator * Tpxt * v))
}

simulate_vasicek_spot_curves <- function(n_paths, T_years, r0 = 0.03, theta = 0.03, k = 0.3, beta = 0.015, seed = 123) {
  set.seed(seed)
  dt <- 1       
  
  spot_paths <- list()
  for(j in 1:n_paths) {
    r_path <- numeric(T_years)
    r_path[1] <- r0
    
    for(i in 2:T_years) {
      r_path[i] <- r_path[i-1] + k * (theta - r_path[i-1]) * dt + beta * sqrt(dt) * rnorm(1, 0, 1)
    }
    
    spot_rates <- exp(-cumsum(r_path))^(-1/(1:T_years)) - 1
    spot_paths[[j]] <- data.frame(T = 1:T_years, i = spot_rates, Scenario = paste("Vasicek Path", j))
  }
  
  return(bind_rows(spot_paths))
}


nld_deaths <- read.table("./HMD Data/Deaths_1x1.txt", skip = 2, header = TRUE)
nld_exposures <- read.table("./HMD Data/Exposures_1x1.txt", skip = 2, header = TRUE)

nld_data <- nld_deaths %>%
  left_join(nld_exposures, by = c("Year", "Age")) %>%
  filter(Age != "110+") %>% mutate(Age = as.numeric(Age), mx = Total.x / Total.y)

Dxt <- nld_data %>% filter(Age %in% 0:90) %>% select(Year, Age, Total.x) %>% pivot_wider(names_from = Year, values_from = Total.x) %>% column_to_rownames(var = "Age") %>% as.matrix()
Ext <- nld_data %>% filter(Age %in% 0:90) %>% select(Year, Age, Total.y) %>% pivot_wider(names_from = Year, values_from = Total.y) %>% column_to_rownames(var = "Age") %>% as.matrix()

LCfit <- fit(lc(), Dxt = Dxt, Ext = Ext, ages = 0:90, years = as.numeric(colnames(Dxt)))
mhat_closed <- kannisto(mhat = t(fitted(LCfit, type = "rates")), est.ages = 80:90, proj.ages = 90:120)

mhat_future <- t(forecast(LCfit, h = 100)$rates)
mhat_all_closed <- rbind(mhat_closed, kannisto(mhat_future, 80:90, 90:120))
mhat_all_long <- rbind(mhat_closed, kannisto(mhat_future * 0.80, 80:90, 90:120))
mhat_all_mort <- rbind(mhat_closed, kannisto(mhat_future * 1.15, 80:90, 90:120))

rfr_files <- list.files(path = "./EIOPA Data", pattern = "\\.xlsx$", full.names = TRUE)
if(length(rfr_files) == 0) stop("ERROR: No EIOPA Excel files found!")

rfr_all_years <- bind_rows(lapply(rfr_files, function(file) {
  read_excel(file, sheet = "RFR_spot_with_VA", skip = 1) %>%
    select(T = 1, i = starts_with("Nether")) %>% 
    mutate(T = suppressWarnings(as.numeric(T)), i = suppressWarnings(as.numeric(i)), RFR_Year = str_extract(file, "\\d{4}")) %>% drop_na()
}))

rfr_curve_base <- rfr_all_years %>% filter(RFR_Year == "2026")

# Note: The paths visualized in the report were generated under an earlier, 
# unseeded global environment. The seed is set to ensure reproducibility of the code.

vasicek_curves <- simulate_vasicek_spot_curves(
  n_paths = 5, 
  T_years = proj_years, 
  r0 = 0.025, 
  theta = 0.035, 
  k = 0.15, 
  beta = 0.015, 
  seed = 1481
)

age_eval <- 30
df_survival <- data.frame(
  Attained_Age = age_eval:120,
  Period = c(1, cumprod(exp(-mhat_closed[as.character(valuation_year), as.character(age_eval:119)]))),
  Cohort = c(1, cumprod(exp(-sapply(1:(120-age_eval), function(k) mhat_all_closed[as.character(valuation_year+k-1), as.character(age_eval+k-1)]))))
) %>% pivot_longer(cols = c(Period, Cohort), names_to = "Method", values_to = "Survival_Prob")

plot_survival_curve <- ggplot(df_survival, aes(x = Attained_Age, y = Survival_Prob, color = Method, linetype = Method)) +
  geom_line(linewidth = 1.2) + scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("Period" = colors_issurance$primary, "Cohort" = colors_issurance$accent)) +
  scale_linetype_manual(values = c("Period" = "dashed", "Cohort" = "solid")) +
  theme_issurance() + labs(title = paste("Survival Probability Curve for a", age_eval, "-year-old"), x = "Attained Age", y = expression(paste("Probability of Survival (", " "[t], p[x], ")")))

df_le_gap <- data.frame(
  Age = all_pf_ages, 
  Period = sapply(all_pf_ages, function(a) LE_period(mhat_closed, a, valuation_year)), 
  Cohort = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_closed, a, valuation_year))
)

plot_le_compare <- ggplot(df_le_gap, aes(x = Age)) +
  geom_ribbon(aes(ymin = Period, ymax = Cohort), fill = colors_issurance$highlight, alpha = 0.2) +
  geom_line(aes(y = Period, color = "Period"), linewidth = 1.2) + geom_line(aes(y = Cohort, color = "Cohort"), linewidth = 1.2) +
  scale_color_manual(values = c("Period" = colors_issurance$primary, "Cohort" = colors_issurance$green_pf)) +
  theme_issurance() + theme(legend.position = "bottom") + labs(title = paste("Period vs. Cohort Life Expectancy in", valuation_year), x = "Participant Age", y = expression(e[x]))

df_le_time <- data.frame(
  Year = 1970:2023,
  Period = sapply(1970:2023, function(y) LE_period(mhat_closed, 65, y)),
  Cohort = sapply(1970:2023, function(y) LE_cohort(mhat_all_closed, 65, y))
)

plot_le_time <- ggplot(df_le_time, aes(x = Year)) +
  geom_ribbon(aes(ymin = Period, ymax = Cohort), fill = colors_issurance$highlight, alpha = 0.2) +
  geom_line(aes(y = Period, color = "Period"), linewidth = 1.2, linetype = "dashed") + geom_line(aes(y = Cohort, color = "Cohort"), linewidth = 1.2, linetype = "solid") +
  scale_color_manual(values = c("Period" = colors_issurance$primary, "Cohort" = colors_issurance$accent)) +
  theme_issurance() + theme(legend.position = "bottom") + labs(title = "Historical Evolution of Life Expectancy at Age 65", x = "Calendar Year", y = expression(e[65]))

df_epv_gap <- data.frame(
  Age = all_pf_ages, 
  Period = sapply(all_pf_ages, function(a) LA_period(mhat_closed, a, valuation_year, amount, rfr_curve_base, 65)), 
  Cohort = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_closed, a, valuation_year, amount, rfr_curve_base, 65))
)

plot_epv_compare <- ggplot(df_epv_gap, aes(x = Age)) +
  geom_ribbon(aes(ymin = Period, ymax = Cohort), fill = colors_issurance$highlight, alpha = 0.25) +
  geom_line(aes(y = Period, color = "Period"), linewidth = 1.2) + geom_line(aes(y = Cohort, color = "Cohort"), linewidth = 1.2) +
  scale_color_manual(values = c("Period" = colors_issurance$primary, "Cohort" = colors_issurance$green_pf)) +
  theme_issurance() + theme(legend.position = "bottom") + labs(title = "Period vs. Cohort EPV Valuation", x = "Participant Age", y = "Expected Present Value (EUR)")

df_error <- do.call(rbind, lapply(c(30, 50, 65), function(a) {
  period_vals <- sapply(1970:2023, function(y) LA_period(mhat_closed, a, y, amount, rfr_curve_base, 65))
  cohort_vals <- sapply(1970:2023, function(y) LA_cohort(mhat_all_closed, a, y, amount, rfr_curve_base, 65))
  data.frame(Year = 1970:2023, Age = as.factor(a), Error_Pct = (cohort_vals - period_vals) / period_vals)
}))

plot_epv_error <- ggplot(df_error, aes(x = Year, y = Error_Pct, color = Age)) +
  geom_line(linewidth = 1.2) + scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("30" = colors_issurance$green_pf, "50" = colors_issurance$primary, "65" = colors_issurance$accent)) +
  theme_issurance() + theme(legend.position = "bottom") + labs(title = "The Cost of Period Thinking: Valuation Deficit", x = "Calendar Year", y = "Under-Reserving Error (%)")

df_le_shocks <- data.frame(
  Age = all_pf_ages, 
  Base = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_closed, a, valuation_year)), 
  Long = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_long, a, valuation_year)), 
  Mort = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_mort, a, valuation_year))
)
plot_s2_le_corridor <- ggplot(df_le_shocks, aes(x = Age)) +
  geom_ribbon(aes(ymin = Mort, ymax = Long), fill = colors_issurance$highlight, alpha = 0.2) +
  geom_line(aes(y = Long, color = "Longevity Shock (0.8x)"), linewidth = 1.2) + geom_line(aes(y = Base, color = "Base Cohort"), linewidth = 1.5) + geom_line(aes(y = Mort, color = "Mortality Shock (1.15x)"), linewidth = 1.2) +
  scale_color_manual(values = c("Longevity Shock (0.8x)" = colors_issurance$accent, "Base Cohort" = colors_issurance$primary, "Mortality Shock (1.15x)" = colors_issurance$green_pf)) +
  theme_issurance() + labs(title = "Solvency II Life Expectancy Corridor", x = "Participant Age", y = "Expected Remaining Years")

df_mort_shocks <- data.frame(
  Age = all_pf_ages, 
  Base = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_closed, a, valuation_year, amount, rfr_curve_base, 65)), 
  Long = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_long, a, valuation_year, amount, rfr_curve_base, 65)), 
  Mort = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_mort, a, valuation_year, amount, rfr_curve_base, 65))
)
plot_s2_corridor <- ggplot(df_mort_shocks, aes(x = Age)) +
  geom_ribbon(aes(ymin = Mort, ymax = Long), fill = colors_issurance$highlight, alpha = 0.2) +
  geom_line(aes(y = Long, color = "Longevity Shock (0.8x)"), linewidth = 1.2) + geom_line(aes(y = Base, color = "Base Cohort"), linewidth = 1.5) + geom_line(aes(y = Mort, color = "Mortality Shock (1.15x)"), linewidth = 1.2) +
  scale_color_manual(values = c("Longevity Shock (0.8x)" = colors_issurance$accent, "Base Cohort" = colors_issurance$primary, "Mortality Shock (1.15x)" = colors_issurance$green_pf)) +
  theme_issurance() + labs(title = "Solvency II Liability Corridor", x = "Participant Age", y = "Expected Present Value")

plot_rfr_history <- ggplot(rfr_all_years, aes(x = T, y = i, color = RFR_Year)) +
  geom_line(linewidth = 1.2, alpha = 0.9) + scale_color_viridis_d(option = "magma", end = 0.85) + scale_y_continuous(labels = label_percent()) +
  theme_issurance() + theme(legend.position = "right") + labs(title = "EIOPA Risk-Free Rate Evolution", x = "Maturity (T) in Years", y = "Spot Interest Rate (i)")

df_rfr_impact <- expand_grid(Age = all_pf_ages, Target_Year = unique(rfr_all_years$RFR_Year)) %>% rowwise() %>%
  mutate(EPV = LA_cohort(mhat_all_closed, Age, valuation_year, amount, rfr_all_years %>% filter(RFR_Year == Target_Year), 65)) %>% ungroup()

plot_epv_rfr_shock <- ggplot(df_rfr_impact, aes(x = Age, y = EPV, color = Target_Year)) +
  geom_line(linewidth = 1.2, alpha = 0.9) + scale_color_viridis_d(option = "magma", end = 0.85) +
  theme_issurance() + theme(legend.position = "right") + labs(title = "Effect of Interest Rate Shock on EPV", x = "Participant Age", y = "Expected Present Value")

plot_vasicek_curves <- ggplot(vasicek_curves, aes(x = T, y = i, group = Scenario, color = Scenario)) +
  geom_line(linewidth = 1.2, alpha = 0.85) + scale_color_brewer(palette = "Set1") + scale_y_continuous(labels = label_percent()) +
  theme_issurance() + theme(legend.position = "right") + labs(title = "Vasicek Interest Rate Scenarios", x = "Maturity (T) in Years", y = "Simulated Spot Rate (i)")

df_vasicek_epv <- expand_grid(Age = all_pf_ages, Sim_Scenario = unique(vasicek_curves$Scenario)) %>% rowwise() %>%
  mutate(EPV = LA_cohort(mhat_all_closed, Age, valuation_year, amount, vasicek_curves %>% filter(Scenario == Sim_Scenario), 65)) %>% ungroup()

plot_vasicek_epv <- ggplot(df_vasicek_epv, aes(x = Age, y = EPV, group = Sim_Scenario, color = Sim_Scenario)) +
  geom_line(linewidth = 1.2, alpha = 0.85) + scale_color_brewer(palette = "Set1") +
  theme_issurance() + theme(legend.position = "right") + labs(title = "Stochastic Liability Distribution", x = "Participant Age", y = "Expected Present Value")

print(plot_survival_curve)
print(plot_le_compare)
print(plot_le_time)
print(plot_epv_compare)
print(plot_epv_error)
print(plot_s2_le_corridor)
print(plot_s2_corridor)
print(plot_rfr_history)
print(plot_epv_rfr_shock)
print(plot_vasicek_curves)
print(plot_vasicek_epv)