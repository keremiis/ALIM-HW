# ==============================================================================
# phase_3.R
# ==============================================================================

# 1. SETUP & CONSTANTS
library(tidyverse)
library(demography)
library(StMoMo) 
library(readxl)
library(scales)
library(plotly)
library(htmlwidgets)

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
all_pf_ages    <- 30:90   
n_sims         <- 1000

# 2. ACTUARIAL FUNCTIONS
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

# 3. DATA PREPARATION
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
rfr_curve_base <- bind_rows(lapply(rfr_files, function(file) {
  read_excel(file, sheet = "RFR_spot_with_VA", skip = 1) %>%
    select(T = 1, i = starts_with("Nether")) %>% 
    mutate(T = suppressWarnings(as.numeric(T)), i = suppressWarnings(as.numeric(i)), RFR_Year = str_extract(file, "\\d{4}")) %>% drop_na()
})) %>% filter(RFR_Year == "2026")

# Calculate Deterministic Baselines
df_le_shocks <- data.frame(
  Age = all_pf_ages, 
  Base = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_closed, a, valuation_year)), 
  Long = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_long, a, valuation_year)), 
  Mort = sapply(all_pf_ages, function(a) LE_cohort(mhat_all_mort, a, valuation_year))
)
df_mort_shocks <- data.frame(
  Age = all_pf_ages, 
  Base = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_closed, a, valuation_year, amount, rfr_curve_base, R)), 
  Long = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_long, a, valuation_year, amount, rfr_curve_base, R)), 
  Mort = sapply(all_pf_ages, function(a) LA_cohort(mhat_all_mort, a, valuation_year, amount, rfr_curve_base, R))
)

# 4. STOCHASTIC SIMULATIONS
set.seed(2026) 
LC_sim <- simulate(LCfit, h = 100, nsim = n_sims)

df_kt_hist <- data.frame(Year = LCfit$years, kt = as.numeric(LCfit$kt), Type = "Historical")
sim_years <- (max(LCfit$years) + 1):(max(LCfit$years) + 100)

df_kt_sims <- bind_rows(lapply(1:n_sims, function(s) {
  data.frame(Year = sim_years, kt = as.numeric(LC_sim$kt.s[[1]][1, , s]), Simulation = s, Type = "Simulated")
}))

mhat_sim_list <- lapply(1:n_sims, function(s) {
  rbind(mhat_closed, kannisto(t(LC_sim$rates[,,s]), est.ages = 80:90, proj.ages = 90:120))
})

df_stoch_all <- expand_grid(Age = all_pf_ages, Simulation = 1:n_sims) %>%
  rowwise() %>%
  mutate(
    LE = LE_cohort(mhat_sim_list[[Simulation]], Age = Age, Year = valuation_year),
    EPV = LA_cohort(mhat_sim_list[[Simulation]], Age = Age, Year = valuation_year, amount = amount, rfr = rfr_curve_base, R = R)
  ) %>% ungroup()

# 5. PHASE 3 PLOTS
plot_sim_kt <- ggplot() +
  geom_line(data = df_kt_sims, aes(x = Year, y = kt, group = Simulation), color = colors_issurance$highlight, alpha = 0.03) +
  geom_line(data = df_kt_hist, aes(x = Year, y = kt), color = colors_issurance$primary, linewidth = 1.5) +
  theme_issurance() + labs(title = expression(paste("Stochastic Evolution of Mortality Trend (", k[t], ")")), subtitle = "1,000 simulations into the future", x = "Year", y = expression(k[t]))

ages_to_plot <- seq(30, 90, by = 5)
le_range <- seq(min(df_stoch_all$LE), max(df_stoch_all$LE), length.out = 200)
plot_sim_le_3d <- plot_ly()
for (current_age in ages_to_plot) {
  dens <- density(df_stoch_all$LE[df_stoch_all$Age == current_age], from = min(le_range), to = max(le_range), n = length(le_range))
  plot_sim_le_3d <- add_paths(plot_sim_le_3d, x = current_age, y = dens$x, z = dens$y, line = list(color = colors_issurance$primary, width = 3), showlegend = FALSE, hoverinfo = "none")
}
plot_sim_le_3d <- plot_sim_le_3d %>% layout(title = "Life Expectancy Distributions Across Ages (3D)", scene = list(xaxis = list(title = "Age", range = c(90, 30)), yaxis = list(title = "Life Expectancy"), zaxis = list(title = "Density", showticklabels = FALSE), camera = list(eye = list(x = -1.5, y = -1.5, z = 0.5))))


df_50 <- df_stoch_all %>% filter(Age == 50)
plot_sim_le_50 <- ggplot() +
  geom_histogram(data = df_50, aes(x = LE, y = after_stat(density)), bins = 30, fill = colors_issurance$silver_pf, color = "white") +
  geom_line(data = data.frame(LE = seq(min(df_50$LE)*0.98, max(df_50$LE)*1.02, length.out=200), Density = dnorm(seq(min(df_50$LE)*0.98, max(df_50$LE)*1.02, length.out=200), mean(df_50$LE), sd(df_50$LE))), aes(x = LE, y = Density), color = colors_issurance$primary, linewidth = 1.5, alpha = 0.8) +
  geom_vline(data = data.frame(Scenario = c("Shock (0.8x)", "Base", "Shock (1.15x)", "VaR (2% & 98%)", "VaR (2% & 98%)"), Value = c(df_le_shocks$Long[df_le_shocks$Age==50], df_le_shocks$Base[df_le_shocks$Age==50], df_le_shocks$Mort[df_le_shocks$Age==50], quantile(df_50$LE, 0.02), quantile(df_50$LE, 0.98))), aes(xintercept = Value, color = Scenario), linewidth = 1.5) +
  scale_color_manual(values = c("Shock (0.8x)" = colors_issurance$accent, "Base" = colors_issurance$highlight, "Shock (1.15x)" = colors_issurance$green_pf, "VaR (2% & 98%)" = colors_issurance$gold_pf)) +
  theme_issurance() + theme(legend.position = "bottom", legend.box = "horizontal") + guides(color = guide_legend(nrow = 1, title = NULL)) + labs(title = "Distribution Fitting: Age 50 Life Expectancy", x = "Expected Remaining Years", y = "Probability Density")

plot_sim_epv_50 <- ggplot() +
  geom_histogram(data = df_50, aes(x = EPV, y = after_stat(density)), bins = 30, fill = colors_issurance$silver_pf, color = "white") +
  geom_line(data = data.frame(EPV = seq(min(df_50$EPV)*0.98, max(df_50$EPV)*1.02, length.out=200), Density = dnorm(seq(min(df_50$EPV)*0.98, max(df_50$EPV)*1.02, length.out=200), mean(df_50$EPV), sd(df_50$EPV))), aes(x = EPV, y = Density), color = colors_issurance$primary, linewidth = 1.5, alpha = 0.8) +
  geom_vline(data = data.frame(Scenario = c("Shock (0.8x)", "Base", "Shock (1.15x)", "VaR (2% & 98%)", "VaR (2% & 98%)"), Value = c(df_mort_shocks$Long[df_mort_shocks$Age==50], df_mort_shocks$Base[df_mort_shocks$Age==50], df_mort_shocks$Mort[df_mort_shocks$Age==50], quantile(df_50$EPV, 0.02), quantile(df_50$EPV, 0.98))), aes(xintercept = Value, color = Scenario), linewidth = 1.5) +
  scale_color_manual(values = c("Shock (0.8x)" = colors_issurance$accent, "Base" = colors_issurance$highlight, "Shock (1.15x)" = colors_issurance$green_pf, "VaR (2% & 98%)" = colors_issurance$gold_pf)) +
  theme_issurance() + theme(legend.position = "bottom", legend.box = "horizontal") + guides(color = guide_legend(nrow = 1, title = NULL)) + labs(title = "Distribution Fitting: Age 50 EPV", x = "Expected Present Value", y = "Probability Density")


epv_range <- seq(min(df_stoch_all$EPV), max(df_stoch_all$EPV), length.out = 200)
plot_sim_epv_3d <- plot_ly()
for (current_age in ages_to_plot) {
  dens <- density(df_stoch_all$EPV[df_stoch_all$Age == current_age], from = min(epv_range), to = max(epv_range), n = length(epv_range))
  plot_sim_epv_3d <- add_paths(plot_sim_epv_3d, x = current_age, y = dens$x, z = dens$y, line = list(color = colors_issurance$highlight, width = 3), showlegend = FALSE)
}
plot_sim_epv_3d <- plot_sim_epv_3d %>% layout(title = "EPV Distributions Across Ages (3D)", scene = list(xaxis = list(title = "Age", range = c(90, 30)), yaxis = list(title = "EPV"), zaxis = list(title = "Density", showticklabels = FALSE), camera = list(eye = list(x = -1.5, y = -1.5, z = 0.5))))


df_var_all <- df_stoch_all %>% group_by(Age) %>% summarize(VaR_02 = quantile(EPV, 0.02), VaR_98 = quantile(EPV, 0.98), .groups = 'drop')
plot_adequacy <- ggplot(df_mort_shocks %>% inner_join(df_var_all, by = "Age") %>% mutate(Longevity = (Long - VaR_98) / VaR_98 * 100, Mortality = (Mort - VaR_02) / VaR_02 * 100) %>% select(Age, Longevity, Mortality) %>% pivot_longer(cols = c(Longevity, Mortality), names_to = "Tail", values_to = "Gap_Pct") %>% mutate(Tail = ifelse(Tail == "Longevity", "Right Tail: Longevity (0.8x vs 98% VaR)", "Left Tail: Early Mortality (1.15x vs 2% VaR)")), aes(x = Age, y = Gap_Pct, color = Tail)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + geom_line(linewidth = 1.5) + geom_point(size = 3) +
  scale_color_manual(values = c("Right Tail: Longevity (0.8x vs 98% VaR)" = colors_issurance$accent, "Left Tail: Early Mortality (1.15x vs 2% VaR)" = colors_issurance$green_pf)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) + theme_issurance() + theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "vertical") + labs(title = "Age-Dependent Adequacy of Regulatory Shocks", x = "Participant Age", y = "Adequacy Gap (%)")

print(plot_sim_kt)
print(plot_sim_le_50)
print(plot_sim_le_3d)
print(plot_sim_epv_50)
print(plot_sim_epv_3d)
print(plot_adequacy)
