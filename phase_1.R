# ==============================================================================
# phase_1.R
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
      plot.caption = element_text(face = "italic", size = 10, hjust = 0, margin = margin(t = 15)),
      axis.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    )
}

valuation_year <- 2023
R              <- 65
amount         <- 1
all_pf_ages    <- 30:90

# 2. ACTUARIAL FUNCTIONS
kannisto <- function(mhat, est.ages, proj.ages) {
  years <- rownames(mhat)
  mhat.proj <- matrix(NA, nrow = nrow(mhat), ncol = length(proj.ages), 
                      dimnames = list(rownames(mhat), proj.ages))
  for (t in years){
    mhat.est       <- mhat[as.character(t), as.character(est.ages)] 
    logit.mhat.est <- log(mhat.est / (1 - mhat.est)) 
    ols.est        <- lm(logit.mhat.est ~ est.ages) 
    phi1 <- coef(ols.est)[1]
    phi2 <- coef(ols.est)[2]
    logit.proj <- phi1 + phi2 * proj.ages
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
    t1 <- prod(exp(-mtx[1:k])) 
    t2 <- (1 - exp(-mtx[k+1])) / mtx[k+1] 
    term2 <- term2 + (t1 * t2)
  }
  return(term1 + term2)
}

LA_period <- function(mhat, Age, Year, amount, rfr, R = 65){
  ages  <- as.integer(colnames(mhat))
  pxt <- exp(-mhat[as.character(Year), as.character(Age:(max(ages)-1))])
  Tpxt <- cumprod(pxt)
  n <- length(Tpxt)
  
  i_rates <- rfr$i[1:n] 
  i_rates[is.na(i_rates)] <- tail(rfr$i[!is.na(rfr$i)], 1) 
  v <- 1 / ((1 + i_rates)^(1:n))
  
  cf_0 <- ifelse(Age >= R, amount, 0)
  future_ages <- Age + (1:n)
  indicator <- ifelse(future_ages >= R, 1, 0)
  return(cf_0 + sum(amount * indicator * Tpxt * v))
}

# 3. DATA PREPARATION
nld_deaths <- read.table("./HMD Data/Deaths_1x1.txt", skip = 2, header = TRUE)
nld_exposures <- read.table("./HMD Data/Exposures_1x1.txt", skip = 2, header = TRUE)

nld_data <- nld_deaths %>%
  left_join(nld_exposures, by = c("Year", "Age")) %>%
  filter(Age != "110+") %>%
  mutate(Age = as.numeric(Age), mx = Total.x / Total.y)

est_ages <- 0:90

Dxt <- nld_data %>%
  filter(Age %in% est_ages) %>%
  select(Year, Age, Total.x) %>%
  pivot_wider(names_from = Year, values_from = Total.x) %>%
  column_to_rownames(var = "Age") %>% as.matrix()

Ext <- nld_data %>%
  filter(Age %in% est_ages) %>%
  select(Year, Age, Total.y) %>%
  pivot_wider(names_from = Year, values_from = Total.y) %>%
  column_to_rownames(var = "Age") %>% as.matrix()

LC_model <- lc() 
LCfit <- fit(LC_model, Dxt = Dxt, Ext = Ext, ages = est_ages, years = as.numeric(colnames(Dxt)))

mhat <- t(fitted(LCfit, type = "rates"))
mhat_closed <- kannisto(mhat = mhat, est.ages = 80:90, proj.ages = 90:120)

rfr_files <- list.files(path = "./EIOPA Data", pattern = "\\.xlsx$", full.names = TRUE)
if(length(rfr_files) == 0) stop("ERROR: No EIOPA Excel files found.")

rfr_list <- lapply(rfr_files, function(file) {
  df <- read_excel(file, sheet = "RFR_spot_with_VA", skip = 1) %>%
    select(T = 1, i = starts_with("Nether")) %>% 
    mutate(T = suppressWarnings(as.numeric(T)), i = suppressWarnings(as.numeric(i)), RFR_Year = str_extract(file, "\\d{4}")) %>%
    drop_na()
  return(df)
})
rfr_all_years <- bind_rows(rfr_list)
rfr_curve_base <- rfr_all_years %>% filter(RFR_Year == "2026")

# 4. PHASE 1 PLOTS
plot_log_mort <- nld_data %>% filter(Age %in% c(30, 40, 50, 60, 70, 80, 90)) %>%
  ggplot(aes(x = Year, y = log(mx), color = as.factor(Age))) +
  geom_line(linewidth = 1.2, alpha = 0.8) + scale_color_viridis_d(option = "plasma", end = 0.8) +
  theme_issurance() + labs(title = "Historical Evolution of Log Mortality Rates", x = "Year", y = expression(paste("Log( ", m[x], " )")), color = "Age")

plot_ax <- ggplot(data.frame(Age = LCfit$ages, ax = as.numeric(LCfit$ax)), aes(x = Age, y = ax)) +
  geom_line(color = colors_issurance$primary, linewidth = 1.2) + theme_issurance() + labs(title = expression(paste("Baseline Mortality Profile (", a[x], ")")), x = "Age", y = expression(a[x]))

plot_bx <- ggplot(data.frame(Age = LCfit$ages, bx = as.numeric(LCfit$bx)), aes(x = Age, y = bx)) +
  geom_line(color = colors_issurance$accent, linewidth = 1.2) + theme_issurance() + labs(title = expression(paste("Improvement Sensitivity (", b[x], ")")), x = "Age", y = expression(b[x]))

plot_kt <- ggplot(data.frame(Year = LCfit$years, kt = as.numeric(LCfit$kt)), aes(x = Year, y = kt)) +
  geom_line(color = colors_issurance$highlight, linewidth = 1.2) + theme_issurance() + labs(title = expression(paste("Mortality Trend Over Time (", k[t], ")")), x = "Year", y = expression(k[t]))

df_combined <- bind_rows(
  data.frame(Age = 0:90, mx = mhat[as.character(valuation_year), ], Type = "Lee-Carter Fit (0-90)"),
  data.frame(Age = 90:120, mx = mhat_closed[as.character(valuation_year), as.character(90:120)], Type = "Kannisto Tail (90-120)"),
  nld_data %>% filter(Year == valuation_year, Age >= 90) %>% select(Age, mx) %>% mutate(Type = "Raw Data")
)

plot_full <- ggplot(df_combined, aes(x = Age, y = mx, color = Type, linetype = Type)) +
  geom_line(data = filter(df_combined, Type != "Raw Data"), linewidth = 1.1) +
  geom_point(data = filter(df_combined, Type == "Raw Data"), size = 2, alpha = 0.5) +
  scale_y_log10() + theme_issurance() + theme(legend.position = "bottom") + labs(title = "Full Mortality Profile (Age 0-120)", x = "Age", y = expression(paste("Log( ", m[x], " )")))

plot_zoomed <- plot_full + coord_cartesian(xlim = c(80, 120), ylim = c(0.02, 1.0)) + labs(title = "The Kannisto Divergence (Age 80-120)")

df_multi_le <- do.call(rbind, lapply(c(35, 45, 55, 65), function(a) {
  data.frame(Year = as.integer(rownames(mhat_closed)), LE = sapply(as.integer(rownames(mhat_closed)), function(y) LE_period(mhat_closed, a, y)), Age = as.factor(a))
}))

plot_multi_le <- ggplot(df_multi_le, aes(x = Year, y = LE, color = Age)) + geom_line(linewidth = 1.1) + theme_issurance() + labs(title = "Historical Evolution of Period Life Expectancy", x = "Year", y = expression(e[x]))

plot_rfr <- ggplot(rfr_curve_base, aes(x = T, y = i)) + geom_line(color = colors_issurance$highlight, linewidth = 1.2) + scale_y_continuous(labels = label_percent()) + theme_issurance() + labs(title = "EIOPA Risk-Free Rate Term Structure", x = "Maturity (T) in Years", y = "Interest Rate (i)")

df_epv_all <- data.frame(Age = all_pf_ages, EPV = sapply(all_pf_ages, function(x) LA_period(mhat_closed, x, valuation_year, amount, rfr_curve_base, R)))

df_epv_decades <- subset(df_epv_all, Age %% 10 == 0)

plot_epv_line <- ggplot(df_epv_all, aes(x = Age, y = EPV)) + 
  geom_line(color = colors_issurance$primary, linewidth = 1.2) +
  geom_point(data = df_epv_decades, color = colors_issurance$accent, size = 4) +
  geom_text(data = df_epv_decades, aes(label = round(EPV, 2)), 
            vjust = -1.2, fontface = "bold", color = colors_issurance$accent, size = 4.5) +
  scale_x_continuous(breaks = seq(30, 90, by = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + 
  theme_issurance() + 
  labs(
    title = "Expected Present Value (EPV) Curve",
    x = "Participant Age", 
    y = "Expected Present Value"
  )

print(plot_log_mort)
print(plot_ax)
print(plot_bx)
print(plot_kt)
print(plot_full)
print(plot_zoomed)
print(plot_multi_le)
print(plot_rfr)
print(plot_epv_line)
