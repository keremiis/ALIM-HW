# ==============================================================================
# app.R (Interactive Dashboard - Classic UI with Data-Rich Charts)
# ==============================================================================

# 1. SETUP & CONSTANTS
library(shiny)
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
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          plot.subtitle = element_text(hjust = 0.5, color = colors_issurance$accent, size = 13, face = "bold"),
          axis.title = element_text(face = "bold", size = 13),
          panel.grid.minor = element_blank())
}

valuation_year <- 2023
R              <- 65
amount         <- 1
all_pf_ages    <- 30:90
portfolio_ages <- seq(30, 90, by = 10) # For point labeling

# 2. ACTUARIAL FUNCTIONS
kannisto <- function(mhat, est.ages, proj.ages) {
  years <- rownames(mhat); mhat.proj <- matrix(NA, nrow = nrow(mhat), ncol = length(proj.ages), dimnames = list(rownames(mhat), proj.ages))
  for (t in years){
    mhat.est <- mhat[as.character(t), as.character(est.ages)] 
    ols.est  <- lm(log(mhat.est / (1 - mhat.est)) ~ est.ages) 
    logit.proj <- coef(ols.est)[1] + coef(ols.est)[2] * proj.ages
    mhat.proj[as.character(t),] <- exp(logit.proj) / (1 + exp(logit.proj))
  }
  cbind(mhat, mhat.proj)
}

LE_cohort <- function(mhat_all, Age, Year){
  n_y <- max(as.integer(colnames(mhat_all))) - Age + 1; target_y <- Year:(Year + n_y - 1)
  mtx <- sapply(1:n_y, function(k) mhat_all[as.character(target_y[k]), as.character(Age + k - 1)])
  term1 <- (1 - exp(-mtx[1])) / mtx[1]; term2 <- 0
  for(k in 1:(length(mtx) - 1)) term2 <- term2 + (prod(exp(-mtx[1:k])) * ((1 - exp(-mtx[k+1])) / mtx[k+1]))
  return(term1 + term2)
}

LA_cohort <- function(mhat_all, Age, Year, amount, rfr, R = 65){
  n_y <- max(as.integer(colnames(mhat_all))) - Age; target_y <- Year:(Year + n_y - 1)
  mtx <- sapply(1:n_y, function(k) mhat_all[as.character(target_y[k]), as.character(Age + k - 1)])
  Tpxt <- cumprod(exp(-mtx)); n <- length(Tpxt); i_r <- rfr$i[1:n]; i_r[is.na(i_r)] <- tail(rfr$i[!is.na(rfr$i)], 1)
  return(ifelse(Age >= R, amount, 0) + sum(amount * ifelse(Age + (1:n) >= R, 1, 0) * Tpxt / ((1 + i_r)^(1:n))))
}

# 3. DATA PREP & BASELINE FORECASTING (Runs once at startup)
cat("Loading data and compiling baseline models for the dashboard...\n")
nld_data <- read.table("./HMD Data/Deaths_1x1.txt", skip=2, header=T) %>% left_join(read.table("./HMD Data/Exposures_1x1.txt", skip=2, header=T), by=c("Year", "Age")) %>% filter(Age != "110+") %>% mutate(Age=as.numeric(Age), mx=Total.x/Total.y)
Dxt <- nld_data %>% filter(Age %in% 0:90) %>% select(Year, Age, Total.x) %>% pivot_wider(names_from = Year, values_from = Total.x) %>% column_to_rownames("Age") %>% as.matrix()
Ext <- nld_data %>% filter(Age %in% 0:90) %>% select(Year, Age, Total.y) %>% pivot_wider(names_from = Year, values_from = Total.y) %>% column_to_rownames("Age") %>% as.matrix()

LCfit <- fit(lc(), Dxt = Dxt, Ext = Ext, ages = 0:90, years = as.numeric(colnames(Dxt)))
mhat_closed <- kannisto(t(fitted(LCfit, type = "rates")), 80:90, 90:120)
mhat_future <- t(forecast(LCfit, h = 100)$rates)
mhat_all_closed <- rbind(mhat_closed, kannisto(mhat_future, 80:90, 90:120))

rfr_files <- list.files(path = "./EIOPA Data", pattern = "\\.xlsx$", full.names = TRUE)
rfr_curve_base <- bind_rows(lapply(rfr_files, function(file) { read_excel(file, sheet = "RFR_spot_with_VA", skip = 1) %>% select(T = 1, i = starts_with("Nether")) %>% mutate(T = as.numeric(T), i = as.numeric(i), RFR_Year = str_extract(file, "\\d{4}")) %>% drop_na() })) %>% filter(RFR_Year == "2026")

# ==============================================================================
# 4. SHINY APP UI & SERVER
# ==============================================================================

ui <- fluidPage(
  
  titlePanel("Dynamic Solvency II Mortality Shock Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Scenario Parameters"),
      p("Adjust the multiplier to see the real-time impact on biometric and financial liabilities."),
      
      sliderInput(inputId = "mort_mult", 
                  label = "Mortality Multiplier:", 
                  min = 0.50, max = 2.00, value = 1.00, step = 0.05),
      
      hr(),
      helpText(strong("Interpretation:")),
      helpText("Values < 1.0x : Longevity risk (participants live longer, higher EPV)."),
      helpText("Values > 1.0x : Mortality risk (participants die faster, lower EPV).")
    ),
    
    mainPanel(
      plotOutput("le_plot", height = "400px"),
      br(),
      plotOutput("epv_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  # Calculate Baseline values ONCE so the app is fast
  base_le <- sapply(all_pf_ages, function(a) LE_cohort(mhat_all_closed, a, valuation_year))
  base_epv <- sapply(all_pf_ages, function(a) LA_cohort(mhat_all_closed, a, valuation_year, amount, rfr_curve_base, R))
  
  # Reactive expression
  shocked_matrix <- reactive({
    mhat_future_shock <- mhat_future * input$mort_mult
    mhat_future_shock_closed <- kannisto(mhat_future_shock, 80:90, 90:120)
    rbind(mhat_closed, mhat_future_shock_closed)
  })
  
  output$le_plot <- renderPlot({
    shock_le <- sapply(all_pf_ages, function(a) LE_cohort(shocked_matrix(), a, valuation_year))
    df_le <- data.frame(Age = all_pf_ages, Base = base_le, Shock = shock_le)
    
    # Calculate real-time impact metric for Age 65
    impact_65 <- df_le$Shock[df_le$Age == 65] - df_le$Base[df_le$Age == 65]
    impact_text <- sprintf("Impact at Age 65: %+.2f Years", impact_65)
    
    ggplot(df_le, aes(x = Age)) +
      geom_ribbon(aes(ymin = pmin(Base, Shock), ymax = pmax(Base, Shock)), fill = colors_issurance$highlight, alpha = 0.2) +
      geom_line(aes(y = Base, color = "Base Scenario (1.0x)"), linewidth = 1.2, linetype = "dashed") +
      geom_line(aes(y = Shock, color = "Shocked Scenario"), linewidth = 1.5) +
      
      # Informative Data Points and Labels
      geom_point(data = filter(df_le, Age %in% portfolio_ages), aes(y = Shock, color = "Shocked Scenario"), size = 3) +
      geom_text(data = filter(df_le, Age %in% portfolio_ages), aes(y = Shock, label = round(Shock, 1)), 
                vjust = ifelse(input$mort_mult < 1, -1.2, 1.8), fontface = "bold", color = colors_issurance$primary) +
      
      scale_color_manual(values = c("Base Scenario (1.0x)" = colors_issurance$primary, "Shocked Scenario" = colors_issurance$accent)) +
      scale_x_continuous(breaks = seq(30, 90, by = 10)) +
      theme_issurance() + theme(legend.position = "bottom") +
      labs(title = paste("Cohort Life Expectancy (Multiplier:", input$mort_mult, "x)"), 
           subtitle = impact_text,
           x = "Participant Age", y = "Remaining Expected Years (e_x)", color = "")
  })
  
  output$epv_plot <- renderPlot({
    shock_epv <- sapply(all_pf_ages, function(a) LA_cohort(shocked_matrix(), a, valuation_year, amount, rfr_curve_base, R))
    df_epv <- data.frame(Age = all_pf_ages, Base = base_epv, Shock = shock_epv)
    
    # Calculate real-time impact metric for Age 65
    impact_pct_65 <- (df_epv$Shock[df_epv$Age == 65] - df_epv$Base[df_epv$Age == 65]) / df_epv$Base[df_epv$Age == 65] * 100
    impact_text <- sprintf("Liability Impact at Age 65: %+.2f%%", impact_pct_65)
    
    ggplot(df_epv, aes(x = Age)) +
      geom_ribbon(aes(ymin = pmin(Base, Shock), ymax = pmax(Base, Shock)), fill = colors_issurance$highlight, alpha = 0.2) +
      geom_line(aes(y = Base, color = "Base Scenario (1.0x)"), linewidth = 1.2, linetype = "dashed") +
      geom_line(aes(y = Shock, color = "Shocked Scenario"), linewidth = 1.5) +
      
      # Informative Data Points and Labels
      geom_point(data = filter(df_epv, Age %in% portfolio_ages), aes(y = Shock, color = "Shocked Scenario"), size = 3) +
      geom_text(data = filter(df_epv, Age %in% portfolio_ages), aes(y = Shock, label = round(Shock, 2)), 
                vjust = ifelse(input$mort_mult < 1, -1.2, 1.8), fontface = "bold", color = colors_issurance$primary) +
      
      scale_color_manual(values = c("Base Scenario (1.0x)" = colors_issurance$primary, "Shocked Scenario" = colors_issurance$green_pf)) +
      scale_x_continuous(breaks = seq(30, 90, by = 10)) +
      theme_issurance() + theme(legend.position = "bottom") +
      labs(title = paste("Cohort Expected Present Value (Multiplier:", input$mort_mult, "x)"), 
           subtitle = impact_text,
           x = "Participant Age", y = "EPV (EUR)", color = "")
  })
}

shinyApp(ui = ui, server = server)