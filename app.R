#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

### BLACK SCHOLES SHINY APP ###
rm(list = ls())

#load packages
library(tidyverse)
library(stargazer)
library(shiny)
library(ggplot2)

### Black-Scholes Option Pricer ###

# Black-Scholes with Greeks 
bs_with_greeks <- function(S, K, T, r, sigma, type = "call") {
  d1 <- (log(S / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  price <- if (type == "call") {
    S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else {
    K * exp(-r * T) * pnorm(-d2) - S * pnorm(-d1)
  }
  
  delta <- if (type == "call") pnorm(d1) else -pnorm(-d1)
  gamma <- dnorm(d1) / (S * sigma * sqrt(T))
  vega <- S * dnorm(d1) * sqrt(T)
  theta <- if (type == "call") {
    - (S * dnorm(d1) * sigma) / (2 * sqrt(T)) - r * K * exp(-r * T) * pnorm(d2)
  } else {
    - (S * dnorm(d1) * sigma) / (2 * sqrt(T)) + r * K * exp(-r * T) * pnorm(-d2)
  }
  
  return(list(price = price, delta = delta, gamma = gamma, vega = vega, theta = theta))
}

# GBM Simulation 
simulate_gbm <- function(S0, mu, sigma, T, steps) {
  dt <- T / steps
  time <- seq(0, T, by = dt)
  n <- length(time)
  Wt <- cumsum(c(0, sqrt(dt) * rnorm(n - 1)))  # Brownian motion
  S <- S0 * exp((mu - 0.5 * sigma^2) * time + sigma * Wt)
  return(data.frame(time = time, S = S))
}

# PnL & Greeks over time 
calc_pnl <- function(asset_path, K, r, sigma, T, type = "call", market_price_spread = 0.02) {
  n <- nrow(asset_path)
  results <- data.frame(
    time = asset_path$time,
    S = asset_path$S,
    price = numeric(n),
    pnl = numeric(n),
    delta = numeric(n),
    gamma = numeric(n),
    vega = numeric(n),
    theta = numeric(n)
  )
  
  init <- bs_with_greeks(S = asset_path$S[1], K, T, r, sigma, type)
  market_price <- init$price * (1 + market_price_spread)
  
  for (i in 1:n) {
    t_remaining <- T - asset_path$time[i]
    if (t_remaining <= 0) {
      results$price[i] <- max(ifelse(type == "call", asset_path$S[i] - K, K - asset_path$S[i]), 0)
    } else {
      bs <- bs_with_greeks(asset_path$S[i], K, t_remaining, r, sigma, type)
      results$price[i] <- bs$price
      results$delta[i] <- bs$delta
      results$gamma[i] <- bs$gamma
      results$vega[i] <- bs$vega
      results$theta[i] <- bs$theta
    }
    results$pnl[i] <- results$price[i] - market_price
  }
  
  return(results)
}

### ui 
ui <- fluidPage(
  titlePanel("Black-Scholes Option PnL Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("S0", "Initial Spot Price (S₀)", 100),
      numericInput("K", "Strike Price (K)", 100),
      numericInput("T", "Time to Maturity (T, in years)", 1),
      numericInput("r", "Risk-Free Rate (r)", 0.05),
      numericInput("sigma", "Volatility (σ)", 0.2),
      numericInput("mu", "Expected Return (μ)", 0.05),
      numericInput("steps", "Simulation Steps", 100),
      selectInput("type", "Option Type", c("call", "put")),
      numericInput("spread", "Market Price Spread (%)", 2, min = 0, step = 0.1)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Price Path", plotOutput("pricePlot")),
        tabPanel("PnL", plotOutput("pnlPlot")),
        tabPanel("Greeks",
                 plotOutput("deltaPlot"),
                 plotOutput("gammaPlot"),
                 plotOutput("vegaPlot"),
                 plotOutput("thetaPlot"))
      )
    )
  )
)

### server
server <- function(input, output) {
  sim_data <- reactive({
    simulate_gbm(
      S0 = input$S0,
      mu = input$mu,
      sigma = input$sigma,
      T = input$T,
      steps = input$steps
    )
  })
  
  pnl_data <- reactive({
    calc_pnl(
      asset_path = sim_data(),
      K = input$K,
      r = input$r,
      sigma = input$sigma,
      T = input$T,
      type = input$type,
      market_price_spread = input$spread / 100
    )
  })
  
  output$pricePlot <- renderPlot({
    ggplot(sim_data(), aes(x = time, y = S)) +
      geom_line(color = "darkgreen") +
      labs(title = "Simulated Asset Price Path", y = "Price", x = "Time") +
      theme_minimal()
  })
  
  output$pnlPlot <- renderPlot({
    ggplot(pnl_data(), aes(x = time, y = pnl)) +
      geom_line(color = "steelblue") +
      labs(title = "Option PnL Over Time", y = "PnL (£)", x = "Time") +
      theme_minimal()
  })
  
  output$deltaPlot <- renderPlot({
    ggplot(pnl_data(), aes(x = time, y = delta)) +
      geom_line(color = "purple") +
      labs(title = "Delta", y = "Delta", x = "Time") +
      theme_minimal()
  })
  
  output$gammaPlot <- renderPlot({
    ggplot(pnl_data(), aes(x = time, y = gamma)) +
      geom_line(color = "red") +
      labs(title = "Gamma", y = "Gamma", x = "Time") +
      theme_minimal()
  })
  
  output$vegaPlot <- renderPlot({
    ggplot(pnl_data(), aes(x = time, y = vega)) +
      geom_line(color = "orange") +
      labs(title = "Vega", y = "Vega", x = "Time") +
      theme_minimal()
  })
  
  output$thetaPlot <- renderPlot({
    ggplot(pnl_data(), aes(x = time, y = theta)) +
      geom_line(color = "brown") +
      labs(title = "Theta", y = "Theta", x = "Time") +
      theme_minimal()
  })
}


# Launch App 
shinyApp(ui = ui, server = server)
