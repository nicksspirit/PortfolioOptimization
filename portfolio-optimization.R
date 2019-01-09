library(shiny)
library(GenSA)
library(ggplot2)
library(GA)
library(magrittr)
library(quantmod)
library(BatchGetSymbols)
library(PerformanceAnalytics)
library(tagsinput)
library(knitr)
library(DEoptim)

css <- "
#plot-container {
  position: relative;
}

#loading-spinner {
  position: absolute;
  left: 50%;
  top: 50%;
  z-index: -1;
  margin-top: -33px;  /* half of the spinner's height */
  margin-left: -33px; /* half of the spinner's width */
}

#plot.recalculating {
  z-index: -2;
}
"

# Downloads price data from yahoo finance for stock tickers
get.ticker_symbols <- function(tickers, first_date, last_date, freq_data) {
    return (
      BatchGetSymbols(
        tickers = tickers,
        first.date = first_date,
        last.date = last_date,
        freq.data = freq_data,
        cache.folder = file.path(tempdir(), "BGS_Cache")
      )
    )
  }

stock_percentage <- function(weight_vals) {
  return (
    weight_vals %>%
    divide_by(sum(weight_vals)) %>%
    multiply_by(100) %>%
    round(2)
  )
}

get.closing_prices <- function(tickers_symbols, ticker_data) {
  df <- as.data.frame(ticker_data$df.tickers)
  available_tickers <- c()
  cl_prices <- NULL

  for (symbol in tickers_symbols) {
    ticker_cl_price <- df %>%
      filter(ticker == symbol) %>%
      select(price.close) %>%
      use_series("price.close")

    if (length(ticker_cl_price) != 0) {
      available_tickers <- c(available_tickers, symbol)
      cl_prices <- cbind(cl_prices, ticker_cl_price)
    }
  }

  colnames(cl_prices) <- available_tickers

  return(cl_prices)
}

get_obj_fn <-  function(cl_prices, threshold = 0.225) {
  log_returns <- diff(log(cl_prices))
  mu <- colMeans(log_returns) # mean of the log-returns series
  sigma <- cov(log_returns) # covariance of the log-returns series

  obj_fn <- function(weights) {
    if (sum(weights) == 0) {
      weights <- weights + 1e-2
    }

    weights <- weights / sum(weights)

    # ES (expected shortfall) computes the conditional value-at-risk
    CVaR <- ES(
      weights = weights,
      method = "gaussian",
      portfolio_method = "component",
      mu = mu,
      sigma = sigma
    )

    return(CVaR$ES + 1e3 * max(CVaR$pct_contrib_ES - threshold, 0))
  }

  return(obj_fn)
}

ui <- fluidPage(
  titlePanel("Optimize Your Portfolio"),
  br(),
  helpText(
    "Mean-risk models were developed in the 1950s for portfolio selection problems. Value-at-Risk (VaR) and Conditional Value-at-Risk (CVaR) are the most popular measures of downside risk.",
    "VaR is the negative value of the portfolio return such that lower returns will only occur with at most a preset probability level, which typically is between one and five percent.",
    "CVaR is the negative value of the mean of all return realizations that are below the VaR."
  ),
  helpText(
    "This app was created to analyze the stock data from a date range and then pick an optimization strategy in order to find the portfolio weights",
    "for which the portfolio has the lowest CVaR and each investment can contribute at most 22.5% or (one set by you) to total portfolio CVaR."
  ),
  sidebarLayout(
    sidebarPanel(
      tagsTextInput(
        inputId = "tickers",
        label = "Add your stock tickers",
        value = "GOOG"
      ),
      dateRangeInput(
        inputId = "date_range",
        label = "Pick a date range",
        start = "2012-10-01",
        end = "2018-10-31"
      ),
      radioButtons(
        inputId = "stock_freq",
        label = "Choose the frequency",
        choices = c("Monthly" = "monthly", "Daily" = "daily")
      ),
      radioButtons(
        inputId = "method",
        label = "Choose a search strategy",
        choices = c(
          "Simulated Anealing" = "SA",
          "Genetic Algorithm" = "GA",
          "Differential Optimization" = "DE"
        )
      ),
      numericInput(
        inputId = "threshold",
        label = "Choose an upper percentage CVaR constraint.",
        value = 0.225,
        min = 0,
        max = 1,
        step = 0.1
      ),
      actionButton(inputId = "generate", label = "Generate")
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.generate",
        tags$head(tags$style(HTML(css))),
        br(),
        div(
          id = "plot-container",
          tags$img(src = "spinner.gif", id = "loading-spinner"),
          conditionalPanel(
            condition = "output.downloaded",
            h4("Stock Download Summary"),
            tableOutput(outputId = "dl_output")
          ),
          br(),
          conditionalPanel(
            condition = "output.plotted",
            h4("Closing Prices"),
            p("A plot of the closing prices given the frequency of the data."),
            plotOutput(outputId = "cl_plot")
          ),
          br(),
          conditionalPanel(
            condition = "output.solved",
            h4("Percentage of Stocks"),
            uiOutput("CVaR"),
            p("The table below computes the percentage CVaR contributions for each asset, showing how much each asset contributes to the total portfolio CVaR."),
            tableOutput(outputId = "search_results")
          )
        )
      )
    )
  )
)



server <- function(input, output) {
  ticker_data <- eventReactive(input$generate, {
    symbols <- toupper(strsplit(input$tickers, ",")[[1]])

    # Downloaded ticker data
    dl_tickers <- get.ticker_symbols(
      tickers = symbols,
      first_date = input$date_range[1],
      last_date = input$date_range[2],
      freq_data = input$stock_freq
    )

    closing_prices <- get.closing_prices(symbols, dl_tickers)

    num_tickers <- ncol(closing_prices)

    set.seed(1234)

    if (input$method == "SA") {
      results <- GenSA(
        fn = get_obj_fn(closing_prices, input$threshold),
        lower = rep(0, num_tickers),
        upper = rep(1, num_tickers),
        control = list(smooth = FALSE, max.call = 3000)
      )
    }
    else if (input$method == "DE") {
      results <- DEoptim(
        fn = get_obj_fn(closing_prices, input$threshold),
        lower = rep(0, num_tickers),
        upper = rep(1, num_tickers)
      )
    }
    else {
      obj_fn <- get_obj_fn(closing_prices, input$threshold)

      results <- ga(
        type = "real-valued",
        fitness = obj_fn,
        lower = rep(0, num_tickers),
        upper = rep(1, num_tickers),
        popSize = 50
      )
    }

    list(
      symbols = colnames(closing_prices),
      first_date = input$date_range[1],
      last_date = input$date_range[2],
      stock_freq = input$stock_freq,
      dl_tickers = dl_tickers,
      search_type = input$method,
      search_method = results
    )
  })

  # Stock Download Summary

  output$downloaded <- reactive({
    dl_summary <- ticker_data()$dl_tickers$df.control
    nrow(dl_summary)
  })

  output$dl_output <- renderTable({
    dl_summary <- ticker_data()$dl_tickers$df.control

    dl_summary %>% select(-total.obs, -perc.benchmark.dates)
  })

  outputOptions(output, "downloaded", suspendWhenHidden = FALSE)

  # Closing Price of Stocks

  output$plotted <- reactive({
    tickers <- ticker_data()$dl_tickers$df.tickers
    nrow(tickers)
  })

  output$cl_plot <- renderPlot({
    tickers <- ticker_data()$dl_tickers$df.tickers

    ggplot(tickers, aes(x = ref.date, y = price.close)) +
      geom_line() +
      facet_wrap( ~ ticker, scales = "free_y")
  })

  outputOptions(output, "plotted", suspendWhenHidden = FALSE)

  # Table to show how much you should invest in each stock

  output$solved <- reactive({
    ticker_data <- ticker_data()
    length(ticker_data$symbols)
  })

  output$search_results <- renderTable({
    ticker_data <- ticker_data()
    model <- ticker_data$search_method

    if (ticker_data$search_type == "SA") {

      output$CVaR <- renderUI({
        h5(sprintf("The CVaR is %g", model$value))
      })

      stock_df <- data.frame(
        Ticker = ticker_data$symbols,
        Weights = model$par,
        Percentage = stock_percentage(model$par)
      )
    }
    else if (input$method == "DE") {
      output$CVaR <- renderUI({
        h5(sprintf("The CVaR is %g", model$optim$bestval))
      })

      stock_df <-  data.frame(
        Ticker = ticker_data$symbols,
        Weights = c(model$optim$bestmem),
        Percentage = c(stock_percentage(model$optim$bestmem))
      )
    }
    else {

      if (length(ticker_data$symbols) == 1) {
        solution <- tail(model@solution, 1)
      }
      else {
        solution <- unclass(model@solution)
      }

      output$CVaR <- renderUI({
        h5(sprintf("The CVaR is %g", model@fitnessValue))
      })

      stock_df <-  data.frame(
        Ticker = ticker_data$symbols,
        Weights = c(solution),
        Percentage = c(stock_percentage(solution))
      )
    }


    stock_df
  })

  outputOptions(output, "solved", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)
