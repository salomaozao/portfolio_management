library(shiny)
library(bslib)
library(tidyquant)
library(tidyverse)
library(slider)
library(jsonlite)
library(lubridate)
library(plotly)
library(quadprog)

source("graphs.r")

ui <- page_sidebar(
  title = "Otimizador de Portfólio Markowitz (Max Sharpe)",
  sidebar = sidebar(
    textInput(
      "tickers",
      "Tickers (separados por vírgula):",
      value = "XINA11.SA, BTC-BRL, VALE3.SA, ITUB4.SA, WEGE3.SA, EGIE3.SA, BRAX11.SA, GOLD11.SA"
    ),
    numericInput(
      "min_weight",
      "Peso Mínimo por Ativo (%):",
      value = 5,
      min = 0,
      max = 20
    ),
    numericInput(
      "ibov_alvo",
      "Ibovespa Alvo (para CAPM):",
      value = 196000,
      step = 1000
    ),
    numericInput(
      "dy_esperado",
      "Dividend Yield Mercado (%):",
      value = 5,
      min = 0,
      step = 0.5
    ),
    actionButton(
      "btn_run",
      "Otimizar Carteira",
      class = "btn-primary w-100 mt-3",
      style = "font-weight: bold;"
    )
  ),

  navset_card_tab(
    nav_panel(
      "Fronteira Eficiente",
      card(
        card_header("Risco vs Retorno (Gráfico Interativo)"),
        plotlyOutput("plot_fronteira", height = "600px")
      )
    ),
    nav_panel(
      "Backtest Histórico",
      card(
        card_header("Performance do Portfólio Tangente vs Ibovespa"),
        plotOutput("plot_backtest", height = "500px")
      )
    ),
    nav_panel(
      "Correlação",
      card(
        card_header("Matriz Linear de Pearson"),
        plotOutput("plot_heatmap", height = "600px")
      )
    ),
    nav_panel(
      "Relatório",
      card(
        card_header("Pesos Finais e Rendimentos"),
        tableOutput("tabela_pesos")
      )
    )
  )
)

server <- function(input, output, session) {
  # Processamento
  modelo <- eventReactive(input$btn_run, {
    showNotification(
      "Baixando dados e calculando fronteira... Pode levar de 15 a 30 segundos.",
      duration = 15,
      type = "message"
    )

    # Processar inputs
    tickers_input <- str_split(input$tickers, ",")[[1]] |> str_trim()
    tickers <- tickers_input[tickers_input != ""]

    peso_min <- input$min_weight / 100
    ibov_alvo_mediana <- input$ibov_alvo
    dy_esperado <- input$dy_esperado / 100

    end <- today()
    start <- end - years(2)

    # 1. Obter dados
    stocks_data <- tq_get(tickers, from = start, to = end, warnings = FALSE) |>
      distinct(symbol, date, .keep_all = TRUE) |>
      group_by(symbol) |>
      mutate(
        log_return = log(adjusted / lag(adjusted)),
        vol_252 = slide_dbl(log_return, sd, .before = 251, .complete = TRUE) *
          sqrt(252)
      ) |>
      drop_na(vol_252)

    ibov_data <- tq_get("^BVSP", from = start, to = end, warnings = FALSE) |>
      distinct(date, .keep_all = TRUE) |>
      mutate(log_return_m = log(adjusted / lag(adjusted))) |>
      drop_na(log_return_m) |>
      select(date, log_return_m, adjusted)

    ibov_atual <- ibov_data |> drop_na(adjusted) |> tail(1) |> pull(adjusted)

    # CDI API
    url_bcb <- "http://api.bcb.gov.br/dados/serie/bcdata.sgs.4391/dados?formato=json"
    cdi <- fromJSON(url_bcb) |>
      as_tibble() |>
      mutate(
        date = as.Date(data, format = "%d/%m/%Y"),
        cdi = as.numeric(valor) / 100
      ) |>
      select(date, cdi)

    r_f <- (1 + tail(cdi$cdi, 1))^12 - 1

    E_Rm <- log(ibov_alvo_mediana / ibov_atual) + dy_esperado

    # CAPM
    capm_results <- stocks_data |>
      select(symbol, date, log_return) |>
      group_by(symbol) |>
      drop_na(log_return) |>
      inner_join(ibov_data, by = "date") |>
      summarise(
        beta = cov(log_return, log_return_m) / var(log_return_m),
        .groups = "drop"
      ) |>
      mutate(
        expected_return = r_f + beta * (E_Rm - r_f)
      )

    # Covariância
    returns_wide <- stocks_data |>
      select(date, symbol, log_return) |>
      drop_na(log_return) |>
      pivot_wider(names_from = symbol, values_from = log_return) |>
      drop_na()

    cov_matrix <- cov(returns_wide[-1])
    cor_matrix <- cov2cor(cov_matrix)

    # Quadprog: Max Sharpe
    mu <- capm_results$expected_return
    n <- length(mu)

    excesso_retorno <- mu - r_f
    A_min_peso <- diag(n) - matrix(peso_min, n, n)
    Amat <- cbind(excesso_retorno, A_min_peso)
    bvec <- c(1, rep(0, n))

    otimizacao <- tryCatch(
      {
        solve.QP(
          Dmat = cov_matrix * 252,
          dvec = rep(0, n),
          Amat = Amat,
          bvec = bvec,
          meq = 1
        )
      },
      error = function(e) {
        showNotification(
          "Erro na Otimização Matemática! Verifique se os ativos permitem esse peso mínimo.",
          type = "error",
          duration = 8
        )
        return(NULL)
      }
    )

    req(otimizacao)

    y_otimo <- otimizacao$solution
    pesos_otimos <- y_otimo / sum(y_otimo)

    df_pesos <- tibble(
      symbol = row.names(cov_matrix),
      alocacao_otima_pct = pesos_otimos * 100
    )

    # Tabela consolidada
    stocks_stats <- capm_results |>
      left_join(
        stocks_data |>
          group_by(symbol) |>
          summarise(
            volatilidade_anual_ultima = last(vol_252),
            .groups = "drop"
          ),
        by = "symbol"
      ) |>
      mutate(
        sharpe_ratio = (expected_return - r_f) / volatilidade_anual_ultima
      ) |>
      left_join(df_pesos, by = "symbol") |>
      mutate(
        alocacao_otima_pct = replace_na(alocacao_otima_pct, 0),
        alocacao_ajustada_pct = ifelse(
          alocacao_otima_pct < 3,
          0,
          alocacao_otima_pct
        ),
        alocacao_ajustada_pct = (alocacao_ajustada_pct /
          sum(alocacao_ajustada_pct)) *
          100,
        across(where(is.numeric), \(x) round(x, 4))
      )

    # Backtest
    w_backtest <- stocks_stats |>
      select(symbol, alocacao_ajustada_pct) |>
      mutate(weight = alocacao_ajustada_pct / 100) |>
      select(symbol, weight)

    portfolio_returns <- stocks_data |>
      filter(symbol %in% w_backtest$symbol) |>
      drop_na(log_return) |>
      tq_portfolio(
        assets_col = symbol,
        returns_col = log_return,
        weights = w_backtest,
        col_rename = "retorno_carteira"
      )

    ibov_benchmark <- ibov_data |>
      mutate(retorno_ibov = log_return_m) |>
      select(date, retorno_ibov)

    comparativo_performance <- portfolio_returns |>
      left_join(ibov_benchmark, by = "date") |>
      drop_na(retorno_carteira, retorno_ibov) |>
      mutate(
        Carteira = exp(cumsum(retorno_carteira)) * 100,
        Ibovespa = exp(cumsum(retorno_ibov)) * 100
      ) |>
      pivot_longer(
        cols = c(Carteira, Ibovespa),
        names_to = "Estrategia",
        values_to = "Valor"
      )

    # Tabela final
    rendimento_historico <- stocks_data |>
      filter(symbol %in% w_backtest$symbol) |>
      group_by(symbol) |>
      summarise(
        Retorno_Real_Total = (exp(sum(log_return)) - 1) * 100,
        .groups = "drop"
      )

    perf_resumo <- comparativo_performance |>
      group_by(Estrategia) |>
      summarise(
        symbol = first(Estrategia),
        Retorno_Real_Total = (last(Valor) - 100),
        .groups = "drop"
      ) |>
      select(symbol, Retorno_Real_Total)

    tabela_final <- stocks_stats |>
      select(
        symbol,
        expected_return,
        beta,
        alocacao_otima_pct,
        alocacao_ajustada_pct
      ) |>
      left_join(rendimento_historico, by = "symbol") |>
      bind_rows(perf_resumo) |>
      mutate(
        expected_return_pct = expected_return * 100,
        across(where(is.numeric), \(x) round(x, 2))
      ) |>
      select(
        Ativo = symbol,
        `Peso Max Sharpe (%)` = alocacao_otima_pct,
        `Peso Ajustado (%)` = alocacao_ajustada_pct,
        `Retorno Esp. (%)` = expected_return_pct,
        `Retorno Real (%)` = Retorno_Real_Total,
        Beta = beta
      )

    names(mu) <- capm_results$symbol

    list(
      tabela_final = tabela_final,
      comparativo_performance = comparativo_performance,
      cor_matrix = cor_matrix,
      mu = mu,
      cov_matrix = cov_matrix,
      r_f = r_f
    )
  })

  # Outputs renderizados com os modulos do graphs.r
  output$plot_fronteira <- renderPlotly({
    req(modelo())
    res <- modelo()
    plot_fronteira_interativa(res$mu, res$cov_matrix, res$r_f)
  })

  output$plot_backtest <- renderPlot({
    req(modelo())
    res <- modelo()
    gera_graf_backtest(res$comparativo_performance)
  })

  output$plot_heatmap <- renderPlot({
    req(modelo())
    res <- modelo()
    gera_heatmap_correlacao(res$cor_matrix)
  })

  output$tabela_pesos <- renderTable(
    {
      req(modelo())
      modelo()$tabela_final
    },
    bordered = TRUE,
    hover = TRUE,
    striped = TRUE
  )
}

shinyApp(ui, server)
