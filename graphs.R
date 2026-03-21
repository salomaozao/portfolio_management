gera_graf_fronteiras = function(a, b) {
  # Seleção dos ativos
  sim_pesos <- seq(0, 1, length.out = 100)
  ativo_a <- a
  ativo_b <- b

  v_a <- stocks_stats$volatilidade_anual_ultima[stocks_stats$symbol == ativo_a]
  v_b <- stocks_stats$volatilidade_anual_ultima[stocks_stats$symbol == ativo_b]
  r_a <- stocks_stats$expected_return[stocks_stats$symbol == ativo_a]
  r_b <- stocks_stats$expected_return[stocks_stats$symbol == ativo_b]
  cor_ab <- cor_matrix[ativo_a, ativo_b]

  fronteira_2 <- tibble(
    w_a = sim_pesos,
    w_b = 1 - sim_pesos,
    retorno = w_a * r_a + w_b * r_b,
    # Fórmula da variância de 2 ativos
    risco = sqrt(
      (w_a^2 * v_a^2) + (w_b^2 * v_b^2) + (2 * w_a * w_b * v_a * v_b * cor_ab)
    )
  )

  ggplot(fronteira_2, aes(x = retorno, y = risco)) +
    # Linha da Fronteira
    geom_path(aes(color = "Fronteira Eficiente"), size = 1) +
    # Ponto do Ativo A (mapeado manualmente para gerar legenda)
    geom_point(aes(x = r_a, y = v_a, color = "Ativo A"), size = 4) +
    # Ponto do Ativo B (mapeado manualmente para gerar legenda)
    geom_point(aes(x = r_b, y = v_b, color = "Ativo B"), size = 4) +
    # Customização das cores e nomes na legenda
    scale_color_manual(
      name = "Legenda",
      values = c(
        "Ativo A" = "red",
        "Ativo B" = "green",
        "Fronteira Eficiente" = "darkblue"
      ),
      labels = c(ativo_a, ativo_b, "Fronteira de Alocação")
    ) +
    labs(
      title = paste("Fronteira de Risco-Retorno:", ativo_a, "vs", ativo_b),
      subtitle = paste("Correlação entre ativos:", round(cor_ab, 4)),
      x = "Volatilidade (Risco)",
      y = "Retorno Esperado (CAPM)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


gera_fronteira_global = function(mu, sigma_mat, stocks_stats, n_sim = 5000) {
  n_ativos <- length(mu)

  # 1. Simulação de Portfólios Aleatórios para preencher o gráfico
  sim_resultados <- matrix(NA, nrow = n_sim, ncol = 2)
  colnames(sim_resultados) <- c("Risco", "Retorno")

  for (i in 1:n_sim) {
    w_random <- runif(n_ativos)
    w_random <- w_random / sum(w_random)

    sim_resultados[i, "Retorno"] <- sum(w_random * mu)
    sim_resultados[i, "Risco"] <- sqrt(
      t(w_random) %*% (sigma_mat * 252) %*% w_random
    )
  }

  df_sim <- as_tibble(sim_resultados)

  # 2. Ponto da Carteira Ótima (já calculada pelo seu solve.QP)
  # Assumindo que 'otimizacao' e 'pesos_otimos' já existem no seu environment
  risco_otimo <- sqrt(t(pesos_otimos) %*% (sigma_mat * 252) %*% pesos_otimos)
  retorno_otimo <- sum(pesos_otimos * mu)

  # 3. Plotagem
  ggplot() +
    # Nuvem de portfólios possíveis
    geom_point(
      data = df_sim,
      aes(x = Risco, y = Retorno),
      alpha = 0.2,
      color = "grey"
    ) +
    # Ativos individuais
    geom_point(
      data = stocks_stats,
      aes(x = volatilidade_anual_ultima, y = expected_return, color = symbol),
      size = 3
    ) +
    geom_text(
      data = stocks_stats,
      aes(x = volatilidade_anual_ultima, y = expected_return, label = symbol),
      vjust = -1,
      size = 3
    ) +
    # Carteira de Variância Mínima (O seu resultado do solve.QP)
    geom_point(
      aes(x = risco_otimo, y = retorno_otimo),
      color = "black",
      shape = 18,
      size = 6
    ) +
    annotate(
      "text",
      x = risco_otimo,
      y = retorno_otimo,
      label = "CARTEIRA ÓTIMA (GMV)",
      vjust = 2,
      fontface = "bold"
    ) +
    labs(
      title = "Fronteira Eficiente Global de Markowitz",
      subtitle = "Comparação entre ativos individuais e a alocação otimizada",
      x = "Risco (Volatilidade Anualizada)",
      y = "Retorno Esperado Anual (CAPM)",
      color = "Ativos"
    ) +
    theme_minimal()
}


plot_fronteira_interativa <- function(mu, cov_matrix, r_f, n_sim = 5000) {
  # Garante que mu seja um vetor nomeado (essencial para o hover)
  if (is.null(names(mu))) {
    names(mu) <- paste0("Ativo_", 1:length(mu))
  }

  n_ativos <- length(mu)
  nomes_ativos <- names(mu)

  # 1. Simulação
  sim_mat <- matrix(NA, nrow = n_sim, ncol = n_ativos + 2)

  for (i in 1:n_sim) {
    w <- runif(n_ativos)
    w <- w / sum(w)

    # Zera os pesos menores que 3%
    w[w < 0.03] <- 0
    w <- w / sum(w)

    ret_sim <- sum(w * mu)
    risco_sim <- sqrt(t(w) %*% (cov_matrix * 252) %*% w)
    sim_mat[i, ] <- c(w, ret_sim, risco_sim)
  }

  # 2. Criação do DF com nomes explícitos e limpos
  df_sim <- as.data.frame(sim_mat)
  # Forçamos nomes válidos para evitar o erro de mutate()
  colnames(df_sim) <- make.names(c(nomes_ativos, "Retorno", "Risco"))

  # Atualizamos a lista de nomes para bater com os novos nomes da coluna
  nomes_colunas_ativos <- colnames(df_sim)[1:n_ativos]

  # 3. Processamento do Hover
  df_sim <- df_sim %>%
    mutate(
      sharpe = (Retorno - r_f) / Risco,
      # Usamos rowwise() para facilitar a iteração por linha
      texto_hover = pmap_chr(
        select(., all_of(nomes_colunas_ativos)),
        function(...) {
          pesos <- c(...)
          txt <- paste0(
            nomes_ativos,
            ": ",
            round(pesos * 100, 1),
            "%",
            collapse = "<br>"
          )
          return(txt)
        }
      )
    )

  # 4. Plotly
  plot_ly(
    data = df_sim,
    x = ~Risco,
    y = ~Retorno,
    type = 'scatter',
    mode = 'markers',
    marker = list(
      color = ~sharpe,
      colorscale = 'Viridis',
      showscale = TRUE,
      size = 5,
      opacity = 0.6
    ),
    text = ~texto_hover,
    hoverinfo = 'text+x+y'
  ) %>%
    plotly::layout(
      title = "Fronteira Eficiente: Passe o mouse para ver os pesos",
      xaxis = list(title = "Risco (Volatilidade Anual)", tickformat = ".1%"),
      yaxis = list(title = "Retorno Esperado (Anual)", tickformat = ".1%"),
      hovermode = "closest"
    )
}

# ================== Gráficos de Backtest e Heatmap ==================

gera_graf_backtest <- function(comparativo_performance) {
  ggplot(comparativo_performance, aes(x = date, y = Valor, color = Estrategia)) +
    geom_line(size = 1) +
    scale_color_manual(
      values = c("Carteira" = "#2c3e50", "Ibovespa" = "#e74c3c")
    ) +
    labs(
      title = "Performance Histórica: Markowitz vs. Ibovespa",
      subtitle = "Simulação baseada na otimização de Máximo Sharpe Ratio",
      x = "Período",
      y = "Patrimônio Acumulado (R$)"
    ) +
    theme_minimal()
}

gera_heatmap_correlacao <- function(cor_matrix) {
  cor_long <- as.data.frame(cor_matrix) |>
    rownames_to_column(var = "Ativo1") |>
    pivot_longer(cols = -Ativo1, names_to = "Ativo2", values_to = "Correlacao")

  ggplot(cor_long, aes(x = Ativo1, y = Ativo2, fill = Correlacao)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlacao, 2)), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "white",
      high = "#4575b4",
      midpoint = 0,
      limit = c(-1, 1),
      name = "Correlação\n(Pearson)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(
      title = "Matriz de Dependência Linear entre Ativos",
      subtitle = "Calculada sobre os log-retornos diários",
      x = "",
      y = ""
    )
}
