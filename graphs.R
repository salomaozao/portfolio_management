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
