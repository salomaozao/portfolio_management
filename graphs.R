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

  ggplot(fronteira_2, aes(x = risco, y = retorno)) +
    # Linha da Fronteira
    geom_path(aes(color = "Fronteira Eficiente"), size = 1) +
    # Ponto do Ativo A (mapeado manualmente para gerar legenda)
    geom_point(aes(x = v_a, y = r_a, color = "Ativo A"), size = 4) +
    # Ponto do Ativo B (mapeado manualmente para gerar legenda)
    geom_point(aes(x = v_b, y = r_b, color = "Ativo B"), size = 4) +
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
