rm(list = ls())


source("graphs.r")


library(tidyquant)
library(tidyverse)
library(slider)
library(jsonlite)
library(lubridate)

library(plotly)

#### ================== 1. Get & Proccess Data
# price, volatility, expected return, correlations (p) para pr+eço e IBOV; Dados do CDI acumulado

### Para as ações do portfólio: Retornos, volatilidade
tickers = c(
  "XINA11.SA",
  "BTC-BRL",
  "BSLV39.SA",
  "EMBR3.SA",
  "BRAX11.SA",
  "GOLD11.SA",
  "ITUB4.SA",
  "VALE3.SA",
  "AXIA3.SA",
  "CMIG4.SA",
  "CYRE4.SA",
  "RDOR3.SA",
  "WEGE3.SA",
  "KLBN11.SA",
  "EGIE3.SA",
  "TTWO",
  "SNDK"
)

end <- today()
start <- end - years(2)

stocks_data = tq_get(tickers, from = start, to = end) |>
  distinct(symbol, date, .keep_all = TRUE) |>
  group_by(symbol) |>
  mutate(
    log_return = log(adjusted / lag(adjusted)),
    # calcula a vol sobre o log_return
    vol_252 = slide_dbl(log_return, sd, .before = 251, .complete = TRUE) *
      sqrt(252)
  ) |>
  drop_na(vol_252) # remove os NAs do início da janela e do primeiro lag


# Checagem:
tickers_puxados <- unique(stocks_data$symbol)
tickers_faltantes <- setdiff(tickers, tickers_puxados)

if (length(tickers_faltantes) > 0) {
  warning(
    "OS SEGUINTES ATIVOS NÃO FORAM INCLUÍDOS NA ANÁLISE: ",
    paste(tickers_faltantes, collapse = ", ")
  )
} else {
  message(
    "SUCESSO: Todos os ",
    length(tickers),
    " ativos foram carregados e processados."
  )
}

### Preço do ibov; Retorno esperado do mercado (R_m)
ibov_data <- tq_get("^BVSP", from = start, to = end) |>
  distinct(date, .keep_all = TRUE) |>
  mutate(log_return_m = log(adjusted / lag(adjusted))) |>
  drop_na(log_return_m) |>
  select(date, log_return_m, adjusted)


ibov_atual <- ibov_data |>
  drop_na(adjusted) |>
  tail(1) |>
  pull(adjusted)


### taxa de juros livre de risco
# cdi <- tq_get("BCB/4391", get = "quandl") |>
#   select(date, cdi = "value") |>
#   mutate(cdi = cdi / 100, date = as.yearmon(date))

# URL direta da API do BCB para a série 4391
url_bcb <- "http://api.bcb.gov.br/dados/serie/bcdata.sgs.4391/dados?formato=json"

cdi <- fromJSON(url_bcb) |>
  as_tibble() |>
  mutate(
    # A API retorna 'data' no formato BR e 'valor' como texto
    date = as.Date(data, format = "%d/%m/%Y"),
    cdi = as.numeric(valor) / 100,
    date = zoo::as.yearmon(date)
  ) |>
  select(date, cdi)

### Para RETORNO ESPERADO:
# CAPM (Capital Asset Pricing Model):

# E[R_p] = R_f + β(E[R_m] - R_f), onde:

## R_f é a **taxa de risco livre** (CDI) (alpha)
# Transformando o último CDI mensal da sua base em uma taxa anualizada
# Fórmula: (1 + taxa_mensal)^12 - 1
r_f <- (1 + tail(cdi$cdi, 1))^12 - 1

## E[R_m] é o retorno esperado do mercado

ibov_alvo_mediana <- 196000 # Média obtida da XP, Safra, Eleven Fnancial, Itaú BBA (8/3/26)
dy_esperado <- 0.05
rf_anual <- 0.12 # Taxa livre de risco (Ajuste conforme a Selic atual)
E_Rm <- log(ibov_alvo_mediana / ibov_atual) + dy_esperado

## R_m - R_f é o prêmio de risco; prêmio por assumir mais risco

## aplicando o CAPM
capm_results <- stocks_data |>
  select(symbol, date, log_return) |>
  group_by(symbol) |>
  drop_na(log_return) |>
  # Alinhamos as datas das ações com as datas do Ibovespa
  inner_join(ibov_data, by = "date") |>
  summarise(
    # Beta = Covariância(Ativo, Mercado) / Variância(Mercado)
    beta = cov(log_return, log_return_m) / var(log_return_m), # Revisar
    .groups = "drop"
  ) |>
  mutate(
    # Aplicação direta da fórmula matemática do CAPM
    expected_return = r_f + beta * (E_Rm - r_f)
  )

## β é o coeficiente angular / de sensibilidade; o quanto o ativo "alavanca" o mercado
## E[R_p] é o retorno esperado do Portfólio

print("RETORNO ANUALIZADO ESPERADO")
print(capm_results)

## (!) O modelo assume NADA DE DIVIDENDOS (que se estivessem sendo somados no R_p estariam junto com o R_f, no alpha - como risco livre?); Isso pode ser implementado no projeto

### Para MATRIZ DE COVARIÂNCIA
# Matriz de Retornos
returns_wide <- stocks_data |>
  select(date, symbol, log_return) |>
  drop_na(log_return) |>
  pivot_wider(names_from = symbol, values_from = log_return) |>
  # Garantir que excluímos dias em que um dos ativos não teve negociação
  drop_na()

# Calcular a Covariância Amostral Diária
# Selecionamos todas as colunas numéricas (excluindo a coluna 1: 'date')
cov_matrix <- cov(returns_wide[-1])
cor_matrix = cov2cor(cov_matrix)

print("MATRIZ DE CORRELAÇÃO ENTRE ATIVOS")
print(cor_matrix)

# ================== . Definir pesos (w_i)
# Os pesos ($w$) são decididos através de três restrições fundamentais:Orçamento:
#  $\sum w_i = 1$ (todo o capital deve ser alocado).
#  Não-negatividade: $w_i \ge 0$ (estratégia long-only, sem venda a descoberto).
#  Retorno Alvo: $\sum w_i E[R_i] = \mu_{alvo}$.
library(quadprog)

# 1. Definir os inputs
mu <- capm_results$expected_return # Vetor E[R] (CAPM)
stocks_processed = row.names(cov_matrix)
sigma_mat <- cov_matrix # Matriz de Covariância (Σ)
n <- length(mu) # Número de ativos (BABA e PETR4)

# 2. Configurar a Matriz de Restrições (Amat)
# Para o Ponto Tangente (Máximo Sharpe Ratio), transformamos os pesos: y = w / (w^T * excesso_retorno)
excesso_retorno <- mu - r_f
A_min_peso <- diag(n) - matrix(0.05, n, n) # Restrição para w_i >= 0.05
Amat <- cbind(excesso_retorno, A_min_peso)

# 3. Configurar o vetor de restrições (bvec)
# A primeira é a igualdade ao retorno 1, o resto é a restrição >= 0
bvec <- c(1, rep(0, n))

# 4. Resolver para a Carteira Tangente (Máximo Sharpe)
# minimize y^T Sigma y  sujeito a y^T excesso_retorno = 1 e y_i - 0.05 * sum(y) >= 0
otimizacao <- solve.QP(
  Dmat = cov_matrix * 252,
  dvec = rep(0, n),
  Amat = Amat,
  bvec = bvec,
  meq = 1
)

# 5. Extrair os pesos ideais
y_otimo <- otimizacao$solution
pesos_otimos <- y_otimo / sum(y_otimo) # Normaliza para recuperar os pesos w que somam 1
df_pesos <- tibble(
  symbol = row.names(sigma_mat),
  alocacao_otima_pct = pesos_otimos * 100
)
print("Pesos Sugeridos pelo Modelo:")
print(df_pesos)

# ================== . Análises e Visualizações

# Criar a tabela base com volatilidade e Sharpe
stocks_stats <- capm_results |>
  left_join(
    stocks_data |>
      group_by(symbol) |>
      summarise(
        volatilidade_anual_ultima = last(vol_252),
        volatilidade_anual_média = mean(vol_252),
        .groups = "drop"
      ),
    by = "symbol"
  ) |>
  mutate(sharpe_ratio = (expected_return - r_f) / volatilidade_anual_ultima)

# Juntar os pesos GARANTINDO que a coluna tenha o nome correto
stocks_stats <- stocks_stats |>
  left_join(df_pesos, by = "symbol") |>
  mutate(
    # Garantir que se algum peso veio como NA (não processado), seja 0
    alocacao_otima_pct = replace_na(alocacao_otima_pct, 0),
    # Zera menores que 3% e recalcula para fechar 100%
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

print("TABELA FINAL CONSOLIDADA:")
print(stocks_stats)


mu <- capm_results$expected_return
names(mu) <- capm_results$symbol


## ====================== BACKTESTING

# 1. Preparar pesos (Garantindo que a soma seja exatamente 1 para evitar erros no tq_portfolio)
w_backtest <- stocks_stats |>
  select(symbol, alocacao_ajustada_pct) |>
  mutate(weight = alocacao_ajustada_pct / 100) |>
  select(symbol, weight)

# 2. Calcular o retorno da carteira
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
# 3. Benchmark e Acúmulo
comparativo_performance <- portfolio_returns |>
  left_join(ibov_benchmark, by = "date") |>
  drop_na(retorno_carteira, retorno_ibov) |>
  mutate(
    # Usando exp(cumsum) para retornos logarítmicos
    Carteira = exp(cumsum(retorno_carteira)) * 100,
    Ibovespa = exp(cumsum(retorno_ibov)) * 100
  ) |>
  pivot_longer(
    cols = c(Carteira, Ibovespa),
    names_to = "Estrategia",
    values_to = "Valor"
  )

# 1. Calcula o retorno acumulado real de cada ativo
rendimento_historico <- stocks_data |>
  filter(symbol %in% w_backtest$symbol) |>
  group_by(symbol) |>
  summarise(
    Retorno_Real_Total = (exp(sum(log_return)) - 1) * 100,
    .groups = "drop"
  )

# 2. Calcula o retorno acumulado da Carteira e do Ibovespa
perf_resumo <- comparativo_performance |>
  group_by(Estrategia) |>
  summarise(
    symbol = first(Estrategia),
    Retorno_Real_Total = (last(Valor) - 100), # Partindo de base 100
    .groups = "drop"
  ) |>
  select(symbol, Retorno_Real_Total)

# 3. Consolida tudo na Stocks Stats
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
    `Peso Máximo Sharpe (%)` = alocacao_otima_pct,
    `Peso Ajustado (%)` = alocacao_ajustada_pct,
    `Retorno Esperado (%)` = expected_return_pct,
    `Retorno Real (%)` = Retorno_Real_Total,
    Beta = beta
  )


# ================== . Visualizações
p_heatmap <- gera_heatmap_correlacao(cor_matrix)
p_heatmap

# gera_graf_fronteiras("TTWO", "VALE3.SA")
gera_fronteira_global(mu, cov_matrix, stocks_stats, n_sim = 50000)

p_backtest <- gera_graf_backtest(comparativo_performance)
p_backtest

p_interativo = plot_fronteira_interativa(mu, cov_matrix, r_f)
print(p_interativo)

print("RELATÓRIO DE RENDIMENTOS:")
print(tabela_final)
  