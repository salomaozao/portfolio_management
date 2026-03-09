rm(list = ls())

source("graphs.r")


library(tidyquant)
library(tidyverse)
library(slider)
library(jsonlite)
library(lubridate)
#### ================== 1. Get & Proccess Data
# price, volatility, expected return, correlations (p) para preço e IBOV; Dados do CDI acumulado

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
  "TTWO",
  "SNDK"
)
end <- today()
start <- end - years(2)

stocks_data = tq_get(tickers, from = start, to = end) |>
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
  drop_na(log_return) |>
  # Alinhamos as datas das ações com as datas do Ibovespa
  inner_join(ibov_data, by = "date") |>
  group_by(symbol) |>
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

# 1. Definir os inputs baseados no que você já calculou
mu <- capm_results$expected_return # Vetor E[R] (CAPM)
stocks_processed = row.names(cov_matrix)
sigma_mat <- cov_matrix # Matriz de Covariância (Σ)
n <- length(mu) # Número de ativos (BABA e PETR4)

# 2. Configurar a Matriz de Restrições (Amat)
# Linha 1: soma dos pesos = 1
# Demais linhas: pesos individuais >= 0 (identidade)
Amat <- cbind(1, diag(n))

# 3. Configurar o vetor de restrições (bvec)
# O primeiro valor é 1 (soma), os demais são 0 (mínimo por ativo)
bvec <- c(1, rep(0, n))

# 4. Resolver para a Carteira de Variância Mínima Global (GMV)
# solve.QP(Dmat, dvec, Amat, bvec, meq)
# meq = 1 indica que a primeira restrição é uma IGUALDADE (=1)
otimizacao <- solve.QP(
  Dmat = cov_matrix * 252,
  dvec = rep(0, n),
  Amat = Amat,
  bvec = bvec,
  meq = 1
)

# 5. Extrair os pesos ideais
pesos_otimos <- otimizacao$solution
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
    across(where(is.numeric), \(x) round(x, 4))
  )

print("TABELA FINAL CONSOLIDADA:")
print(stocks_stats)

# Rodar os gráficos (Agora com a tabela correta)
gera_graf_fronteiras("TTWO", "VALE3.SA")
gera_fronteira_global(mu, cov_matrix, stocks_stats)

## ====================== BACKTESTING

# 1. Preparar pesos (Garantindo que a soma seja exatamente 1 para evitar erros no tq_portfolio)
w_backtest <- stocks_stats |>
  select(symbol, alocacao_otima_pct) |>
  mutate(weight = alocacao_otima_pct / 100) |>
  select(symbol, weight)

# 2. Calcular o retorno da carteira
portfolio_returns <- stocks_data |>
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

# 5. Plot Final
ggplot(comparativo_performance, aes(x = date, y = Valor, color = Estrategia)) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("Carteira" = "#2c3e50", "Ibovespa" = "#e74c3c")
  ) +
  labs(
    title = "Performance Histórica: Markowitz vs. Ibovespa",
    subtitle = "Simulação baseada na otimização de Variância Mínima",
    x = "Período",
    y = "Patrimônio Acumulado (R$)"
  ) +
  theme_minimal()
