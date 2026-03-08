library(tidyquant)
library(tidyverse)
library(slider)
library(jsonlite)

#### ================== 1. Get & Proccess Data
# price, volatility, expected return, correlations (p) para preço e IBOV; Dados do CDI acumulado

### Para as ações do portfólio: Retornos, volatilidade
tickers = c("BABA", "PETR4.SA")

start = "2018-01-01"
end = "2022-08-15"

stocks_data = tq_get(tickers, from = start, to = end) |>
  group_by(symbol) |>
  mutate(
    log_return = log(adjusted / lag(adjusted)),
    # calcula a vol sobre o log_return
    vol_252 = slide_dbl(log_return, sd, .before = 251, .complete = TRUE) *
      sqrt(252)
  ) |>
  drop_na(vol_252) # remove os NAs do início da janela e do primeiro lag


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
  inner_join(ibov_returns, by = "date") |>
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
  Dmat = sigma_mat,
  dvec = rep(0, n),
  Amat = Amat,
  bvec = bvec,
  meq = 1
)

# 5. Extrair os pesos ideais
pesos_otimos <- otimizacao$solution
names(pesos_otimos) <- names(mu)

print("Pesos Sugeridos pelo Modelo:")
print(round(pesos_otimos, 4))
