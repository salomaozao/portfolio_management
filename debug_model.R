library(tidyquant)
library(tidyverse)
library(slider)
library(jsonlite)
library(lubridate)
library(quadprog)

tickers_input <- "XINA11.SA, BTC-BRL, VALE3.SA, ITUB4.SA, WEGE3.SA, EGIE3.SA, BRAX11.SA, GOLD11.SA, XOM, LMT"
tickers <- str_split(tickers_input, ",")[[1]] |> str_trim()
peso_min <- 5 / 100
ibov_alvo_mediana <- 196000
dy_esperado <- 5 / 100

end <- today()
start <- end - years(2)

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

returns_wide <- stocks_data |>
  select(date, symbol, log_return) |>
  drop_na(log_return) |>
  pivot_wider(names_from = symbol, values_from = log_return) |>
  drop_na()

cat("Dimensions of returns_wide: ", dim(returns_wide), "\n")

cov_matrix <- cov(returns_wide[-1])
cor_matrix <- cov2cor(cov_matrix)

cat("Eigenvalues of cov_matrix: ", min(eigen(cov_matrix)$values), "\n")

mu <- capm_results$expected_return
n <- length(mu)

excesso_retorno <- mu - r_f
A_min_peso <- diag(n) - matrix(peso_min, n, n)
Amat <- cbind(excesso_retorno, A_min_peso)
bvec <- c(1, rep(0, n))

cat("Solving QP...\n")
tryCatch({
  otimizacao <- solve.QP(
    Dmat = cov_matrix * 252,
    dvec = rep(0, n),
    Amat = Amat,
    bvec = bvec,
    meq = 1
  )
  print("QP Solved!")
  print(otimizacao$solution)
}, error = function(e) {
  print("Error in solve.QP:")
  print(e)
})
