library(httr)
library(jsonlite)

#' Função para buscar dados de investimentos via BRAPI
#' @param ticker Ex: "WEGE3", "BTC", "BRAX11"
#' @param tipo Ex: "preco", "var24h", "var7d", "var1mes", "pl", "lpa", "nome"
#' @param token Seu token da API BRAPI
brapi_dados <- function(
  ticker,
  tipo = "preco",
  token = "tLki4EByiYvzUj9VTpT1Xe"
) {
  if (missing(ticker) || ticker == "") {
    return("Ticker vazio")
  }
  if (missing(tipo) || tipo == "") {
    return("Tipo não informado")
  }

  ticker <- toupper(trimws(ticker))
  tipo_desejado <- tolower(tipo)

  # 1. Identifica se é Crypto ou Ação
  is_crypto <- ticker %in% c("BTC", "ETH", "SOL", "USDT", "ADA", "DOT")

  if (is_crypto) {
    url <- paste0(
      "https://brapi.dev/api/v2/crypto?coin=",
      ticker,
      "&currency=BRL&token=",
      token
    )
  } else {
    url <- paste0(
      "https://brapi.dev/api/quote/",
      ticker,
      "?token=",
      token,
      "&range=1mo&interval=1d"
    )
  }

  # 2. Requisição
  res <- GET(url)

  if (status_code(res) != 200) {
    return("Erro na API ou Ticker inválido")
  }

  json <- fromJSON(content(res, as = "text", encoding = "UTF-8"))

  # 3. Extração dos dados
  # CORREÇÃO AQUI: Pegamos a primeira linha do data.frame, não a primeira coluna
  dados <- if (is_crypto) json$coins[1, ] else json$results[1, ]

  # 4. Lógica de retorno
  resultado <- switch(
    tipo_desejado,
    "preco" = dados$regularMarketPrice,
    "var24h" = (if (is_crypto) {
      dados$regularMarketChangePercent
    } else {
      dados$regularMarketChangePercent
    }) /
      100,
    "var7d" = {
      # CORREÇÃO AQUI: O histórico vira uma lista contendo um data.frame, precisamos do [[1]]
      hist <- if (!is_crypto) dados$historicalDataPrice[[1]] else NULL
      if (is.null(hist) || nrow(hist) < 5) {
        return("N/A")
      }
      (dados$regularMarketPrice / hist$close[1]) - 1
    },
    "var1mes" = {
      hist <- if (!is_crypto) dados$historicalDataPrice[[1]] else NULL
      if (is.null(hist) || nrow(hist) < 15) {
        return("N/A")
      }
      (dados$regularMarketPrice / hist$close[1]) - 1
    },
    "pl" = ifelse(
      !is.null(dados$priceEarnings) && !is.na(dados$priceEarnings),
      dados$priceEarnings,
      "N/A"
    ),
    "lpa" = ifelse(
      !is.null(dados$earningsPerShare) && !is.na(dados$earningsPerShare),
      dados$earningsPerShare,
      "N/A"
    ),
    "nome" = if (is_crypto) dados$coinName else dados$longName,
    "Tipo inválido"
  )

  return(resultado)
}

library(httr)
library(jsonlite)

#' Função para retornar matriz de High/Low de uma ação
#' @param ticker Ex: "PETR4", "VALE3"
#' @param periodo Ex: "5d", "1mo", "3mo", "6mo", "1y", "2y", "5y", "max"
#' @param token Seu token da API BRAPI
brapi_high_low <- function(ticker, periodo = "1mo", token = "tLki4EByiYvzUj9VTpT1Xe") {
  
  if (missing(ticker) || ticker == "") stop("Ticker não pode ser vazio.")
  
  ticker <- toupper(trimws(ticker))
  
  # Montando a URL com o range dinâmico
  url <- paste0("https://brapi.dev/api/quote/", ticker, 
                "?token=", token, 
                "&range=", periodo, 
                "&interval=1d")
  
  # Requisição
  res <- GET(url)
  
  if (status_code(res) != 200) {
    stop(paste("Erro na API:", status_code(res), "- Verifique o ticker ou o token."))
  }
  
  json <- fromJSON(content(res, as = "text", encoding = "UTF-8"))
  
  # Extrair o histórico de preços (desempacotando o data.frame interno)
  hist_df <- json$results$historicalDataPrice[[1]]
  
  if (is.null(hist_df) || nrow(hist_df) == 0) {
    stop("Nenhum dado histórico encontrado para este período.")
  }
  
  # A BRAPI retorna a data em formato Timestamp Unix. Vamos converter para formato de Data legível.
  datas_legiveis <- as.Date(as.POSIXct(hist_df$date, origin = "1970-01-01"))
  
  # Criando a matriz apenas com as colunas High e Low
  matriz_hl <- cbind(High = hist_df$high, Low = hist_df$low)
  
  # Nomeando as linhas da matriz com as datas para facilitar a identificação
  rownames(matriz_hl) <- as.character(datas_legiveis)
  
  return(matriz_hl)
}