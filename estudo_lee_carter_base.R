# ESTUDO LEE-CARTER — IMPLEMENTACAO DIDATICA, R BASE
# Referencia principal: Lee & Carter (1992), JASA 87(419), 659-671.
# Finalidade: reconstruir a aproximacao SVD e projetar o indice temporal.
# Dados: SIMULADOS; nenhuma estimativa corresponde aos artigos.
# Requisitos: R >= 4.1; sem dependencias adicionais.

set.seed(20261012)
idades <- seq(0, 90, by = 5)
anos <- 1980:2024
nx <- length(idades)
nt <- length(anos)

# E[x,t] = exposicao, medida em pessoas-ano.
E <- outer(
  95000 * exp(-idades / 58),
  1 + 0.004 * (anos - min(anos)),
  FUN = "*"
)

# Parametros generativos de uma populacao hipotetica.
alpha_real <- -9.0 + 0.075 * idades
beta_real <- (1.2 - 0.005 * idades)
beta_real <- beta_real / sum(beta_real)
kappa_real <- -0.55 * (anos - mean(anos))

log_m <- outer(alpha_real, rep(1, nt)) + outer(beta_real, kappa_real)
mu <- E * exp(log_m)
D <- matrix(rpois(length(mu), lambda = as.vector(mu)), nrow = nx)
rownames(D) <- rownames(E) <- as.character(idades)
colnames(D) <- colnames(E) <- as.character(anos)

# Ajuste log-linear via SVD (1 componente).
ajustar_lee_carter <- function(obitos, exposicao, idades, anos,
                                recalibrar_obitos = FALSE) {
  if (!is.matrix(obitos) || !is.matrix(exposicao) ||
      !identical(dim(obitos), dim(exposicao))) {
    stop("'obitos' e 'exposicao' devem ser matrizes da mesma dimensao.")
  }
  if (length(idades) != nrow(obitos) || length(anos) != ncol(obitos)) {
    stop("Numero de idades ou anos incompativel com as matrizes.")
  }
  if (anyNA(obitos) || anyNA(exposicao) ||
      any(obitos < 0) || any(exposicao <= 0)) {
    stop("Dados ausentes, obitos negativos ou exposicoes nao positivas.")
  }
  if (ncol(obitos) < 3L || nrow(obitos) < 2L) {
    stop("Necessarias pelo menos 2 idades e 3 anos.")
  }

  # Correcao de 0.5 somente para log de taxas com zero obito.
  # Este passo nao e igual ao ajuste Poisson por maxima verossimilhanca.
  log_taxas <- log((obitos + 0.5) / exposicao)
  alpha <- rowMeans(log_taxas)
  centralizada <- sweep(log_taxas, 1, alpha, FUN = "-")
  decomp <- svd(centralizada, nu = 1, nv = 1)
  beta_bruto <- decomp$u[, 1]
  s <- sum(beta_bruto)
  if (abs(s) < 1e-10) {
    stop("Primeiro vetor singular tem soma proxima de zero.")
  }
  beta <- beta_bruto / s
  kappa <- as.numeric(decomp$d[1] * decomp$v[, 1] * s)

  # Impor soma(kappa)=0, preservando o produto alpha + beta*kappa.
  centro <- mean(kappa)
  alpha <- alpha + beta * centro
  kappa <- kappa - centro

  # Etapa original opcional: ajustar kappa_t para igualar totais de obitos.
  # Aqui requer beta_x estritamente positivo para assegurar monotonicidade.
  if (isTRUE(recalibrar_obitos)) {
    if (any(beta <= 0)) {
      stop("Recalibracao simples exige todos os beta_x > 0.")
    }
    for (j in seq_along(anos)) {
      alvo <- sum(obitos[, j])
      f <- function(k) {
        sum(exposicao[, j] * exp(alpha + beta * k)) - alvo
      }
      kappa[j] <- uniroot(f, interval = kappa[j] + c(-300, 300))$root
    }
    centro <- mean(kappa)
    alpha <- alpha + beta * centro
    kappa <- kappa - centro
  }

  taxas_ajustadas <- exp(outer(alpha, rep(1, length(anos))) +
                         outer(beta, kappa))
  dimnames(taxas_ajustadas) <- list(as.character(idades),
                                    as.character(anos))
  list(alpha = alpha, beta = beta, kappa = kappa,
       taxas_ajustadas = taxas_ajustadas,
       fracao_singular = decomp$d[1]^2 / sum(decomp$d^2),
       idades = idades, anos = anos)
}

# Passeio aleatorio com deriva, horizonte h.
# Limites APROXIMADOS: variacao futura de k e incerteza na deriva.
# NAO incluem incerteza na estimacao de alpha/beta ou observacao Poisson.
prever_lee_carter <- function(ajuste, h = 10L, nivel = 0.95) {
  if (length(h) != 1L || h < 1 || h != as.integer(h)) {
    stop("h precisa ser um inteiro positivo.")
  }
  if (nivel <= 0 || nivel >= 1) stop("nivel deve estar entre 0 e 1.")
  delta <- diff(ajuste$kappa)
  deriva <- mean(delta)
  sigma <- sd(delta)
  passos <- seq_len(h)
  k_futuro <- tail(ajuste$kappa, 1) + passos * deriva
  # Variancia de inovacao h*sigma^2 e de deriva h^2*sigma^2/n_delta.
  se_k <- sigma * sqrt(passos + passos^2 / length(delta))
  z <- qnorm((1 + nivel) / 2)
  log_m <- outer(ajuste$alpha, rep(1, h)) +
    outer(ajuste$beta, k_futuro)
  amplitude <- abs(outer(ajuste$beta, z * se_k))
  anos_futuros <- max(ajuste$anos) + passos
  nomes <- list(as.character(ajuste$idades), as.character(anos_futuros))
  taxas <- exp(log_m)
  inferior <- exp(log_m - amplitude)
  superior <- exp(log_m + amplitude)
  dimnames(taxas) <- dimnames(inferior) <- dimnames(superior) <- nomes
  list(anos = anos_futuros, kappa = k_futuro, se_k = se_k,
       deriva = deriva, sigma_inovacao = sigma,
       taxas = taxas, inferior = inferior, superior = superior)
}

# VALIDACAO TEMPORAL: treino ate 2017; teste de 2018 ate 2024.
idx_treino <- which(anos <= 2017)
idx_teste <- which(anos > 2017)
modelo <- ajustar_lee_carter(
  D[, idx_treino, drop = FALSE], E[, idx_treino, drop = FALSE],
  idades, anos[idx_treino], recalibrar_obitos = FALSE
)
previsao <- prever_lee_carter(modelo, h = length(idx_teste))
stopifnot(identical(as.integer(colnames(previsao$taxas)), anos[idx_teste]))

obitos_teste <- D[, idx_teste, drop = FALSE]
expo_teste <- E[, idx_teste, drop = FALSE]
taxa_teste <- (obitos_teste + 0.5) / expo_teste
rmse_log <- sqrt(mean((log(taxa_teste) - log(previsao$taxas))^2))
lambda_prevista <- expo_teste * previsao$taxas
positivos <- obitos_teste > 0
poisson_deviance <- 2 * (sum(lambda_prevista - obitos_teste) +
    sum(obitos_teste[positivos] *
          log(obitos_teste[positivos] / lambda_prevista[positivos])))

cat("Soma(beta) =", sum(modelo$beta), "\n")
cat("Soma(kappa) =", sum(modelo$kappa), "\n")
cat("Fracao singular (SVD) =", modelo$fracao_singular, "\n")
cat("Deriva estimada =", previsao$deriva, "\n")
cat("RMSE em log-taxa no teste =", rmse_log, "\n")
cat("Deviance Poisson no teste =", poisson_deviance, "\n")

# Ajustar com todo o periodo, depois projetar os 10 anos seguintes.
modelo_completo <- ajustar_lee_carter(D, E, idades, anos)
projecao_2034 <- prever_lee_carter(modelo_completo, h = 10)

# Figuras. Os graficos utilizam idades quinquenais no exemplo simulado.
png("01_parametros_lee_carter.png", width = 1100, height = 400)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(idades, modelo_completo$alpha, type = "b", pch = 16,
     xlab = "Idade", ylab = "alpha(x)", main = "Padrao por idade")
plot(idades, modelo_completo$beta, type = "b", pch = 16,
     xlab = "Idade", ylab = "beta(x)", main = "Sensibilidade temporal")
plot(anos, modelo_completo$kappa, type = "l", lwd = 2,
     xlab = "Ano", ylab = "kappa(t)", main = "Indice de mortalidade")
dev.off()

idade_alvo <- "60"
png("02_previsao_idade_60.png", width = 850, height = 540)
serie <- modelo_completo$taxas_ajustadas[idade_alvo, ]
limites_y <- range(c(serie, projecao_2034$inferior[idade_alvo, ],
                     projecao_2034$superior[idade_alvo, ]))
plot(anos, serie, type = "l", lwd = 2, ylim = limites_y,
     xlim = range(c(anos, projecao_2034$anos)),
     xlab = "Ano", ylab = "Taxa central (por pessoa-ano)",
     main = "Lee-Carter: previsao para idade 60")
lines(projecao_2034$anos, projecao_2034$taxas[idade_alvo, ], lwd = 2,
      lty = 2)
lines(projecao_2034$anos, projecao_2034$inferior[idade_alvo, ],
      lty = 3)
lines(projecao_2034$anos, projecao_2034$superior[idade_alvo, ],
      lty = 3)
legend("topright", c("Ajuste", "Previsao", "Limites aproximados 95%"),
       lty = c(1, 2, 3), bty = "n")
dev.off()

# Objetos 'D', 'E', 'modelo_completo' e 'projecao_2034' podem
# ser reutilizados no segundo script, caso ele seja executado depois.
