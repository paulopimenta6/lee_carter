# EXTENSOES LEE-CARTER — EXEMPLOS ILUSTRATIVOS EM R
# Dependencia sugerida para aplicacoes profissionais: StMoMo (CRAN).
# As simulacoes abaixo NAO reproduzem os algoritmos MCMC dos artigos.

# (A) Czado et al. (2005): verossimilhanca Poisson log-bilinear.
logverossimilhanca_poisson <- function(D, E, alpha, beta, kappa) {
  stopifnot(is.matrix(D), is.matrix(E), identical(dim(D), dim(E)),
            nrow(D) == length(alpha), nrow(D) == length(beta),
            ncol(D) == length(kappa), all(E > 0), all(D >= 0))
  eta <- outer(alpha, rep(1, ncol(D))) + outer(beta, kappa)
  sum(dpois(as.vector(D), lambda = as.vector(E * exp(eta)), log = TRUE))
}

# Executar este bloco depois de 'source("estudo_lee_carter_base.R")'.
if (exists("modelo_completo") && exists("D") && exists("E")) {
  ll <- logverossimilhanca_poisson(
    D, E, modelo_completo$alpha, modelo_completo$beta,
    modelo_completo$kappa
  )
  cat("Log-verossimilhanca Poisson do ajuste SVD:", ll, "\n")
}

# Para ajuste Poisson profissional por maxima verossimilhanca:
# install.packages("StMoMo")
# library(StMoMo)
# data("EWMaleData", package = "StMoMo")
# ajuste_pois <- StMoMo::fit(StMoMo::lc(), data = EWMaleData,
#                            ages.fit = 55:89, verbose = FALSE)
# proj_pois <- forecast(ajuste_pois, h = 10)
# plot(ajuste_pois)
# plot(proj_pois)
# NOTA: essa estimacao Poisson por MV NAO e uma analise bayesiana MCMC.

# (B) Pedroza (2006): geracao didatica da equacao de estado.
set.seed(2206)
Tn <- 35L
p <- 6L
alpha <- seq(-8, -4, length.out = p)
beta <- rep(1 / p, p)
theta <- -0.15
sigma_omega <- 0.10
sigma_epsilon <- 0.08
k <- numeric(Tn)
for (t in 2:Tn) {
  k[t] <- k[t - 1] + theta + rnorm(1, 0, sigma_omega)
}
Y <- outer(alpha, rep(1, Tn)) + outer(beta, k) +
  matrix(rnorm(p * Tn, 0, sigma_epsilon), nrow = p)
cat("Dimensoes Y (idades x anos) no espaco de estados:", dim(Y), "\n")
# Aqui apenas SIMULAMOS y_t = alpha + beta*k_t + epsilon_t.
# Inferencia bayesiana exigiria prioris, Kalman + FFBS e MCMC.

# (C) Liu & Li (2015): J0, J1 e J2 diferenciam a idade do choque.
set.seed(2015)
idades <- seq(0, 90, by = 5)
p <- length(idades)
alpha <- -9 + 0.075 * idades
beta <- rep(1 / p, p)
beta_choque <- 0.5 + 2 * exp(-((idades - 25) / 15)^2)
beta_choque <- beta_choque / sum(beta_choque)
k_atual <- -5
choque <- 10  # uma severidade positiva, apenas ilustrativa
rho <- 0.65  # correlacao entre faixas etarias: 0 <= rho < 1
w_comum <- rnorm(1)
w_local <- rnorm(p)
w <- sqrt(rho) * w_comum + sqrt(1 - rho) * w_local
J2 <- beta_choque * (choque + 1.5 * w)
log_base <- alpha + beta * k_atual
m_base <- exp(log_base)
m_J0 <- exp(log_base + beta * choque)
m_J1 <- exp(log_base + beta_choque * choque)
m_J2 <- exp(log_base + J2)

png("03_choques_mortalidade.png", width = 850, height = 530)
plot(idades, m_J0 / m_base, type = "l", lwd = 2,
     ylim = range(c(m_J0, m_J1, m_J2) / m_base),
     xlab = "Idade", ylab = "Razao de taxas em relacao ao caso sem choque",
     main = "Choques de mortalidade: padroes etarios")
lines(idades, m_J1 / m_base, lty = 2, lwd = 2)
lines(idades, m_J2 / m_base, lty = 3, lwd = 2)
legend("topright", c("J0: mesmo beta", "J1: beta do choque",
                     "J2: padrao aleatorio"), lty = 1:3, bty = "n")
dev.off()

# (D) Wisniowski et al. (2015): esqueleto de coortes.
# APENAS DIDATICO: taxas constantes no ano e migracao simplificada;
# NAO reproduz o modelo conjunto nem o ajuste bayesiano do artigo.
projetar_coorte <- function(pop, morte, emigracao, imigrantes,
                            nascimentos_sobreviventes) {
  n <- length(pop)
  stopifnot(n >= 2L, length(morte) == n,
            length(emigracao) == n, length(imigrantes) == n,
            all(pop >= 0), all(morte >= 0), all(emigracao >= 0),
            all(imigrantes >= 0), nascimentos_sobreviventes >= 0)
  sobrev <- pop * exp(-(morte + emigracao))
  proxima <- numeric(n)
  proxima[1] <- nascimentos_sobreviventes + imigrantes[1]
  proxima[2:n] <- sobrev[1:(n-1)] + imigrantes[2:n]
  # A ultima faixa etaria e aberta (ex.: 90+).
  proxima[n] <- proxima[n] + sobrev[n]
  proxima
}

pop <- c(8000, 7500, 7200, 6400, 5200)
morte <- c(0.002, 0.0005, 0.003, 0.025, 0.08)
emig <- c(0.002, 0.006, 0.01, 0.003, 0.001)
immig <- c(50, 100, 80, 25, 5)
print(projetar_coorte(pop, morte, emig, immig, 600))
