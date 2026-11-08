# Guia de estudo de Lee Carter

##### [**Undermind**](https://undermind.ai)

---


## Table of Contents

- [Guia de estudo — Lee–Carter e extensões](#guia-de-estudo-leecarter-e-extensões)
  - [Questão de pesquisa](#questão-de-pesquisa)
  - [Percurso dos cinco artigos](#percurso-dos-cinco-artigos)
  - [Prática em R](#prática-em-r)
  - [Pontos de rigor](#pontos-de-rigor)
  - [References](#references)

# Guia de estudo — Lee–Carter e extensões

## Questão de pesquisa

Como decompor log-taxas de mortalidade por idade e ano, prever a evolução temporal com incerteza e incorporar contagens, efeitos de coorte, migração e choques?

## Percurso dos cinco artigos

1.  Modelo original \[Lee92\]: $`\log m_{x,t}=\alpha_x+\beta_x\kappa_t+\varepsilon_{x,t}`$, restrições $`\sum_x\beta_x=1`$, $`\sum_t\kappa_t=0`$. Estimação SVD de primeira componente e projeção de $`\kappa_t`$ via passeio aleatório com deriva. Verificar a recalibração do índice para os totais de óbitos.
2.  Alternativa Poisson bayesiana \[Cza05\]: $`D_{x,t}\sim\text{Poisson}(E_{x,t}\exp(\alpha_x+\beta_x\kappa_t))`$. MCMC conjunto com priori autoregressiva de tendência para o índice temporal, normal para beta, gama para exponencial de alfa e hiperprioris.
3.  Espaço de estados bayesiano \[Ped06\]: $`y_t=\alpha+\beta\kappa_t+\varepsilon_t`$; $`\kappa_t=\kappa_{t-1}+\theta+\omega_t`$. Filtro de Kalman, simulação posterior e dados faltantes.
4.  Projeção populacional \[Wis15\]: Poisson para mortalidade, fecundidade, emigrantes e imigrantes; imigrantes modelados como contagens com exposição igual a 1. Alguns modelos incluem efeito de coorte. Projeção demográfica por coortes com sobreviventes, nascimentos e fluxos migratórios.
5.  Choques transitórios \[Liu15\]: $`\log m_{x,t}=\alpha_x+\beta_x\kappa_t+N_tJ_{x,t}+\varepsilon_{x,t}`$. J0 usa mesmo perfil beta; J1 usa perfil etário de choque separado; J2 permite variação aleatória com correlação entre idades. Aplicação em títulos vinculados à mortalidade catastrófica.

## Prática em R

- Ajustar SVD de posto 1 à matriz de log-taxas, aplicar normalização e explorar os parâmetros.
- Projetar kappa usando passeio aleatório com deriva; avaliar RMSE em dados separados cronologicamente e construir limites condicionais.
- Comparar com ajuste Poisson por máxima verossimilhança usando StMoMo (que não equivale ao modelo bayesiano de Czado).
- Gerar amostras das equações de estado de Pedroza, experimentar choques J0/J1/J2 e construir uma recursão de coortes simplificada.

## Pontos de rigor

- Uma taxa central por pessoa-ano não equivale à probabilidade de morrer em um ano.
- Primeira componente SVD pode explicar grande parte da variação e ainda prever mal; separar treino e teste por ano.
- Intervalos condicionais a alfa e beta não capturam toda a incerteza de parâmetros, resíduos ou escolha do modelo.
- Os cinco artigos não propõem todos o mesmo estimador nem o mesmo modelo probabilístico.
- A leitura aqui se baseia nos cinco PDFs enviados; os programas R produzidos são didáticos e simulados, não reproduções numéricas das tabelas dos autores.

---

## References

\[Lee92\] R. Lee and L. Carter, “Modeling and forecasting U. S. mortality,” Sep. 01, 1992. doi: [10.1080/01621459.1992.10475265](https://doi.org/10.1080/01621459.1992.10475265).

\[Cza05\] C. Czado, A. Delwarde, and M. Denuit, “Bayesian Poisson log-bilinear mortality projections,” Jun. 24, 2005. doi: [10.1016/J.INSMATHECO.2005.01.001](https://doi.org/10.1016/J.INSMATHECO.2005.01.001).

\[Ped06\] C. Pedroza, “A Bayesian forecasting model: predicting U.S. male mortality.” *Biostatistics*, vol. 7 4, pp. 530–50, Oct. 2006, doi: [10.1093/BIOSTATISTICS/KXJ024](https://doi.org/10.1093/BIOSTATISTICS/KXJ024).

\[Wis15\] A. Wiśniowski, P. W. F. Smith, J. Bijak, J. Raymer, and J. J. Forster, “Bayesian Population Forecasting: Extending the Lee-Carter Method,” *Demography*, vol. 52, pp. 1035–1059, Jun. 2015, doi: [10.1007/s13524-015-0389-y](https://doi.org/10.1007/s13524-015-0389-y).

\[Liu15\] Y. Liu and J. S. Li, “The age pattern of transitory mortality jumps and its impact on the pricing of catastrophic mortality bonds,” Sep. 01, 2015. doi: [10.1016/J.INSMATHECO.2015.05.005](https://doi.org/10.1016/J.INSMATHECO.2015.05.005).
