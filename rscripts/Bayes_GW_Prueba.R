## Setup------------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(scales)
library(ggplot2)
library(tibble)
library(rstan)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(coda)

## Cambia el default del tamaño de fuente --------------------------------------
theme_set(theme_linedraw(base_size = 14))

##Cambia el número de decimales ------------------------------------------------
options(digits = 4)

##Problemas con la Emacs
options(pillar.subtle = FALSE)
options(rlang_backtrace_on_error = "none")
options(crayon.enabled = FALSE)

##Para el tema del ggplot
sin_lineas <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())
sin_leyenda <- theme(legend.position = "none")
sin_ejes <- theme(axis.ticks = element_blank(), axis.text = element_blank())

## Llamamos el modelo ----------------------------------------------------------

mod_prueba <- cmdstan_model("models/GW_Bayesian_Prueba.stan")
#Parece que compila bien :)

##Prueba con datos sintéticos --------------------------------------------------
#Si esto nos sale ya estamos del otros lado, solo habría que ajustar la manera en la que aparecen los datos


#P's verdaderas
p1_true <- c(0.7, 0.3)
p2_true <- c(0.4, 0.6)

#Num generaciones
N <- 10

#Población inicial
initial_population <- c(Type1 = 100, Type2 = 100)

#Matriz para establecer los counts de cada tipo en cada generación

Z <- matrix(nrow = 2, ncol = N)
Z[, 1] <- initial_population

#Generamos datos sintéticos
set.seed(182120) #Resultados reproducibles
for (n in 1:(N - 1)) {
  #Counts en la genreación actual
  current_counts <- Z[,n]
  
  #Counts esperados de acuerdo con las probs de transición
  expected_counts <- matrix(c(p1_true[1] * current_counts[1] + p2_true[1] * current_counts[2],
                              p1_true[2] * current_counts[1] + p2_true[2] * current_counts[2]),
                            nrow = 2, byrow = TRUE)
  
  #Distribución multinomial para generar los counts en la sig generación
  Z[, n+1] <- rmultinom(1, sum(current_counts), c(expected_counts[1, ] / sum(expected_counts), expected_counts[2, ] / sum(expected_counts)))
}

#Transponemos para el modelo de Stan acepte los datos
Z <- base::t(Z)

##Fit del modelo (QUE SALGA POR FAVOR) ------------------------------

#Prep de los datos
stan_data_prueba <- list(N = N, Z = t(Z))

#Fit
fit_prueba <- mod_prueba$sample(data = stan_data_prueba, chains = 4,
                                parallel_chains = 4, iter_sampling = 10000,
                                iter_warmup = 5000)
#Resumen del fit
print(fit_prueba$summary())
print(fit_prueba$diagnostic_summary())

#Veamos las posteriores marginales
draws <- as_draws_df(fit_prueba)

color_scheme_set("purple")
mcmc_dens(draws, pars = c("p1[1]","p1[2]"))
color_scheme_set("teal")
mcmc_dens(draws, pars = c("p2[1]","p2[2]"))

#Análisis de las p's
draws_mcmc <- as.mcmc(draws)
mcmc_list_prueba <- as.mcmc.list(draws_mcmc)

#Resumen
summary(mcmc_list_prueba)

#Plots importantes #OJO: Solo nos fijamos en las p's
plot(draws_mcmc, trace = TRUE, density = TRUE, smooth = TRUE,
     auto.layout = TRUE, ask = dev.interactive())
