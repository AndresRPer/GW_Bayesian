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




 












