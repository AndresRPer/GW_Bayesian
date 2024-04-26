## Prueba de modelo con 5 tipos ------------------------------------------------

## Compilamos el modelos -------------------------------------------------------
mod_5types <- cmdstan_model("models/GW_Bayesian_5Types_Tests.stan")
#Ya quedo bien

## Datos sintéticos para hacer la prueba ---------------------------------------

#Probs de transición reales
transition_probs <- matrix(c(
  0.2, 0.2, 0.2, 0.2, 0.2,  
  0.2, 0.2, 0.2, 0.2, 0.2,  
  0.2, 0.2, 0.2, 0.2, 0.2,  
  0.2, 0.2, 0.2, 0.2, 0.2,  
  0.2, 0.2, 0.2, 0.2, 0.2  
), byrow = TRUE, nrow = 5)

#Poblaciones iniciales para cada tipo

initial_counts <- c(100,100,100,100,100)

#Generaciones
N <- 10

#Matriz de poblacion
Z <- matrix(nrow = 5, ncol = N)
Z[,1] <- initial_counts

#Generamos los datos
set.seed(182120) #Datos reproducibles
for (n in 1:(N-1)) {
  for (i in 1:5) {
    #Calculamos el total de incividuos del tipo i en la n-ésima generación
    total_individuals_type_i <- Z[i,n]
    
    #El count de cada tipo en la sig. generación es muestreado de una miltinomial
    Z[,n+1] <- rmultinom(1, size = total_individuals_type_i,
                         prob = transition_probs[i,])
  }
}
Z <- t(Z)

## Fit del modelo --------------------------------------------------------------
fit_5tipos <- mod_5types$sample(data = list(N = N, Z = t(Z)),
                                seed = 182120, chains = 1,
                                parallel_chains = 1,
                                iter_warmup = 5000,
                                iter_sampling = 10000)









