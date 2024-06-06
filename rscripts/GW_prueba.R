##Necessary libraries
library(cmdstanr)
library(posterior)
library(MCMCpack)

#Paremeters of synthethic data
set.seed(182120)
K <- 2 #Types
N <- 10 #Generations
initial_counts <- c(100,100) #Initial population
alpha <- matrix(c(0.5, 0.5,
                0.5, 0.5), nrow = K, ncol = K, byrow = TRUE)

#Function to generate synthetic data
generetate_synthethic_data <- function(K, N, initial_counts, alpha) {
  offspring_distribution <- matrix(NA, nrow = K, ncol = K)
  for (k in 1:K) {
    offspring_distribution[k, ] <- MCMCpack::rdirichlet(1, alpha[k, ])
  }
  
  counts <- matrix(0, nrow = N, ncol = K)
  counts[1, ] <- initial_counts
  
  for (n in 2:N) {
   new_counts <- rep(0, K)
   for (k in 1:K){
     new_counts <- new_counts + rmultinom(1, counts[n - 1, k], offspring_distribution[k, ])
   }
   counts[n, ] <- new_counts
  }
  
  list(
    K = K,
    N = N,
    initial_counts = initial_counts,
    counts = counts
  )
}

#Generate synthetic data 
synthethic_data <- generetate_synthethic_data(K, N, initial_counts, alpha)

#Compile Stan model
mod_test <- cmdstan_model("models/GW_BayesianTest.stan")

#Fit Stan model to synthetic data
fit_test <- mod_test$sample(data = synthethic_data, seed = 182120,
                            chains = 4, parallel_chains = 4,
                            iter_warmup = 5000, iter_sampling = 100000,
                            adapt_delta = 0.9)
#Print summary of the fit
print(fit_test)

#Extract posterior samples
posterior_samples <- fit_test$draws()
posterior_summary <- summarise_draws(posterior_samples)

#Print summary of posterior distribution
print(posterior_summary)

# Check trace plots
mcmc_trace(posterior_samples, pars = c("offspring_distribution[1,1]",
                                       "offspring_distribution[1,2]", 
                                       "offspring_distribution[2,1]", 
                                       "offspring_distribution[2,2]"))
