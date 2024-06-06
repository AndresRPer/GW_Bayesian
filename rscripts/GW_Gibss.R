# Load necessary packages
library(MCMCpack)

# Step 1: Simulate the branching process to generate data
simulate_branching <- function(population, P, generations) {
  M <- nrow(P)
  data <- array(0, dim = c(generations, M, M))
  
  for (gen in 1:generations) {
    new_population <- numeric(M)
    for (i in 1:M) {
      offspring <- rmultinom(1, population[i], P[i, ])
      new_population <- new_population + offspring
      data[gen, i, ] <- offspring
    }
    population <- new_population
  }
  
  return(data)
}

# Parameters
M <- 3
generations <- 10
population <- c(10, 5, 5)
P_true <- matrix(c(0.3, 0.5, 0.2,
                   0.2, 0.3, 0.5,
                   0.4, 0.1, 0.5), 
                 nrow = M, byrow = TRUE)

data <- simulate_branching(population, P_true, generations)

# Step 2: Prepare the data for MCMC
counts <- apply(data, 2:3, sum)  # Summarize counts over generations
alpha_prior <- matrix(1, nrow = M, ncol = M)  # Prior parameters for Dirichlet distribution

# Step 3: Run the MCMC using MCmultinomialdirichlet
# Run MCMC for each row of P
samples_list <- lapply(1:M, function(i) {
  MCmultinomdirichlet(counts[i, ], alpha_prior[i, ], mcmc = 10000, burnin = 1000, thin = 10)
})

# Combine samples into an array
samples_array <- array(0, dim = c(1000, M, M))
for (i in 1:M) {
  samples_array[, i, ] <- samples_list[[i]]
}

# Step 4: Analyze the results
# Calculate the mean of the posterior samples
P_estimated <- apply(samples_array, c(2, 3), mean)

cat("Estimated Offspring Distribution Matrix P:\n")
print(P_estimated)

