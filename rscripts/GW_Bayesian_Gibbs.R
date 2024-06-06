# Load necessary packages
library(MCMCpack)
library(reshape2)
library(ggplot2)

# Step 1: Simulate the branching process to generate data
simulate_branching <- function(population, P, generations, kappa) {
  M <- nrow(P)
  max_kappa <- max(kappa)
  data <- array(0, dim = c(generations, M, M, max_kappa + 1))
  
  for (gen in 1:generations) {
    new_population <- numeric(M)
    for (i in 1:M) {
      for (j in 1:population[i]) {
        offspring <- rmultinom(1, 1, P[i, ])
        for (k in 1:M) {
          data[gen, i, k, offspring[k] + 1] <- data[gen, i, k, offspring[k] + 1] + 1
        }
        new_population <- new_population + offspring
      }
    }
    population <- new_population
  }
  
  return(data)
}

# Parameters
M <- 3
generations <- 10
population <- c(10, 5, 5)
kappa <- matrix(c(5, 5, 5,
                  5, 5, 5,
                  5, 5, 5), 
                nrow = M, byrow = TRUE)
P_true <- matrix(c(0.3, 0.5, 0.2,
                   0.2, 0.3, 0.5,
                   0.4, 0.1, 0.5), 
                 nrow = M, byrow = TRUE)

data <- simulate_branching(population, P_true, generations, kappa)

# Step 2: Prepare the data for MCMC
counts <- array(0, dim = c(M, M, max(kappa) + 1))
for (gen in 1:generations) {
  for (i in 1:M) {
    for (j in 1:M) {
      counts[i, j, ] <- counts[i, j, ] + data[gen, i, j, ]
    }
  }
}

alpha_prior <- array(1, dim = c(M, M, max(kappa) + 1))  # Prior parameters for Dirichlet distribution

# Step 3: Run the MCMC using MCmultinomialdirichlet
# Run MCMC for each row of P
samples_list <- list()
for (i in 1:M) {
  for (j in 1:M) {
    samples_list[[paste(i, j, sep = "_")]] <- MCmultinomdirichlet(
      counts[i, j, 1:(kappa[i, j] + 1)], 
      alpha_prior[i, j, 1:(kappa[i, j] + 1)], 
      mcmc = 10000, burnin = 1000, thin = 10
    )
  }
}

# Combine samples into an array
samples_array <- array(0, dim = c(1000, M, M, max(kappa) + 1))
for (i in 1:M) {
  for (j in 1:M) {
    samples_array[, i, j, 1:(kappa[i, j] + 1)] <- samples_list[[paste(i, j, sep = "_")]]
  }
}

# Step 4: Analyze the results
# Calculate the mean of the posterior samples
P_estimated <- apply(samples_array, c(2, 3, 4), mean)

cat("Estimated Offspring Distribution Matrix P:\n")
print(P_estimated)

# Obtain marginal distributions and their parameters
marginals <- list()
marginal_parameters <- list()

for (i in 1:M) {
  for (j in 1:M) {
    for (k in 1:(kappa[i, j] + 1)) {
      marginal_name <- paste("P[", i, ",", j, ",", k - 1, "]", sep = "")
      marginals[[marginal_name]] <- samples_array[, i, j, k]
      alpha_post <- alpha_prior[i, j, k] + counts[i, j, k]
      alpha_0 <- sum(alpha_prior[i, j, 1:(kappa[i, j] + 1)]) + sum(counts[i, j, 1:(kappa[i, j] + 1)])
      beta_post <- alpha_0 - alpha_post
      marginal_parameters[[marginal_name]] <- c(alpha_post, beta_post)
    }
  }
}

# Summary statistics for each marginal distribution
marginal_summaries <- lapply(marginals, summary)
marginal_summaries

# Parameters for the marginal Beta distributions
marginal_parameters

# Prepare P_true for comparison
P_true_extended <- array(0, dim = c(M, M, max(kappa) + 1))
for (i in 1:M) {
  for (j in 1:M) {
    P_true_extended[i, j, 1:2] <- c(P_true[i, j], 1 - P_true[i, j])
  }
}

# Compute MAE and MSE
mae <- mean(abs(P_estimated - P_true_extended))
mse <- mean((P_estimated - P_true_extended)^2)
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Mean Squared Error (MSE):", mse, "\n")
