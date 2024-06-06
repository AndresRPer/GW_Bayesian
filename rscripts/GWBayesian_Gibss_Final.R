# Load necessary packages
library(MCMCpack)
library(reshape2)
library(ggplot2)

# Define three P_true matrices for each parent type with probabilities for up to 5 offspring
M <- 3
max_kappa <- 5
P_true_list <- list(
  matrix(c(0.3, 0.1, 0.1, 0.1, 0.1, 0.3,
           0.2, 0.1, 0.1, 0.1, 0.1, 0.4,
           0.4, 0.1, 0.1, 0.1, 0.1, 0.2), nrow = M, byrow = TRUE),
  matrix(c(0.5, 0.1, 0.1, 0.1, 0.1, 0.1,
           0.3, 0.1, 0.1, 0.1, 0.1, 0.3,
           0.1, 0.1, 0.1, 0.1, 0.1, 0.5), nrow = M, byrow = TRUE),
  matrix(c(0.2, 0.1, 0.1, 0.1, 0.1, 0.4,
           0.5, 0.1, 0.1, 0.1, 0.1, 0.1,
           0.3, 0.1, 0.1, 0.1, 0.1, 0.3), nrow = M, byrow = TRUE)
)

# Simulate the branching process to generate data
simulate_branching <- function(population, P_true_list, generations, max_kappa) {
  M <- length(population)
  data <- array(0, dim = c(generations, M, M, max_kappa + 1))
  
  for (gen in 1:generations) {
    new_population <- numeric(M)
    for (i in 1:M) {
      for (j in 1:population[i]) {
        for (k in 1:M) {
          offspring <- rmultinom(1, 1, P_true_list[[i]][k, ])
          data[gen, i, k, ] <- data[gen, i, k, ] + offspring
        }
        new_population <- new_population + colSums(offspring)
      }
    }
    population <- new_population
  }
  
  return(data)
}

# Parameters
generations <- 10
population <- c(10, 5, 5)

# Generate synthetic data
data <- simulate_branching(population, P_true_list, generations, max_kappa)

# Prepare the data for MCMC
counts <- array(0, dim = c(M, M, max_kappa + 1))
for (gen in 1:generations) {
  for (i in 1:M) {
    for (j in 1:M) {
      counts[i, j, ] <- counts[i, j, ] + data[gen, i, j, ]
    }
  }
}

alpha_prior <- array(0.2, dim = c(M, M, max_kappa + 1))  # Prior parameters for Dirichlet distribution

# Run the MCMC using MCmultinomialdirichlet
samples_list <- list()
for (i in 1:M) {
  for (j in 1:M) {
    samples_list[[paste(i, j, sep = "_")]] <- MCmultinomdirichlet(
      counts[i, j, 1:(max_kappa + 1)], 
      alpha_prior[i, j, 1:(max_kappa + 1)], 
      mcmc = 10000, burnin = 1000, thin = 10
    )
  }
}

# Combine samples into an array
samples_array <- array(0, dim = c(1000, M, M, max_kappa + 1))
for (i in 1:M) {
  for (j in 1:M) {
    samples_array[, i, j, 1:(max_kappa + 1)] <- samples_list[[paste(i, j, sep = "_")]]
  }
}

# Calculate the mean of the posterior samples
P_estimated <- apply(samples_array, c(2, 3, 4), mean)

cat("Estimated Offspring Distribution Matrix P:\n")


# Rearrange P_estimated into three matrices of size M x (max_kappa + 1)
P_rearranged <- list()
for (i in 1:M) {
  P_rearranged[[i]] <- matrix(0, nrow = M, ncol = max_kappa + 1)
  for (j in 1:M) {
    for (k in 1:(max_kappa + 1)) {
      P_rearranged[[i]][j, k] <- P_estimated[i, j, k]
    }
  }
}
print(P_rearranged)

# Calculate the MSE for each parent type
mse_values <- numeric(M)
for (i in 1:M) {
  mse_values[i] <- mean((P_rearranged[[i]] - P_true_list[[i]])^2)
}

# Print the MSE values
for (i in 1:M) {
  cat(sprintf("MSE for Parent Type %d: %f\n", i, mse_values[i]))
}

# Visualize the rearranged P_estimated
for (i in 1:M) {
  P_melted <- melt(P_rearranged[[i]])
  colnames(P_melted) <- c("OffspringType", "OffspringCount", "Probability")
  P_melted$ParentType <- i
  ggplot(P_melted, aes(x = OffspringCount - 1, y = Probability, fill = factor(OffspringType))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = sprintf("Offspring Distribution for Parent Type %d", i),
         x = "Number of Offspring", y = "Probability", fill = "Offspring Type") +
    theme_minimal() +
    print()
}
