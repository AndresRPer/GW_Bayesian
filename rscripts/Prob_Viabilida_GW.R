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

alpha_prior <- array(1, dim = c(M, M, max_kappa + 1))  # Prior parameters for Dirichlet distribution

# Run the MCMC using MCmultinomialdirichlet
samples_list <- list()
for (i in 1:M) {
  for (j in 1:M) {
    samples <- MCmultinomdirichlet(
      counts[i, j, 1:(max_kappa + 1)], 
      alpha_prior[i, j, 1:(max_kappa + 1)], 
      mcmc = 10000, burnin = 1000, thin = 10
    )
    samples_list[[paste(i, j, sep = "_")]] <- samples
  }
}

# Combine samples into an array
num_samples <- nrow(samples_list[[1]])
samples_array <- array(0, dim = c(num_samples, M, M, max_kappa + 1))
for (i in 1:M) {
  for (j in 1:M) {
    samples_array[, i, j, 1:(max_kappa + 1)] <- samples_list[[paste(i, j, sep = "_")]]
  }
}

# Calculate the mean of the posterior samples
P_estimated <- apply(samples_array, c(2, 3, 4), mean)

cat("Estimated Offspring Distribution Matrix P:\n")
print(P_estimated)

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

# Calculate the mean matrix
calculate_mean_matrix <- function(P_matrix, max_kappa) {
  mean_matrix <- matrix(0, nrow = nrow(P_matrix), ncol = nrow(P_matrix))
  for (i in 1:nrow(P_matrix)) {
    for (j in 1:nrow(P_matrix)) {
      mean_matrix[i, j] <- sum((0:max_kappa) * P_matrix[j, ])
    }
  }
  return(mean_matrix)
}

# Simulate the process and estimate super-criticality probability
m <- 1000  # Number of simulations
super_critical_count <- 0

for (sim in 1:m) {
  # Simulate d^2 posterior offspring laws
  P_simulated <- list()
  for (i in 1:M) {
    P_simulated[[i]] <- matrix(0, nrow = M, ncol = max_kappa + 1)
    for (j in 1:M) {
      sample_idx <- sample(1:num_samples, 1)  # Sample from posterior
      P_simulated[[i]][j, ] <- samples_array[sample_idx, i, j, ]
    }
  }
  
  # Calculate the mean matrix for the simulated P
  mean_matrix_sim <- matrix(0, nrow = M, ncol = M)
  for (i in 1:M) {
    mean_matrix_sim <- mean_matrix_sim + calculate_mean_matrix(P_simulated[[i]], max_kappa)
  }
  mean_matrix_sim <- mean_matrix_sim / M
  
  # Calculate the largest eigenvalue
  eigenvalues_sim <- eigen(mean_matrix_sim)$values
  largest_eigenvalue_sim <- max(Re(eigenvalues_sim))
  
  # Check if the largest eigenvalue is greater than 1
  if (largest_eigenvalue_sim > 1) {
    super_critical_count <- super_critical_count + 1
  }
}

# Calculate the probability of super-criticality
prob_super_critical <- super_critical_count / m
cat("Estimated Probability of Super-Criticality: ", prob_super_critical, "\n")

# Prepare P_true for comparison
P_true_extended <- list()
for (i in 1:M) {
  P_true_extended[[i]] <- matrix(0, nrow = M, ncol = max_kappa + 1)
  for (j in 1:M) {
    P_true_extended[[i]][j, 1:(max_kappa + 1)] <- P_true_list[[i]][j, 1:(max_kappa + 1)]
  }
}

# Convert to data frames for plotting
P_estimated_df <- data.frame()
P_true_df <- data.frame()

for (i in 1:M) {
  for (j in 1:M) {
    for (k in 0:max_kappa) {
      P_estimated_df <- rbind(P_estimated_df, data.frame(
        ParentType = i,
        OffspringType = j,
        OffspringCount = k,
        Probability = P_rearranged[[i]][j, k + 1],
        Type = "Estimated"
      ))
      P_true_df <- rbind(P_true_df, data.frame(
        ParentType = i,
        OffspringType = j,
        OffspringCount = k,
        Probability = P_true_extended[[i]][j, k + 1],
        Type = "True"
      ))
    }
  }
}

# Combine the data frames
comparison_df <- rbind(P_estimated_df, P_true_df)

# Plot the comparison
ggplot(comparison_df, aes(x = factor(OffspringCount), y = Probability, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ ParentType, labeller = label_both) +
  labs(x = "Número de Descendientes", y = "Probabilidad", 
       title = "Comparación de la Distribución de Descendencia Estimada y Verdadera") +
  theme_minimal()
