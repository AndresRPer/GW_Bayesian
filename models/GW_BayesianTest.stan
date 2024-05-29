data {
  int<lower=1> K; // Number of types (2 in this case)
  int<lower=1> N; // Number of generations
  array[K] int<lower=0> initial_counts; // Initial counts of individuals of each type
  array[N,K] int<lower=0> counts; // Observed counts of individuals at each generation
}

parameters {
  array[K] simplex[K] offspring_distribution; // Offspring distribution for each type
}

model {
  // Non-informative priors for the offspring distribution
  for (k in 1:K) {
    // Non-informative prior for Dirichlet (effectively, Dirichlet(1))
    offspring_distribution[k] ~ dirichlet(rep_vector(1.0, K));
  }

  // Likelihood: simulate the counts of individuals at the next generation
  for (n in 2:N) {
    for (k in 1:K) {
      vector[K] expected_counts = offspring_distribution[k] * counts[n-1, k];
      counts[n, ] ~ multinomial(expected_counts / sum(expected_counts));
    }
  }
}

generated quantities {
  array[K] int final_counts;
  for (k in 1:K) {
    final_counts[k] = counts[N, k];
  }
}
