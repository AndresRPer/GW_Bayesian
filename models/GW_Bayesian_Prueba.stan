data {
  int<lower=1> K;  // Number of types
  int<lower=1> T;  // Number of time points
  int<lower=1> max_kappa;  // Maximum number of offspring any type can produce
  int<lower=0> n[K, K, T+1, max_kappa+1];  // Life table data using max_kappa as a fixed size
  int<lower=1> kappa[K, K];  // Max number of offspring for each pair (i, j), not used for array sizes
}

parameters {
  simplex[max_kappa+1] p[K, K]; // Offspring distribution probabilities for each type pair
}

model {
  for (i in 1:K) {
    for (j in 1:K) {
      vector[max_kappa+1] alpha = rep_vector(1.0, max_kappa + 1);  // Initialize Dirichlet alphas
      for (k in 0:kappa[i,j]) {
        for (t in 0:T) {
          alpha[k+1] += n[i, j, t, k];  // Sum life table counts into alpha
        }
      }
      p[i,j] ~ dirichlet(alpha);  // Posterior Dirichlet distribution
    }
  }
}
