data {
  int<lower = 0> N; //Number of generations
  array[5,N] int<lower=0> Z; // Observed counts per generation per type
}

parameters {
  simplex[5] p1;
  simplex[5] p2;
  simplex[5] p3;
  simplex[5] p4;
  simplex[5] p5;  //Transition probabilities from type i to type j
}

model {
  //Priors
  p1 ~ dirichlet([1,1,1,1,1]);
  p2 ~ dirichlet([1,1,1,1,1]);
  p3 ~ dirichlet([1,1,1,1,1]);
  p4 ~ dirichlet([1,1,1,1,1]);
  p5 ~ dirichlet([1,1,1,1,1]);
  
  //Likelihood
  for (n in 1:N-1) {
    vector[5] probabilities_for_next_gen;
    array[5] int counts_for_next_gen;
    
    //Calculate probs for the next generation
    for (i in 1:5) {
      probabilities_for_next_gen[i] = p1[i] * Z[1,n] + p2[i] * Z[2,n] +
                                      p3[i] * Z[3,n] + p4[i] * Z[4,n] +
                                      p5[i] * Z[5,n];
      
    } 
    
    //Normalize to make them probablities
    probabilities_for_next_gen /= sum(probabilities_for_next_gen);
    
    //Setting up the counts for the multinomial distribution
    for (i in 1:5) {
      counts_for_next_gen[i] = Z[i, n+1];
    }
    
    //Multinomial distribution with observed counts and calculated probs
    counts_for_next_gen ~ multinomial(probabilities_for_next_gen);
    
  }
}

generated quantities {
  matrix[5,5] P;  //Transition probability matrix
  P[,1] = p1;
  P[,2] = p2;
  P[,3] = p3;
  P[,4] = p4;
  P[,5] = p5;
}
