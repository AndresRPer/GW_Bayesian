data {
    int<lower=0> N; // Number of generations
    array[2,N] int<lower=1> Z; // Observed counts of each type at each generation
}

parameters {
    simplex[2] p1; // Transition probabilities for type 1
    simplex[2] p2; // Transition probabilities for type 2
}

model {
    // Priors
    p1 ~ dirichlet([1, 1]); // Adjust these parameters based on your prior knowledge
    p2 ~ dirichlet([1, 1]); // Adjust these parameters based on your prior knowledge

    // Likelihood
    for (n in 1:N-1) {
        vector[2] probabilities_for_next_gen;
        array[2] int counts_for_next_gen;

        probabilities_for_next_gen[1] = p1[1] * Z[1, n] + p2[1] * Z[2, n];
        probabilities_for_next_gen[2] = p1[2] * Z[1, n] + p2[2] * Z[2, n];
        
        // Normalize to make them probabilities
        probabilities_for_next_gen /= sum(probabilities_for_next_gen);

        // Setting up the counts for the multinomial distribution
        counts_for_next_gen[1] = Z[1, n+1];
        counts_for_next_gen[2] = Z[2, n+1];

        // Multinomial distribution with observed counts and calculated probabilities
        counts_for_next_gen ~ multinomial(probabilities_for_next_gen);
    }
}

generated quantities {
    matrix[2,2] P;
    P[1,1] = p1[1];
    P[1,2] = p1[2];
    P[2,1] = p2[1];
    P[2,2] = p2[2];
}
