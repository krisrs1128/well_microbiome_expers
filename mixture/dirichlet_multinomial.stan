/*
 * Dirichlet Multinomial Clustering
 *
 * This approach to clustering counts models variation in depths while
 * essentially focusing on relative abundances, so might be a good compromise
 * between using relative abundances vs. raw data.
 *
 * reference: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030126
 *
 * author: sankaran.kris@gmail.com
 * date: 09/20/2017
 */

data {
  int<lower=0> N; // Number of samples
  int<lower=0> V; // Number of species
  int<lower=0> K; // Number of clusters
  int<lower=0> x[N, V]; // species abundances

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] beta;
}

parameters {
  simplex[V] p[K]; // cluster probabilities
  simplex[K] thetas; // overall mixture proportions
}

model {
  // priors
  thetas ~ dirichlet(alpha);
  for (k in 1:K) {
    p[k] ~ dirichlet(beta);
  }

  // likelihood
  for (i in 1:N) {
    real probs[K];
    for (k in 1:K) {
      probs[k] = log(thetas[k]) + multinomial_lpmf(x[i] | p[k]);
    }

    target += log_sum_exp(probs);
  }
}
