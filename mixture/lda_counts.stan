/*
 * Latent Dirichlet Allocation
 *
 * This is a mixed-membership analysis, so somewhere in between (discrete)
 * clustering and (continuous) ordination.
 *
 * It's also my favorite method lately, so might as well try it out.
 * reference: https://arxiv.org/abs/1706.04969
 *
 * author: sankaran.kris@gmail.com
 * date: 09/21/2017
 */
data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num species
  int<lower=0> N; // num samples
  int<lower=0> x[N, V]; // species counts for each sample

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma;
}

parameters {
  simplex[K] theta[N]; // topic mixtures
  simplex[V] beta[K]; // word dist for k^th topic
}

model {
  for (i in 1:N) {
    theta[i] ~ dirichlet(alpha);
  }

  for (k in 1:K) {
    beta[k] ~ dirichlet(gamma);
  }

  for (i in 1:N) {
    vector[V] eta;
    eta = beta[1] * theta[i, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[i, k];
    }
    x[i] ~ multinomial(eta);
  }

}

generated quantities {
  int<lower=0> x_sim[N, V]; // simulated species counts, for posterior checking

  for (i in 1:N) {
    vector[V] eta;
    eta = beta[1] * theta[i, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[i, k];
    }
    x_sim[i] = multinomial_rng(eta, sum(x[i]));
  }
}
