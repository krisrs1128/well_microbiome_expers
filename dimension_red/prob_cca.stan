/*
 * Fit a probabilistic CCA model
 *
 * author: sankaran.kris@gmail.com
 * date: 09/28/2017
 */

data {
  int<lower=1> K;
  int<lower=1> L1;
  int<lower=1> L2;
  int<lower=1> p1;
  int<lower=1> p2;
  int<lower=1> n;
  real<lower=0> tau;
  matrix[n, p1] x;
  matrix[n, p2] y;
  matrix[p1, p1] id_x;
  matrix[p2, p2] id_y;
  matrix[K, K] id_k;
  vector<lower=0, upper=0>[K] zeros_k;
}


parameters {
  real<lower=0> sigma_x;
  real<lower=0> sigma_y;
  matrix[n, K] xi_s;
  matrix[n, L1] xi_x;
  matrix[n, L2] xi_y;
  matrix[p1, K] Bx;
  matrix[p2, K] By;
  matrix[p1, L1] Wx;
  matrix[p2, L2] Wy;
}

transformed parameters {
  matrix[p1, p1] cov_x;
  matrix[p2, p2] cov_y;
  cov_x = sigma_x ^ 2 * id_x;
  cov_y = sigma_y ^ 2 * id_y;
}

model {
  // prior
  for (i in 1:n) {
    xi_s[i] ~ multi_normal(zeros_k, tau * id_k);
    xi_x[i] ~ multi_normal(zeros_k, tau * id_k);
    xi_y[i] ~ multi_normal(zeros_k, tau * id_k);
  }

  // likelihood
  for (i in 1:n) {
    x[i] ~ multi_normal(Bx * xi_s[i]' + Wx * xi_x[i]', cov_x);
    y[i] ~ multi_normal(By * xi_s[i]' + Wy * xi_y[i]', cov_y);
  }
}
