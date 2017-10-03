/*
 * Fit a LDA + CCA model
 *
 * author: sankaran.kris@gmail.com
 * date: 10/03/2017
 */

data {
  int<lower=1> n;
  int<lower=1> p1;
  int<lower=1> p2;
  int<lower=1> K;
  int<lower=1> L1;
  int<lower=1> L2;
  real<lower=0> sigma;
  real<lower=0> tau;

  matrix[n, p1] x;
  matrix[n, p2] y;
  matrix[p2, p2] id_y;
  matrix[K, K] id_k;
  matrix[L1, L1] id_l1;
  matrix[L2, L2] id_l2;
}


parameters {
  real<lower=0> sigma;
  matrix[n, K] xi_s;
  matrix[n, K] xi_x;
  matrix[n, K] xi_y;
  matrix[p1, K] Wx;
  matrix[p2, K] Wy;
  matrix[p1, K] Bx;
  matrix[p2, K] By;
}

transformed parameters {
  matrix[p2, p2] cov_y;
  cov_y = sigma ^ 2 * id_y;
}

model {
  // prior
  for (i in 1:n) {
    xi_s[i] ~ multi_normal(zeros_k, tau * id_k);
    xi_x[i] ~ multi_normal(zeros_k, tau * id_l1);
    xi_y[i] ~ multi_normal(zeros_k, tau * id_l2);
  }

  // likelihood
  for (i in 1:n) {
    x[i] ~ multinomial(softmax(Bx * xi_s[i]' + Wx * xi_x[i]') sum(x[i]));
    y[i] ~ multi_normal(By * xi_s[i]' + Wy * xi_y[i]', cov_y);
  }
}
