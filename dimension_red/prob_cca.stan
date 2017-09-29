/*
 * Fit a probabilistic CCA model
 *
 * author: sankaran.kris@gmail.com
 * date: 09/28/2017
 */

data {
  int<lower=1> K;
  int<lower=1> p1;
  int<lower=1> p2;
  int<lower=1> n;
  matrix[p1, p1] id1;
  matrix[p2, p2] id2;
  matrix[p2, p2] idk;
  matrix[n, p1] x;
  matrix[n, p2] y;
}


parameters {
  real<lower=0> sigma_x;
  real<lower=0> sigma_y;
  matrix[n, K] xi_s;
  matrix[n, K] xi_x;
  matrix[n, K] xi_y;
  matrix[p1, K] Wx;
  matrix[p2, K] Wy;
  matrix[p1, K] Bx;
  matrix[p2, K] By;
}

transformed parameters {
  vector[K] zeros;
  matrix[p1, p1] cov_x;
  matrix[p2, p2] cov_y;

  for (i in 1:p1) {
    for (j in 1:p1) {
      cov_x[i, j] = sigma_x ^ 2;
    }
  }

  for (i in 1:p2) {
    for (j in 1:p2) {
      cov_y[i, j] = sigma_y ^ 2;
    }
  }
}

model {
  // prior
  for (i in 1:n) {
    xi_s[i] ~ multi_normal(zeros, idk);
    xi_x[i] ~ multi_normal(zeros, idk);
    xi_y[i] ~ multi_normal(zeros, idk);
  }

  // likelihood
  for (i in 1:n) {
    x[i] ~ multi_normal(Wx * xi_s[i]' + Bx * xi_x[i]', cov_x);
    y[i] ~ multi_normal(Wy * xi_s[i]' + By * xi_y[i]', cov_y);
  }
}
