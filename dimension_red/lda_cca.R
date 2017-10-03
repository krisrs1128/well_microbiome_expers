#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using Latent Dirichlet Allocation + Canonical Correlation
## Analysis simultaneously across count and gaussian tables.
##
## author: sankaran.kris@gmail.com
## date: 10/03/2017

###############################################################################
## Libraries and setup
###############################################################################
library("rstan")
library("tidyverse")
library("reshape2")
set.seed(1032017)

#' i.i.d. Gaussian Matrix
matnorm <- function(n, p, mu = 0, sigma = 1) {
  matrix(
    rnorm(n * p, mu, sigma),
    n, p
  )
}

#' Simulate parameters used in LDA + CCA model
simulate_parameters <- function(n = 100, N = 1000, p1 = 200, p2 = 30, K = 1,
                                L1 = 2, L2 = 2) {
  list(
    "n" = n,
    "p1" = p1,
    "p2" = p2,
    "N" = 1000,
    "z" = matnorm(n, K),
    "xi_x" = matnorm(n, L1),
    "xi_y" = matnorm(n, L2),
    "A" = matnorm(p1, K),
    "B" = matnorm(p2, K),
    "C" = matnorm(p1, L2),
    "D" = matnorm(p2, L1),
    "sigma" = 0.2
  )
}

#' Softmax function
softmax <- function(x) {
  exp(x) / sum(exp(x))
}

#' Simulate from LDA + CCA model
simulate_data <- function(theta) {
  mu_x <- theta$z %*% t(theta$A) + theta$xi_x %*% t(theta$C)
  mu_y <- theta$z %*% t(theta$B) + theta$xi_y %*% t(theta$D)

  X <- matrix(nrow = theta$n, ncol = theta$p1)
  Y <- matrix(nrow = theta$n, ncol = theta$p2)

  for (i in seq_len(theta$n)) {
    X[i, ] <- rmultinom(1, theta$N, softmax(mu_x[i, ]))
  }

  for (i in seq_len(theta$n)) {
    for (j in seq_len(theta$p2)) {
      Y[i, j] <- rnorm(1, mu_y[i, j], theta$sigma)
    }
  }

  list("X" = X, "Y" = Y)
}

###############################################################################
## simulate toy data
###############################################################################
theta <- simulate_parameters()
sim <- simulate_data(theta)
plot(sim$X[, 121], sim$Y[, 11])
cormat <- cor(sim$X, sim$Y)
heatmap(cormat)
