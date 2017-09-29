#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## A study of using probabilistic CCA on simulated and real data.
##
## author: sankaran.kris@gmail.com
## date: 09/28/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("rstan")
set.seed(9282017)

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
  panel.border = element_rect(size = 0.5),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

matnorm <- function(n, p, mu = 0, sigma = 1) {
  matrix(
    rnorm(n * p, mu, sigma),
    n, p
  )
}

###############################################################################
## Simulation study to validate approach
###############################################################################

K <- 5
ps <- c(20, 50)
n <- 100

## simulate data
theta <- list(
  "Wx" = matnorm(20, K, 0, 2),
  "Wy" = matnorm(50, K, 0, 2),
  "Bx" = matnorm(20, K),
  "By" = matnorm(50, K),
  "xi_s" = matnorm(n, K),
  "xi_x" = matnorm(n, K),
  "xi_y" = matnorm(n, K)
)

x <- theta$xi_s %*% t(theta$Wx) +
  theta$xi_x %*% t(theta$Bx) +
  matnorm(n, ps[1], 0, 0.1)
y <- theta$xi_s %*% t(theta$Wy) +
  theta$xi_y %*% t(theta$By) +
  matnorm(n, ps[2], 0, 0.1)

## fit model
m <- stan_model("prob_cca.stan")
stan_data <- list(
  "K" = K,
  "p1" = ps[1],
  "p2" = ps[2],
  "n" = n,
  "id1" = diag(ps[1]),
  "id2" = diag(ps[2]),
  "idk" = diag(K),
  "x" = x,
  "y" = y
)

fit <- vb(m, stan_data)

###############################################################################
## Real data application
###############################################################################
