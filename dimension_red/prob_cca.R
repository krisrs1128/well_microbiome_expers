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
library("viridis")
source("prep_tables.R")
source("plot.R")
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

#' Average across first dimension
slice_mean <- function(x) {
  apply(x, c(2, 3), mean)
}

#' Means from posterior samples
parameter_means <- function(theta_samples) {
  theta_hat <- list()
  for (i in seq_along(theta_samples)) {
    if (length(dim(theta_samples[[i]])) > 2) {
      theta_hat[[i]] <- slice_mean(theta_samples[[i]])
    } else {
      theta_hat[[i]] <- mean(theta_samples[[i]])
    }
  }

  names(theta_hat) <- names(theta_samples)
  theta_hat
}

###############################################################################
## Simulation study to validate approach
###############################################################################
K <- 3
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
  "tau" = 5,
  "x" = x,
  "y" = y,
  "id_x" = diag(ps[1]),
  "id_y" = diag(ps[2]),
  "id_k" = diag(K),
  "zeros_k" = rep(0, K)
)

fit <- vb(m, stan_data, iter = 5000)
theta_samples <- rstan::extract(fit)
rm(fit)
theta_hat <- parameter_means(theta_samples)
plot(theta_hat$xi_s[, 1], theta$xi_s[, 2], asp = 1)

###############################################################################
## Real data application
###############################################################################

## read and prepare data
raw <- read_data()
processed <- process_data(raw$seqtab, raw$bc, raw$taxa)

stan_data <- list(
  "K" = K,
  "p1" = ncol(processed$bc),
  "p2" = ncol(processed$x_seq),
  "n" = nrow(processed$bc),
  "tau" = 5,
  "x" = scale(processed$bc),
  "y" = scale(processed$x_seq),
  "id_x" = diag(ncol(processed$bc)),
  "id_y" = diag(ncol(processed$x_seq)),
  "id_k" = diag(K),
  "zeros_k" = rep(0, K)
)

micro_fit <- vb(m, stan_data, iter = 100)
micro_samples <- rstan::extract(micro_fit)
rm(micro_fit)

micro_hat <- parameter_means(micro_samples)
shared_scores <- cbind(micro_hat$xi_s, processed$bc)
colnames(shared_scores)[1:K] <- paste0("Axis", 1:K)

## should try to plot posteriors, not just means
ggplot(shared_scores) +
  geom_point(
    aes(x = Axis1, y = Axis2, col = weight_dxa)
  ) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  )

rownames(micro_hat$Wx) <- colnames(processed$bc)
rownames(micro_hat$Wy) <- colnames(processed$x_seq)
rownames(micro_hat$Bx) <- colnames(processed$bc)
rownames(micro_hat$By) <- colnames(processed$x_seq)

loadings <- prepare_loadings(
  list(micro_hat$Wx, micro_hat$Wy),
  c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)
plot_loadings(loadings, c(1, 1))

loadings_distinct <- prepare_loadings(
  list(micro_hat$Bx, micro_hat$By),
  c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)
plot_loadings(loadings_distinct, c(1, 1))
