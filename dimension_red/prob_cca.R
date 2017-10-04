#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Applying probabilistic CCA to real data.
##
## author: sankaran.kris@gmail.com
## date: 10/04/2017

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

m <- stan_model("prob_cca.stan")
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
