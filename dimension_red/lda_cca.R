#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using Latent Dirichlet Allocation + Canonical Correlation
## Analysis simultaneously across count and gaussian tables.
##
## author: sankaran.kris@gmail.com
## date: 10/04/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("rstan")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
library("vegan")
source("prep_tables.R")
source("plot.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## cleaner ggplot theme
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
## read and prepare the data
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.02, "rlog" = FALSE)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
bc_mat <- scale(processed$bc)

K <- 3
L1 <- 3
L2 <- 3

stan_data <- list(
  "n" = nrow(bc_mat),
  "p1" = ncol(processed$x_seq),
  "p2" = ncol(bc_mat),
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "sigma" = 1,
  "tau" = 8,
  "x" = processed$x_seq,
  "y" = bc_mat,
  "id_y" = diag(ncol(bc_mat)),
  "id_k" = diag(K),
  "id_l1" = diag(L1),
  "id_l2" = diag(L2),
  "zeros_k" = rep(0, K),
  "zeros_l1" = rep(0, L1),
  "zeros_l2" = rep(0, L2)
)

m <- stan_model("lda_cca.stan")
vb_fit <- vb(m, data = stan_data)

posterior <- rstan::extract(vb_fit)
pmeans <- parameter_means(posterior)

###############################################################################
## Plot the results
###############################################################################
scv <- scale_color_viridis(guide = guide_colorbar(barwidth = 0.15, ticks = FALSE))

## Plot some of the shared scores
plot_scores_wrapper(pmeans$xi_s, raw, processed, scv)
rownames(pmeans$By) <- colnames(processed$bc)
rownames(pmeans$Bx) <- colnames(processed$x_seq)
loadings <- prepare_loadings(
  list(cbind(100 * pmeans$By, 5 * 1), cbind(pmeans$Bx, 1)),
  c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)
loadings[loadings$type == "body_comp", "family"] <- "Body Comp."

plot_loadings(loadings, c(1, 1)) +
  facet_wrap(~family, ncol = 4) +
  scale_color_brewer(palette = "Set2", guide = FALSE) +
  scale_size_continuous(range = c(0, 2), guide = FALSE)

## now plot unshared scores (species abundances first)
plot_scores_wrapper(pmeans$xi_x, raw, processed, scv)
plot_scores_wrapper(pmeans$xi_y, raw, processed, scv)

## Now plot unshared loadings
rownames(pmeans$Wx) <- colnames(processed$x_seq)
loadings_x <- prepare_loadings(list(cbind(pmeans$Wx, 1)), "seq") %>%
  left_join(processed$mseqtab)
plot_loadings(loadings_x, c(1, 1)) +
  scale_size_continuous(range = c(0, 2), guide = FALSE)

rownames(pmeans$Wy) <- colnames(processed$bc)
loadings_y <- prepare_loadings(list(cbind(pmeans$Wy, 1)), "body_comp") %>%
  mutate(seq_num = "NA") %>%
  left_join(processed$mseqtab)
plot_loadings(loadings_y, c(1, 1)) +
  scale_size_continuous(range = c(1, 4), guide = FALSE)
