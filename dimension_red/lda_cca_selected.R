#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Variant of lda_cca.R that runs the shared LDA scores model only only CCA.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

library("phyloseq")
library("rstan")
library("tidyverse")
library("PMD")
library("reshape2")
library("ggrepel")
library("viridis")
source("prep_tables.R")
source("plot.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## cleaner ggplot theme
scale_colour_discrete <- function(...)
  scale_color_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#464646"),
    na.value = "black"
  )
scale_fill_discrete <- function(...)
  scale_fill_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#464646"),
    na.value = "black"
  )

theme_set(theme_bw())
theme_update(
  panel.background = element_rect(fill = "#F8F8F8"),
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
opts <- list(
  "filt_k" = 0.02,
  "stan_file" = "lda_cca.stan",
  "outdir" = "../chapter/figure/lda_cca_selected/"
)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
bc_mat <- scale(processed$bc)

###############################################################################
## Identify subset of features to run on using CCA
###############################################################################
K <- 2
cca_res <- CCA(bc_mat, scale(processed$x_seq), penaltyx = 0.6, penaltyz = 0.3, K = K)

opts$rlog <- FALSE
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
keep_ix <- rowSums(cca_res$v) == 0

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
L1 <- 3
L2 <- 3

stan_data <- list(
  "n" = nrow(bc_mat),
  "p1" = sum(keep_ix),
  "p2" = ncol(bc_mat),
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "sigma" = 1,
  "a0" = 1,
  "b0" = 1,
  "x" = processed$x_seq[, keep_ix],
  "y" = bc_mat,
  "id_y" = diag(ncol(bc_mat)),
  "id_k" = diag(K),
  "id_l1" = diag(L1),
  "id_l2" = diag(L2),
  "zeros_k" = rep(0, K),
  "zeros_l1" = rep(0, L1),
  "zeros_l2" = rep(0, L2)
)

m <- stan_model(opts$stan_file)
vb_fit <- vb(m, data = stan_data)
