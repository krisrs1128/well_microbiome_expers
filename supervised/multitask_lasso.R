#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using lasso to model multiple responses in the WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("viridis")
library("glmnet")
source("../dimension_red/prep_tables.R")

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
## Apply parallel lassos
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.02, "rlog" = FALSE)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
y <- scale(processed$bc)
x <- scale(processed$x_seq)

n_lambda <- 50
beta_hats <- array(dim = c(ncol(y), 1 + ncol(x), n_lambda))
for (r in seq_len(ncol(y))) {
  fit <- glmnet(x, y[, r], nlambda = n_lambda, alpha = 0.7)
  beta_hats[r,, ] <- as.matrix(coef(fit))
}

