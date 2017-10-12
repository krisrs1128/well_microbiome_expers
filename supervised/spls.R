#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Sparse partial least squares applied to the WELL microbiome data.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("reshape2")
source("../dimension_red/prep_tables.R")
library("spls")

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
opts <- list("filt_k" = 0.5)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

y <- scale(processed$bc)
x <- scale(processed$x_seq)
cv_eval <- cv.spls(x, y, K = 1:6, eta = seq(0, 0.9, 0.05), scale.x = FALSE, fold = 5)
cv_eval

train_ix <- sample(1:nrow(x), 80)
fit <- spls(x[train_ix, ], y[train_ix, ], cv_eval$K.opt, cv_eval$eta.opt)
y_hat <- x %*% fit$betahat
plot(y[train_ix, 24], y_hat[train_ix, 24])
points(y[-train_ix, 24], y_hat[-train_ix, 24], col = "blue")
abline(a = 0, b = 1, col = "red")
