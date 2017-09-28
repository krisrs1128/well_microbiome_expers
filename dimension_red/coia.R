#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Example using Co-Inertia analysis.
##
## author: sankaran.kris@gmail.com
## date: 09/27/2017

library("phyloseq")
library("DESeq2")
library("tidyverse")
library("forcats")
library("reshape2")
library("ggrepel")
library("viridis")
library("ade4")

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

cca_perc <- function(cca_res, i) {
  round(100 * cca_res$CanCorr[i] / sum(cca_res$CanCorr), 2)
}

perc_label <- function(cca_res, i) {
  sprintf("CC%s [%s%%]", i, cca_perc(cca_res, i))
}

###############################################################################
## Load data
###############################################################################
raw <- read_data()
opts <- list(filt_k = 0.5, filt_a = 0)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## Run CoIA on the two (scaled) tables
###############################################################################
dudi1 <- dudi.pca(x_seq, scan = FALSE, nf = 3)
dudi2 <- dudi.pca(bc_mat, scan = FALSE, nf = 3)
coin1 <- coinertia(dudi1, dudi2, scan = FALSE, nf = 3)
