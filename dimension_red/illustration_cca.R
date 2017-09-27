#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An illustration of the geometric view of CCA using data from the WELL study.
##
## author: sankaran.kris@gmail.com
## date: 09/27/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("plotly")
library("plot3D")

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
## Read and prepare data
###############################################################################
bc <- readRDS("../data/sample_data_bc.rds")
counts <- readRDS("../data/seqtab.rds")

bc_sub <- scale(bc[2:4, c("weight_dxa", "Android_LM", "Total_FM")])
counts_sub <- scale(counts[2:4, 1:3])

x_grid <- seq(-0.5, 0.5, 0.1)
xyz_grid <- as.matrix(expand.grid(x_grid, x_grid, x_grid))
bc_plane <- matrix(nrow = nrow(xyz_grid), ncol = 3)
counts_plane <- matrix(nrow = nrow(xyz_grid), ncol = 3)
for (i in seq_len(nrow(xyz_grid))) {
  bc_plane[i, ] <- as.numeric(t(bc_sub) %*% xyz_grid[i, ])
  counts_plane[i, ] <- as.numeric(t(counts_sub) %*% xyz_grid[i, ])
}

## span of the body measurement data
points3D(
  x = bc_plane[, 1],
  y = bc_plane[, 2],
  z = bc_plane[, 3],
  ticktype = "detailed",
  xlab = rownames(bc_sub)[1],
  ylab = rownames(bc_sub)[2],
  zlab = rownames(bc_sub)[3],
  col = "#ff9966",
  alpha = 0.1,
  theta = 40,
  phi = 20,
  cex = 0.5
)

## span of the counts
points3D(
  x = counts_plane[, 1],
  y = counts_plane[, 2],
  z = counts_plane[, 3],
  ticktype = "detailed",
  xlab = rownames(bc_sub)[1],
  ylab = rownames(bc_sub)[2],
  zlab = rownames(bc_sub)[3],
  col = "#ccb3ff",
  alpha = 0.1,
  theta = 40,
  phi = 20,
  cex = 0.5,
  add = TRUE
)

## labels for the body measurements
text3D(
  bc_sub[, 1],
  bc_sub[, 2],
  bc_sub[, 3],
  col = "#ff9966",
  labels = colnames(bc_sub),
  add = TRUE
)

## labels for the otu abundances
text3D(
  counts_sub[, 1],
  counts_sub[, 2],
  counts_sub[, 3],
  col = "#661aff",
  labels = colnames(counts_sub),
  add = TRUE
)
