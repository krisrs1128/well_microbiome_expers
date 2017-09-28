#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An example studying multiple tables by concatenating and then using PCA.
##
## author: sankaran.kris@gmail.com
## date: 09/25/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("phyloseq")
library("ggrepel")
library("reshape2")
library("viridis")
source("prep_tables.R")
source("plot.R")

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
## Load data
###############################################################################
raw <- read_data()
opts <- list(filt_k = 0.07, filt_a = 0)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## run and visualize PCA
###############################################################################
combined <- cbind(processed$bc, processed$x_seq)
pc_res <- prcomp(scale(combined))

## extract scores and join in sample data
scores <- prepare_scores(list(pc_res$x), c("combined")) %>%
  left_join(
    processed$bc %>%
    rownames_to_column("Number")
  )

## extract loadings and join taxa information
loadings_list <- list(
  pc_res$rotation[rownames(pc_res$rotation) %in% colnames(processed$bc), ],
  pc_res$rotation[rownames(pc_res$rotation) %in% colnames(processed$x_seq), ]
)

loadings <- prepare_loadings(loadings_list, c("body_comp", "seq")) %>%
  left_join(processed$mseqtab)

plot_loadings(loadings, pc_res$sdev) +
  ylim(-0.15, 0.35) +
  xlim(-0.15, 0.2)
ggsave("../chapter/figure/pca/loadings.png", width = 4.56, height = 3)

## and study the scores
plot_scores(scores, "weight_dxa", "Weight", pc_res$sdev) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca/scores_weight.png", width = 3.56, height = 2.6)

##  also study scores in relation to overall ruminoccocus / lachospiraceae ratio
scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "Rum. / Lach ratio", pc_res$sdev) +
  scale_color_viridis(
    "Rum. / Lach. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca/scores_rl_ratio.png", width = 3.56, height = 2.6)
