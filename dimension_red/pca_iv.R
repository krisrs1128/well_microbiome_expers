#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Illustration of PCA-IV on the body composition and microbiome elements of the
## WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 09/29/2017

library("phyloseq")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
library("ade4")
source("prep_tables.R")
source("plot.R")

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
## Run PCA-IV
###############################################################################
raw <- read_data()
processed <- process_data(raw$seqtab, raw$bc, raw$taxa)

K <- 3
pca_micro <- dudi.pca(
  scale(processed$x_seq),
  center = FALSE,
  scale = FALSE,
  scan = FALSE,
  nf = 3
)
pcaiv_res <- pcaiv(
  pca_micro, scale(processed$bc),
  scan = FALSE,
  nf = 3
)

###############################################################################
## study equivalent of scores (projections of x and y onto principal axes)
###############################################################################
pcaiv_res
summary(pcaiv_res)
plot(pcaiv_res)

## correlate columns from both data frames with the principal axes
loadings <- prepare_loadings(
  list(0.001 * pcaiv_res$fa, pcaiv_res$c1),
    c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)
plot_loadings(loadings, pcaiv_res$eig) +
  scale_size_continuous(range = c(0.3, 1))
ggsave("../chapter/figure/pca_iv/loadings.png", width = 4.56, height = 3)

## project the samples onto the principal axes
scores <- prepare_scores(
  list(pcaiv_res$li, pcaiv_res$ls),
  c("body_comp", "seq")
) %>%
  left_join(raw$bc)

mscores <- melt_scores(scores)
plot_scores(scores, "Total_LM", "Total LM", pcaiv_res$eig) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca_iv/scores_total_lm.png", width = 3.56, height = 1.8)

scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "Rum. / Lach. ratio", pcaiv_res$eig) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca_iv/scores_rl_ratio.png", width = 3.56, height = 1.8)
