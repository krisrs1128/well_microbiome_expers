#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Illustration of canonical correlation analysis on the body composition and
## microbiome elements of the WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 09/26/2017

library("phyloseq")
library("tidyverse")
library("reshape2")
library("ggrepel")
library("viridis")
library("vegan")
source("prep_tables.R")

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
###############################################################################
## read and prepare the data
###############################################################################
raw <- read_data()
processed <- process_data(raw$seqtab, raw$bc, raw$taxa)

###############################################################################
## Run and plot CCA on the two (scaled) tables
###############################################################################
bc_mat <- scale(processed$bc)
x_seq <- scale(processed$x_seq)
cca_res <- CCorA(bc_mat, x_seq)

loadings <- prepare_loadings(
  list(cca_res$corr.Y.Cy, cca_res$corr.X.Cx),
  c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)

plot_loadings(loadings, cca_res$Eigenvalues) +
  ylim(-0.5, 0.4) +
  xlim(-0.9, 0.3)
ggsave("../chapter/figure/cca/loadings.png", width = 4.56, height = 2.3)

###############################################################################
## Plot the scores
###############################################################################

## extract scores and join in sample data
scores <- prepare_scores(
  list(cca_res$Cx, cca_res$Cy),
  c("body_comp", "seq")
) %>%
  left_join(
    processed$bc %>%
    add_rownames("Number")
  )

## color by data type
mscores <- melt_scores(scores)
plot_scores(scores, "type", "Meas. Type", cca_res$Eigenvalues) +
  link_scores(mscores) +
  scale_color_brewer(palette = "Set1")
ggsave("../chapter/figure/cca/scores_linked.png", width = 3.56, height = 2.6)

## color by weight
plot_scores(scores, "weight_dxa", "Weight", cca_res$Eigenvalues) +
  link_scores(mscores) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/cca/scores_weight.png", width = 3.56, height = 2.6)

## color by ruminoccocus / lachnospiraceae ratios
scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "Rum. / Lach ratio", cca_res$Eigenvalues) +
  link_scores(mscores) +
  scale_color_viridis(
    "Rum. / Lach. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/cca/scores_rl_ratio.png", width = 3.56, height = 2.6)
