#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Application of the canonical correspondence analysis to the body composition
## data.
##
## author: sankaran.kris@gmail.com
## date: 10/06/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("vegan")
library("PMA")
library("viridis")
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
## Apply canonical correspondence analysis
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.02, "rlog" = FALSE)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
x <- scale(processed$bc)
y <- processed$x_seq
cca_res <- cca(
  y ~ age + weight_dxa + Total_FM + Total_LM + aoi,
  data.frame(x)
)

###############################################################################
## Plot loadings and scores
###############################################################################
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

## Plot the scores
cc_scores <- prepare_scores(
  scores(cca_res, choices = 1:3)[2],
  "body_comp"
) %>%
  left_join(raw$bc)

plot_scores(cc_scores, "Total_LM", "Total LM", cca_res$CCA$eig) +
  scale_color_viridis(
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/ccpna/scores_total_lm.png", width = 5.47, height = 3.19)

loadings <- prepare_loadings(
  list(4 * cca_res$CCA$biplot, cca_res$CCA$v),
  c("body_comp", "seq")
) %>%
  left_join(seq_families)
plot_loadings(loadings, cca_res$CCA$eig, c(-8, 4)) +
  scale_size_continuous(range = c(0.1, 5))
ggsave("../chapter/figure/ccpna/loadings.png", width = 9.81, height = 4.62)
