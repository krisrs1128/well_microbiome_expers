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
## Load data
###############################################################################
raw <- read_data()
processed <- process_data(raw$seqtab, raw$bc, raw$bc_full, raw$taxa, raw$tree)

###############################################################################
## run and visualize PCA
###############################################################################
samp <- processed$ps %>%
  sample_data() %>%
  select(-id, -number, -gender, -batch, -operator)
combined <- cbind(samp, get_taxa(processed$ps))
pc_res <- prcomp(scale(combined))

## extract scores and join in sample data
scores <- prepare_scores(list(pc_res$x), c("combined")) %>%
  left_join(
    processed$bc %>%
    rownames_to_column("Number")
  )

## extract loadings and join taxa information
loadings_list <- list(
  pc_res$rotation[rownames(pc_res$rotation) %in% names(sample_data(processed$ps)), ],
  pc_res$rotation[rownames(pc_res$rotation) %in% taxa_names(processed$ps), ]
)

seq_fam <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

loadings <- prepare_loadings(loadings_list, c("body_comp", "seq")) %>%
  left_join(seq_fam)

plot_loadings(loadings, pc_res$sdev)
ggsave("../chapter/figure/pca/loadings.png", width = 7.69, height = 5.17)

## and study the scores
plot_scores(scores, "android_lm", "Android FM", pc_res$sdev) +
  scale_color_viridis(
    "Android FM ",
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca/scores_total_lm.png", width = 4.45, height = 2.63)

##  also study scores in relation to overall bacteroides / ruminoccocus ratio
scores <- scores %>%
  left_join(family_means(processed$mseqtab)) %>%
  mutate(
    rl_ratio = tanh(rl_ratio)
  )
plot_scores(scores, "rl_ratio", "Bact. / Rumino. ratio", pc_res$sdev) +
  scale_color_viridis(
    "Bact. / Rumino. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  )
ggsave("../chapter/figure/pca/scores_rl_ratio.png", width = 3.56, height = 2.6)
