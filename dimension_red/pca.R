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

pc_perc <- function(pc_res, i) {
  round(100 * pc_res$sdev[i] / sum(pc_res$sdev), 2)
}

perc_label <- function(pc_res, i) {
  sprintf("PC%s [%s%%]", i, pc_perc(pc_res, i))
}

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
scores <- data.frame(
  Number = rownames(combined),
  pc_res$x
) %>%
  left_join(
    data.frame(
      Number = rownames(processed$bc),
      processed$bc
    )
  )

## extract loadings and join taxa information
loadings <- data.frame(
  "variable" = colnames(combined),
  pc_res$rotation
  ) %>%
  mutate(
    type = ifelse(variable %in% colnames(processed$x_seq), "seq", "body_comp"),
    seq_num = variable
  ) %>%
  left_join(processed$mseqtab)
loadings[loadings$type != "seq", "seq_num"] <- NA

## how to species data relate to body composition?
asp_ratio <- sqrt(pc_res$sdev[2] / pc_res$sdev[1])
ggplot(loadings) +
  geom_hline(yintercept = 0, size = 0.5) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = loadings %>%
      filter(type == "seq"),
    aes(x = PC1, y = PC2, size = PC3, col = family),
    alpha = 0.6
  ) +
  geom_text_repel(
    data = loadings %>%
      filter(type == "body_comp"),
    aes(x = PC1, y = PC2, label = variable, size = PC3),
    segment.size = 0.3,
    segment.alpha = 0.5
  ) +
  labs(
    "x" = perc_label(pc_res, 1),
    "y" = perc_label(pc_res, 2),
    "col" = "Family"
  ) +
  scale_size_continuous(range = c(0, 2.5), breaks = c(-0.1, 0.1)) +
  ylim(-0.1, 0.2) +
  xlim(-0.15, 0.15) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/pca/loadings.png", width = 4.56, height = 3.78)

## and study the scores
ggplot(scores) +
  geom_point(
    aes(x = PC1, y = PC2, size = PC3, col = weight_dxa)
  ) +
  labs(
    "x" = perc_label(pc_res, 1),
    "y" = perc_label(pc_res, 2),
    "col" = "Family"
  ) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  ) +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(ratio = asp_ratio)
ggsave("../chapter/figure/pca/scores_weight.png", width = 3.56, height = 2.6)

##  also study scores in relation to overall ruminoccocus / lachospiraceae ratio
family_means <- processed$mseqtab %>%
  group_by(family, Number) %>%
  summarise(family_mean = mean(value)) %>%
  spread(family, family_mean) %>%
  group_by(Number) %>%
  summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)

scores <- scores %>%
  left_join(family_means)
ggplot(scores) +
  geom_point(
    aes(x = PC1, y = PC2, size = PC3, col = rl_ratio)
  ) +
  labs(
    "x" = perc_label(pc_res, 1),
    "y" = perc_label(pc_res, 2)
  ) +
  scale_color_viridis(
    "Rum. / Lach. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  ) +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(ratio = asp_ratio)
ggsave("../chapter/figure/pca/scores_rl_ratio.png", width = 3.56, height = 2.6)
