#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Application of the penalized matrix decomposition to the body composition
## data.
##
## author: sankaran.kris@gmail.com
## date: 10/06/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("PMA")
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
## read and prepare the data
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.02)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
x <- scale(processed$bc)
y <- scale(processed$x_seq)

cca_res <- CCA(x, y, penaltyx = 0.6, penaltyz = 0.3, K = 3)

## Plot the loadings
rownames(cca_res$u) <- colnames(x)
rownames(cca_res$v) <- colnames(y)
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

loadings <- prepare_loadings(
  list(data.frame(cca_res$u), data.frame(cca_res$v)),
  c("body_comp", "seq")
) %>%
  left_join(seq_families)

p <- plot_loadings(
  loadings %>%
  filter((Axis.1 != 0 & Axis.2 != 0) | type == "body_comp"),
  cca_res$d
) +
  geom_point(
    data = loadings %>%
      filter(!((Axis.1 != 0 & Axis.2 != 0) | type == "body_comp")),
    aes(x = Axis.1, y = Axis.2, col = family),
    alpha = 0.8, size = 0.9
  ) +
  scale_size(range = c(1.5, 4), breaks = c(-5, 5))

p$layers[c(4, 5)] <- p$layers[c(5, 4)]
p

## Plot the scores
scores <- prepare_scores(
  list(x %*% cca_res$u, y %*% cca_res$v),
  c("body_comp", "seq")
) %>%
  left_join(raw$bc)

mscores <- melt_scores(scores)
plot_scores(scores, "type", "Meas. Type", cca_res$d) +
  link_scores(mscores) +
  scale_color_brewer(palette = "Set1")

## color by weight
plot_scores(scores, "Total_FM", "Total FM", cca_res$d) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )

## color by ruminoccocus / lachnospiraceae ratios
scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "Rum. / Lach. ratio", cca_res$d) +
  link_scores(mscores) +
  scale_color_viridis(
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  )
