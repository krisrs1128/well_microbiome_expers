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

plot_loadings(loadings, cca_res$CanCorr) +
  ylim(-0.5, 0.4) +
  xlim(-0.9, 0.3)
ggsave("../chapter/figure/cca/loadings.png", width = 4.56, height = 2.3)

###############################################################################
## Plot the scores
###############################################################################

## extract scores and join in sample data
scores <- rbind(
  data.frame(
    type = "body_comp",
    Number = rownames(bc_mat),
    cca_res$Cx[, 1:3]
  ),
  data.frame(
    type = "seq",
    Number = rownames(x_seq),
    cca_res$Cy[, 1:3]
  )
) %>%
  left_join(
    data.frame(
      Number = rownames(processed$bc),
      processed$bc
    )
  )

mscores <- scores %>%
  select(Number, type, starts_with("CanAxis")) %>%
  gather(comp, value, starts_with("CanAxis")) %>%
  unite(comp_type, comp, type) %>%
  spread(comp_type, value)

ggplot() +
  geom_segment(
    data = mscores,
    aes(
      x = CanAxis1_body_comp, xend = CanAxis1_seq,
      y = CanAxis2_body_comp, yend = CanAxis2_seq,
      size = (CanAxis3_body_comp + CanAxis3_seq) / 2
    ),
    alpha = 0.1
  ) +
  geom_point(
    data = scores,
    aes(
      x = CanAxis1,
      y = CanAxis2,
      size = CanAxis3,
      col = type
    )
  ) +
  labs(
    "col" = "Meas. Type",
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2)
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/scores_linked.png", width = 3.56, height = 2.6)

## same figure but shaded by weight
ggplot() +
  geom_segment(
    data = mscores,
    aes(
      x = CanAxis1_body_comp, xend = CanAxis1_seq,
      y = CanAxis2_body_comp, yend = CanAxis2_seq,
      size = (CanAxis3_body_comp + CanAxis3_seq) / 2
    ),
    alpha = 0.1
  ) +
  geom_point(
    data = scores,
    aes(
      x = CanAxis1,
      y = CanAxis2,
      size = CanAxis3,
      col = weight_dxa
    )
  ) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  ) +
  labs(
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2)
  ) +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/scores_weight.png", width = 3.56, height = 2.6)

##  also study scores in relation to overall ruminoccocus / lachospiraceae ratio
family_means <- processed$mseqtab %>%
  group_by(family, Number) %>%
  summarise(family_mean = mean(value)) %>%
  spread(family, family_mean) %>%
  group_by(Number) %>%
  summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)

scores <- scores %>%
  left_join(family_means)
ggplot() +
  geom_segment(
    data = mscores,
    aes(
      x = CanAxis1_body_comp, xend = CanAxis1_seq,
      y = CanAxis2_body_comp, yend = CanAxis2_seq,
      size = (CanAxis3_body_comp + CanAxis3_seq) / 2
    ),
    alpha = 0.1
  ) +
  geom_point(
    data = scores,
    aes(
      x = CanAxis1,
      y = CanAxis2,
      size = CanAxis3,
      col = rl_ratio
    )
  ) +
  scale_color_viridis(
    "Rum. / Lach. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  ) +
  labs(
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2)
  ) +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/scores_rl_ratio.png", width = 3.56, height = 2.6)
