#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Graph-fused lasso applied to the body composition data.
##
## author: sankaran.kris@gmail.com
## date: 10/16/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("reshape2")
source("../dimension_red/prep_tables.R")
library("gflasso")
library("glasso")

## cleaner ggplot theme
scale_colour_discrete <- function(...)
  scale_color_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#ffc100"),
    na.value = "#464646"
  )
scale_fill_discrete <- function(...)
  scale_fill_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',"#ffc100"),
    na.value = "#464646"
  )

theme_set(theme_bw())
theme_update(
  panel.background = element_rect(fill = "#F8F8F8"),
  panel.border = element_rect(size = 0.5),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 10),
  strip.background = element_blank(),
  strip.text = element_text(size = 10),
  legend.key = element_blank()
)

###############################################################################
## read and prepare the data
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.07)
processed <- process_data(raw$seqtab, raw$bc, raw$bc_full, raw$taxa, opts)

y <- scale(processed$bc)
x <- scale(processed$x_seq)

R <- glasso(cov(y), rho = 1e-1)
opts <- list(
  gamma = 0.01,
  lambda = 0.1,
  eps = 0.01,
  verbose = TRUE
)

fit <- gflasso(y, x, R$w, opts)

###############################################################################
## plot fitted coefficients
###############################################################################
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

site_ordered <- c(
  "aoi", "age", "height_dxa", "weight_dxa",
  "bmi", "android_fm", "android_lm", "gynoid_fm", "gynoid_lm", "l_trunk_fm",
  "l_trunk_lm", "r_trunk_fm", "r_trunk_lm", "trunk_fm", "trunk_lm",
  "l_total_fm", "l_total_lm", "r_total_fm", "r_total_lm", "total_fm",
  "total_lm", "l_leg_fm", "l_leg_lm", "r_leg_fm", "r_leg_lm", "legs_fm",
  "legs_lm", "l_arm_fm", "l_arm_lm", "r_arm_fm", "r_arm_lm", "arms_fm",
  "arms_lm"
)
mass_type_ordered <- c(
  site_ordered[!grepl("fm|lm", site_ordered)],
  site_ordered[grepl("fm", site_ordered)],
  site_ordered[grepl("lm", site_ordered)]
)

mR <- melt(
  R$w,
  varnames = c("x", "y"),
  value.name = "fill"
)
mR$x <- colnames(y)[mR$x]
mR$y <- colnames(y)[mR$y]
mR$x <- factor(mR$x, mass_type_ordered)
mR$y <- factor(mR$y, mass_type_ordered)
ggscaffold::ggheatmap(mR) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0),
    axis.title = element_blank(),
    legend.position = "bottom"
  )
ggsave(
  "../chapter/figure/graph_lasso/graph.png",
  width = 3.17,
  height = 3.07
)

species_order <- colnames(x)[hclust(dist(fit$B))$order]
write.csv(species_order, file = "species_order.csv", row.names = FALSE)

colnames(fit$B) <- colnames(y)
mbeta <- fit$B %>%
  melt(varnames = c("seq_num", "feature")) %>%
  left_join(seq_families) %>%
  mutate(
    feature = factor(feature, mass_type_ordered),
    seq_num = factor(seq_num, species_order)
  )

ggplot(mbeta) +
  geom_tile(
    aes(x = seq_num, y = feature, fill = value)
  ) +
  geom_rect(
    aes(col = family),
    fill = "transparent", size = 2,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  scale_fill_gradient2(
    guide = guide_colorbar(ticks = FALSE, barheight = 0.6),
    mid = "#F8F8F8", low = "#40004b", high = "#00441b"
  ) +
  scale_colour_discrete() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ family, scale = "free", space = "free") +
  theme(
    axis.text = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0, "cm"),
    axis.text.y = element_text(angle = 0, hjust = 0),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/graph_lasso/coef_heatmap.png",
  width = 8,
  height = 3.78
)
