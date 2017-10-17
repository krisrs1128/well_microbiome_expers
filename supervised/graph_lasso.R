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
## read and prepare the data
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.07)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

y <- scale(processed$bc)
x <- scale(processed$x_seq)

R <- cor(y)
opts <- list(
  gamma = 0.01,
  lambda = 0.1,
  eps = 0.01,
  verbose = TRUE
)

fit <- gflasso(y, x, R, opts)

###############################################################################
## plot fitted coefficients
###############################################################################
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

site_ordered <- c(
  "aoi", "age", "height_dxa", "weight_dxa",
  "bmi", "Android_FM", "Android_LM", "Gynoid_FM", "Gynoid_LM", "L_Trunk_FM",
  "L_Trunk_LM", "R_Trunk_FM", "R_Trunk_LM", "Trunk_FM", "Trunk_LM",
  "L_Total_FM", "L_Total_LM", "R_Total_FM", "R_Total_LM", "Total_FM",
  "Total_LM", "L_Leg_FM", "L_Leg_LM", "R_Leg_FM", "R_Leg_LM", "Legs_FM",
  "Legs_LM", "L_Arm_FM", "L_Arm_LM", "R_Arm_FM", "R_Arm_LM", "Arms_FM",
  "Arms_LM"
)
mass_type_ordered <- c(
  site_ordered[!grepl("FM|LM", site_ordered)],
  site_ordered[grepl("FM", site_ordered)],
  site_ordered[grepl("LM", site_ordered)]
)

species_order <- colnames(x)[hclust(dist(t(x)))$order]

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
  scale_fill_gradient2(
    guide = guide_colorbar(ticks = FALSE, barheight = 0.6),
    low = "#40004b",
    high = "#00441b"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ family, scale = "free", space = "free") +
  theme(
    axis.text = element_blank(),
    panel.spacing = unit(0, "cm"),
    axis.text.y = element_text(size = 5, angle = 0, hjust = 0),
    strip.text.x = element_text(size = 7, angle = 90, hjust = 0),
    legend.position = "bottom"
  )
