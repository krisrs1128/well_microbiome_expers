#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Sparse partial least squares applied to the WELL microbiome data.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("reshape2")
source("../dimension_red/prep_tables.R")
library("spls")

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
opts <- list("filt_k" = 0.02, "filt_a" = 0)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

y <- scale(processed$bc)
x <- scale(processed$x_seq)
cv_eval <- cv.spls(x, y, K = 1:6, eta = seq(0, 0.9, 0.05), scale.x = FALSE, fold = 5)
cv_eval

train_ix <- sample(1:nrow(x), 80)
#fit <- spls(x[train_ix, ], y[train_ix, ], cv_eval$K.opt, cv_eval$eta.opt)
fit <- spls(x[train_ix, ], y[train_ix, ], 4, 0.6)
y_hat <- x %*% fit$betahat
plot(y[train_ix, 24], y_hat[train_ix, 24])
points(y[-train_ix, 24], y_hat[-train_ix, 24], col = "blue")
abline(a = 0, b = 1, col = "red")

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
species_order <- colnames(x)[hclust(dist(fit$betahat))$order]

mbeta <- fit$betahat %>%
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
    guide = guide_colorbar(ticks = FALSE, keyheight = 0.5),
    low = "#40004b",
    high = "#00441b"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ family, scale = "free", space = "free") +
  theme(
    axis.text = element_blank(),
    panel.spacing = unit(0, "cm"),
    axis.text.y = element_text(size = 7, angle = 0, hjust = 0),
    strip.text.x = element_text(size = 7, angle = 90, hjust = 0),
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/spls/coef_heatmap.png",
  width = 8.77,
  height = 5.17
)

large_species <- mbeta %>%
  filter(
    feature %in% c("Total_LM", "Total_FM")
  ) %>%
  group_by(seq_num) %>%
  mutate(norm = sqrt(sum(value ^ 2))) %>%
  filter(norm > 0.1) %>%
  arrange(desc(norm))

mlarge_species <- melt(
  data.frame(
    "Number" = rownames(x),
    x[, unique(large_species$seq_num)],
    y[, c("Total_FM", "Total_LM")]
  ),
  id.vars = c("Number", "Total_FM", "Total_LM"),
  variable.name = "seq_num"
) %>%
  left_join(seq_families)

mlarge_species$seq_num <- factor(
  mlarge_species$seq_num,
  mbeta %>%
    filter(feature == "Total_LM") %>%
    arrange(value) %>%
    .[["seq_num"]]
)

ggplot(mlarge_species) +
  geom_hline(yintercept = 0, size = 0.1, alpha = 0.8) +
  geom_vline(xintercept = 0, size = 0.1, alpha = 0.8) +
  geom_point(
    aes(x = value, y = Total_LM, col = family),
    size = 0.7, alpha = 0.8
  ) +
  facet_wrap(~seq_num, ncol = 8) +
  theme(
    legend.position = "bottom"
  )
ggsave(
  "../chapter/figure/spls/total_lm_species.png",
  width = 7.4,
  height = 4.3
)

## same plot for Total FM
mlarge_species$seq_num <- factor(
  mlarge_species$seq_num,
  mbeta %>%
    filter(feature == "Total_FM") %>%
    arrange(value) %>%
    .[["seq_num"]]
)

ggplot(mlarge_species) +
  geom_hline(yintercept = 0, size = 0.1, alpha = 0.8) +
  geom_vline(xintercept = 0, size = 0.1, alpha = 0.8) +
  geom_point(
    aes(x = value, y = Total_FM, col = family),
    size = 0.7, alpha = 0.8
  ) +
  facet_wrap(~seq_num, ncol = 8) +
  theme(
    legend.position = "bottom"
  )
ggsave(
  "../chapter/figure/spls/total_fm_species.png",
  width = 7.4,
  height = 4.3
)
