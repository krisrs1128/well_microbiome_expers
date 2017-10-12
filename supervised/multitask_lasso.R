#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using lasso to model multiple responses in the WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("viridis")
library("glmnet")
source("../dimension_red/prep_tables.R")

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
## Apply parallel lassos
###############################################################################
raw <- read_data()
opts <- list("filt_k" = 0.02, "rlog" = FALSE)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
y <- scale(processed$bc)
x <- scale(processed$x_seq)

n_lambda <- 50
beta_hats <- array(dim = c(ncol(y), 1 + ncol(x), n_lambda))
for (r in seq_len(ncol(y))) {
  fit <- glmnet(x, y[, r], nlambda = n_lambda, alpha = 0.7)
  beta_hats[r,, ] <- as.matrix(coef(fit))
}
plot(y[, r], predict(fit, x)[, 20]) ## seems like some real associations (though not test set)

###############################################################################
## Plot the results
###############################################################################
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

rownames(beta_hats) <- colnames(y)
colnames(beta_hats) <- c("intercept", colnames(x))
mbeta <- beta_hats %>%
  melt(
    varnames = c("feature", "seq_num", "lambda")
  ) %>%
  left_join(seq_families)

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

mbeta$feature <- factor(
  mbeta$feature,
  levels = mass_type_ordered
)

ggplot(mbeta) +
  geom_tile(
    aes(x = seq_num, y = lambda, fill = value)
  ) +
  scale_fill_gradient2(
    guide = guide_colorbar(ticks = FALSE, keyheight = 0.5),
    low = "#40004b",
    high = "#00441b"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(feature ~ family, scale = "free", space = "free") +
  theme(
    axis.text = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
    strip.text.x = element_text(size = 7, angle = 90, hjust = 0),
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/multitask_lasso_hm.png",
  height = 11.3,
  width = 8.01
)
