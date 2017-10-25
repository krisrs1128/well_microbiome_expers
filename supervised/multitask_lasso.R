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
opts <- list("filt_k" = 0.07)
processed <- process_data(raw$seqtab, raw$bc, raw$bc_full, raw$taxa, opts)
y <- scale(processed$bc)
x <- scale(processed$x_seq)

n_lambda <- 30
lambdas <- seq(0.001, 0.7, length.out = n_lambda)
beta_hats <- array(dim = c(ncol(y), 1 + ncol(x), n_lambda))
fits <- list()

for (r in seq_len(ncol(y))) {
  message("tuning response ", r)
  fit[[r]] <- cv.glmnet(x, y[, r], lambda = lambdas, alpha = 0.3)
  beta_hats[r,, ] <- as.matrix(coef(fit[[r]]$glmnet.fit))
}

###############################################################################
## Plot the cross validation errors
###############################################################################
cv_err <- do.call(cbind, sapply(fits, function(x) x$cvm))
cv_err[cv_err > 1.4] <- NA
colnames(cv_err) <- colnames(y)
rownames(cv_err) <- lambdas

ggplot(melt(cv_err)) +
  geom_tile(
    aes(x = Var1, y = Var2, fill = value)
  ) +
  scale_fill_gradient2(midpoint = 1, low = "blue", high = "red") +
  scale_x_continuous(expand = c(0, 0))
  scale_y_discrete(expand = c(0, 0))

###############################################################################
## Plot the results
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
  site_ordered[!grepl("FM|LM", site_ordered)],
  site_ordered[grepl("FM", site_ordered)],
  site_ordered[grepl("LM", site_ordered)]
)

rownames(beta_hats) <- colnames(y)
colnames(beta_hats) <- c("intercept", colnames(x))
mbeta <- beta_hats %>%
  melt(
    varnames = c("feature", "seq_num", "lambda")
  ) %>%
  left_join(seq_families) %>%
  mutate(feature = factor(feature, mass_type_ordered))

ggplot(mbeta) +
  geom_tile(
    aes(x = seq_num, y = lambda, fill = value)
  ) +
  geom_rect(
    aes(col = family),
    fill = "transparent", size = 1,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
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
  "../chapter/figure/graph_lasso/multitask_lasso_hm.png",
  width = 6.86,
  height = 3.78
)
