#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Variant of lda_cca.R that runs the shared LDA scores model only only CCA.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

library("phyloseq")
library("rstan")
library("tidyverse")
library("PMA")
library("reshape2")
library("ggrepel")
library("viridis")
source("prep_tables.R")
source("plot.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
opts <- list(
  "filt_k" = 0.02,
  "stan_file" = "lda_cca.stan",
  "outdir" = "../chapter/figure/lda_cca_selected/"
)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
bc_mat <- scale(processed$bc)

###############################################################################
## Identify subset of features to run on using CCA
###############################################################################
K <- 2
cca_res <- CCA(
  bc_mat,
  scale(processed$x_seq),
  penaltyx = 0.6,
  penaltyz = 0.3,
  K = K
)

opts$rlog <- FALSE
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)
keep_ix <- rowSums(cca_res$v) != 0

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
L1 <- 3
L2 <- 3

stan_data <- list(
  "n" = nrow(bc_mat),
  "p1" = sum(keep_ix),
  "p2" = ncol(bc_mat),
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "sigma" = 1,
  "a0" = 1,
  "b0" = 1,
  "x" = processed$x_seq[, keep_ix],
  "y" = bc_mat,
  "id_y" = diag(ncol(bc_mat)),
  "id_k" = diag(K),
  "id_l1" = diag(L1),
  "id_l2" = diag(L2),
  "zeros_k" = rep(0, K),
  "zeros_l1" = rep(0, L1),
  "zeros_l2" = rep(0, L2)
)

m <- stan_model(opts$stan_file)
vb_fit <- vb(m, data = stan_data)
posterior <- rstan::extract(vb_fit)

###############################################################################
## Plot the results
###############################################################################
scv <- scale_color_viridis(
  guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
)
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

mdist <- melt_parameters(posterior)
mdist$xi_s <- reshape_posterior_score(mdist$xi_s, raw$bc)
mdist$xi_x <- reshape_posterior_score(mdist$xi_x, raw$bc)
mdist$xi_y <- reshape_posterior_score(mdist$xi_y, raw$bc)

ggplot() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point(
    data = mdist$xi_s,
    aes(
      x = axis_1,
      y = axis_2,
      col = Total_LM
    ),
    size = 0.1,
    alpha = 0.01
  ) +
  coord_equal() +
  scv +
  labs(x = "Axis 1", y = "Axis 2", col = "Total LM")

ggsave(
  sprintf("%s/shared_scores_lm_posterior.png", opts$outdir),
  width = 5.63, height = 3.42
)

## using the scores from just the bc table
ggplot() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point(
    data = mdist$xi_y,
    aes(
      x = axis_1,
      y = axis_2,
      col = Total_LM
    ),
    size = 0.1,
    alpha = 0.1
  ) +
  coord_equal() +
  scv +
  labs(x = "Axis 1", y = "Axis 2", col = "Total LM")
ggsave(
  sprintf("%s/unshared_scores_lm_posterior.png", opts$outdir),
  width = 9.24,
  height = 2.5
)

## Now plotting loadings boxplots
mdist$Wx$seq_num <- colnames(processed$x_seq)[mdist$Wx$row]
mdist$Wx <- mdist$Wx %>%
  left_join(seq_families)
mdist$Wx$seq_num <- factor(
  mdist$Wx$seq_num,
  levels = names(sort(colSums(processed$x_seq), decreasing = TRUE))
)

wx_summary <- mdist$Wx %>%
  group_by(col, seq_num) %>%
  dplyr::summarise(
    family = family[1],
    lower = quantile(value, 0.25),
    upper = quantile(value, 0.75),
    med = median(value)
  )

ggplot(wx_summary) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_pointrange(
    aes(x = seq_num, ymin = lower, ymax = upper, y = med, col = family, fill = family),
    fatten = 1.2,
    size = 0.1
  ) +
  facet_grid(col ~ family, scale = "free_x", space = "free_x") +
  theme(
    panel.spacing.x = unit(0, "cm"),
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave(
  sprintf("%s/within_loadings_seq_boxplots.png", opts$outdir),
  width = 7.45,
  height = 3.37
)

mdist$Bx$seq_num <- colnames(processed$x_seq)[mdist$Bx$row]
mdist$Bx <- mdist$Bx %>%
  left_join(seq_families)
mdist$Bx$seq_num <- factor(
  mdist$Bx$seq_num,
  levels = names(sort(colSums(processed$x_seq), decreasing = TRUE))
)

bx_summary <- mdist$Bx %>%
  group_by(col, seq_num) %>%
  dplyr::summarise(
    family = family[1],
    lower = quantile(value, 0.25),
    upper = quantile(value, 0.75),
    med = median(value)
  )

ggplot(bx_summary) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_pointrange(
    aes(x = seq_num, ymin = lower, ymax = upper, y = med, col = family, fill = family),
    fatten = 1.2,
    size = 0.1
  ) +
  facet_grid(col ~ family, scale = "free_x", space = "free_x") +
  theme(
    panel.spacing.x = unit(0, "cm"),
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave(
  sprintf("%s/between_loadings_seq_boxplots.png", opts$outdir),
  width = 7.45,
  height = 3.37
)

## and finally loadings boxplot for body composition variables
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

mdist$Wy$variable <- tolower(colnames(processed$bc)[mdist$Wy$row])
mdist$Wy$variable <- factor(
  mdist$Wy$variable,
  levels = mass_type_ordered
)

ggplot(mdist$Wy) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_boxplot(
    aes(
      x = variable,
      y = value
    ),
    outlier.size = 1
  ) +
  facet_grid(col ~ .) +
  theme(axis.text.x = element_text(angle = -90))

ggsave(
  sprintf("%s/within_loadings_body_comp_boxplots.png", opts$outdir),
  width = 6.44,
  height = 3.95
)

mdist$By$variable <- tolower(colnames(processed$bc)[mdist$By$row])
mdist$By$variable <- factor(
  mdist$By$variable,
  levels = mass_type_ordered
)

ggplot(mdist$By) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_boxplot(
    aes(
      x = variable,
      y = value
    ),
    outlier.size = 1
  ) +
  facet_grid(col ~ .) +
  theme(axis.text.x = element_text(angle = -90))

ggsave(
  sprintf("%s/between_loadings_body_comp_boxplots.png", opts$outdir),
  width = 6.44,
  height = 3.95
)
