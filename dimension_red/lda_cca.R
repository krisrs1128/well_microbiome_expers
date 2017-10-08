#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An experiment using Latent Dirichlet Allocation + Canonical Correlation
## Analysis simultaneously across count and gaussian tables.
##
## author: sankaran.kris@gmail.com
## date: 10/04/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("rstan")
library("tidyverse")
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
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c', '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', '#6a3d9a'),
    na.value = "black"
  )
scale_fill_discrete <- function(...)
  scale_fill_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c', '#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', '#6a3d9a'),
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
opts <- list("filt_k" = 0.02, "rlog" = FALSE)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
bc_mat <- scale(processed$bc)

K <- 3
L1 <- 4
L2 <- 3

stan_data <- list(
  "n" = nrow(bc_mat),
  "p1" = ncol(processed$x_seq),
  "p2" = ncol(bc_mat),
  "K" = K,
  "L1" = L1,
  "L2" = L2,
  "sigma" = 1,
  "a0" = 1,
  "b0" = 1,
  "x" = processed$x_seq,
  "y" = bc_mat,
  "id_y" = diag(ncol(bc_mat)),
  "id_k" = diag(K),
  "id_l1" = diag(L1),
  "id_l2" = diag(L2),
  "zeros_k" = rep(0, K),
  "zeros_l1" = rep(0, L1),
  "zeros_l2" = rep(0, L2)
)

m <- stan_model("lda_cca.stan")
vb_fit <- vb(m, data = stan_data)

posterior <- rstan::extract(vb_fit)
rm(vb_fit)
pmeans <- parameter_means(posterior)
for (i in seq_along(pmeans)) {
  pmeans[[i]] <- cbind(pmeans[[i]], 1) # in case only used 2 dimensions
}

###############################################################################
## Plot the results
###############################################################################
scv <- scale_color_viridis(
  guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
)
seq_families <- processed$mseqtab %>%
  select(seq_num, family) %>%
  unique()

## Plot some of the shared scores and loadings
p <- plot_scores_wrapper(pmeans$xi_s, raw, processed, scv)
for (i in seq_along(p)) {
  p[[i]] +
    scale_size_continuous(range = c(0, 1.5), guide = FALSE) +
    theme(axis.title = element_blank())
  ggsave(
    sprintf("../chapter/figure/lda_cca/shared_scores_%s.png", i),
    width = 3.56, height = 2.6
  )
}

rownames(pmeans$By) <- colnames(processed$bc)
rownames(pmeans$Bx) <- colnames(processed$x_seq)
loadings <- prepare_loadings(
  list(pmeans$By, pmeans$Bx),
  c("body_comp", "seq")
) %>%
  left_join(seq_families)

plot_loadings(
  loadings %>% filter(type == "seq", !is.na(family)),
  c(1, 1),
  a = 0.9
) +
  facet_wrap(~family, ncol = 5) +
  scale_size_continuous(range = c(0, 2), guide = FALSE) +
  theme(
    axis.title = element_blank(),
    legend.position="none"
  )

ggsave(
  sprintf("../chapter/figure/lda_cca/shared_loadings_seq.png"),
  width = 6.56, height = 3.9
  )
plot_loadings(loadings %>% filter(type == "body_comp"), c(1, 1)) +
  scale_size_continuous(range = c(1.3, 2), guide = FALSE) +
  theme(axis.title = element_blank())
ggsave(
  "../chapter/figure/lda_cca/shared_loadings_body_comp.png",
  width = 4.0, height = 4.0
)

## now plot unshared scores (species abundances first)
p <- plot_scores_wrapper(pmeans$xi_x, raw, processed, scv)
p <- c(p, plot_scores_wrapper(pmeans$xi_y, raw, processed, scv))

for (i in seq_along(p)) {
  p[[i]] + scale_size_continuous(range = c(0, 1.5), guide = FALSE)
  ggsave(
    sprintf("../chapter/figure/lda_cca/unshared_scores_%s.png", i),
    width = 3.56, height = 2.6
  )
}

## Now plot unshared loadings
rownames(pmeans$Wx) <- colnames(processed$x_seq)
loadings_x <- prepare_loadings(list(pmeans$Wx), "seq") %>%
  left_join(seq_families) %>%
  filter(!is.na(family))
plot_loadings(loadings_x, c(1, 1), a = 0.9) +
  facet_wrap(~family, ncol = 5) +
  scale_size_continuous(range = c(0, 2), guide = FALSE) +
  theme(
    axis.title = element_blank(),
    legend.position = "none"
  )
plot_topics(loadings_x)

ggsave(
  "../chapter/figure/lda_cca/loadings_seq.png",
  width = 5.56, height = 3.5
)

rownames(pmeans$Wy) <- colnames(processed$bc)
loadings_y <- prepare_loadings(list(pmeans$Wy), "body_comp") %>%
  mutate(seq_num = "NA") %>%
  left_join(seq_families)
plot_loadings(loadings_y, c(1, 1)) +
  scale_size_continuous(range = c(1.5, 3), guide = FALSE) +
  theme(axis.title = element_blank())
ggsave(
  "../chapter/figure/lda_cca/loadings_body_comp.png",
  width = 6.0, height = 1.5
)

## Plot posteriors associated with scores / topics
mdist <- melt_parameters(posterior)
mdist$xi_s <-reshape_posterior_score(mdist$xi_s, bc)
mdist$xi_x <-reshape_posterior_score(mdist$xi_x, bc)
mdist$xi_y <-reshape_posterior_score(mdist$xi_y, bc)

ggplot() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point(
    data = mdist$xi_s %>%
      mutate(
        Total_FM_cut = cut(round(Total_FM / 1e4, 2), 8)
      ),
    aes(
      x = topic_1,
      y = topic_2,
      col = topic_3
    ),
    size = 0.1,
    alpha = 0.01
  ) +
  coord_equal() +
  facet_wrap(~Total_FM_cut, ncol = 4) +
  scale_color_viridis(
    option = "magma",
    guide = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  labs(x = "Axis 1", y = "Axis 2", col = "Axis 3")

ggplot() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point(
    data = mdist$xi_s %>%
      mutate(
        weight_dxa_cut = cut(weight_dxa, 8)
      ),
    aes(
      x = topic_1,
      y = topic_2,
      col = topic_3
    ),
    size = 0.1,
    alpha = 0.01
  ) +
  coord_equal() +
  facet_wrap(~weight_dxa_cut, ncol = 4) +
  scale_color_viridis(
    option = "magma",
    guide = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  labs(x = "Axis 1", y = "Axis 2", col = "Axis 3")

## using the scores from just the bc table
ggplot() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point(
    data = mdist$xi_y %>%
      mutate(
        Total_FM_cut = cut(round(Total_FM / 1e4, 2), 4)
      ),
    aes(
      x = topic_1,
      y = topic_2,
      col = topic_3
    ),
    size = 0.1,
    alpha = 0.01
  ) +
  coord_equal() +
  facet_wrap(~Total_FM_cut, nrow = 1) +
  scale_color_viridis(
    option = "magma",
    guide = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  labs(x = "Axis 1", y = "Axis 2", col = "Axis 3")

## Now plotting loadings boxplots
mdist$Wx$seq_num <- colnames(processed$x_seq)[mdist$Wx$row]
mdist$Wx <- mdist$Wx %>%
  left_join(seq_families)
mdist$Wx$seq_num <- factor(
  mdist$Wx$seq_num,
  levels = colnames(processed$x_seq)[taxa_order(loadings)]
)

ggplot(mdist$Wx) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_boxplot(
    aes(x = seq_num, y = value, col = family, fill = family),
    outlier.size = 0.05,
    width = 0.1
  ) +
  facet_grid(col ~ family, scale = "free_x", space = "free_x") +
  theme(
    panel.spacing.x = unit(0, "cm"),
    axis.text.x = element_blank(),
    strip.text.x = element_blank()
  )
ggsave(
  "../chapter/figure/lda_cca/within_loadings_boxplots.png",
  width = 8, height = 3
)

mdist$Bx$seq_num <- colnames(processed$x_seq)[mdist$Bx$row]
mdist$Bx <- mdist$Bx %>%
  left_join(seq_families)
mdist$Bx$seq_num <- factor(
  mdist$Bx$seq_num,
  levels = colnames(processed$x_seq)[taxa_order(loadings)]
)

ggplot(mdist$Bx) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_boxplot(
    aes(x = seq_num, y = value, col = family, fill = family),
    outlier.size = 0.05,
    width = 0.1
  ) +
  facet_grid(col ~ family, scale = "free_x", space = "free_x") +
  theme(
    panel.spacing.x = unit(0, "cm"),
    axis.text.x = element_blank(),
    strip.text.x = element_blank()
  )
ggsave(
  "../chapter/figure/lda_cca/between_loadings_boxplots.png",
  width = 8, height = 3
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
