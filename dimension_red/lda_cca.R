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
source("lda_cca_plot.R")

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
  "rlog" = FALSE,
  ## "stan_file" = "lda_cca.stan",
  ## "outdir" = "../chapter/figure/lda_cca/"
  "stan_file" = "lda_cca_exp.stan",
  "outdir" = "../chapter/figure/lda_cca_exp/"
)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## Run and plot LDA / CCA on the two (scaled) tables
###############################################################################
bc_mat <- scale(processed$bc)

K <- 2
L1 <- 3
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

m <- stan_model(opts$stan_file)
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
    sprintf("%s/shared_scores_%s.png", opts$outdir, i),
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
  loadings %>% filter(type == "seq"),
  c(1, 1),
  a = 0.9
) +
  facet_wrap(~family, ncol = 5) +
  scale_size_continuous(range = c(0, 2), guide = FALSE) +
  theme(
    axis.title = element_blank(),
    legend.position="none"
  )
plot_topics(loadings %>% filter(type == "seq"))

ggsave(
  sprintf("%s/shared_loadings_seq.png", opts$outdir),
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
    sprintf("%s/unshared_scores_%s.png", opts$outdir, i),
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
ggsave(
  sprintf("%s/loadings_seq_scatter.png", opts$outdir),
  width = 5.56, height = 3.5
)
plot_topics(loadings_x)
ggsave(
  "../chapter/figure/lda_cca/loadings_seq_rows.png",
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
  sprintf("%s/loadings_body_comp.png", opts$outdir),
  width = 6.0, height = 1.5
)

###############################################################################
## Instead of providing point estimate posterior means, we can display the
## entire posteriors for different scores and loadings
###############################################################################
mdist <- melt_parameters(posterior)
mdist$xi_s <- reshape_posterior_score(mdist$xi_s, raw$bc)
mdist$xi_x <- reshape_posterior_score(mdist$xi_x, raw$bc)
mdist$xi_y <- reshape_posterior_score(mdist$xi_y, raw$bc)
lda_cca_plots(mdist, seq_families, processed, opts)
