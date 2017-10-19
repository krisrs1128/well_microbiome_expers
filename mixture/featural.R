#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Instead of simulataneously identifying scores, fit separate models to each
## data set and then study correlations between cluster memberships across
## tables. This implements this idea for pca on the body composition data and
## latent dirichlet allocation on the counts table.
##
## author: sankaran.kris@gmail.com
## date: 10/18/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("rstan")
library("tidyverse")
library("reshape2")
library("viridis")
source("../dimension_red/prep_tables.R")

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
  "outdir" = "../chapter/figure/featural/"
)
dir.create(opts$outdir)
processed <- process_data(
  raw$seqtab,
  raw$bc,
  raw$bc_full,
  raw$taxa,
  opts
)

###############################################################################
## Fit separate models
###############################################################################
## pca on body composition data
pc_res <- princomp(scale(processed$bc_full))

## lda on counts data
m <- stan_model("lda_counts.stan")

stan_data <- list(
  K = 5,
  V = ncol(processed$x_seq),
  N = nrow(processed$x_seq),
  x = processed$x_seq,
  alpha = rep(1, 5),
  gamma = rep(0.1, ncol(processed$x_seq))
)


lda_fit <- vb(m, stan_data, adapt_engaged = FALSE, eta = 0.1)
lda_res <- rstan::extract(lda_fit)
rm(lda_fit)

###############################################################################
## Cluster scores across tables
###############################################################################
theta_hat <- lda_res$theta %>%
  apply(c(2, 3), mean)

id_map <- setNames(as.character(raw$bc$id), raw$bc$number)
subset_ids <- id_map[rownames(processed$bc)]

cross_scores <- cbind(
  pc_res$scores[subset_ids, 1:5],
  theta_hat
) %>%
  scale()

hc <- hclust(dist(cross_scores))
cross_df <- data.frame(cross_scores) %>%
  rownames_to_column("id") %>%
  melt(
    id.vars = "id",
    variable.name = "axis",
    value.name = "score"
  ) %>%
  mutate(
    id = factor(id, rownames(cross_scores)[hc$order])
  )

ggplot(cross_df) +
  geom_tile(
    aes(x = id, y = axis, fill = score)
  ) +
  scale_fill_gradient2(low = "darkpurple", high = "darkgreen")

pairs(cross_scores[, 1:5], cross_scores[, 5:10])

###############################################################################
## Can we interpret loadings / clusters associated with these scores?
###############################################################################
