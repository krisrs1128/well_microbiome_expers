#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Example using Co-Inertia analysis.
##
## author: sankaran.kris@gmail.com
## date: 09/27/2017

library("phyloseq")
library("DESeq2")
library("tidyverse")
library("forcats")
library("reshape2")
library("ggrepel")
library("viridis")
library("ade4")

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

out_path <- "../chapter/figure/coia"
dir.create(out_path)

###############################################################################
## Load data
###############################################################################
raw <- read_data()
opts <- list(filt_k = 0.5, filt_a = 0)
processed <- process_data(raw$seqtab, raw$bc, raw$taxa, opts)

###############################################################################
## Run CoIA on the two (scaled) tables
###############################################################################
dudi1 <- dudi.pca(x_seq, scan = FALSE, nf = 3)
dudi2 <- dudi.pca(bc_mat, scan = FALSE, nf = 3)
coia_res <- coinertia(dudi1, dudi2, scan = FALSE, nf = 3)

loadings <- prepare_loadings(
  list(coia_res$l1, coia_res$c1),
  c("body_comp", "seq")
) %>%
  left_join(processed$mseqtab)
plot_loadings(loadings, coia_res$eig) +
  ylim(-0.6, 0.5) +
  xlim(-0.45, 0.35) +
  scale_size_continuous(range = c(0, 3), breaks = c(-8, -8))
ggsave(file.path(out_path, "loadings.png"), width = 4.56, height = 2.3)

scores <- prepare_scores(
  list(coia_res$lX, coia_res$lY),
  c("body_comp", "seq")
) %>%
  left_join(
    processed$bc %>%
    rownames_to_column("Number")
  )

mscores <- melt_scores(scores)
plot_scores(scores, "type", "Meas. Type", coia_res$eig) +
  link_scores(mscores) +
  scale_color_brewer(palette = "Set1")
ggsave(file.path(out_path, "scores_linked.png"), width = 4.7, height = 1.7)

plot_scores(scores, "Total_FM", "Trunk FM", coia_res$eig) +
  link_scores(mscores) +
  scale_color_viridis(
    "Total FM ",
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
ggsave(file.path(out_path, "scores_total_fm.png"), width = 4.7, height = 1.7)

scores <- scores %>%
  left_join(family_means(processed$mseqtab))
plot_scores(scores, "rl_ratio", "Rum. / Lach ratio", coia_res$eig) +
  link_scores(mscores) +
  scale_color_viridis(
    "Rum. / Lach. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  )
ggsave(file.path(out_path, "scores_rl_ratio.png"), width = 4.7, height = 1.7)
