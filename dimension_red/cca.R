#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Illustration of canonical correlation analysis on the body composition and
## microbiome elements of the WELL-China data.
##
## author: sankaran.kris@gmail.com
## date: 09/26/2017

library("phyloseq")
library("DESeq2")
library("tidyverse")
library("forcats")
library("reshape2")
library("ggrepel")
library("viridis")
library("vegan")

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

cca_perc <- function(cca_res, i) {
  round(100 * cca_res$CanCorr[i] / sum(cca_res$CanCorr), 2)
}

perc_label <- function(cca_res, i) {
  sprintf("CC%s [%s%%]", i, cca_perc(cca_res, i))
}

###############################################################################
## Load data
###############################################################################
opts <- list(
  gender = "Male",
  sf_quantile = 0.95,
  filt_k = 0.5,
  filt_a = 0
)

seqtab <- readRDS("../data/seqtab.rds")
bc <- readRDS("../data/sample_data_bc.rds")
colnames(seqtab) <- paste0("species_", seq_len(ntaxa(seqtab)))

taxa <- readRDS("../data/taxa.rds") %>%
  data.frame() %>%
  rownames_to_column("seq")
taxa$seq_num <- colnames(seqtab)
taxa$family <- fct_lump(taxa$Family, n = 7)

mseqtab <- seqtab %>%
  melt(varnames = c("Number", "seq_num")) %>%
  left_join(bc) %>%
  left_join(taxa)

###############################################################################
## normalize with DESeq2's varianceStabilizingTransformation
###############################################################################
seqtab <- seqtab %>%
  filter_taxa(function(x) mean(x > opts$filt_a) > opts$filt_k, prune = TRUE)
dds <- DESeqDataSetFromMatrix(
  countData = t(get_taxa(seqtab)),
  colData = data.frame("unused" = rep(1, nrow(seqtab))),
  design = ~1
)

## are quantiles for one sample systematically larger than those for others (if
## so, give it a large size factor). Basically related to sequencing depth.
##
## code copied from here: https://support.bioconductor.org/p/76548/
qs <- apply(counts(dds), 2, quantile, opts$sf_quantile)
sizeFactors(dds) <- qs / exp(mean(log(qs)))

## and the regularized data
vsd <- varianceStabilizingTransformation(dds, fitType = "local")
x_seq <- t(assay(vsd))
x_seq <- x_seq[bc$Number, ]
x_seq <- scale(x_seq[bc$gender == opts$gender, ])

###############################################################################
## Run CCA on the two (scaled) tables
###############################################################################
bc_mat <- data.frame(bc) %>%
  filter(gender == opts$gender) %>%
  select(-id, -gender)
rownames(bc_mat) <- bc_mat$Number
bc_mat$Number <- NULL
bc_mat <- scale(bc_mat)

cca_res <- CCorA(bc_mat, x_seq)

###############################################################################
## Plot the loadings
###############################################################################
loadings <- rbind(
  data.frame(
    variable = rownames(cca_res$corr.Y.Cy),
    seq_num = NA,
    type = "body_comp",
    cca_res$corr.Y.Cy[, 1:3]
  ),
  data.frame(
    variable = rownames(cca_res$corr.X.Cx),
    seq_num = rownames(cca_res$corr.X.Cx),
    type = "seq",
    cca_res$corr.X.Cx[, 1:3]
  )
) %>%
  left_join(taxa)

asp_ratio <- sqrt(cca_res$CanCorr[2] / cca_res$CanCorr[1])
ggplot(loadings) +
  geom_hline(yintercept = 0, size = 0.5) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    data = loadings %>%
      filter(type == "seq"),
    aes(x = CanAxis1, y = CanAxis2, size = CanAxis3, col = family),
    alpha = 1
  ) +
  geom_text_repel(
    data = loadings %>%
      filter(type == "body_comp"),
    aes(x = CanAxis1, y = CanAxis2, label = variable, size = CanAxis3),
    segment.size = 0.3,
    segment.alpha = 0.5
  ) +
  labs(
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2),
    "col" = "Family"
  ) +
  scale_size_continuous(
    range = c(0, 4),
    breaks = c(-5, 5)
  ) +
  ylim(-0.5, 0.4) +
  xlim(-0.9, 0.3) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/loadings.png", width = 4.56, height = 2.3)

###############################################################################
## Plot the scores
###############################################################################

## extract scores and join in sample data
scores <- rbind(
  data.frame(
    type = "body_comp",
    Number = rownames(bc_mat),
    cca_res$Cx[, 1:3]
  ),
  data.frame(
    type = "seq",
    Number = rownames(x_seq),
    cca_res$Cy[, 1:3]
  )
) %>%
  left_join(bc)

mscores <- scores %>%
  select(Number, type, starts_with("CanAxis")) %>%
  gather(comp, value, starts_with("CanAxis")) %>%
  unite(comp_type, comp, type) %>%
  spread(comp_type, value)

ggplot() +
  geom_segment(
    data = mscores,
    aes(
      x = CanAxis1_body_comp, xend = CanAxis1_seq,
      y = CanAxis2_body_comp, yend = CanAxis2_seq,
      size = (CanAxis3_body_comp + CanAxis3_seq) / 2
    ),
    alpha = 0.1
  ) +
  geom_point(
    data = scores,
    aes(
      x = CanAxis1,
      y = CanAxis2,
      size = CanAxis3,
      col = type
    )
  ) +
  labs(
    "col" = "Meas. Type",
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2)
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/scores_linked.png", width = 3.56, height = 2.6)

## same figure but shaded by weight
ggplot() +
  geom_segment(
    data = mscores,
    aes(
      x = CanAxis1_body_comp, xend = CanAxis1_seq,
      y = CanAxis2_body_comp, yend = CanAxis2_seq,
      size = (CanAxis3_body_comp + CanAxis3_seq) / 2
    ),
    alpha = 0.1
  ) +
  geom_point(
    data = scores,
    aes(
      x = CanAxis1,
      y = CanAxis2,
      size = CanAxis3,
      col = weight_dxa
    )
  ) +
  scale_color_viridis(
    "Weight ",
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  ) +
  labs(
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2)
  ) +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/scores_weight.png", width = 3.56, height = 2.6)

##  also study scores in relation to overall ruminoccocus / lachospiraceae ratio
family_means <- mseqtab %>%
  group_by(family, Number) %>%
  summarise(family_mean = mean(value)) %>%
  spread(family, family_mean) %>%
  group_by(Number) %>%
  summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)

scores <- scores %>%
  left_join(family_means)
ggplot() +
  geom_segment(
    data = mscores,
    aes(
      x = CanAxis1_body_comp, xend = CanAxis1_seq,
      y = CanAxis2_body_comp, yend = CanAxis2_seq,
      size = (CanAxis3_body_comp + CanAxis3_seq) / 2
    ),
    alpha = 0.1
  ) +
  geom_point(
    data = scores,
    aes(
      x = CanAxis1,
      y = CanAxis2,
      size = CanAxis3,
      col = rl_ratio
    )
  ) +
  scale_color_viridis(
    "Rum. / Lach. Ratio  ",
    guide = guide_colorbar(barwidth= 0.15, ticks = FALSE)
  ) +
  labs(
    "x" = perc_label(cca_res, 1),
    "y" = perc_label(cca_res, 2)
  ) +
  scale_size_continuous(range = c(0, 1.5), breaks = c(-8, 8)) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/cca/scores_rl_ratio.png", width = 3.56, height = 2.6)
