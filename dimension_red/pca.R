#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An example studying multiple tables by concatenating and then using PCA.
##
## author: sankaran.kris@gmail.com
## date: 09/25/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("phyloseq")
library("forcats")
library("DESeq2")
library("ggrepel")
library("reshape2")
library("viridis")

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

pc_perc <- function(pc_res, i) {
  round(100 * pc_res$sdev[i] / sum(pc_res$sdev), 2)
}

perc_label <- function(pc_res, i) {
  sprintf("PC%s [%s%%]", i, pc_perc(pc_res, i))
}

###############################################################################
## Load data
###############################################################################
seqtab <- readRDS("../data/seqtab.rds")
bc <- readRDS("../data/sample_data_bc.rds")
colnames(seqtab) <- paste0("species_", seq_len(ntaxa(seqtab)))

taxa <- readRDS("../data/taxa.rds") %>%
  data.frame() %>%
  rownames_to_column("seq")
taxa$seq_num <- colnames(seqtab)

###############################################################################
## normalize with DESeq2's varianceStabilizingTransformation
###############################################################################
seqtab <- seqtab[, 1:700]
dds <- DESeqDataSetFromMatrix(
  countData = t(get_taxa(seqtab)),
  colData = data.frame("unused" = rep(1, nrow(seqtab))),
  design = ~1
)

## are quantiles for one sample systematically larger than those for others (if
## so, give it a large size factor). Basically related to sequencing depth.
##
## code copied from here: https://support.bioconductor.org/p/76548/
qs <- apply(counts(dds), 2, quantile, .95)
sizeFactors(dds) <- qs / exp(mean(log(qs)))

## and the regularized data
vsd <- varianceStabilizingTransformation(dds, fitType = "local")
x_seq <- assay(vsd)

###############################################################################
## concatenate quantitative body composition (split by gender) and (somewhat
## filtered) microbiome counts, and transformed
###############################################################################
bc[, -c(1:3)] <- scale(bc[, -c(1:3)])
combined_df <- data.frame(bc, scale(t(x_seq[, bc$Number])))
combined <- combined_df %>%
  select(-Number, -id) %>%
  split(.$gender, drop = TRUE) %>%
  lapply(function(x) as.matrix(x[, colnames(x) != "gender"]))

###############################################################################
## combining into single melted data set (mainly for plotting)
###############################################################################
taxa$family <- fct_lump(taxa$Family, n = 7)
mseqtab <- seqtab %>%
  melt(varnames = c("Number", "seq_num")) %>%
  left_join(bc) %>%
  left_join(taxa)

###############################################################################
## run and visualize PCA
###############################################################################
pc_res <- lapply(combined, prcomp)

## extract scores and join in sample data
scores <- data.frame(
  Number = rownames(combined$Male),
  pc_res$Male$x
) %>%
  left_join(bc)

## extract loadings and join taxa information
loadings <- data.frame(
  "variable" = colnames(combined[[1]]),
  pc_res$Male$rotation
  ) %>%
  mutate(
    type = ifelse(variable %in% taxa_names(seqtab), "seq", "body_comp"),
    seq_num = variable
  ) %>%
  left_join(taxa)
loadings[loadings$type != "seq", "seq_num"] <- NA

## how to species data relate to body composition?
asp_ratio <- sqrt(pc_res$Male$sdev[2] / pc_res$Male$sdev[1])
ggplot(loadings) +
  geom_point(
    data = loadings %>%
      filter(type == "seq"),
    aes(x = PC1, y = PC2, col = family),
    size = 1, alpha = 0.6
  ) +
  geom_text_repel(
    data = loadings %>%
      filter(type == "body_comp"),
    aes(x = PC1, y = PC2, label = variable),
    size = 2
  ) +
  labs(
    "x" = perc_label(pc_res$Male, 1),
    "y" = perc_label(pc_res$Male, 2)
  ) +
  coord_fixed(asp_ratio)
ggsave("../chapter/figure/pca/loadings.png")

## and study the scores
ggplot(scores) +
  geom_point(
    aes(x = PC1, y = PC2, col = weight_dxa)
  ) +
  labs(
    "x" = perc_label(pc_res$Male, 1),
    "y" = perc_label(pc_res$Male, 2)
  ) +
  scale_color_viridis() +
  coord_fixed(ratio = asp_ratio)
ggsave("../chapter/figure/pca/scores_weight.png")

##  also study scores in relation to overall ruminoccocus / lachospiraceae ratio
family_means <- mseqtab %>%
  group_by(family, Number) %>%
  summarise(family_mean = mean(value)) %>%
  spread(family, family_mean) %>%
  group_by(Number) %>%
  summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)

scores <- scores %>%
  left_join(family_means)
ggplot(scores) +
  geom_point(
    aes(x = PC1, y = PC2, col = rl_ratio)
  ) +
  labs(
    "x" = perc_label(pc_res$Male, 1),
    "y" = perc_label(pc_res$Male, 2)
  ) +
  scale_color_viridis() +
  coord_fixed(ratio = asp_ratio)
ggsave("../chapter/figure/pca/scores_rl_ratio.png")
