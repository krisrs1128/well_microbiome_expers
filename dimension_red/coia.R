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
x_seq <- x_seq[bc$gender == opts$gender, ]

###############################################################################
## Run CoIA on the two (scaled) tables
###############################################################################
bc_mat <- data.frame(bc) %>%
  filter(gender == opts$gender) %>%
  select(-id, -gender)
rownames(bc_mat) <- bc_mat$Number
bc_mat$Number <- NULL

dudi1 <- dudi.pca(x_seq, scan = FALSE, nf = 3)
dudi2 <- dudi.pca(bc_mat, scan = FALSE, nf = 3)
coin1 <- coinertia(dudi1, dudi2, scan = FALSE, nf = 3)
