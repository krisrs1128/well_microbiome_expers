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

###############################################################################
## concatenate quantitative body composition (split by gender) and (somewhat
## filtered) microbiome counts, and transformed
###############################################################################
bc[, -c(1:3)] <- bc[, -c(1:3)]
bc_mat <- data.frame(bc) %>%
  filter(gender == "Male") %>%
  select(-id, -Number, -gender) %>%
  as.matrix() %>%
  scale()

cca_res <- cancor(bc_mat, x_seq[bc$gender == "Male", ])
cca_res$xcoef[, 1:2]
cca_res$ycoef[, 1:2]
