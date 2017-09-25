#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## An example studying multiple tables by concatenating and then using PCA.
##
## author: sankaran.kris@gmail.com
## date: //2017

###############################################################################
## Libraries and setup
###############################################################################

library("tidyverse")
library("phyloseq")
library("DESeq2")

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

seqtab <- readRDS("../data/seqtab.rds")
taxa <- readRDS("../data/taxa.rds")
bc <- readRDS("../data/sample_data_bc.rds")

###############################################################################
## Concatenate quantitative body composition (split by gender) and (somewhat
## filtered) microbiome counts, and transformed
###############################################################################
bc_vars <- setdiff(colnames(bc), c("Number", "id"))

combined <- data.frame(bc, seqtab[bc$Number, ])
dim(combined)
cbind(rownames(combined), combined$Number)

library("DESeq2")
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
x <- assay(vsd)
## I would prefer rlog, but I haven't been patient enough to let it finish
## (vst is much faster)
## rlog(dds, fitType = "local")

## for the record, the difference between asinh and rlog
hist(x, breaks = 100)
hist(asinh(get_taxa(seqtab)), breaks = 100)
