#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Preprocessing steps to apply before applying dimensionality reduction
## methods.
##
## author: sankaran.kris@gmail.com
## date: 09/27/2017

library("phyloseq")
library("DESeq2")
library("tidyverse")
library("forcats")

###############################################################################
## Prepare opts and read data
###############################################################################
merge_default_opts <- function(opts = list()) {
  default_opts <- list(
    gender = "Male",
    sf_quantile = 0.95,
    filt_k = 0.5,
    filt_a = 0
  )
  modifyList(default_opts, opts)
}

read_data <- function(data_dir = "../data/") {
  seqtab <- readRDS(file.path(data_dir, "seqtab.rds"))
  bc <- readRDS(file.path(data_dir, "sample_data_bc.rds"))
  colnames(seqtab) <- paste0("species_", seq_len(ntaxa(seqtab)))

  taxa <- readRDS(file.path(data_dir, "taxa.rds")) %>%
    data.frame() %>%
    rownames_to_column("seq")
  taxa$seq_num <- colnames(seqtab)
  taxa$family <- fct_lump(taxa$Family, n = 7)
  list("taxa" = taxa, "seqtab" = seqtab, "bc" = bc)
}

process_data <- function(seqtab, bc, taxa, opts = list()) {
  opts <- merge_default_opts(opts)
  mseqtab <- seqtab %>%
    melt(varnames = c("Number", "seq_num")) %>%
    left_join(bc) %>%
    left_join(taxa)

  ## normalize with DESeq2
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

  bc_mat <- data.frame(bc) %>%
    filter(gender == opts$gender) %>%
    select(-id, -gender)
  rownames(bc_mat) <- bc_mat$Number
  bc_mat$Number <- NULL
  list("x_seq" = x_seq, "bc" = bc_mat, "mseqtab" = mseqtab)
}

family_means <- function(mseqtab) {
  processed$mseqtab %>%
    group_by(family, Number) %>%
    summarise(family_mean = mean(value)) %>%
    spread(family, family_mean) %>%
    group_by(Number) %>%
    summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)
}
