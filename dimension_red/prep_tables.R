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
library("reshape2")

###############################################################################
## Prepare opts and read data
###############################################################################
merge_default_opts <- function(opts = list()) {
  default_opts <- list(
    gender = "Female", ## has more data
    sf_quantile = 0.95,
    filt_k = 0.5,
    filt_a = 0,
    rlog = TRUE
  )
  modifyList(default_opts, opts)
}

read_data <- function(data_dir = "../data/") {
  seqtab <- readRDS(file.path(data_dir, "seqtab.rds"))
  bc <- readRDS(file.path(data_dir, "sample_data_bc.rds"))
  colnames(bc) <- tolower(colnames(bc))
  colnames(seqtab) <- paste0("species_", seq_len(ntaxa(seqtab)))
  bc_full <- read_csv(file.path(data_dir, "WELL_China_1969_7.25.2017.csv")) %>%
    rename(
      gender_it = "gender",
      age_it = "age",
      height = "height_dxa",
      weight = "weight_dxa"
    ) %>%
    mutate(
      aoi = android_fm / gynoid_fm
    ) %>%
    select_at(
      .vars = setdiff(tolower(colnames(bc)), c("number", "aoi"))
    )

  taxa <- readRDS(file.path(data_dir, "taxa.rds")) %>%
    data.frame() %>%
    rownames_to_column("seq")
  taxa$seq_num <- colnames(seqtab)
  list("taxa" = taxa, "seqtab" = seqtab, "bc" = bc, "bc_full" = bc_full)
}

process_data <- function(seqtab, bc, bc_full, taxa, opts = list()) {
  opts <- merge_default_opts(opts)

  ## filter counts
  seqtab <- seqtab[bc$number, ]
  seqtab <- seqtab[bc$gender == opts$gender, ]
  seqtab <- seqtab %>%
    filter_taxa(function(x) mean(x > opts$filt_a) > opts$filt_k, prune = TRUE)
  taxa <- taxa %>%
    filter(seq_num %in% colnames(seqtab)) %>%
    mutate(
      family = fct_lump(Family, n = 9, ties.method = "first")
    )
  taxa$family <- factor(
    taxa$family,
    levels = names(sort(table(taxa$family), decreasing = TRUE))
  )

  mseqtab <- seqtab %>%
    melt(varnames = c("number", "seq_num")) %>%
    left_join(bc) %>%
    left_join(taxa)

  dds <- DESeqDataSetFromMatrix(
    countData = t(get_taxa(seqtab)),
    colData = data.frame("unused" = rep(1, nrow(seqtab))),
    design = ~1
  )

  ## are quantiles for one sample systematically larger than those for others (if
  ## so, give it a large size factor). Basically related to sequencing depth.
  ##
  ## code copied from here: https://support.bioconductor.org/p/76548/
  if (opts$rlog) {
    fname <- paste0(as.character(opts), collapse = "")
    fpath <- file.path("..", "data", digest::digest(fname))
    if (file.exists(fpath)) {
      dds <- get(load(fpath))
    } else {
      qs <- apply(counts(dds), 2, quantile, opts$sf_quantile)
      sizeFactors(dds) <- qs / exp(mean(log(qs)))
      dds <- rlog(dds, fitType = "local", betaPriorVar = 0.5)
      save(dds, file = fpath)
    }
  }

  x_seq <- t(assay(dds))
  bc_mat <- data.frame(bc) %>%
    filter(gender == opts$gender) %>%
    select(-id, -gender)
  rownames(bc_mat) <- bc_mat$number
  bc_mat$number <- NULL

  bc_full <- data.frame(bc_full) %>%
    mutate(
      gender = recode(gender, `1` = "Male", `2` = "Female")
    ) %>%
    filter(gender == opts$gender) %>%
    select( -gender) %>%
    na.omit()
  rownames(bc_full) <- bc_full$id
  bc_full$id <- NULL

  list(
    "x_seq" = x_seq,
    "bc" = bc_mat,
    "bc_full" = bc_full,
    "mseqtab" = mseqtab
  )
}

family_means <- function(mseqtab) {
  processed$mseqtab %>%
    group_by(family, Number) %>%
    dplyr::summarise(family_mean = mean(value)) %>%
    spread(family, family_mean) %>%
    group_by(Number) %>%
    dplyr::summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)
}
