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
    rlog = c("rlog" = TRUE, "sf_quantile" = 0.95),
    filt_k = 0.5,
    filt_a = 0,
    outdir = "../data"
  )
  modifyList(default_opts, opts)
}

read_data <- function(data_dir = "../data/") {
  seqtab <- readRDS(file.path(data_dir, "seqtab.rds"))
  bc <- readRDS(file.path(data_dir, "sample_data_bc.rds")) %>%
    mutate(id = as.character(id))
  colnames(bc) <- tolower(colnames(bc))
  colnames(seqtab) <- paste0("species_", seq_len(ntaxa(seqtab)))
  bc_full <- read_csv(file.path(data_dir, "WELL_China_1969_7.25.2017.csv")) %>%
    rename(
      gender = gender_it,
      age = age_it,
      height_dxa = height,
      weight_dxa = weight
    ) %>%
    mutate(
      aoi = android_fm / gynoid_fm,
      id = as.character(id)
    ) %>%
    select_at(
      .vars = setdiff(tolower(colnames(bc)), "number")
    )

  taxa <- readRDS(file.path(data_dir, "taxa.rds")) %>%
    data.frame() %>%
    rownames_to_column("seq")
  taxa$seq_num <- colnames(seqtab)
  tree <- readRDS(file.path(data_dir, "phylo_tree.rds"))
  list("taxa" = taxa, "seqtab" = seqtab, "bc" = bc, "bc_full" = bc_full, "tree" = tree)
}

#' Prepare sample data
#'
#' reorder columns and convert to ratios.
prep_sample <- function(survey, scale_sample_data = FALSE) {
  sample <- data.frame(survey) %>%
    mutate(
      id = as.character(id),
      gender = as.factor(gender)
    ) %>%
    select( ## reorder columns
      id, age, gender, height_dxa, weight_dxa,
      bmi, aoi, ends_with("_fm"), ends_with("_lm"),
      starts_with("diet_")
    )%>%
    mutate_at( ## fm ratios
      vars(ends_with("fm"), -total_fm),
      .funs = funs(. / total_fm)
    ) %>%
    mutate_at( ## lm ratios
      vars(ends_with("lm"), -total_lm),
      .funs = funs(. / total_lm)
    ) %>%
    mutate_at(
      vars(starts_with("diet_")),
      .funs = log
    ) %>%
    mutate( ## new variables
      fat_lean_ratio = total_fm / total_lm,
      total_fm = total_fm / (weight_dxa * 1000),
      total_lm = total_lm / (weight_dxa * 1000)
    )

  ## scale body measurements by gender
  if (scale_sample_data) {
    sample <- sample %>%
      group_by(gender) %>%
      mutate_if(is.numeric, safe_scale) %>%
      ungroup()
  }

  as.data.frame(sample)
}

prepare_taxa <- function(taxa) {
  taxa <- taxa %>%
    mutate(
      family = fct_lump(Family, n = 9, ties.method = "first")
    )
  taxa$family <- factor(
    taxa$family,
    levels = names(sort(table(taxa$family), decreasing = TRUE))
  )
  tax_cols <- colnames(taxa)
  taxa <- tax_table(taxa)
  colnames(taxa) <- tax_cols
  taxa_names(taxa) <- taxa[, "seq_num"]
  taxa
}


#' scaling that ignores NAs
safe_scale <- function(x) {
  x[!is.finite(x)] <- NA
  z <- x - mean(x, na.rm = TRUE)
  z / sd(x, na.rm = TRUE)
}

process_data <- function(seqtab, bc, bc_full, taxa, tree, opts = list()) {
  opts <- merge_default_opts(opts)

  ## preparing taxa and survey data
  taxa <- prepare_taxa(taxa) %>%
    filter(seq_num %in% colnames(seqtab))
  taxa_names(tree) <- taxa[, "seq_num"]

  bc_full <- prep_sample(bc_full, opts$scale_sample_data) %>%
    left_join(bc %>% select(id, number))
  rownames(bc) <- bc$number

  ## combine into ps
  ps <- phyloseq(
    otu_table(seqtab, taxa_are_rows = FALSE),
    sample_data(bc),
    taxa,
    phy_tree(tree)
  ) %>%
    subset_samples(gender == opts$gender) %>%
    filter_taxa(
      function(x) {
        mean(x > opts$filter_a) > opts$filter_k
      },
      prune = TRUE
    )

  if (opts$rlog$rlog) {
    fname <- paste0(as.character(opts), collapse = "")
    dir.create("..", "data")
    fpath <- file.path(opts$outdir, digest::digest(fname))
    if (file.exists(fpath)) {
      rlog_dds <- get(load(fpath))
    } else {
      rlog_dds <- do.call(rlog_ps, c(fpath, ps, opts$rlog))
    }
    otu_table(ps)@.Data <- t(assay(rlog_dds))
  }

  list(
    "ps" = ps,
    "bc_full" = bc_full,
    "mseqtab" = melt_ps(ps)
  )
}

melt_ps <- function(ps) {
  get_taxa(ps) %>%
    melt(varnames = c("number", "seq_num")) %>%
    left_join(sample_data(ps)) %>%
    left_join(data.frame(tax_table(ps)))
}

family_means <- function(mseqtab) {
  processed$mseqtab %>%
    group_by(family, Number) %>%
    dplyr::summarise(family_mean = mean(value)) %>%
    spread(family, family_mean) %>%
    group_by(Number) %>%
    dplyr::summarise(rl_ratio = Ruminococcaceae / Lachnospiraceae)
}
