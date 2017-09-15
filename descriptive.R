#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Descriptive statistics on the WELL microbiome dataset.
##
## author: sankaran.kris@gmail.com
## date: 09/14/2017

## load packages
library("tidyverse")
library("reshape2")

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


## load data
bc <- readRDS("data/sample_data_bc.rds")
survey <- read_csv("data/WELL_China_1969_7.25.2017.csv")
seqtab <- readRDS("data/seqtab.rds")
taxa <- readRDS("data/taxa.rds")

## getting names to agree across bc and survey -- not sure where aoi came from
bc_survey <- survey %>%
  rename(
    age = age_it,
    gender = gender_it,
    height_dxa = height,
    weight_dxa = weight
  )
bc_cols <- intersect(colnames(bc_survey), tolower(colnames(bc)))
id_map <- setNames(sample_names(bc), bc$id)

## subsetting to the body composition fields
bc_survey <- bc_survey %>%
  select_(.dots = bc_cols) %>%
  mutate(
    id = as.character(id),
    Number = id_map[id]
  ) %>%
  select_(.dots = c("id", "Number", bc_cols)) # reorder columns

## make histograms
histo_plot <- function(msurvey) {
  ggplot(msurvey) +
    geom_histogram(
      aes(x = value, fill = as.factor(gender)),
      bins = 60, position = "identity", alpha = 0.6
    ) +
    facet_wrap(~measurement, scales = "free") +
    theme(
      axis.text.y = element_blank()
    )
}

melt_bc_survey <- bc_survey %>%
  gather(measurement, value, -id, -Number, -gender)
histo_plot(melt_bc_survey)

## reorder columns according to a clustering
clust_order <- function(x) {
  D <- dist(scale(x))
  rownames(x)[hclust(D)$order]
}

q_levels <- clust_order(t(bc_survey[, -c(1, 2)]))
melt_bc_survey <- melt_bc_survey %>%
  mutate(
    measurement = factor(measurement, levels = q_levels)
  )

histo_plot(melt_bc_survey)

## let's make some heatmaps
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

melt_bc_survey <- melt_bc_survey %>%
  group_by(measurement, gender) %>%
  mutate(
    scaled_meas = scale(value)
  )

## reorder rows according to clustering also
id_levels <- clust_order(bc_survey[, -c(1:2)])
id_levels <- bc_survey$id[as.numeric(id_levels)]
melt_bc_survey <- melt_bc_survey %>%
  mutate(
    id = factor(id, levels = id_levels)
  )

ggplot(melt_bc_survey) +
  geom_tile(
    aes(x = measurement, y = id, fill = scaled_meas)
  ) +
  scale_fill_gradient2() +
  facet_grid(gender ~ ., scale = "free", space = "free") +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_blank()
  )

## pairs plot
q_levels
for (i in seq(1, 29, 5)) {
  pairs(bc_survey[, q_levels[i:(i + 5)]])
}

###############################################################################
## And now the microbiome!
###############################################################################
dim(seqtab)
seq_map <- setNames(
  paste0("species_", seq_along(seqtab)),
  colnames(seqtab)
)

colnames(seqtab) <- seq_map[colnames(seqtab)]
mseqtab <- seqtab %>%
  melt(varnames = c("Number", "species")) %>%
  left_join()

ggplot(seqtab) +
  geom_histogram()
