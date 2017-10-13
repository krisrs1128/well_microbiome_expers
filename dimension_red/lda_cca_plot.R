#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Not generic enough to be included in plot.R, but don't want to copy code
## everywhere.
##
## author: sankaran.kris@gmail.com
## date: //2017

lda_cca_plots <- function(mdist, seq_families, processed, opts) {
  scv <- scale_color_viridis(
    guide = guide_colorbar(barwidth = 0.15, ticks = FALSE)
  )
  ggplot() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_point(
      data = mdist$xi_s,
      aes(
        x = axis_1,
        y = axis_2,
        col = Total_LM
      ),
      size = 0.1,
      alpha = 0.01
    ) +
    coord_equal() +
    scv +
    labs(x = "Axis 1", y = "Axis 2", col = "Total LM")

  ggsave(
    sprintf("%s/shared_scores_lm_posterior.png", opts$outdir),
    width = 5.63, height = 3.42
  )

  ## using the scores from just the bc table
  ggplot() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_point(
      data = mdist$xi_y,
      aes(
        x = axis_1,
        y = axis_2,
        col = Total_LM
      ),
      size = 0.1,
      alpha = 0.1
    ) +
    coord_equal() +
    scv +
    labs(x = "Axis 1", y = "Axis 2", col = "Total LM")
  ggsave(
    sprintf("%s/unshared_scores_lm_posterior.png", opts$outdir),
    width = 9.24,
    height = 2.5
  )

  ## Now plotting loadings boxplots
  mdist$Wx$seq_num <- colnames(processed$x_seq)[mdist$Wx$row]
  mdist$Wx <- mdist$Wx %>%
    left_join(seq_families)
  mdist$Wx$seq_num <- factor(
    mdist$Wx$seq_num,
    levels = names(sort(colSums(processed$x_seq), decreasing = TRUE))
  )

  wx_summary <- mdist$Wx %>%
    group_by(col, seq_num) %>%
    dplyr::summarise(
      family = family[1],
      lower = quantile(value, 0.25),
      upper = quantile(value, 0.75),
      med = median(value)
    )

  ggplot(wx_summary) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_pointrange(
      aes(x = seq_num, ymin = lower, ymax = upper, y = med, col = family, fill = family),
      fatten = 1.2,
      size = 0.1
    ) +
    facet_grid(col ~ family, scale = "free_x", space = "free_x") +
    theme(
      panel.spacing.x = unit(0, "cm"),
      axis.text.x = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom"
    )
  ggsave(
    sprintf("%s/within_loadings_seq_boxplots.png", opts$outdir),
    width = 7.45,
    height = 3.37
  )

  mdist$Bx$seq_num <- colnames(processed$x_seq)[mdist$Bx$row]
  mdist$Bx <- mdist$Bx %>%
    left_join(seq_families)
  mdist$Bx$seq_num <- factor(
    mdist$Bx$seq_num,
    levels = names(sort(colSums(processed$x_seq), decreasing = TRUE))
  )

  bx_summary <- mdist$Bx %>%
    group_by(col, seq_num) %>%
    dplyr::summarise(
      family = family[1],
      lower = quantile(value, 0.25),
      upper = quantile(value, 0.75),
      med = median(value)
    )

  ggplot(bx_summary) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_pointrange(
      aes(x = seq_num, ymin = lower, ymax = upper, y = med, col = family, fill = family),
      fatten = 1.2,
      size = 0.1
    ) +
    facet_grid(col ~ family, scale = "free_x", space = "free_x") +
    theme(
      panel.spacing.x = unit(0, "cm"),
      axis.text.x = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom"
    )
  ggsave(
    sprintf("%s/between_loadings_seq_boxplots.png", opts$outdir),
    width = 7.45,
    height = 3.37
  )

  ## and finally loadings boxplot for body composition variables
  site_ordered <- c(
    "aoi", "age", "height_dxa", "weight_dxa",
    "bmi", "android_fm", "android_lm", "gynoid_fm", "gynoid_lm", "l_trunk_fm",
    "l_trunk_lm", "r_trunk_fm", "r_trunk_lm", "trunk_fm", "trunk_lm",
    "l_total_fm", "l_total_lm", "r_total_fm", "r_total_lm", "total_fm",
    "total_lm", "l_leg_fm", "l_leg_lm", "r_leg_fm", "r_leg_lm", "legs_fm",
    "legs_lm", "l_arm_fm", "l_arm_lm", "r_arm_fm", "r_arm_lm", "arms_fm",
    "arms_lm"
  )
  mass_type_ordered <- c(
    site_ordered[!grepl("fm|lm", site_ordered)],
    site_ordered[grepl("fm", site_ordered)],
    site_ordered[grepl("lm", site_ordered)]
  )

  mdist$Wy$variable <- tolower(colnames(processed$bc)[mdist$Wy$row])
  mdist$Wy$variable <- factor(
    mdist$Wy$variable,
    levels = mass_type_ordered
  )

  ggplot(mdist$Wy) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_boxplot(
      aes(
        x = variable,
        y = value
      ),
      outlier.size = 1
    ) +
    facet_grid(col ~ .) +
    theme(axis.text.x = element_text(angle = -90))

  ggsave(
    sprintf("%s/within_loadings_body_comp_boxplots.png", opts$outdir),
    width = 6.44,
    height = 3.95
  )

  mdist$By$variable <- tolower(colnames(processed$bc)[mdist$By$row])
  mdist$By$variable <- factor(
    mdist$By$variable,
    levels = mass_type_ordered
  )

  ggplot(mdist$By) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_boxplot(
      aes(
        x = variable,
        y = value
      ),
      outlier.size = 1
    ) +
    facet_grid(col ~ .) +
    theme(axis.text.x = element_text(angle = -90))

  ggsave(
    sprintf("%s/between_loadings_body_comp_boxplots.png", opts$outdir),
    width = 6.44,
    height = 3.95
  )
}
