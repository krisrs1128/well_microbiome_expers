#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Helper functions for plotting multitable scores
##
## author: sankaran.kris@gmail.com
## date: 09/26/2017

library("ggplot2")
library("ggrepel")

#' Combine Loadings
#'
#' This combines loadings across data types. It's just a helper for the
#' WELL-China dataset
prepare_loadings <- function(loadings_list, types, K = 3) {
  df_list <- list()
  for (m in seq_along(loadings_list)) {
    colnames(loadings_list[[m]]) <- NULL
    df_list[[m]] <- data.frame(
      "variable" = rownames(loadings_list[[m]]),
      "seq_num" = NA,
      "type" = types[m],
      "Axis" = loadings_list[[m]][, 1:K]
    )
  }

  seq_ix <- which(types == "seq")
  if (length(seq_ix) > 0) {
    df_list[[seq_ix]]$seq_num <- df_list[[seq_ix]]$variable
  }

  do.call(rbind, df_list)
}

prepare_scores <- function(scores_list, types, K = 3) {
  df_list <- list()
  for (m in seq_along(scores_list)) {
    colnames(scores_list[[m]]) <- NULL
    df_list[[m]] <- data.frame(
      "Number" = rownames(scores_list[[m]]),
      "type" = types[m],
      "Axis" = scores_list[[m]][, 1:K]
    )
  }

  do.call(rbind, df_list)
}

melt_scores <- function(scores) {
  scores %>%
    select(Number, type, starts_with("Axis")) %>%
    gather(comp, value, starts_with("Axis")) %>%
    unite(comp_type, comp, type) %>%
    spread(comp_type, value)
}

reshape_posterior_score <- function(xi, bc) {
  xi %>%
    mutate(
      Number = bc$Number[row],
      col = paste0("topic_", col)
    ) %>%
    spread(col, value) %>%
    left_join(bc)
}

perc_label <- function(eigs, i) {
  perc <- 100 * eigs[i] / sum(eigs)
  sprintf("Axis %s [%s%%]", i, round(perc, 2))
}

plot_topics <- function(loadings) {
  taxa_hc <- loadings %>%
    select(starts_with("Axis")) %>%
    as.matrix() %>%
    dist() %>%
    hclust()

  mloadings <- loadings %>%
    filter(type == "seq") %>%
    select(variable, starts_with("Axis"), family) %>%
    gather(topic, loading, -variable, -family)

  mloadings$variable <- factor(
    mloadings$variable,
    levels = loadings$seq_num[taxa_hc$order]
  )
  mloadings$topic <- gsub("Axis\\.", "Topic ", mloadings$topic)

  ggplot(mloadings) +
    geom_hline(yintercept = 0) +
    geom_point(
      aes(
        x = variable,
        col = family,
        y = loading
      )
    ) +
    facet_grid(topic ~ family, scale = "free_x", space = "free") +
    theme(
      panel.spacing.x = unit(0, "cm"),
      axis.text.x = element_blank(),
      strip.text.x = element_blank()
    )
}


plot_loadings <- function(loadings, eigs, size_breaks = c(-5, 5), a = 1) {
  ggplot(loadings) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_point(
      data = loadings %>%
        filter(type == "seq"),
      aes(x = Axis.1, y = Axis.2, size = Axis.3, col = family),
      alpha = a
    ) +
    geom_text_repel(
      data = loadings %>%
        filter(type == "body_comp"),
      aes(x = Axis.1, y = Axis.2, label = variable, size = Axis.3),
      segment.size = 0.3,
      segment.alpha = 0.5,
      force = 0.05
    ) +
    labs(
      "x" = perc_label(eigs, 1),
      "y" = perc_label(eigs, 2),
      "col" = "Family"
    ) +
    scale_size_continuous(
      range = c(0, 4),
      breaks = size_breaks
    ) +
    coord_fixed(sqrt(eigs[2] / eigs[1]))
}

plot_scores <- function(scores, col_var, col_label, eigs, size_breaks = c(-8, 8)) {
  ggplot() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_point(
      data = scores,
      aes_string(
        x = "Axis.1",
        y = "Axis.2",
        size = "Axis.3",
        col = col_var
      )
    ) +
    labs(
      "col" = col_label,
      "x" = perc_label(eigs, 1),
      "y" = perc_label(eigs, 2)
    ) +
    scale_size_continuous(range = c(0, 1.5), breaks = size_breaks) +
    coord_fixed(sqrt(eigs[2] / eigs[1]))
}

link_scores <- function(mscores, alpha = 0.1) {
  geom_segment(
    data = mscores,
    aes(
      x = Axis.1_body_comp, xend = Axis.1_seq,
      y = Axis.2_body_comp, yend = Axis.2_seq,
      size = (Axis.3_body_comp + Axis.3_seq) / 2
    ),
    alpha = alpha
  )
}

#' Useful for plotting scores for several features
plot_scores_wrapper <- function(xi, raw, processed, scv) {
  xi_df <- data.frame(
    "Axis" = xi,
    "Number" = rownames(processed$bc)
  ) %>%
    left_join(raw$bc) %>%
    left_join(family_means(processed$mseqtab))
  list(
    plot_scores(xi_df, "age", "age", c(1, 1)) + scv,
    plot_scores(xi_df, "bmi", "BMI", c(1, 1)) + scv,
    plot_scores(xi_df, "Trunk_LM", "Trunk LM", c(1, 1)) + scv,
    plot_scores(xi_df, "rl_ratio", "Rum. / Lachn.", c(1, 1)) + scv
  )
}

#' Average across first dimension
slice_mean <- function(x) {
  apply(x, c(2, 3), mean)
}

#' Means from posterior samples
parameter_means <- function(theta_samples) {
  theta_hat <- list()
  for (i in seq_along(theta_samples)) {
    if (length(dim(theta_samples[[i]])) > 2) {
      theta_hat[[i]] <- slice_mean(theta_samples[[i]])
    } else {
      theta_hat[[i]] <- mean(theta_samples[[i]])
    }
  }

  names(theta_hat) <- names(theta_samples)
  theta_hat
}

melt_parameters <- function(theta_samples) {
  mtheta <- list()
  for (i in seq_along(theta_samples)) {
    if (length(dim(theta_samples[[i]])) > 2) {
      mtheta[[i]] <- theta_samples[[i]] %>%
        melt(
          varnames = c("iteration", "row", "col")
        )
    } else {
      mtheta[[i]] <- theta_samples[[i]] %>%
        melt(value.name = names(theta_samples)[i])
    }
  }
  names(mtheta) <- names(theta_samples)

  mtheta
}
