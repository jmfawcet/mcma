# S3 plot methods
# ================

#' Plot mcma_priors
#'
#' Shows implied prior densities on the probability scale for Se, Sp,
#' and prevalence.
#'
#' @param x An mcma_priors object.
#' @param bounded Logical; if TRUE, shows bounded prior densities.
#' @param prev_center Prevalence prior centre (probability scale). If NULL,
#'   prevalence panel is omitted.
#' @param prev_sd Prevalence prior SD on logit scale.
#' @param n_samples Number of samples for density estimation.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.mcma_priors <- function(x,
                             bounded     = FALSE,
                             prev_center = NULL,
                             prev_sd     = 1.5,
                             n_samples   = 10000,
                             ...) {

  # Simulate from the same transformed-normal accuracy priors used for
  # fitting, so their shapes can be inspected on the probability scale.
  screen <- x[!x$is_gold, ]
  if (nrow(screen) == 0) {
    rlang::warn("No screening measures to plot.")
    return(ggplot2::ggplot())
  }

  plot_data <- list()

  for (i in seq_len(nrow(screen))) {
    mid <- screen$measure_id[i]
    se_c <- screen$se[i]
    sp_c <- screen$sp[i]
    kap  <- screen$kappa[i]

    # Compute the same approximate logit-scale SDs used by the prior builder,
    # including the bounded chain-rule adjustment when requested.
    se_sd <- prior_sd_from_kappa(se_c, kap, bounded = bounded)
    sp_sd <- prior_sd_from_kappa(sp_c, kap, bounded = bounded)

    # Draw on the fitted normal scale and back-transform so the graph shows
    # prior uncertainty in interpretable probabilities.
    if (bounded) {
      se_mu <- stats::qlogis(map_to_inner(se_c))
      sp_mu <- stats::qlogis(map_to_inner(sp_c))
      se_samp <- 0.5 + 0.5 * stats::plogis(stats::rnorm(n_samples, se_mu, se_sd))
      sp_samp <- 0.5 + 0.5 * stats::plogis(stats::rnorm(n_samples, sp_mu, sp_sd))
    } else {
      se_mu <- stats::qlogis(se_c)
      sp_mu <- stats::qlogis(sp_c)
      se_samp <- stats::plogis(stats::rnorm(n_samples, se_mu, se_sd))
      sp_samp <- stats::plogis(stats::rnorm(n_samples, sp_mu, sp_sd))
    }

    plot_data <- c(plot_data, list(
      data.frame(measure = mid, param = "Se", value = se_samp),
      data.frame(measure = mid, param = "Sp", value = sp_samp)
    ))
  }

  # Optionally add the prevalence prior as a separate panel; it always uses
  # the ordinary logit transformation.
  # Prevalence
  if (!is.null(prev_center)) {
    prev_mu <- stats::qlogis(prev_center)
    prev_samp <- stats::plogis(stats::rnorm(n_samples, prev_mu, prev_sd))
    plot_data <- c(plot_data, list(
      data.frame(measure = "prevalence", param = "Prevalence", value = prev_samp)
    ))
  }

  # Combine simulated values from all instruments, then estimate smooth
  # densities for display.
  df <- do.call(rbind, plot_data)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value, fill = .data$measure)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::facet_wrap(~ .data$param, scales = "free") +
    ggplot2::labs(x = "Probability", y = "Density",
                  title = "Implied Prior Densities") +
    ggplot2::theme_minimal()
}


#' Plot mcma_kwga results
#'
#' Produces a heatmap of Se/Sp posterior weights from the kernel-weighted grid averaging (KWGA) procedure.
#'
#' @param x An mcma_kwga object.
#' @param show_top_n Number of top-weighted cells to label.
#' @param ... Ignored.
#' @param plot_estimate Accuracy estimate to outline on the heatmap:
#'   `"joint_map"` selects the highest-weight grid pair and
#'   `"marginal_map"` selects the separate Se and Sp marginal modes.
#'   NULL (the default) draws no outline.
#' @return A ggplot object.
#' @export
plot.mcma_kwga <- function(x, show_top_n = 10, plot_estimate = NULL, ...) {

  # Draw one tile per candidate accuracy pair, with colour showing its share
  # of the total kernel weight.

  grid <- x$grid

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$sp, y = .data$se,
                                          fill = .data$weight)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Weight", option = "C") +
    ggplot2::labs(x = "Specificity", y = "Sensitivity",
                  title = "KWGA Weight Surface") +
    ggplot2::theme_minimal()

  # Annotate only the highest-weight candidates, keeping the rest of the
  # heatmap readable.
  # Label top cells
  if (show_top_n > 0 && show_top_n <= nrow(grid)) {
    top <- utils::head(grid, show_top_n)
    p <- p + ggplot2::geom_text(
      data = top,
      ggplot2::aes(label = sprintf("%.1f%%", .data$weight * 100)),
      size = 2.5, colour = "white"
    )
  }

  # Mark estimates
  if(!is.null(plot_estimate) && plot_estimate %in% c('joint_map', 'marginal_map'))
  {
    est <- x$estimates
    if (nrow(est) > 0) {
      wm <- est[est$method == plot_estimate, ]
      if (nrow(wm) == 1) {
        # p <- p + ggplot2::geom_point(
        #   data = data.frame(sp = wm$sp, se = wm$se),
        #   inherit.aes = FALSE,
        #   ggplot2::aes(x = .data$sp, y = .data$se),
        #   colour = "red", size = 4, shape = 4, stroke = 2
        # )
        # Determine tile size from grid spacing

        # Use the grid spacing to size the outline around the selected
        # estimate.
        sp_vals <- sort(unique(grid$sp))
        se_vals <- sort(unique(grid$se))
        tile_w <- if (length(sp_vals) > 1) diff(sp_vals)[1] else 0.025
        tile_h <- if (length(se_vals) > 1) diff(se_vals)[1] else 0.025

        p <- p + ggplot2::annotate(
          "rect",
          xmin = wm$sp - tile_w / 2, xmax = wm$sp + tile_w / 2,
          ymin = wm$se - tile_h / 2, ymax = wm$se + tile_h / 2,
          fill = NA, colour = "red", linewidth = 1.5
        )
      }
    }
  }

  p
}
