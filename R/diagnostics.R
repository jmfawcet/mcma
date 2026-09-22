# Model diagnostics
# =================

#' Posterior predictive checks
#'
#' Compares replicated data from the posterior predictive distribution to the
#' observed data.
#'
#' @param fit A brmsfit object.
#' @param n_draws Number of posterior predictive draws.
#' @param summary Logical; if TRUE, returns summary statistics.
#' @param plot Logical; if TRUE, also produces a visual PPC comparison.
#' @return A tibble with PPC summary statistics, or the raw yrep matrix.
#' @export
mcma_ppc <- function(fit, n_draws = 500, summary = TRUE, plot = FALSE) {

  # Generate possible data sets under the fitted model and compare their study
  # proportions with the observed proportions. Each study contributes equally
  # to these summaries.

  # brms::posterior_predict() for nonlinear models evaluates the formula
  # in R and needs inv_logit defined. Define it in the calling environment
  # if not already available.
  if (!exists("inv_logit", envir = globalenv())) {
    assign("inv_logit", stats::plogis, envir = globalenv())
  }

  # Each row is one replicated dataset, with one positive count per observed
  # study. Divide by the original sample sizes to compare prevalence
  # proportions.
  yrep <- brms::posterior_predict(fit, ndraws = n_draws)
  y    <- fit$data$y
  n    <- fit$data$n

  if (!summary) return(yrep)

  # Observed proportions
  p_obs <- y / n

  # Replicated proportions
  p_rep <- sweep(yrep, 2, n, "/")

  # Compare the unweighted mean prevalence across studies in the replicated
  # and observed datasets.
  # Mean discrepancy
  mean_obs <- mean(p_obs)
  mean_rep <- rowMeans(p_rep)
  ppc_md   <- mean(mean_rep) - mean_obs

  # A ratio near one means the replicated data have a similar spread of study
  # proportions to the observed data.
  # Variance ratio
  var_obs  <- stats::var(p_obs)
  var_rep  <- apply(p_rep, 1, stats::var)
  ppc_var_ratio <- mean(var_rep) / var_obs

  # Calculate the fraction of replicated datasets whose mean or variance
  # exceeds the observed value. Posterior predictive tail areas.
  # Bayesian p-values
  ppc_p_mean <- mean(mean_rep >= mean_obs)
  ppc_p_var  <- mean(var_rep >= var_obs)

  # Return the discrepancy summaries together; summary = FALSE instead
  # returns the replicated counts above.
  out <- tibble::tibble(
    ppc_md        = ppc_md,
    ppc_var_ratio = ppc_var_ratio,
    ppc_p_mean    = ppc_p_mean,
    ppc_p_var     = ppc_p_var,
    n_obs         = length(y)
  )

  if (plot) {
    message("PPC plotting not yet implemented. Use brms::pp_check() directly.")
  }

  out
}


#' Quick convergence summary
#'
#' @param ... Named or unnamed brmsfit objects. Unnamed inputs receive
#'   `model_` followed by their position in the argument list.
#' @return A tibble with one row per model.
#' @export
mcma_convergence <- function(...) {

  # Collect sampling diagnostics across models. Rhat measures chain
  # agreement; effective sample sizes describe how much independent
  # information the correlated draws contain.

  fits <- list(...)
  fit_names <- names(fits)
  if (is.null(fit_names)) fit_names <- rep("", length(fits))
  unnamed <- is.na(fit_names) | !nzchar(fit_names)
  fit_names[unnamed] <- paste0("model_", which(unnamed))

  # Inspect each supplied fit independently so one invalid object does not
  # prevent summaries of the remaining models.
  # Index by position so missing or repeated names never select the wrong fit.
  rows <- lapply(seq_along(fits), function(i) {
    nm <- fit_names[i]
    fit <- fits[[i]]
    if (!inherits(fit, "brmsfit")) {
      rlang::warn(sprintf("'%s' is not a brmsfit, skipping.", nm))
      return(NULL)
    }

    # Rhat should be close to one when chains agree. Report the largest value
    # to expose the least well-mixed parameter.
    # Rhat
    rhat_vals <- tryCatch({
      rh <- posterior::summarise_draws(posterior::as_draws(fit), "rhat")
      rh$rhat
    }, error = function(e) NA_real_)

    # Bulk ESS concerns the centre of the posterior; tail ESS concerns its
    # extremes. The smallest values identify the weakest sampling precision.
    # ESS
    ess <- tryCatch({
      posterior::summarise_draws(
        posterior::as_draws(fit),
        "ess_bulk", "ess_tail"
      )
    }, error = function(e) NULL)

    # Divergences flag trajectories where the sampler had difficulty
    # exploring the posterior geometry.
    # Divergences
    np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    ndiv <- if (!is.null(np)) {
      sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
    } else {
      NA_integer_
    }

    # Tree depth measures the sampler work used for a transition. This
    # reports the observed maximum, not a count of transitions that hit the
    # configured limit.
    # Max treedepth
    max_td <- if (!is.null(np)) {
      max(np$Value[np$Parameter == "treedepth__"], na.rm = TRUE)
    } else {
      NA_integer_
    }

    # Reduce the parameter-level diagnostics to one row per fitted model,
    # preserving missing values when extraction fails.
    tibble::tibble(
      model         = nm,
      max_rhat      = if (all(is.na(rhat_vals))) NA_real_ else max(rhat_vals, na.rm = TRUE),
      min_ess_bulk  = if (is.null(ess)) NA_real_ else min(ess$ess_bulk, na.rm = TRUE),
      min_ess_tail  = if (is.null(ess)) NA_real_ else min(ess$ess_tail, na.rm = TRUE),
      n_divergences = ndiv,
      max_treedepth = max_td
    )
  })

  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}


#' Plot joint Se/Sp posterior
#'
#' Visualises the joint posterior distribution of Se and Sp from a joint
#' model (M8).
#'
#' @param fit A brmsfit from `mcma_fit_joint()`.
#' @param bounded Logical; auto-detected if NULL.
#' @param true_se True Se for overlay (simulation validation).
#' @param true_sp True Sp for overlay (simulation validation).
#' @param show_iso Logical; show the iso-correction line.
#' @param p_obs_screen Observed screening prevalence for iso line (computed from posterior if NULL).
#' @param show_marginals Logical; show marginal density plots.
#' @param show_border Logical; if TRUE, draw a box around the plot panel.
#' @param plot_title Title for the plot (NULL removes title).
#' @param other_options A list of ggplot elements to add to the plot.
#' @param n_draws Number of draws to subsample for plotting.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot_joint_sesp <- function(fit,
                            bounded        = NULL,
                            true_se        = NULL,
                            true_sp        = NULL,
                            show_iso       = TRUE,
                            p_obs_screen   = NULL,
                            show_marginals = FALSE,
                            show_border    = TRUE,
                            plot_title     = "Joint Posterior for Se and Sp",
                            other_options  = NULL,
                            n_draws        = 2000,
                            ...) {

  # Plot paired sensitivity and specificity draws so their posterior
  # dependence remains visible. Optional overlays show known truth and a
  # constant-correction curve.
  if (show_marginals && !requireNamespace("ggExtra", quietly = TRUE)) {
    stop("Package 'ggExtra' is required for marginal density plots. ",
         "Install with: install.packages('ggExtra')")
  }

  # Extract the shared screening accuracy on its probability scale, using the
  # bounded setting saved with the fit unless overridden.
  # Extract draws
  sesp <- extract_sesp(fit, bounded = bounded, level = "global", summary = FALSE)
  se_draws <- sesp$se
  sp_draws <- sesp$sp

  # Prevalence draws for iso line
  prev_draws <- plogis(posterior::as_draws_df(fit)$b_pi_Intercept)

  # Apply the same random indices to Se, Sp, and prevalence so each plotted
  # point still represents one joint posterior draw.
  # Subsample
  if (length(se_draws) > n_draws) {
    idx <- sample(length(se_draws), n_draws)
    se_draws   <- se_draws[idx]
    sp_draws   <- sp_draws[idx]
    prev_draws <- prev_draws[idx]
  }

  plot_df <- data.frame(se = se_draws, sp = sp_draws, prevalence = prev_draws)

  # Combine a lightly shaded scatterplot with filled density contours to show
  # where posterior draws concentrate.
  # Core plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$sp, y = .data$se)) +
    ggplot2::geom_point(alpha = 0.05, size = 0.5, color = "steelblue") +
    ggplot2::stat_density_2d(
      ggplot2::aes(fill = ggplot2::after_stat(.data$level)),
      geom = "polygon", alpha = 0.4
    ) +
    ggplot2::scale_fill_viridis_c(option = "D", guide = "none") +
    ggplot2::labs(
      x = "Specificity (Sp)",
      y = "Sensitivity (Se)",
      title = plot_title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::xlim(c(0.5, 1)) +
    ggplot2::ylim(c(0.5, 1))

  # Panel border
  if (show_border) {
    p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.7
      )
    )
  }

  # True values overlay
  if (!is.null(true_se) && !is.null(true_sp)) {
    p <- p + ggplot2::geom_point(
      data = data.frame(sp = true_sp, se = true_se),
      ggplot2::aes(x = .data$sp, y = .data$se),
      colour = "red", size = 3, shape = 4, stroke = 2
    )
  }

  # Draw the ordinary Rogan-Gladen curve for a fixed observed and true
  # prevalence. This overlay currently uses baseline Se/Sp without c/o
  # adjustments.
  # Iso-correction line
  if (show_iso) {
    theta <- median(plot_df$prevalence, na.rm = TRUE)

    if (!is.na(theta) && theta > 0 && theta < 1) {
      if (is.null(p_obs_screen)) {
        p_obs_draws <- plot_df$se * plot_df$prevalence +
          (1 - plot_df$sp) * (1 - plot_df$prevalence)
        p_obs_screen <- median(p_obs_draws, na.rm = TRUE)
      }

      # Solve the observation equation for Sp at each Se value, then keep the
      # part of the curve inside the displayed accuracy range.
      se_seq <- seq(0.5, 1, by = 0.005)
      sp_iso <- 1 - (p_obs_screen - theta * se_seq) / (1 - theta)
      iso_df <- data.frame(se = se_seq, sp = sp_iso)
      iso_df <- iso_df[iso_df$sp > 0.5 & iso_df$sp < 1 & (iso_df$se + iso_df$sp) > 1, ]

      if (nrow(iso_df) > 0) {
        p <- p + ggplot2::geom_line(
          data = iso_df,
          ggplot2::aes(x = .data$sp, y = .data$se),
          inherit.aes = FALSE,
          colour = "red", linetype = "dashed", linewidth = 0.8
        )
      }
    }
  }

  # Add caller-supplied plot elements before optional marginal densities
  # change the returned object class.
  # Additional ggplot layers
  if (!is.null(other_options)) {
    p <- purrr::reduce(other_options, `+`, .init = p)
  }

  # Marginal densities (must be last — changes object class)
  if (show_marginals) {
    p <- ggExtra::ggMarginal(p, type = "density", fill = "grey80",
                             color = "grey40", alpha = 0.6)
  }

  p
}
