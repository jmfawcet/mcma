# Sensitivity analysis across Se/Sp grid
# ========================================

#' Sensitivity analysis over an Se x Sp grid
#'
#' Fits the corrected model at each point in an Se x Sp grid, updating the
#' mcma_priors object at each point. Models are cached to `model_dir`.
#' Cache names retain up to 15 significant digits while preserving existing
#' two-decimal names for ordinary grids.
#'
#' @param data Data frame.
#' @param priors Base mcma_priors object (Se/Sp centres will be overridden).
#' @param se_grid Numeric vector of Se values.
#' @param sp_grid Numeric vector of Sp values.
#' @param prev_center Prevalence prior centre (probability scale).
#' @param bounded Logical; whether to use bounded parameterization.
#' @param prev_re Prevalence random effects formula.
#' @param sesp_re Se/Sp random effects formula.
#' @param correlated_re Logical; whether to correlate random effects.
#' @param measure_col Measure identifier column.
#' @param model_dir Directory for caching fitted models.
#' @param rerun Logical; if `TRUE`, refit and overwrite cached grid models even
#'   when a cached fit already exists (sets brms `file_refit = "always"`).
#'   Default `FALSE` reuses any cached fits.
#' @param ... Additional arguments passed to `mcma_fit()`.
#' @return An S3 object of class `mcma_sensitivity`.
#' @export
mcma_sensitivity <- function(data,
                             priors,
                             se_grid,
                             sp_grid,
                             prev_center,
                             bounded       = FALSE,
                             prev_re       = ~ (1 | es_id),
                             sesp_re       = ~ (1 | es_id),
                             correlated_re = TRUE,
                             measure_col   = "measure_id",
                             model_dir     = "models/sensitivity/",
                             rerun         = FALSE,
                             ...) {

  # Repeat the corrected fit over candidate prior centres. Each grid row
  # records how the prevalence estimate changes under a different accuracy
  # assumption.

  # Subsets need priors only for the instruments they actually contain.
  priors <- .mcma_match_priors(priors, data[[measure_col]])

  if (!is.null(model_dir) && !dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }

  # Create all prior-centre pairs and omit combinations with a nonpositive
  # correction denominator.
  grid <- tidyr::expand_grid(prior_se = se_grid, prior_sp = sp_grid)
  grid <- grid[grid$prior_se + grid$prior_sp > 1, ]
  n_grid <- nrow(grid)

  results <- vector("list", n_grid)

  for (i in seq_len(n_grid)) {
    c_se <- grid$prior_se[i]
    c_sp <- grid$prior_sp[i]

    # Clamp grid centres to the admissible range used by as_brms_prior() so
    # boundary cells (e.g. Se/Sp = 1.0) pass validation and remain fittable.
    lower <- if (bounded) 0.5 else 1e-4
    c_se  <- clamp_probability(c_se, lower, 1 - 1e-4)
    c_sp  <- clamp_probability(c_sp, lower, 1 - 1e-4)

    # Use a separate fit file for each pair. rerun controls whether brms
    # reuses that cached result or fits it again.
    # Cache path
    file_path <- NULL
    if (!is.null(model_dir)) {
      file_path <- file.path(model_dir,
                              .mcma_grid_key(c_se, c_sp))
    }

    tryCatch({
      # Update priors with grid-point centres
      updated_priors <- update(priors, se = c_se, sp = c_sp)

      fit <- mcma_fit(
        data          = data,
        priors        = updated_priors,
        prev_center   = prev_center,
        bounded       = bounded,
        prev_re       = prev_re,
        sesp_re       = sesp_re,
        correlated_re = correlated_re,
        measure_col   = measure_col,
        file          = file_path,
        file_refit    = if (rerun) "always" else "never",
        ...
      )

      # Extract prevalence and heterogeneity from the fitted model;
      # unavailable heterogeneity is recorded as missing.
      prev <- extract_prevalence(fit)
      tau  <- tryCatch(extract_tau(fit), error = function(e) {
        tibble::tibble(mean = NA_real_, median = NA_real_,
                       sd = NA_real_, l95 = NA_real_, u95 = NA_real_)
      })

      # Add the spread expected for a new study, which includes heterogeneity
      # as well as uncertainty in the pooled intercept.
      # Prediction interval
      pi_bounds <- tryCatch({
        extract_prevalence(fit, prediction_interval = TRUE)
      }, error = function(e) {
        tibble::tibble(pi_l95 = NA_real_, pi_u95 = NA_real_)
      })

      # Record chain mixing and divergent transitions alongside each estimate
      # so problematic grid fits can be identified.
      # Diagnostics
      dr <- posterior::as_draws_df(fit)
      rhat_vals <- tryCatch(
        posterior::rhat(fit),
        error = function(e) NA_real_
      )
      max_rhat <- if (is.numeric(rhat_vals)) max(rhat_vals, na.rm = TRUE) else NA_real_

      np <- tryCatch(
        brms::nuts_params(fit),
        error = function(e) NULL
      )
      ndiv <- if (!is.null(np)) {
        sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
      } else {
        NA_integer_
      }

      results[[i]] <- tibble::tibble(
        prior_se    = c_se,
        prior_sp    = c_sp,
        prev_mean   = prev$mean,
        prev_median = prev$median,
        prev_sd     = prev$sd,
        ci_lb       = prev$l95,
        ci_ub       = prev$u95,
        ci_width    = prev$u95 - prev$l95,
        pi_lb       = if ("pi_l95" %in% names(pi_bounds)) pi_bounds$pi_l95 else NA_real_,
        pi_ub       = if ("pi_u95" %in% names(pi_bounds)) pi_bounds$pi_u95 else NA_real_,
        tau_mean    = tau$mean,
        maxrhat     = max_rhat,
        ndiv        = ndiv
      )

    # Keep the grid complete after a failed fit: warn and return a labelled
    # row of missing estimates instead of aborting all remaining pairs.
    }, error = function(e) {
      rlang::warn(sprintf("Failed at Se=%.3f, Sp=%.3f: %s", c_se, c_sp, e$message))
      results[[i]] <<- tibble::tibble(
        prior_se = c_se, prior_sp = c_sp,
        prev_mean = NA_real_, prev_median = NA_real_, prev_sd = NA_real_,
        ci_lb = NA_real_, ci_ub = NA_real_, ci_width = NA_real_,
        pi_lb = NA_real_, pi_ub = NA_real_, tau_mean = NA_real_,
        maxrhat = NA_real_, ndiv = NA_integer_
      )
    })

    message(sprintf("[%d/%d] Se=%.3f, Sp=%.3f done", i, n_grid, c_se, c_sp))
  }

  # Combine the grid rows and attach the settings needed by the print and
  # heatmap methods.
  result_tbl <- do.call(rbind, results)

  out <- list(
    results     = result_tbl,
    se_grid     = se_grid,
    sp_grid     = sp_grid,
    prev_center = prev_center,
    bounded     = bounded,
    model_dir   = model_dir
  )
  class(out) <- c("mcma_sensitivity", "list")
  out
}


#' Sensitivity comparison with separate gold/screen prevalence
#'
#' Same as `mcma_sensitivity()`, but uses a comparison formula with separate
#' gold/screen prevalence intercepts at each grid point.
#'
#' @inheritParams mcma_sensitivity
#' @param gold_column Name of the gold indicator column.
#' @param study_col Column used by the default prevalence and accuracy
#'   random-effects formulas. Explicitly supplied formulas are used unchanged.
#' @param refresh Iteration interval for progress printing (0 = silent).
#' @inheritParams mcma_fit
#' @param ... Additional arguments passed to `brms::brm()` at each grid point.
#' @return An `mcma_sensitivity` object with additional gold/screen columns.
#' @export
mcma_sensitivity_comparison <- function(data,
                                        priors,
                                        se_grid,
                                        sp_grid,
                                        prev_center,
                                        bounded       = FALSE,
                                        prev_re       = ~ (1 | es_id),
                                        sesp_re       = ~ (1 | es_id),
                                        correlated_re = TRUE,
                                        measure_col   = "measure_id",
                                        gold_column   = "is_gold",
                                        study_col     = "es_id",
                                        model_dir     = "models/sensitivity_comp/",
                                        adapt_delta = 0.9999,
                                        max_treedepth = 20,
                                        iter = 6000,
                                        warmup = floor(iter / 2),
                                        refresh = 50,
                                        rerun = FALSE,
                                        ...) {

  # At each prior setting, estimate separate corrected interview and
  # screening prevalences and their posterior difference.
  if (missing(prev_re)) prev_re <- .mcma_study_re(study_col)
  if (missing(sesp_re)) sesp_re <- .mcma_study_re(study_col)

  # Subsets need priors only for the instruments they actually contain.
  priors <- .mcma_match_priors(priors, data[[measure_col]])

  if (!is.null(model_dir) && !dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }

  # Evaluate the same accuracy-prior combinations used by the ordinary
  # sensitivity grid.
  grid <- tidyr::expand_grid(prior_se = se_grid, prior_sp = sp_grid)
  grid <- grid[grid$prior_se + grid$prior_sp > 1, ]
  n_grid <- nrow(grid)

  # Ensure measure_col is character
  data[[measure_col]] <- as.character(data[[measure_col]])

  results <- vector("list", n_grid)

  for (i in seq_len(n_grid)) {
    c_se <- grid$prior_se[i]
    c_sp <- grid$prior_sp[i]

    # Clamp grid centres to the admissible range used by as_brms_prior() so
    # boundary cells (e.g. Se/Sp = 1.0) pass validation and remain fittable.
    lower <- if (bounded) 0.5 else 1e-4
    c_se  <- clamp_probability(c_se, lower, 1 - 1e-4)
    c_sp  <- clamp_probability(c_sp, lower, 1 - 1e-4)

    file_path <- NULL
    if (!is.null(model_dir)) {
      file_path <- file.path(model_dir,
                              .mcma_grid_key(c_se, c_sp))
    }

    tryCatch({
      updated_priors <- update(priors, se = c_se, sp = c_sp)

      # Give gold and screening studies distinct prevalence coefficients
      # while retaining the accuracy correction in the likelihood.
      # Build comparison formula (pi ~ 0 + is_gold)
      formula <- mcma_formula(
        bounded       = bounded,
        prev_re       = prev_re,
        sesp_re       = sesp_re,
        correlated_re = correlated_re,
        moderators    = NULL,
        measure_col   = measure_col,
        gold_column   = gold_column,
        comparison    = TRUE
      )

      # Translate the updated accuracy table into brms priors for this grid
      # point.
      brms_prior <- as_brms_prior(
        priors      = updated_priors,
        prev_center = prev_center,
        bounded     = bounded,
        has_prev_re = !is.null(prev_re),
        has_sesp_re = !is.null(sesp_re),
        measure_col = measure_col # Match prior names to the formula's column.
      )

      # Ensure is_gold is a factor
      data[[gold_column]] <- factor(data[[gold_column]])

      # Fit or reload this grid point using an identity link, since the
      # nonlinear formula already yields a probability.
      fit <- brms::brm(
        formula = formula,
        data    = data,
        family  = stats::binomial(link = "identity"),
        prior   = brms_prior,
        file    = file_path,
        file_refit = if (rerun) "always" else "never",
        refresh = refresh,
        control = list(
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth
        ),
        iter = iter,
        warmup=warmup,
        ...
      )

      # Retain the accuracy scale and priors when grid caches are read directly.
      fit <- .mcma_finalize_fit(fit, list(
        bounded = bounded, comparison = TRUE, measure_col = measure_col,
        gold_column = gold_column, study_col = study_col, prev_center = prev_center
      ), priors = updated_priors, file = file_path, ...)

      # Extract gold/screen prevalence
      dr   <- posterior::as_draws_df(fit)
      nms  <- names(dr)
      cols <- .find_gold_screen_cols(nms, gold_column = gold_column)

      # Transform each group to the probability scale, then subtract paired
      # draws to preserve their posterior dependence.
      gold_prev   <- stats::plogis(as.numeric(dr[[cols$gold]]))
      screen_prev <- stats::plogis(as.numeric(dr[[cols$screen]]))
      diff_draws  <- gold_prev - screen_prev

      results[[i]] <- tibble::tibble(
        prior_se    = c_se,
        prior_sp    = c_sp,
        gold_mean   = mean(gold_prev),
        gold_ci_lb  = stats::quantile(gold_prev, 0.025, names = FALSE),
        gold_ci_ub  = stats::quantile(gold_prev, 0.975, names = FALSE),
        screen_mean = mean(screen_prev),
        screen_ci_lb = stats::quantile(screen_prev, 0.025, names = FALSE),
        screen_ci_ub = stats::quantile(screen_prev, 0.975, names = FALSE),
        diff_mean   = mean(diff_draws),
        diff_ci_lb  = stats::quantile(diff_draws, 0.025, names = FALSE),
        diff_ci_ub  = stats::quantile(diff_draws, 0.975, names = FALSE)
      )

    # Report failed grid points as missing estimates and continue with the
    # remaining settings.
    }, error = function(e) {
      rlang::warn(sprintf("Failed at Se=%.3f, Sp=%.3f: %s", c_se, c_sp, e$message))
      results[[i]] <<- tibble::tibble(
        prior_se = c_se, prior_sp = c_sp,
        gold_mean = NA_real_, gold_ci_lb = NA_real_, gold_ci_ub = NA_real_,
        screen_mean = NA_real_, screen_ci_lb = NA_real_, screen_ci_ub = NA_real_,
        diff_mean = NA_real_, diff_ci_lb = NA_real_, diff_ci_ub = NA_real_
      )
    })

    message(sprintf("[%d/%d] Se=%.3f, Sp=%.3f done", i, n_grid, c_se, c_sp))
  }

  # Return one row per prior pair with the group estimates and their
  # difference.
  result_tbl <- do.call(rbind, results)

  out <- list(
    results     = result_tbl,
    se_grid     = se_grid,
    sp_grid     = sp_grid,
    prev_center = prev_center,
    bounded     = bounded,
    model_dir   = model_dir
  )
  class(out) <- c("mcma_sensitivity", "list")
  out
}


#' Print an mcma_sensitivity object
#' @param x An mcma_sensitivity object.
#' @param ... Ignored.
#' @export
print.mcma_sensitivity <- function(x, ...) {

  # Display the grid ranges and, when present, the range of posterior median
  # prevalence estimates.

  cat(sprintf("mcma_sensitivity: %d grid points\n", nrow(x$results)))
  cat(sprintf("  Se range: [%.3f, %.3f]\n", min(x$se_grid), max(x$se_grid)))
  cat(sprintf("  Sp range: [%.3f, %.3f]\n", min(x$sp_grid), max(x$sp_grid)))
  if (!all(is.na(x$results$prev_median))) {
    cat(sprintf("  Prevalence range: [%.4f, %.4f]\n",
                min(x$results$prev_median, na.rm = TRUE),
                max(x$results$prev_median, na.rm = TRUE)))
  }
  invisible(x)
}


#' Plot sensitivity analysis results
#'
#' @param x An mcma_sensitivity object.
#' @param type `"heatmap"` for prevalence heatmap, `"difference"` for
#'   gold-screen difference heatmap.
#' @param metric Which metric to colour the heatmap by.
#' @param as_pct Logical; display as percentages.
#' @param divergent_palette Logical; use divergent blue-red palette for heatmap.
#' @param fill_limits Optional limits for the fill scale.
#' @param text_size Size of cell label text.
#' @param digits Number of decimal places in labels.
#' @param include_cell_estimates Logical; for difference heatmap, show gold/screen estimates.
#' @param title Optional plot title (auto-generated if NULL).
#' @param highlight_cells Data frame with `prior_se` and `prior_sp` columns, or
#'   list of lists with `se` and `sp` elements.
#' @param highlight_color Colour for highlight rectangles.
#' @param highlight_size Linewidth for highlight rectangles.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.mcma_sensitivity <- function(x,
                                  type                   = c("heatmap", "difference"),
                                  metric                 = "prev_median",
                                  as_pct                 = TRUE,
                                  divergent_palette      = FALSE,
                                  fill_limits            = NULL,
                                  text_size              = 2.5,
                                  digits                 = 2,
                                  include_cell_estimates = FALSE,
                                  title                  = NULL,
                                  highlight_cells        = NULL,
                                  highlight_color        = "red",
                                  highlight_size         = 1.5,
                                  ...) {

  # Turn grid results into a prevalence or group-difference heatmap.
  # Multiplying probabilities by 100 expresses values in percent or
  # percentage points.

  type <- match.arg(type)
  dat  <- x$results
  mult <- if (as_pct) 100 else 1
  fmt  <- paste0("%.", digits, "f")

  if (type == "heatmap") {

    # Select the quantity to colour by and convert the numerical prior
    # centres to evenly spaced plotting categories.
    fill_var <- if (metric %in% names(dat)) metric else "prev_median"
    is_prev  <- fill_var %in% c("prev_mean", "prev_median")

    plot_data <- dat
    plot_data$fill_value <- plot_data[[fill_var]] * mult
    plot_data$x_fac <- factor(plot_data$prior_sp)
    plot_data$y_fac <- factor(plot_data$prior_se)

    # Prevalence cells display a point estimate, credible interval for the
    # typical study, and prediction interval for a new study.
    if (is_prev) {
      plot_data$label <- paste0(
        sprintf(fmt, plot_data$fill_value), "\n",
        "CI [", sprintf(fmt, plot_data$ci_lb * mult), ", ",
        sprintf(fmt, plot_data$ci_ub * mult), "]\n",
        "PI [", sprintf(fmt, plot_data$pi_lb * mult), ", ",
        sprintf(fmt, plot_data$pi_ub * mult), "]"
      )
    } else {
      plot_data$label <- sprintf(fmt, plot_data$fill_value)
    }

    metric_labels <- c(
      prev_mean   = "Posterior Mean Prevalence",
      prev_median = "Posterior Median Prevalence",
      ci_width    = "95% CI Width",
      pi_lb       = "95% PI Lower Bound",
      pi_ub       = "95% PI Upper Bound",
      prev_sd     = "Posterior SD"
    )

    if (is.null(title)) {
      title <- metric_labels[fill_var]
      if (is.na(title)) title <- fill_var
      if (as_pct) title <- paste0(title, " (%)")
    }

    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = .data$x_fac, y = .data$y_fac,
                                      fill = .data$fill_value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$label),
                         size = text_size, color = "black", lineheight = 0.9) +
      ggplot2::labs(x = "Prior Sp Centre", y = "Prior Se Centre",
                    fill = if (as_pct) "%" else fill_var, title = title) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text  = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

    # Choose a sequential palette for magnitudes or a two-sided palette
    # centred at the median plotted value.
    if (divergent_palette) {
      p <- p + ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = stats::median(plot_data$fill_value), limits = fill_limits
      )
    } else {
      p <- p + ggplot2::scale_fill_gradient(
        low = "#F7FBFF", high = "#08519C", limits = fill_limits
      )
    }

  } else {

    # For group comparisons, centre the colour scale at zero so direction and
    # size of the difference are visible.
    # Difference heatmap
    if (!"diff_mean" %in% names(dat)) {
      rlang::abort("Difference heatmap requires gold/screen comparison results.")
    }

    if (is.null(title)) {
      title <- "Gold Standard \u2013 Corrected Screening Prevalence"
      if (as_pct) title <- paste0(title, " (%)")
    }

    plot_data <- dat
    plot_data$fill_value <- plot_data$diff_mean * mult
    plot_data$x_fac <- factor(plot_data$prior_sp)
    plot_data$y_fac <- factor(plot_data$prior_se)

    # Optionally show both group estimates in addition to their difference
    # and credible limits.
    if (include_cell_estimates) {
      plot_data$label <- paste0(
        "Gold: ", sprintf(fmt, plot_data$gold_mean * mult), "\n",
        "Scr:  ", sprintf(fmt, plot_data$screen_mean * mult), "\n",
        "\u0394: ", sprintf(fmt, plot_data$diff_mean * mult), "\n",
        "[", sprintf(fmt, plot_data$diff_ci_lb * mult), ", ",
        sprintf(fmt, plot_data$diff_ci_ub * mult), "]"
      )
    } else {
      plot_data$label <- paste0(
        sprintf(fmt, plot_data$diff_mean * mult), "\n",
        "[", sprintf(fmt, plot_data$diff_ci_lb * mult), ", ",
        sprintf(fmt, plot_data$diff_ci_ub * mult), "]"
      )
    }

    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = .data$x_fac, y = .data$y_fac,
                                      fill = .data$fill_value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$label),
                         size = text_size, color = "black", lineheight = 0.9) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = fill_limits
      ) +
      ggplot2::labs(x = "Prior Sp Centre", y = "Prior Se Centre",
                    fill = if (as_pct) "%" else "Diff", title = title) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text  = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }

  # Add outlines for selected prior combinations after constructing either
  # kind of heatmap.
  # Highlight cells
  p <- .add_highlight_cells(p, plot_data, highlight_cells,
                            highlight_color, highlight_size)
  p
}


# Internal helper for highlight rectangles
.add_highlight_cells <- function(p, plot_data, highlight_cells,
                                 highlight_color = "red",
                                 highlight_size  = 1.5) {

  # Add outlines around requested prior combinations. Work with the plotted
  # factor positions so boxes align with the centres of the heatmap cells.
  if (is.null(highlight_cells)) return(p)

  if (is.list(highlight_cells) && !is.data.frame(highlight_cells)) {
    highlight_cells <- dplyr::bind_rows(lapply(highlight_cells, function(x) {
      tibble::tibble(prior_se = x$se, prior_sp = x$sp)
    }))
  }

  # Read the axis category order rather than assuming the grid values start
  # at a particular coordinate.
  sp_levels <- levels(plot_data$x_fac)
  se_levels <- levels(plot_data$y_fac)

  highlight_cells <- highlight_cells %>%
    dplyr::mutate(
      x_num = match(as.character(.data$prior_sp), sp_levels),
      y_num = match(as.character(.data$prior_se), se_levels)
    ) %>%
    dplyr::filter(!is.na(.data$x_num), !is.na(.data$y_num))

  # Draw a one-cell-wide rectangle around each matched combination;
  # combinations absent from the grid have already been removed.
  if (nrow(highlight_cells) > 0) {
    for (r in seq_len(nrow(highlight_cells))) {
      p <- p + ggplot2::annotate(
        "rect",
        xmin = highlight_cells$x_num[r] - 0.5,
        xmax = highlight_cells$x_num[r] + 0.5,
        ymin = highlight_cells$y_num[r] - 0.5,
        ymax = highlight_cells$y_num[r] + 0.5,
        fill = NA, color = highlight_color, linewidth = highlight_size
      )
    }
  }

  p
}


.mcma_grid_key <- function(se, sp) {

  # Keep useful grid precision so nearby prior settings get distinct caches.
  # nsmall preserves the familiar two-decimal filenames for existing grids.
  label <- function(x) format(x, digits = 15, nsmall = 2,
                              scientific = FALSE, trim = TRUE)
  paste0("Se", label(se), "_Sp", label(sp))
}
