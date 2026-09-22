# Kernel-weighted grid averaging (KWGA) for Se/Sp estimation
# ================================================

#' Discrete kernel-weighted grid averaging (KWGA) over a grid of Se/Sp candidates
#'
#' Scores each (Se, Sp) pair by the agreement between the Rogan-Gladen-
#' corrected screening posterior and the gold-standard posterior from the
#' comparison model. Uses a draw-paired Gaussian kernel.
#'
#' By default the Rogan-Gladen-corrected screening prevalence is clamped to the
#' unit interval before it is compared with the gold-standard draws, so an
#' accuracy pair that implies a negative prevalence is scored as if it implied
#' zero. With `clamp_scoring = FALSE` the unclamped value is scored instead, so
#' such pairs are penalised by their full distance from the gold-standard
#' draws. Pair that setting with joint resampling in [mcma_kwga_prevalence()]
#' (its default `resample = "auto"` does so automatically).
#'
#' @param fit_comparison A brmsfit from `mcma_fit_comparison()`, or a list
#'   with `$gold_draws` and `$screen_draws` on the probability scale.
#' @param se_grid Candidate Se values.
#' @param sp_grid Candidate Sp values.
#' @param bandwidth Kernel bandwidth. Defaults to the SD of gold draws.
#' @param prior_weights Optional prior weights for the grid (same length as
#'   the number of valid grid points, or NULL for uniform).
#' @param gold_column Name of the gold indicator column.
#' @param clamp_scoring Logical. `TRUE` (default) scores the clamped
#'   Rogan-Gladen correction (original construction); `FALSE` scores the
#'   unclamped correction. The per-grid-point summary columns are always
#'   computed from the clamped values.
#' @return An S3 object of class `mcma_kwga`.
#' @export
mcma_kwga <- function(fit_comparison,
                     se_grid     = seq(0.60, 0.95, by = 0.025),
                     sp_grid     = seq(0.60, 0.95, by = 0.025),
                     bandwidth   = NULL,
                     prior_weights = NULL,
                     gold_column = "is_gold",
                     clamp_scoring = TRUE) {

  # Compare corrected screening prevalence with interview prevalence over
  # candidate accuracy pairs. Agreement determines the relative kernel
  # weights.

  # --- Extract draws ---
  if (inherits(fit_comparison, "brmsfit")) {
    dr  <- posterior::as_draws_df(fit_comparison)
    nms <- names(dr)
    cols <- .find_gold_screen_cols(nms, gold_column = gold_column)
    gold_draws   <- stats::plogis(as.numeric(dr[[cols$gold]]))
    screen_draws <- stats::plogis(as.numeric(dr[[cols$screen]]))
  } else if (is.list(fit_comparison)) {
    gold_draws   <- as.numeric(fit_comparison$gold_draws)
    screen_draws <- as.numeric(fit_comparison$screen_draws)
  } else {
    rlang::abort("fit_comparison must be a brmsfit or a list with $gold_draws and $screen_draws.")
  }

  n_draws <- length(gold_draws)
  stopifnot(length(screen_draws) == n_draws, n_draws > 0)

  # The kernel bandwidth controls how strongly a mismatch is penalized. By
  # default it follows uncertainty in the interview prevalence.
  # --- Bandwidth ---
  if (is.null(bandwidth)) {
    bandwidth <- stats::sd(gold_draws)
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 0.02
  }

  # Evaluate every requested Se/Sp combination whose Rogan-Gladen denominator
  # is positive.
  # --- Build grid ---
  grid <- tidyr::expand_grid(se = se_grid, sp = sp_grid)
  grid <- grid[grid$se + grid$sp > 1, ]
  n_grid <- nrow(grid)
  if (n_grid == 0) rlang::abort("No valid (se, sp) pairs with se + sp > 1.")

  # Allow initial preferences over accuracy pairs; equal preferences leave
  # the agreement scores to determine relative weights.
  # --- Prior weights ---
  if (is.null(prior_weights)) {
    log_prior <- rep(0, n_grid)
  } else {
    if (length(prior_weights) != n_grid) {
      rlang::warn("prior_weights length mismatch; using uniform.")
      log_prior <- rep(0, n_grid)
    } else {
      pw <- prior_weights / sum(prior_weights)
      log_prior <- log(pw)
    }
  }

  # Arrange posterior draws down rows and candidate accuracy pairs across
  # columns, allowing every correction to be computed together.
  # --- Vectorized scoring ---
  # Matrices: rows = draws, cols = grid points
  screen_mat <- matrix(screen_draws, nrow = n_draws, ncol = n_grid)
  sp_minus1  <- matrix(grid$sp - 1,  nrow = n_draws, ncol = n_grid, byrow = TRUE)
  denom_mat  <- matrix(grid$se + grid$sp - 1, nrow = n_draws, ncol = n_grid, byrow = TRUE)
  gold_mat   <- matrix(gold_draws, nrow = n_draws, ncol = n_grid)

  # Rogan-Gladen correction per draw per grid point
  corrected_raw <- (screen_mat + sp_minus1) / denom_mat
  corrected_raw <- matrix(corrected_raw, nrow = n_draws, ncol = n_grid)
  corrected <- pmin(1, pmax(0, corrected_raw))
  corrected <- matrix(corrected, nrow = n_draws, ncol = n_grid)

  # Small draw-paired differences receive high normal-kernel scores. These
  # scores measure agreement, not a fitted binomial likelihood.
  # Discrepancy and scoring
  # With clamp_scoring = TRUE a pair implying a negative prevalence is scored as
  # zero prevalence; with FALSE it is scored at its full (negative) distance.
  discrepancy  <- gold_mat - (if (isTRUE(clamp_scoring)) corrected else corrected_raw)
  discrepancy  <- matrix(discrepancy, nrow = n_draws, ncol = n_grid)
  log_lik_mat  <- stats::dnorm(discrepancy, mean = 0, sd = bandwidth, log = TRUE)
  log_lik_mat  <- matrix(log_lik_mat, nrow = n_draws, ncol = n_grid)

  # Subtract the largest log score before exponentiating, preventing
  # numerical underflow when agreement is poor.
  # Stable log-mean-exp per column
  max_ll <- vapply(seq_len(n_grid), function(k) max(log_lik_mat[, k]), numeric(1))
  log_ml <- vapply(seq_len(n_grid), function(k) {
    m <- max_ll[k]
    m + log(mean(exp(log_lik_mat[, k] - m)))
  }, numeric(1))

  # Summary stats per grid point
  grid$corrected_mean <- vapply(seq_len(n_grid), function(k) mean(corrected[, k]), numeric(1))
  grid$corrected_median <- vapply(seq_len(n_grid), function(k) stats::median(corrected[, k]), numeric(1))
  grid$corrected_ci_lb <- vapply(seq_len(n_grid), function(k)
    stats::quantile(corrected[, k], 0.025, names = FALSE), numeric(1))
  grid$corrected_ci_ub <- vapply(seq_len(n_grid), function(k)
    stats::quantile(corrected[, k], 0.975, names = FALSE), numeric(1))

  # Combine agreement scores with initial grid weights and normalize to a
  # total of one; sort highest-weight pairs first.
  # --- Posterior weights ---
  grid$log_ml <- log_ml
  grid$log_prior <- log_prior
  log_post <- log_ml + log_prior
  max_lp <- max(log_post)
  grid$weight <- exp(log_post - max_lp)
  grid$weight <- grid$weight / sum(grid$weight)

  grid <- grid[order(-grid$weight), ]
  grid$cumweight <- cumsum(grid$weight)

  # Sum over the other accuracy parameter to obtain separate Se and Sp weight
  # distributions.
  # --- Marginals ---
  se_marg <- stats::aggregate(weight ~ se, data = grid, FUN = sum)
  sp_marg <- stats::aggregate(weight ~ sp, data = grid, FUN = sum)
  names(se_marg)[2] <- "w"
  names(sp_marg)[2] <- "w"

  # Keep several summaries because the highest-weight pair, marginal modes,
  # and weighted averages need not coincide.
  # --- Point estimates ---
  joint_map <- grid[1, ]
  se_map_marg <- se_marg$se[which.max(se_marg$w)]
  sp_map_marg <- sp_marg$sp[which.max(sp_marg$w)]
  se_mean <- sum(grid$se * grid$weight)
  sp_mean <- sum(grid$sp * grid$weight)
  se_median <- .weighted_quantile(se_marg$se, se_marg$w, 0.5)
  sp_median <- .weighted_quantile(sp_marg$sp, sp_marg$w, 0.5)

  estimates <- tibble::tibble(
    method = c("joint_map", "marginal_map", "weighted_mean", "weighted_median"),
    se = c(joint_map$se, se_map_marg, se_mean, se_median),
    sp = c(joint_map$sp, sp_map_marg, sp_mean, sp_median)
  )

  # Return the scores, summaries, and original paired draws together for
  # later prevalence correction and plotting.
  out <- list(
    grid         = tibble::as_tibble(grid),
    estimates    = estimates,
    gold_draws   = gold_draws,
    screen_draws = screen_draws,
    bandwidth    = bandwidth,
    gold_column  = gold_column,
    clamp_scoring = isTRUE(clamp_scoring)
  )
  class(out) <- c("mcma_kwga", "list")
  out
}


#' Kernel-weighted grid-averaged prevalence estimate
#'
#' Computes a kernel-weighted grid-averaged prevalence by sampling Se/Sp pairs according to
#' their posterior weights and applying the Rogan-Gladen correction.
#'
#' @param kwga An `mcma_kwga` object from `mcma_kwga()`.
#' @param fit_comparison Optional comparison fit, or a list with `screen_draws`
#'   (and, for joint resampling, `gold_draws`) on the probability scale. When
#'   supplied, its draws replace those stored in `kwga`; NULL uses the stored
#'   draws. Grid weights remain those in `kwga`, so normally supply the same fit
#'   used to construct that object.
#' @param mode `"analytic"` applies the Rogan-Gladen correction post hoc.
#' @param resample `"independent"` (original construction) samples an accuracy
#'   pair by its grid weight and, independently, a screening draw.
#'   `"joint"` samples the (posterior draw, accuracy pair) combination with
#'   probability proportional to prior weight times kernel score, so every
#'   prevalence draw is conditioned on its agreement with the gold-standard
#'   draws. `"auto"` (default) uses `"joint"` when `kwga` was built with
#'   `clamp_scoring = FALSE` and `"independent"` otherwise.
#' @param n_draws Number of Monte Carlo draws.
#' @param seed Optional integer seed for reproducible Monte Carlo draws; the
#'   global RNG state is restored on exit. Defaults to NULL (unseeded).
#' @return A list with `$summary`, `$draws`, `$weights`, and `$resample` (the
#'   scheme actually used).
#' @export
mcma_kwga_prevalence <- function(kwga,
                                fit_comparison = NULL,
                                mode     = c("analytic"),
                                n_draws  = 4000,
                                seed     = NULL,
                                resample = c("auto", "independent", "joint")) {

  # Generate corrected prevalence draws from the grid-weighted mixture. Each
  # draw combines one sampled accuracy pair with a screening prevalence draw.

  mode <- match.arg(mode)
  resample <- match.arg(resample)

  # Use an optional local reproducibility seed and restore the previously
  # existing global random-number state on exit.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  grid <- kwga$grid

  if (mode == "analytic") {
    screen_draws <- kwga$screen_draws
    gold_draws   <- kwga$gold_draws

    if (!is.null(fit_comparison)) {
      # Use the supplied comparison rather than silently retaining old draws.
      if (inherits(fit_comparison, "brmsfit")) {
        cfg <- attr(fit_comparison, "mcma_config")
        gold_column <- cfg$gold_column
        if (is.null(gold_column)) gold_column <- kwga$gold_column
        if (is.null(gold_column)) gold_column <- "is_gold"
        gs <- extract_gold_screen_diff(
          fit_comparison, summary = FALSE, gold_column = gold_column
        )
        screen_draws <- gs$screen
        gold_draws   <- gs$gold
      } else if (is.list(fit_comparison)) {
        screen_draws <- fit_comparison$screen_draws
        gold_draws   <- fit_comparison$gold_draws
      } else {
        stop("fit_comparison must be a brmsfit or a list with screen_draws.", call. = FALSE)
      }
    }

    .mcma_validate_probability(screen_draws, "screen_draws", length(screen_draws))

    # "auto" pairs joint resampling with unclamped scoring and the original
    # independent scheme with clamped scoring.
    if (resample == "auto") {
      resample <- if (isFALSE(kwga$clamp_scoring)) "joint" else "independent"
    }

    if (resample == "independent") {
      # Sample grid rows by weight
      idx <- sample(
        seq_len(nrow(grid)),
        size    = n_draws,
        replace = TRUE,
        prob    = grid$weight
      )
      sampled <- grid[idx, ]

      # For each sampled (Se, Sp), pick a random screening draw and correct
      draw_idx <- sample(seq_along(screen_draws), size = n_draws, replace = TRUE)
      p_screen <- screen_draws[draw_idx]
      corrected <- (p_screen + sampled$sp - 1) / (sampled$se + sampled$sp - 1)
      corrected <- pmin(1, pmax(0, corrected))
    } else {
      # Joint resampling: pick (posterior draw, accuracy pair) together with
      # probability proportional to prior weight x kernel score, scored the same
      # way the grid weights were (clamped or unclamped), so each prevalence
      # draw is conditioned on its agreement with the gold-standard draws.
      if (is.null(gold_draws) || length(gold_draws) != length(screen_draws)) {
        stop("Joint resampling needs gold_draws paired with screen_draws (same length).",
             call. = FALSE)
      }
      .mcma_validate_probability(gold_draws, "gold_draws", length(gold_draws))
      n_post <- length(screen_draws)
      n_grid <- nrow(grid)
      log_prior <- if (is.null(grid$log_prior)) rep(0, n_grid) else grid$log_prior
      raw <- outer(screen_draws, grid$sp - 1, "+") /
        matrix(grid$se + grid$sp - 1, nrow = n_post, ncol = n_grid, byrow = TRUE)
      scored <- if (isFALSE(kwga$clamp_scoring)) raw else pmin(1, pmax(0, raw))
      score <- stats::dnorm(gold_draws - scored, mean = 0, sd = kwga$bandwidth, log = TRUE) +
        matrix(log_prior, nrow = n_post, ncol = n_grid, byrow = TRUE)
      idx <- sample.int(n_post * n_grid, size = n_draws, replace = TRUE,
                        prob = exp(score - max(score)))
      corrected <- pmin(1, pmax(0, raw[idx]))
    }

    # Measure weight concentration: equal weights over K pairs give K
    # effective pairs; one dominant pair gives about one.
    n_eff <- 1 / sum(grid$weight^2)

    # Summarise the simulated mixture, including uncertainty from both the
    # sampled accuracy pairs and screening prevalence draws.
    out_summary <- tibble::tibble(
      mean    = mean(corrected),
      median  = stats::median(corrected),
      ci_lb   = stats::quantile(corrected, 0.025, names = FALSE),
      ci_ub   = stats::quantile(corrected, 0.975, names = FALSE),
      sd      = stats::sd(corrected),
      n_effective_models = n_eff
    )

    return(list(
      summary = out_summary,
      draws   = corrected,
      weights = grid,
      resample = resample
    ))
  }
}


#' Print an mcma_kwga object
#' @param x An mcma_kwga object.
#' @param ... Ignored.
#' @export
print.mcma_kwga <- function(x, ...) {

  # Display the grid size, agreement bandwidth, highest-weight pair, and
  # alternative summaries without modifying the stored results.

  cat("Kernel-weighted grid averaging (KWGA) over the Se/Sp grid\n")
  cat(sprintf("  Grid size: %d valid points\n", nrow(x$grid)))
  cat(sprintf("  Bandwidth: %.4f\n", x$bandwidth))
  cat(sprintf("  Scoring: %s\n",
              if (isFALSE(x$clamp_scoring)) "unclamped correction" else "clamped correction"))
  cat(sprintf("  Top weight: %.4f at Se=%.3f, Sp=%.3f\n",
              x$grid$weight[1], x$grid$se[1], x$grid$sp[1]))
  cat("\nEstimates:\n")
  print(x$estimates, n = nrow(x$estimates))
  invisible(x)
}


#' Summary for mcma_kwga
#' @param object An mcma_kwga object.
#' @param ... Ignored.
#' @export
summary.mcma_kwga <- function(object, ...) {

  # Reuse the same display when the user calls summary() on a grid-averaging
  # result.
  print(object)
}


# --- Internal helpers ---

.weighted_quantile <- function(x, w, probs = 0.5) {

  # Sort values and accumulate their weights. Return the first value reaching
  # the requested cumulative probability; no interpolation is used.

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= probs)[1]]
}
