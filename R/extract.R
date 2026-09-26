# Posterior extraction functions
# ==============================

#' Extract or predict prevalence
#'
#' Extracts pooled prevalence or predicts latent prevalence for specified data
#' rows and random effects. Handles naive models
#' (`b_Intercept`), corrected models (`b_pi_Intercept`), and comparison
#' models (separate gold/screen intercepts).
#'
#' @param fit A brmsfit object.
#' @param summary Logical; if TRUE, returns a summary tibble. If FALSE,
#'   returns raw posterior draws on the scale selected by `transform`.
#' @param prediction_interval Logical; if TRUE, adds `pi_l95` and `pi_u95`.
#'   For a brms-style prediction request, these are the same target-draw
#'   quantiles as `l95` and `u95`. Otherwise, the original new-study PI is used.
#' @param probs Quantile probabilities for CI/PI bounds.
#' @param transform Logical; if TRUE, returns probability scale.
#' @param prediction_groups Legacy argument: names of grouping factors whose levels are new.
#'   NULL includes every prevalence grouping factor not in `condition_on`
#'   or `zero_groups`.
#' @param condition_on Legacy argument: named list of existing group levels, for example
#'   `list(refid = "Smith2020")`. Their fitted prevalence deviations are
#'   retained in predictions, including posterior uncertainty. Every grouping
#'   factor must be new, conditioned on, or explicitly held at zero when
#'   requesting a PI.
#' @param gold_column Gold indicator column; NULL uses the saved configuration
#'   or `"is_gold"` for older fits.
#' @param zero_groups Legacy argument: names of grouping factors whose random effects are fixed
#'   at zero, for example `"refid"`. These contribute neither a fitted deviation
#'   nor heterogeneity to the prediction interval. NULL holds no groups at zero.
#' @param newdata Optional data frame for prediction. Explicitly supplying
#'   `newdata` (including `NULL`), `re_formula`, or `re.form` selects brms-style
#'   prediction. `NULL` uses the original model data. Include all columns
#'   required by the fitted brms model, including `n` and `measure_id` or
#'   `acc_group` where applicable; copying rows of `fit$data` is a convenient
#'   starting point. Set moderators to the values of interest. Existing group
#'   labels use fitted effects; unseen labels require `allow_new_levels = TRUE`.
#' @param re_formula Formula selecting group-level terms, as in brms. `NULL`
#'   includes all terms, `NA` or `~0` excludes all, and `~(1 | es_id)` includes
#'   study effects while holding omitted report effects at zero. Excluding an
#'   effect does not average over its heterogeneity.
#' @param re.form Alias of `re_formula`. Supply only one of these arguments.
#' @param allow_new_levels Logical; allow unseen group labels in `newdata`.
#'   Defaults to `FALSE`, as in brms.
#' @param sample_new_levels How brms generates new group effects. Defaults to
#'   `"gaussian"`: draw from the fitted normal distribution of group effects.
#'   Alternatives are `"uncertainty"` (select an existing group's effect at
#'   each posterior draw) and `"old_levels"` (use one existing group's full
#'   set of draws per new level). These alternatives require a brms-style
#'   prediction request.
#' @return A summary tibble or numeric vector of draws. Comparison models
#'   return one summary row per `population`, or a list with `gold` and `screen`
#'   draw vectors. Population summaries remain at reference moderator values;
#'   `condition_on` and `zero_groups` affect only the legacy prediction interval.
#'   A brms-style request instead returns one summary row per data row, identified
#'   by `.row`, or a draws-by-row matrix when `summary = FALSE`. All summaries
#'   and draws describe the requested moderators and random effects.
#' @details With all groups new, random-intercept variances are added within
#'   each posterior draw. For a new estimate within an existing report, use
#'   `prediction_groups = "es_id"` and `condition_on = list(refid = "reportID")`.
#'   To include study heterogeneity while holding the report effect at zero,
#'   use `prediction_groups = "es_id"` and `zero_groups = "refid"`. This fixes
#'   the report deviation at its population centre on the log-odds scale; it
#'   does not average over report heterogeneity. Intercept and selected SD
#'   uncertainty remain in the interval. Prediction intervals describe latent prevalence, 
#'   without binomial noise.
#' @examples
#' \dontrun{
#' nd <- fit$data[1, , drop = FALSE]
#' nd$es_id <- "new_study"
#' extract_prevalence(fit, newdata = nd, re_formula = ~(1 | es_id),
#'                    allow_new_levels = TRUE, sample_new_levels = "gaussian",
#'                    prediction_interval = TRUE)
#' }
#' @export
extract_prevalence <- function(fit,
                               summary             = TRUE,
                               prediction_interval = FALSE,
                               probs               = c(0.025, 0.975),
                               transform           = TRUE,
                               prediction_groups   = NULL,
                               condition_on        = NULL,
                               gold_column         = NULL,
                               zero_groups         = NULL,
                               newdata             = NULL,
                               re_formula          = NULL,
                               re.form             = NULL,
                               allow_new_levels    = FALSE,
                               sample_new_levels   = "gaussian") {

  # Explicit prediction inputs select row-wise predictions.
  predict_rows <- !missing(newdata) || !missing(re_formula) || !missing(re.form)
  if (!missing(re.form)) {
    if (!missing(re_formula)) stop("Supply only one of re_formula and re.form.", call. = FALSE)
    re_formula <- re.form
  }
  sample_new_levels <- match.arg(sample_new_levels, c("gaussian", "uncertainty", "old_levels"))

  if (predict_rows) {
    if (!is.null(prediction_groups) || !is.null(condition_on) || !is.null(zero_groups)) {
      stop("Use newdata/re_formula or the legacy grouping arguments, not both.", call. = FALSE)
    }

    # Extract true prevalence before the observation-error equation.
    # brms supplies moderator effects and existing or newly sampled group effects.
    nlpar <- if ("pi" %in% names(fit$formula$pforms)) "pi" else NULL
    draws <- .mcma_prediction_draws(fit, newdata, re_formula, allow_new_levels,
                                    sample_new_levels, nlpars = nlpar)[[1]]
    if (transform) draws <- stats::plogis(draws)
    if (!summary) return(draws)
    return(.mcma_prediction_summary(draws, probs, prediction_interval))
  }

  # The legacy PI samples Gaussian heterogeneity directly; a different
  # sampling method must be requested through the brms prediction interface.
  if (sample_new_levels != "gaussian") {
    stop("Supply newdata or re_formula to use sample_new_levels other than 'gaussian'.", call. = FALSE)
  }

  # Summarise the population intercept on the probability scale by default,
  # or on the log-odds scale if requested. It represents the typical study at
  # the reference moderator values.

  # Keep paired posterior draws together, including comparison groups.
  dr <- posterior::as_draws_df(fit)
  cfg <- attr(fit, "mcma_config")
  if (is.null(gold_column)) {
    gold_column <- if (!is.null(cfg$gold_column)) cfg$gold_column else "is_gold"
  }

  if ("b_pi_Intercept" %in% names(dr)) {
    columns <- c(prevalence = "b_pi_Intercept")
  } else if ("b_Intercept" %in% names(dr)) {
    columns <- c(prevalence = "b_Intercept")
  } else {
    columns <- unlist(.find_gold_screen_cols(names(dr), gold_column))
  }

  on_scale <- if (transform) stats::plogis else identity
  draws <- lapply(columns, function(column) on_scale(as.numeric(dr[[column]])))
  if (!summary) return(if (length(draws) == 1L) draws[[1]] else draws)

  # Existing reports contribute their fitted deviations; new levels
  # receive fresh random deviations. Groups explicitly held at zero contribute
  # neither. Intercept and included-effect uncertainty are still propagated.
  if (prediction_interval) {
    target <- .mcma_prediction_target(dr, prediction_groups, condition_on, zero_groups)
  }

  rows <- lapply(names(columns), function(population) {
    out <- .mcma_draw_summary(draws[[population]], probs)

    if (prediction_interval) {
      pred_logit <- stats::rnorm(nrow(dr),
        mean = as.numeric(dr[[columns[[population]]]]) + target$offset,
        sd = target$tau
      )
      pred <- on_scale(pred_logit)
      out$pi_l95 <- stats::quantile(pred, probs[1], names = FALSE)
      out$pi_u95 <- stats::quantile(pred, probs[2], names = FALSE)
    }

    if (length(columns) > 1L) out <- tibble::add_column(out, population = population, .before = 1)
    out
  })

  do.call(rbind, rows)
}


#' Extract between-study heterogeneity SD (tau)
#'
#' @param fit A brmsfit object.
#' @param summary Logical; if TRUE, returns summary tibble.
#' @param probs Quantile probabilities.
#' @param groups Grouping factors to include, for example `"es_id"` or
#'   `c("es_id", "refid")`. NULL includes all prevalence random-intercept
#'   variances. Accuracy SDs are never included. With no random intercepts,
#'   the combined SD is zero.
#' @param re_formula Formula selecting prevalence random intercepts. `NULL`
#'   includes all, `NA` or `~0` includes none, and `~(1 | es_id)` selects study
#'   heterogeneity. This function summarises SDs, not existing group deviations.
#'   Random slopes are not supported by this combined-intercept SD.
#' @param re.form Alias of `re_formula`. Supply only one of these arguments.
#' @return A tibble or numeric vector.
#' @details `groups` remains available for compatibility. Do not combine a
#'   non-NULL `groups` with `re_formula` or `re.form`.
#' @export
extract_tau <- function(fit, summary = TRUE, probs = c(0.025, 0.975),
                        groups = NULL, re_formula = NULL, re.form = NULL) {

  # Translate brms-style formulas into the existing variance selection.
  select_formula <- !missing(re_formula) || !missing(re.form)
  if (!missing(re.form)) {
    if (!missing(re_formula)) stop("Supply only one of re_formula and re.form.", call. = FALSE)
    re_formula <- re.form
  }
  if (select_formula) {
    if (!is.null(groups)) stop("Supply groups or re_formula, not both.", call. = FALSE)
    groups <- .mcma_re_groups(re_formula)
  }

  # Add selected variances within each posterior draw before taking the
  # square root. Selecting one group returns that component alone.
  draws <- .tau_draws(fit, groups = groups)
  if (!summary) return(as.numeric(draws))

  # Summarise posterior uncertainty in the heterogeneity parameter itself.
  tibble::tibble(
    mean   = mean(draws),
    median = stats::median(draws),
    sd     = stats::sd(draws),
    l95    = stats::quantile(draws, probs[1], names = FALSE),
    u95    = stats::quantile(draws, probs[2], names = FALSE)
  )
}


#' Extract Se/Sp posterior summaries
#'
#' Extracts sensitivity and specificity posteriors from corrected or joint
#' models.
#'
#' @param fit A brmsfit object (corrected or joint model).
#' @param bounded Logical; if NULL, auto-detected from the fit.
#' @param level `"global"` returns a single estimate (joint models);
#'   `"per_measure"` returns one row per screening measure. When omitted,
#'   joint models use `"global"` and other models use `"per_measure"`.
#' @param summary Logical; if TRUE, returns summary tibble.
#' @param probs Quantile probabilities.
#' @param measure_col Name of the instrument identifier column used when
#'   fitting a per-measure model. This selects the matching Se/Sp coefficients;
#'   it is not used for a joint model's global accuracy summary.
#' @inheritParams extract_prevalence
#' @return A tibble or list of draw vectors. With explicit `newdata`,
#'   `re_formula`, or `re.form`, returns one summary row per data row, identified
#'   by `.row`, or a list containing `se` and `sp` draws-by-row matrices.
#' @details Without prediction inputs, `level` selects the original coefficient
#'   summaries. With prediction inputs, rows determine the accuracy target and
#'   `level` is not used. These predictions describe baseline Se/Sp, before any
#'   additional c/o adjustment. For joint models, `acc_group` identifies gold
#'   or screening rows. Both accuracy parameters use the same sampled group
#'   effects so their posterior dependence is retained.
#' @export
extract_sesp <- function(fit,
                         bounded = NULL,
                         level   = c("per_measure", "global"),
                         summary = TRUE,
                         probs   = c(0.025, 0.975),
                         measure_col = "measure_id",
                         newdata = NULL, re_formula = NULL, re.form = NULL, 
                         allow_new_levels = FALSE,
                         sample_new_levels = "gaussian") {

  # Keep the original coefficient summaries unless prediction inputs
  # explicitly request accuracy for particular rows and random effects.
  predict_rows <- !missing(newdata) || !missing(re_formula) || !missing(re.form)
  if (!missing(re.form)) {
    if (!missing(re_formula)) stop("Supply only one of re_formula and re.form.", call. = FALSE)
    re_formula <- re.form
  }
  sample_new_levels <- match.arg(sample_new_levels, c("gaussian", "uncertainty", "old_levels"))
  if (!predict_rows && sample_new_levels != "gaussian") {
    stop("Supply newdata or re_formula to use sample_new_levels other than 'gaussian'.", call. = FALSE)
  }

  # Select the accuracy coefficients and transform them back to
  # probabilities. These are baseline accuracy estimates, before any
  # additional c/o adjustment.
  auto_level <- missing(level)
  level <- match.arg(level)

  # Read the fit metadata before transforming accuracy coefficients; the
  # bounded and unbounded models use different inverse maps.
  # Auto-detect bounded
  if (is.null(bounded)) {
    cfg <- attr(fit, "mcma_config")
    if (!is.null(cfg)) {
      bounded <- isTRUE(cfg$bounded)
    } else {
      rlang::warn("mcma_config attribute not found; assuming unbounded. Set bounded = TRUE/FALSE explicitly.")
      bounded <- FALSE
    }
  }

  # Prepare Se and Sp together to retain their correlated new effects.
  # Back-transform using the accuracy scale saved with the fitted model.
  if (predict_rows) {
    draws <- .mcma_prediction_draws(fit, newdata, re_formula, allow_new_levels,
                                    sample_new_levels, nlpars = c("Se", "Sp"))
    se <- .backtransform_sesp(draws$Se, bounded)
    sp <- .backtransform_sesp(draws$Sp, bounded)
    if (!summary) return(list(se = se, sp = sp))

    rows <- lapply(seq_len(ncol(se)), function(i) {
      out <- .summarise_sesp_pair(as.character(i), se[, i], sp[, i], probs)
      out$measure_id <- NULL
      tibble::add_column(out, .row = i, .before = 1)
    })
    return(do.call(rbind, rows))
  }

  # Coefficient draws are needed only by the original summary route.
  dr <- posterior::as_draws_df(fit)
  nms <- names(dr)
  se_prefix <- paste0("^b_Se_", measure_col)
  sp_prefix <- paste0("^b_Sp_", measure_col)

  # Detect model type: joint (acc_group) vs per-measure (measure_id)
  is_joint <- any(grepl("Se_acc_group", nms))
  if (auto_level && is_joint) level <- "global"

  if (level == "global") {
    if (!is_joint) {
      rlang::warn(paste0(
        "Per-measure model detected. 'level = \"global\"' is only meaningful ",
        "for joint models. Falling back to 'per_measure'."
      ))
      level <- "per_measure"
    }
  }

  # A joint model has one screening coefficient for each accuracy parameter,
  # so return a single shared Se/Sp summary.
  if (level == "global" && is_joint) {
    se_col <- grep("Se_acc_groupscreen", nms, value = TRUE)[1]
    sp_col <- grep("Sp_acc_groupscreen", nms, value = TRUE)[1]
    if (is.na(se_col) || is.na(sp_col)) {
      rlang::abort("Could not find Se/Sp screening coefficients in joint model.")
    }
    se_raw <- as.numeric(dr[[se_col]])
    sp_raw <- as.numeric(dr[[sp_col]])
    se_draws <- .backtransform_sesp(se_raw, bounded)
    sp_draws <- .backtransform_sesp(sp_raw, bounded)

    if (!summary) {
      return(list(se = se_draws, sp = sp_draws))
    }
    return(.summarise_sesp_pair("global", se_draws, sp_draws, probs))
  }

  # For instrument-specific models, locate coefficients by their name prefix
  # and omit constants used for gold-standard accuracy.
  # Per-measure extraction
  se_cols <- grep(se_prefix, nms, value = TRUE)
  sp_cols <- grep(sp_prefix, nms, value = TRUE)

  # Filter out gold-standard (constant) columns by checking variance
  se_cols <- se_cols[vapply(se_cols, function(c) stats::sd(dr[[c]]) > 1e-6, logical(1))]
  sp_cols <- sp_cols[vapply(sp_cols, function(c) stats::sd(dr[[c]]) > 1e-6, logical(1))]

  if (length(se_cols) == 0) {
    rlang::abort("No Se measure coefficients found in the model draws.")
  }

  # Pair Se and Sp by the instrument label, rather than relying on their
  # column order.
  # Extract measure IDs from column names
  se_ids <- sub(se_prefix, "", se_cols)
  sp_ids <- sub(sp_prefix, "", sp_cols)
  measure_ids <- intersect(se_ids, sp_ids)

  # Translate brms-safe coefficient names back to the user's instrument labels.
  measure_labels <- stats::setNames(measure_ids, measure_ids)
  original_ids <- attr(fit, "mcma_priors")$measure_id
  if (is.null(original_ids) && is.data.frame(fit$data)) {
    original_ids <- unique(as.character(fit$data[[measure_col]]))
  }
  if (length(original_ids)) {
    coefficients <- .mcma_measure_coefficients(original_ids, measure_col)
    labels <- names(coefficients)[match(paste0(measure_col, measure_ids), coefficients)]
    measure_labels[!is.na(labels)] <- labels[!is.na(labels)]
  }

  if (!summary) {
    out <- lapply(measure_ids, function(mid) {
      se_raw <- dr[[paste0("b_Se_", measure_col, mid)]]
      sp_raw <- dr[[paste0("b_Sp_", measure_col, mid)]]
      list(se = .backtransform_sesp(se_raw, bounded),
           sp = .backtransform_sesp(sp_raw, bounded))
    })
    names(out) <- unname(measure_labels)
    return(out)
  }

  # Transform and summarise each instrument separately, then stack the
  # labelled rows into a table.
  rows <- lapply(measure_ids, function(mid) {
    se_raw <- as.numeric(dr[[paste0("b_Se_", measure_col, mid)]])
    sp_raw <- as.numeric(dr[[paste0("b_Sp_", measure_col, mid)]])
    .summarise_sesp_pair(unname(measure_labels[mid]),
                         .backtransform_sesp(se_raw, bounded),
                         .backtransform_sesp(sp_raw, bounded),
                         probs)
  })
  do.call(rbind, rows)
}


#' Extract moderator slope posteriors
#'
#' @param fit A brmsfit object with moderator terms.
#' @param variable Character vector of moderator names. If NULL, extracts all.
#' @param summary Logical.
#' @param probs Quantile probabilities.
#' @param gold_column Name of the gold-standard indicator used in the model.
#'   Its group intercepts are excluded from the moderator results. If omitted,
#'   the saved model configuration is used when available; otherwise the
#'   default is `"is_gold"`.
#' @return A tibble.
#' @export
extract_slope <- function(fit,
                          variable = NULL,
                          summary  = TRUE,
                          probs    = c(0.025, 0.975),
                          gold_column = "is_gold") {

  # Corrected models prefix prevalence coefficients with b_pi_; naive
  # models use b_. Select one family so accuracy coefficients cannot enter.

  dr  <- posterior::as_draws_df(fit)
  nms <- names(dr)

  prefix <- if (any(startsWith(nms, "b_pi_"))) "b_pi_" else "b_"
  mod_cols <- nms[startsWith(nms, prefix)]
  mod_names <- substring(mod_cols, nchar(prefix) + 1L)

  # Exclude the population intercept and gold/screen group intercepts from
  # the moderator list.
  cfg <- attr(fit, "mcma_config")
  
  if (missing(gold_column) && !is.null(cfg$gold_column)) gold_column <- cfg$gold_column
  
  group_intercepts <- paste0(gold_column, c("TRUE", "FALSE", "true", "false", "1", "0"))
  keep <- !mod_names %in% c("Intercept", group_intercepts)
  mod_cols <- mod_cols[keep]
  mod_names <- mod_names[keep]

  if (length(mod_cols) == 0) {
    rlang::abort("No moderator coefficients found in the model.")
  }

  # Optionally retain only the moderator coefficients requested by the
  # caller.
  if (!is.null(variable)) {
    keep <- mod_names %in% variable
    mod_cols  <- mod_cols[keep]
    mod_names <- mod_names[keep]
  }

  if (!summary) {
    out <- lapply(mod_cols, function(c) as.numeric(dr[[c]]))
    names(out) <- mod_names
    return(out)
  }

  # Summarise each retained log-odds coefficient; with summary = FALSE the
  # earlier branch returns the full draw vectors.
  rows <- lapply(seq_along(mod_cols), function(i) {
    draws <- as.numeric(dr[[mod_cols[i]]])
    tibble::tibble(
      variable = mod_names[i],
      mean     = mean(draws),
      median   = stats::median(draws),
      sd       = stats::sd(draws),
      ci_lb    = stats::quantile(draws, probs[1], names = FALSE),
      ci_ub    = stats::quantile(draws, probs[2], names = FALSE)
    )
  })
  do.call(rbind, rows)
}


#' Extract gold vs screen prevalence difference
#'
#' Extracts the posterior difference between gold-standard and screening
#' prevalence from a comparison model.
#'
#' @param fit A brmsfit from a comparison model.
#' @param summary Logical.
#' @param probs Quantile probabilities.
#' @param gold_column Name of the gold-standard indicator used when fitting
#'   the comparison model. This identifies its gold and screening intercepts.
#' @return A tibble.
#' @export
extract_gold_screen_diff <- function(fit,
                                     summary = TRUE,
                                     probs   = c(0.025, 0.975),
                                     gold_column = "is_gold") {

  # Subtract paired posterior draws to preserve the modelled dependence
  # between the two group estimates. Positive differences mean higher
  # interview prevalence.
  dr   <- posterior::as_draws_df(fit)
  nms  <- names(dr)
  cols <- .find_gold_screen_cols(nms, gold_column = gold_column)

  gold_draws   <- stats::plogis(as.numeric(dr[[cols$gold]]))
  screen_draws <- stats::plogis(as.numeric(dr[[cols$screen]]))
  diff_draws   <- gold_draws - screen_draws

  if (!summary) {
    return(list(gold = gold_draws, screen = screen_draws, diff = diff_draws))
  }

  # Summarise both group means and the paired difference, including the
  # posterior probability that the difference exceeds zero.
  tibble::tibble(
    gold_mean   = mean(gold_draws),
    screen_mean = mean(screen_draws),
    diff_mean   = mean(diff_draws),
    diff_ci_lb  = stats::quantile(diff_draws, probs[1], names = FALSE),
    diff_ci_ub  = stats::quantile(diff_draws, probs[2], names = FALSE),
    prob_gt_zero = mean(diff_draws > 0)
  )
}


# --- Internal helpers ---

# Prepare predictions once so nonlinear parameters share the same
# posterior draws and the same sampled effects for each new study or report.
.mcma_prediction_draws <- function(fit, newdata, re_formula, allow_new_levels,
                                   sample_new_levels, nlpars = NULL) {
  prep <- brms::prepare_predictions(
    fit, newdata = newdata, re_formula = re_formula,
    allow_new_levels = allow_new_levels, sample_new_levels = sample_new_levels,
    check_response = FALSE
  )

  # The linear scale gives logits for prevalence and inner logits for
  # bounded accuracy. The caller applies the appropriate inverse transform.
  parameters <- if (is.null(nlpars)) list(NULL) else as.list(nlpars)
  draws <- lapply(parameters, function(nlpar) {
    brms::posterior_epred(prep, dpar = NULL, nlpar = nlpar, sort = FALSE,
                          scale = "linear", summary = FALSE)
  })
  if (!is.null(nlpars)) names(draws) <- nlpars
  draws
}

# Each prediction row gets a summary of exactly the draws returned by
# summary = FALSE. For new groups those draws already include heterogeneity;
# do not add a second random deviation when reporting a prediction interval.
.mcma_prediction_summary <- function(draws, probs, prediction_interval = FALSE) {
  rows <- lapply(seq_len(ncol(draws)), function(i) {
    out <- .mcma_draw_summary(draws[, i], probs)
    if (prediction_interval) {
      out$pi_l95 <- out$l95
      out$pi_u95 <- out$u95
    }
    tibble::add_column(out, .row = i, .before = 1)
  })
  do.call(rbind, rows)
}

# Let brms parse grouping syntax, including nested groups. Tau here is
# the combined random-intercept SD; slopes would need covariate-specific
# variances and must not silently be treated as intercepts.
.mcma_re_groups <- function(re_formula) {
  if (is.null(re_formula)) return(NULL)
  if (is.atomic(re_formula) && length(re_formula) == 1L && is.na(re_formula)) {
    return(character())
  }
  if (!inherits(re_formula, "formula") || length(re_formula) != 2L) {
    stop("re_formula must be NULL, NA, or a one-sided random-effects formula.", call. = FALSE)
  }

  terms <- brms::brmsterms(stats::update(re_formula, .mcma_response ~ .))$dpars$mu
  if (length(all.vars(terms$fe))) {
    stop("re_formula must contain only random-intercept terms.", call. = FALSE)
  }
  if (is.null(terms$re)) return(character())
  intercept_only <- vapply(terms$re$form, function(form) {
    term <- stats::terms(form)
    !length(attr(term, "term.labels")) && attr(term, "intercept") == 1L
  }, logical(1))
  if (!all(intercept_only)) {
    stop("extract_tau() supports random-intercept selection only, not random slopes.", call. = FALSE)
  }
  unique(terms$re$group)
}

.prev_draws <- function(fit) {

  # Recognise the intercept names used by corrected and naive models, then
  # convert log odds to probabilities. Comparison fits have separate group
  # intercepts.
  dr <- posterior::as_draws_df(fit)
  nms <- names(dr)
  
  if ("b_pi_Intercept" %in% nms) {
    stats::plogis(as.numeric(dr$b_pi_Intercept))
  } else if ("b_Intercept" %in% nms) {
    stats::plogis(as.numeric(dr$b_Intercept))
  } else {
    rlang::abort("Could not find population intercept in draws.")
  }
}

.prev_draws_logit <- function(fit) {

  # Read the same population intercept without back-transformation, for
  # calculations that must take place on the log-odds scale.
  dr <- posterior::as_draws_df(fit)
  nms <- names(dr)
  
  if ("b_pi_Intercept" %in% nms) {
    as.numeric(dr$b_pi_Intercept)
  } else if ("b_Intercept" %in% nms) {
    as.numeric(dr$b_Intercept)
  } else {
    rlang::abort("Could not find population intercept in draws.")
  }
}

.tau_draws <- function(fit, groups = NULL) {

  # Reuse the same component selection as the prediction helper.
  dr <- posterior::as_draws_df(fit)
  components <- .mcma_tau_components(dr)
  .mcma_combine_tau(dr, components, groups)
}

.mcma_tau_components <- function(dr) {
  # Match prevalence intercept SDs exactly. Se/Sp variability belongs to
  # the measurement process and must not inflate latent-prevalence intervals.
  corrected <- any(startsWith(names(dr), "b_pi_"))
  suffix <- if (corrected) "__pi_+Intercept$" else "__Intercept$"
  columns <- grep(paste0("^sd_.+", suffix), names(dr), value = TRUE)
  groups <- sub(suffix, "", substring(columns, 4L))
  stats::setNames(columns, groups)
}

.mcma_combine_tau <- function(dr, components, groups = NULL) {
  # Each grouping factor is an independent source of variation. Their
  # variances add; their SDs do not. No selected groups means no added noise.
  if (is.null(groups)) groups <- names(components)
  if (anyNA(groups) || anyDuplicated(groups) || any(!groups %in% names(components))) {
    stop("groups must name distinct prevalence random-intercept grouping factors.", call. = FALSE)
  }
  if (!length(groups)) return(rep(0, nrow(dr)))

  variances <- lapply(unname(components[groups]), function(column) as.numeric(dr[[column]])^2)
  sqrt(Reduce(`+`, variances))
}

.mcma_prediction_target <- function(dr, prediction_groups, condition_on,
                                     zero_groups = NULL) {
  # Require a complete target: each group is new, an identified existing
  # level, or explicitly fixed at zero. An omitted group is never silently removed.
  components <- .mcma_tau_components(dr)
  if (is.null(condition_on)) condition_on <- list()
  
  fixed_groups <- names(condition_on)
  if (!is.list(condition_on) ||
      (length(condition_on) && (is.null(fixed_groups) || anyNA(fixed_groups) ||
        any(!nzchar(fixed_groups)) || anyDuplicated(fixed_groups)))) {
    stop("condition_on must be a named list of existing group levels.", call. = FALSE)
  }

  # Holding a group at zero needs no saved effect for a particular report.
  # Check names so a typo cannot accidentally change the prediction target.
  if (is.null(zero_groups)) zero_groups <- character()
  if (!is.character(zero_groups) || anyNA(zero_groups) ||
      anyDuplicated(zero_groups) || any(!zero_groups %in% names(components))) {
    stop("zero_groups must name distinct prevalence random-intercept grouping factors.", call. = FALSE)
  }
  if (is.null(prediction_groups)) {
    prediction_groups <- setdiff(names(components), c(fixed_groups, zero_groups))
  }

  # Each group has exactly one role. Only new groups supply fresh
  # heterogeneity; zero-valued groups need no addition to the offset below.
  assigned_groups <- c(prediction_groups, fixed_groups, zero_groups)
  if (anyDuplicated(assigned_groups) ||
      !setequal(assigned_groups, names(components))) {
    stop("Each prevalence grouping factor must be in exactly one of prediction_groups, condition_on, or zero_groups.", call. = FALSE)
  }
  tau <- .mcma_combine_tau(dr, components, prediction_groups)

  # Match brms draw names directly, preserving the pairing between each
  # existing report effect, the population intercept, and the sampled SDs.
  offset <- rep(0, nrow(dr))
  corrected <- any(startsWith(names(dr), "b_pi_"))
  for (group in fixed_groups) {
    level <- condition_on[[group]]
    if (length(level) != 1L || is.na(level)) {
      stop("Each condition_on entry must identify one existing level.", call. = FALSE)
    }
    column <- if (corrected) {
      paste0("r_", group, "__pi[", level, ",Intercept]")
    } else {
      paste0("r_", group, "[", level, ",Intercept]")
    }
    if (!column %in% names(dr)) {
      stop("No saved prevalence effect for ", group, " = ", level, ".", call. = FALSE)
    }
    offset <- offset + as.numeric(dr[[column]])
  }

  list(tau = tau, offset = offset)
}

.mcma_draw_summary <- function(draws, probs) {

  # Use the same summary columns for ordinary and comparison prevalence.
  tibble::tibble(
    mean = mean(draws), median = stats::median(draws), sd = stats::sd(draws),
    l95 = stats::quantile(draws, probs[1], names = FALSE),
    u95 = stats::quantile(draws, probs[2], names = FALSE)
  )
}

.backtransform_sesp <- function(raw, bounded) {

  # Undo the fitted accuracy parameterization. The bounded model maps the
  # inner probability to the interval from 0.5 to 1.

  if (bounded) {
    0.5 + 0.5 * stats::plogis(raw)
  } else {
    stats::plogis(raw)
  }
}

.summarise_sesp_pair <- function(measure_id, se_draws, sp_draws, probs) {

  # Create one labelled row containing posterior means, medians, and
  # equal-tailed credible limits for sensitivity and specificity.
  tibble::tibble(
    measure_id = measure_id,
    se_mean    = mean(se_draws),
    se_median  = stats::median(se_draws),
    se_ci_lb   = stats::quantile(se_draws, probs[1], names = FALSE),
    se_ci_ub   = stats::quantile(se_draws, probs[2], names = FALSE),
    sp_mean    = mean(sp_draws),
    sp_median  = stats::median(sp_draws),
    sp_ci_lb   = stats::quantile(sp_draws, probs[1], names = FALSE),
    sp_ci_ub   = stats::quantile(sp_draws, probs[2], names = FALSE)
  )
}

.find_gold_screen_cols <- function(nms,
                                   gold_column = "is_gold") {

  # Locate the interview and screening coefficients across the naming
  # conventions produced by logical, numeric, and factor gold indicators.
  gold_candidates <- c(
    paste0("b_pi_", gold_column, "TRUE"),
    paste0("b_pi_", gold_column, "1"),
    paste0("b_pi_", gold_column, "true"),
    paste0("b_", gold_column, "TRUE"),
    paste0("b_", gold_column, "1"),
    paste0("b_", gold_column, "true")
  )

  screen_candidates <- c(
    paste0("b_pi_", gold_column, "FALSE"),
    paste0("b_pi_", gold_column, "0"),
    paste0("b_pi_", gold_column, "false"),
    paste0("b_", gold_column, "FALSE"),
    paste0("b_", gold_column, "0"),
    paste0("b_", gold_column, "false")
  )

  # Prefer explicit known coefficient names before trying broader pattern
  # matches.
  gold_col   <- intersect(gold_candidates, nms)[1]
  screen_col <- intersect(screen_candidates, nms)[1]

  if (is.na(gold_col)) {
    gold_pattern <- paste0(
      gold_column,
      ".*TRUE|",
      gold_column,
      ".*1$"
    )

    gold_col <- grep(
      gold_pattern,
      nms,
      value = TRUE,
      ignore.case = TRUE
    )[1]
  }

  if (is.na(screen_col)) {
    scr_pattern <- paste0(
      gold_column,
      ".*FALSE|",
      gold_column,
      ".*0$"
    )

    screen_col <- grep(
      scr_pattern,
      nms,
      value = TRUE,
      ignore.case = TRUE
    )[1]
  }

  # Stop if either group cannot be located; returning the wrong coefficient
  # would reverse or invalidate the comparison.
  if (is.na(gold_col) || is.na(screen_col)) {
    rlang::abort(paste0(
      "Could not find gold/screen intercept columns. ",
      "Available: ", paste(utils::head(nms, 30), collapse = ", ")
    ))
  }

  list(gold = gold_col, screen = screen_col)
}


#' Extract the concealment and over-diagnosis terms from a joint fit
#'
#' Summarises the concealment probability `c` and the over-diagnosis
#' (additional false-positive) probability `o` used by [mcma_fit_joint()].
#' A term given a prior is summarised from its posterior draws; a fixed term
#' is reported at its fixed value (the mean of the row-specific values, with
#' their range as bounds, when it varies by row); an absent term is reported
#' as zero.
#'
#' @param fit A brmsfit object returned by [mcma_fit_joint()].
#' @param summary Logical; if FALSE, returns a list with the raw posterior
#'   draws of each term (`NULL` for a term without a prior).
#' @param probs Quantile probabilities for the CI bounds.
#' @return A tibble with one row per term (`c`, `o`) giving `mode`
#'   ("none", "fixed" or "prior"), the prior string, and summary
#'   statistics, or a list of draws when `summary = FALSE`.
#' @export
extract_bias <- function(fit, summary = TRUE, probs = c(0.025, 0.975)) {

  # Read how each bias term was specified so fixed inputs are kept distinct
  # from parameters that have a posterior distribution.
  bias <- attr(fit, "mcma_config")$bias
  if (is.null(bias)) {
    stop("This fit carries no concealment / over-diagnosis settings; ",
         "fit it with mcma_fit_joint().", call. = FALSE)
  }
  
  # Fits from before the prior route record only the fixed values
  mode_of <- function(value, prior, mode) {
    if (!is.null(mode)) return(mode)
    if (!is.null(prior)) return("prior")
    if (any(value != 0)) "fixed" else "none"
  }

  # Map the user-facing c and o terms to their saved mode and the internal
  # brms coefficient names.
  spec <- list(
    c = list(value = bias$c, prior = bias$c_prior, par = "b_conc_Intercept",
             mode = mode_of(bias$c, bias$c_prior, bias$c_mode)),
    o = list(value = bias$o, prior = bias$o_prior, par = "b_overdx_Intercept",
             mode = mode_of(bias$o, bias$o_prior, bias$o_mode))
  )

  # Read posterior draws only if at least one bias term was estimated; fixed
  # inputs require no sampling output.
  needs_draws <- any(vapply(spec, function(s) s$mode == "prior", logical(1)))
  dr <- if (needs_draws) posterior::as_draws_df(fit) else NULL
  draws_of <- function(s) if (s$mode == "prior") as.numeric(dr[[s$par]]) else NULL

  if (!summary) return(lapply(spec, draws_of))

  # Estimated terms receive credible intervals. For row-specific fixed
  # inputs, the reported bounds are their range, not an uncertainty interval.
  rows <- lapply(names(spec), function(term) {
    s <- spec[[term]]
    if (s$mode == "prior") {
      d <- draws_of(s)
      tibble::tibble(
        term = term, mode = s$mode, prior = s$prior,
        mean = mean(d), median = stats::median(d), sd = stats::sd(d),
        l95 = stats::quantile(d, probs[1], names = FALSE),
        u95 = stats::quantile(d, probs[2], names = FALSE)
      )
    } else {
      v <- if (s$mode == "fixed") as.numeric(s$value) else 0
      tibble::tibble(
        term = term, mode = s$mode, prior = NA_character_,
        mean = mean(v), median = stats::median(v), sd = 0,
        l95 = min(v), u95 = max(v)
      )
    }
  })
  dplyr::bind_rows(rows)
}
