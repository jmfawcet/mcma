# Prior specification for misclassification-corrected models
# ===========================================================

#' Map probability to inner scale
#'
#' Maps an external probability (>= 0.5) to the inner probability scale used
#' by the bounded parameterization, where Se = 0.5 + 0.5 * inv_logit(eta).
#' The inner probability is t = 2p - 1.
#'
#' @param p Probability value(s) on the external scale (>= 0.5).
#' @param eps Small constant to avoid 0/1 on the inner scale.
#' @return Numeric vector of inner probabilities.
#' @export
map_to_inner <- function(p, eps = 1e-4) {

  # Invert p = 0.5 + 0.5*t. Keep the inner probability away from 0 and 1 so
  # its logit stays finite.

  clamp_probability(2 * p - 1, eps, 1 - eps)
}

#' Compute prior SD from effective sample size
#'
#' Computes the logit-normal prior SD from a centre and effective sample size
#' via the delta-method approximation to a Beta distribution.
#'
#' @param p Prior centre on the probability scale (0, 1).
#' @param kappa Effective prior sample size.
#' @param bounded Logical; if TRUE, propagates uncertainty through the full
#'   bounded transformation (Se -> t = 2Se - 1 -> logit(t)) using the chain rule.
#' @param sd_method Method used when bounded = TRUE. "simple" applies the
#'   approximation directly to t = 2p - 1; "chain_rule" propagates uncertainty
#'   from p through the bounded transformation. "simple" should not be used
#'   and is included for legacy purposes.
#' @return Numeric scalar: the prior SD on the (inner) logit scale.
#' @export
prior_sd_from_kappa <- function(p, kappa, bounded = TRUE,
                                sd_method = c("chain_rule", "simple")) {

  # Translate a probability-scale prior and its effective sample size into an
  # approximate normal SD on the fitted logit scale.

  sd_method <- match.arg(sd_method)
  stopifnot(all(p > 0 & p < 1), all(kappa > 0))

  if (!bounded) {
    # Unbounded: delta method on Se directly
    sqrt(1 / ((kappa + 1) * p * (1 - p)))

  } else if (sd_method == "simple") {
    # Bounded, simple: delta method on inner prob t = 2p-1
    t <- map_to_inner(p)
    sqrt(1 / ((kappa + 1) * t * (1 - t)))

  } else {
    # Bounded, chain rule: propagate Se variance through Se -> t -> logit(t)
    t <- map_to_inner(p)

    # Start with the Beta variance on the original probability scale. The
    # derivative converts that SD through p -> 2p - 1 -> logit(2p - 1).
    sd_se <- sqrt(p * (1 - p) / (kappa + 1))
    deriv <- 2 / (t * (1 - t))
    sd_se * deriv
  }
}


#' Construct an mcma_priors object
#'
#' Creates an S3 object of class `mcma_priors` from user-supplied accuracy
#' information. Each measure in the dataset gets one row. Gold-standard
#' measures should be flagged via `is_gold = TRUE`.
#'
#' @param measure_id Character vector of measure identifiers matching the data.
#' @param se Per-measure Se prior centres (probability scale). Recycled if
#'   length 1. NA for gold-standard measures.
#' @param sp Per-measure Sp prior centres (probability scale). Recycled if
#'   length 1. NA for gold-standard measures.
#' @param kappa Effective prior sample size. Recycled if length 1. NA for
#'   gold-standard measures.
#' @param is_gold Logical vector identifying gold-standard measures. Recycled
#'   if length 1.
#' @param gold_logit_constant Logit-scale constant for gold-standard Se/Sp.
#' @return An S3 object of class `mcma_priors`.
#' @export
mcma_priors <- function(measure_id,
                        se,
                        sp,
                        kappa = 200,
                        is_gold = FALSE,
                        gold_logit_constant = 10) {

  # Build one row of prior settings per measurement instrument. Shared inputs
  # are repeated across instruments to simplify common-prior specifications.

  n <- length(measure_id)
  # Recycle scalars
  se       <- rep_len(se, n)
  sp       <- rep_len(sp, n)
  kappa    <- rep_len(kappa, n)
  is_gold  <- rep_len(is_gold, n)

  # Construct, validate and return
  tbl <- tibble::tibble(
    measure_id          = as.character(measure_id),
    se                  = se,
    sp                  = sp,
    kappa               = kappa,
    is_gold             = is_gold,
    gold_logit_constant = gold_logit_constant
  )

  # Give the table its package-specific class so print(), update(), and
  # plot() dispatch to the appropriate methods.
  tbl <- structure(
    tbl,
    class = c("mcma_priors", class(tbl))
  )

  .validate_mcma_priors(tbl)

  tbl
}


#' Validate an mcma_priors object
#'
#' Internal validation routine used by both `mcma_priors()` and
#' `update.mcma_priors()` to ensure that sensitivity and specificity priors
#' remain mathematically admissible.
#'
#' For screening measures, sensitivity and specificity must lie strictly
#' within (0, 1), and must satisfy Se + Sp > 1. The latter condition ensures
#' that the implied misclassification model is identifiable and that the
#' Rogan–Gladen correction denominator remains positive.
#'
#' Gold-standard measures are exempt from these checks because their Se/Sp
#' values are represented internally as fixed near-perfect constants.
#'
#' @param x An object of class `mcma_priors`.
#'
#' @return Invisibly returns `x` if validation succeeds.
#'
#' @details
#' This function is called automatically whenever an `mcma_priors` object is
#' created or modified. Users will normally not call it directly.
#'
#' Validation checks include:
#' \itemize{
#'   \item Sensitivity values are strictly between 0 and 1.
#'   \item Specificity values are strictly between 0 and 1.
#'   \item Sensitivity plus specificity exceeds 1.
#'   \item Required columns are present and correctly formatted.
#' }
#'
#' @keywords internal
.validate_mcma_priors <- function(x) {

  # Check the table structure before checking its values. Interview rows use
  # fixed accuracy coefficients and are exempt from screening-prior checks.
  required <- c(
    "measure_id",
    "se",
    "sp",
    "kappa",
    "is_gold",
    "gold_logit_constant"
  )

  missing <- setdiff(required, names(x))

  if (length(missing) > 0) {
    rlang::abort(
      paste(
        "Missing required column(s):",
        paste(missing, collapse = ", ")
      )
    )
  }

  # Validate only screening rows, since gold rows can legitimately have
  # missing Se, Sp, and kappa inputs.
  screen <- !as.logical(x$is_gold)

  bad_range <- screen & (
    is.na(x$se) | is.na(x$sp) |
      x$se <= 0 | x$se >= 1 |
      x$sp <= 0 | x$sp >= 1
  )

  if (any(bad_range)) {
    rlang::abort("Screening Se/Sp priors must be inside (0, 1).")
  }

  # Require the prior centres to describe a test with a positive Rogan-Gladen
  # denominator.
  bad_sum <- screen & (x$se + x$sp <= 1)

  if (any(bad_sum)) {
    rlang::abort("Screening Se + Sp must be greater than 1.")
  }

  # A positive effective sample size is needed to calculate a prior variance.
  bad_kappa <- screen & (is.na(x$kappa) | x$kappa <= 0)

  if (any(bad_kappa)) {
    rlang::abort("Kappa must be greater than 0.")
  }

  invisible(x)
}



#' Update an mcma_priors object
#'
#' S3 method for modifying a priors object without rebuilding from scratch.
#'
#' @param object An mcma_priors object to modify.
#' @param se New Se centre(s). Scalar updates all screening; named vector
#'   updates specific measures.
#' @param sp New Sp centre(s). Same recycling logic as se.
#' @param kappa New kappa value(s).
#' @param ... Additional arguments (ignored).
#' @return A new mcma_priors object with the updated values.
#' @export
update.mcma_priors <- function(object, se = NULL, sp = NULL, kappa = NULL, ...) {

  # Work on a copy of the prior table. Scalar changes apply to screening
  # rows; named Se/Sp changes select instruments by their identifiers.

  tbl <- object
  screen <- !tbl$is_gold

  if (!is.null(se)) {
    if (is.null(names(se)) && length(se) == 1L) {
      tbl$se[screen] <- se
    } else if (!is.null(names(se))) {
      idx <- match(names(se), tbl$measure_id)
      idx <- idx[!is.na(idx)]
      tbl$se[idx] <- se[tbl$measure_id[idx]]
    } else {
      tbl$se[screen] <- rep_len(se, sum(screen))
    }
  }

  # Apply the corresponding update to specificity while retaining unchanged
  # columns.
  if (!is.null(sp)) {
    if (is.null(names(sp)) && length(sp) == 1L) {
      tbl$sp[screen] <- sp
    } else if (!is.null(names(sp))) {
      idx <- match(names(sp), tbl$measure_id)
      idx <- idx[!is.na(idx)]
      tbl$sp[idx] <- sp[tbl$measure_id[idx]]
    } else {
      tbl$sp[screen] <- rep_len(sp, sum(screen))
    }
  }

  # Adjust prior certainty separately from prior centres.
  if (!is.null(kappa)) {
    if (length(kappa) == 1L) {
      tbl$kappa[screen] <- kappa
    } else {
      tbl$kappa[screen] <- rep_len(kappa, sum(screen))
    }
  }

  # Restore the specialised table class and validate the updated values
  # before returning them.
  tbl <- structure(
    tbl,
    class = c("mcma_priors", setdiff(class(tbl), "mcma_priors"))
  )

  .validate_mcma_priors(tbl)

  tbl
}


#' Print an mcma_priors object
#'
#' @param x An mcma_priors object.
#' @param ... Additional arguments (ignored).
#' @export
print.mcma_priors <- function(x, ...) {

  # Show how many instruments are treated as screening or gold-standard
  # measures, followed by their stored prior settings.

  n_gold   <- sum(x$is_gold)
  n_screen <- sum(!x$is_gold)
  cat(sprintf("mcma_priors: %d measure(s) (%d screening, %d gold-standard)\n",
              nrow(x), n_screen, n_gold))
  cat("\n")
  print(tibble::as_tibble(x), n = nrow(x))
  invisible(x)
}


#' Convert mcma_priors to brms prior
#'
#' Translates an `mcma_priors` object into a `brms` prior object. This is
#' called internally by `mcma_fit()` but is exported for debugging.
#'
#' @param priors An mcma_priors object.
#' @param prev_center Prevalence prior centre on the probability scale.
#' @param bounded Logical; if TRUE, uses bounded parameterization.
#' @param prev_prior_sd SD of the prevalence prior on the logit scale.
#' @param tau_prior_sd SD of the prior on the prevalence RE SD.
#' @param sesp_re_sd SD of the prior on Se/Sp RE SDs.
#' @param moderator_terms Character vector of moderator variable names.
#' @param moderator_prior_sd SD of the prior on moderator coefficients.
#' @param has_prev_re Logical; whether the model includes prevalence REs.
#' @param has_sesp_re Logical; whether the model includes Se/Sp REs.
#' @param sigma_bias Random perturbation SD for prior centres on logit scale.
#' @param delta_se Systematic shift to Se prior centres on probability scale.
#' @param delta_sp Systematic shift to Sp prior centres on probability scale.
#' @param rho_bias Correlation of Se/Sp bias perturbations.
#' @inheritParams prior_sd_from_kappa
#' @param measure_col Name of the instrument identifier column used in the
#'   model formula. It determines the coefficient names to which the
#'   instrument-specific accuracy priors are attached.
#' @return A brms prior object.
#' @export
as_brms_prior <- function(priors,
                          prev_center,
                          bounded          = TRUE,
                          prev_prior_sd    = 1.5,
                          tau_prior_sd     = 1.0,
                          sesp_re_sd       = 0.5,
                          moderator_terms  = NULL,
                          moderator_prior_sd = 1.0,
                          has_prev_re      = TRUE,
                          has_sesp_re      = TRUE,
                          sigma_bias       = 0,
                          delta_se         = 0,
                          delta_sp         = 0,
                          rho_bias         = -0.5,
                          sd_method = "chain_rule",
                          measure_col = "measure_id") {

  # Translate the human-readable prior table into priors attached to specific
  # brms coefficients. Probability-scale inputs are converted to the fitted
  # parameter scale.

  prev_logit <- stats::qlogis(prev_center)

  # --- Prevalence priors ---
  pri <- brms::set_prior(
    sprintf("normal(%0.6f, %0.6f)", prev_logit, prev_prior_sd),
    nlpar = "pi"
  )

  # Add a heterogeneity prior only when the prevalence model actually
  # contains random effects.
  if (has_prev_re) {
    pri <- c(pri, brms::set_prior(
      sprintf("normal(0, %0.6f)", tau_prior_sd),
      nlpar = "pi", class = "sd"
    ))
  }

  # Regularise how much instrument accuracy can vary between studies; brms
  # constrains these SD parameters to be nonnegative.
  # --- Se/Sp RE SD priors ---
  if (has_sesp_re) {
    pri <- c(
      pri,
      brms::set_prior(sprintf("normal(0, %0.6f)", sesp_re_sd),
                       nlpar = "Se", class = "sd"),
      brms::set_prior(sprintf("normal(0, %0.6f)", sesp_re_sd),
                       nlpar = "Sp", class = "sd")
    )
  }

  # Attach zero-centred priors to the expanded moderator coefficients
  # supplied by mcma_fit().
  # --- Moderator priors ---
  if (!is.null(moderator_terms) && length(moderator_terms) > 0) {
    for (mod in moderator_terms) {
      pri <- c(pri, brms::set_prior(
        sprintf("normal(0, %0.6f)", moderator_prior_sd),
        nlpar = "pi", class = "b", coef = mod
      ))
    }
  }

  # Let brms name the instrument coefficients, including spaces and hyphens.
  coefficients <- .mcma_measure_coefficients(priors$measure_id, measure_col)

  # --- Per-measure Se/Sp priors ---
  gold   <- priors[priors$is_gold, , drop = FALSE]
  screen <- priors[!priors$is_gold, , drop = FALSE]

  # Fix the gold instrument coefficients at the requested near-perfect value
  # on the fitted scale.
  # Gold-standard measures: constant(10) (or user-specified constant)
  for (i in seq_len(nrow(gold))) {
    const <- gold$gold_logit_constant[i]
    coef <- coefficients[[gold$measure_id[i]]]
    pri <- c(
      pri,
      brms::prior_string(sprintf("constant(%g)", const),
                          class = "b", coef = coef, nlpar = "Se"),
      brms::prior_string(sprintf("constant(%g)", const),
                          class = "b", coef = coef, nlpar = "Sp")
    )
  }

  # Screening measures: logit-normal priors
  n_screen <- nrow(screen)
  if (n_screen > 0) {
    se_centres <- screen$se + delta_se
    sp_centres <- screen$sp + delta_sp

    # Clamping
    lower <- if (bounded) 0.5 else 1e-4
    upper <- 1 - 1e-4

    se_centres <- clamp_probability(se_centres, lower, upper)
    sp_centres <- clamp_probability(sp_centres, lower, upper)

    # For simulation studies, optionally introduce correlated errors into the
    # assumed accuracy centres on the logit scale.
    # Random perturbation (for simulation use)
    if (sigma_bias > 0 && n_screen > 0) {
      Sigma <- matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
      # Preserve the Se/Sp columns when only one screening instrument exists.
      E <- matrix(MASS::mvrnorm(n = n_screen, mu = c(0, 0), Sigma = Sigma), ncol = 2)
      se_logit <- stats::qlogis(se_centres) + E[, 1]
      sp_logit <- stats::qlogis(sp_centres) + E[, 2]
      se_centres <- stats::plogis(se_logit)
      sp_centres <- stats::plogis(sp_logit)

      # Clamping
      lower <- if (bounded) 0.5 else 1e-4
      upper <- 1 - 1e-4

      se_centres <- clamp_probability(se_centres, lower, upper)
      sp_centres <- clamp_probability(sp_centres, lower, upper)
    }

    # Build separate Se and Sp priors for each instrument, using its
    # effective sample size to determine their spreads.
    for (i in seq_len(n_screen)) {
      mid   <- screen$measure_id[i]
      kap   <- screen$kappa[i]
      coef  <- coefficients[[mid]]

      if (bounded) {
        # Bounded: map to inner scale, compute prior on inner logit
        se_t <- map_to_inner(se_centres[i])
        sp_t <- map_to_inner(sp_centres[i])
        se_sd <- prior_sd_from_kappa(se_centres[i], kap, bounded = TRUE, sd_method = sd_method)
        sp_sd <- prior_sd_from_kappa(sp_centres[i], kap, bounded = TRUE, sd_method = sd_method)
        se_mu <- stats::qlogis(se_t)
        sp_mu <- stats::qlogis(sp_t)
      } else {
        # Unbounded: standard logit-normal
        se_sd <- prior_sd_from_kappa(se_centres[i], kap, bounded = FALSE, sd_method = sd_method)
        sp_sd <- prior_sd_from_kappa(sp_centres[i], kap, bounded = FALSE, sd_method = sd_method)
        se_mu <- stats::qlogis(se_centres[i])
        sp_mu <- stats::qlogis(sp_centres[i])
      }

      # Attach each normal distribution to the exact coefficient name
      # generated from the instrument column and label.
      pri <- c(
        pri,
        brms::prior_string(sprintf("normal(%0.6f, %0.6f)", se_mu, se_sd),
                            class = "b", coef = coef, nlpar = "Se"),
        brms::prior_string(sprintf("normal(%0.6f, %0.6f)", sp_mu, sp_sd),
                            class = "b", coef = coef, nlpar = "Sp")
      )
    }
  }

  pri
}


.mcma_match_priors <- function(priors, measures) {

  # A shared prior table can cover more instruments than a particular subset.
  # Keep only the rows that have corresponding model coefficients.
  measures <- unique(as.character(measures))
  missing <- setdiff(measures, priors$measure_id)
  unused <- setdiff(priors$measure_id, measures)
  if (length(missing)) {
    stop("Missing Se/Sp priors for measure(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (length(unused)) {
    warning("Ignoring priors for measure(s) not present in data: ",
            paste(unused, collapse = ", "), call. = FALSE)
  }

  priors[priors$measure_id %in% measures, , drop = FALSE]
}

.mcma_measure_coefficients <- function(measure_ids, measure_col) {

  # A one-row-per-instrument design lets brms provide its own safe names.
  # Each row selects exactly one coefficient, so label order cannot mix priors.
  labels <- unique(as.character(measure_ids))
  levels <- labels
  if (length(levels) == 1L) levels <- c(levels, paste0(levels, "_mcma_unused"))
  data <- data.frame(.mcma_response = rep(0, length(levels)))
  data[[measure_col]] <- levels
  formula <- stats::reformulate(paste0("`", measure_col, "`"),
                                response = ".mcma_response", intercept = FALSE)
  design <- brms::make_standata(formula, data = data)$X
  coefficients <- colnames(design)[max.col(design, ties.method = "first")]

  stats::setNames(coefficients[seq_along(labels)], labels)
}
