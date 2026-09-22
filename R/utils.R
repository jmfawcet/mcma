# Algebraic utilities for the Rogan-Gladen correction
# ====================================================

#' Rogan-Gladen correction
#'
#' Corrects an observed prevalence for misclassification given known Se and Sp.
#'
#' @param p_obs Observed prevalence (numeric vector).
#' @param se Sensitivity (numeric vector).
#' @param sp Specificity (numeric vector).
#' @param clamp Logical; if TRUE, result is bounded to \[0, 1\].
#' @return Corrected prevalence estimate(s).
#' @export
rg_correct <- function(p_obs, se, sp, clamp = TRUE) {

  # Remove expected false positives and rescale by the test discrimination,
  # Se + Sp - 1. This correction treats the supplied accuracies as known.
  denom <- se + sp - 1

  # Stop when the inverse is undefined because the test provides no
  # discrimination between cases and non-cases.
  if (any(abs(denom) < 1e-12)) {
    rlang::abort("Se + Sp must not equal 1 (denominator is zero).")
  }
  out <- (p_obs + sp - 1) / denom

  # Optionally truncate impossible prevalence estimates to [0, 1]. This makes
  # them admissible but can introduce boundary bias.
  if (clamp) out <- pmin(1, pmax(0, out))
  out
}

#' Forward misclassification model
#'
#' Returns the expected observed prevalence given true prevalence, Se, and Sp.
#'
#' @param theta True prevalence (numeric vector).
#' @param se Sensitivity (numeric vector).
#' @param sp Specificity (numeric vector).
#' @return Expected observed prevalence.
#' @export
rg_forward <- function(theta, se, sp) {
  # Add the probability of a true case testing positive to the probability of
  # a non-case testing positive.
  se * theta + (1 - sp) * (1 - theta)
}

#' Implied sensitivity
#'
#' Solves the forward misclassification model algebraically for Se.
#'
#' @param p_obs Observed prevalence.
#' @param theta True prevalence.
#' @param sp Specificity.
#' @return Implied sensitivity.
#' @export
implied_se <- function(p_obs, theta, sp) {

  # Rearrange the observation equation to find the sensitivity compatible
  # with the supplied prevalence and specificity.
  stopifnot(all(theta > 0 & theta < 1))

  # Subtract the expected false-positive contribution, then divide by the
  # fraction of true cases.
  (p_obs - (1 - sp) * (1 - theta)) / theta
}

#' Implied specificity
#'
#' Solves the forward misclassification model algebraically for Sp.
#'
#' @param p_obs Observed prevalence.
#' @param theta True prevalence.
#' @param se Sensitivity.
#' @return Implied specificity.
#' @export
implied_sp <- function(p_obs, theta, se) {

  # Rearrange the observation equation to find the specificity compatible
  # with the supplied prevalence and sensitivity.

  stopifnot(all(theta > 0 & theta < 1))

  # Solve for the true-negative classification probability; prevalence must
  # be below one to avoid a zero denominator.
  (p_obs - 1 + theta * (1 - se)) / (theta - 1)
}

#' Iso-correction line
#'
#' Returns (Se, Sp) pairs jointly consistent with the given true and observed
#' prevalence. Useful for overlay on joint Se/Sp posterior plots.
#'
#' @param theta True prevalence.
#' @param p_obs Observed prevalence.
#' @param se_range Numeric vector of Se values to evaluate.
#' @param lb Lower bound for filtering (default 0.5).
#' @param ub Upper bound for filtering (default 1.0).
#' @return A tibble with columns `se` and `sp`.
#' @export
iso_correction_line <- function(theta, p_obs,
                                se_range = seq(0.50, 1.00, by = 0.005),
                                lb = 0.5, ub = 1.0) {

  # For a fixed true and observed prevalence, calculate all requested Se/Sp
  # pairs that produce that same observation probability.
  stopifnot(theta > 0, theta < 1, p_obs > 0, p_obs < 1)
  sp_vals <- 1 - (p_obs - theta * se_range) / (1 - theta)
  out <- tibble::tibble(se = se_range, sp = sp_vals)

  # Retain only accuracy pairs within the displayed limits and with positive
  # test discrimination.
  out <- out[out$se >= lb & out$se <= ub & out$sp >= lb & out$sp <= ub &
               out$se + out$sp > 1, ]
  out
}


#' Clamp values to optional limits
#'
#' Replaces values below a lower limit or above an upper limit with that limit.
#' This internal helper keeps accuracy prior centres within their allowed ranges.
#'
#' @param p Numeric vector of values to clamp.
#' @param lower Numeric lower limit, or `NULL` to omit the lower bound.
#' @param upper Numeric upper limit, or `NULL` to omit the upper bound.
#' @details The lower bound is applied first with [base::pmax()], followed by
#'   the upper bound with [base::pmin()]. Scalar limits apply to every element;
#'   vector limits follow the usual R recycling rules. Callers must supply
#'   compatible limits with `lower <= upper` when both bounds are present.
#'   Inputs are not validated and missing values propagate.
#'
#'   With both limits set to `NULL`, `p` is returned unchanged. The function
#'   only enforces the supplied limits; it does not automatically restrict
#'   values to the probability range `0` to `1`.
#' @return A numeric vector with the requested limits applied. With scalar
#'   limits, the result has the same length as `p`.
#' @seealso [safe_qlogis()], [map_to_inner()]
#' @examples
#' mcma:::clamp_probability(c(-0.1, 0.25, 1.2), lower = 0, upper = 1)
#' mcma:::clamp_probability(c(0, 0.2, 1), lower = 0.01)
#' @keywords internal
clamp_probability = function(p, lower = NULL, upper = NULL) {
  # Replace values outside the requested limits with the nearest limit.
  # Either limit can be omitted independently.
  if (!is.null(lower))
    p = pmax(p, lower)

  # Apply the upper bound after the lower bound; pmin/pmax work element by
  # element on a vector.
  if (!is.null(upper))
    p = pmin(p, upper)
  p
}

#' Logit transformation with finite boundary values
#'
#' Clamps probabilities away from zero and one before converting them to log
#' odds with [stats::qlogis()].
#'
#' @param p Numeric vector of probabilities. Values outside the probability
#'   range are clipped to the specified limits as well.
#' @param eps Positive numeric scalar smaller than `0.5`. Probabilities are
#'   clipped to `eps` and `1 - eps` before transformation; defaults to `1e-4`.
#' @details This internal helper calls [clamp_probability()] and then applies
#'   [stats::qlogis()]. With the default `eps`, inputs of zero and one have
#'   finite log odds of approximately -9.21 and 9.21. Missing inputs remain
#'   missing. The boundary values depend on `eps`. Callers
#'   are responsible for supplying a valid `eps`; inputs are not validated.
#' @return A numeric vector of log odds with the same length as `p`.
#' @seealso [clamp_probability()], [stats::qlogis()]
#' @examples
#' mcma:::safe_qlogis(c(0, 0.5, 1))
#' mcma:::safe_qlogis(c(0.01, 0.25, 0.99))
#' @keywords internal
safe_qlogis = function(p, eps = 1e-4) {

  # Keep probabilities away from 0 and 1 before taking log odds, avoiding
  # infinite values in prior specifications.

  qlogis(clamp_probability(p, eps, 1 - eps))
}
