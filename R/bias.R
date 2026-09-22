#' Rogan-Gladen correction with concealment and additional false positives
#'
#' Applies the correction using effective sensitivity `(1 - c * gamma) * se`
#' and effective specificity `(1 - o * delta) * sp`.
#'
#' @param p_obs Observed prevalence.
#' @param se,sp Baseline sensitivity and specificity, before the additional
#'   concealment and false-positive processes.
#' @param c Concealment probability in diagnostic interviews.
#' @param o Additional false-positive probability in diagnostic interviews.
#' @param gamma,delta Screening weights for concealment and additional false
#'   positives, respectively. Each is between zero and one. In particular,
#'   `o * delta` is the probability that the additional process makes an
#'   otherwise correctly classified non-case screen positive.
#' @param clamp Logical; if TRUE, truncate estimates to `[0, 1]`. Defaults to
#'   FALSE because truncation can introduce bias.
#' @details All probability inputs must be finite and in `[0, 1]`, and have
#'   length one or a common vector length. The effective sensitivity plus
#'   effective specificity must differ from one. Negative denominators are
#'   allowed. For diagnostic interviews, use `se = sp = 1` and
#'   `gamma = delta = 1`, giving `(p_obs - o) / (1 - c - o)`.
#'   With `c = o = 0`, this is the ordinary Rogan-Gladen correction.
#'   Inputs are treated as known; this function does not propagate their
#'   uncertainty. Do not apply these additional effects to accuracy estimates
#'   that already include the same effects.
#' @return Corrected prevalence estimate(s).
#' @seealso [rg_correct()], [mcma_fit_joint()]
#' @examples
#' rg_correct_bias(0.137008, se = 0.85, sp = 0.89,
#'                 c = 0.30, o = 0.02, gamma = 0.5, delta = 0.5)
#' @export
rg_correct_bias <- function(p_obs, se, sp, c = 0, o = 0,
                            gamma = 1, delta = 1, clamp = FALSE) {

  # Treat all supplied error rates as known. A single value applies to every
  # observation; vectors allow observation-specific corrections.
  inputs <- list(p_obs = p_obs, se = se, sp = sp,
                 c = c, o = o, gamma = gamma, delta = delta)

  # Use the longest input as the intended output length; every other input
  # must be either scalar or that length.
  n <- max(lengths(inputs))
  for (nm in names(inputs)) {
    .mcma_validate_probability(inputs[[nm]], nm, n)
  }

  if (!is.logical(clamp) || length(clamp) != 1L || is.na(clamp)) {
    stop("clamp must be TRUE or FALSE.", call. = FALSE)
  }

  # Apply the ordinary Rogan-Gladen inverse using the effective accuracy
  # after the extra bias processes.
  mult <- .mcma_bias_multipliers(n, c, o, gamma, delta)

  rg_correct(p_obs, se * mult$se, sp * mult$sp, clamp = clamp)
}

.mcma_validate_probability <- function(x, name, n) {

  # Reject missing values, invalid probabilities, and incompatible lengths
  # before R can silently recycle mismatched inputs.
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      !length(x) || !length(x) %in% c(1L, n) ||
      any(!is.finite(x)) || any(x < 0 | x > 1)) {
    stop(name, " must contain finite probabilities in [0, 1], ",
         "with length one or ", n, ".", call. = FALSE)
  }

  # Return the validated input without printing it, allowing callers to use
  # this purely as a check.
  invisible(x)
}

.mcma_bias_multipliers <- function(n, c, o, gamma, delta, is_gold = FALSE) {

  # Calculate the fractions of baseline sensitivity and specificity that
  # remain after the additional reporting errors.
  inputs <- list(c = c, o = o, gamma = gamma, delta = delta)

  for (nm in names(inputs)) {
    .mcma_validate_probability(inputs[[nm]], nm, n)
  }

  if (!(is.logical(is_gold) || is.numeric(is_gold)) ||
      !length(is_gold) || !length(is_gold) %in% c(1L, n) ||
      anyNA(is_gold) || any(!is_gold %in% c(0, 1))) {
    stop("The gold indicator must contain TRUE/FALSE or 0/1, ",
         "with length one or ", n, ".", call. = FALSE)
  }

  # The weight is one for interviews and the supplied attenuation for
  # screens; subtracting the weighted error gives the retained fraction.
  w <- .mcma_bias_weights(n, gamma, delta, is_gold)   # weights shared with the prior route in mcma_fit_joint()

  list(
    se = 1 - c * w$gwt,
    sp = 1 - o * w$dwt
  )
}

# per-row weights for the concealment (gamma) and over-diagnosis (delta) terms;
# gold rows always 1, screening rows gamma / delta. Used by the fixed multipliers
# above and, as data columns, by the prior-driven terms in mcma_fit_joint().
.mcma_bias_weights <- function(n, gamma, delta, is_gold = FALSE) {

  # Give diagnostic interviews the full reporting-error rates. For screens,
  # gamma and delta set how strongly each error process operates.
  .mcma_validate_probability(gamma, "gamma", n)
  .mcma_validate_probability(delta, "delta", n)

  if (!(is.logical(is_gold) || is.numeric(is_gold)) ||
      !length(is_gold) || !length(is_gold) %in% c(1L, n) ||
      anyNA(is_gold) || any(!is_gold %in% c(0, 1))) {
    stop("The gold indicator must contain TRUE/FALSE or 0/1, ",
         "with length one or ", n, ".", call. = FALSE)
  }

  # Convert the group flag to zero/one values and repeat scalar inputs to
  # obtain one weight per data row.
  gold <- rep_len(as.numeric(is_gold), n)

  list(
    gwt = gold + (1 - gold) * rep_len(gamma, n),
    dwt = gold + (1 - gold) * rep_len(delta, n)
  )
}

# a prior for the concealment / over-diagnosis parameters must be NULL or a single
# Stan distribution string (as accepted by brms::set_prior), e.g. "beta(4, 12)"
.mcma_validate_bias_prior <- function(x, name) {

  # NULL means no estimated bias parameter. Otherwise require one nonempty
  # distribution specification; brms checks the distribution itself later.
  if (is.null(x)) return(invisible(NULL))

  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(name, " must be NULL or a single Stan distribution string, e.g. \"beta(4, 12)\".",
         call. = FALSE)
  }

  invisible(x)
}
