# Check the new prediction interface without running Stan sampling.
# These draws make the expected summaries and variance components explicit.

testthat::test_that("random-effect formulas select the intended heterogeneity", {
  dr <- posterior::as_draws_df(data.frame(
    b_pi_Intercept = c(-3, -2, -1),
    sd_es_id__pi_Intercept = c(3, 3, 3),
    sd_refid__pi_Intercept = c(4, 4, 4),
    sd_es_id__Se_Intercept = c(100, 100, 100)
  ))
  testthat::expect_equal(extract_tau(dr, re_formula = NULL)$mean, 5)
  testthat::expect_equal(extract_tau(dr, re_formula = ~(1 | es_id))$mean, 3)
  testthat::expect_equal(extract_tau(dr, re.form = ~(1 | refid))$mean, 4)
  testthat::expect_equal(extract_tau(dr, re_formula = NA)$mean, 0)
  testthat::expect_equal(extract_tau(dr, re_formula = ~0)$mean, 0)
  testthat::expect_equal(extract_tau(dr), extract_tau(dr, re_formula = NULL))
  testthat::expect_error(extract_tau(dr, re_formula = ~(1 | typo)), "grouping factors")
  testthat::expect_error(extract_tau(dr, re_formula = ~(1 + x | es_id)), "random slopes")
  testthat::expect_error(extract_tau(dr, re_formula = ~x), "random-intercept")
  testthat::expect_error(extract_tau(dr, groups = "es_id", re_formula = NULL), "not both")
  testthat::expect_error(extract_tau(dr, re.form = NA, re_formula = NA), "only one")

  expected <- plogis(dr$b_pi_Intercept)
  testthat::expect_equal(extract_prevalence(dr, summary = FALSE), expected)
  out <- extract_prevalence(dr, prediction_interval = TRUE,
                            zero_groups = c("es_id", "refid"))
  testthat::expect_equal(out$mean, mean(expected))
  testthat::expect_equal(out$pi_l95, out$l95)
  testthat::expect_equal(out$pi_u95, out$u95)
})

testthat::test_that("prevalence summaries and raw draws use the same prediction target", {
  eta <- matrix(c(-3, -2, -1, 0, -2, -1, 0, 1), ncol = 2)
  fit <- list(formula = list(pforms = list(pi = pi ~ 1)))
  calls <- list()
  testthat::local_mocked_bindings(
    prepare_predictions = function(...) { calls$prepare <<- list(...); list() },
    posterior_epred = function(object, dpar, nlpar, sort, scale, summary) {
      calls$nlpar <<- nlpar
      testthat::expect_equal(scale, "linear")
      eta
    }, .package = "brms"
  )
  nd <- data.frame(es_id = c("new1", "new2"))
  re <- ~(1 | es_id)
  raw <- extract_prevalence(fit, newdata = nd, re.form = re,
                            allow_new_levels = TRUE, summary = FALSE)
  testthat::expect_equal(raw, plogis(eta))
  testthat::expect_equal(calls$nlpar, "pi")
  testthat::expect_equal(calls$prepare$newdata, nd)
  testthat::expect_equal(calls$prepare$re_formula, re)
  testthat::expect_true(calls$prepare$allow_new_levels)
  testthat::expect_equal(calls$prepare$sample_new_levels, "gaussian")

  out <- extract_prevalence(fit, newdata = nd, re_formula = re,
                            prediction_interval = TRUE, probs = c(.1, .9))
  testthat::expect_equal(out$.row, 1:2)
  testthat::expect_equal(out$mean, colMeans(raw))
  testthat::expect_equal(out$l95, apply(raw, 2, quantile, .1, names = FALSE))
  testthat::expect_equal(out$u95, apply(raw, 2, quantile, .9, names = FALSE))
  testthat::expect_equal(out$pi_l95, out$l95)
  testthat::expect_equal(out$pi_u95, out$u95)
  testthat::expect_equal(extract_prevalence(fit, re_formula = NA,
                           transform = FALSE, summary = FALSE), eta)

  fit$formula$pforms <- list()
  extract_prevalence(fit, newdata = NULL, sample_new_levels = "old_levels")
  testthat::expect_null(calls$nlpar)
  testthat::expect_equal(calls$prepare$sample_new_levels, "old_levels")
  extract_prevalence(fit, re_formula = NA, sample_new_levels = "uncertainty")
  testthat::expect_equal(calls$prepare$sample_new_levels, "uncertainty")
  testthat::expect_error(extract_prevalence(fit, re.form = NA, re_formula = NA), "only one")
  testthat::expect_error(extract_prevalence(fit, newdata = nd, zero_groups = "refid"), "not both")
  testthat::expect_error(extract_prevalence(fit, sample_new_levels = "typo"), "arg")
  testthat::expect_error(extract_prevalence(fit, sample_new_levels = "old_levels"), "Supply newdata")
})

testthat::test_that("accuracy predictions share preparation and retain the fitted scale", {
  eta <- matrix(c(-2, -1, 0, 1, -1, 0, 1, 2), ncol = 2)
  n_prepare <- 0L
  testthat::local_mocked_bindings(
    prepare_predictions = function(...) { n_prepare <<- n_prepare + 1L; list() },
    posterior_epred = function(object, dpar, nlpar, sort, scale, summary) {
      if (nlpar == "Se") eta else -eta
    }, .package = "brms"
  )
  fit <- structure(list(), mcma_config = list(bounded = TRUE))
  raw <- extract_sesp(fit, newdata = NULL, re_formula = NA, summary = FALSE)
  testthat::expect_equal(n_prepare, 1L)
  testthat::expect_equal(raw$se, .5 + .5 * plogis(eta))
  testthat::expect_equal(raw$sp, .5 + .5 * plogis(-eta))
  out <- extract_sesp(fit, re.form = NA)
  testthat::expect_equal(out$.row, 1:2)
  testthat::expect_equal(out$se_mean, colMeans(raw$se))
  testthat::expect_equal(out$sp_ci_lb, apply(raw$sp, 2, quantile, .025, names = FALSE))
  testthat::expect_equal(extract_sesp(fit, newdata = NULL, bounded = FALSE,
                                     summary = FALSE)$se, plogis(eta))
  testthat::expect_error(extract_sesp(fit, re.form = NA, re_formula = NA), "only one")
  testthat::expect_error(extract_sesp(fit, sample_new_levels = "old_levels"), "Supply newdata")
})

testthat::test_that("brms evaluates categorical moderators and interactions at new rows", {
  # Use known coefficient draws with brms' prediction machinery.
  # Empty fits build model structure without compiling or sampling a Stan model.
  dat <- data.frame(y = c(1, 2, 3, 4), n = 100, es_id = 1:4,
                    measure_id = c("A", "B", "A", "B"),
                    setting = factor(c("a", "a", "b", "b")), x = c(0, 1, 0, 1))
  for (corrected in c(FALSE, TRUE)) {
    priors <- if (corrected) mcma_priors(c("A", "B"), .85, .9) else NULL
    fit <- mcma_fit(dat, priors = priors, prev_center = .06, prev_re = NULL,
                    sesp_re = NULL, moderators = ~setting * x, bounded = TRUE,
                    backend = "rstan", empty = TRUE)
    m <- cbind(Intercept = c(-3, -2), settingb = c(.2, .3), x = c(.4, .5),
                `settingb:x` = c(.1, .2))
    expected <- m[, 1] + outer(m[, 2], c(0, 0, 1, 1)) +
      outer(m[, 3], dat$x) + outer(m[, 4], c(0, 0, 0, 1))
    colnames(m) <- paste0(if (corrected) "b_pi_" else "b_", colnames(m))
    if (corrected) {
      m <- cbind(m, b_Se_measure_idA = c(1, 2), b_Sp_measure_idA = c(2, 3),
                  b_Se_measure_idB = c(1, 2), b_Sp_measure_idB = c(2, 3))
    }
    dr <- posterior::as_draws_matrix(m)
    testthat::with_mocked_bindings({
      actual <- extract_prevalence(fit, newdata = dat, re_formula = NA,
                                    summary = FALSE, transform = FALSE)
      testthat::expect_equal(actual, unname(expected))
      tab <- extract_prevalence(fit, newdata = dat, re_formula = NA)
      testthat::expect_equal(tab$mean, colMeans(plogis(expected)))
    }, as_draws_matrix = function(...) dr, ndraws = function(...) nrow(dr),
       .package = "brms")
  }
})
