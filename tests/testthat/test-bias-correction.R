library(mcma)

testthat::test_that("the correction recovers prevalence and preserves ordinary RG", {
  testthat::expect_equal(rg_correct_bias(.137008, .85, .89,
                                        c = .30, o = .02, gamma = .5, delta = .5), .03)
  testthat::expect_equal(rg_correct_bias(.0307, 1, 1, c = .30, o = .01), .03)
  p <- c(.01, .15, .30)
  testthat::expect_equal(rg_correct_bias(p, .85, .89),
                        rg_correct(p, .85, .89, clamp = FALSE))
  testthat::expect_lt(rg_correct_bias(.01, .85, .89), 0)
  testthat::expect_equal(rg_correct_bias(.01, .85, .89, clamp = TRUE), 0)
  testthat::expect_error(rg_correct_bias(.1, .8, .6, c = .5), "denominator")
  testthat::expect_error(rg_correct_bias(.1, .8, .9, delta = 1.1), "delta")
  testthat::expect_error(rg_correct_bias(.1, .8, .9, c = NA_real_), "c must")
  testthat::expect_error(rg_correct_bias(c(.1, .2, .3), c(.8, .9), .9), "length")
})

testthat::test_that("joint fitting uses the intended observation model", {
  # Capture the actual brm arguments, then validate them using brms itself.
  testthat::local_mocked_bindings(brm = function(...) list(...), .package = "brms")
  dat <- data.frame(y = c(3, 7, 12, 16), n = 100,
                    es_id = 1:4, is_gold = c(TRUE, TRUE, FALSE, FALSE))
  inv_logit <- stats::plogis
  theta <- c(.03, .06, .03, .06)
  se <- ifelse(dat$is_gold, 1, .85)
  sp <- ifelse(dat$is_gold, 1, .89)
  for (bounded in c(TRUE, FALSE)) {
    ordinary <- mcma_fit_joint(dat, prev_center = .06, bounded = bounded)
    adjusted <- mcma_fit_joint(dat, prev_center = .06, bounded = bounded,
                              c = .30, o = .02, gamma = .5, delta = .5)
    testthat::expect_false("mcma_cmult" %in% names(ordinary$data))
    testthat::expect_equal(adjusted$data$mcma_cmult, c(.7, .7, .85, .85))
    testthat::expect_equal(adjusted$data$mcma_omult, c(.98, .98, .99, .99))
    pi <- qlogis(theta)
    Se <- qlogis(if (bounded) 2 * se - 1 else se)
    Sp <- qlogis(if (bounded) 2 * sp - 1 else sp)
    mcma_cmult <- adjusted$data$mcma_cmult
    mcma_omult <- adjusted$data$mcma_omult
    p0 <- eval(ordinary$formula$formula[[3]])
    p1 <- eval(adjusted$formula$formula[[3]])
    expected <- ifelse(dat$is_gold, .7 * theta + .02 * (1 - theta),
                       .85 * .85 * theta + (1 - .99 * .89) * (1 - theta))
    testthat::expect_equal(p0, se * theta + (1 - sp) * (1 - theta))
    testthat::expect_equal(p1, expected)
    testthat::expect_equal(rg_correct_bias(p1, se, sp, c = .3, o = .02,
                           gamma = ifelse(dat$is_gold, 1, .5),
                           delta = ifelse(dat$is_gold, 1, .5)), theta)
    testthat::expect_equal(attr(adjusted, "mcma_config")$bias$c, .3)
    args <- adjusted[c("formula", "data", "family", "prior")]
    testthat::expect_type(do.call(brms::make_stancode, args), "character")
    testthat::expect_type(do.call(brms::make_standata, args), "list")
  }
  # Optional random effects must not leave priors for nonexistent parameters.
  for (prev in list(NULL, ~ (1 | es_id))) {
    for (accuracy in list(NULL, ~ (1 | es_id))) {
      fit <- mcma_fit_joint(dat, prev_center = .06, prev_re = prev,
                            sesp_re = accuracy, c = .3, o = .02)
      args <- fit[c("formula", "data", "family", "prior")]
      testthat::expect_type(do.call(brms::make_standata, args), "list")
    }
  }
  testthat::expect_error(mcma_fit_joint(dat, prev_center = .06, c = c(.1, .2)),
                         "length")
})

testthat::test_that("priors can be placed on the concealment / over-diagnosis terms", {
  testthat::local_mocked_bindings(brm = function(...) list(...), .package = "brms")
  dat <- data.frame(y = c(3, 7, 12, 16), n = 100,
                    es_id = 1:4, is_gold = c(TRUE, TRUE, FALSE, FALSE))
  as_args <- function(f) f[c("formula", "data", "family", "prior")]

  # Concealment with a prior, over-diagnosis absent (the default)
  f  <- mcma_fit_joint(dat, prev_center = .06, c_prior = "beta(4, 12)", gamma = .5)
  fs <- deparse1(f$formula$formula)
  testthat::expect_true(grepl("(1 - conc * mcma_gwt) *", fs, fixed = TRUE))
  testthat::expect_false(grepl("overdx|mcma_omult|mcma_cmult", fs))
  testthat::expect_equal(f$data$mcma_gwt, c(1, 1, .5, .5))
  testthat::expect_true("conc" %in% names(f$formula$pforms))
  row <- f$prior[f$prior$nlpar == "conc", ]
  testthat::expect_equal(nrow(row), 1L)
  testthat::expect_equal(row$prior, "beta(4, 12)")
  testthat::expect_equal(as.numeric(row$lb), 0)
  testthat::expect_equal(as.numeric(row$ub), 1)
  code <- do.call(brms::make_stancode, as_args(f))
  testthat::expect_true(grepl("b_conc", code))
  testthat::expect_type(do.call(brms::make_standata, as_args(f)), "list")
  testthat::expect_equal(attr(f, "mcma_config")$bias$c_mode, "prior")
  testthat::expect_equal(attr(f, "mcma_config")$bias$o_mode, "none")

  # Both terms with priors
  g  <- mcma_fit_joint(dat, prev_center = .06, c_prior = "beta(4, 12)",
                       o_prior = "beta(1, 49)", delta = .5)
  gs <- deparse1(g$formula$formula)
  testthat::expect_true(grepl("(1 - overdx * mcma_dwt) *", gs, fixed = TRUE))
  testthat::expect_equal(g$data$mcma_dwt, c(1, 1, .5, .5))
  testthat::expect_true(all(c("conc", "overdx") %in% names(g$formula$pforms)))
  testthat::expect_type(do.call(brms::make_standata, as_args(g)), "list")

  # Mixed: prior on c, fixed o
  h  <- mcma_fit_joint(dat, prev_center = .06, c_prior = "beta(4, 12)", o = .02)
  hs <- deparse1(h$formula$formula)
  testthat::expect_true(grepl("conc", hs) && grepl("mcma_omult", hs))
  testthat::expect_equal(h$data$mcma_omult, c(.98, .98, .98, .98))

  # Zeroing: no prior and a zero value drop the term, giving the original model
  z <- mcma_fit_joint(dat, prev_center = .06)
  testthat::expect_identical(
    deparse1(z$formula$formula),
    deparse1(mcma_fit_joint(dat, prev_center = .06, c = 0, o = 0,
                            c_prior = NULL, o_prior = NULL)$formula$formula))
  testthat::expect_false(any(c("mcma_gwt", "mcma_dwt", "mcma_cmult", "mcma_omult") %in%
                               names(z$data)))
  testthat::expect_false(any(c("conc", "overdx") %in% names(z$formula$pforms)))

  # extract_bias() needs no draws for fixed / absent terms
  eb <- extract_bias(mcma_fit_joint(dat, prev_center = .06, c = .3))
  testthat::expect_equal(eb$mode, c("fixed", "none"))
  testthat::expect_equal(eb$mean, c(.3, 0))
  testthat::expect_equal(extract_bias(z)$mean, c(0, 0))

  # Invalid input
  testthat::expect_error(mcma_fit_joint(dat, prev_center = .06, c = .2, c_prior = "beta(4, 12)"),
                         "not both")
  testthat::expect_error(mcma_fit_joint(dat, prev_center = .06, o_prior = 3), "o_prior")
})
