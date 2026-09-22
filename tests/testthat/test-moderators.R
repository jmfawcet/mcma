library(mcma)

testthat::test_that("categorical moderators and interactions receive valid priors", {
  testthat::local_mocked_bindings(brm = function(...) list(...), .package = "brms")
  dat <- data.frame(
    y = c(8, 12, 16, 9, 13, 17, 10, 14, 18, 11, 15, 19),
    n = 100, es_id = 1:12, measure_id = rep(c("A", "B"), 6),
    setting = factor(rep(c("field", "hospital", "community"), 4),
                     levels = c("field", "hospital", "community")),
    x_mod = rep(c(-1.5, -0.5, 0.5, 1.5), each = 3)
  )
  pri <- mcma_priors(c("A", "B"), se = .85, sp = .9)
  cases <- list(
    list(formula = ~ setting, expected = c("settinghospital", "settingcommunity")),
    list(formula = ~ setting * x_mod,
         expected = c("settinghospital", "settingcommunity", "x_mod",
                      "settinghospital:x_mod", "settingcommunity:x_mod")),
    list(formula = ~ x_mod, expected = "x_mod")
  )
  for (bounded in c(FALSE, TRUE)) {
    for (case in cases) {
      f <- mcma_fit(dat, priors = pri, prev_center = .06,
                    moderators = case$formula, bounded = bounded)
      mp <- as.data.frame(f$prior)
      mp <- mp[mp$class == "b" & mp$nlpar == "pi" & nzchar(mp$coef), ]
      testthat::expect_setequal(mp$coef, case$expected)
      testthat::expect_true(all(mp$prior == "normal(0, 1.000000)"))
      args <- f[c("formula", "data", "family", "prior")]
      testthat::expect_type(do.call(brms::make_standata, args), "list")
      testthat::expect_type(do.call(brms::make_stancode, args), "character")
    }
  }

  # Custom contrasts must also use the coefficient names recognized by brms.
  stats::contrasts(dat$setting) <- stats::contr.sum(3)
  f <- mcma_fit(dat, priors = pri, prev_center = .06, moderators = ~ setting)
  mp <- as.data.frame(f$prior)
  testthat::expect_setequal(mp$coef[mp$nlpar == "pi" & nzchar(mp$coef)],
                           c("setting1", "setting2"))
  testthat::expect_type(do.call(brms::make_standata,
                                f[c("formula", "data", "family", "prior")]), "list")
})
