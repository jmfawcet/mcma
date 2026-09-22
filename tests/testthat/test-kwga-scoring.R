# Synthetic paired draws: gold prevalence around 2.5%, screening around 14.6%,
# so any Sp at or below about 0.85 implies a negative Rogan-Gladen correction,
# while Sp near 0.875 reproduces the gold prevalence.
kwga_test_draws <- function(n = 500) {
  set.seed(20260922)
  list(gold_draws   = plogis(rnorm(n, qlogis(0.025), 0.35)),
       screen_draws = plogis(rnorm(n, qlogis(0.146), 0.08)))
}
kwga_test_grid <- list(se = seq(0.70, 0.975, by = 0.025), sp = seq(0.70, 0.975, by = 0.025))

testthat::test_that("clamp_scoring = FALSE penalises pairs that imply a negative prevalence", {
  draws <- kwga_test_draws()
  clamped   <- mcma_kwga(draws, se_grid = kwga_test_grid$se, sp_grid = kwga_test_grid$sp)
  unclamped <- mcma_kwga(draws, se_grid = kwga_test_grid$se, sp_grid = kwga_test_grid$sp,
                         clamp_scoring = FALSE)
  testthat::expect_true(clamped$clamp_scoring)
  testthat::expect_false(unclamped$clamp_scoring)
  testthat::expect_true("log_prior" %in% names(unclamped$grid))

  # Per-cell display summaries are always computed from clamped values
  key <- function(g) g[order(g$se, g$sp), c("se", "sp", "corrected_mean", "corrected_ci_lb")]
  testthat::expect_equal(key(clamped$grid), key(unclamped$grid))

  # Cells whose correction clamps to zero keep most of the weight under clamped
  # scoring and lose the majority of it under unclamped scoring
  low_sp <- function(g) sum(g$weight[g$sp <= 0.85])
  testthat::expect_gt(low_sp(clamped$grid), 0.3)
  testthat::expect_lt(low_sp(unclamped$grid), 0.3)
  testthat::expect_lt(low_sp(unclamped$grid), 0.5 * low_sp(clamped$grid))
  testthat::expect_output(print(unclamped), "unclamped correction")
  testthat::expect_output(print(clamped), "clamped correction")
})

testthat::test_that("resample = 'auto' follows clamp_scoring and joint draws are conditioned", {
  draws <- kwga_test_draws()
  clamped   <- mcma_kwga(draws, se_grid = kwga_test_grid$se, sp_grid = kwga_test_grid$sp)
  unclamped <- mcma_kwga(draws, se_grid = kwga_test_grid$se, sp_grid = kwga_test_grid$sp,
                         clamp_scoring = FALSE)
  ind <- mcma_kwga_prevalence(clamped, seed = 1)
  jnt <- mcma_kwga_prevalence(unclamped, seed = 1)
  testthat::expect_equal(ind$resample, "independent")
  testthat::expect_equal(jnt$resample, "joint")
  testthat::expect_length(jnt$draws, 4000)
  testthat::expect_true(all(jnt$draws >= 0 & jnt$draws <= 1))

  # The original construction pins a large share of draws at exactly zero; the
  # corrected one does not, and its interval sits near the gold-standard posterior
  testthat::expect_gt(mean(ind$draws == 0), 0.2)
  testthat::expect_lt(mean(jnt$draws == 0), 0.05)
  testthat::expect_gt(jnt$summary$ci_lb, 0)
  testthat::expect_lt(abs(jnt$summary$mean - mean(draws$gold_draws)), 0.02)

  # Explicit overrides work in both directions and are reproducible with a seed
  forced <- mcma_kwga_prevalence(clamped, seed = 1, resample = "joint")
  testthat::expect_equal(forced$resample, "joint")
  testthat::expect_identical(mcma_kwga_prevalence(unclamped, seed = 7)$draws,
                             mcma_kwga_prevalence(unclamped, seed = 7)$draws)
  testthat::expect_equal(mcma_kwga_prevalence(unclamped, seed = 1, resample = "independent")$resample,
                         "independent")

  # Joint resampling needs paired gold draws
  testthat::expect_error(
    mcma_kwga_prevalence(unclamped, fit_comparison = list(screen_draws = rep(.3, 3))),
    "gold_draws")
})
