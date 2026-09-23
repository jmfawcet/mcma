# Regression checks for instrument naming, caches, grouping and diagnostics.
# Model formulas are validated by brms without running Stan sampling.

review_data <- function() {
  data.frame(y = c(3, 5, 12, 14, 18, 10), n = 100, es_id = 1:6,
             measure_id = rep(c("A", "gold"), 3),
             is_gold = rep(c(FALSE, TRUE), 3))
}

review_priors <- function() {
  mcma_priors(c("A", "gold"), c(.85, NA), c(.9, NA),
              is_gold = c(FALSE, TRUE))
}

validate_review_model <- function(fit) {
  args <- fit[c("formula", "data", "family", "prior")]
  testthat::expect_type(do.call(brms::make_standata, args), "list")
  testthat::expect_type(do.call(brms::make_stancode, args), "character")
}

testthat::test_that("brms caches retain mcma metadata after a direct reload", {
  testthat::local_mocked_bindings(brm = function(..., file = NULL) {
    args <- list(...)
    path <- if (endsWith(file, ".rds")) file else paste0(file, ".rds")
    fit <- structure(list(file = path, data = args$data), class = "brmsfit")
    # brms writes its object before the package wrapper annotates it.
    saveRDS(fit, path, compress = FALSE)
    fit
  }, .package = "brms")
  dat <- review_data()
  fits <- list(
    mcma_fit(dat, prev_center = .06, file = tempfile(), file_compress = FALSE),
    mcma_fit(dat, review_priors(), prev_center = .06, bounded = TRUE,
              file = tempfile(fileext = ".rds"), file_compress = FALSE),
    mcma_fit_joint(dat, prev_center = .06, c = .3, o = .02,
                    file = tempfile(), file_compress = FALSE),
    mcma_fit_comparison(dat, gold_prev_center = .06, scr_prev_center = .2,
                         file = tempfile(), file_compress = FALSE)
  )
  for (fit in fits) {
    cached <- readRDS(fit$file)
    testthat::expect_identical(attr(cached, "mcma_config"), attr(fit, "mcma_config"))
    testthat::expect_identical(attr(cached, "mcma_priors"), attr(fit, "mcma_priors"))
  }
  cached <- readRDS(fits[[3]]$file)
  testthat::expect_equal(extract_bias(cached)$mean, c(.3, .02))
  dr <- posterior::as_draws_df(data.frame(
    b_Se_acc_groupscreen = c(1, 2, 3), b_Sp_acc_groupscreen = c(2, 3, 4)))
  testthat::with_mocked_bindings({
    testthat::expect_equal(extract_sesp(cached)$se_mean, mean(.5 + .5 * plogis(c(1, 2, 3))))
  }, as_draws_df = function(...) dr, .package = "posterior")

  # Loading an already annotated cache should not rewrite the large file.
  Sys.setFileTime(cached$file, as.POSIXct("2000-01-01", tz = "UTC"))
  before <- file.info(cached$file)$mtime
  mcma:::.mcma_finalize_fit(cached, attr(cached, "mcma_config"), file = cached$file)
  testthat::expect_identical(file.info(cached$file)$mtime, before)
})

testthat::test_that("instrument labels and unused prior rows are handled consistently", {
  testthat::local_mocked_bindings(brm = function(...) list(...), .package = "brms")
  dat <- review_data()
  dat$instrument <- rep(c("PHQ-9", "CES D", "gold"), 2)
  pri <- mcma_priors(c("unused", "CES D", "gold", "PHQ-9"),
                     c(.75, .8, NA, .85), c(.8, .9, NA, .95),
                     is_gold = c(FALSE, FALSE, TRUE, FALSE))
  for (bounded in c(FALSE, TRUE)) {
    testthat::expect_warning(
      fit <- mcma_fit(dat, pri, prev_center = .06, measure_col = "instrument",
                      bounded = bounded), "Ignoring priors.*unused")
    validate_review_model(fit)
    testthat::expect_setequal(attr(fit, "mcma_priors")$measure_id,
                              c("PHQ-9", "CES D", "gold"))
    testthat::expect_true(all(c("instrumentPHQM9", "instrumentCESD") %in% fit$prior$coef))
  }
  dr <- posterior::as_draws_df(data.frame(
    b_Se_instrumentPHQM9 = c(1, 2, 3), b_Sp_instrumentPHQM9 = c(2, 3, 4),
    b_Se_instrumentCESD = c(.5, 1, 1.5), b_Sp_instrumentCESD = c(1, 1.5, 2),
    b_Se_instrumentgold = 10, b_Sp_instrumentgold = 10))
  testthat::with_mocked_bindings({
    testthat::expect_setequal(extract_sesp(fit, measure_col = "instrument")$measure_id,
                              c("PHQ-9", "CES D"))
    testthat::expect_setequal(names(extract_sesp(fit, measure_col = "instrument", summary = FALSE)),
                              c("PHQ-9", "CES D"))
  }, as_draws_df = function(...) dr, .package = "posterior")
  testthat::expect_error(mcma_fit(dat, review_priors(), prev_center = .06,
                                  measure_col = "instrument"), "Missing Se/Sp priors")
  # The exported prior builder also supports a one-instrument prior table.
  single <- as_brms_prior(mcma_priors("PHQ-9", .85, .9), prev_center = .06)
  testthat::expect_true("measure_idPHQM9" %in% single$coef)
})

testthat::test_that("study_col controls defaults but preserves explicit formulas", {
  testthat::local_mocked_bindings(brm = function(...) list(...), .package = "brms")
  dat <- review_data()
  names(dat)[names(dat) == "es_id"] <- "study"
  joint <- mcma_fit_joint(dat, prev_center = .06, study_col = "study")
  comparison <- mcma_fit_comparison(dat, gold_prev_center = .06, scr_prev_center = .2,
                                     study_col = "study")
  validate_review_model(joint)
  validate_review_model(comparison)
  testthat::expect_true("study" %in% all.vars(joint$formula$pforms$pi))
  testthat::expect_true("study" %in% all.vars(joint$formula$pforms$Se))
  testthat::expect_true("study" %in% all.vars(comparison$formula))
  dat$report <- rep(1:3, each = 2)
  explicit <- mcma_fit_joint(dat, prev_center = .06, study_col = "study",
                              prev_re = ~(1 | report), sesp_re = NULL)
  validate_review_model(explicit)
  testthat::expect_true("report" %in% all.vars(explicit$formula$pforms$pi))
  testthat::expect_false("study" %in% all.vars(explicit$formula$pforms$pi))
  testthat::expect_false("study" %in% all.vars(explicit$formula$pforms$Se))
})

testthat::test_that("sensitivity grids retain precise cache names and subset priors", {
  paths <- character()
  calls <- list()
  dr <- posterior::as_draws_df(data.frame(
    b_pi_Intercept = c(-3, -2, -1),
    b_pi_is_goldTRUE = c(-3, -2, -1), b_pi_is_goldFALSE = c(-2, -1, 0)))
  testthat::local_mocked_bindings(brm = function(...) {
    args <- list(...)
    validate_review_model(args)
    paths <<- c(paths, args$file)
    calls[[length(calls) + 1L]] <<- args
    dr
  }, nuts_params = function(...) data.frame(Parameter = "divergent__", Value = 0),
  .package = "brms")
  testthat::local_mocked_bindings(rhat = function(...) 1, .package = "posterior")
  dat <- review_data()
  dat$study <- dat$es_id
  pri <- mcma_priors(c("A", "gold", "unused"),c(.85, NA, .8),c(.9, NA, .9),
                     is_gold = c(FALSE, TRUE, FALSE))
  testthat::expect_warning(
    ordinary <- mcma_sensitivity(dat, pri, se_grid = c(.8, .804), sp_grid = .9,
                                 prev_center = .06, model_dir = tempfile()), "Ignoring priors")
  testthat::expect_equal(length(unique(paths)), 2L)
  testthat::expect_equal(basename(paths), c("Se0.80_Sp0.90", "Se0.804_Sp0.90"))
  testthat::expect_true(all(is.finite(ordinary$results$prev_mean)))
  paths <- character()
  names(dat)[names(dat) == "measure_id"] <- "instrument"
  dat$es_id <- NULL
  testthat::expect_warning(
    comparison <- mcma_sensitivity_comparison(dat, pri, se_grid = c(.8, .804), sp_grid = .9,
      prev_center = .06, study_col = "study", measure_col = "instrument", model_dir = tempfile()),
    "Ignoring priors")
  testthat::expect_equal(length(unique(paths)), 2L)
  testthat::expect_true(all(is.finite(comparison$results$diff_mean)))
  testthat::expect_true("study" %in% all.vars(calls[[3]]$formula$pforms$pi))
  testthat::expect_equal(mcma:::.mcma_grid_key(.8, .9), "Se0.80_Sp0.90")
  testthat::expect_false(identical(mcma:::.mcma_grid_key(.800001, .9),
                                  mcma:::.mcma_grid_key(.800002, .9)))
})

testthat::test_that("convergence summaries retain mixed and repeated model names", {
  testthat::local_mocked_bindings(
    as_draws = function(x, ...) x$id,
    summarise_draws = function(x, ...) data.frame(rhat = 1 + x/1000, ess_bulk = 1000, ess_tail = 900),
    .package = "posterior")
  testthat::local_mocked_bindings(nuts_params = function(...) {
    data.frame(Parameter = c("divergent__", "treedepth__"), Value = c(0, 5))
  }, .package = "brms")
  one <- structure(list(id = 1), class = "brmsfit")
  two <- structure(list(id = 2), class = "brmsfit")
  out <- mcma_convergence(first = one, two)
  testthat::expect_equal(out$model, c("first", "model_2"))
  testthat::expect_equal(out$max_rhat, c(1.001, 1.002))
  repeated <- mcma_convergence(same = one, same = two)
  testthat::expect_equal(repeated$max_rhat, c(1.001, 1.002))
  testthat::expect_equal(mcma_convergence(one, two)$model, c("model_1", "model_2"))
})

testthat::test_that("KWGA uses an explicitly supplied comparison without changing weights", {
  original <- list(gold_draws = c(.04, .05, .06), screen_draws = c(.20, .21, .22))
  kwga <- mcma_kwga(original, se_grid = .85, sp_grid = .9, clamp_scoring = TRUE)
  stored <- mcma_kwga_prevalence(kwga, seed = 1)
  same <- mcma_kwga_prevalence(kwga, fit_comparison = original, seed = 1)
  testthat::expect_identical(stored, same)
  replaced <- mcma_kwga_prevalence(kwga, fit_comparison = list(screen_draws = rep(.3, 3)), seed = 1)
  testthat::expect_equal(replaced$summary$mean, (.3 + .9 - 1) / (.85 + .9 - 1))
  testthat::expect_identical(replaced$weights, kwga$grid)
  testthat::expect_gt(replaced$summary$mean, stored$summary$mean)
  testthat::expect_error(mcma_kwga_prevalence(kwga, fit_comparison = "not a fit"), "fit_comparison")
  testthat::expect_error(mcma_kwga_prevalence(kwga, fit_comparison = list()), "screen_draws")
  dr <- posterior::as_draws_df(data.frame(
    b_diagnosticTRUE = qlogis(c(.04, .05, .06)), b_diagnosticFALSE = qlogis(rep(.3, 3))))
  fit <- structure(list(), class = "brmsfit", mcma_config = list(gold_column = "diagnostic"))
  testthat::with_mocked_bindings({
    result <- mcma_kwga_prevalence(kwga, fit_comparison = fit, seed = 1)
    testthat::expect_equal(result$draws, replaced$draws)
  }, as_draws_df = function(...) dr, .package = "posterior")
})
