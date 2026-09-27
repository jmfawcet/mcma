# Exercise brms' real cache validation, replacing only compilation/sampling.
# This keeps the regression checks fast without duplicating its cache rules.

cache_test_data <- function() {
  data.frame(y = c(3, 5, 12, 14, 18, 10), n = 100, es_id = 1:6,
             measure_id = rep(c("A", "gold"), 3),
             is_gold = rep(c(FALSE, TRUE), 3))
}

cache_test_priors <- function() {
  mcma_priors(c("A", "gold"), c(.85, NA), c(.9, NA),
              is_gold = c(FALSE, TRUE))
}

cache_test_sampler <- function(env = parent.frame()) {
  state <- new.env(parent = emptyenv())
  state$runs <- 0L
  testthat::local_mocked_bindings(
    compile_model = function(...) list(),
    fit_model = function(...) {
      state$runs <- state$runs + 1L
      NULL
    }, .package = "brms", .env = env)
  state
}

cache_test_fit <- function(kind, file, ...) {
  args <- list(data = cache_test_data(), prev_re = NULL, backend = "rstan",
               file = file, rename = FALSE, refresh = 0)
  fitter <- switch(kind, joint = mcma_fit_joint, comparison = mcma_fit_comparison,
                   mcma_fit)
  if (kind == "comparison") {
    args$gold_prev_center <- .06
    args$scr_prev_center <- .2
  } else {
    args$prev_center <- .06
  }
  if (kind %in% c("corrected", "joint")) args$sesp_re <- NULL
  if (kind == "corrected") args$priors <- cache_test_priors()
  # Preserve explicit NULL overrides, such as removing random effects.
  overrides <- list(...)
  args[names(overrides)] <- overrides
  do.call(fitter, args)
}

testthat::test_that("unchanged caches preserve metadata and are not rewritten", {
  state <- cache_test_sampler()
  for (kind in c("naive", "corrected", "joint", "comparison")) {
    path <- tempfile(fileext = ".rds")
    first <- cache_test_fit(kind, path)
    runs <- state$runs
    Sys.setFileTime(path, as.POSIXct("2000-01-01", tz = "UTC"))
    before <- list(hash = tools::md5sum(path), mtime = file.info(path)$mtime)
    second <- cache_test_fit(kind, path)
    testthat::expect_identical(state$runs, runs)
    testthat::expect_identical(attr(second, "mcma_config"), attr(first, "mcma_config"))
    testthat::expect_identical(attr(second, "mcma_priors"), attr(first, "mcma_priors"))
    testthat::expect_identical(second$prior, first$prior)
    testthat::expect_identical(tools::md5sum(path), before$hash)
    testthat::expect_identical(file.info(path)$mtime, before$mtime)
  }
})

testthat::test_that("changed data priors and model structure refit cached models", {
  state <- cache_test_sampler()
  for (kind in c("naive", "corrected", "joint", "comparison")) {
    changes <- list(
      list(data = transform(cache_test_data(), y = y + 1)),
      if (kind == "comparison") list(gold_prev_center = .12) else list(prev_center = .12),
      list(prev_re = ~ (1 | es_id))
    )
    for (change in changes) {
      path <- tempfile(fileext = ".rds")
      cache_test_fit(kind, path)
      runs <- state$runs
      updated <- do.call(cache_test_fit, c(list(kind = kind, file = path), change))
      testthat::expect_identical(state$runs, runs + 1L)
      testthat::expect_identical(attr(readRDS(path), "mcma_config"),
                                attr(updated, "mcma_config"))
    }
  }
  # Accuracy-prior uncertainty also matters even at the same Se/Sp centres.
  path <- tempfile(fileext = ".rds")
  cache_test_fit("corrected", path)
  runs <- state$runs
  pri <- update(cache_test_priors(), kappa = 50)
  updated <- cache_test_fit("corrected", path, priors = pri)
  testthat::expect_identical(state$runs, runs + 1L)
  testthat::expect_identical(attr(updated, "mcma_priors"), pri)
})

testthat::test_that("never preserves stored metadata and cache bytes", {
  state <- cache_test_sampler()
  for (kind in c("naive", "corrected", "joint", "comparison")) {
    path <- tempfile(fileext = ".RDS")
    first <- cache_test_fit(kind, path)
    runs <- state$runs
    Sys.setFileTime(path, as.POSIXct("2000-01-01", tz = "UTC"))
    before <- list(hash = tools::md5sum(path), mtime = file.info(path)$mtime)
    change <- if (kind == "comparison") list(gold_prev_center = .12) else list(prev_center = .12)
    if (kind %in% c("corrected", "joint")) change$bounded <- FALSE
    second <- do.call(cache_test_fit, c(list(kind = kind, file = path,
      file_refit = "never", data = transform(cache_test_data(), y = y + 1)), change))
    testthat::expect_identical(state$runs, runs)
    testthat::expect_identical(attr(second, "mcma_config"), attr(first, "mcma_config"))
    testthat::expect_identical(attr(second, "mcma_priors"), attr(first, "mcma_priors"))
    testthat::expect_identical(second$prior, first$prior)
    testthat::expect_identical(tools::md5sum(path), before$hash)
    testthat::expect_identical(file.info(path)$mtime, before$mtime)
  }
})

testthat::test_that("always refits while sampling-only changes do not", {
  state <- cache_test_sampler()
  path <- tempfile(fileext = ".rds")
  cache_test_fit("naive", path)
  testthat::expect_identical(state$runs, 1L)
  cache_test_fit("naive", path, iter = 200, warmup = 100, adapt_delta = .9)
  testthat::expect_identical(state$runs, 1L)
  cache_test_fit("naive", path, file_refit = "always", iter = 200, warmup = 100)
  testthat::expect_identical(state$runs, 2L)
})

testthat::test_that("legacy caches require validation before annotation", {
  state <- cache_test_sampler()
  path <- tempfile(fileext = ".rds")
  first <- cache_test_fit("corrected", path)
  legacy <- first
  attr(legacy, "mcma_config") <- NULL
  attr(legacy, "mcma_priors") <- NULL
  saveRDS(legacy, path)
  before <- tools::md5sum(path)
  testthat::expect_error(cache_test_fit("corrected", path, file_refit = "never"),
                        "no mcma_config")
  testthat::expect_null(attr(readRDS(path), "mcma_config"))
  testthat::expect_identical(tools::md5sum(path), before)
  validated <- cache_test_fit("corrected", path)
  testthat::expect_identical(state$runs, 1L)
  testthat::expect_identical(attr(validated, "mcma_config"), attr(first, "mcma_config"))
  testthat::expect_identical(attr(readRDS(path), "mcma_priors"), attr(first, "mcma_priors"))

  saveRDS(legacy, path)
  updated <- cache_test_fit("corrected", path, bounded = FALSE)
  testthat::expect_identical(state$runs, 2L)
  testthat::expect_false(attr(updated, "mcma_config")$bounded)
})

testthat::test_that("cached bounded accuracy draws retain their correct transformation", {
  cache_test_sampler()
  path <- tempfile(fileext = ".rds")
  cache_test_fit("corrected", path)
  cached <- cache_test_fit("corrected", path, bounded = FALSE, file_refit = "never")
  dr <- posterior::as_draws_df(data.frame(b_Se_measure_idA = c(1, 2, 3),
                                         b_Sp_measure_idA = c(2, 3, 4)))
  testthat::local_mocked_bindings(as_draws_df = function(...) dr, .package = "posterior")
  accuracy <- extract_sesp(cached)
  testthat::expect_equal(accuracy$se_mean, mean(.5 + .5 * plogis(c(1, 2, 3))))
  testthat::expect_equal(accuracy$sp_mean, mean(.5 + .5 * plogis(c(2, 3, 4))))
})

testthat::test_that("both sensitivity grids validate caches and honour rerun", {
  state <- cache_test_sampler()
  dr <- posterior::as_draws_df(data.frame(b_pi_Intercept = c(-3, -2, -1),
    b_pi_is_goldTRUE = c(-3, -2, -1), b_pi_is_goldFALSE = c(-2, -1, 0)))
  testthat::local_mocked_bindings(as_draws_df = function(...) dr,
    rhat = function(...) 1, .package = "posterior")
  testthat::local_mocked_bindings(nuts_params = function(...) {
    data.frame(Parameter = "divergent__", Value = 0)
  }, .package = "brms")

  for (fitter in list(mcma_sensitivity, mcma_sensitivity_comparison)) {
    args <- list(data = cache_test_data(), priors = cache_test_priors(),
      se_grid = .8, sp_grid = .9, prev_center = .06, prev_re = NULL,
      sesp_re = NULL, model_dir = tempfile(), backend = "rstan", rename = FALSE)
    do.call(fitter, args)
    runs <- state$runs
    path <- list.files(args$model_dir, full.names = TRUE)
    before <- tools::md5sum(path)
    do.call(fitter, args)
    testthat::expect_identical(state$runs, runs)
    testthat::expect_identical(tools::md5sum(path), before)
    args$data$y <- args$data$y + 1
    do.call(fitter, args)
    testthat::expect_identical(state$runs, runs + 1L)
    do.call(fitter, c(args, list(rerun = TRUE)))
    testthat::expect_identical(state$runs, runs + 2L)
  }
})
