# Small saved records exercise repair without running a sampler.
make_repair_fixture <- function() {
  dir <- tempfile("mcma-repair-")
  cor_dir <- file.path(dir, "backups")
  dir.create(cor_dir, recursive = TRUE)
  file <- file.path(dir, "Sim1_1.rds")
  backup <- file.path(cor_dir, basename(file))
  original <- list(data = data.frame(y = 3, n = 10), result = list(value = 0.2))
  previous <- list(data = original$data, result = list(value = 0.1))
  saveRDS(original, file)
  saveRDS(previous, backup)
  list(dir = dir, cor_dir = cor_dir, file = file, backup = backup,
       original = original, previous = previous)
}

testthat::test_that("repeated failed refits preserve the current result and backup", {
  for (failure in c("error", "NULL")) {
    x <- make_repair_fixture()
    attempts <- 0L
    fitter <- function(...) {
      attempts <<- attempts + 1L
      if (failure == "error") stop("controlled fitting failure")
      NULL
    }

    suppressWarnings(mcma_sim_repair(
      fit_fn = fitter, dir = x$dir, cor_dir = x$cor_dir, cycles = 2,
      bad_criteria = function(x) TRUE
    ))

    testthat::expect_identical(attempts, 2L)
    testthat::expect_identical(readRDS(x$file), x$original)
    testthat::expect_identical(readRDS(x$backup), x$previous)
    unlink(x$dir, recursive = TRUE)
  }
})

testthat::test_that("a later successful refit replaces the result and preserves its input", {
  x <- make_repair_fixture()
  attempts <- 0L
  repaired <- list(value = 0.3)
  fitter <- function(data, ...) {
    attempts <<- attempts + 1L
    testthat::expect_identical(data, x$original$data)
    if (attempts == 1L) return(NULL)
    repaired
  }

  mcma_sim_repair(
    fit_fn = fitter, dir = x$dir, cor_dir = x$cor_dir, cycles = 2,
    bad_criteria = function(x) TRUE
  )

  testthat::expect_identical(attempts, 2L)
  testthat::expect_identical(readRDS(x$file),
                            list(data = x$original$data, result = repaired))
  testthat::expect_identical(readRDS(x$backup), x$original)
  testthat::expect_length(list.files(x$dir, pattern = "^\\.mcma_repair_", all.files = TRUE), 0L)
  unlink(x$dir, recursive = TRUE)
})

testthat::test_that("storage failures leave the current result readable", {
  for (failure in c("write", "backup", "rename")) {
    x <- make_repair_fixture()

    testthat::with_mocked_bindings({
      testthat::expect_error(mcma_sim_repair(
        fit_fn = function(...) list(value = 0.3),
        dir = x$dir, cor_dir = x$cor_dir, cycles = 1,
        bad_criteria = function(x) TRUE
      ), "controlled write failure|Could not back up|Could not replace")
    },
    saveRDS = if (failure == "write") function(object, file, ...) {
      writeLines("incomplete replacement", file)
      stop("controlled write failure")
    } else base::saveRDS,
    file.copy = if (failure == "backup") function(...) FALSE else base::file.copy,
    file.rename = if (failure == "rename") function(...) FALSE else base::file.rename,
    .package = "base")

    testthat::expect_identical(readRDS(x$file), x$original)
    testthat::expect_identical(readRDS(x$backup),
                              if (failure == "rename") x$original else x$previous)
    testthat::expect_length(list.files(x$dir, pattern = "^\\.mcma_repair_", all.files = TRUE), 0L)
    unlink(x$dir, recursive = TRUE)
  }
})

