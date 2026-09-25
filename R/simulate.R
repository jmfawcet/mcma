# Simulation functions
# ====================

#' Generate simulated prevalence data
#'
#' Generates simulated datasets for power analysis and methods validation.
#' When arguments are vectors, creates a grid of scenarios. Each scenario
#' is replicated `n_reps` times.
#'
#' @param k Number of studies per dataset.
#' @param di Number of distinct screening measures.
#' @param mean_se Population mean Se (probability scale). Can be a vector for
#'   grid construction.
#' @param mean_sp Population mean Sp. Can be a vector.
#' @param sd_se_logit Between-measure SD of Se on logit scale.
#' @param sd_sp_logit Between-measure SD of Sp on logit scale.
#' @param rho_sesp Correlation between logit(Se) and logit(Sp).
#' @param omega_se Study-level jitter SD for Se on logit scale.
#' @param omega_sp Study-level jitter SD for Sp on logit scale.
#' @param mu_prevalence True population prevalence. Can be a vector.
#' @param tau_logit Between-study heterogeneity SD on logit scale.
#' @param prop_gold Proportion of gold-standard studies (0 = none). Can be
#'   a vector.
#' @param moderator If not NULL, a list with `slope` and `sd`.
#' @param meanlog_n Log-normal mean for sample sizes.
#' @param sdlog_n Log-normal SD for sample sizes.
#' @param add_const_n Constant added to sample sizes.
#' @param floor_prob Minimum Se/Sp on probability scale.
#' @param ceiling_prob Maximum Se/Sp on probability scale.
#' @param n_reps Number of replications per scenario.
#' @param seed Random seed for reproducibility.
#' @return An S3 object of class `mcma_sim_data`.
#' @export
mcma_sim_data <- function(k             = 20,
                          di            = 10,
                          mean_se       = 0.85,
                          mean_sp       = 0.80,
                          sd_se_logit   = 0.2,
                          sd_sp_logit   = 0.2,
                          rho_sesp      = -0.5,
                          omega_se      = 0.2,
                          omega_sp      = 0.2,
                          mu_prevalence = 0.06,
                          tau_logit     = 1.0,
                          prop_gold     = 0,
                          moderator     = NULL,
                          meanlog_n     = 5.50,
                          sdlog_n       = 0.55,
                          add_const_n   = 30,
                          floor_prob    = 0.5,
                          ceiling_prob  = 0.999,
                          n_reps        = 1,
                          seed          = NULL) {

  # Generate repeated datasets for each combination of design settings.
  # Setting a seed makes the sequence of random draws reproducible.

  if (!is.null(seed)) set.seed(seed)

  # Each vector-valued design argument contributes a dimension to the grid;
  # one row is one experimental condition.
  # Build scenario grid from vectorized arguments
  scenario_grid <- tidyr::expand_grid(
    k             = k,
    di            = di,
    mean_se       = mean_se,
    mean_sp       = mean_sp,
    mu_prevalence = mu_prevalence,
    prop_gold     = prop_gold,
    tau_logit     = tau_logit
  )
  scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))

  # Generate datasets
  all_data <- vector("list", nrow(scenario_grid))

  for (s in seq_len(nrow(scenario_grid))) {
    sc <- scenario_grid[s, ]
    reps <- vector("list", n_reps)

    # Generate independent replications under the current condition while
    # retaining all generation settings.
    for (r in seq_len(n_reps)) {
      reps[[r]] <- .generate_one_dataset(
        k           = sc$k,
        di          = sc$di,
        mean_se     = sc$mean_se,
        mean_sp     = sc$mean_sp,
        sd_se_logit = sd_se_logit,
        sd_sp_logit = sd_sp_logit,
        rho_sesp    = rho_sesp,
        omega_se    = omega_se,
        omega_sp    = omega_sp,
        mu_prevalence = sc$mu_prevalence,
        tau_logit   = sc$tau_logit,
        prop_gold   = sc$prop_gold,
        moderator   = moderator,
        meanlog_n   = meanlog_n,
        sdlog_n     = sdlog_n,
        add_const_n = add_const_n,
        floor_prob  = floor_prob,
        ceiling_prob = ceiling_prob
      )
    }

    all_data[[s]] <- reps
  }

  # Keep the design table, generated datasets, and shared settings together
  # for reproducible fitting and inspection.
  out <- list(
    scenarios = scenario_grid,
    data      = all_data,
    n_reps    = n_reps,
    params    = list(
      sd_se_logit = sd_se_logit, sd_sp_logit = sd_sp_logit,
      rho_sesp = rho_sesp, omega_se = omega_se, omega_sp = omega_sp,
      meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
      floor_prob = floor_prob, ceiling_prob = ceiling_prob,
      moderator = moderator
    )
  )
  class(out) <- c("mcma_sim_data", "list")
  out
}


#' Run simulation fitting
#'
#' Fits models to each replicate in a simulation object. Results are saved
#' as per-replicate RDS files.
#'
#' @param sim_data An `mcma_sim_data` object.
#' @param fit_fn Fitting function (e.g., `mcma_fit`).
#' @param dir Output directory for RDS files.
#' @param skip_num Number of initial replicates to skip.
#' @param ... Additional arguments passed to `fit_fn`.
#' @export
mcma_sim_run <- function(sim_data, fit_fn, dir, skip_num = 0, ...) {

  # Fit each generated dataset and save it with its scenario metadata.
  # Existing files are skipped so a long simulation can resume after
  # interruption.

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  n_scenarios <- nrow(sim_data$scenarios)
  n_reps      <- sim_data$n_reps

  for (s in seq_len(n_scenarios)) {
    for (r in seq_len(n_reps)) {
      if (s == 1 && r <= skip_num) next

      out_file <- file.path(dir, sprintf("Sim%d_%d.rds", s, r))
      if (file.exists(out_file)) {
        message(sprintf("Skipping Sim%d_%d (exists)", s, r))
        next
      }

      dat <- sim_data$data[[s]][[r]]
      message(sprintf("Fitting Sim%d_%d...", s, r))

      # Capture fitting errors as NULL results so later replications can
      # continue and failed runs can be revisited.
      result <- tryCatch({
        fit_fn(data = dat, ...)
      }, error = function(e) {
        rlang::warn(sprintf("Sim%d_%d failed: %s", s, r, e$message))
        NULL
      })

      # Save both the result and its exact input data, avoiding the need to
      # regenerate random data during a repair.
      saveRDS(
        list(
          scenario_id = s,
          rep_id      = r,
          factors     = as.list(sim_data$scenarios[s, ]),
          result      = result,
          data        = dat
        ),
        file = out_file
      )
    }
  }
}


#' Re-fit failed simulation replicates
#'
#' Identifies replicates with convergence issues and re-fits with tighter
#' sampler settings.
#'
#' @param sim_data An `mcma_sim_data` object.
#' @param fit_fn Fitting function.
#' @param dir Directory containing simulation results.
#' @param cor_dir Directory for storing replaced (bad) files.
#' @param iter Iterations for re-fitting.
#' @param adapt_delta Acceptance rate for re-fitting.
#' @param step_size Step size for re-fitting.
#' @param cycles Number of repair passes.
#' @param bad_criteria Function taking a result list, returning TRUE if bad.
#'   NULL checks failed fits, Rhat, effective sample sizes, and divergences.
#' @param rhat_limit Largest acceptable Rhat (default 1.01).
#' @param min_ess Smallest acceptable bulk or tail effective sample size
#'   (default 400). Set to zero to disable this check.
#' @param max_divergences Largest acceptable number of divergent transitions.
#' @param ... Additional arguments passed to `fit_fn`.
#' @export
mcma_sim_repair <- function(sim_data,
                            fit_fn,
                            dir,
                            cor_dir     = paste0(dir, "_corrected/"),
                            iter        = 12000,
                            adapt_delta = 0.9999,
                            step_size   = 0.002,
                            cycles      = 3,
                            bad_criteria = NULL,
                            rhat_limit  = 1.01,
                            min_ess     = 400,
                            max_divergences = 0,
                            ...) {

  # Apply the same convergence checks as simulation summaries. A caller
  # can replace these checks with a criterion suited to a custom fitter.

  if (is.null(bad_criteria)) {
    bad_criteria <- function(x) {
      nzchar(.mcma_sim_bad_reason(x, rhat_limit, min_ess, max_divergences))
    }
  }

  if (!dir.exists(cor_dir)) dir.create(cor_dir, recursive = TRUE)

  for (cycle in seq_len(cycles)) {
    files <- list.files(dir, pattern = "^Sim\\d+_\\d+\\.rds$", full.names = TRUE)
    bad_files <- character(0)

    # Inspect every saved record with the selected criterion and collect the
    # files that need another fit.
    for (f in files) {
      res <- readRDS(f)
      if (bad_criteria(res)) bad_files <- c(bad_files, f)
    }

    if (length(bad_files) == 0) {
      message(sprintf("Cycle %d: no bad replicates found. Done.", cycle))
      break
    }

    message(sprintf("Cycle %d: re-fitting %d bad replicates", cycle, length(bad_files)))

    for (f in bad_files) {

      res <- readRDS(f)
      dat <- res$data

      # Forward all requested tuning settings. The package fitters map
      # step_size to the control name required by the selected Stan backend.
      new_result <- tryCatch({
        fit_fn(data = dat, iter = iter, adapt_delta = adapt_delta,
               step_size = step_size, ...)
      }, error = function(e) {
        rlang::warn(sprintf("Re-fit of %s failed: %s", basename(f), e$message))
        NULL
      })

      # Keep the current result and backup after a failed fit, so the next
      # cycle can still find and retry this replicate.
      if (is.null(new_result)) next

      res$result <- new_result

      # Finish writing beside the original before backing it up and replacing
      # it. A failed write or rename leaves the current result in place.
      tmp_file <- tempfile(pattern = ".mcma_repair_", tmpdir = dirname(f))
      tryCatch({
        saveRDS(res, file = tmp_file)
        if (!file.copy(f, file.path(cor_dir, basename(f)), overwrite = TRUE)) {
          stop("Could not back up ", f, "; the current result was kept.")
        }
        if (!file.rename(tmp_file, f)) {
          stop("Could not replace ", f, "; the current result was kept.")
        }
      }, finally = {
        if (file.exists(tmp_file)) unlink(tmp_file)
      })
    }
  }
}


#' Read simulation results
#'
#' Reads and combines all per-replicate RDS result files from a directory.
#'
#' @param dir Directory containing simulation results.
#' @param pattern File pattern to match.
#' @return A list of result objects.
#' @export
mcma_sim_read <- function(dir, pattern = "^Sim\\d+_\\d+\\.rds$") {
  # Read the numbered per-replicate files and retain their filenames as list
  # names, so individual records can be traced back to disk.
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    rlang::warn("No simulation result files found.")
    return(list())
  }

  # Return a named list rather than combining fitted objects into a single
  # data frame.
  results <- lapply(files, readRDS)
  names(results) <- basename(files)
  results
}


#' Summarise simulation results
#'
#' Groups replicate estimates by design factors, with bias, RMSE and coverage
#' when generation truth is available. Failed fits and extraction failures are
#' always recorded; extractable fits with poor convergence are retained only
#' when `exclude_bad = FALSE`.
#'
#' @param results List of simulation results from `mcma_sim_read()`.
#' @param group_vars Variables to group by; NULL uses `"scenario_id"`.
#' @param summary_type `"prevalence"` or `"moderator"`.
#' @param exclude_bad Logical; exclude fits flagged by the convergence criterion.
#' @param bad_criteria Optional function taking a saved replicate record and
#'   returning TRUE when it should be flagged. NULL uses the default checks.
#' @param rhat_limit Largest acceptable Rhat (default 1.01).
#' @param min_ess Smallest acceptable bulk or tail ESS (default 400);
#'   zero disables this check.
#' @param max_divergences Largest acceptable divergence count (default zero).
#' @return A list containing grouped `$summary`, all extractable `$replicates`
#'   (with an `included` flag), `$bad_reps` with reasons, and `$n_bad` counting
#'   replicate records rather than coefficients. Summary interval endpoints are
#'   averages of replicate endpoints, not an interval for the simulation mean.
#'   `rejection_rate` is the fraction of moderator intervals excluding zero;
#'   it estimates power for a nonzero slope and type-I error for a zero slope.
#' @export
mcma_sim_summarise <- function(results,
                               group_vars   = NULL,
                               summary_type = c("prevalence", "moderator"),
                               exclude_bad  = FALSE,
                               bad_criteria = NULL,
                               rhat_limit   = 1.01,
                               min_ess      = 400,
                               max_divergences = 0) {

  # Keep an audit record for every failed replicate before aggregating
  # estimates. Filtering must not erase the evidence that a fit failed.
  summary_type <- match.arg(summary_type)
  rows <- vector("list", length(results))
  bad <- vector("list", length(results))

  for (i in seq_along(results)) {
    res <- results[[i]]
    metadata <- tibble::tibble(scenario_id = res$scenario_id, rep_id = res$rep_id)
    for (nm in setdiff(names(res$factors), names(metadata))) {
      if (length(res$factors[[nm]]) == 1L) metadata[[nm]] <- res$factors[[nm]]
    }

    reason <- if (is.null(res$result) || is.null(bad_criteria)) {
      .mcma_sim_bad_reason(res, rhat_limit, min_ess, max_divergences)
    } else if (isTRUE(bad_criteria(res))) {
      "Flagged by bad_criteria"
    } else {
      ""
    }

    # Extract the requested estimand. Comparison prevalence keeps one row
    # per population; moderator summaries keep one row per coefficient.
    metrics <- NULL
    if (!is.null(res$result)) {
      metrics <- tryCatch({
        if (summary_type == "prevalence") {
          prev <- extract_prevalence(res$result)
          tau <- extract_tau(res$result)$mean
          out <- tibble::tibble(
            prev_mean = prev$mean, prev_median = prev$median,
            prev_sd = prev$sd, ci_lb = prev$l95, ci_ub = prev$u95,
            tau_mean = tau, estimate = prev$mean
          )
          if ("population" %in% names(prev)) out$population <- prev$population
          truth <- res$factors$mu_prevalence
          out$truth <- if (length(truth) == 1L) truth else NA_real_
          out$reject_zero <- NA
        } else {
          slope <- extract_slope(res$result)
          out <- tibble::tibble(
            variable = slope$variable, slope_mean = slope$mean,
            slope_median = slope$median, slope_sd = slope$sd,
            ci_lb = slope$ci_lb, ci_ub = slope$ci_ub,
            estimate = slope$mean, truth = NA_real_
          )

          # The package generator records the true x_mod slope in each
          # dataset. Other coefficients remain unscored when truth is unknown.
          truth <- unique(res$data$slope_true)
          if (length(truth) == 1L) out$truth[out$variable == "x_mod"] <- truth
          out$reject_zero <- out$ci_lb > 0 | out$ci_ub < 0
        }

        if (!nrow(out) || any(!is.finite(out$estimate)) ||
            any(!is.finite(out$ci_lb)) || any(!is.finite(out$ci_ub))) {
          stop("No finite estimates and intervals could be extracted.")
        }
        out
      }, error = function(e) {
        reason <<- paste(c(reason[nzchar(reason)], paste("Extraction failed:", conditionMessage(e))),
                         collapse = "; ")
        NULL
      })
    }

    if (nzchar(reason)) {
      bad[[i]] <- tibble::add_column(metadata, reason = reason)
    }
    if (!is.null(metrics)) {
      out <- cbind(metadata[rep(1L, nrow(metrics)), , drop = FALSE], metrics)
      out$bad <- nzchar(reason)
      out$included <- !(exclude_bad && nzchar(reason))
      rows[[i]] <- out
    }
  }

  # Preserve usable replicate rows for inspection, including excluded
  # ones. n_bad counts each saved record once even if it has several slopes.
  replicates <- dplyr::bind_rows(rows)
  bad_reps <- dplyr::bind_rows(bad)
  included <- if (nrow(replicates)) replicates[replicates$included, , drop = FALSE] else replicates
  if (!nrow(included)) {
    rlang::warn("No included results to summarise.")
    return(list(summary = tibble::tibble(), replicates = replicates,
                bad_reps = bad_reps, n_bad = nrow(bad_reps)))
  }

  if (is.null(group_vars)) group_vars <- "scenario_id"
  group_vars <- unique(c(group_vars, intersect(c("population", "variable"), names(included))))
  if (any(!group_vars %in% names(included))) {
    stop("Unknown group_vars: ", paste(setdiff(group_vars, names(included)), collapse = ", "), call. = FALSE)
  }

  # Bias and coverage compare each estimate with its own generation
  # truth before averaging, so grouping over different truths remains valid.
  included$error <- included$estimate - included$truth
  included$covered <- included$ci_lb <= included$truth & included$ci_ub >= included$truth
  included$ci_width <- included$ci_ub - included$ci_lb
  average_columns <- intersect(
    c("prev_mean", "prev_median", "prev_sd", "tau_mean",
      "slope_mean", "slope_median", "slope_sd", "ci_lb", "ci_ub", "ci_width"),
    names(included)
  )
  average_columns <- setdiff(average_columns, group_vars)

  grouped <- dplyr::group_by(included, dplyr::across(dplyr::all_of(group_vars)))
  summary_tbl <- dplyr::summarise(
    grouped,
    n_reps = dplyr::n(),
    dplyr::across(dplyr::all_of(average_columns), .mcma_mean_available),
    empirical_sd = stats::sd(.data$estimate),
    bias = .mcma_mean_available(.data$error),
    rmse = sqrt(.mcma_mean_available(.data$error^2)),
    coverage = .mcma_mean_available(.data$covered),
    rejection_rate = .mcma_mean_available(.data$reject_zero),
    .groups = "drop"
  )

  list(summary = summary_tbl, replicates = replicates,
       bad_reps = bad_reps, n_bad = nrow(bad_reps))
}

.mcma_mean_available <- function(x) {

  # Missing truth or unavailable diagnostics should stay missing rather
  # than becoming NaN or being interpreted as a zero error rate.
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

.mcma_sim_bad_reason <- function(record, rhat_limit, min_ess, max_divergences) {

  # Use one convergence policy for repair and summary. Missing diagnostics
  # cannot establish convergence, so flag them unless a custom criterion is used.
  fit <- record$result
  if (is.null(fit)) return("Fit failed (NULL result)")

  diagnostics <- tryCatch({
    if (inherits(fit, "brmsfit")) {
      mcma_convergence(model = fit)
    } else if (!is.null(fit$diagnostics)) {
      fit$diagnostics
    } else {
      fit
    }
  }, error = function(e) NULL)

  # Unsupported custom results have unavailable diagnostics; they must
  # not be mistaken for successful convergence or crash the repair loop.
  if (!is.list(diagnostics)) diagnostics <- list()

  # Accept both package diagnostic names and the shorter names used in
  # saved simulation summaries from custom fitting functions.
  value <- function(names) {
    for (nm in names) {
      x <- diagnostics[[nm]]
      if (is.numeric(x) && length(x) == 1L) return(x)
    }
    NA_real_
  }
  rh <- value(c("max_rhat", "maxrhat"))
  nd <- value(c("n_divergences", "ndiv"))
  bulk <- value(c("min_ess_bulk", "ess_bulk", "min_ess"))
  tail <- value(c("min_ess_tail", "ess_tail", "min_ess"))

  reasons <- character()
  if (!is.finite(rh)) reasons <- c(reasons, "Rhat unavailable")
  else if (rh > rhat_limit) reasons <- c(reasons, "Rhat exceeds limit")
  if (!is.finite(nd)) reasons <- c(reasons, "Divergence count unavailable")
  else if (nd > max_divergences) reasons <- c(reasons, "Divergences exceed limit")

  if (min_ess > 0) {
    if (any(!is.finite(c(bulk, tail)))) reasons <- c(reasons, "ESS unavailable")
    else if (min(bulk, tail) < min_ess) reasons <- c(reasons, "ESS below limit")
  }

  paste(reasons, collapse = "; ")
}


# --- Internal helpers ---

.generate_one_dataset <- function(k, di, mean_se, mean_sp,
                                   sd_se_logit, sd_sp_logit, rho_sesp,
                                   omega_se, omega_sp,
                                   mu_prevalence, tau_logit, prop_gold,
                                   moderator, meanlog_n, sdlog_n, add_const_n,
                                   floor_prob, ceiling_prob) {

  # Construct one meta-analysis with known study prevalences and instrument
  # accuracies, then simulate the positive counts that an analyst would
  # observe.

  # Gold-standard studies
  n_gold   <- max(0, round(k * prop_gold))
  n_screen <- k - n_gold

  dat_parts <- list()

  n_all     <- .sample_sizes_logn(k, meanlog_n, sdlog_n, add_const_n)
  theta_all <- .simulate_true_prevalence(k, mu_prevalence, tau_logit)

  # Add the moderator contribution to each study logit prevalence, retaining
  # the true slope for later recovery checks.
  if (!is.null(moderator)) {

    x_mod_all <- stats::rnorm(k, mean = 0, sd = moderator$sd)
    theta_all <- stats::plogis(
      stats::qlogis(theta_all) + moderator$slope * x_mod_all
    )

    slope_true_all <- rep(moderator$slope, k)
  }

  # Generate interview observations using the configured near-perfect ceiling
  # probability for both Se and Sp.
  # --- Gold standard ---
  gold_idx <- seq_len(n_gold)
  if (n_gold > 0) {
    n_i   <- n_all[gold_idx]
    theta <- theta_all[gold_idx]
    Se_g  <- rep(ceiling_prob, n_gold)
    Sp_g  <- rep(ceiling_prob, n_gold)
    fm    <- .forward_misclassification(n_i, theta, Se_g, Sp_g)

    dat_gold <- data.frame(
      study      = seq_len(n_gold),
      n          = n_i,
      y          = fm$y,
      p_obs      = fm$p_obs,
      theta_true = theta,
      Se         = Se_g,
      Sp         = Sp_g,
      measure_id = "gold",
      se_base    = ceiling_prob,
      sp_base    = ceiling_prob,
      is_gold    = TRUE,
      stringsAsFactors = FALSE
    )

    if (!is.null(moderator)) {
      dat_gold$x_mod      <- x_mod_all[gold_idx]
      dat_gold$slope_true <- slope_true_all[gold_idx]
    }

    dat_parts <- c(dat_parts, list(dat_gold))
  }

  # Assign a randomly selected screening instrument to each remaining study,
  # then add study-specific accuracy variation.
  # --- Screening studies ---
  screen_idx <- seq.int(n_gold + 1, k)
  if (n_screen > 0) {
    measures    <- .draw_measures(di, mean_se, mean_sp, sd_se_logit, sd_sp_logit,
                                  rho_sesp, floor_prob, ceiling_prob)
    measure_idx <- sample(measures$measure_id, size = n_screen, replace = TRUE)
    acc         <- .jitter_study_accuracy(measures, measure_idx,
                                           omega_se, omega_sp,
                                           floor_prob, ceiling_prob)

    n_i   <- n_all[screen_idx]
    theta <- theta_all[screen_idx]
    fm    <- .forward_misclassification(n_i, theta, acc$Se, acc$Sp)

    # Keep latent truth and baseline accuracy alongside the observed counts
    # so simulation performance can be evaluated.
    dat_screen <- data.frame(
      study      = seq_len(n_screen) + n_gold,
      n          = n_i,
      y          = fm$y,
      p_obs      = fm$p_obs,
      theta_true = theta,
      Se         = acc$Se,
      Sp         = acc$Sp,
      measure_id = as.character(measure_idx),
      se_base    = measures$se_base[measure_idx],
      sp_base    = measures$sp_base[measure_idx],
      is_gold    = FALSE,
      stringsAsFactors = FALSE
    )

    # Moderator
    if (!is.null(moderator)) {
      dat_screen$x_mod      <- x_mod_all[screen_idx]
      dat_screen$slope_true <- slope_true_all[screen_idx]
    }

    dat_parts <- c(dat_parts, list(dat_screen))
  }

  # Stack interview and screening studies and supply the factor ID expected
  # by the package fitting defaults.
  dat <- do.call(rbind, dat_parts)
  dat$es_id <- factor(dat$study)
  dat
}


.draw_measures <- function(di, mean_se, mean_sp, sd_se_logit, sd_sp_logit,
                            rho, floor_prob, ceiling_prob) {

  # Draw correlated instrument-level sensitivity and specificity on the logit
  # scale. Sigma combines their individual SDs and their correlation.

  mu <- c(stats::qlogis(mean_se), stats::qlogis(mean_sp))
  sds <- c(sd_se_logit, sd_sp_logit)
  R <- matrix(c(1, rho, rho, 1), 2, 2)
  Sigma <- diag(sds) %*% R %*% diag(sds)

  # Draw one Se/Sp pair per instrument, transform to probabilities, and
  # enforce the requested accuracy limits.
  # mvrnorm returns a vector for one instrument; preserve two columns.
  etas <- matrix(MASS::mvrnorm(n = di, mu = mu, Sigma = Sigma), ncol = 2)
  se_base <- stats::plogis(etas[, 1])
  sp_base <- stats::plogis(etas[, 2])
  se_base <- pmin(pmax(se_base, floor_prob), ceiling_prob)
  sp_base <- pmin(pmax(sp_base, floor_prob), ceiling_prob)

  # Store both probability and logit representations, since later steps use
  # them for different parts of generation.
  data.frame(
    measure_id  = seq_len(di),
    eta_se_base = stats::qlogis(se_base),
    eta_sp_base = stats::qlogis(sp_base),
    se_base     = se_base,
    sp_base     = sp_base
  )
}

.jitter_study_accuracy <- function(measures, measure_idx,
                                    omega_se, omega_sp,
                                    floor_prob, ceiling_prob) {

  # Allow studies using the same instrument to have different accuracies by
  # adding independent logit-scale variation around its baseline Se and Sp.

  eta_se <- stats::rnorm(length(measure_idx),
                          mean = measures$eta_se_base[measure_idx],
                          sd = omega_se)
  eta_sp <- stats::rnorm(length(measure_idx),
                          mean = measures$eta_sp_base[measure_idx],
                          sd = omega_sp)

  # Transform the perturbed logits back to probabilities and keep them within
  # the chosen simulation limits.
  Se <- pmin(pmax(stats::plogis(eta_se), floor_prob), ceiling_prob)
  Sp <- pmin(pmax(stats::plogis(eta_sp), floor_prob), ceiling_prob)
  list(Se = Se, Sp = Sp)
}

.sample_sizes_logn <- function(n, meanlog, sdlog, add_const) {

  # Create right-skewed sample sizes, add a minimum-size offset, and round to
  # whole participants.

  round(stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog) + add_const)
}

.simulate_true_prevalence <- function(n, mu, tau) {

  # Draw study log odds around qlogis(mu), with between-study SD tau, then
  # transform to probabilities. mu is the distribution median on the
  # probability scale.

  stats::plogis(stats::rnorm(n, mean = stats::qlogis(mu), sd = tau))
}

.forward_misclassification <- function(n_i, theta_i, Se_i, Sp_i) {

  # Add true positives and false positives to obtain apparent prevalence,
  # then draw an observed count with ordinary binomial sampling variation.
  p_obs <- Se_i * theta_i + (1 - Sp_i) * (1 - theta_i)

  # Use the apparent probability to generate sampling noise; p_obs itself is
  # the expectation, not the realised proportion y/n.
  y <- stats::rbinom(length(n_i), size = n_i, prob = p_obs)
  list(p_obs = p_obs, y = y)
}
