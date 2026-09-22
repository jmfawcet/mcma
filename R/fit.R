# Model fitting functions
# =======================

#' Build the brms nonlinear formula
#'
#' Returns the `brms::bf()` object that `mcma_fit()` would construct, without
#' fitting the model. Useful for inspection, debugging, or manual modification.
#'
#' @param bounded Logical; if TRUE, uses Se = 0.5 + 0.5 * inv_logit(eta).
#' @param prev_re One-sided formula for prevalence random effects, e.g.,
#'   `~ (1 | es_id)`. Set to NULL to omit.
#' @param sesp_re One-sided formula for Se/Sp random effects. Set to NULL
#'   to omit (M5-equivalent).
#' @param correlated_re Logical; if TRUE, uses `(1 | s | ...)` to correlate
#'   pi, Se, and Sp random effects.
#' @param moderators One-sided formula for moderator fixed effects, e.g.,
#'   `~ setting + age_centered`. Set to NULL for no moderators.
#' @param measure_col Column name in data identifying the measure.
#' @param gold_column Column name flagging gold-standard studies; used only when
#'   `comparison = TRUE`. Defaults to `"is_gold"`.
#' @param comparison Logical; if TRUE, builds comparison formula with
#'   `pi ~ 0 + is_gold` for separate gold/screen intercepts.
#' @return A brms `bf()` formula object.
#' @export
mcma_formula <- function(bounded       = FALSE,
                         prev_re       = ~ (1 | es_id),
                         sesp_re       = ~ (1 | es_id),
                         correlated_re = TRUE,
                         moderators    = NULL,
                         measure_col   = "measure_id",
                         gold_column   = "is_gold",
                         comparison    = FALSE) {

  # Build the observation model and its three linked regressions. pi, Se, and
  # Sp are linear predictors here; inverse-logit transformations turn them
  # into probabilities.

  # The binomial probability is true positives plus false positives. Bounding
  # Se/Sp changes their transformation, not the count likelihood.
  
  #
  # Main likelihood
  # 
  if (bounded) {
    main <- paste0(
      "y | trials(n) ~ inv_logit(pi) * (0.5 + 0.5 * inv_logit(Se)) + ",
      "(1 - inv_logit(pi)) * (1 - (0.5 + 0.5 * inv_logit(Sp)))"
    )
  } else {
    main <- paste0(
      "y | trials(n) ~ inv_logit(pi) * inv_logit(Se) + ",
      "(1 - inv_logit(pi)) * (1 - inv_logit(Sp))"
    )
  }

  # An ordinary model has one prevalence intercept; a comparison model has
  # separate intercepts for the gold and screening groups.
  
  #
  # Pi submodel
  # 
  if (comparison) {
    pi_fixed <- paste0("0 + ", gold_column)
  } else {
    pi_fixed <- "1"
  }

  # Add requested moderator terms to the prevalence regression, where their
  # effects act on log odds.
  if (!is.null(moderators)) {
    mod_terms <- labels(stats::terms(moderators))
    pi_fixed <- paste(c(pi_fixed, mod_terms), collapse = " + ")
  }

  # Add study/report deviations to prevalence. Matching correlation
  # identifiers link these deviations across the accuracy submodels.
  pi_re <- .build_re_string(prev_re, correlated_re, has_partner_re = !is.null(sesp_re))
  if (nzchar(pi_re)) {
    pi_rhs <- paste(pi_fixed, pi_re, sep = " + ")
  } else {
    pi_rhs <- pi_fixed
  }
  pi_formula <- stats::as.formula(paste("pi ~", pi_rhs))

  # Give each instrument its own accuracy coefficient. Optional random
  # effects allow the same instrument to vary between studies.
  
  #
  # Se/Sp submodels
  # 
  sesp_fixed <- paste0("0 + ", measure_col)
  sesp_re_str <- .build_re_string(sesp_re, correlated_re, has_partner_re = !is.null(prev_re))
  if (nzchar(sesp_re_str)) {
    sesp_rhs <- paste(sesp_fixed, sesp_re_str, sep = " + ")
  } else {
    sesp_rhs <- sesp_fixed
  }
  se_formula <- stats::as.formula(paste("Se ~", sesp_rhs))
  sp_formula <- stats::as.formula(paste("Sp ~", sesp_rhs))

  # Assemble the count model and submodels into one nonlinear brms formula;
  # this step does not fit or sample the model.
  
  #
  # Build bf()
  # 
  brms::bf(
    stats::as.formula(main),
    pi_formula,
    se_formula,
    sp_formula,
    nl = TRUE
  )
}


#' Fit a misclassification-corrected meta-analysis model
#'
#' The primary fitting function. When `priors` is NULL, fits a naive
#' (uncorrected) binomial model. When priors are supplied, constructs the
#' nonlinear formula and fits the corrected model.
#'
#' @param data Data frame with columns: y, n, study/effect-size IDs, measure
#'   IDs, and any moderators.
#' @param priors An mcma_priors object, or NULL for a naive model. Prior rows
#'   for instruments absent from data are ignored with a warning.
#' @param prev_center Prevalence prior centre on the probability scale.
#' @param prev_re One-sided formula for prevalence random effects.
#' @param moderators One-sided formula for moderator fixed effects.
#' @param sesp_re One-sided formula for Se/Sp random effects.
#' @param correlated_re Logical; if TRUE, correlates RE across submodels.
#' @param bounded Logical; if TRUE, uses bounded parameterization.
#' @param measure_col Column name for the measure identifier.
#' @param iter Total iterations per chain.
#' @param warmup Warmup iterations per chain.
#' @param adapt_delta Target acceptance rate.
#' @param max_treedepth Maximum tree depth.
#' @param chains Number of MCMC chains.
#' @param cores Number of cores.
#' @param backend brms backend.
#' @param file File path for model caching. Cached RDS files retain the package
#'   configuration and priors needed by the extraction functions.
#' @param refresh Iteration interval for progress printing (0 = silent).
#' @param step_size Optional initial sampler step size; NULL lets Stan choose it.
#' @param ... Additional arguments passed to `brms::brm()`.
#' @param sd_method Method for converting accuracy-prior uncertainty to the
#'   fitted scale when `bounded = TRUE`. `"chain_rule"` propagates the
#'   probability-scale variance through the full bounded transformation.
#'   `"simple"` is retained for legacy comparisons; see
#'   [prior_sd_from_kappa()]. Ignored for a naive fit.
#' @return A brmsfit object with the mcma_priors attached as an attribute.
#' @export
mcma_fit <- function(data,
                     priors        = NULL,
                     prev_center,
                     prev_re       = ~ (1 | es_id),
                     moderators    = NULL,
                     sesp_re       = ~ (1 | es_id),
                     correlated_re = TRUE,
                     bounded       = FALSE,
                     sd_method = "chain_rule",
                     measure_col   = "measure_id",
                     iter          = 6000,
                     warmup        = floor(iter / 2),
                     adapt_delta   = 0.9999,
                     max_treedepth = 20,
                     chains        = 4,
                     cores         = 4,
                     backend       = "cmdstanr",
                     file          = NULL,
                     refresh       = 50,
                     step_size     = NULL, 
                     ...) {

  # Convert the prevalence prior centre to log odds, then choose the
  # uncorrected or misclassification-adjusted fitting route.

  prev_logit <- stats::qlogis(prev_center)

  if (is.null(priors)) {
    
    #
    # Naive model
    # 
    re_str <- .build_re_string(prev_re, correlated_re = FALSE, has_partner_re = FALSE)
    if (!is.null(moderators)) {
      mod_terms <- labels(stats::terms(moderators))
      fixed <- paste(c("1", mod_terms), collapse = " + ")
    } else {
      fixed <- "1"
    }
    if (nzchar(re_str)) {
      rhs <- paste(fixed, re_str, sep = " + ")
    } else {
      rhs <- fixed
    }
    formula <- stats::as.formula(paste("y | trials(n) ~", rhs))

    # Specify the uncorrected prevalence and heterogeneity priors. This route
    # uses the usual binomial logit link.
    naive_prior <- c(
      brms::set_prior(sprintf("normal(%0.6f, 1.5)", prev_logit),
                       class = "Intercept"),
      # A model without random effects has no SD parameter to assign a prior.
      if (!is.null(prev_re)) brms::set_prior("normal(0, 1)", class = "sd")
    )

    # Pass the model to brms/Stan. Warmup tunes the sampler; later iterations
    # supply the posterior draws.
    fit <- brms::brm(
      formula = formula,
      data    = data,
      family  = stats::binomial(),
      prior   = naive_prior,
      iter    = iter,
      warmup  = warmup,
      chains  = chains,
      cores   = cores,
      backend = backend,
      control = .mcma_sampler_control(adapt_delta, max_treedepth, step_size, backend), 
      file    = file,
      refresh = refresh,
      ...
    )

    # Attach settings used by the extraction and plotting helpers, then
    # return the uncorrected fit.
    config <- list(
      bounded = FALSE, naive = TRUE, measure_col = measure_col,
      prev_center = prev_center
    )
    return(.mcma_finalize_fit(fit, config, file = file, ...))
  }

  # Match priors to this dataset; shared tables may include unused instruments.
  priors <- .mcma_match_priors(priors, data[[measure_col]])

  # --- Corrected model ---
  # Ensure measure_col is character

  # Treat instrument identifiers as labels, even when they look numeric,
  # rather than fitting an accuracy trend across their numeric values.
  data[[measure_col]] <- as.character(data[[measure_col]])

  # Build formula
  formula <- mcma_formula(
    bounded       = bounded,
    prev_re       = prev_re,
    sesp_re       = sesp_re,
    correlated_re = correlated_re,
    moderators    = moderators,
    measure_col   = measure_col
  )

  # Build priors
  # Use brms coefficient names after expanding factor contrasts/interactions.
  mod_terms <- NULL
  if (!is.null(moderators)) {
    prior_table <- brms::get_prior(
      formula = formula, data = data,
      family = stats::binomial(link = "identity")
    )
    mod_terms <- unique(prior_table$coef[
      prior_table$class == "b" & prior_table$nlpar == "pi" &
        nzchar(prior_table$coef) & prior_table$coef != "Intercept"
    ])
  }

  # Build accuracy, prevalence, heterogeneity, and moderator priors on their
  # fitted scales.
  brms_prior <- as_brms_prior(
    priors          = priors,
    prev_center     = prev_center,
    bounded         = bounded,
    sd_method = sd_method,
    has_prev_re     = !is.null(prev_re),
    has_sesp_re     = !is.null(sesp_re),
    moderator_terms = mod_terms,
    measure_col = measure_col
  )

  # The nonlinear formula already calculates a probability, so use the
  # identity link to avoid applying a second logistic transformation.
  fit <- brms::brm(
    formula = formula,
    data    = data,
    family  = stats::binomial(link = "identity"),
    prior   = brms_prior,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    cores   = cores,
    backend = backend,
    control = .mcma_sampler_control(adapt_delta, max_treedepth, step_size, backend), 
    file    = file,
    refresh = refresh,
    ...
  )

  # Retain the supplied priors and parameterization alongside the fit so
  # downstream summaries can interpret its coefficients.
  config <- list(
    bounded = bounded, naive = FALSE, measure_col = measure_col,
    prev_center = prev_center, correlated_re = correlated_re
  )
  .mcma_finalize_fit(fit, config, priors = priors, file = file, ...)
}


#' Fit a joint model (M8) estimating Se/Sp from data
#'
#' Fits the joint model that simultaneously estimates prevalence, Se, and Sp.
#' Uses an accuracy-group parameterization (gold vs screen) rather than
#' per-measure Se/Sp coefficients.
#'
#' @details Concealment (`c`) and over-diagnosis / additional false positives
#'   (`o`) enter the observation model as multipliers on sensitivity and
#'   specificity: gold rows use `1 - c` and `1 - o`, screening rows use
#'   `1 - c * gamma` and `1 - o * delta`. Each term can be handled in one of
#'   three ways: absent (`c = 0` / `o = 0` with no prior, the default, which
#'   reproduces the original model), fixed at a known value (`c` / `o`, scalar
#'   or per row), or treated as a single model parameter in \[0, 1\] with a prior
#'   (`c_prior` / `o_prior`, e.g. `"beta(4, 12)"`). These parameters are
#'   estimated jointly with prevalence and accuracy. The likelihood can update
#'   their distributions even when they are not separately identifiable, so
#'   their posteriors need not reproduce their priors. This route propagates
#'   their joint posterior uncertainty into the prevalence estimate within a
#'   single fit. `gamma` and `delta` are always fixed. Baseline gold accuracy
#'   retains the existing near-perfect specification. Screening priors and
#'   extracted Se/Sp refer to baseline accuracy before these additional
#'   effects; use [extract_bias()] to summarise the terms themselves. Bounding
#'   applies to baseline accuracy, not the effective accuracy after adjustment.
#'   When reusing a cache path across scenarios, pass
#'   `file_refit = "on_change"` through `...`.
#' @param c,o Fixed interview concealment and additional false-positive
#'   probabilities, respectively. Each is a scalar or a vector of length
#'   `nrow(data)`, with values in \[0, 1\]. Zero (the default) removes the term.
#' @param c_prior,o_prior Optional prior for the concealment and
#'   over-diagnosis parameters, as a Stan distribution string on the
#'   probability scale (e.g. `"beta(4, 12)"`); the parameter is bounded to
#'   \[0, 1\]. `NULL` (the default) keeps the term fixed at `c` / `o`. A prior
#'   and a non-zero fixed value cannot both be given for the same term.
#' @param gamma,delta Fixed screening weights in \[0, 1\], scalar or length
#'   `nrow(data)`. Gold rows always use weight one. `o * delta` is the
#'   probability of making an otherwise correctly classified non-case screen
#'   positive. See [rg_correct_bias()].
#'
#' @param data Data frame with y, n, study/es_id, and a gold indicator column.
#' @param gold_column Name of the logical column identifying gold-standard
#'   studies.
#' @param prev_center Prevalence prior centre (probability scale).
#' @param se_prior_center Prior centre for screening Se.
#' @param sp_prior_center Prior centre for screening Sp.
#' @param se_prior_sd Prior SD for screening Se on the (inner) logit scale.
#' @param sp_prior_sd Prior SD for screening Sp on the (inner) logit scale.
#' @param prev_re One-sided formula for prevalence random effects.
#' @param sesp_re One-sided formula for Se/Sp random effects.
#' @param correlated_re Logical; if TRUE, REs are correlated across submodels.
#' @param bounded Logical; whether to use the bounded parameterization.
#' @param study_col Column used by the default prevalence and accuracy
#'   random-effects formulas. Explicitly supplied formulas are used unchanged.
#' @param file File path for caching. Cached RDS files retain the package
#'   configuration needed by the extraction functions.
#' @param refresh Iteration interval for progress printing (0 = silent).
#' @param step_size Optional initial sampler step size; NULL lets Stan choose it.
#' @param ... Additional arguments passed to `brms::brm()`.
#' @inheritParams mcma_fit
#' @return A brmsfit object.
#' @export
mcma_fit_joint <- function(data,
                           gold_column    = "is_gold",
                           prev_center,
                           se_prior_center = 0.80,
                           sp_prior_center = 0.80,
                           se_prior_sd     = 1.0,
                           sp_prior_sd     = 1.0,
                           prev_re         = ~ (1 | es_id),
                           sesp_re         = ~ (1 | es_id),
                           correlated_re   = TRUE,
                           bounded         = TRUE,
                           study_col       = "es_id",
                           file            = NULL,
                           iter            = 6000,
                           warmup          = floor(iter / 2),
                           adapt_delta     = 0.9999,
                           max_treedepth   = 20,
                           chains          = 4,
                           cores           = 4,
                           backend         = "cmdstanr",
                           refresh         = 50,
                           c = 0, o = 0, gamma = 1, delta = 1,
                           c_prior = NULL, o_prior = NULL,      # priors on the concealment / over-diagnosis parameters (NULL = fixed at c / o)
                           step_size = NULL, # optional initial sampler step size
                           ...) {

  # Use interview studies to anchor prevalence while estimating a shared
  # screening sensitivity and specificity. Optional bias terms alter the
  # observation probabilities.

  # Use the selected study column only when the caller has not supplied a
  # random-effects formula of their own.
  if (missing(prev_re)) prev_re <- .mcma_study_re(study_col)
  if (missing(sesp_re)) sesp_re <- .mcma_study_re(study_col)

  # Fixed bias multipliers; full strength for gold, weighted for screens.
  .mcma_validate_bias_prior(c_prior, "c_prior")                                             
  .mcma_validate_bias_prior(o_prior, "o_prior")                                             
  if (!is.null(c_prior) && any(c != 0)) stop("Give either a fixed c or c_prior, not both.", call. = FALSE)   
  if (!is.null(o_prior) && any(o != 0)) stop("Give either a fixed o or o_prior, not both.", call. = FALSE)  

  # Validate the fixed probabilities and construct interview/screen weights
  # shared by the fixed and estimated-bias routes.
  mult <- .mcma_bias_multipliers(nrow(data), c, o, gamma, delta,
                                  is_gold = data[[gold_column]])
  wts  <- .mcma_bias_weights(nrow(data), gamma, delta, is_gold = data[[gold_column]])    

  # Choose each bias term independently: absent, fixed by the analyst, or
  # estimated with a supplied prior.
  c_mode <- if (!is.null(c_prior)) "prior" else if (any(mult$se != 1)) "fixed" else "none" 
  o_mode <- if (!is.null(o_prior)) "prior" else if (any(mult$sp != 1)) "fixed" else "none" 
  has_bias <- c_mode != "none" || o_mode != "none"    

  # Include only the row-level multipliers or weights needed by the chosen
  # likelihood.
  if (c_mode == "fixed") data$mcma_cmult <- mult$se    # columns are added only for the terms in use
  if (o_mode == "fixed") data$mcma_omult <- mult$sp    
  if (c_mode == "prior") data$mcma_gwt   <- wts$gwt    # per-row concealment weight (gold 1, screen gamma)
  if (o_mode == "prior") data$mcma_dwt   <- wts$dwt    # per-row over-diagnosis weight (gold 1, screen delta)

  # Use two accuracy groups instead of instrument-specific coefficients.
  # Their fixed level order makes prior names predictable.
  # Create acc_group factor
  data$acc_group <- factor(
    ifelse(data[[gold_column]], "gold", "screen"),
    levels = c("gold", "screen")
  )

  prev_logit <- stats::qlogis(prev_center)

  # The original likelihood is unchanged when both bias rates are zero.

  # First construct baseline accuracy, then multiply it by the fraction
  # retained after concealment or additional false positives.
  se_expr <- if (bounded) "(0.5 + 0.5 * inv_logit(Se))" else "inv_logit(Se)"
  sp_expr <- if (bounded) "(0.5 + 0.5 * inv_logit(Sp))" else "inv_logit(Sp)"
  se_expr <- switch(c_mode,                                          
    none  = se_expr,
    fixed = paste("mcma_cmult *", se_expr),
    prior = paste("(1 - conc * mcma_gwt) *", se_expr))
  sp_expr <- switch(o_mode,                                          
    none  = sp_expr,
    fixed = paste("mcma_omult *", sp_expr),
    prior = paste("(1 - overdx * mcma_dwt) *", sp_expr))

  # Combine adjusted sensitivity, adjusted specificity, and true prevalence
  # into the probability of an observed positive count.
  main <- paste0(
    "y | trials(n) ~ inv_logit(pi) * ", se_expr, " + ",
    "(1 - inv_logit(pi)) * (1 - ", sp_expr, ")"
  )

  # Build the requested study/report structure for prevalence and accuracy;
  # matching group IDs allow their deviations to correlate.
  pi_re <- .build_re_string(prev_re, correlated_re, has_partner_re = !is.null(sesp_re))
  sesp_re_str <- .build_re_string(sesp_re, correlated_re, has_partner_re = !is.null(prev_re))

  pi_rhs <- if (nzchar(pi_re)) paste("1 +", pi_re) else "1"
  sesp_rhs <- if (nzchar(sesp_re_str)) {
    paste("0 + acc_group +", sesp_re_str)
  } else {
    "0 + acc_group"
  }

  # A prior-driven bias term adds a single shared coefficient on the
  # probability scale. Fixed terms need only data columns.
  nl_forms <- list(                                                 
    stats::as.formula(main),
    stats::as.formula(paste("pi ~", pi_rhs)),
    stats::as.formula(paste("Se ~", sesp_rhs)),
    stats::as.formula(paste("Sp ~", sesp_rhs))
  )
  if (c_mode == "prior") nl_forms <- c(nl_forms, list(stats::as.formula("conc ~ 1")))     
  if (o_mode == "prior") nl_forms <- c(nl_forms, list(stats::as.formula("overdx ~ 1")))   
  formula <- do.call(brms::bf, c(nl_forms, list(nl = TRUE)))        

  # --- Priors ---
  # Se/Sp inner logit centres for screening
  if (bounded) {
    se_inner <- stats::qlogis(map_to_inner(se_prior_center))
    sp_inner <- stats::qlogis(map_to_inner(sp_prior_center))
  } else {
    se_inner <- stats::qlogis(se_prior_center)
    sp_inner <- stats::qlogis(sp_prior_center)
  }

  # Set population and heterogeneity priors. Gold coefficients are fixed near
  # perfect; screening coefficients remain estimated from data and priors.
  pri <- c(
    # Prevalence
    brms::set_prior(sprintf("normal(%0.6f, 1.5)", prev_logit), nlpar = "pi"),
    # Include SD priors only for random effects actually requested.
    if (!is.null(prev_re))
      brms::set_prior("normal(0, 1)", nlpar = "pi", class = "sd"),
    if (!is.null(sesp_re))
      brms::set_prior("normal(0, 0.5)", nlpar = "Se", class = "sd"),
    if (!is.null(sesp_re))
      brms::set_prior("normal(0, 0.5)", nlpar = "Sp", class = "sd"),
    # Gold: fixed at near-perfect
    brms::prior_string("constant(10)", class = "b",
                        coef = "acc_groupgold", nlpar = "Se"),
    brms::prior_string("constant(10)", class = "b",
                        coef = "acc_groupgold", nlpar = "Sp"),
    # Screening: wide priors
    brms::prior_string(sprintf("normal(%0.6f, %0.6f)", se_inner, se_prior_sd),
                        class = "b", coef = "acc_groupscreen", nlpar = "Se"),
    brms::prior_string(sprintf("normal(%0.6f, %0.6f)", sp_inner, sp_prior_sd),
                        class = "b", coef = "acc_groupscreen", nlpar = "Sp")
  )
  # concealment / over-diagnosis parameters live on the probability scale in [0, 1]
  if (c_mode == "prior") pri <- c(pri, brms::set_prior(c_prior, class = "b", nlpar = "conc",   lb = 0, ub = 1))   
  if (o_mode == "prior") pri <- c(pri, brms::set_prior(o_prior, class = "b", nlpar = "overdx", lb = 0, ub = 1))   

  # Fit all active parameters jointly; brms applies the binomial likelihood
  # to the observation probability constructed above.
  fit <- brms::brm(
    formula = formula,
    data    = data,
    family  = stats::binomial(link = "identity"),
    prior   = pri,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    cores   = cores,
    backend = backend,
    control = .mcma_sampler_control(adapt_delta, max_treedepth, step_size, backend), 
    file    = file,
    refresh = refresh,
    ...
  )

  # Store the bias specifications as well as the accuracy transformation for
  # later extraction and interpretation.
  config <- list(
    bounded = bounded, joint = TRUE, gold_column = gold_column,
    study_col = study_col, prev_center = prev_center,
    bias = list(c = c, o = o, gamma = gamma, delta = delta, 
                c_prior = c_prior, o_prior = o_prior,       
                c_mode = c_mode, o_mode = o_mode)           
  )
  .mcma_finalize_fit(fit, config, file = file, ...)
}


#' Fit a comparison model with separate gold/screen intercepts
#'
#' Fits an uncorrected binomial model with separate intercepts for
#' gold-standard and screening studies. Used as input to `mcma_kwga()`.
#'
#' @param data Data frame with y, n, study/es_id, and a gold indicator column.
#' @param gold_column Name of the logical column identifying gold-standard
#'   studies.
#' @param gold_prev_center Prevalence prior centre (probability scale).
#' @param scr_prev_center Prevalence prior centre (probability scale) for screeners.
#' @param gold_prev_sd_logit SD prior (logit scale).
#' @param scr_prev_sd_logit SD prior (logit scale) for screeners.
#' @param prev_re One-sided formula for random effects.
#' @param study_col Column used by the default prevalence random-effects
#'   formula. An explicitly supplied formula is used unchanged.
#' @param file File path for caching. Cached RDS files retain the package
#'   configuration needed by the extraction functions.
#' @param refresh Iteration interval for progress printing (0 = silent).
#' @param step_size Optional initial sampler step size; NULL lets Stan choose it.
#' @param ... Additional arguments passed to `brms::brm()`.
#' @inheritParams mcma_fit
#' @return A brmsfit object.
#' @export
mcma_fit_comparison <- function(data,
                                gold_column   = "is_gold",
                                gold_prev_center,
                                scr_prev_center,
                                gold_prev_sd_logit = 1.5,
                                scr_prev_sd_logit = 1.5,
                                prev_re       = ~ (1 | es_id),
                                study_col     = "es_id",
                                file          = NULL,
                                iter          = 6000,
                                warmup        = floor(iter / 2),
                                adapt_delta   = 0.9999,
                                max_treedepth = 20,
                                chains        = 4,
                                cores         = 4,
                                backend       = "cmdstanr",
                                refresh       = 50,
                                step_size     = NULL, 
                                ...) {

  # Estimate separate apparent prevalence levels for interviews and screens.
  # These group estimates supply the inputs for the grid-based accuracy
  # calibration.

  # An explicit random-effects formula takes precedence over study_col.
  if (missing(prev_re)) prev_re <- .mcma_study_re(study_col)

  gold_prev_center_logit <- stats::qlogis(gold_prev_center)
  scr_prev_center_logit <- stats::qlogis(scr_prev_center)

  # Ensure is_gold is a factor for 0 + is_gold formula
  data[[gold_column]] <- factor(data[[gold_column]])

  # Use one coefficient per group, with no shared intercept, plus the
  # requested between-study/report variation.
  re_str <- .build_re_string(prev_re, correlated_re = FALSE, has_partner_re = FALSE)
  rhs <- paste0("0 + ", gold_column)
  if (nzchar(re_str)) rhs <- paste(rhs, re_str, sep = " + ")
  formula <- stats::as.formula(paste("y | trials(n) ~", rhs))

  # Match each factor level to the intended interview or screening prior,
  # using the centres on the log-odds scale.
  # Priors for both gold and screen intercepts
  levels_gold <- levels(data[[gold_column]])
  pri_list <- lapply(levels_gold, function(lev) {
    is_gold_level <- toupper(as.character(lev)) %in% c("TRUE", "1")
    center <- if (is_gold_level) gold_prev_center_logit else scr_prev_center_logit
    sd_l   <- if (is_gold_level) gold_prev_sd_logit     else scr_prev_sd_logit
    brms::set_prior(
      sprintf("normal(%0.6f, %0.6f)", center, sd_l),
      class = "b", coef = paste0(gold_column, lev)
    )
  })
  pri <- do.call(c, pri_list)
  # Only specify heterogeneity when the comparison includes random effects.
  if (!is.null(prev_re)) {
    pri <- c(pri, brms::set_prior("normal(0, 1)", class = "sd"))
  }

  # Fit apparent positive-count prevalence with the standard binomial logit
  # link; no accuracy correction is applied in this comparison model.
  fit <- brms::brm(
    formula = formula,
    data    = data,
    family  = stats::binomial(),
    prior   = pri,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    cores   = cores,
    backend = backend,
    control = .mcma_sampler_control(adapt_delta, max_treedepth, step_size, backend), 
    file    = file,
    refresh = refresh,
    ...
  )

  # Record the column names and prior centres so the fitted comparison can be
  # identified later.
  config <- list(
    comparison = TRUE, gold_column = gold_column, study_col = study_col,
    prev_center = gold_prev_center,  scr_prev_center = scr_prev_center
  )
  .mcma_finalize_fit(fit, config, file = file, ...)
}


.mcma_sampler_control <- function(adapt_delta, max_treedepth, step_size, backend) {

  # brms forwards control names to the selected backend. CmdStanR calls
  # the option step_size; RStan calls it stepsize. Keep this distinction here
  # so simulation repair can use one argument for every package fitter.
  control <- list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  if (!is.null(step_size)) {
    if (!is.numeric(step_size) || length(step_size) != 1L ||
        !is.finite(step_size) || step_size <= 0) {
      stop("step_size must be a positive finite number or NULL.", call. = FALSE)
    }
    key <- if (backend == "cmdstanr") "step_size" else "stepsize"
    control[[key]] <- step_size
  }

  control
}


# --- Internal helper: build RE string from one-sided formula ---

.build_re_string <- function(re_formula, correlated_re, has_partner_re) {

  # Turn a random-effects formula into text that can be inserted into each
  # submodel. NULL means that this submodel has no random effects.

  if (is.null(re_formula)) return("")

  # Extract the right-hand side of the one-sided formula and join any wrapped
  # lines into one string.
  re_text <- paste(deparse(re_formula[[2]]), collapse = " ")

  # If correlated_re is TRUE and there's a partner submodel with REs, add a
  # correlation ID per grouping factor so that matching groups across
  # submodels share a correlation block while distinct groups do not.
  if (correlated_re && has_partner_re) {

    # Find ordinary grouping terms, then give each grouping factor its own
    # shared correlation label.
    pattern <- "\\(([^|]+)\\|\\s*([^|)]+)\\)"
    m  <- gregexpr(pattern, re_text)[[1]]
    if (m[1] != -1) {
      hits <- regmatches(re_text, list(m))[[1]]
      rebuilt <- vapply(hits, function(h) {
        lhs <- trimws(sub(pattern, "\\1", h))
        grp <- trimws(sub(pattern, "\\2", h))
        id  <- paste0("mcma_", gsub("[^A-Za-z0-9_]", "_", grp))
        sprintf("(%s | %s | %s)", lhs, id, grp)
      }, character(1))

      # Replace the matched terms while leaving the rest of the requested
      # formula intact.
      regmatches(re_text, list(m)) <- list(rebuilt)
    }
  }

  re_text
}


.mcma_study_re <- function(study_col) {

  # Build the usual study intercept using the requested column name.
  stats::reformulate(sprintf("(1 | `%s`)", study_col))
}

.mcma_finalize_fit <- function(fit, config, priors = NULL, file = NULL, ...) {

  # brms saves before returning, so add our metadata and update that cache.
  # Avoid rewriting a large cached fit when its metadata is already current.
  changed <- !identical(attr(fit, "mcma_config"), config) ||
    !identical(attr(fit, "mcma_priors"), priors)
  attr(fit, "mcma_config") <- config
  attr(fit, "mcma_priors") <- priors

  fit_options <- list(...)
  if (changed && !is.null(file) && !is.null(fit$file) && !isTRUE(fit_options$empty)) {
    compress <- if (is.null(fit_options$file_compress)) TRUE else fit_options$file_compress
    saveRDS(fit, file = fit$file, compress = compress)
  }

  fit
}
