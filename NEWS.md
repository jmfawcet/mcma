# mcma 0.1.1

* Kernel-weighted grid averaging (KWGA) is now corrected. `mcma_kwga()` scores the
  unclamped Rogan-Gladen correction of each screening draw against the
  gold-standard draws, so an accuracy pair that implies a negative prevalence is
  penalised by its full value rather than being scored as zero;
  `mcma_kwga_prevalence()` resamples posterior draws and accuracy pairs jointly
  (`resample = "joint"`, chosen automatically). Old version remains available for comparison via 
  `mcma_kwga(clamp_scoring = TRUE)` and `mcma_kwga_prevalence(resample = "independent")`.
* `mcma_kwga()` records the scoring version in its result and in `print()`;
  `mcma_kwga_prevalence()` returns the resampling used and accepts an
  optional `seed`.
* `mcma_kwga_prevalence()` drops the unused `mode` argument and reports
  `n_effective_pairs` (previously `n_effective_models`).
* `bounded = TRUE` is now the default in `mcma_fit()`, `mcma_sensitivity()`,
  `mcma_sensitivity_comparison()` and `mcma_formula()`, matching
  `mcma_fit_joint()` and the analyses in the paper.
* `prior_sd_from_kappa()`, `as_brms_prior()` and `plot()` for `mcma_priors`
  objects also default to `bounded = TRUE`, so the priors they report or draw
  are the ones `mcma_fit()` uses.
* `mcma_sensitivity_comparison()` uses the same sampler defaults (`chains`,
  `cores`, `backend`, `step_size`) as the other fitters. Both sensitivity
  functions cache to a temporary directory unless `model_dir` is supplied, and
  expose `file_refit`.
* `mcma_sim_repair()` writes each refit to a temporary file and keeps the
  current result if the refit fails or cannot be saved.
* New tests for the KWGA scoring and resampling options; case-study vignette
  updated to the corrected estimator.

# mcma 0.1.0

* Initial release.
