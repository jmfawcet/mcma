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
* New tests for the KWGA scoring and resampling options; case-study vignette
  updated to the corrected estimator.

# mcma 0.1.0

* Initial release.
