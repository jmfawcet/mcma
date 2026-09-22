# Misclassification Corrected Meta-analysis (mcma)

<!-- badges: start -->
![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)
![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-1f72aa.svg)
<!-- badges: end -->

Bayesian misclassification-corrected meta-analysis of disorder prevalence in R. Built on 
[**brms**](https://paul-buerkner.github.io/brms/) and Stan (via `cmdstanr`), `mcma` pools 
prevalence estimates across studies while modelling the imperfect sensitivity (*Se*) and 
specificity (*Sp*) of the instruments that produced them, so that screening-based and 
diagnostic interview-based ("gold standard" for our purposes) studies can be synthesized 
together without the upward bias of naive pooling.


**Companion repository (paper, simulations and case study):** https://github.com/jmfawcet/mcma_case

The correction is built on the Rogan–Gladen relationship between the **observed** (apparent) screening 
rate and the **true** prevalence θ:

```
p_obs = Se · θ + (1 − Sp) · (1 − θ)
```

`mcma` embeds this relationship inside hierarchical Bayesian models and provides tools for specifying 
accuracy priors from external evidence, estimating *Se*/*Sp* from the data when gold-standard studies are 
available, running prior-sensitivity analyses, and validating the workflow by simulation.

## Companion paper

The methods implemented here are described in:

> Fawcett, J. M., Whitridge, J., Bartoš, F., & Fawcett, E. J. *Bayesian meta-analysis of disorder prevalence with misclassification correction.* Manuscript submitted for publication.

A preprint link will be added here when available. All simulation code, simulation outputs and the case-study analysis reported in the paper are in the companion repository, [**mcma_case**](https://github.com/jmfawcet/mcma_case).

## Installation

`mcma` is not on CRAN. Install from GitHub with [**remotes**](https://remotes.r-lib.org/):

```r
# install.packages("remotes")
remotes::install_github("jmfawcet/mcma")
```

To build the case-study vignette as well (requires [Quarto](https://quarto.org/) and an internet connection, since the vignette downloads its data from OSF):

```r
remotes::install_github("jmfawcet/mcma", build_vignettes = TRUE)
```

### Stan backend

`mcma` fits models through **brms** using the **cmdstanr** backend by default, so you will also need
[`cmdstanr`](https://mc-stan.org/cmdstanr/) and a working CmdStan installation:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

(Pass `backend = "rstan"` to the fitting functions if you prefer rstan)

## Quick start

A minimal corrected meta-analysis. Your data should be one row per study × measure, with columns `y` 
(positives), `n` (sample size), an effect-size/study id (`es_id`), a `measure_id`, and a logical `is_gold` flag.

```r
library(mcma)

# 1. Specify sensitivity/specificity priors from external evidence.
#    One row per measure; gold-standard measures are flagged is_gold = TRUE
#    (Se/Sp are fixed, so se/sp/kappa are NA).
priors <- mcma_priors(
  measure_id = c("screen_A", "screen_B", "gold"),
  se         = c(0.85, 0.80, NA),
  sp         = c(0.75, 0.85, NA),
  kappa      = c(100,  100,  NA),   # prior "effective sample size"
  is_gold    = c(FALSE, FALSE, TRUE)
)

# 2. Fit the misclassification-corrected hierarchical model.
fit <- mcma_fit(
  data        = my_data,
  priors      = priors,
  prev_center = 0.06,               # prior centre for population prevalence
  prev_re     = ~ (1 | es_id),      # study-level random effects on prevalence
  sesp_re     = ~ (1 | es_id)
)

# 3. Summarise the corrected pooled prevalence.
extract_prevalence(fit)

# Diagnostics
mcma_convergence(fit)
mcma_ppc(fit)
```

Passing `priors = NULL` to `mcma_fit()` fits a naive (uncorrected) binomial model for comparison.

### When you don't have external *Se*/*Sp*

If some studies use a gold standard, you can **estimate** *Se*/*Sp* from the data instead of supplying priors:

```r
# Discrete kernel-weighted grid averaging (KWGA) over a grid of (Se, Sp) candidates
fit_comp <- mcma_fit_comparison(data = my_data, gold_column = "is_gold",
                                gold_prev_center = 0.06, scr_prev_center = 0.20)
kwga      <- mcma_kwga(fit_comp, se_grid = seq(0.6, 0.95, 0.025),
                                sp_grid = seq(0.6, 0.95, 0.025))
mcma_kwga_prevalence(kwga, fit_comparison = fit_comp)

# or, better yet, a joint model that estimates prevalence, Se, and Sp simultaneously
fit_joint <- mcma_fit_joint(data = my_data, gold_column = "is_gold",
                            prev_center = 0.06, bounded = TRUE)
```

### Prior sensitivity

```r
sens <- mcma_sensitivity(data = my_data, priors = priors, prev_center = 0.06,
                         se_grid = seq(0.70, 0.95, 0.05),
                         sp_grid = seq(0.70, 0.95, 0.05))
plot(sens, type = "heatmap")
```

## Function reference

| Area | Functions |
|------|-----------|
| **Model fitting** | `mcma_fit()`, `mcma_fit_joint()`, `mcma_fit_comparison()`, `mcma_formula()` |
| **Accuracy priors** | `mcma_priors()`, `as_brms_prior()`, `prior_sd_from_kappa()`, `update()` / `plot()` methods |
| **Estimating Se/Sp from data** | `mcma_kwga()`, `mcma_kwga_prevalence()` |
| **Posterior summaries** | `extract_prevalence()`, `extract_sesp()`, `extract_tau()`, `extract_slope()`, `extract_gold_screen_diff()` |
| **Sensitivity analysis** | `mcma_sensitivity()`, `mcma_sensitivity_comparison()` |
| **Diagnostics** | `mcma_convergence()`, `mcma_ppc()` |
| **Misclassification algebra** | `rg_correct()`, `rg_forward()`, `implied_se()`, `implied_sp()`, `iso_correction_line()` |
| **Plotting** | `plot_joint_sesp()`, `plot()` for `mcma_priors` / `mcma_kwga` / `mcma_sensitivity` |
| **Simulation** | `mcma_sim_data()`, `mcma_sim_run()`, `mcma_sim_repair()`, `mcma_sim_read()`, `mcma_sim_summarise()` |

A worked end-to-end analysis, reproducing the depression case study from the paper, is provided in the package vignette:

```r
vignette("case-study", package = "mcma")
```

## Citation

If you use `mcma` in your work, please cite both the companion paper and the package:

- Fawcett, J. M., Whitridge, J., Bartoš, F., & Fawcett, E. J. *Bayesian meta-analysis of disorder prevalence with misclassification correction.* Manuscript submitted for publication.

- Fawcett, J. M., & Whitridge, J. (2026). *mcma: Bayesian Misclassification-Corrected Meta-Analysis* (Version 0.1.0) [R package]. https://github.com/jmfawcet/mcma

## License

MIT. See [`LICENSE.md`](LICENSE.md).
