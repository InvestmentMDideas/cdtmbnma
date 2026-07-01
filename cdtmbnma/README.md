# cdtmbnma

Component dose-response network meta-analysis with a dose-dependent interaction
surface.

`cdtmbnma` fits an arm-based component network meta-analysis in which each
treatment arm is decomposed into component-specific dose-response curves, and the
interaction between two components is modelled as a function of the doses of
both. That surface collapses the combinatorial set of configuration-level
interaction terms into a few structural parameters that extrapolate across dose
pairs, including pairs that no trial observed. The package implements the methods
in the accompanying paper.

## What it does

- Arm-based component network meta-analysis for any number of components.
- Continuous (normal, known standard error) and binary (binomial-logit) outcomes.
- Three interaction surfaces: additive (`none`), `bilinear`, and a saturating
  general-pharmacodynamic-interaction surface (`gpdi`).
- Study-level random effects on relative effects, with free per-study reference
  levels that preserve randomisation.
- Prediction of relative effects at arbitrary dose combinations.
- Estimation in Stan through `cmdstanr` (recommended) or `rstan`.

## Status

This is version 0.1.0, a source release. The Stan model generalises a
two-component implementation that was validated by parameter recovery against a
known generating process, and the R interface and tests are written against that
model. The package has not yet been compiled or `R CMD check`ed in the authoring
environment, which lacked an R and Stan toolchain. Build and check it on a
machine with R and a Stan backend before relying on results.

## Install

```r
# 1. A Stan backend (cmdstanr is recommended)
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()        # one-time toolchain setup

# 2. Documentation, checks, and install from source
install.packages(c("devtools", "roxygen2", "posterior"))
devtools::document("cdtmbnma")     # regenerate man/ from the roxygen comments
devtools::install("cdtmbnma", build_vignettes = TRUE)
devtools::check("cdtmbnma")
```

`man/` is intentionally empty in this source tree; `devtools::document()`
generates the help pages from the roxygen comments in `R/`.

## Quick start

```r
library(cdtmbnma)

sv <- sacval_example()             # the paper's sacubitril/valsartan dose plane

fit <- cdtmbnma(
  sv, study = "study", components = c("d_sac", "d_val"),
  outcome = "continuous", y = "y", se = "se",
  dstar = c(d_sac = 200, d_val = 320),
  interaction = "bilinear")

summary(fit)                       # component effects, interaction, diagnostics
plot(fit)                          # component dose-response curves
predict(fit, data.frame(d_sac = 100, d_val = 160))   # an unobserved combination
```

The data frame has one row per arm, with a dose column for each component and
zero meaning the component is absent. For a binary outcome, pass
`outcome = "binary"`, `events`, and `n_binary` instead of `y` and `se`.

See `vignette("cdtmbnma")` for the full sacubitril and valsartan walkthrough.

## Model

For arm \(i\) in study \(s\), with component doses \(d_{ic}\), the linear
predictor is

\[
\mu_i \;=\; m_s \;+\; \sum_c \frac{a_c\, d_{ic}}{e_c + d_{ic}}
\;+\; \sum_{(c,c')} \psi(d_{ic}, d_{ic'}) \;+\; u_i,
\]

with a free per-study reference level \(m_s\), a per-component Emax
dose-response, a pairwise interaction surface \(\psi\), and a random effect
\(u_i\) on active arms. The bilinear surface is
\(\psi = \eta\,(d_c/d_c^\star)(d_{c'}/d_{c'}^\star)\); the saturating surface is
\(\psi = \mathrm{INT}\,\frac{d_c}{\kappa_c + d_c}\,\frac{d_{c'}}{\kappa_{c'} + d_{c'}}\).
The likelihood is normal with known standard error for a continuous outcome and
binomial-logit for a binary outcome.

## Citation

Cite the methods paper. A `CITATION` file will accompany the first tagged
release.
