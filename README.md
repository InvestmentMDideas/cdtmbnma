# cdtmbnma

**Component, dose, and time network meta-analysis with dose-dependent interaction surfaces.**

This repository is the research compendium for the methods paper *Component, dose,
and time network meta-analysis with dose-dependent interaction surfaces*. It holds
the `cdtmbnma` R package and the worked data-extraction artifacts for the paper's
empirical applications. The manuscript and its supplement are distributed separately
and are not included here.

`cdtmbnma` fits an arm-based component network meta-analysis in which each treatment
arm is decomposed into component-specific dose-response curves, and the interaction
between two components is modelled as a function of the doses of both. That surface
collapses the combinatorial set of configuration-level interaction terms into a few
structural parameters that extrapolate across dose pairs, including pairs that no
trial observed.

## Contents

| Path | Description |
|------|-------------|
| `cdtmbnma/` | Source of the `cdtmbnma` R package (arm-based component dose-response NMA with a dose-dependent interaction surface, estimated in Stan). |
| `cdtmbnma_0.2.1.tar.gz` | Built source tarball of the package, version 0.2.1. |
| `data/sacval_msSBP_doseplane.csv` | Sacubitril–valsartan sitting-SBP dose plane analysed in the paper (one row per trial arm). |
| `data/copd_bgf_extraction_form.csv` | Arm-level extraction skeleton for the COPD triple-therapy (budesonide / glycopyrronium / formoterol) network. |
| `data/copd_bgf_data_dictionary.csv` | Field definitions and extraction rules for the COPD extraction form. |
| `data/antihtn_factorial_extraction_form.csv` | Arm-level extraction skeleton for the antihypertensive factorial (amlodipine + ARB) network. |
| `data/antihtn_factorial_data_dictionary.csv` | Field definitions and extraction rules for the antihypertensive extraction form. |

## The `cdtmbnma` package

- Arm-based component network meta-analysis for any number of components.
- Continuous (normal, known standard error) and binary (binomial–logit) outcomes.
- Three interaction surfaces: additive (`none`), `bilinear`, and a saturating
  general pharmacodynamic interaction surface (`gpdi`).
- Study-level random effects on relative effects, with free per-study reference
  levels that preserve randomisation.
- Prediction of relative effects at arbitrary dose combinations, including
  combinations that no trial observed.
- Estimation in Stan through `cmdstanr` (recommended) or `rstan`.

### Install

```r
# 1. A Stan backend (cmdstanr recommended)
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

# 2. Document, install, and check from source (run from the repository root)
install.packages(c("devtools", "roxygen2", "posterior"))
devtools::document("cdtmbnma")
devtools::install("cdtmbnma", build_vignettes = TRUE)
devtools::check("cdtmbnma")
```

`man/` is intentionally empty in the source tree; `devtools::document()` generates
the help pages from the roxygen comments in `cdtmbnma/R/`.

### Quick start

```r
library(cdtmbnma)

sv <- sacval_example()             # the paper's sacubitril/valsartan dose plane

fit <- cdtmbnma(
  sv, study = "study", components = c("d_sac", "d_val"),
  outcome = "continuous", y = "y", se = "se",
  dstar = c(d_sac = 200, d_val = 320),
  interaction = "bilinear")

summary(fit)                                       # effects, interaction, diagnostics
plot(fit)                                          # component dose-response curves
predict(fit, data.frame(d_sac = 100, d_val = 160)) # an unobserved combination
```

See `vignette("cdtmbnma")` for the full sacubitril and valsartan walkthrough.

## Data

Each `data/` file is an arm-level table for one empirical application, with a dose
column per component and `0` encoding an absent component:

- **`sacval_msSBP_doseplane.csv`** — the analysed dataset for the sacubitril–valsartan
  application. Columns: `study`, `arm`, `d_sac`, `d_val`, `n`, `y` (per-arm change
  from baseline in sitting systolic blood pressure, mmHg), `sd`, `se`, and `sd_source`
  (how each arm's dispersion was obtained: reported/network-table, imputed, or
  reconstructed from published least-squares-mean standard errors).
- **`copd_bgf_extraction_form.csv`** / **`copd_bgf_data_dictionary.csv`** — a worked
  arm-level extraction instrument for the inhaled triple-therapy programme (ETHOS,
  KRONOS), with the dose and design cells pre-filled and the outcome cells left for
  extraction, plus the accompanying field definitions.
- **`antihtn_factorial_extraction_form.csv`** / **`antihtn_factorial_data_dictionary.csv`**
  — a worked arm-level extraction template for a factorial antihypertensive design
  (amlodipine crossed with an angiotensin-receptor blocker), the off-diagonal dose
  design that identifies a dose-dependent interaction.

## Status

Version 0.2.1 is a source release. The general Stan model generalises a two-component
implementation that was validated by parameter recovery against a known generating
process, and the R interface and tests are written against that model. Build and check
the package on a machine with R and a Stan backend before relying on results.

## License

Released under the MIT License. Copyright © 2026 Tyler Pitre. See [LICENSE](LICENSE).

## Citation

Please cite the methods paper. A `CITATION` file will accompany the first tagged release.
