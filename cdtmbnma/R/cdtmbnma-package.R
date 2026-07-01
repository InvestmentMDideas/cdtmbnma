#' cdtmbnma: component dose-response network meta-analysis with dose-dependent interaction
#'
#' The package fits an arm-based component network meta-analysis in which each
#' treatment arm is decomposed into component-specific dose-response curves, and
#' the interaction between two components is modelled as a function of the doses
#' of both. That dose-dependent interaction surface collapses the combinatorial
#' set of configuration-level interaction terms into a small number of structural
#' parameters that extrapolate across dose pairs.
#'
#' The workflow is \code{\link{cdt_data}} to build the design, then
#' \code{\link{cdt_fit}} to sample, or \code{\link{cdtmbnma}} to do both in one
#' call. Use \code{summary}, \code{coef}, \code{plot}, and
#' \code{\link{predict.cdtmbnma}} on the result. Estimation uses Stan through
#' \pkg{cmdstanr} or \pkg{rstan}.
#'
#' @section Interaction surfaces:
#' \describe{
#'   \item{none}{additive component effects.}
#'   \item{bilinear}{one parameter per pair; the interaction grows with the
#'     product of the normalised doses.}
#'   \item{gpdi}{a saturating surface adapted from the general pharmacodynamic
#'     interaction model, bilinear at low dose and plateauing as both doses grow.}
#' }
#'
#' @keywords internal
"_PACKAGE"

#' Example dose plane: sacubitril and valsartan blood pressure
#'
#' Reads the eight-week sitting-systolic-blood-pressure dose plane used in the
#' package vignette, eighteen arms across five trials, with a valsartan dose
#' axis, a fixed-ratio diagonal, and off-diagonal arms that vary sacubitril at
#' fixed valsartan.
#'
#' @return A data frame with columns \code{study}, \code{arm}, \code{d_sac},
#'   \code{d_val}, \code{n}, \code{y} (mean change in mmHg), \code{sd},
#'   \code{se}.
#' @export
sacval_example <- function() {
  f <- system.file("extdata", "sacval_msSBP.csv", package = "cdtmbnma")
  if (f == "") stop("Example data not found in the installed package.")
  utils::read.csv(f, stringsAsFactors = FALSE)
}
