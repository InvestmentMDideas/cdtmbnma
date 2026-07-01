#' Prior settings for cdtmbnma
#'
#' Weakly-informative defaults on the model's structural parameters. All effects
#' live on the outcome scale for a continuous outcome and on the log-odds scale
#' for a binary outcome, so scale the component-effect prior to the outcome.
#'
#' @param emax_sd Prior SD for each component's Emax (maximal effect).
#' @param logED50_mean,logED50_sd Prior mean and SD for each component's log
#'   half-maximal dose. Also used for the saturating-surface half-doses.
#' @param int_sd Prior SD for the interaction parameter (bilinear \eqn{\eta} or
#'   saturating asymptote), centred at additivity.
#' @param ref_sd Prior SD for the free per-study reference level.
#' @param omega_sd Scale of the half-normal prior on the random-effect SD.
#'
#' @return A named list of prior hyperparameters.
#' @export
cdt_priors <- function(emax_sd = 10, logED50_mean = log(50), logED50_sd = 1,
                       int_sd = 5, ref_sd = 10, omega_sd = 2) {
  list(prior_emax_sd = emax_sd,
       prior_logED50_mean = logED50_mean,
       prior_logED50_sd = logED50_sd,
       prior_int_sd = int_sd,
       prior_ref_sd = ref_sd,
       prior_omega_sd = omega_sd)
}
