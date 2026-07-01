## Internal: locate and compile the Stan model, cached per backend ------------
.cdt_cache <- new.env(parent = emptyenv())

.cdt_backend <- function(backend = c("auto", "cmdstanr", "rstan")) {
  backend <- match.arg(backend)
  has_cmd <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rst <- requireNamespace("rstan", quietly = TRUE)
  if (backend == "cmdstanr" && !has_cmd) stop("cmdstanr is not installed.")
  if (backend == "rstan" && !has_rst) stop("rstan is not installed.")
  if (backend == "auto") {
    backend <- if (has_cmd) "cmdstanr" else if (has_rst) "rstan" else
      stop("Install 'cmdstanr' (recommended) or 'rstan' to fit models.")
  }
  backend
}

.cdt_model <- function(backend) {
  key <- paste0("model_", backend)
  if (!is.null(.cdt_cache[[key]])) return(.cdt_cache[[key]])
  stan_file <- system.file("stan", "cdtmbnma.stan", package = "cdtmbnma")
  if (stan_file == "") stop("Stan model file not found in the installed package.")
  mod <- if (backend == "cmdstanr")
    cmdstanr::cmdstan_model(stan_file)
  else
    rstan::stan_model(file = stan_file)
  .cdt_cache[[key]] <- mod
  mod
}

#' Fit a component dose-response network meta-analysis
#'
#' Compiles the Stan model on first use and samples the posterior. The
#' interaction surface is chosen here; the component and study structure comes
#' from the \code{cdt_data} object.
#'
#' @param data A \code{\link{cdt_data}} object.
#' @param interaction Interaction surface: \code{"bilinear"} (default),
#'   \code{"none"} (additive), or \code{"gpdi"} (saturating, general
#'   pharmacodynamic interaction).
#' @param priors A list from \code{\link{cdt_priors}}.
#' @param newdata Optional data frame of component-dose combinations at which to
#'   predict the relative effect against the all-zero reference. Must carry the
#'   component columns.
#' @param backend One of \code{"auto"}, \code{"cmdstanr"}, or \code{"rstan"}.
#' @param chains,iter_warmup,iter_sampling Sampler settings.
#' @param adapt_delta Target acceptance probability.
#' @param seed Random seed.
#' @param refresh Console refresh interval (0 silences progress).
#' @param ... Passed to the backend sampler.
#'
#' @return An object of class \code{cdtmbnma}.
#' @export
cdt_fit <- function(data, interaction = c("bilinear", "none", "gpdi"),
                    priors = cdt_priors(), newdata = NULL,
                    backend = "auto", chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, seed = 1, refresh = 0, ...) {
  stopifnot(inherits(data, "cdt_data"))
  interaction <- match.arg(interaction)
  icode <- c(none = 0L, bilinear = 1L, gpdi = 2L)[[interaction]]
  if (icode > 0 && data$standata$P == 0)
    warning("No component pairs co-occur, so the interaction is not identifiable; ",
            "fitting the additive model.")

  C <- data$standata$C
  if (!is.null(newdata)) {
    Dp <- as.matrix(newdata[, data$components, drop = FALSE]); storage.mode(Dp) <- "double"
    G <- nrow(Dp)
  } else {
    Dp <- matrix(0, 1, C); G <- 0L
  }

  standata <- c(data$standata, priors,
                list(interaction = icode, G = as.integer(G), Dpred = Dp))

  be <- .cdt_backend(backend)
  mod <- .cdt_model(be)

  if (be == "cmdstanr") {
    fit <- mod$sample(data = standata, chains = chains,
                      parallel_chains = chains,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                      adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
  } else {
    fit <- rstan::sampling(mod, data = standata, chains = chains,
                           warmup = iter_warmup, iter = iter_warmup + iter_sampling,
                           control = list(adapt_delta = adapt_delta),
                           seed = seed, refresh = refresh, ...)
  }

  draws <- posterior::as_draws_df(if (be == "cmdstanr") fit$draws() else as.array(fit))
  structure(list(
    draws = draws, fit = fit, data = data, backend = be,
    spec = list(interaction = interaction, priors = priors),
    newdata = newdata
  ), class = "cdtmbnma")
}

#' Component dose-response network meta-analysis (one-call interface)
#'
#' Convenience wrapper that builds the design with \code{\link{cdt_data}} and
#' fits it with \code{\link{cdt_fit}}.
#'
#' @inheritParams cdt_data
#' @param interaction,priors,newdata,backend,chains,iter_warmup,iter_sampling,adapt_delta,seed,refresh,... Passed to \code{\link{cdt_fit}}.
#' @return An object of class \code{cdtmbnma}.
#' @export
cdtmbnma <- function(data, study, components,
                     outcome = c("continuous", "binary"),
                     y = NULL, se = NULL, sd = NULL, n = NULL,
                     events = NULL, n_binary = NULL,
                     ref = NULL, dstar = NULL, interactions = NULL,
                     interaction = c("bilinear", "none", "gpdi"),
                     priors = cdt_priors(), newdata = NULL, backend = "auto",
                     chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                     adapt_delta = 0.95, seed = 1, refresh = 0, ...) {
  d <- cdt_data(data, study = study, components = components,
                outcome = match.arg(outcome), y = y, se = se, sd = sd, n = n,
                events = events, n_binary = n_binary, ref = ref,
                dstar = dstar, interactions = interactions)
  cdt_fit(d, interaction = match.arg(interaction), priors = priors,
          newdata = newdata, backend = backend, chains = chains,
          iter_warmup = iter_warmup, iter_sampling = iter_sampling,
          adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
}
