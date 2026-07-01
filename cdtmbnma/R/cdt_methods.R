#' @export
print.cdtmbnma <- function(x, ...) {
  d <- x$data
  cat("<cdtmbnma fit>\n")
  cat(sprintf("  %d arms, %d studies, %d components\n",
              d$n_arms, d$n_studies, length(d$components)))
  cat(sprintf("  outcome: %s   interaction: %s   backend: %s\n",
              d$outcome, x$spec$interaction, x$backend))
  cat("\nComponent effects (Emax):\n")
  print(coef(x), row.names = FALSE)
  invisible(x)
}

#' Posterior summary of a cdtmbnma fit
#'
#' @param object A fitted \code{cdtmbnma} object.
#' @param ... Unused.
#' @return A data frame summarising the structural parameters, with convergence
#'   diagnostics.
#' @export
summary.cdtmbnma <- function(object, ...) {
  d <- object$data; C <- length(d$components); P <- nrow(d$pairs)
  itype <- object$spec$interaction
  vars <- c(sprintf("emax[%d]", seq_len(C)),
            sprintf("ED50[%d]", seq_len(C)),
            if (itype == "bilinear" && P) sprintf("eta[%d]", seq_len(P)),
            if (itype == "gpdi" && P) sprintf("INT[%d]", seq_len(P)),
            if (itype == "gpdi") sprintf("kappa[%d]", seq_len(C)),
            "omega")
  s <- posterior::summarise_draws(
    object$draws, "mean", "sd",
    ~ stats::quantile(.x, probs = c(0.025, 0.975)),
    "rhat", "ess_bulk")
  s <- s[s$variable %in% vars, , drop = FALSE]

  lab <- s$variable
  for (c in seq_len(C)) {
    lab <- sub(sprintf("^emax\\[%d\\]$", c), paste0("emax: ", d$components[c]), lab)
    lab <- sub(sprintf("^ED50\\[%d\\]$", c), paste0("ED50: ", d$components[c]), lab)
    lab <- sub(sprintf("^kappa\\[%d\\]$", c), paste0("kappa: ", d$components[c]), lab)
  }
  if (P) for (p in seq_len(P)) {
    nm <- paste(d$components[d$pairs[p, ]], collapse = " x ")
    lab <- sub(sprintf("^eta\\[%d\\]$", p), paste0("interaction: ", nm), lab)
    lab <- sub(sprintf("^INT\\[%d\\]$", p), paste0("interaction: ", nm), lab)
  }
  s$variable <- lab
  as.data.frame(s)
}

#' @export
coef.cdtmbnma <- function(object, ...) {
  d <- object$data; C <- length(d$components)
  em <- .draw_mat(object$draws, "emax", C)
  ed <- .draw_mat(object$draws, "ED50", C)
  data.frame(
    component = d$components,
    emax = colMeans(em),
    emax_q2.5 = apply(em, 2, stats::quantile, 0.025),
    emax_q97.5 = apply(em, 2, stats::quantile, 0.975),
    ED50 = colMeans(ed),
    row.names = NULL, check.names = FALSE)
}

#' Plot component dose-response curves
#'
#' Draws each component's marginal dose-response curve (the named component
#' varied, the others held at zero) with a posterior credible band, on the model
#' scale.
#'
#' @param x A fitted \code{cdtmbnma} object.
#' @param ngrid Number of dose points per curve.
#' @param probs Quantiles for the band.
#' @param ... Passed to \code{plot}.
#' @return Invisibly, a list of the plotted grids.
#' @export
plot.cdtmbnma <- function(x, ngrid = 60, probs = c(0.025, 0.975), ...) {
  d <- x$data; comp <- d$components; C <- length(comp)
  Dobs <- d$standata$D
  draws <- x$draws
  emax <- .draw_mat(draws, "emax", C); ED50 <- .draw_mat(draws, "ED50", C)

  op <- graphics::par(mfrow = c(1, C), mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
  grids <- vector("list", C)
  for (c in seq_len(C)) {
    dmax <- max(Dobs[, c]); if (dmax == 0) dmax <- 1
    g <- seq(0, dmax, length.out = ngrid)
    M <- sapply(g, function(dose) emax[, c] * dose / (ED50[, c] + dose))
    m <- colMeans(M); lo <- apply(M, 2, stats::quantile, probs[1]); hi <- apply(M, 2, stats::quantile, probs[2])
    graphics::plot(g, m, type = "n", ylim = range(lo, hi, 0),
                   xlab = paste0(comp[c], " dose"), ylab = "effect", main = comp[c])
    graphics::polygon(c(g, rev(g)), c(lo, rev(hi)), col = grDevices::adjustcolor("steelblue", 0.2), border = NA)
    graphics::lines(g, m, col = "steelblue", lwd = 2)
    graphics::abline(h = 0, col = "grey60", lty = 2)
    grids[[c]] <- data.frame(dose = g, mean = m, lo = lo, hi = hi)
  }
  invisible(stats::setNames(grids, comp))
}
