## Internal: pull a named parameter block as a draws-by-index matrix ----------
.draw_mat <- function(draws, base, k) {
  if (k == 0) return(matrix(0, nrow(draws), 0))
  cols <- sprintf("%s[%d]", base, seq_len(k))
  as.matrix(draws[, cols, drop = FALSE])
}

#' Predict relative effects at dose combinations
#'
#' Computes the posterior of the relative effect against the all-zero reference
#' (placebo) at one or more component-dose combinations, including combinations
#' and dose pairs that no trial observed. The random effect is marginalised at
#' zero, so the prediction is the population-average structural effect rather
#' than a new-study predictive interval. The effect is on the outcome scale for
#' a continuous outcome and on the log-odds scale for a binary outcome.
#'
#' @param object A fitted \code{cdtmbnma} object.
#' @param newdata Data frame of dose combinations, carrying the component
#'   columns used in fitting.
#' @param probs Quantiles for the credible interval.
#' @param ... Unused.
#'
#' @return A data frame with the posterior mean and quantiles for each row of
#'   \code{newdata}.
#' @export
predict.cdtmbnma <- function(object, newdata, probs = c(0.025, 0.975), ...) {
  d <- object$data
  comp <- d$components; C <- length(comp)
  miss <- setdiff(comp, names(newdata))
  if (length(miss)) stop("newdata is missing component column(s): ", paste(miss, collapse = ", "))
  Dp <- as.matrix(newdata[, comp, drop = FALSE]); storage.mode(Dp) <- "double"

  draws <- object$draws
  emax <- .draw_mat(draws, "emax", C)
  ED50 <- .draw_mat(draws, "ED50", C)
  itype <- object$spec$interaction
  P <- nrow(d$pairs)
  dstar <- d$dstar
  eta <- if (itype == "bilinear") .draw_mat(draws, "eta", P) else NULL
  INT <- if (itype == "gpdi") .draw_mat(draws, "INT", P) else NULL
  kappa <- if (itype == "gpdi") .draw_mat(draws, "kappa", C) else NULL

  out <- vector("list", nrow(Dp))
  for (g in seq_len(nrow(Dp))) {
    dose <- Dp[g, ]
    lp <- rep(0, nrow(draws))
    for (c in seq_len(C))
      lp <- lp + emax[, c] * dose[c] / (ED50[, c] + dose[c])
    if (itype == "bilinear" && P > 0) {
      for (p in seq_len(P)) {
        a <- d$pairs[p, 1]; b <- d$pairs[p, 2]
        lp <- lp + eta[, p] * (dose[a] / dstar[a]) * (dose[b] / dstar[b])
      }
    } else if (itype == "gpdi" && P > 0) {
      for (p in seq_len(P)) {
        a <- d$pairs[p, 1]; b <- d$pairs[p, 2]
        lp <- lp + INT[, p] * (dose[a] / (kappa[, a] + dose[a])) *
          (dose[b] / (kappa[, b] + dose[b]))
      }
    }
    qs <- stats::quantile(lp, probs)
    out[[g]] <- c(mean = mean(lp), sd = stats::sd(lp), qs)
  }
  res <- as.data.frame(do.call(rbind, out))
  cbind(newdata[, comp, drop = FALSE], res)
}
