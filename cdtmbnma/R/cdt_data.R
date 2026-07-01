#' Assemble a component dose-response network for cdtmbnma
#'
#' Turns an arm-level data frame into the design that \code{\link{cdt_fit}}
#' consumes. Each row is one treatment arm. Every component named in
#' \code{components} must have a dose column, with zero meaning the component is
#' absent from that arm.
#'
#' @param data A data frame with one row per study arm.
#' @param study Name of the column identifying the study.
#' @param components Character vector of component dose column names. A value of
#'   zero marks the component absent.
#' @param outcome Either \code{"continuous"} or \code{"binary"}.
#' @param y For a continuous outcome, the column of arm mean changes.
#' @param se For a continuous outcome, the column of standard errors. Supply this
#'   or both \code{sd} and \code{n}.
#' @param sd,n For a continuous outcome, columns of arm standard deviation and
#'   sample size, used as \code{se = sd / sqrt(n)} when \code{se} is absent.
#' @param events,n_binary For a binary outcome, the event-count and sample-size
#'   columns.
#' @param ref Optional column flagging the study-baseline arm (logical or 0/1).
#'   When absent, the all-zero-dose arm of each study is used, and a study with
#'   no such arm takes its lowest-total-dose arm as baseline.
#' @param dstar Optional named numeric vector of dose-normalisation references
#'   used by the bilinear surface. Defaults to the maximum observed dose of each
#'   component.
#' @param interactions Pairs of components to give an interaction term. Either
#'   \code{NULL} (default: all component pairs that co-occur in some arm), or a
#'   list of length-2 character vectors naming the pairs.
#'
#' @return An object of class \code{cdt_data}: a list holding the partial Stan
#'   data, the component and study labels, and the interaction pairs.
#' @export
cdt_data <- function(data, study, components,
                     outcome = c("continuous", "binary"),
                     y = NULL, se = NULL, sd = NULL, n = NULL,
                     events = NULL, n_binary = NULL,
                     ref = NULL, dstar = NULL, interactions = NULL) {
  outcome <- match.arg(outcome)
  stopifnot(is.data.frame(data))
  if (!study %in% names(data)) stop("Study column '", study, "' not found.")
  miss <- setdiff(components, names(data))
  if (length(miss)) stop("Component column(s) not found: ", paste(miss, collapse = ", "))

  D <- as.matrix(data[, components, drop = FALSE])
  storage.mode(D) <- "double"
  if (anyNA(D)) stop("Component dose columns must not contain NA (use 0 for absent).")
  if (any(D < 0)) stop("Component doses must be non-negative.")
  C <- ncol(D)
  N <- nrow(D)

  sf <- factor(data[[study]])
  study_idx <- as.integer(sf)
  S <- nlevels(sf)

  # outcome vectors
  yv <- rep(0, N); sev <- rep(1, N); rv <- rep(0L, N); nv <- rep(0L, N)
  if (outcome == "continuous") {
    if (is.null(y) || !y %in% names(data)) stop("Provide 'y' for a continuous outcome.")
    yv <- as.double(data[[y]])
    if (!is.null(se) && se %in% names(data)) {
      sev <- as.double(data[[se]])
    } else if (!is.null(sd) && !is.null(n) && all(c(sd, n) %in% names(data))) {
      sev <- as.double(data[[sd]]) / sqrt(as.double(data[[n]]))
    } else {
      stop("Provide 'se', or both 'sd' and 'n', for a continuous outcome.")
    }
    if (anyNA(yv) || anyNA(sev)) stop("Continuous outcome contains NA.")
    if (any(sev <= 0)) stop("Standard errors must be positive.")
  } else {
    if (is.null(events) || is.null(n_binary)) stop("Provide 'events' and 'n_binary' for a binary outcome.")
    rv <- as.integer(data[[events]])
    nv <- as.integer(data[[n_binary]])
    if (anyNA(rv) || anyNA(nv)) stop("Binary outcome contains NA.")
    if (any(rv < 0) || any(nv < rv)) stop("Need 0 <= events <= n_binary.")
  }

  # study-baseline (reference) flag: exactly one per study, no random effect
  if (!is.null(ref) && ref %in% names(data)) {
    is_ref <- as.integer(as.logical(data[[ref]]))
  } else {
    is_ref <- integer(N)
    tot <- rowSums(D)
    for (s in seq_len(S)) {
      rows <- which(study_idx == s)
      zero <- rows[tot[rows] == 0]
      base <- if (length(zero)) zero[1] else rows[which.min(tot[rows])]
      if (!length(zero))
        warning("Study '", levels(sf)[s], "' has no all-zero arm; using its lowest-dose arm as baseline.")
      is_ref[base] <- 1L
    }
  }

  # dose-normalisation references
  if (is.null(dstar)) {
    dstar <- apply(D, 2, max)
    dstar[dstar == 0] <- 1
  } else {
    dstar <- dstar[components]
    if (anyNA(dstar)) stop("'dstar' must name every component.")
  }
  names(dstar) <- components

  # interaction pairs (component index pairs that co-occur unless overridden)
  if (is.null(interactions)) {
    pa <- integer(0); pb <- integer(0)
    if (C >= 2) for (i in 1:(C - 1)) for (j in (i + 1):C) {
      if (any(D[, i] > 0 & D[, j] > 0)) { pa <- c(pa, i); pb <- c(pb, j) }
    }
  } else {
    pa <- integer(0); pb <- integer(0)
    for (pr in interactions) {
      i <- match(pr[1], components); j <- match(pr[2], components)
      if (anyNA(c(i, j))) stop("Interaction pair names must match 'components'.")
      pa <- c(pa, i); pb <- c(pb, j)
    }
  }
  P <- length(pa)

  standata <- list(
    N = N, S = S, C = C, P = P,
    study = study_idx, D = D, is_ref = is_ref,
    pair_a = as.array(pa), pair_b = as.array(pb),
    dstar = as.array(unname(dstar)),
    outcome = if (outcome == "continuous") 1L else 2L,
    y = yv, se = sev, r = rv, nn = nv
  )
  structure(list(
    standata = standata,
    components = components, studies = levels(sf),
    outcome = outcome, pairs = cbind(pair_a = pa, pair_b = pb),
    dstar = dstar, n_arms = N, n_studies = S
  ), class = "cdt_data")
}

#' @export
print.cdt_data <- function(x, ...) {
  cat("<cdt_data>\n")
  cat(sprintf("  arms:       %d across %d studies\n", x$n_arms, x$n_studies))
  cat(sprintf("  components: %s\n", paste(x$components, collapse = ", ")))
  cat(sprintf("  outcome:    %s\n", x$outcome))
  np <- nrow(x$pairs)
  if (np)
    cat(sprintf("  interaction pairs: %s\n",
                paste(apply(x$pairs, 1, function(p)
                  paste(x$components[p], collapse = " x ")), collapse = "; ")))
  else
    cat("  interaction pairs: none identifiable (no components co-occur)\n")
  invisible(x)
}
