#' Numerically determine optimal time lags
#'
#' @description Numerically searches a grid of time lags for the largest
#' absolute cross-effect implied by a continuous-time drift matrix. The OTL in
#' a standard CoTiMA summary uses the mean estimated drift matrix only. With
#' \code{CI = TRUE}, the search is repeated for transformed drift draws.
#'
#' @param ctmaFitFit A fit object created with \code{ctmaFit()}.
#' @param timeRange Numeric vector of candidate time lags. For a fit, the
#' default comprises 101 points from zero to three times the largest observed
#' lag.
#' @param driftMat A square drift matrix, list of such matrices, or a
#' three-dimensional array in draw-by-row-by-column order. Supply this or
#' \code{ctmaFitFit}, not both.
#' @param undoTimeScaling Use original-time units for a fit (default TRUE).
#' @param digits Number of digits used for rounding.
#' @param CI If TRUE, calculate draw-wise median OTLs and intervals.
#' @param confidence Confidence level; defaults to 0.95.
#'
#' @details For ML/MAP fits, intervals use the approximate transformed draws
#' stored by ctsem (normally Hessian-based); for Bayesian fits they are
#' posterior credible intervals. Precision is limited by \code{timeRange}.
#' An optimum on its boundary suggests using a wider grid.
#'
#' @examples
#' \dontrun{
#' OTL <- ctmaOTL(ctmaFitFit = CoTiMA::CoTiMAFullFit_6_new)
#' OTL_CI <- ctmaOTL(ctmaFitFit = CoTiMA::CoTiMAFullFit_6_new,
#'                   CI = TRUE, confidence = .95)
#' }
#'
#' @importFrom OpenMx expm
#' @export ctmaOTL
#'
#' @return A list with the grid, mean drift matrix, maximum absolute
#' cross-effects, and OTLs. With \code{CI = TRUE}, it also contains
#' \code{optimalCrossLagCI}, \code{maxCrossEffectCI}, \code{confidence},
#' \code{nDriftDraws}, and \code{uncertaintySource}.
ctmaOTL <- function(ctmaFitFit = NULL, timeRange = NULL, driftMat = NULL,
                    undoTimeScaling = NULL, digits = 4, CI = FALSE,
                    confidence = .95) {
  if (is.null(ctmaFitFit) == is.null(driftMat)) {
    stop("Supply exactly one of 'ctmaFitFit' or 'driftMat'.", call. = FALSE)
  }
  if (!is.logical(CI) || length(CI) != 1L || is.na(CI)) {
    stop("'CI' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      is.na(confidence) || confidence <= 0 || confidence >= 1) {
    stop("'confidence' must be strictly between zero and one.", call. = FALSE)
  }

  as_drift_array <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      if (!length(x) || !all(vapply(x, is.matrix, logical(1)))) {
        stop("Every element of 'driftMat' must be a matrix.", call. = FALSE)
      }
      dims <- vapply(x, function(z) paste(dim(z), collapse = "x"), character(1))
      if (length(unique(dims)) != 1L) {
        stop("All drift matrices must have identical dimensions.", call. = FALSE)
      }
      nr <- nrow(x[[1L]]); nc <- ncol(x[[1L]])
      out <- array(NA_real_, c(length(x), nr, nc))
      for (i in seq_along(x)) out[i, , ] <- x[[i]]
      return(out)
    }
    if (is.array(x) && length(dim(x)) == 3L) return(x)
    NULL
  }

  drift_draws <- NULL
  uncertainty_source <- NULL
  if (!is.null(ctmaFitFit)) {
    if (!inherits(ctmaFitFit, "CoTiMAFit")) {
      stop("'ctmaFitFit' must have been created with ctmaFit().", call. = FALSE)
    }
    n.latent <- ctmaFitFit$n.latent
    if (is.null(undoTimeScaling)) {
      undoTimeScaling <- TRUE
      message("'undoTimeScaling' was omitted; original-time units are used.")
    }
    scaleTime <- ctmaFitFit$argumentList$scaleTime
    if (is.null(scaleTime) || !is.finite(scaleTime)) scaleTime <- 1
    if (is.null(timeRange)) {
      maxDelta <- max(ctmaFitFit$statisticsList$allDeltas, na.rm = TRUE)
      if (maxDelta > 1) maxDelta <- round(maxDelta)
      timeRange <- seq(0, 3 * maxDelta, length.out = 101L)
      message("'timeRange' was omitted; using 0 to ",
              round(3 * maxDelta, digits), " in 101 grid points.")
    }
    driftMat <- if (isTRUE(undoTimeScaling)) {
      matrix(ctmaFitFit$modelResults$DRIFToriginal_time_scale,
             n.latent, n.latent, byrow = TRUE)
    } else {
      matrix(ctmaFitFit$modelResults$DRIFT,
             n.latent, n.latent, byrow = TRUE)
    }
    if (isTRUE(CI)) {
      drift_draws <- ctmaFitFit$studyFitList$stanfit$transformedpars$pop_DRIFT
      if (is.null(drift_draws) || length(dim(drift_draws)) != 3L) {
        stop("CI = TRUE requires transformed pop_DRIFT draws. Refit with uncertainty/transformed draws enabled, or supply multiple drift matrices.",
             call. = FALSE)
      }
      if (isTRUE(undoTimeScaling)) drift_draws <- drift_draws * scaleTime
      uncertainty_source <- if (isTRUE(ctmaFitFit$argumentList$optimize)) {
        "approximate transformed parameter draws"
      } else "posterior transformed parameter draws"
    }
  } else {
    drift_draws <- as_drift_array(driftMat)
    if (!is.null(drift_draws)) {
      if (dim(drift_draws)[2L] != dim(drift_draws)[3L]) {
        stop("Every drift matrix must be square.", call. = FALSE)
      }
      driftMat <- apply(drift_draws, c(2, 3), mean)
      uncertainty_source <- "supplied drift-matrix draws"
    }
    if (!is.matrix(driftMat) || nrow(driftMat) != ncol(driftMat)) {
      stop("'driftMat' must contain one or more square matrices.", call. = FALSE)
    }
    n.latent <- ncol(driftMat); scaleTime <- 1
    if (is.null(timeRange)) stop("Supply 'timeRange' with 'driftMat'.", call. = FALSE)
    if (isTRUE(CI) && is.null(drift_draws)) {
      stop("CI = TRUE needs multiple drift matrices.", call. = FALSE)
    }
  }
  if (!is.numeric(timeRange) || length(timeRange) < 2L ||
      any(!is.finite(timeRange)) || any(timeRange < 0)) {
    stop("'timeRange' needs at least two finite, non-negative lags.", call. = FALSE)
  }
  timeRange <- sort(unique(timeRange))

  calculate_otl <- function(drift) {
    lag <- effect <- matrix(NA_real_, n.latent, n.latent)
    for (row in seq_len(n.latent)) for (column in seq_len(n.latent)) {
      if (row != column && drift[row, column] != 0) {
        values <- vapply(timeRange, function(time) {
          OpenMx::expm(drift * time)[row, column]
        }, numeric(1))
        index <- which.max(abs(values))
        lag[row, column] <- timeRange[index]
        effect[row, column] <- abs(values[index])
      }
    }
    list(lag = lag, effect = effect)
  }

  point <- calculate_otl(driftMat)
  result <- list(
    timeRangeUsed = timeRange, driftMatrixUsed = driftMat,
    timeScaleUsed = scaleTime,
    maxCrossEffect = round(point$effect, digits),
    optimalCrossLag = round(point$lag, digits),
    pointEstimateNote = paste(
      "The point OTL is based on the mean estimated drift matrix.",
      "Use CI = TRUE to incorporate parameter uncertainty."
    )
  )
  if (isTRUE(CI)) {
    n_draws <- dim(drift_draws)[1L]
    lag_draws <- effect_draws <- array(NA_real_, c(n_draws, n.latent, n.latent))
    for (draw in seq_len(n_draws)) {
      one <- calculate_otl(drift_draws[draw, , ])
      lag_draws[draw, , ] <- one$lag
      effect_draws[draw, , ] <- one$effect
    }
    probs <- c((1 - confidence) / 2, .5, 1 - (1 - confidence) / 2)
    interval <- function(x) {
      out <- apply(x, c(2, 3), stats::quantile, probs = probs,
                   na.rm = TRUE, names = FALSE)
      aperm(out, c(2, 3, 1))
    }
    lag_ci <- interval(lag_draws); effect_ci <- interval(effect_draws)
    dimnames(lag_ci)[[3L]] <- dimnames(effect_ci)[[3L]] <-
      c("lower", "median", "upper")
    result$optimalCrossLagCI <- round(lag_ci, digits)
    result$maxCrossEffectCI <- round(effect_ci, digits)
    result$confidence <- confidence
    result$nDriftDraws <- n_draws
    result$uncertaintySource <- uncertainty_source
  }
  message("Drift matrix used for the point OTL (mean estimate):")
  print(driftMat)
  result
}
