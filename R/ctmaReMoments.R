#' Estimate study-level random-effects moments from a CoTiMA-FE fit
#'
#' @description
#' Reconstructs transformed study-specific parameters from the raw posterior
#' draws of a CoTiMA-FE model with study indicators and summarizes their
#' empirical between-study distribution. This is a post-hoc method-of-moments
#' analysis: it does not refit a hierarchical model and does not shrink study
#' estimates toward the population mean.
#'
#' @param fit A CoTiMA full-fit object containing one pooled ctsem fit, or that
#'   pooled ctsem fit itself.
#' @param study_design Optional numeric matrix with one row per study and one
#'   column for every time-independent predictor in the ctsem model. Column
#'   names should match \code{fit$ctstanmodelbase$TIpredNames}. If omitted,
#'   unique predictor rows are taken from the fitted data. Supply this argument
#'   when additional moderators produce more than one predictor row per study.
#' @param probs Two probabilities defining the lower and upper uncertainty
#'   limits.
#' @param weighting Study weighting for empirical population means and
#'   between-study SDs: \code{"equal"} or \code{"inverse_variance"}.
#' @param SEweighting Optional logical compatibility shortcut. \code{TRUE}
#'   selects inverse-variance weighting and \code{FALSE} selects equal
#'   weighting. When supplied, it overrides \code{weighting}.
#' @param SEsource Source label and validation policy for study-specific SEs:
#'   \code{"auto"}, \code{"posterior"}, or \code{"hessian"}. All SEs are
#'   calculated as the SD across the transformed study draws in
#'   \code{stanfit$rawposterior}. With \code{"auto"}, the function uses the
#'   original fit's optimization flag to identify posterior versus
#'   Hessian/optimization-based draws when possible.
#' @param cores Number of cores used by \code{ctsem::ctExtract()} when
#'   summaries of individually varying parameters are requested.
#' @param include_matrices Reconstruct and summarize study-specific model
#'   matrices.
#' @param include_individual Summarize person-level random-effect covariance
#'   parameters when the fitted model contains individually varying parameters.
#'
#' @return An object inheriting from \code{"CoTiMAFit"} and
#'   \code{"CoTiMAReMoments"}. Its main elements are
#'   \code{studyParameters}, \code{distributionParameters},
#'   \code{studySE}, \code{studyWeights}, \code{effectiveNumberOfStudies},
#'   \code{covarianceDraws}, \code{correlationDraws},
#'   \code{studyMatrices}, and \code{matrixSummary}.
#'
#' @details
#' For posterior draw \eqn{d}, study \eqn{s}, and parameter \eqn{p}, the
#' function first calculates the raw study value
#' \deqn{\eta_{sp}^{(d)}=\mu_p^{(d)}+\sum_j x_{sj}\beta_{pj}^{(d)}}
#' and then applies the parameter-specific ctsem transformation
#' \eqn{\theta_{sp}^{(d)}=g_p(\eta_{sp}^{(d)})}. With equal weighting, every
#' study receives the same weight and the between-study variance uses the
#' ordinary sample-variance denominator \eqn{K-1}.
#' With inverse-variance weighting, parameter-specific fixed weights are
#' estimated from the across-draw SD of each transformed study parameter. The
#' weighted sample-variance denominator is
#' \eqn{\sum w-(\sum w^2)/(\sum w)}. Precision weighting changes the estimand
#' and is best treated as a sensitivity analysis.
#'
#' @export
#' @seealso \code{\link{ctmaParMap}}, \code{\link{ctmaRePrep}},
#'   \code{\link{ctmaReFit}}
#'
#' @examples
#' \dontrun{
#' moments <- ctmaReMoments(CoTiMAFullFit, cores = 4)
#' summary(moments)
#' moments$studyParameters
#' moments$distributionParameters
#' }
ctmaReMoments <- function(
    fit,
    study_design = NULL,
    probs = c(.025, .975),
    weighting = c("equal", "inverse_variance"),
    SEweighting = NULL,
    SEsource = c("auto", "posterior", "hessian"),
    cores = 2,
    include_matrices = TRUE,
    include_individual = TRUE) {

  if (!is.numeric(probs) || length(probs) != 2L || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1) || probs[1L] >= probs[2L]) {
    stop("probs must contain two increasing probabilities between zero and one.",
      call. = FALSE)
  }
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) stop("cores must be a positive integer.",
    call. = FALSE)

  weighting <- match.arg(weighting)
  SEsource <- match.arg(SEsource)
  if (!is.null(SEweighting)) {
    if (!is.logical(SEweighting) || length(SEweighting) != 1L ||
        is.na(SEweighting)) stop(
      "SEweighting must be TRUE, FALSE, or NULL.", call. = FALSE)
    weighting <- if (isTRUE(SEweighting)) "inverse_variance" else "equal"
  }

  ctfit <- .ctmaParMapFit(fit)
  if (is.null(ctfit$stanfit$rawposterior)) stop(
    "The pooled ctsem fit has no fit$stanfit$rawposterior matrix. Moment-based ",
    "estimation requires posterior draws from the fitted CoTiMA-FE model.",
    call. = FALSE
  )

  .ctma_print_message(paste0(
    "Computing post-hoc study-level moments from the CoTiMA-FE posterior. ",
    if (weighting == "equal") "Studies are equally weighted. " else
      "Study effects use parameter-specific inverse-variance weights. ",
    "No hierarchical shrinkage is applied."
  ))

  parameter_map <- ctmaParMap(ctfit)
  design <- .ctmaReMomentsDesign(ctfit, study_design)
  draws <- .ctmaReMomentsStudyDraws(ctfit, parameter_map, design)
  se_source <- .ctmaReMomentsSESource(ctfit, SEsource)
  moments <- .ctmaReMomentsDistribution(
    draws$transformed, probs, weighting = weighting
  )
  if (weighting == "inverse_variance") .ctma_print_message(paste0(
    "Inverse-variance study weighting uses parameter-specific SEs from ",
    se_source, ". This changes the estimand and is recommended as a ",
    "sensitivity analysis alongside equal weighting. Covariance and ",
    "correlation draws remain equally weighted because parameter-specific ",
    "weights do not define a unique joint covariance weight."
  ))

  study_matrices <- matrix_summary <- NULL
  if (isTRUE(include_matrices)) {
    study_matrices <- .ctmaReMomentsMatrices(ctfit, draws$transformed,
      parameter_map$rawpopmeans)
    matrix_summary <- .ctmaReMomentsMatrixSummary(study_matrices, probs)
  }

  individual <- NULL
  if (isTRUE(include_individual) &&
      nrow(parameter_map$individually_varying)) {
    .ctma_print_message(paste0(
      "Extracting ctsem constrained draws for the individually varying ",
      "parameter summary."
    ))
    extracted <- ctsem::ctExtract(ctfit, cores = cores)
    individual <- .ctmaReMomentsIndividual(extracted, parameter_map, probs)
  }

  out <- structure(list(
    method = "post-hoc empirical moments",
    weighting = weighting,
    covarianceWeighting = "equal study weights",
    SEsource = se_source,
    pooledFit = ctfit,
    parameterMap = parameter_map,
    studyDesign = design,
    rawStudyDraws = draws$raw,
    transformedStudyDraws = draws$transformed,
    studyParameters = moments$study,
    distributionParameters = moments$distribution,
    studySE = moments$study_se,
    studyWeights = moments$weights,
    effectiveNumberOfStudies = moments$effective_studies,
    covarianceDraws = moments$covariance,
    correlationDraws = moments$correlation,
    studyMatrices = study_matrices,
    matrixSummary = matrix_summary,
    individualRandomEffects = individual,
    probs = probs
  ), class = c("CoTiMAFit", "CoTiMAReMoments"))
  out$summary <- .ctmaReMomentsSummary(out, probs = probs, print_note = FALSE)
  out
}

.ctmaReMomentsData <- function(fit) {
  .ctmaParMapFirst(fit$data, fit$standata)
}

.ctmaReMomentsDesign <- function(fit, study_design = NULL) {
  ti_names <- fit$ctstanmodelbase$TIpredNames
  if (!length(ti_names)) stop(
    "The fitted CoTiMA-FE model has no time-independent predictors from which ",
    "study-specific parameters can be reconstructed.", call. = FALSE)
  if (!is.null(study_design)) {
    x <- as.matrix(study_design)
    storage.mode(x) <- "double"
    if (!is.null(colnames(x))) {
      missing <- setdiff(ti_names, colnames(x))
      if (length(missing)) stop(
        "study_design lacks predictor columns: ",
        paste(missing, collapse = ", "), ".", call. = FALSE)
      x <- x[, ti_names, drop = FALSE]
    }
    if (ncol(x) != length(ti_names)) stop(
      "study_design must have one column for every time-independent predictor.",
      call. = FALSE)
  } else {
    d <- .ctmaReMomentsData(fit)
    if (is.null(d$tipredsdata)) stop(
      "The fitted object has no tipredsdata. Supply study_design explicitly.",
      call. = FALSE)
    x <- as.matrix(unique(as.data.frame(d$tipredsdata[, ti_names,
      drop = FALSE])))
    storage.mode(x) <- "double"
  }
  if (!nrow(x)) stop("study_design contains no studies.", call. = FALSE)
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) {
    rownames(x) <- paste0("Study_", seq_len(nrow(x)))
  }
  if (anyNA(x) || any(!is.finite(x))) stop(
    "study_design must contain only finite, non-missing values.",
    call. = FALSE)
  x
}

.ctmaReMomentsSoftplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

.ctmaReMomentsTransform <- function(x, transform) {
  if (is.na(transform) || !nzchar(transform) || identical(transform, "param"))
    return(x)
  eval(parse(text = transform), envir = list(
    param = x, log1p_exp = .ctmaReMomentsSoftplus,
    exp = exp, log = log, log1p = log1p, sqrt = sqrt
  ))
}

.ctmaReMomentsStudyDraws <- function(fit, map, design) {
  posterior <- as.matrix(fit$stanfit$rawposterior)
  n_draw <- nrow(posterior)
  n_par <- nrow(map$rawpopmeans)
  n_study <- nrow(design)
  n_ti <- ncol(design)
  d <- .ctmaReMomentsData(fit)

  beta <- array(0, c(n_draw, n_par, n_ti))
  for (j in seq_len(n_ti)) {
    idx <- d$TIPREDEFFECTsetup[seq_len(n_par), j]
    affected <- which(idx > 0L)
    if (length(affected)) {
      posterior_cols <- map$tipredeffectparams$rawposterior_col[
        match(idx[affected], map$tipredeffectparams$tipredeffectparams_col)]
      if (anyNA(posterior_cols)) stop(
        "A TI-effect parameter could not be matched to rawposterior.",
        call. = FALSE)
      beta[, affected, j] <- posterior[, posterior_cols, drop = FALSE]
    }
  }

  dn <- list(
    draw = NULL,
    parameter = map$rawpopmeans$parameter,
    study = rownames(design)
  )
  raw <- transformed <- array(NA_real_, c(n_draw, n_par, n_study),
    dimnames = dn)
  raw_means <- posterior[, map$rawpopmeans$rawposterior_col, drop = FALSE]
  for (s in seq_len(n_study)) {
    raw[, , s] <- raw_means
    for (j in seq_len(n_ti)) {
      raw[, , s] <- raw[, , s] + beta[, , j] * design[s, j]
    }
    for (p in seq_len(n_par)) {
      transformed[, p, s] <- .ctmaReMomentsTransform(
        raw[, p, s], map$rawpopmeans$transform[p]
      )
    }
  }
  list(raw = raw, transformed = transformed)
}

.ctmaReMomentsStats <- function(x, probs) {
  c(mean = mean(x), median = stats::median(x),
    lower = stats::quantile(x, probs[1L], names = FALSE),
    upper = stats::quantile(x, probs[2L], names = FALSE))
}

.ctmaReMomentsSESource <- function(fit, requested) {
  optimized <- fit$args$optimize
  identified <- if (isTRUE(optimized)) {
    "Hessian/optimization-based rawposterior draws"
  } else if (identical(optimized, FALSE)) {
    "Bayesian posterior rawposterior draws"
  } else {
    "rawposterior draws (origin not identified in the fit object)"
  }
  if (requested == "posterior" && isTRUE(optimized)) stop(
    "SEsource='posterior' was requested, but the fit records optimize=TRUE. ",
    "Use SEsource='hessian' or 'auto', or provide a Bayesian CoTiMA-FE fit.",
    call. = FALSE
  )
  if (requested == "hessian" && identical(optimized, FALSE)) stop(
    "SEsource='hessian' was requested, but the fit records optimize=FALSE. ",
    "Use SEsource='posterior' or 'auto', or provide an optimization-based fit.",
    call. = FALSE
  )
  if (requested == "posterior") "Bayesian posterior rawposterior draws" else
    if (requested == "hessian")
      "Hessian/optimization-based rawposterior draws" else identified
}

.ctmaReMomentsDistribution <- function(x, probs,
    weighting = c("equal", "inverse_variance")) {
  weighting <- match.arg(weighting)
  dn <- dimnames(x)
  k <- dim(x)[3L]
  p_count <- dim(x)[2L]
  study_se <- matrix(NA_real_, nrow = k, ncol = p_count,
    dimnames = list(study = dn[[3L]], parameter = dn[[2L]]))
  for (s in seq_len(k)) for (p in seq_len(p_count)) {
    study_se[s, p] <- stats::sd(x[, p, s])
  }
  if (any(!is.finite(study_se))) stop(
    "At least one study-specific SE is non-finite. Precision weights cannot ",
    "be calculated from these draws.", call. = FALSE)
  if (weighting == "inverse_variance" && any(study_se <= 0)) stop(
    "At least one study-specific SE is zero. Inverse-variance weighting would ",
    "give it infinite weight; use equal weighting or inspect the fitted draws.",
    call. = FALSE)
  weights <- if (weighting == "equal") {
    matrix(1 / k, nrow = k, ncol = p_count, dimnames = dimnames(study_se))
  } else {
    w <- 1 / study_se^2
    sweep(w, 2L, colSums(w), "/")
  }
  effective_studies <- 1 / colSums(weights^2)

  study <- do.call(rbind, lapply(seq_len(dim(x)[3L]), function(s) {
    do.call(rbind, lapply(seq_len(dim(x)[2L]), function(p) {
      z <- .ctmaReMomentsStats(x[, p, s], probs)
      data.frame(study = dn[[3L]][s], parameter = dn[[2L]][p],
        standard_error = study_se[s, p], weight = weights[s, p],
        mean = z[1L], median = z[2L], lower = z[3L], upper = z[4L],
        row.names = NULL)
    }))
  }))

  re_mean <- re_sd <- matrix(NA_real_, nrow = dim(x)[1L], ncol = p_count)
  for (p in seq_len(p_count)) {
    wp <- weights[, p]
    xp <- matrix(x[, p, ], nrow = dim(x)[1L], ncol = k)
    re_mean[, p] <- as.vector(xp %*% wp)
    denominator <- 1 - sum(wp^2)
    if (k > 1L && denominator > .Machine$double.eps) {
      centered <- sweep(xp, 1L, re_mean[, p], "-")
      re_sd[, p] <- sqrt(as.vector(centered^2 %*% wp) / denominator)
    }
  }
  distribution <- do.call(rbind, lapply(seq_len(dim(x)[2L]), function(p) {
    a <- .ctmaReMomentsStats(re_mean[, p], probs)
    b <- .ctmaReMomentsStats(re_sd[, p], probs)
    rbind(
      data.frame(parameter = dn[[2L]][p],
        distribution_parameter = "mean", mean = a[1L], median = a[2L],
        lower = a[3L], upper = a[4L], weighting = weighting,
        effective_studies = effective_studies[p], row.names = NULL),
      data.frame(parameter = dn[[2L]][p],
        distribution_parameter = "between_study_sd", mean = b[1L],
        median = b[2L], lower = b[3L], upper = b[4L],
        weighting = weighting, effective_studies = effective_studies[p],
        row.names = NULL)
    )
  }))

  covariance <- correlation <- array(NA_real_, c(
    dim(x)[1L], dim(x)[2L], dim(x)[2L]
  ), dimnames = list(NULL, dn[[2L]], dn[[2L]]))
  if (dim(x)[3L] > 1L) for (i in seq_len(dim(x)[1L])) {
    study_by_parameter <- t(matrix(
      x[i, , ], nrow = dim(x)[2L], ncol = dim(x)[3L]
    ))
    covariance[i, , ] <- stats::cov(study_by_parameter)
    correlation[i, , ] <- suppressWarnings(stats::cor(study_by_parameter))
  }
  list(study = study, distribution = distribution,
    study_se = study_se, weights = weights,
    effective_studies = effective_studies,
    covariance = covariance, correlation = correlation)
}

.ctmaReMomentsMatrices <- function(fit, draws, map) {
  ans <- list()
  d <- .ctmaReMomentsData(fit)
  for (m in unique(map$matrix)) {
    ii <- which(map$matrix == m)
    nr <- max(map$row[ii]); nc <- max(map$col[ii])
    a <- array(0, c(dim(draws)[1L], nr, nc, dim(draws)[3L]),
      dimnames = list(NULL, seq_len(nr), seq_len(nc), dimnames(draws)[[3L]]))
    for (k in ii) a[, map$row[k], map$col[k], ] <- draws[, k, ]
    ans[[m]] <- a

    if (m %in% c("DIFFUSION", "T0VAR", "MANIFESTVAR") && nr == nc) {
      cv <- a
      for (s in seq_len(dim(draws)[3L])) for (draw in seq_len(dim(draws)[1L])) {
        mat <- matrix(a[draw, , , s], nr, nc)
        if (isTRUE(d$choleskymats >= 1)) {
          cv[draw, , , s] <- tcrossprod(mat)
        } else {
          ss <- sv <- numeric(nr)
          for (i in seq_len(nr)) {
            other <- vapply(setdiff(seq_len(nr), i), function(j)
              if (j > i) mat[j, i] else mat[i, j], numeric(1L))
            ss[i] <- sum(other^2) + 1e-5
            sv[i] <- sum(other) + 1e-5
          }
          o <- matrix(0, nr, nr)
          for (i in seq_len(nr)) {
            r1 <- sqrt(ss[i])
            r3 <- abs(sv[i]) / r1 - 1
            r4 <- sqrt(.ctmaReMomentsSoftplus(
              2 * (abs(sv[i]) - sv[i] - 1) - 4
            ))
            rr <- sqrt(ss[i] + ((r4 * r3 + 1) * r4 + 1))
            for (j in seq_len(nr)) {
              if (j > i) o[i, j] <- mat[j, i] / rr
              if (j < i) o[i, j] <- mat[i, j] / rr
            }
            o[i, i] <- sqrt(1 - sum(o[i, ]^2) + 1e-5)
          }
          cv[draw, , , s] <- tcrossprod(diag(diag(mat), nr, nr) %*% o)
        }
      }
      ans[[paste0(m, "cov")]] <- cv
    }
  }
  ans
}

.ctmaReMomentsMatrixSummary <- function(matrices, probs) {
  if (!length(matrices)) return(data.frame())
  do.call(rbind, lapply(names(matrices), function(nm) {
    a <- matrices[[nm]]
    do.call(rbind, lapply(seq_len(dim(a)[4L]), function(s) {
      do.call(rbind, lapply(seq_len(dim(a)[2L]), function(r) {
        do.call(rbind, lapply(seq_len(dim(a)[3L]), function(cc) {
          z <- .ctmaReMomentsStats(a[, r, cc, s], probs)
          data.frame(matrix = nm, study = dimnames(a)[[4L]][s], row = r,
            col = cc, mean = z[1L], median = z[2L], lower = z[3L],
            upper = z[4L], row.names = NULL)
        }))
      }))
    }))
  }))
}

.ctmaReMomentsIndividual <- function(est, map, probs) {
  iv <- map$individually_varying
  labels <- iv$parameter
  covariance <- est$popcov
  if (is.null(covariance)) return(list(map = iv, summary = data.frame()))
  dimnames(covariance)[2:3] <- list(labels, labels)
  correlation <- covariance
  for (i in seq_len(dim(covariance)[1L])) {
    correlation[i, , ] <- suppressWarnings(stats::cov2cor(covariance[i, , ]))
  }
  rows <- do.call(rbind, lapply(seq_len(nrow(iv)), function(r) {
    do.call(rbind, lapply(seq_len(r), function(cc) {
      value <- if (r == cc) sqrt(covariance[, r, r]) else covariance[, r, cc]
      z <- .ctmaReMomentsStats(value, probs)
      data.frame(quantity = if (r == cc) "sd" else "covariance",
        row_parameter = labels[r], col_parameter = labels[cc],
        mean = z[1L], median = z[2L], lower = z[3L], upper = z[4L])
    }))
  }))
  list(map = iv, covariance = covariance, correlation = correlation,
    T0VAR = est$pop_T0VAR, T0cov = est$pop_T0cov, summary = rows)
}

.ctmaReMomentsSummary <- function(object, probs = object$probs,
    print_note = TRUE) {
  if (!inherits(object, "CoTiMAReMoments")) stop(
    "object must be returned by ctmaReMoments().", call. = FALSE)
  if (isTRUE(print_note)) .ctma_print_message(paste0(
    "This is a post-hoc empirical moment summary with ",
    if (object$weighting == "equal") "equal study weights. " else
      "parameter-specific inverse-variance study weights. ",
    "It is not a hierarchical fit and contains no partial pooling."
  ))
  list(
    method = object$method,
    weighting = object$weighting,
    SEsource = object$SEsource,
    covarianceWeighting = object$covarianceWeighting,
    numberOfStudies = nrow(object$studyDesign),
    studyParameters = object$studyParameters,
    distributionParameters = object$distributionParameters,
    studySE = object$studySE,
    studyWeights = object$studyWeights,
    effectiveNumberOfStudies = object$effectiveNumberOfStudies,
    matrixSummary = object$matrixSummary,
    individualRandomEffects = object$individualRandomEffects,
    probs = probs
  )
}
