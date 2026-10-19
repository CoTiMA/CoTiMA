#' Map ctsem parameters to posterior columns
#'
#' @description
#' Creates labels and column mappings for ctsem population means,
#' time-independent predictor effects, individually varying parameters, and
#' their covariance parameters. The input may be either a fitted ctsem object
#' or a CoTiMA fit containing one pooled ctsem fit in \code{studyFitList}.
#'
#' @param fit A fitted ctsem object or a \code{"CoTiMAFit"} object. A CoTiMA
#'   object must contain exactly one pooled ctsem fit. CoTiMA initial-fit
#'   objects containing several separately fitted studies are ambiguous and
#'   must first be reduced to the desired ctsem fit.
#' @param est Optional result returned by \code{ctsem::ctExtract()}. If
#'   supplied, the dimensions of \code{rawpopmeans} and
#'   \code{tipredeffectparams} are checked against the model mapping.
#'
#' @return A list with data frames named \code{rawpopmeans},
#'   \code{tipredeffectparams}, \code{individually_varying}, and
#'   \code{random_effect_covariance}. The mapping includes the corresponding
#'   columns of \code{fit$stanfit$rawposterior} where applicable.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # A CoTiMA full-fit object:
#' parameterMap <- ctmaParMap(CoTiMAFullFit)
#'
#' # The pooled ctsem fit extracted from the same object:
#' ctsemFit <- CoTiMAFullFit$studyFitList
#' parameterMap <- ctmaParMap(ctsemFit)
#'
#' # Optional consistency checks against ctExtract():
#' estimates <- ctsem::ctExtract(ctsemFit)
#' parameterMap <- ctmaParMap(ctsemFit, est = estimates)
#' }
ctmaParMap <- function(fit, est = NULL) {
  fit <- .ctmaParMapFit(fit)
  pars <- fit$ctstanmodelbase$pars
  setup <- .ctmaParMapFirst(
    fit$data$matsetup,
    fit$setup$matsetup,
    fit$standata$matsetup
  )
  model_data <- .ctmaParMapFirst(fit$data, fit$standata)
  npar <- .ctmaParMapFirst(
    fit$data$nparams,
    fit$standata$nparams
  )

  if (is.null(pars) || is.null(setup) || is.null(npar)) {
    stop(
      "The selected ctsem fit does not contain the parameter setup ",
      "information required by ctmaParMap().",
      call. = FALSE
    )
  }
  npar <- as.integer(npar)[1L]
  posterior_layout <- .ctmaParMapRawposteriorLayout(model_data)

  rows <- lapply(seq_len(npar), function(i) {
    k <- which(setup[, "param"] == i)[1L]
    if (is.na(k)) stop(
      "No model-matrix entry found for raw parameter ", i, ".",
      call. = FALSE
    )

    matrix_name <- rownames(setup)[k]
    matrix_name <- sub("[0-9]+$", "", matrix_name)
    r <- unname(setup[k, "row"])
    cc <- unname(setup[k, "col"])
    internal_matrix <- matrix_name
    internal_row <- r
    internal_col <- cc
    iv_index <- unname(setup[k, "indvarying"])

    z <- which(pars$matrix == matrix_name & pars$row == r & pars$col == cc)
    z <- z[1L]
    label <- if (length(z) && !is.na(pars$param[z]) &&
        nzchar(pars$param[z])) {
      pars$param[z]
    } else {
      expanded <- fit$ctstanmodel$pars
      ez <- if (!is.null(expanded)) which(
        expanded$matrix == matrix_name & expanded$row == r &
          expanded$col == cc
      ) else integer()
      candidate <- if (length(ez)) expanded$param[ez[1L]] else NA_character_
      if (!is.na(candidate) && nzchar(candidate)) candidate else
        paste0(matrix_name, "_r", r, "_c", cc)
    }
    transform <- if (length(z)) pars$transform[z] else NA_character_

    # ctsem represents an individually varying manifest mean internally as an
    # augmented T0MEANS row. Map it back to its substantive MANIFESTMEANS row.
    original_iv <- which(!is.na(pars$indvarying) & pars$indvarying)
    if (iv_index > 0L && iv_index <= length(original_iv)) {
      iz <- original_iv[iv_index]
      matrix_name <- pars$matrix[iz]
      r <- pars$row[iz]
      cc <- pars$col[iz]
      if (!is.na(pars$param[iz]) && nzchar(pars$param[iz])) {
        label <- pars$param[iz]
      }
      transform <- pars$transform[iz]
    }

    data.frame(
      rawpopmeans_col = i,
      rawposterior_col = posterior_layout$rawpopmeans[i],
      matrix = matrix_name,
      row = r,
      col = cc,
      internal_matrix = internal_matrix,
      internal_row = internal_row,
      internal_col = internal_col,
      parameter = label,
      transform = transform,
      indvarying_index = iv_index,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  ti_names <- fit$ctstanmodelbase$TIpredNames
  ti_setup <- .ctmaParMapFirst(
    fit$data$TIPREDEFFECTsetup,
    fit$standata$TIPREDEFFECTsetup
  )
  if (length(ti_names) && length(ti_setup)) {
    ti <- do.call(rbind, lapply(seq_along(ti_names), function(j) {
      idx <- ti_setup[seq_len(npar), j]
      affected <- which(idx > 0L)
      data.frame(
        tipredeffectparams_col = idx[affected],
        rawposterior_col =
          posterior_layout$tipredeffectparams[idx[affected]],
        TI_predictor = ti_names[j],
        rawpopmeans_col = affected,
        parameter = out$parameter[affected],
        effect_label = paste0(
          ti_names[j], "_TO_", out$parameter[affected]
        ),
        stringsAsFactors = FALSE
      )
    }))
    ti <- ti[order(ti$tipredeffectparams_col), , drop = FALSE]
  } else {
    ti <- data.frame()
  }

  if (!is.null(est)) {
    if (is.null(est$rawpopmeans) || ncol(est$rawpopmeans) != nrow(out)) {
      stop(
        "est$rawpopmeans has a different number of columns than the model ",
        "setup.",
        call. = FALSE
      )
    }
    if (nrow(ti) && (is.null(est$tipredeffectparams) ||
        ncol(est$tipredeffectparams) != nrow(ti))) {
      stop(
        "est$tipredeffectparams has a different number of columns than its ",
        "model setup.",
        call. = FALSE
      )
    }
  }

  iv <- out[out$indvarying_index > 0L, c(
    "indvarying_index", "rawpopmeans_col", "matrix", "row", "col",
    "parameter"
  )]
  iv <- iv[order(iv$indvarying_index), , drop = FALSE]
  if (nrow(iv)) {
    iv$rawpopsdbase_col <- seq_len(nrow(iv))
    iv$rawposterior_mean_col <-
      posterior_layout$rawpopmeans[iv$rawpopmeans_col]
    iv$rawposterior_sd_col <-
      posterior_layout$rawpopsdbase[iv$rawpopsdbase_col]
    iv$sd_label <- paste0("SD_", iv$parameter)

    sqrt_index <- matrix(NA_integer_, nrow(iv), nrow(iv))
    counter <- 0L
    if (nrow(iv) > 1L) {
      for (cc in seq_len(nrow(iv) - 1L)) {
        for (r in seq.int(cc + 1L, nrow(iv))) {
          counter <- counter + 1L
          sqrt_index[r, cc] <- counter
        }
      }
    }
    covmap <- do.call(rbind, lapply(seq_len(nrow(iv)), function(r) {
      do.call(rbind, lapply(seq_len(r), function(cc) {
        sc <- sqrt_index[r, cc]
        data.frame(
          row = r,
          col = cc,
          row_indvarying_index = iv$indvarying_index[r],
          col_indvarying_index = iv$indvarying_index[cc],
          row_rawpopmeans_col = iv$rawpopmeans_col[r],
          col_rawpopmeans_col = iv$rawpopmeans_col[cc],
          row_rawpopsdbase_col = iv$rawpopsdbase_col[r],
          col_rawpopsdbase_col = iv$rawpopsdbase_col[cc],
          sqrtpcov_col = sc,
          rawposterior_sqrtpcov_col = if (is.na(sc)) NA_integer_ else
            posterior_layout$sqrtpcov[sc],
          row_parameter = iv$parameter[r],
          col_parameter = iv$parameter[cc],
          covariance_label = if (r == cc) {
            paste0("VAR_", iv$parameter[r])
          } else {
            paste0("COV_", iv$parameter[r], "_WITH_", iv$parameter[cc])
          },
          correlation_label = if (r == cc) {
            paste0("COR_", iv$parameter[r], "_WITH_SELF")
          } else {
            paste0("COR_", iv$parameter[r], "_WITH_", iv$parameter[cc])
          },
          stringsAsFactors = FALSE
        )
      }))
    }))
  } else {
    covmap <- data.frame()
  }

  structure(
    list(
      rawpopmeans = out,
      tipredeffectparams = ti,
      individually_varying = iv,
      random_effect_covariance = covmap
    ),
    ctsem_fit_source = attr(fit, "ctmaParMap_source", exact = TRUE)
  )
}

.ctmaParMapFirst <- function(...) {
  values <- list(...)
  for (value in values) if (!is.null(value)) return(value)
  NULL
}

.ctmaParMapIsCtsemFit <- function(x) {
  is.list(x) && !is.null(x$ctstanmodelbase$pars) &&
    (!is.null(x$data$matsetup) || !is.null(x$setup$matsetup) ||
      !is.null(x$standata$matsetup))
}

.ctmaParMapFit <- function(object) {
  if (.ctmaParMapIsCtsemFit(object)) {
    attr(object, "ctmaParMap_source") <- "ctsem fit"
    return(object)
  }

  # A hierarchical CoTiMA-RE fit retains its originating ctsem fit here.
  if (inherits(object, "CoTiMAStudyREFit") &&
      .ctmaParMapIsCtsemFit(object$prepared$original)) {
    fit <- object$prepared$original
    attr(fit, "ctmaParMap_source") <- "CoTiMA-RE fit"
    return(fit)
  }

  candidates <- object$studyFitList
  if (.ctmaParMapIsCtsemFit(candidates)) {
    attr(candidates, "ctmaParMap_source") <- "CoTiMA fit"
    return(candidates)
  }
  if (is.list(candidates)) {
    is_fit <- vapply(candidates, .ctmaParMapIsCtsemFit, logical(1L))
    fits <- candidates[is_fit]
    if (length(fits) == 1L) {
      fit <- fits[[1L]]
      attr(fit, "ctmaParMap_source") <- "CoTiMA fit"
      return(fit)
    }
    if (length(fits) > 1L) stop(
      "This CoTiMA object contains ", length(fits),
      " separately fitted ctsem models. ctmaParMap() requires one pooled ",
      "ctsem fit; supply the desired studyFitList element explicitly.",
      call. = FALSE
    )
  }

  stop(
    "fit must be either a fitted ctsem object or a CoTiMA fit object ",
    "containing exactly one pooled ctsem fit in $studyFitList.",
    call. = FALSE
  )
}

.ctmaParMapRawposteriorLayout <- function(data) {
  required <- c(
    "nparams", "nindvarying", "nindvaryingoffdiagonals", "intoverpop",
    "nsubjects", "ntipredeffects"
  )
  missing <- required[vapply(required, function(x) is.null(data[[x]]),
    logical(1L))]
  if (length(missing)) stop(
    "The ctsem fit lacks required data fields: ",
    paste(missing, collapse = ", "), ".",
    call. = FALSE
  )

  npar <- as.integer(data$nparams)[1L]
  niv <- as.integer(data$nindvarying)[1L]
  noff <- as.integer(data$nindvaryingoffdiagonals)[1L]
  nbase <- if (isTRUE(as.integer(data$intoverpop)[1L] == 0L)) {
    as.integer(data$nsubjects)[1L] * niv
  } else {
    0L
  }
  nti <- as.integer(data$ntipredeffects)[1L]

  start_sd <- npar + 1L
  start_cor <- npar + niv + 1L
  start_base <- npar + niv + noff + 1L
  start_ti <- npar + niv + noff + nbase + 1L
  list(
    rawpopmeans = seq_len(npar),
    rawpopsdbase = if (niv) seq.int(start_sd, length.out = niv) else integer(),
    sqrtpcov = if (noff) seq.int(start_cor, length.out = noff) else integer(),
    baseindparams = if (nbase) {
      seq.int(start_base, length.out = nbase)
    } else {
      integer()
    },
    tipredeffectparams = if (nti) {
      seq.int(start_ti, length.out = nti)
    } else {
      integer()
    }
  )
}
