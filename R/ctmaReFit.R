# Study-level random effects for ctsem models returned with fit = FALSE.
#
# This hierarchical fitting script intentionally defines only
# ctmaRePrep(), ctmaReFit(), and summary().
# Functions beginning with ct_ (for example ct_random_effect_summary()) belong
# to the separate post-hoc FI reconstruction script
# CoTiMA_RE_moments.R; they do not fit the hierarchical RE model.
#
# CRAN/API NOTE:
# MAP and conditional ML use the exported rstan::optimizing() interface. No
# unexported ctsem optimizer is required. A future exported ctsem wrapper could
# still be considered if ctsem-specific optimizer uncertainty is restored.
# Historical maintainer reminder requested for future development:
# "Ask the ctsem maintainer to export a supported optimization function or
# provide a public wrapper around stanoptimis()."
#
# This prototype modifies the generated ctsem Stan program. Study effects are
# added on ctsem's unconstrained/raw scale and are shared by every person in a
# study. Person-level random effects remain a separate level.

.ctma_print_message <- function(text, width = 81L,
    prefix = "Just a note: ") {
  width <- max(40L, as.integer(width)[1L])
  left <- "######## "
  right <- " #########"
  content_width <- width - nchar(left) - nchar(right)
  text <- paste0(prefix, paste(as.character(text), collapse = ""))
  lines <- unlist(strwrap(text, width = content_width, simplify = FALSE),
    use.names = FALSE)
  if (!length(lines)) lines <- ""
  border <- strrep("#", width)
  print(border)
  for (line in lines) {
    padding <- content_width - nchar(line, type = "width")
    print(paste0(left, line, strrep(" ", max(0L, padding)), right))
  }
  print(border)
  invisible(NULL)
}

.ctma_replace_once <- function(text, pattern, replacement, label = pattern) {
  hit <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  n <- if (identical(hit[1L], -1L)) 0L else length(hit)
  if (n != 1L) stop("Expected exactly one Stan insertion point for ", label,
                    "; found ", n, ". The ctsem template may have changed.")
  sub(pattern, replacement, text, fixed = TRUE)
}

.ctma_progress <- function(verbose, stage, total, text, started = NULL) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  elapsed <- if (is.null(started)) "" else
    paste0("; elapsed ", format(round(as.numeric(difftime(
      Sys.time(), started, units = "secs")), 1), trim = TRUE), " s")
  .ctma_print_message(sprintf(
    "CoTiMA study RE %d/%d: %s%s", stage, total, text, elapsed
  ))
  invisible(NULL)
}

.ctma_add_draw_dimension <- function(x) {
  lapply(x, function(z) {
    olddim <- dim(z)
    oldnames <- dimnames(z)
    if (is.null(olddim)) {
      array(z, dim = c(1L, length(z)),
        dimnames = list(draw = NULL, names(z)))
    } else {
      array(z, dim = c(1L, olddim),
        dimnames = c(list(draw = NULL), oldnames))
    }
  })
}

.ctma_optimizer_data <- function(fit_data, estimation) {
  if (!estimation %in% c("map", "ml")) stop(
    "Optimizer data can only be configured for estimation='map' or 'ml'."
  )
  if (is.null(fit_data$priors)) stop(
    "The ctsem Stan data do not contain the expected 'priors' switch."
  )
  enabled <- as.integer(estimation == "map")
  fit_data$priors <- enabled
  fit_data$studyre_use_hyperpriors <- enabled
  fit_data
}

# rstan distinguishes a scalar from a Stan array/vector by the R object's
# dimension attribute.  A plain atomic vector of length one can therefore be
# interpreted as a scalar, which breaks models with exactly one study-level
# random effect.  Always retain an explicit one-dimensional shape.
.ctma_stan_1d <- function(x, mode = c("integer", "double")) {
  mode <- match.arg(mode)
  x <- if (mode == "integer") as.integer(x) else as.numeric(x)
  array(x, dim = length(x))
}

.ctma_bool_matrix <- function(x, n, name) {
  if (is.null(x)) return(matrix(FALSE, n, n))
  x <- as.matrix(x)
  if (!identical(dim(x), c(n, n)))
    stop(name, " must be a ", n, " x ", n, " matrix.")
  if (anyNA(x)) stop(name, " may not contain NA.")
  matrix(as.logical(x), n, n, dimnames = dimnames(x))
}

.ctma_named_ci <- function(x, name) {
  if (is.null(x) || is.null(names(x))) return(NULL)
  j <- match(tolower(name), tolower(names(x)))
  if (is.na(j)) NULL else x[[j]]
}

.ctma_infer_study_id <- function(model, study_dummy_ti = NULL,
    verbose = TRUE) {
  ti_names <- model$ctstanmodelbase$TIpredNames
  x <- as.data.frame(model$standata$tipredsdata)
  if (ncol(x) != length(ti_names))
    stop("The TI predictor data and TIpredNames have incompatible dimensions.")
  names(x) <- ti_names

  if (!is.null(study_dummy_ti)) {
    if (is.numeric(study_dummy_ti)) {
      if (anyNA(study_dummy_ti) || any(study_dummy_ti < 1L) ||
          any(study_dummy_ti > length(ti_names)))
        stop("Numeric study_dummy_ti contains an invalid TI column index.")
      study_dummy_ti <- ti_names[as.integer(study_dummy_ti)]
    }
    missing <- setdiff(study_dummy_ti, ti_names)
    if (length(missing)) stop("Unknown study dummy TI column(s): ",
      paste(missing, collapse = ", "))
    groups <- list(explicit = study_dummy_ti)
  } else {
    pars <- model$ctstanmodelbase$pars
    eligible <- vapply(ti_names, function(nm) {
      effect <- paste0(nm, "_effect")
      effect %in% names(pars) && any(pars[[effect]], na.rm = TRUE)
    }, logical(1))
    candidates <- ti_names[eligible]
    if (!length(candidates)) stop(
      "No 'study' TI predictor was found and no TI predictors with enabled ",
      "parameter effects were available for inferring study dummies. Supply ",
      "study_id or study_dummy_ti explicitly."
    )

    # Unscaled dummy columns contain two levels (usually 0/1), while weighted
    # effect coding can contain three. Centring/scaling applies an affine
    # transformation but preserves the number and row pattern of these levels,
    # so recognize arbitrary finite two- or three-level values. `signif` avoids
    # treating negligible floating-point differences as additional levels.
    # Group columns by their complete parameter-effect pattern so a binary
    # substantive moderator with a different target is not silently absorbed.
    candidates <- candidates[vapply(x[candidates], function(z) {
      if (anyNA(z) || any(!is.finite(z))) return(FALSE)
      nlevels <- length(unique(signif(as.numeric(z), 12L)))
      nlevels %in% c(2L, 3L)
    }, logical(1))]
    if (!length(candidates)) stop(
      "No 'study' TI predictor was found and the enabled TI predictors do not ",
      "look like two- or three-level (possibly scaled) study dummies. ",
      "Supply study_id or ",
      "study_dummy_ti explicitly."
    )
    signatures <- vapply(candidates, function(nm)
      paste(as.integer(pars[[paste0(nm, "_effect")]]), collapse = ""), "")
    groups <- split(candidates, signatures)
  }

  valid <- lapply(groups, function(cols) {
    profiles <- unique(as.data.frame(lapply(
      x[cols], function(z) signif(as.numeric(z), 12L)
    )))
    if (nrow(profiles) == length(cols) + 1L) cols else NULL
  })
  valid <- Filter(Negate(is.null), valid)
  if (!length(valid)) stop(
    "TI dummy inference failed: no candidate group had K - 1 columns and K ",
    "distinct row profiles. Supply study_id or study_dummy_ti explicitly."
  )
  max_size <- max(lengths(valid))
  valid <- valid[lengths(valid) == max_size]
  if (length(valid) != 1L) stop(
    "TI dummy inference is ambiguous. Candidate sets: ",
    paste(vapply(valid, paste, collapse = ",", FUN.VALUE = ""),
      collapse = " | "),
    ". Supply study_dummy_ti explicitly."
  )

  cols <- valid[[1L]]
  profile_data <- lapply(x[cols], function(z) signif(as.numeric(z), 12L))
  profile_key <- do.call(paste, c(profile_data, sep = "\r"))
  study_id <- match(profile_key, unique(profile_key))
  if (isTRUE(verbose)) .ctma_print_message(paste0(
    "No TI predictor named 'study' was found. Constructed study membership ",
    "for ", length(unique(study_id)), " studies from TI dummy columns: ",
    paste(cols, collapse = ", "),
    ". Scaled and unscaled dummy values are both supported."
  ))
  list(id = study_id, dummy_names = cols)
}

.ctma_stored_re <- function(model) {
  al <- model$argumentList
  block <- .ctma_named_ci(al, "randomEffect")
  if (is.null(block)) block <- .ctma_named_ci(al, "randomEffects")
  list(
    reDRIFT = .ctma_named_ci(block, "reDRIFT"),
    reDIFFUSION = .ctma_named_ci(block, "reDIFFUSION"),
    reT0VAR = .ctma_named_ci(block, "reT0VAR")
  )
}

.ctma_study_id <- function(model, study_id = NULL, study_id_ti = NULL) {
  ns <- as.integer(model$standata$nsubjects)
  subject <- model$standata$subject

  if (is.null(study_id)) {
    ti_names <- model$ctstanmodelbase$TIpredNames
    if (is.null(study_id_ti))
      stop("Supply study_id or the name/index of its TI column in study_id_ti.")
    j <- if (is.character(study_id_ti)) match(study_id_ti, ti_names) else
      as.integer(study_id_ti)
    if (length(j) != 1L || is.na(j) || j < 1L || j > length(ti_names))
      stop("study_id_ti does not identify a time-independent predictor.")
    study_id <- model$standata$tipredsdata[, j]
  }

  if (length(study_id) == length(subject)) {
    by_subject <- split(study_id, subject)
    bad <- vapply(by_subject, function(z) length(unique(z)) != 1L, logical(1))
    if (any(bad)) stop("study_id must be constant within person-level id.")
    study_id <- vapply(by_subject, `[`, numeric(1), 1L)
  }
  if (length(study_id) != ns)
    stop("study_id must have length nsubjects or ndatapoints.")
  if (anyNA(study_id)) stop("study_id may not contain missing values.")

  study_labels <- unique(as.character(study_id))
  encoded <- match(as.character(study_id), study_labels)
  if (length(study_labels) < 2L)
    stop("Study-level random effects require at least two distinct studies.")
  list(id = as.integer(encoded), labels = study_labels,
       nstudies = length(study_labels))
}

.ctma_raw_layout <- function(model) {
  d <- model$standata
  npar <- as.integer(d$nparams); niv <- as.integer(d$nindvarying)
  noff <- as.integer(d$nindvaryingoffdiagonals)
  nbase <- if (as.integer(d$intoverpop) == 0L)
    as.integer(d$nsubjects) * niv else 0L
  nti <- as.integer(d$ntipredeffects)
  list(
    rawpopmeans = seq_len(npar),
    tipredeffectparams = if (nti) seq.int(
      npar + niv + noff + nbase + 1L, length.out = nti
    ) else integer()
  )
}

.ctma_raw_point <- function(model) {
  rp <- model$stanfit$rawposterior
  source <- NULL
  if (!is.null(rp) && length(rp)) {
    point <- apply(as.matrix(rp), 2L, stats::median, na.rm = TRUE)
    source <- "posterior medians from fit$stanfit$rawposterior"
  } else if (!is.null(model$stanfit$rawest) && length(model$stanfit$rawest)) {
    point <- as.numeric(model$stanfit$rawest)
    source <- "optimized estimates from fit$stanfit$rawest"
  } else return(NULL)
  if (any(!is.finite(point))) return(NULL)
  list(point = point, source = source, layout = .ctma_raw_layout(model))
}

.ctma_auto_init_template <- function(model, spec, sid, dummy_names,
    remove_study_fixed_effects, verbose) {
  src <- .ctma_raw_point(model)
  if (is.null(src)) {
    if (isTRUE(verbose)) .ctma_print_message(
      "Automatic initialization was requested, but the supplied object has no usable raw fitted estimates. Stan/ctsem random initialization will be used."
    )
    return(NULL)
  }
  d <- model$standata
  rawmean <- as.numeric(src$point[src$layout$rawpopmeans])
  old_setup <- d$TIPREDEFFECTsetup
  old_nti <- as.integer(d$ntipredeffects)
  beta <- if (old_nti) as.numeric(src$point[src$layout$tipredeffectparams]) else numeric()
  if (length(beta) != old_nti) return(NULL)

  ti_names <- model$ctstanmodelbase$TIpredNames
  dummy_cols <- match(dummy_names, ti_names)
  study_rows <- match(seq_len(sid$nstudies), sid$id)
  xstudy <- if (length(dummy_cols))
    as.matrix(d$tipredsdata[study_rows, dummy_cols, drop = FALSE]) else
    matrix(0, sid$nstudies, 0L)
  effects <- matrix(0, sid$nstudies, nrow(spec))
  means <- rawmean
  for (q in seq_len(nrow(spec))) if (spec$type[q] == 1L) {
    theta <- rep(rawmean[spec$index[q]], sid$nstudies)
    if (length(dummy_cols)) for (j in seq_along(dummy_cols)) {
      e <- old_setup[spec$index[q], dummy_cols[j]]
      if (e > 0L) theta <- theta + xstudy[, j] * beta[e]
    }
    means[spec$index[q]] <- mean(theta)
    effects[, q] <- theta - mean(theta)
  }
  sds <- apply(effects, 2L, stats::sd)
  sds[!is.finite(sds) | sds < 0.05] <- 0.05

  # Determine which old TI coefficients remain after study-indicator removal.
  setup_new <- old_setup
  if (isTRUE(remove_study_fixed_effects) && length(dummy_cols))
    setup_new[unique(spec$index[spec$type == 1L]), dummy_cols] <- 0L
  used_old <- sort(unique(as.integer(setup_new[setup_new > 0L])))
  retained_beta <- if (length(used_old)) beta[used_old] else numeric()
  ti_pos <- src$layout$tipredeffectparams
  prefix <- if (length(ti_pos)) src$point[seq_len(min(ti_pos) - 1L)] else src$point
  suffix <- if (length(ti_pos) && max(ti_pos) < length(src$point))
    src$point[seq.int(max(ti_pos) + 1L, length(src$point))] else numeric()
  prefix[seq_along(means)] <- means
  unconstrained_base <- c(prefix, retained_beta, suffix)

  if (isTRUE(verbose)) .ctma_print_message(paste0(
    "Automatic initialization constructed from ", src$source, ": ",
    length(means), " population parameters, ", length(retained_beta),
    " retained TI effects, and ", sid$nstudies, " study-effect rows."
  ))
  list(rawpopmeans = means, tipredeffectparams = retained_beta,
       studyre_effect = effects, studyre_sd = sds,
       unconstrained_base = unconstrained_base,
       setup_new = setup_new, used_old = used_old, source = src$source)
}

.ctma_init_function <- function(prepared, jitter = 0.02) {
  template <- prepared$init_template
  if (is.null(template)) return("random")
  force(prepared); force(jitter)
  function(chain_id = 1L) {
    set.seed(as.integer(104729L + chain_id))
    q <- prepared$standata$nstudyre; k <- prepared$standata$nstudies
    out <- list(
      rawpopmeans = array(template$rawpopmeans +
        stats::rnorm(length(template$rawpopmeans), 0, jitter),
        dim = length(template$rawpopmeans)),
      tipredeffectparams = array(template$tipredeffectparams +
        stats::rnorm(length(template$tipredeffectparams), 0, jitter),
        dim = length(template$tipredeffectparams)),
      studyre_sd = array(pmax(0.01, template$studyre_sd *
        exp(stats::rnorm(q, 0, jitter))), dim = q)
    )
    if (prepared$parameterization == "centered")
      out$studyre_effect <- matrix(template$studyre_effect +
        stats::rnorm(k * q, 0, jitter), k, q)
    else out$studyre_z <- t(template$studyre_effect /
      matrix(template$studyre_sd, k, q, byrow = TRUE)) +
      matrix(stats::rnorm(k * q, 0, jitter), q, k)
    blocks <- prepared$covariance_blocks
    for (i in seq_len(nrow(blocks))) {
      suffix <- paste0("b", blocks$block_id[i])
      if (blocks$structure[i] == "full")
        out[[paste0("studyre_Lcorr_", suffix)]] <- diag(blocks$size[i])
      if (blocks$structure[i] == "equicorrelated")
        out[[paste0("studyre_rho_", suffix)]] <- 0
    }
    out
  }
}

.ctma_sqrtpcov_index <- function(r, cc, n) {
  if (r <= cc) stop("A T0 correlation must be specified in the lower triangle.")
  counter <- 0L
  for (j in seq_len(n - 1L)) for (i in seq.int(j + 1L, n)) {
    counter <- counter + 1L
    if (i == r && j == cc) return(counter)
  }
  stop("Could not map the T0 correlation parameter.")
}

.ctma_re_spec <- function(model, reDRIFT, reDIFFUSION, reT0VAR) {
  nl <- as.integer(model$ctstanmodelbase$n.latent)
  masks <- list(
    DRIFT = .ctma_bool_matrix(reDRIFT, nl, "reDRIFT"),
    DIFFUSION = .ctma_bool_matrix(reDIFFUSION, nl, "reDIFFUSION"),
    T0VAR = .ctma_bool_matrix(reT0VAR, nl, "reT0VAR")
  )
  basepars <- model$ctstanmodelbase$pars
  ms <- model$setup$matsetup
  matrix_code <- c(DRIFT = 3L, DIFFUSION = 4L, T0VAR = 8L)
  niv <- as.integer(model$standata$nindvarying)
  state_index <- as.integer(model$standata$intoverpopindvaryingindex)

  ans <- list()
  q <- 0L
  for (mat in names(masks)) {
    selected <- which(masks[[mat]], arr.ind = TRUE)
    if (!nrow(selected)) next
    for (a in seq_len(nrow(selected))) {
      r <- selected[a, "row"]; cc <- selected[a, "col"]
      bp <- which(basepars$matrix == mat & basepars$row == r & basepars$col == cc)
      if (!length(bp) || is.na(basepars$param[bp[1L]]) ||
          (!is.na(basepars$value[bp[1L]]) && !is.nan(basepars$value[bp[1L]]))) {
        stop(mat, "[", r, ",", cc, "] is not a free parameter.")
      }
      q <- q + 1L
      label <- basepars$param[bp[1L]]

      if (mat != "T0VAR") {
        # popsetup maps the substantive parameter to rawpopmeans even when
        # ctsem implements an individually varying parameter as an augmented
        # latent state. Coordinate-based matsetup is only a fallback.
        ps <- model$setup$popsetup
        zps <- which(!is.na(ps$parname) & ps$parname == label & ps$param > 0L)
        z <- which(ms$matrix == matrix_code[mat] & ms$row == r & ms$col == cc &
                     ms$param > 0L)
        raw_index <- if (length(zps)) ps$param[zps[1L]] else if (length(z))
          ms$param[z[1L]] else NA_integer_
        if (is.na(raw_index))
          stop("No raw population-mean parameter index for ", mat, "[", r, ",", cc,
               "]. This may be a state-dependent expression rather than a free parameter.")
        person_varying <- isTRUE(basepars$indvarying[bp[1L]])
        ans[[q]] <- data.frame(
          re = q, matrix = mat, row = r, col = cc, parameter = label,
          type = 1L, index = raw_index, t0_row = 0L, t0_col = 0L,
          person_varying = person_varying,
          interpretation = if (person_varying)
            "study shift in the mean of a person-varying parameter" else
            "study-varying parameter without person-level variation",
          transform = basepars$transform[bp[1L]], stringsAsFactors = FALSE
        )
      } else {
        # With no individually varying initial states, a free T0VAR element is
        # an ordinary rawpopmeans parameter. If ctsem has absorbed T0VAR into
        # its individual initial-state covariance, matsetup$param is zero and
        # we instead modify rawpopsdbase/sqrtpcov below.
        mz <- which(ms$matrix == matrix_code[mat] & ms$row == r & ms$col == cc &
                      ms$param > 0L)
        if (length(mz)) {
          ps <- model$setup$popsetup
          zps <- which(!is.na(ps$parname) & ps$parname == label & ps$param > 0L)
          raw_index <- if (length(zps)) ps$param[zps[1L]] else ms$param[mz[1L]]
          ans[[q]] <- data.frame(
            re = q, matrix = mat, row = r, col = cc, parameter = label,
            type = 1L, index = raw_index, t0_row = 0L, t0_col = 0L,
            person_varying = FALSE,
            interpretation = "study-varying T0 parameter without person-level T0 variation",
            transform = basepars$transform[bp[1L]], stringsAsFactors = FALSE
          )
          next
        }
        if (cc > r) stop("Select T0VAR elements in the lower triangle only.")
        ir <- match(r, state_index); ic <- match(cc, state_index)
        if (is.na(ir) || is.na(ic))
          stop("T0VAR[", r, ",", cc,
               "] is not represented in ctsem's individually varying initial states.")
        if (ir < ic) { tmp <- ir; ir <- ic; ic <- tmp }
        if (ir == ic) {
          type <- 2L; index <- ir
        } else {
          type <- 3L; index <- .ctma_sqrtpcov_index(ir, ic, niv)
        }
        ans[[q]] <- data.frame(
          re = q, matrix = mat, row = r, col = cc, parameter = label,
          type = type, index = index, t0_row = ir, t0_col = ic,
          person_varying = TRUE,
          interpretation = if (type == 2L)
            "study-varying person-level initial-state SD" else
            "study-varying person-level initial-state correlation basis",
          transform = basepars$transform[bp[1L]], stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(ans)) stop("At least one study-level random effect must be selected.")
  do.call(rbind, ans)
}

.ctma_covariance_name <- function(x) {
  aliases <- c(
    independent = "independent", ure = "independent",
    equicorrelated = "equicorrelated", equal = "equicorrelated",
    ere = "equicorrelated", full = "full", cre = "full",
    block = "block"
  )
  if (length(x) != 1L || is.na(x)) stop(
    "A covariance structure must be one of 'independent', ",
    "'equicorrelated', 'full', or 'block'."
  )
  out <- unname(aliases[tolower(as.character(x))])
  if (is.na(out)) stop("Unknown covariance structure: ", x, ".")
  out
}

.ctma_re_blocks <- function(spec, covariance, covarianceBlocks = NULL,
    reBlocks = NULL, blockStructures = NULL, nlatent, verbose = TRUE) {
  nlatent <- as.integer(nlatent)[1L]
  if (is.na(nlatent) || nlatent < 1L) stop("nlatent must be a positive integer.")
  covariance <- .ctma_covariance_name(covariance)
  q <- nrow(spec)

  if (covariance != "block") {
    if (!is.null(covarianceBlocks) || !is.null(reBlocks) ||
        !is.null(blockStructures)) warning(
      "covarianceBlocks, reBlocks, and blockStructures are ignored unless ",
      "covariance='block'."
    )
    block_name <- "all"
    structure_name <- covariance
    spec$block <- block_name
  } else {
    if (is.null(reBlocks)) {
      spec$block <- spec$matrix
      if (is.null(covarianceBlocks)) stop(
        "With covariance='block' and automatic matrix blocks, supply ",
        "covarianceBlocks, for example list(DRIFT='equicorrelated', ",
        "DIFFUSION='independent', T0VAR='independent')."
      )
      structures <- unlist(covarianceBlocks, use.names = TRUE)
    } else {
      if (!is.list(reBlocks)) stop(
        "reBlocks must be a named list of DRIFT, DIFFUSION, and T0VAR matrices."
      )
      for (mat in intersect(names(reBlocks), c("DRIFT", "DIFFUSION", "T0VAR"))) {
        bm <- as.matrix(reBlocks[[mat]])
        if (!identical(dim(bm), c(nlatent, nlatent))) stop(
          "reBlocks$", mat, " must be a ", nlatent, " x ", nlatent, " matrix."
        )
        assigned <- !is.na(bm) & nzchar(as.character(bm))
        selected <- matrix(FALSE, nlatent, nlatent)
        z <- which(spec$matrix == mat)
        if (length(z)) selected[cbind(spec$row[z], spec$col[z])] <- TRUE
        if (any(assigned & !selected)) stop(
          "reBlocks$", mat, " assigns block names to cells that are not ",
          "selected by the corresponding random-effect mask."
        )
      }
      assignments <- character(q)
      for (i in seq_len(q)) {
        mat <- spec$matrix[i]
        bm <- reBlocks[[mat]]
        if (is.null(bm)) stop("reBlocks has no ", mat, " assignment matrix.")
        bm <- as.matrix(bm)
        if (!identical(dim(bm), c(nlatent, nlatent))) stop(
          "reBlocks$", mat, " must be a ", nlatent, " x ", nlatent, " matrix."
        )
        value <- bm[spec$row[i], spec$col[i]]
        if (length(value) != 1L || is.na(value) || !nzchar(as.character(value)))
          stop(mat, "[", spec$row[i], ",", spec$col[i],
            "] is selected as random but has no reBlocks assignment.")
        assignments[i] <- as.character(value)
      }
      spec$block <- assignments
      structures <- unlist(blockStructures, use.names = TRUE)
      if (is.null(blockStructures) || is.null(names(structures))) stop(
        "Custom reBlocks require a named blockStructures vector or list."
      )
    }

    if (is.null(names(structures)) || any(!nzchar(names(structures)))) stop(
      "Block structures must be named by block."
    )
    missing <- setdiff(unique(spec$block), names(structures))
    extra <- setdiff(names(structures), unique(spec$block))
    if (length(missing)) stop("No covariance structure supplied for block(s): ",
      paste(missing, collapse = ", "), ".")
    if (length(extra)) warning("Unused covariance block definition(s): ",
      paste(extra, collapse = ", "), ".")
  }

  block_order <- unique(spec$block)
  if (covariance != "block") structures <- stats::setNames(
    structure_name, block_name
  )
  structures <- vapply(block_order, function(b)
    .ctma_covariance_name(structures[[b]]), character(1L))
  if (any(structures == "block")) stop(
    "A block's own structure cannot be 'block'. Use independent, ",
    "equicorrelated, or full."
  )

  # Keep members of each block contiguous so Stan can use block-diagonal
  # Cholesky factors efficiently. The original matrix/row/column coordinates
  # remain in the specification and all RE indices are regenerated.
  spec <- spec[order(match(spec$block, block_order), spec$re), , drop = FALSE]
  rownames(spec) <- NULL
  spec$re <- seq_len(nrow(spec))
  sizes <- as.integer(table(factor(spec$block, levels = block_order)))
  starts <- cumsum(c(1L, head(sizes, -1L)))
  effective <- structures
  singleton <- sizes == 1L & structures != "independent"
  if (any(singleton)) {
    if (isTRUE(verbose)) .ctma_print_message(paste0(
      "Singleton covariance block(s) ",
      paste(block_order[singleton], collapse = ", "),
      " contain no correlation; their effective structure is independent."
    ))
    effective[singleton] <- "independent"
  }
  info <- data.frame(
    block_id = seq_along(block_order), block = block_order,
    start = starts, size = sizes, requested = unname(structures),
    structure = unname(effective), stringsAsFactors = FALSE
  )
  spec$block_id <- match(spec$block, info$block)

  if (isTRUE(verbose)) {
    lines <- paste0(info$block, ": ", info$size, " parameter(s), ",
      info$structure, "; ", ifelse(info$structure == "full",
        info$size * (info$size - 1L) / 2L,
        ifelse(info$structure == "equicorrelated", 1L, 0L)),
      " correlation parameter(s)")
    .ctma_print_message(paste0(
      "Study-effect covariance blocks: ", paste(lines, collapse = "; "),
      ". Correlations between blocks are fixed to zero."
    ))
  }
  list(specification = spec, blocks = info, covariance = covariance)
}

.ctma_resolve_parameterization <- function(parameterization, study_id,
    verbose = TRUE) {
  parameterization <- match.arg(parameterization,
    c("auto", "centered", "noncentered"))
  if (parameterization != "auto") return(parameterization)

  # A centred parameterisation is generally preferable when each group effect
  # is informed by many independent subjects; a non-centred parameterisation
  # is generally preferable for sparse/weakly informed groups. This observable
  # rule is deliberately simple and is reported rather than hidden.
  group_sizes <- as.integer(table(study_id))
  resolved <- if (stats::median(group_sizes) >= 20) "centered" else "noncentered"
  if (isTRUE(verbose)) .ctma_print_message(paste0(
    "parameterization='auto' selected '", resolved,
    "' (study sizes: min = ", min(group_sizes),
    ", median = ", format(stats::median(group_sizes), trim = TRUE),
    ", max = ", max(group_sizes), ")."
  ))
  resolved
}

.ctma_stan_block_code <- function(blocks) {
  correlated <- blocks$structure != "independent"
  has_corr <- any(correlated)
  par <- tp_decl <- tp_build <- prior <- character()
  ncorr <- 0L
  if (has_corr) {
    tp_decl <- c(tp_decl, "  matrix[nstudyre,nstudyre] studyre_Lcorr;")
    tp_build <- c(tp_build,
      "  studyre_Lcorr = diag_matrix(rep_vector(1.0,nstudyre));")
  }
  for (i in seq_len(nrow(blocks))) {
    size <- blocks$size[i]; start <- blocks$start[i]
    end <- start + size - 1L
    suffix <- paste0("b", blocks$block_id[i])
    if (blocks$structure[i] == "full") {
      par <- c(par, sprintf(
        "  cholesky_factor_corr[%d] studyre_Lcorr_%s;", size, suffix))
      tp_build <- c(tp_build, sprintf(
        "  studyre_Lcorr[%d:%d,%d:%d] = studyre_Lcorr_%s;",
        start, end, start, end, suffix))
      prior <- c(prior, sprintf(paste0(
        "  if(studyre_use_hyperpriors == 1)\n",
        "    target += lkj_corr_cholesky_lpdf(studyre_Lcorr_%s | studyre_lkj_eta);"),
        suffix))
      ncorr <- ncorr + size * (size - 1L) / 2L
    } else if (blocks$structure[i] == "equicorrelated") {
      lower <- -1 / (size - 1)
      par <- c(par, sprintf(
        "  real<lower=%.17g,upper=1> studyre_rho_%s;", lower, suffix))
      tp_decl <- c(tp_decl, sprintf("  matrix[%d,%d] studyre_R_%s;",
        size, size, suffix))
      tp_build <- c(tp_build,
        sprintf("  studyre_R_%s = rep_matrix(studyre_rho_%s,%d,%d);",
          suffix, suffix, size, size),
        sprintf("  for(j in 1:%d) studyre_R_%s[j,j] = 1;", size, suffix),
        sprintf(paste0("  studyre_Lcorr[%d:%d,%d:%d] = ",
          "cholesky_decompose(studyre_R_%s);"),
          start, end, start, end, suffix))
      # This is the LKJ kernel restricted to the one-dimensional
      # equicorrelation family and regularizes the common rho toward zero.
      prior <- c(prior, sprintf(paste0(
        "  if(studyre_use_hyperpriors == 1)\n",
        "    target += lkj_corr_lpdf(studyre_R_%s | studyre_lkj_eta);"),
        suffix))
      ncorr <- ncorr + 1L
    }
  }
  list(has_corr = has_corr, par = paste(par, collapse = "\n"),
    tp_decl = paste(tp_decl, collapse = "\n"),
    tp_build = paste(tp_build, collapse = "\n"),
    prior = paste(prior, collapse = "\n"), ncorr = ncorr)
}

.ctma_patch_stan <- function(text, covariance, parameterization, blocks) {
  block_code <- .ctma_stan_block_code(blocks)
  data_anchor <- "  int nindvaryingoffdiagonals; //number of off diagonal parameters needed for popcov matrix"
  data_add <- paste0(data_anchor, "\n",
    "  // CoTiMA study-level random effects\n",
    "  int<lower=1> nstudies;\n",
    "  int<lower=1> nstudyre;\n",
    "  array[nsubjects] int<lower=1,upper=nstudies> study_id;\n",
    "  array[nstudyre] int<lower=1,upper=3> studyre_type;\n",
    "  array[nstudyre] int<lower=1> studyre_index;\n",
    "  array[nstudyre] int<lower=0> studyre_row;\n",
    "  array[nstudyre] int<lower=0> studyre_col;\n",
    "  vector<lower=0>[nstudyre] studyre_sd_prior;\n",
    "  real<lower=1> studyre_lkj_eta;\n",
    "  int<lower=0,upper=1> studyre_use_hyperpriors;")
  text <- .ctma_replace_once(text, data_anchor, data_add, "data block")

  par_anchor <- "  vector[(nsubsets > 1) ? 1 : 0] subsetpar;"
  par_add <- paste0(par_anchor, "\n",
    "  // CoTiMA study-effect distribution parameters\n",
    if (parameterization == "noncentered")
      "  matrix[nstudyre,nstudies] studyre_z;\n" else
      "  matrix[nstudies,nstudyre] studyre_effect;\n",
    "  vector<lower=0>[nstudyre] studyre_sd;\n",
    block_code$par)
  text <- .ctma_replace_once(text, par_anchor, par_add, "parameters block")

  tp_anchor <- "  matrix[nindvarying, nindvarying] rawpopcorr;"
  effect_declaration <- if (parameterization == "noncentered")
    "  matrix[nstudies,nstudyre] studyre_effect;\n" else ""
  effect_formula <- if (parameterization == "centered") "" else
    if (block_code$has_corr)
      "  studyre_effect = (diag_pre_multiply(studyre_sd, studyre_Lcorr) * studyre_z)';"
    else
      "  studyre_effect = (diag_matrix(studyre_sd) * studyre_z)';"
  tp_add <- paste0(tp_anchor, "\n",
    effect_declaration,
    if (nzchar(block_code$tp_decl)) paste0(block_code$tp_decl, "\n") else "",
    "  array[nstudies] matrix[nindvarying,nindvarying] studyre_rawpopcovbase;\n",
    "  array[nstudies] matrix[nindvarying,nindvarying] studyre_T0VAR;\n",
    "  array[nstudies] matrix[nindvarying,nindvarying] studyre_T0cov;\n",
    if (nzchar(block_code$tp_build)) paste0(block_code$tp_build, "\n") else "",
    effect_formula)
  text <- .ctma_replace_once(text, tp_anchor, tp_add, "transformed parameters declarations")

  cov_anchor <- "  }//end indvarying par setup"
  cov_add <- paste0(cov_anchor, "\n",
    "\n  // Construct study-specific raw T0 SD/correlation matrices.\n",
    "  for(s in 1:nstudies){\n",
    "    studyre_rawpopcovbase[s] = rawpopcovbase;\n",
    "    for(q in 1:nstudyre){\n",
    "      if(studyre_type[q] == 2){\n",
    "        int j = studyre_index[q];\n",
    "        studyre_rawpopcovbase[s][j,j] = log1p_exp(2*(rawpopsdbase[j] + studyre_effect[s,q])-1) * sdscale[j] + 1e-10;\n",
    "      }\n",
    "      if(studyre_type[q] == 3){\n",
    "        int r = studyre_row[q];\n",
    "        int c = studyre_col[q];\n",
    "        studyre_rawpopcovbase[s][r,c] = inv_logit(sqrtpcov[studyre_index[q]] + studyre_effect[s,q])*2-1;\n",
    "        studyre_rawpopcovbase[s][c,r] = 0;\n",
    "      }\n",
    "    }\n",
    "    studyre_T0VAR[s] = studyre_rawpopcovbase[s];\n",
    "    studyre_T0cov[s] = sdcovsqrt2cov(studyre_T0VAR[s],choleskymats);\n",
    "    // Apply the same T0-mean scaling used by ctsem.\n",
    "    if(nindvarying > 0){\n",
    "      for(q in 1:nindvarying){\n",
    "        for(ri in 1:size(matsetup)){\n",
    "          if(matsetup[ri,7]==1 && matsetup[ri,5] && matsetup[ri,1]==intoverpopindvaryingindex[q]){\n",
    "            studyre_T0cov[s][q,] *= matvalues[ri,2] * matvalues[ri,3];\n",
    "            studyre_T0cov[s][,q] *= matvalues[ri,2] * matvalues[ri,3];\n",
    "          }\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "  }")
  text <- .ctma_replace_once(text, cov_anchor, cov_add, "T0 covariance construction")

  raw_anchor <- "  if(si > 0 &&  ntieffects > 0){"
  raw_add <- paste0(
    "  if(si > 0){\n",
    "    for(q in 1:nstudyre) if(studyre_type[q] == 1)\n",
    "      rawindparams[studyre_index[q]] += studyre_effect[study_id[si],q];\n",
    "  }\n\n", raw_anchor)
  text <- .ctma_replace_once(text, raw_anchor, raw_add, "subject raw parameters")

  t0_anchor <- "if(intoverpop && nindvarying > 0) T0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovbase;"
  t0_add <- paste0(
    "    if(intoverpop && nindvarying > 0){\n",
    "      if(si > 0) T0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = studyre_rawpopcovbase[study_id[si]];\n",
    "      else T0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovbase;\n",
    "    }")
  text <- .ctma_replace_once(text, t0_anchor, t0_add, "subject T0 covariance")

  model_anchor <- "  real priormod2 = priormod / nsubsets;"
  study_distribution <- if (parameterization == "noncentered")
    "  target += std_normal_lpdf(to_vector(studyre_z));\n" else if (!block_code$has_corr)
    paste0("  for(q in 1:nstudyre)\n",
      "    target += normal_lpdf(studyre_effect[,q] | 0, studyre_sd[q]);\n") else
    paste0("  for(s in 1:nstudies)\n",
      "    target += multi_normal_cholesky_lpdf(studyre_effect[s]' | rep_vector(0,nstudyre), diag_pre_multiply(studyre_sd, studyre_Lcorr));\n")
  model_add <- paste0(model_anchor, "\n",
    "  // This term defines the study-effect distribution and\n",
    "  // is retained for Bayes, MAP, and joint/conditional ML. Hyperpriors\n",
    "  // on its scale and correlation can be disabled for ML optimization.\n",
    study_distribution,
    "  if(studyre_use_hyperpriors == 1)\n",
    "    target += normal_lpdf(studyre_sd | 0, studyre_sd_prior);\n",
    block_code$prior)
  .ctma_replace_once(text, model_anchor, model_add, "model priors")
}

#' Prepare a hierarchical study-level random-effects CoTiMA model
#'
#' @description
#' Adds study-level random effects to selected ctsem DRIFT, DIFFUSION, and
#' T0VAR parameters. Effects are inserted on ctsem's unconstrained parameter
#' scale before parameter-specific transformations are applied.
#'
#' @param model A pooled fitted ctsem model extracted from a CoTiMA full fit.
#' @param study_id Optional study identifier at the subject or observation row
#'   level.
#' @param study_id_ti Optional name or index of a time-independent predictor
#'   containing study membership.
#' @param study_dummy_ti Optional names or indices of the K - 1 study
#'   indicator predictors. Ordinary, scaled, and weighted-effect coding are
#'   supported.
#' @param reDRIFT Logical \code{n.latent} by \code{n.latent} matrix selecting
#'   study-varying DRIFT cells. If \code{NULL}, the stored CoTiMA setting is
#'   used.
#' @param reDIFFUSION Logical matrix selecting free lower-triangular
#'   DIFFUSION cells.
#' @param reT0VAR Logical matrix selecting free lower-triangular T0VAR cells.
#' @param covariance Study-effect correlation structure: \code{"independent"},
#'   \code{"equicorrelated"}, \code{"full"}, or \code{"block"}. Aliases
#'   \code{"ure"}, \code{"ere"}, \code{"cre"}, and \code{"equal"} are
#'   accepted.
#' @param covarianceBlocks Named list assigning an independent,
#'   equicorrelated, or full structure to automatic DRIFT, DIFFUSION, and
#'   T0VAR blocks when \code{covariance = "block"}.
#' @param reBlocks Optional named list of character matrices assigning selected
#'   cells to custom blocks.
#' @param blockStructures Named vector or list assigning a covariance structure
#'   to each custom block in \code{reBlocks}.
#' @param parameterization \code{"centered"}, \code{"noncentered"}, or
#'   \code{"auto"}.
#' @param init Initialization policy. \code{"auto"} reconstructs starting
#'   values from the pooled fit; \code{"random"} delegates to Stan. Expert
#'   users may supply a list or function.
#' @param init_jitter Non-negative raw-scale initialization jitter.
#' @param remove_study_fixed_effects Remove study-indicator effects for
#'   parameters selected as hierarchical, while retaining substantive
#'   moderator effects.
#' @param sd_prior Positive half-normal prior scale for study-effect standard
#'   deviations; scalar or one value per selected parameter.
#' @param lkj_eta LKJ shape used for full blocks and as a regularizing kernel
#'   for equicorrelated blocks.
#' @param verbose Print CoTiMA-style preparation messages.
#'
#' @return An object of class \code{"CoTiMARePrep"} and
#'   \code{"ctmaStudyREModel"}, containing patched Stan code, Stan data,
#'   parameter mappings, study labels, covariance blocks, and initialization.
#' @export
#' @seealso \code{\link{ctmaReFit}}, \code{\link{ctmaReMLFit}}
#'
#' @examples
#' \dontrun{
#' prepared <- ctmaRePrep(
#'   pooledFit,
#'   study_dummy_ti = paste0("TI", 1:9),
#'   reDRIFT = matrix(TRUE, 2, 2),
#'   reDIFFUSION = matrix(FALSE, 2, 2),
#'   reT0VAR = matrix(FALSE, 2, 2),
#'   covariance = "block",
#'   covarianceBlocks = list(DRIFT = "equicorrelated")
#' )
#' }
ctmaRePrep <- function(
    model,
    study_id = NULL,
    study_id_ti = NULL,
    study_dummy_ti = NULL,
    reDRIFT = NULL,
    reDIFFUSION = NULL,
    reT0VAR = NULL,
    covariance = c("independent", "equicorrelated", "full", "block"),
    covarianceBlocks = NULL,
    reBlocks = NULL,
    blockStructures = NULL,
    parameterization = c("auto", "centered", "noncentered"),
    init = c("auto", "random"),
    init_jitter = 0.02,
    remove_study_fixed_effects = TRUE,
    sd_prior = 0.25,
    lkj_eta = 2,
    verbose = TRUE) {

  covariance <- .ctma_covariance_name(covariance[1L])
  parameterization_requested <- match.arg(parameterization)
  if (is.character(init)) init <- match.arg(init)
  if (!is.numeric(init_jitter) || length(init_jitter) != 1L ||
      !is.finite(init_jitter) || init_jitter < 0)
    stop("init_jitter must be one finite non-negative number.")
  if (!is.list(model) || is.null(model$stanmodeltext) || is.null(model$standata))
    stop("model must be the object returned by CoTiMA/ctsem with fit = FALSE.")

  stored <- .ctma_stored_re(model)
  if (is.null(reDRIFT)) reDRIFT <- stored$reDRIFT
  if (is.null(reDIFFUSION)) reDIFFUSION <- stored$reDIFFUSION
  if (is.null(reT0VAR)) reT0VAR <- stored$reT0VAR
  if (is.null(study_id) && is.null(study_id_ti) &&
      "study" %in% model$ctstanmodelbase$TIpredNames)
    study_id_ti <- "study"
  dummy_names <- character()
  if (is.null(study_id) && is.null(study_id_ti)) {
    inferred <- .ctma_infer_study_id(
      model, study_dummy_ti = study_dummy_ti, verbose = verbose
    )
    study_id <- inferred$id
    dummy_names <- inferred$dummy_names
  } else if (!is.null(study_dummy_ti)) {
    dummy_names <- if (is.numeric(study_dummy_ti))
      model$ctstanmodelbase$TIpredNames[as.integer(study_dummy_ti)] else
      as.character(study_dummy_ti)
  }

  # The identifier may be transported as a TI column, but it must not act as
  # a numeric moderator in the original ctsem model.
  if (!is.null(study_id_ti) && is.character(study_id_ti)) {
    effect_col <- paste0(study_id_ti, "_effect")
    if (effect_col %in% names(model$ctstanmodelbase$pars) &&
        any(model$ctstanmodelbase$pars[[effect_col]], na.rm = TRUE))
      stop("The study identifier has ordinary TI effects enabled. Set all ",
           effect_col, " entries to FALSE; it is an identifier, not a moderator.")
  }
  sid <- .ctma_study_id(model, study_id, study_id_ti)
  parameterization <- .ctma_resolve_parameterization(
    parameterization_requested, sid$id, verbose = verbose
  )
  spec <- .ctma_re_spec(model, reDRIFT, reDIFFUSION, reT0VAR)
  block_result <- .ctma_re_blocks(
    spec = spec, covariance = covariance,
    covarianceBlocks = covarianceBlocks, reBlocks = reBlocks,
    blockStructures = blockStructures,
    nlatent = as.integer(model$ctstanmodelbase$n.latent), verbose = verbose
  )
  spec <- block_result$specification
  covariance_blocks <- block_result$blocks
  q <- nrow(spec)

  init_template <- if (is.character(init) && identical(init, "auto"))
    .ctma_auto_init_template(model, spec, sid, dummy_names,
      remove_study_fixed_effects, verbose) else NULL

  # Remove only the study-indicator effects for parameters selected as random;
  # later TI columns representing substantive moderators remain untouched.
  if (isTRUE(remove_study_fixed_effects) && length(dummy_names)) {
    dummy_cols <- match(dummy_names, model$ctstanmodelbase$TIpredNames)
    selected_raw <- unique(spec$index[spec$type == 1L])
    setup <- model$standata$TIPREDEFFECTsetup
    setup[selected_raw, dummy_cols] <- 0L
    used <- sort(unique(as.integer(setup[setup > 0L])))
    active <- setup > 0L
    if (any(active)) setup[active] <- match(setup[active], used)
    model$standata$TIPREDEFFECTsetup <- setup
    model$standata$ntipredeffects <- length(used)
    if (!is.null(model$data$TIPREDEFFECTsetup)) {
      model$data$TIPREDEFFECTsetup <- setup
      model$data$ntipredeffects <- length(used)
    }
    effect_cols <- paste0(dummy_names, "_effect")
    parameter_rows <- model$ctstanmodelbase$pars$param %in%
      spec$parameter[spec$type == 1L]
    model$ctstanmodelbase$pars[parameter_rows, effect_cols] <- FALSE
    if (isTRUE(verbose)) .ctma_print_message(paste0(
      "Removed fixed study-indicator effects from ", length(selected_raw),
      " selected raw parameter row(s); retained substantive moderator effects remain active."
    ))
  }

  mean_and_person <- spec$type == 1L & spec$person_varying
  if (any(mean_and_person)) warning(
    "The following parameters vary at both levels: ",
    paste(spec$parameter[mean_and_person], collapse = ", "),
    ". The study effect shifts the study-specific parameter mean; ctsem's existing ",
    "individual random effect describes person deviations around that study mean."
  )
  t0_person_cov <- spec$type %in% c(2L, 3L)
  if (any(t0_person_cov) && as.integer(model$standata$intoverpop) != 1L)
    stop(
      "Study-varying T0 covariance parameters currently require intoverpop = TRUE. ",
      "With intoverpop = FALSE, ctsem generates person deviations using one common ",
      "rawpopcovchol before the subject loop, so the present patch cannot safely apply ",
      "a study-specific covariance."
    )
  if (any(t0_person_cov)) warning(
    "Selected T0VAR elements are part of ctsem's person-level initial-state ",
    "covariance. Their study effects therefore model study differences in the ",
    "within-study distribution of individual initial states, not shifts in the ",
    "study mean initial state."
  )
  full_blocks <- covariance_blocks$structure == "full"
  if (any(full_blocks & sid$nstudies <= covariance_blocks$size + 1L)) warning(
    "At least one full study-effect correlation block is large relative to ",
    "the ", sid$nstudies, " studies. Consider equicorrelated or independent ",
    "blocks, fewer random parameters, or stronger priors."
  )
  if (length(sd_prior) == 1L) sd_prior <- rep(sd_prior, q)
  if (length(sd_prior) != q || any(!is.finite(sd_prior)) || any(sd_prior <= 0))
    stop("sd_prior must be positive and have length one or nstudyre.")

  standata <- model$standata
  standata$nstudies <- sid$nstudies
  standata$nstudyre <- q
  standata$study_id <- sid$id
  standata$studyre_type <- .ctma_stan_1d(spec$type, "integer")
  standata$studyre_index <- .ctma_stan_1d(spec$index, "integer")
  standata$studyre_row <- .ctma_stan_1d(spec$t0_row, "integer")
  standata$studyre_col <- .ctma_stan_1d(spec$t0_col, "integer")
  standata$studyre_sd_prior <- .ctma_stan_1d(sd_prior, "double")
  standata$studyre_lkj_eta <- as.numeric(lkj_eta)
  standata$studyre_use_hyperpriors <- 1L

  # ctsem's generated Stan program uses whenvecp/whenmat to decide which raw
  # parameters and system matrices must be recomputed for every subject. A
  # parameter that was invariant in the source FI model is otherwise computed
  # only for si == 0, so adding a study shift to rawindparams has no effect on
  # the likelihood. Mark every selected study-RE parameter as subject-varying
  # in the scheduling data. DRIFT also requires its analytic Jacobian (JAx) to
  # be refreshed because ctsem may use JAxDRIFTequiv for the transition matrix.
  matrix_code <- c(DRIFT = 3L, DIFFUSION = 4L, T0VAR = 8L)
  raw_re <- which(spec$type == 1L)
  if (length(raw_re)) {
    # ctsem calls whichequals(whenvecp[2, ], 0, 0), where comparison=0
    # means "not equal". Consequently, a zero DISABLES subject-level
    # transformation; an active entry contains its parameter index. Earlier
    # versions of this script set zero here and therefore silently left the
    # population parameter in every subject's likelihood.
    raw_indices <- unique(spec$index[raw_re])
    standata$whenvecp[2L, raw_indices] <- raw_indices
    # ctsem's parvectform() and mcalc() contain a second subject-level gate:
    # the relevant matsetup row must be marked as individually varying,
    # TI-varying, or state dependent. Merely changing whenvecp/whenmat is not
    # sufficient for a parameter that was invariant in the source model.
    # Column 6 is ctsem's TI-varying scheduling flag. Using it here activates
    # subject-level transformation and matrix construction without falsely
    # adding the parameter to ctsem's person-level random-effect covariance.
    matsetup_rows <- which(standata$matsetup[, 3L] %in% raw_indices)
    if (!length(matsetup_rows)) stop(
      "Could not activate subject-level matrix construction for the selected ",
      "study random effect(s): no matching matsetup rows were found."
    )
    standata$matsetup[matsetup_rows, 6L] <- 1L
    for (mat in unique(spec$matrix[raw_re]))
      standata$whenmat[matrix_code[[mat]], 5L] <- 1L
  }
  if (any(spec$matrix == "DRIFT")) {
    standata$whenmat[matrix_code[["DRIFT"]], 5L] <- 1L
    if (nrow(standata$whenmat) >= 52L)
      standata$whenmat[52L, 5L] <- 1L
  }
  if (any(spec$matrix == "DIFFUSION"))
    standata$whenmat[matrix_code[["DIFFUSION"]], 5L] <- 1L
  if (any(spec$matrix == "T0VAR"))
    standata$whenmat[matrix_code[["T0VAR"]], 5L] <- 1L

  # The Stan program has been patched and is therefore not the cached ctsem
  # model represented by the inherited recompile flag. Retain this marker so
  # the prepared object cannot be mistaken for ctsem's original cached model.
  standata$recompile <- 1L

  structure(list(
    original = model,
    stanmodeltext = .ctma_patch_stan(
      model$stanmodeltext, covariance, parameterization, covariance_blocks
    ),
    standata = standata,
    specification = spec,
    study_labels = sid$labels,
    covariance = covariance,
    covariance_blocks = covariance_blocks,
    parameterization = parameterization,
    parameterization_requested = parameterization_requested,
    init = init,
    init_jitter = init_jitter,
    init_template = init_template,
    study_dummy_names = dummy_names,
    remove_study_fixed_effects = remove_study_fixed_effects,
    priors = list(sd_prior = sd_prior, lkj_eta = lkj_eta)
  ), class = c("CoTiMARePrep", "ctmaStudyREModel"))
}

#' Fit a hierarchical study-level random-effects CoTiMA model
#'
#' @description
#' Fits an object returned by \code{\link{ctmaRePrep}} using Bayesian sampling, maximum
#' a posteriori optimization, or joint/conditional maximum likelihood. MAP and
#' conditional ML use the exported \code{rstan::optimizing()} interface.
#'
#' @param prepared Object returned by \code{\link{ctmaRePrep}}.
#' @param estimation \code{"bayes"}, \code{"map"}, or \code{"ml"}. The ML
#'   option optimizes study effects jointly and is not marginal ML.
#' @param method Deprecated compatibility alias: \code{"sample"} selects Bayes
#'   and \code{"optimize"} selects MAP.
#' @param cores Number of parallel Bayesian chains. A single optimizer run uses
#'   one core; use \code{\link{ctmaReMLFit}} for parallel ML starts.
#' @param chains Number of Bayesian chains.
#' @param iter Total sampling iterations per chain or maximum optimizer
#'   iterations.
#' @param seed Random-number seed.
#' @param init Optional fitting-stage initialization override.
#' @param init_jitter Optional override of the preparation-stage jitter.
#' @param verbose Print CoTiMA-style progress messages.
#' @param refresh Backend progress refresh interval.
#' @param optimizer_uncertainty Optimization uncertainty method. Currently only
#'   \code{"none"} is supported; use Bayesian estimation for intervals.
#' @param compile_args Named arguments passed to \code{rstan::stan_model()}.
#' @param fit_args Named arguments passed to \code{rstan::sampling()} or
#'   \code{rstan::optimizing()}.
#'
#' @return A \code{"CoTiMAFit"} object with subtype
#'   \code{"CoTiMAStudyREFit"}. Apply \code{summary()} for population,
#'   heterogeneity, correlation, and study-specific estimates.
#' @export
#' @seealso \code{\link{ctmaRePrep}}, \code{\link{ctmaReMLFit}}
#'
#' @examples
#' \dontrun{
#' bayesFit <- ctmaReFit(prepared, estimation = "bayes", chains = 4,
#'   iter = 2000, cores = 4)
#' mapFit <- ctmaReFit(prepared, estimation = "map", cores = 1)
#' }
ctmaReFit <- function(
    prepared,
    estimation = c("bayes", "map", "ml"),
    method = NULL,
    cores = 2,
    chains = 4,
    iter = 2000,
    seed = 1234,
    init = NULL,
    init_jitter = NULL,
    verbose = TRUE,
    refresh = 100,
    optimizer_uncertainty = c("none", "hessian", "surrogate", "is",
      "bootstrap", "fullbootstrap", "sandwich", "opg"),
    compile_args = list(),
    fit_args = list()) {

  # `method` is retained as a backward-compatible alias. Previously,
  # method="sample" meant Bayesian sampling and method="optimize" meant MAP.
  if (!is.null(method)) {
    method <- match.arg(method, c("optimize", "sample"))
    estimation <- if (method == "sample") "bayes" else "map"
  } else estimation <- match.arg(estimation)
  optimizer_uncertainty <- match.arg(optimizer_uncertainty)
  if (is.null(init)) init <- prepared$init
  if (is.null(init_jitter)) init_jitter <- prepared$init_jitter
  resolved_init <- if (is.character(init) && identical(init, "auto"))
    .ctma_init_function(prepared, init_jitter) else init
  if (estimation != "bayes" && optimizer_uncertainty != "none")
    stop(
      "optimizer_uncertainty = '", optimizer_uncertainty,
      "' is not currently supported for MAP/ML study-RE fits. ctsem's ",
      "post-optimization uncertainty and transformed-sample machinery either ",
      "subsets data without subsetting study_id or cannot restore the ",
      "dynamically compiled model, leading to dimension, invalid-model, or ",
      "'No admissible samples' errors. Use optimizer_uncertainty = 'none' for ",
      "point estimates, or estimation = 'bayes' for uncertainty intervals."
    )
  if (!inherits(prepared, "ctmaStudyREModel"))
    stop("prepared must be returned by ctmaRePrep().")
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Package 'rstan' is required.")

  started <- Sys.time()
  if (isTRUE(verbose)) .ctma_print_message(paste0(
    "Using CoTiMA_RE_hierarchical.R version ",
    attr(ctmaReFit, "script_version"), "; estimation = ", estimation,
    "; parameterization = ", prepared$parameterization,
    if (estimation == "bayes") "." else
      paste0("; optimizer_uncertainty = ", optimizer_uncertainty, ".")
  ))
  .ctma_progress(verbose, 1L, 4L,
    paste0("Validating and compiling the Stan model (", estimation, ")."), started)

  compile_call <- c(list(
    model_code = prepared$stanmodeltext,
    model_name = paste0("cotima_study_re_", prepared$covariance, "_",
      prepared$parameterization)
  ), compile_args)
  sm <- do.call(rstan::stan_model, compile_call)

  .ctma_progress(verbose, 2L, 4L,
    if (estimation == "bayes")
      paste0("Sampling ", chains, " chain(s), ", iter,
        " iterations per chain. Stan progress follows.")
    else paste0("Running ", toupper(estimation),
      " optimization with rstan. Optimizer progress follows."), started)

  fit_data <- prepared$standata
  fit_data$studyre_use_hyperpriors <- as.integer(estimation != "ml")
  custom_1d <- c("studyre_type", "studyre_index", "studyre_row",
                 "studyre_col", "studyre_sd_prior")
  bad_shape <- custom_1d[vapply(fit_data[custom_1d], function(x)
    !identical(dim(x), as.integer(length(x))), logical(1L))]
  if (length(bad_shape)) stop(
    "The following Stan data fields lost their one-dimensional array shape: ",
    paste(bad_shape, collapse = ", "),
    ". Recreate prepared with ctmaRePrep() from this script version."
  )
  if (estimation != "bayes" && cores > 1L && isTRUE(verbose)) .ctma_print_message(paste0(
    "A single rstan optimization uses one core. The cores argument remains ",
    "active for Bayesian chains; use ctmaReMLFit() for parallel independent ",
    "ML starts."
  ))

  if (estimation == "bayes") {
    if (!is.null(fit_args$init)) resolved_init <- fit_args$init
    call <- utils::modifyList(
      list(object = sm, data = fit_data, chains = chains,
           iter = iter, cores = cores, seed = seed, refresh = refresh,
           init = resolved_init),
      fit_args
    )
    engine_fit <- do.call(rstan::sampling, call)
    if (!methods::is(engine_fit, "stanfit") ||
        length(methods::slot(engine_fit, "sim")$samples) == 0L) stop(
      "Stan sampling did not create any samples. Review the Stan messages above; ",
      "no ctmaStudyREFit object was created.", call. = FALSE
    )
    .ctma_progress(verbose, 3L, 4L,
      "Sampling finished; extracting posterior draws.", started)
    draws <- rstan::extract(engine_fit, permuted = TRUE)
  } else {
    if (estimation == "ml") warning(
      "estimation='ml' maximizes the joint/conditional likelihood, retaining ",
      "the Gaussian study-effect density but omitting ctsem priors and study ",
      "hyperpriors. It does not integrate study effects and is therefore not ",
      "marginal ML; variance estimates can be biased or degenerate."
    )
    fit_data <- .ctma_optimizer_data(fit_data, estimation)

    # Use rstan's exported optimizer directly. This avoids relying on ctsem's
    # unexported stanoptimis() implementation and keeps the dynamically patched
    # Stan model and its custom study data together.
    if (!is.null(fit_args$init)) resolved_init <- fit_args$init
    forbidden <- intersect(names(fit_args), c(
      "object", "data", "seed", "as_vector", "hessian", "verbose"
    ))
    if (length(forbidden)) stop(
      "These fit_args are managed by ctmaReFit() for MAP/ML and may not be ",
      "supplied: ", paste(forbidden, collapse = ", "), "."
    )
    defaults <- list(
      object = sm, data = fit_data, seed = as.integer(seed),
      init = resolved_init, algorithm = "LBFGS", iter = as.integer(iter),
      refresh = as.integer(refresh), as_vector = FALSE, hessian = FALSE,
      verbose = isTRUE(verbose)
    )
    optim_args <- utils::modifyList(defaults, fit_args)
    engine_fit <- do.call(rstan::optimizing, optim_args)
    if (!is.list(engine_fit) || !is.list(engine_fit$par) ||
        length(engine_fit$value) != 1L || !is.finite(engine_fit$value) ||
        (!is.null(engine_fit$return_code) &&
         as.integer(engine_fit$return_code) != 0L)) stop(
      "Stan optimization returned no usable finite solution. Review the ",
      "optimizer output and try different initialization, more iterations, ",
      "or parallel multistart fitting with ctmaReMLFit().",
      call. = FALSE
    )
    .ctma_progress(verbose, 3L, 4L,
      "Optimization finished; collecting constrained point estimates.",
      started)
    # as_vector = FALSE returns all constrained parameters, transformed
    # parameters, and generated quantities as a named list. Add the one-draw
    # dimension expected by the common CoTiMA-RE summary code.
    draws <- .ctma_add_draw_dimension(engine_fit$par)
  }

  out <- structure(list(prepared = prepared, stanmodel = sm,
                 estimation = estimation,
                 optimizer_uncertainty = if (estimation == "bayes")
                   "posterior" else optimizer_uncertainty,
                 method = if (estimation == "bayes") "sample" else "optimize",
                 optimizer = if (estimation == "bayes")
                   "rstan::sampling" else "rstan::optimizing",
                 objective = if (estimation == "bayes") NULL else
                   engine_fit$value,
                 prior_switches = if (estimation == "bayes") NULL else list(
                   ctsem_priors = fit_data$priors,
                   study_hyperpriors = fit_data$studyre_use_hyperpriors
                 ),
                 engine_fit = engine_fit, draws = draws),
            class = c("CoTiMAFit", "CoTiMAStudyREFit", "ctmaStudyREFit"))
  .ctma_progress(verbose, 4L, 4L, "Fit object is ready.", started)
  out
}

# Diagnostic marker: useful for detecting an older function definition that
# is still present in an interactive R session after this file was updated.
attr(ctmaReFit, "script_version") <- "2026-08-11.20"

.ctma_draw_summary <- function(x, probs = c(.025, .5, .975)) {
  c(mean = mean(x), stats::quantile(x, probs, names = FALSE))
}

.ctma_transform <- function(param, transform) {
  if (is.na(transform) || !nzchar(transform) || transform == "param") return(param)
  softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  eval(parse(text = transform), envir = list(
    param = param, log1p_exp = softplus, exp = exp, log = log,
    log1p = log1p, sqrt = sqrt
  ))
}

# Model-implied SD of g(alpha + u), where u ~ N(0, tau^2).  Gauss-Hermite
# quadrature avoids the common but incorrect shortcut g(tau), which is not a
# standard deviation when g() is nonlinear.  The calculation is performed for
# every posterior draw and is therefore subsequently summarized like any other
# derived posterior quantity.
.ctma_transformed_re_sd <- function(alpha, tau, transform, nodes = 31L) {
  alpha <- as.numeric(alpha)
  tau <- as.numeric(tau)
  if (length(alpha) != length(tau)) stop(
    "alpha and tau must have equal lengths.", call. = FALSE
  )
  if (!length(alpha)) return(numeric())

  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 3L) stop(
    "nodes must be one integer of at least 3.", call. = FALSE
  )

  # Golub-Welsch construction for physicists' Gauss-Hermite quadrature.
  jacobi <- matrix(0, nodes, nodes)
  off_diagonal <- sqrt(seq_len(nodes - 1L) / 2)
  jacobi[cbind(seq_len(nodes - 1L), 2:nodes)] <- off_diagonal
  jacobi[cbind(2:nodes, seq_len(nodes - 1L))] <- off_diagonal
  eig <- eigen(jacobi, symmetric = TRUE)
  order_nodes <- order(eig$values)
  normal_nodes <- sqrt(2) * eig$values[order_nodes]
  normal_weights <- eig$vectors[1L, order_nodes]^2

  raw_values <- outer(tau, normal_nodes, `*`) + alpha
  transformed <- .ctma_transform(raw_values, transform)
  transformed <- matrix(transformed, nrow = length(alpha), ncol = nodes)
  first_moment <- as.vector(transformed %*% normal_weights)
  second_moment <- as.vector((transformed^2) %*% normal_weights)
  sqrt(pmax(0, second_moment - first_moment^2))
}

.ctmaReSummary <- function(fit, probs = c(.025, .5, .975)) {
  if (!inherits(fit, "ctmaStudyREFit")) stop("fit must be a ctmaStudyREFit.")
  if (identical(fit$optimizer_uncertainty, "none")) .ctma_print_message(paste0(
    "This MAP/ML fit contains point estimates only. The lower, median, and ",
    "upper summary columns therefore equal the optimized estimate."
  ))
  d <- fit$draws
  spec <- fit$prepared$specification
  studies <- fit$prepared$study_labels
  effect <- d$studyre_effect
  if (is.null(effect)) stop("studyre_effect was not returned by the Stan fit.")

  effect_summary <- do.call(rbind, lapply(seq_along(studies), function(s)
    do.call(rbind, lapply(seq_len(nrow(spec)), function(q) {
      z <- .ctma_draw_summary(effect[, s, q], probs)
      data.frame(study = studies[s], re = q, matrix = spec$matrix[q],
                 row = spec$row[q], col = spec$col[q],
                 parameter = spec$parameter[q], block = spec$block[q],
                 scale = "raw deviation",
                 mean = z[1], lower = z[2], median = z[3], upper = z[4])
    }))))

  sd_summary <- do.call(rbind, lapply(seq_len(nrow(spec)), function(q) {
    z <- .ctma_draw_summary(d$studyre_sd[, q], probs)
    data.frame(re = q, matrix = spec$matrix[q], row = spec$row[q],
               col = spec$col[q], parameter = spec$parameter[q],
               block = spec$block[q],
               mean = z[1], lower = z[2], median = z[3], upper = z[4])
  }))

  # For scalar ctsem parameters, convert the complete N(alpha, tau^2)
  # distribution to the final parameter scale.  T0VAR types 2 and 3 are
  # excluded because their final variances/covariances depend jointly on the
  # complete Cholesky matrix rather than on a scalar transformation.
  transformed_sd_summary <- data.frame()
  transformed_re <- which(spec$type == 1L)
  if (length(transformed_re)) {
    transformed_sd_summary <- do.call(rbind, lapply(transformed_re, function(q) {
      transformed_sd_draws <- .ctma_transformed_re_sd(
        alpha = d$rawpopmeans[, spec$index[q]],
        tau = d$studyre_sd[, q],
        transform = spec$transform[q]
      )
      z <- .ctma_draw_summary(transformed_sd_draws, probs)
      data.frame(re = q, matrix = spec$matrix[q], row = spec$row[q],
                 col = spec$col[q], parameter = spec$parameter[q],
                 block = spec$block[q],
                 scale = "transformed parameter distribution",
                 mean = z[1], lower = z[2], median = z[3], upper = z[4])
    }))
  }

  rawmean_re <- which(spec$type == 1L)
  study_parameter_summary <- population_parameter_summary <- data.frame()
  if (length(rawmean_re)) {
    study_parameter_summary <- do.call(rbind, lapply(seq_along(studies), function(s)
      do.call(rbind, lapply(rawmean_re, function(q) {
        raw <- d$rawpopmeans[, spec$index[q]] + effect[, s, q]
        value <- .ctma_transform(raw, spec$transform[q])
        z <- .ctma_draw_summary(value, probs)
        data.frame(study = studies[s], re = q, matrix = spec$matrix[q],
                   row = spec$row[q], col = spec$col[q],
                   parameter = spec$parameter[q], block = spec$block[q],
                   scale = "transformed absolute",
                   mean = z[1], lower = z[2], median = z[3], upper = z[4])
      }))))
    population_parameter_summary <- do.call(rbind, lapply(rawmean_re, function(q) {
      value <- .ctma_transform(d$rawpopmeans[, spec$index[q]], spec$transform[q])
      z <- .ctma_draw_summary(value, probs)
      data.frame(re = q, matrix = spec$matrix[q], row = spec$row[q],
                 col = spec$col[q], parameter = spec$parameter[q],
                 block = spec$block[q],
                 scale = "transformed population mean",
                 mean = z[1], lower = z[2], median = z[3], upper = z[4])
    }))
  }

  correlation_summary <- data.frame()
  if (any(fit$prepared$covariance_blocks$structure != "independent")) {
    nd <- dim(d$studyre_Lcorr)[1L]
    cor_draws <- array(NA_real_, c(nd, nrow(spec), nrow(spec)))
    for (i in seq_len(nd)) cor_draws[i,,] <- tcrossprod(d$studyre_Lcorr[i,,])
    pairs <- do.call(rbind, lapply(seq_len(nrow(spec)), function(r) {
      cc <- seq_len(r - 1L)
      if (!length(cc)) return(NULL)
      cc <- cc[spec$block_id[cc] == spec$block_id[r]]
      if (!length(cc)) return(NULL)
      cbind(r = r, cc = cc)
    }))
    if (!is.null(pairs) && nrow(pairs)) correlation_summary <- do.call(rbind,
      lapply(seq_len(nrow(pairs)), function(i) {
        r <- pairs[i, "r"]; cc <- pairs[i, "cc"]
        z <- .ctma_draw_summary(cor_draws[, r, cc], probs)
        data.frame(block = spec$block[r],
                   structure = fit$prepared$covariance_blocks$structure[
                     spec$block_id[r]],
                   row_re = r, col_re = cc,
                   row_parameter = spec$parameter[r], col_parameter = spec$parameter[cc],
                   mean = z[1], lower = z[2], median = z[3], upper = z[4])
      }))
  }

  t0var_summary <- t0cov_summary <- data.frame()
  if (!is.null(d$studyre_T0VAR)) {
    selected_t0 <- which(spec$matrix == "T0VAR" & spec$type %in% c(2L, 3L))
    if (length(selected_t0)) {
      t0var_summary <- do.call(rbind, lapply(seq_along(studies), function(s)
        do.call(rbind, lapply(selected_t0, function(q) {
          z <- .ctma_draw_summary(d$studyre_T0VAR[, s, spec$t0_row[q], spec$t0_col[q]], probs)
          data.frame(study = studies[s], parameter = spec$parameter[q],
                     row = spec$row[q], col = spec$col[q],
                     mean = z[1], lower = z[2], median = z[3], upper = z[4])
        }))))
      t0cov_summary <- do.call(rbind, lapply(seq_along(studies), function(s)
        do.call(rbind, lapply(selected_t0, function(q) {
          z <- .ctma_draw_summary(d$studyre_T0cov[, s, spec$t0_row[q], spec$t0_col[q]], probs)
          data.frame(study = studies[s], parameter = spec$parameter[q],
                     row = spec$row[q], col = spec$col[q],
                     mean = z[1], lower = z[2], median = z[3], upper = z[4])
        }))))
    }
  }

  list(specification = spec,
       covariance_blocks = fit$prepared$covariance_blocks,
       study_effects_raw = effect_summary,
       study_sd_raw = sd_summary,
       study_sd_transformed = transformed_sd_summary,
       study_correlations = correlation_summary,
       study_parameters = study_parameter_summary,
       population_parameters = population_parameter_summary,
       study_T0VAR = t0var_summary, study_T0cov = t0cov_summary,
       ctsem_draws = d)
}

# Backward-compatible aliases. New code should use ctmaRePrep(), ctmaReFit(),
# and summary().
ctma_prepare_study_re <- ctmaRePrep
ctma_fit_study_re <- ctmaReFit
ctma_summary_study_re <- .ctmaReSummary

# Example for the revised mf2.rds ------------------------------------------
# mf2 <- readRDS("mf2.rds")
#
# # Automatically reads mf2$argumentList$randomEffect and TI column "study":
# prepared <- ctmaRePrep(mf2, covariance = "independent")
#
# # Explicit matrices override the settings stored in argumentList:
# # prepared <- ctmaRePrep(mf2, reT0VAR = matrix(TRUE, 1, 1))
#
# fit <- ctmaReFit(prepared, estimation = "bayes", cores = 4)
# mapfit <- ctmaReFit(prepared, estimation = "map", cores = 4)
# mlfit <- ctmaReFit(prepared, estimation = "ml", cores = 4)
# results <- summary(fit)
