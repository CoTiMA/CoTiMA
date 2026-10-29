# Internal, shared power-analysis pipeline. Statistical calculations are kept
# separate from fitting, file I/O and presentation. Matrix coordinates, not
# vector positions or parameter labels, identify outcomes and predictors.

.ctmaPowerMessage <- function(text) {
  message("#################################################################################\n",
    "CoTiMA power analysis: ", text, "\n",
    "#################################################################################")
}

.ctmaPowerDesign <- function(init, timeRange, statisticalPower, failSafeN, failSafeP,
    activeDirectory, skipScaling, useSampleFraction, digits) {
  if (is.null(init)) stop("A fitted CoTiMA object (ctmaInitFit) has to be supplied.", call. = FALSE)
  n <- init$n.latent
  if (length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n))
    stop("Power analysis requires at least two latent variables.", call. = FALSE)
  if (!length(statisticalPower)) statisticalPower <- c(.9, .8)
  if (!is.numeric(statisticalPower) || any(!is.finite(statisticalPower)) ||
      any(statisticalPower <= 0 | statisticalPower >= 1) || anyDuplicated(statisticalPower))
    stop("statisticalPower must contain distinct finite probabilities strictly between 0 and 1.", call. = FALSE)
  if (!isTRUE(skipScaling) || !is.null(useSampleFraction)) stop(
    "skipScaling=FALSE and useSampleFraction are not supported by the invariant fitting interface. Prepare/rescale or subset the input studies explicitly instead.", call. = FALSE)
  if (length(digits) != 1L || !is.finite(digits) || digits < 0 || digits > 15 || digits != floor(digits))
    stop("digits must be an integer between 0 and 15.", call. = FALSE)
  sizes <- as.numeric(unlist(init$statisticsList$allSampleSizes))
  waves <- as.numeric(unlist(init$statisticsList$allTpoints))
  if (!length(sizes) || length(waves) != length(sizes) || any(!is.finite(sizes)) ||
      any(sizes <= 0) || any(!is.finite(waves)) || any(waves < 2 | waves != floor(waves)))
    stop("Each study must supply a positive sample size and at least two measurement occasions.", call. = FALSE)
  intervals <- init$statisticsList$allDeltas
  if (is.list(intervals) && length(intervals) == length(waves) &&
      identical(as.integer(lengths(intervals)), as.integer(waves - 1))) {
    intervals <- lapply(intervals, as.numeric)
  } else {
    flat <- as.numeric(unlist(intervals))
    if (length(flat) != sum(waves - 1)) stop(
      "Study-specific time intervals do not match allTpoints; their mapping is ambiguous.", call. = FALSE)
    intervals <- split(flat, rep(seq_along(waves), waves - 1))
  }
  observed <- as.numeric(unlist(intervals))
  if (any(!is.finite(observed)) || any(observed <= 0)) stop("Study time intervals must be finite and positive.", call. = FALSE)
  if (is.null(timeRange)) timeRange <- c(0, 1.5 * max(observed), 1)
  if (!is.numeric(timeRange) || length(timeRange) != 3L || any(!is.finite(timeRange)) ||
      timeRange[1] < 0 || timeRange[2] < timeRange[1] || timeRange[3] <= 0)
    stop("timeRange must be a finite nonnegative start/end/positive-step triple.", call. = FALSE)
  if ((timeRange[2] - timeRange[1]) / timeRange[3] > 1e5)
    stop("timeRange requests more than 100000 grid steps; choose a coarser step.", call. = FALSE)
  times <- sort(unique(c(seq(timeRange[1], timeRange[2], by = timeRange[3]),
    observed[observed >= timeRange[1] & observed <= timeRange[2]])))
  times <- times[times > 0]
  if (!length(times)) stop("timeRange must include at least one positive lag.", call. = FALSE)
  if (is.null(failSafeN)) failSafeN <- round(mean(sizes))
  if (is.null(failSafeP)) failSafeP <- .01
  if (length(failSafeN) != 1L || !is.finite(failSafeN) || failSafeN <= n + 1)
    stop("failSafeN must exceed the number of predictors plus one.", call. = FALSE)
  if (length(failSafeP) != 1L || !is.finite(failSafeP) || failSafeP <= 0 || failSafeP >= 1)
    stop("failSafeP must be strictly between 0 and 1.", call. = FALSE)
  if (is.null(activeDirectory)) activeDirectory <- init$activeDirectory
  if (!is.character(activeDirectory) || length(activeDirectory) != 1L ||
      is.na(activeDirectory) || !dir.exists(activeDirectory))
    stop("activeDirectory must name an existing directory.", call. = FALSE)
  labels <- outer(seq_len(n), seq_len(n), function(row, col) paste0("V", col, "toV", row))
  targets <- expand.grid(row = seq_len(n), col = seq_len(n))
  targets <- targets[targets$row != targets$col, , drop = FALSE]
  targets$parameter <- labels[as.matrix(targets)]
  list(n = as.integer(n), sizes = sizes, intervals = intervals, times = times,
    probabilities = statisticalPower, failSafeN = failSafeN, failSafeP = failSafeP,
    directory = normalizePath(activeDirectory), labels = labels, targets = targets)
}

.ctmaPowerGaussian <- function(init, pooled) {
  if (!is.null(init$n.manifest) && init$n.manifest > 0 && init$n.manifest != init$n.latent)
    stop("This state-regression power projection requires one directly observed variable per state; latent measurement/augmented-state models require simulation-based power analysis.", call. = FALSE)
  objects <- c(list(init, pooled), init$studyFitList, pooled$studyFitList)
  types <- list()
  for (object in objects) {
    if (!is.list(object)) next
    for (model in list(object$ctModel, object$ctstanmodel, object$ctstanmodelbase)) {
      if (is.list(model) && !is.null(model$manifesttype)) types <- c(types, list(model$manifesttype))
      if (is.list(model) && is.matrix(model$LAMBDA)) {
        loading <- suppressWarnings(matrix(as.numeric(model$LAMBDA), nrow(model$LAMBDA)))
        if (!identical(dim(loading), c(as.integer(init$n.latent), as.integer(init$n.latent))) ||
            any(!is.finite(loading)) || max(abs(loading - diag(init$n.latent))) > 1e-10)
          stop("Non-identity or freely estimated measurement loadings require simulation-based power analysis.", call. = FALSE)
      }
      if (is.list(model) && !is.null(model$MANIFESTVAR)) {
        measurement_error <- suppressWarnings(as.numeric(model$MANIFESTVAR))
        if (any(!is.finite(measurement_error)) || any(abs(measurement_error) > 1e-10))
          stop("This state-regression power projection does not include measurement-error variance; use simulation-based power analysis.", call. = FALSE)
      }
    }
    if (!is.null(object$argumentList$binaries)) types <- c(types, list(object$argumentList$binaries))
  }
  codes <- suppressWarnings(as.numeric(unlist(types)))
  if (length(codes) && (any(!is.finite(codes)) || any(codes != 0))) stop(
    "ctmaPower's Gaussian covariance/regression approximation does not support binary or ordinal indicators. Use a simulation-based power analysis for this design.", call. = FALSE)
  list(verified = length(codes) > 0L,
    reason = if (length(codes)) "Continuous Gaussian indicators recorded in model metadata."
      else "Indicator-type metadata is unavailable; a Gaussian covariance model is assumed.")
}

.ctmaPowerMatrices <- function(fit, n) {
  summary <- fit$summary
  if (is.null(summary$parmatrices)) summary <- fit$resultsSummary
  if (is.null(summary$parmatrices) && is.null(fit$modelResults)) summary <- summary(fit)
  table <- summary$parmatrices
  read_matrix <- function(names, fallback, symmetric) {
    if (!is.null(table) && all(c("matrix", "row", "col") %in% colnames(table))) {
      mean_col <- intersect(c("Mean", "mean"), colnames(table))
      for (name in names) {
        rows <- which(table[, "matrix"] == name)
        if (!length(rows) || !length(mean_col)) next
        index <- cbind(as.integer(table[rows, "row"]), as.integer(table[rows, "col"]))
        if (anyNA(index) || any(index < 1 | index > n) || anyDuplicated(data.frame(index)))
          stop("Ambiguous matrix coordinates for ", name, ".", call. = FALSE)
        value <- matrix(NA_real_, n, n)
        value[index] <- as.numeric(table[rows, mean_col[1]])
        if (symmetric) {
          missing <- is.na(value); value[missing] <- t(value)[missing]
        }
        if (!anyNA(value)) return(value)
      }
    }
    value <- fit$modelResults[[fallback]]
    if (!is.matrix(value) || !identical(dim(value), c(n, n))) stop(
      "Cannot extract a complete ", fallback, " matrix with explicit coordinates or matrix dimensions.", call. = FALSE)
    value
  }
  drift <- read_matrix("DRIFT", "DRIFT", FALSE)
  diffusion <- read_matrix("DIFFUSIONcov", "DIFFUSION", TRUE)
  initial <- read_matrix(c("T0cov", "T0VAR"), "T0VAR", TRUE)
  if (any(!is.finite(c(drift, diffusion, initial)))) stop("Model matrices must be finite.", call. = FALSE)
  for (name in c("diffusion", "initial")) {
    value <- get(name)
    if (!isTRUE(all.equal(value, t(value), tolerance = 1e-8))) stop(name, " must be a covariance matrix, not a Cholesky factor.", call. = FALSE)
    if (min(eigen(value, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 * max(1, abs(value)))
      stop(name, " is not positive semidefinite.", call. = FALSE)
  }
  scale <- fit$CoTiMAStanctArgs$scaleTime
  if (is.null(scale)) scale <- fit$argumentList$CoTiMAStanctArgs$scaleTime
  if (is.null(scale)) scale <- fit$argumentList$scaleTime
  if (is.null(scale)) scale <- 1
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0)
    stop("The invariant fit has invalid time-scaling metadata.", call. = FALSE)
  # Input times were multiplied by scale during invariant fitting. Restore
  # original units before evaluating the original study lags.
  list(DRIFT = drift * scale, DIFFUSION = diffusion * scale, T0VAR = initial,
    time_scale = scale, summary = summary)
}

.ctmaPowerImplied <- function(matrices, lag) {
  # The augmented exponential integrates exp(A*s) Q exp(A'*s) without
  # inverting a Kronecker sum; zero/singular drift matrices are supported.
  a <- matrices$DRIFT; n <- nrow(a)
  transition <- OpenMx::expm(a * lag)
  generator <- kronecker(a, diag(n)) + kronecker(diag(n), a)
  augmented <- rbind(cbind(generator, as.numeric(matrices$DIFFUSION)), rep(0, n^2 + 1L))
  noise <- matrix(OpenMx::expm(augmented * lag)[seq_len(n^2), n^2 + 1L], n, n)
  initial <- matrices$T0VAR
  covariance <- rbind(cbind(initial, initial %*% t(transition)),
    cbind(transition %*% initial, transition %*% initial %*% t(transition) + noise))
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0) ||
      inherits(try(chol(covariance), silent = TRUE), "try-error")) stop(
    "The implied covariance is not positive definite at lag ", lag, ". Power is not identified for this projection.", call. = FALSE)
  correlation <- stats::cov2cor(covariance)
  dimnames(correlation) <- rep(list(c(paste0("V", seq_len(n), "T0"), paste0("V", seq_len(n), "T1"))), 2)
  correlation
}

.ctmaPowerRegressionModel <- function(n, syntax = NULL) {
  x <- paste0("V", seq_len(n), "T0"); y <- paste0("V", seq_len(n), "T1")
  if (is.null(syntax)) {
    syntax <- c(vapply(y, function(outcome) paste(outcome, "~", paste(x, collapse = " + ")), character(1)),
      unlist(lapply(seq_len(n), function(i) if (i < n) paste(y[i], "~~", y[seq.int(i + 1L, n)]))))
  }
  if (!is.character(syntax) || !length(syntax)) stop("regressionModel must be lavaan syntax or NULL.", call. = FALSE)
  table <- lavaan::lavaanify(syntax, auto = TRUE, model.type = "sem", fixed.x = FALSE,
    meanstructure = FALSE)
  if (!setequal(lavaan::lavNames(table, "ov"), c(x, y))) stop(
    "regressionModel must use exactly V1T0,...,VpT0 and V1T1,...,VpT1 as observed variables.", call. = FALSE)
  if (any(!table$op %in% c("~", "~~"))) stop(
    "Only observed Gaussian regressions and covariance restrictions are supported; latent variables, means, equality/inequality constraints and defined parameters require simulation-based power analysis.", call. = FALSE)
  slopes <- table$op == "~"
  if (any(!table$lhs[slopes] %in% y | !table$rhs[slopes] %in% x)) stop(
    "Regression equations must predict T1 outcomes from T0 predictors only.", call. = FALSE)
  # Shared labels induce multi-parameter hypotheses: the one-coefficient
  # noncentral-t power formula cannot silently be applied to those hypotheses.
  labels <- table$label[nzchar(table$label)]
  if (anyDuplicated(labels)) stop("Shared parameter labels are not supported by the one-coefficient power approximation.", call. = FALSE)
  covariance <- table$op == "~~"
  if (any((table$lhs[covariance] %in% x) != (table$rhs[covariance] %in% x))) stop(
    "Cross-time residual covariances are not supported in the regression power approximation.", call. = FALSE)
  required <- expand.grid(lhs = y, rhs = x, stringsAsFactors = FALSE)
  actual <- paste(table$lhs[slopes], table$rhs[slopes])
  complete_slopes <- length(actual) == n^2 && !anyDuplicated(actual) &&
    setequal(actual, paste(required$lhs, required$rhs)) && all(table$free[slopes] > 0)
  covariance_free <- all(table$free[covariance] > 0)
  complete_covariance <- sum(covariance) == n * (n + 1)
  list(syntax = syntax, table = table, unrestricted = complete_slopes && covariance_free && complete_covariance,
    reason = if (complete_slopes && covariance_free && complete_covariance)
      "All regression coefficients and both within-time covariance matrices are unrestricted."
      else "The regression syntax contains fixed/omitted coefficients or covariance restrictions.")
}

.ctmaPowerDirect <- function(correlation, targets, nobs) {
  n <- ncol(correlation) / 2L
  xx <- correlation[seq_len(n), seq_len(n), drop = FALSE]
  xy <- correlation[seq_len(n), n + seq_len(n), drop = FALSE]
  yy <- correlation[n + seq_len(n), n + seq_len(n), drop = FALSE]
  inverse <- solve(xx)
  coefficients <- t(inverse %*% xy)
  residual <- yy - t(xy) %*% inverse %*% xy
  result <- lapply(seq_len(nrow(targets)), function(k) {
    i <- targets$row[k]; j <- targets$col[k]; retained <- setdiff(seq_len(n), j)
    without <- yy[i, i] - xy[retained, i] %*% solve(xx[retained, retained, drop = FALSE], xy[retained, i])
    variance <- residual[i, i]
    if (variance <= 0) stop("A regression has nonpositive residual variance.", call. = FALSE)
    estimate <- coefficients[i, j]
    # ML expected-information SE, matching lavaan's Gaussian covariance fit.
    se <- sqrt(variance * inverse[j, j] / nobs)
    data.frame(estimate = estimate, se = se, pvalue = 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE),
      R2 = 1 - variance / yy[i, i], R2_without = 1 - as.numeric(without) / yy[i, i], predictors = n, tested = TRUE)
  })
  do.call(rbind, result)
}

.ctmaPowerProbability <- function(N, f2, predictors, alpha = .05, useMBESS = FALSE) {
  if (!is.finite(N) || N <= predictors + 1) return(NA_real_)
  if (is.infinite(f2) && f2 > 0) return(1)
  if (!is.finite(f2) || f2 < 0) stop("The partial effect size must be nonnegative.", call. = FALSE)
  if (useMBESS) return(as.numeric(MBESS::ss.power.reg.coef(Cohen.f2 = f2, p = predictors,
    Specified.N = N, alpha.level = alpha)$Actual.Power))
  # T^2 has a noncentral F distribution: includes BOTH rejection tails of
  # the two-sided t test. This corrects the legacy alpha/2 null-power limit.
  df <- N - predictors - 1
  stats::pf(stats::qf(1 - alpha, 1, df), 1, df, ncp = N * f2, lower.tail = FALSE)
}

.ctmaPowerRequiredN <- function(f2, predictors, target, useMBESS = FALSE) {
  if (is.na(f2)) return(NA_real_)
  probability <- function(N) .ctmaPowerProbability(N, f2, predictors, useMBESS = useMBESS)
  lower <- predictors + 2
  if (probability(lower) >= target) return(lower)
  if (f2 == 0) return(Inf)
  upper <- lower
  maximum <- 2^52
  while (probability(upper) < target) {
    if (upper >= maximum) return(Inf)
    upper <- min(maximum, 2 * upper)
  }
  # Integer bisection avoids scanning successive N, loose root rounding and
  # assuming a sign change when the minimum admissible N already suffices.
  while (upper - lower > 1) {
    middle <- floor(lower + (upper - lower) / 2)
    if (probability(middle) >= target) upper <- middle else lower <- middle
  }
  upper
}

.ctmaPowerF2 <- function(result) {
  increment <- result$R2 - result$R2_without
  tested <- result$tested
  if (any(!is.finite(increment[tested])) || any(increment[tested] < -1e-7) || any(result$R2[tested] >= 1)) stop(
    "Full/reduced regression fits do not define a valid nonnegative partial effect size. Inspect convergence and model restrictions.", call. = FALSE)
  pmax(0, increment) / (1 - result$R2)
}

.ctmaPowerSave <- function(object, path, overwrite) {
  if (file.exists(path) && !isTRUE(overwrite)) stop("File exists: ", path, ". Set silentOverwrite=TRUE to replace it.", call. = FALSE)
  saveRDS(object, path)
}

.ctmaPowerOptima <- function(values, times) {
  valid <- !is.na(values)
  minimum <- if (any(valid)) min(values[valid]) else NA_real_
  tied <- if (is.finite(minimum)) times[valid & values == minimum] else numeric()
  list(minimum = minimum, first = if (length(tied)) tied[1] else NA_real_, all = tied)
}

.ctmaPowerPath <- function(directory, prefix) {
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) || !nzchar(prefix))
    stop("A cache prefix must be one nonempty string.", call. = FALSE)
  file.path(directory, if (grepl("[.]rds$", prefix)) prefix else paste0(prefix, ".rds"))
}

.ctmaPowerLavaan <- function(correlation, model, targets, nobs, directory,
    loadPrefix, savePrefix, lagIndex, overwrite) {
  reference_n <- 1000
  fit_one <- function(table) {
    fit <- lavaan::sem(table, sample.cov = correlation, sample.nobs = reference_n,
      sample.cov.rescale = FALSE, fixed.x = FALSE, meanstructure = FALSE,
      optim.method = "BFGS", control = list(reltol = 1e-12))
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) stop("A lavaan regression did not converge at lag index ", lagIndex, ".", call. = FALSE)
    fit
  }
  full <- fit_one(model$table)
  estimates <- lavaan::parameterEstimates(full)
  residual <- function(fit, outcome) {
    table <- lavaan::parameterEstimates(fit)
    value <- table$est[table$op == "~~" & table$lhs == outcome & table$rhs == outcome]
    if (length(value) != 1L || !is.finite(value) || value <= 0)
      stop("A lavaan regression has invalid residual variance.", call. = FALSE)
    value
  }
  rows <- lapply(seq_len(nrow(targets)), function(k) {
    outcome <- paste0("V", targets$row[k], "T1")
    predictor <- paste0("V", targets$col[k], "T0")
    index <- which(model$table$op == "~" & model$table$lhs == outcome & model$table$rhs == predictor)
    predictors <- sum(model$table$op == "~" & model$table$lhs == outcome & model$table$free > 0)
    if (length(index) != 1L || model$table$free[index] == 0) return(data.frame(
      estimate = NA_real_, se = NA_real_, pvalue = NA_real_, R2 = NA_real_,
      R2_without = NA_real_, predictors = predictors, tested = FALSE))
    reduced_table <- model$table
    keys <- function(table) paste(table$lhs, table$op, table$rhs)
    matched <- match(keys(reduced_table), keys(estimates))
    warm <- reduced_table$free > 0 & !is.na(matched)
    reduced_table$ustart[warm] <- estimates$est[matched[warm]]
    reduced_table$free[index] <- 0L
    reduced_table$ustart[index] <- 0
    reduced_table$label[index] <- ""
    reduced_table$free[reduced_table$free > 0] <- seq_len(sum(reduced_table$free > 0))
    # Build each reduction independently. Never delete predictors with regexes
    # or carry the preceding target's reduced model into the next hypothesis.
    metadata <- list(version = 1L, covariance = correlation, model = reduced_table, nobs = reference_n)
    if (length(loadPrefix)) {
      path <- file.path(directory, paste0(loadPrefix, "_", k, "_", lagIndex, ".rds"))
      cached <- readRDS(path)
      if (!is.list(cached) || !identical(cached$metadata, metadata)) stop(
        "Reduced-fit cache does not match this design/model, or lacks metadata: ", path,
        ". Regenerate it with saveAllInvWOSingFit.", call. = FALSE)
      reduced <- cached$fit
      if (!inherits(reduced, "lavaan") || !isTRUE(lavaan::lavInspect(reduced, "converged")))
        stop("Cached lavaan fit is invalid.", call. = FALSE)
    } else reduced <- fit_one(reduced_table)
    if (length(savePrefix)) .ctmaPowerSave(list(metadata = metadata, fit = reduced),
      file.path(directory, paste0(savePrefix, "_", k, "_", lagIndex, ".rds")), overwrite)
    coefficient <- estimates[estimates$op == "~" & estimates$lhs == outcome & estimates$rhs == predictor, , drop = FALSE]
    se <- coefficient$se * sqrt(reference_n / nobs)
    if (nrow(coefficient) != 1L || !is.finite(se) || se <= 0)
      stop("The target regression coefficient has no valid lavaan SE.", call. = FALSE)
    variance <- correlation[outcome, outcome]
    data.frame(estimate = coefficient$est, se = se,
      pvalue = 2 * stats::pnorm(abs(coefficient$est / se), lower.tail = FALSE),
      R2 = 1 - residual(full, outcome) / variance,
      R2_without = 1 - residual(reduced, outcome) / variance,
      predictors = predictors, tested = TRUE)
  })
  do.call(rbind, rows)
}

.ctmaPowerAnalyse <- function(init, pooled, design, useMBESS, powerMethod,
    regressionModel, loadReduced, saveReduced, overwrite, digits) {
  matrices <- .ctmaPowerMatrices(pooled, design$n)
  gaussian <- .ctmaPowerGaussian(init, pooled)
  model <- .ctmaPowerRegressionModel(design$n, regressionModel)
  cache_requested <- length(loadReduced) > 0L || length(saveReduced) > 0L
  eligible <- gaussian$verified && model$unrestricted && !cache_requested
  if (powerMethod == "covariance" && !eligible) stop(
    "The direct covariance route is unavailable: ", gaussian$reason, " ", model$reason,
    if (cache_requested) " Reduced lavaan caches were requested." else "", call. = FALSE)
  route <- if (powerMethod == "auto") if (eligible) "covariance" else "lavaan" else powerMethod
  note <- if (route == "covariance") paste(
    "Direct covariance-matrix calculations used: the generated discrete-time model is an unrestricted Gaussian multiple regression.",
    "No lavaan regression fits were required.") else paste(
    "Lavaan regression fits used (once per lag and hypothesis, reused across requested power levels).",
    if (powerMethod == "lavaan") "The lavaan route was explicitly requested." else paste(gaussian$reason, model$reason),
    if (cache_requested) "Reduced-fit caching requires the lavaan route." else "")
  .ctmaPowerMessage(note)
  n_times <- length(design$times); n_targets <- nrow(design$targets)
  regression <- vector("list", n_times)
  for (i in seq_len(n_times)) {
    correlation <- .ctmaPowerImplied(matrices, design$times[i])
    regression[[i]] <- if (route == "covariance")
      .ctmaPowerDirect(correlation, design$targets, design$failSafeN) else
      .ctmaPowerLavaan(correlation, model, design$targets, design$failSafeN,
        design$directory, loadReduced, saveReduced, i, overwrite)
    regression[[i]]$f2 <- .ctmaPowerF2(regression[[i]])
  }
  names(regression) <- format(design$times, digits = 15, trim = TRUE, scientific = FALSE)
  if (anyDuplicated(names(regression))) stop("Time lags are too close to label unambiguously.", call. = FALSE)
  sample_sizes <- matrix(NA_real_, n_times, n_targets * length(design$probabilities))
  coefficients <- matrix(NA_real_, n_times, n_targets)
  all_optimal <- list(); significance <- character(n_targets)
  column <- 0L
  for (target in seq_len(n_targets)) {
    effects <- vapply(regression, function(x) x$f2[target], numeric(1))
    predictors <- vapply(regression, function(x) x$predictors[target], numeric(1))
    coefficients[, target] <- vapply(regression, function(x) x$estimate[target], numeric(1))
    for (probability in design$probabilities) {
      column <- column + 1L
      sample_sizes[, column] <- vapply(seq_len(n_times), function(i)
        .ctmaPowerRequiredN(effects[i], predictors[i], probability, useMBESS), numeric(1))
    }
    pvalues <- vapply(regression, function(x) x$pvalue[target], numeric(1))
    significant <- design$times[is.finite(pvalues) & pvalues < design$failSafeP]
    significance[target] <- if (length(significant)) paste0(design$targets$parameter[target],
      ": shortest significant evaluated lag = ", min(significant), "; longest = ", max(significant),
      " (two-sided ML Wald p < ", design$failSafeP, ", N = ", design$failSafeN,
      "). Significance need not hold at every lag between these endpoints.") else paste0(
      design$targets$parameter[target], ": no evaluated lag is significant (or the target is fixed) at p < ",
      design$failSafeP, ", N = ", design$failSafeN, ".")
  }
  colnames(sample_sizes) <- unlist(lapply(design$targets$parameter,
    function(name) paste0(name, " Power=", design$probabilities)))
  rownames(sample_sizes) <- names(regression)
  colnames(coefficients) <- design$targets$parameter; rownames(coefficients) <- names(regression)
  minimum <- optimal <- rep(NA_real_, ncol(sample_sizes))
  for (j in seq_len(ncol(sample_sizes))) {
    selected <- .ctmaPowerOptima(sample_sizes[, j], design$times)
    minimum[j] <- selected$minimum; optimal[j] <- selected$first
    all_optimal[[colnames(sample_sizes)[j]]] <- selected$all
  }
  required <- rbind(cbind(sample_sizes, round(coefficients, digits)),
    "Min N" = c(minimum, rep(NA_real_, n_targets)),
    "Opt. Lag" = c(optimal, rep(NA_real_, n_targets)))
  posthoc <- lapply(seq_len(n_targets), function(target) {
    table <- matrix(NA_real_, length(design$sizes), 1 + 3 * max(lengths(design$intervals)))
    table[, 1] <- design$sizes
    for (study in seq_along(design$sizes)) for (occasion in seq_along(design$intervals[[study]])) {
      lag <- design$intervals[[study]][occasion]
      index <- match(lag, design$times)
      offset <- 2 + (occasion - 1) * 3
      table[study, offset] <- lag
      if (!is.na(index) && regression[[index]]$tested[target]) {
        for (a in seq_along(c(.05, .01))) table[study, offset + a] <- .ctmaPowerProbability(
          design$sizes[study], regression[[index]]$f2[target], regression[[index]]$predictors[target], c(.05, .01)[a], useMBESS)
      }
    }
    colnames(table) <- c("N", rep(c("Time Lag", "Power (\u03b1=.05)", "Power (\u03b1=.01)"), max(lengths(design$intervals))))
    # Unavailable lags/designs remain NA, not zero. This avoids claiming zero
    # power for a study outside the evaluated range or without residual df.
    summary_rows <- matrix(NA_real_, 2, ncol(table), dimnames = list(c("Mean", "Median"), NULL))
    for (j in c(3, 4)) {
      values <- table[, j]; values <- values[is.finite(values)]
      if (length(values)) summary_rows[, j] <- c(mean(values), stats::median(values))
    }
    rownames(table) <- paste0("Study_No_", seq_along(design$sizes))
    round(rbind(table, summary_rows), digits)
  })
  names(posthoc) <- design$targets$parameter
  interpretation <- paste(
    "Plug-in power for discrete-time regression coefficients at evaluated lags; not the power of the pooled continuous-time drift test.",
    "CT parameter uncertainty and clustered/multiwave sampling uncertainty are not propagated.",
    if (!model$unrestricted) "For restricted regressions, partial-R-squared power remains an approximation to the target's sampling distribution; it is not a general SEM power test." else "",
    "CT coefficient constraints are not converted into discrete-time zero constraints.",
    "Fixed/omitted coefficients in an explicit regressionModel have NA power.",
    "Inf required N denotes an unattainable target or a requirement beyond exact-integer numerical resolution.",
    "Opt. Lag reports the first actual minimizing lag; all ties are retained separately.",
    "Post-hoc Mean/Median use available first-interval powers; missing values are not replaced by zero.",
    paste0("Lag calculations use original input time units (invariant-fit scaleTime = ", matrices$time_scale,
      "); legacy caches without scaling metadata are assumed to use original units."))
  probability_note <- if (useMBESS) "MBESS legacy upper-tail probability with a two-sided critical value; integer bisection replaces the sequential N search."
    else "Complete two-sided noncentral-t rejection probability, evaluated through the equivalent noncentral F distribution; integer bisection finds the minimum N."
  .ctmaPowerMessage("Power and sample-size calculations completed.")
  list(matrices = matrices, route = route, check = list(gaussian = gaussian, regression = model$reason,
      unrestricted = model$unrestricted, eligible = eligible), note = note,
    probability_note = probability_note, interpretation = interpretation,
    details = regression, required = required, coefficients = coefficients, posthoc = posthoc,
    all_optimal = all_optimal, significance = significance)
}

.ctmaPowerEstimateTable <- function(matrices, labels) {
  n <- nrow(labels)
  table <- matrices$summary$parmatrices
  families <- list(DRIFT = "DRIFT", DIFFUSION = "DIFFUSIONcov", T0VAR = c("T0cov", "T0VAR"))
  output <- lapply(names(families), function(family) {
    se <- matrix(NA_real_, n, n)
    if (!is.null(table) && all(c("matrix", "row", "col") %in% colnames(table))) {
      se_name <- intersect(c("Sd", "sd", "SE", "se"), colnames(table))
      for (matrix_name in families[[family]]) {
        rows <- which(table[, "matrix"] == matrix_name)
        if (!length(rows) || !length(se_name)) next
        se[cbind(as.integer(table[rows, "row"]), as.integer(table[rows, "col"]))] <- as.numeric(table[rows, se_name[1]])
        if (family != "DRIFT") { missing <- is.na(se); se[missing] <- t(se)[missing] }
        break
      }
    }
    if (family != "T0VAR") se <- se * matrices$time_scale
    estimate <- as.numeric(t(matrices[[family]]))
    se <- as.numeric(t(se))
    as.numeric(t(cbind(estimate, se, estimate / se)))
  })
  output <- do.call(rbind, output)
  rownames(output) <- c("Fixed Effects Drift", "Fixed Effects Diffusion", "Fixed Effects T0Var")
  colnames(output) <- as.vector(rbind(as.vector(t(labels)), "(SE)", "Tvalue"))
  output
}

.ctmaPowerResult <- function(init, pooled, design, analysis, controls, digits) {
  detail <- do.call(rbind, lapply(seq_along(analysis$details), function(i)
    cbind(lag = design$times[i], design$targets, analysis$details[[i]])))
  pooled$resultsSummary <- analysis$matrices$summary
  parameter_names <- init$parameterNames
  parameter_names$DRIFT <- as.vector(design$labels)
  structure(list(activeDirectory = file.path(design$directory, ""),
    backend = controls$backend, backendArgs = controls$backendArgs,
    resource_resolution = controls$resource_resolution, coresToUse = controls$cores,
    plot.type = "power", model.type = "ct", n.studies = 1L, n.latent = design$n,
    studyList = init$studyList, studyFitList = pooled, statisticsList = init$statisticsList,
    emprawList = NULL, modelResults = analysis$matrices[c("DRIFT", "DIFFUSION", "T0VAR")],
    parameterNames = parameter_names,
    scaleTime = controls$scaleTime, argumentList = controls,
    power_method = analysis$route, regression_check = analysis$check,
    regression_details = detail,
    summary = list(model = "Analysis of Statistical Power and Required Sample Sizes",
      calculation_note = analysis$note, probability_note = analysis$probability_note,
      interpretation = analysis$interpretation,
      estimates = list("Estimates of Model with all Effects Invariant" = round(.ctmaPowerEstimateTable(analysis$matrices, design$labels), digits),
        "Requested Statistical Power" = design$probabilities,
        "Power (post hoc) for Drift Effects" = analysis$posthoc,
        "Required Sample Sizes" = analysis$required,
        "Effect Sizes (based on discrete-time calcs; used for power calcs.)" = round(analysis$coefficients, digits),
        "All optimal lags" = analysis$all_optimal,
        "Range of significant effects" = analysis$significance))), class = "CoTiMAFit")
}
