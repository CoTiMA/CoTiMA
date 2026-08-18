# Experimental parallel multistart ML fitting for CoTiMA study random effects
#
# This file intentionally does not modify CoTiMA_RE_hierarchical.R.
# It requires that file because model preparation and result transformation are
# shared with the Bayesian implementation.
#
# IMPORTANT STATISTICAL QUALIFICATION
# -----------------------------------
# This maximizes the joint/conditional Stan objective over both ordinary model
# parameters and latent study effects.  It does not integrate the study effects
# out and is therefore not marginal maximum likelihood (FIML).  In addition,
# the chosen study-effect parameterisation affects the joint mode. Treat this
# as an experimental optimizer
# and compare its result with the Bayesian fit.
#
# Multiple cores are used for independent starts.  A single rstan::optimizing()
# run is not parallel.  Parallel starts can nevertheless improve robustness to
# local optima and make useful use of several cores.

# Minimum preparation version required by this optimizer. Earlier prepared
# objects may contain the reversed whenvecp coding that silently prevented
# study-level raw parameter effects from entering ctsem's subject matrices.
.ctma_parallel_required_re_version <- "2026-08-10.17"

.ctma_parallel_validate_prepared <- function(prepared) {
  spec <- prepared$specification
  standata <- prepared$standata
  raw_re <- which(spec$type == 1L)

  if (length(raw_re)) {
    raw_indices <- unique(as.integer(spec$index[raw_re]))
    observed <- as.integer(standata$whenvecp[2L, raw_indices])
    if (!identical(observed, raw_indices)) stop(
      "The prepared object uses stale or invalid ctsem scheduling: selected ",
      "study-RE raw parameters must contain their own indices in ",
      "standata$whenvecp[2, ]. Re-source CoTiMA_RE_hierarchical.R version ",
      .ctma_parallel_required_re_version,
      " or later and recreate the object with ctmaRePrep()."
    )

    mat_rows <- which(standata$matsetup[, 3L] %in% raw_indices)
    if (!length(mat_rows) || any(standata$matsetup[mat_rows, 6L] == 0L)) stop(
      "The prepared object does not activate ctsem's subject-level matsetup ",
      "gate for every selected study random effect. Recreate it with ",
      "ctmaRePrep() from CoTiMA_RE_hierarchical.R version ",
      .ctma_parallel_required_re_version, " or later."
    )
  }
  invisible(TRUE)
}

.ctma_parallel_validate_init <- function(init, prepared) {
  if (is.character(init) && length(init) == 1L) return(init)
  if (is.numeric(init) && length(init) == 1L) return(init)
  if (!is.list(init)) stop(
    "init must be 'random', a numeric radius, a list, or a function returning ",
    "one of these."
  )

  expected_ti <- as.integer(prepared$standata$ntipredeffects)
  if (!is.null(init$tipredeffectparams) &&
      length(init$tipredeffectparams) != expected_ti) stop(
    "The initialization contains ", length(init$tipredeffectparams),
    " tipredeffectparams, but the prepared model requires ", expected_ti,
    ". This commonly occurs after removing fixed study-indicator effects."
  )

  q <- as.integer(prepared$standata$nstudyre)
  k <- as.integer(prepared$standata$nstudies)
  if (!is.null(init$studyre_sd)) {
    if (length(init$studyre_sd) != q) stop(
      "studyre_sd initialization must have length nstudyre = ", q, "."
    )
    # Preserve the dimension required by rstan for vector[1].
    init$studyre_sd <- array(as.numeric(init$studyre_sd), dim = q)
  }
  if (!is.null(init$studyre_effect) &&
      !identical(dim(init$studyre_effect), c(k, q))) stop(
    "For a centered model, studyre_effect initialization must have dimensions ",
    "nstudies x nstudyre = ", k, " x ", q, "."
  )
  if (!is.null(init$studyre_z) &&
      !identical(dim(init$studyre_z), c(q, k))) stop(
    "For a noncentered model, studyre_z initialization must have dimensions ",
    "nstudyre x nstudies = ", q, " x ", k, "."
  )
  init
}

.ctma_try_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    text <- sprintf(...)
    if (exists(".ctma_print_message", mode = "function")) {
      .ctma_print_message(text)
    } else {
      width <- 81L
      left <- "######## "
      right <- " #########"
      content_width <- width - nchar(left) - nchar(right)
      lines <- strwrap(paste0("Just a note: ", text), width = content_width)
      border <- strrep("#", width)
      print(border)
      for (line in lines) print(paste0(
        left, line,
        strrep(" ", max(0L, content_width - nchar(line, type = "width"))),
        right
      ))
      print(border)
    }
  }
  invisible(NULL)
}

.ctma_try_one_start <- function(i, sm, standata, seed, init, optim_args,
    prepared) {
  started <- Sys.time()
  warnings <- character()
  ans <- tryCatch(
    withCallingHandlers({
      this_init <- if (is.function(init)) {
        fml <- names(formals(init))
        if (length(fml)) init(i) else init()
      } else init
      this_init <- .ctma_parallel_validate_init(this_init, prepared)
      call <- c(list(
        object = sm,
        data = standata,
        seed = as.integer(seed + i - 1L),
        init = this_init,
        as_vector = FALSE,
        hessian = FALSE,
        verbose = FALSE
      ), optim_args)
      do.call(rstan::optimizing, call)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(ans, "error")) return(list(
    start = i, ok = FALSE, value = -Inf, elapsed = elapsed,
    error = conditionMessage(ans), warnings = warnings, result = NULL
  ))

  value <- ans$value
  return_code <- ans$return_code
  ok <- length(value) == 1L && is.finite(value) &&
    (is.null(return_code) || identical(as.integer(return_code), 0L)) &&
    is.list(ans$par)
  list(
    start = i, ok = ok,
    value = if (length(value) == 1L && is.finite(value)) value else -Inf,
    elapsed = elapsed,
    error = if (ok) NULL else paste0(
      "Optimization returned no usable finite solution",
      if (!is.null(return_code)) paste0(" (return code ", return_code, ")") else "",
      "."
    ),
    warnings = warnings,
    result = if (ok) ans else NULL
  )
}

.ctma_try_parallel_fork <- function(n_starts, cores, worker,
    poll_seconds, verbose) {
  pending_ids <- seq_len(n_starts)
  completed <- vector("list", n_starts)

  while (length(pending_ids)) {
    batch_ids <- head(pending_ids, cores)
    pending_ids <- setdiff(pending_ids, batch_ids)
    jobs <- lapply(batch_ids, function(i)
      parallel::mcparallel(worker(i), silent = TRUE, mc.set.seed = FALSE))
    names(jobs) <- vapply(jobs, function(x) as.character(x$pid), character(1L))
    active <- jobs

    on.exit({
      if (length(active)) invisible(lapply(active, function(job)
        try(tools::pskill(job$pid), silent = TRUE)))
    }, add = TRUE)

    while (length(active)) {
      got <- parallel::mccollect(active, wait = FALSE)
      if (is.null(got)) {
        Sys.sleep(poll_seconds)
        next
      }
      for (pid in names(got)) {
        z <- got[[pid]]
        if (inherits(z, "try-error")) {
          i <- batch_ids[match(pid, names(jobs))]
          z <- list(start = i, ok = FALSE, value = -Inf, elapsed = NA_real_,
            error = as.character(z), warnings = character(), result = NULL)
        }
        completed[[z$start]] <- z
        .ctma_try_message(verbose,
          "[parallel ML] start %d/%d finished in %.1f s; %s%s",
          z$start, n_starts, z$elapsed,
          if (z$ok) "objective = " else "FAILED",
          if (z$ok) format(z$value, digits = 10) else
            paste0(": ", z$error))
      }
      active <- active[setdiff(names(active), names(got))]
    }
  }
  completed
}

#' Parallel multistart conditional ML for CoTiMA study random effects
#'
#' @description
#' Compiles the prepared Stan model once and runs independent
#' \code{rstan::optimizing()} starts on forked workers. The best finite
#' objective is retained. Study effects are optimized rather than integrated
#' out, so this is joint/conditional and not marginal maximum likelihood.
#'
#' @param prepared Object returned by \code{\link{ctmaRePrep}}.
#' @param cores Maximum number of concurrent optimization workers.
#' @param n_starts Number of independent optimizer starts.
#' @param parameterization Centered, noncentered, or automatic study-effect
#'   parameterization. The prepared Stan text is regenerated if necessary.
#' @param seed Base random-number seed.
#' @param init Initialization policy, list, or function. \code{NULL} uses the
#'   policy stored in \code{prepared}.
#' @param init_jitter Optional override of the stored automatic-init jitter.
#' @param algorithm Optimization algorithm passed to rstan.
#' @param iter Maximum iterations per start.
#' @param refresh Optimizer progress refresh interval.
#' @param poll_seconds Positive interval for checking forked workers.
#' @param compile_args Named arguments passed to \code{rstan::stan_model()}.
#' @param optim_args Additional named arguments passed to
#'   \code{rstan::optimizing()}.
#' @param verbose Print CoTiMA-style progress messages.
#'
#' @return A \code{"CoTiMAFit"} object with subtype
#'   \code{"CoTiMAStudyREFit"}, plus all optimization attempts, the selected
#'   start, and its objective value.
#' @export
#' @seealso \code{\link{ctmaRePrep}}, \code{\link{ctmaReFit}}
#'
#' @examples
#' \dontrun{
#' mlFit <- ctmaReMLFit(prepared, cores = 4, n_starts = 8,
#'   parameterization = "centered", iter = 2000)
#' }
ctmaReMLFit <- function(
    prepared,
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
    n_starts = max(4L, cores),
    parameterization = c("auto", "centered", "noncentered"),
    seed = 1234L,
    init = NULL,
    init_jitter = NULL,
    algorithm = "LBFGS",
    iter = 2000L,
    refresh = 100L,
    poll_seconds = 2,
    compile_args = list(),
    optim_args = list(),
    verbose = TRUE) {

  if (!inherits(prepared, "ctmaStudyREModel")) stop(
    "prepared must be returned by ctmaRePrep()."
  )
  .ctma_parallel_validate_prepared(prepared)
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required.")
  if (is.null(init)) init <- prepared$init
  if (is.null(init_jitter)) init_jitter <- prepared$init_jitter
  .ctma_try_message(verbose,
    "Using CoTiMA_RE_hierarchical_parallel.R version %s.",
    attr(ctmaReMLFit, "script_version"))
  if (!is.null(parameterization)) {
    parameterization_requested <- match.arg(parameterization)
    resolved <- .ctma_resolve_parameterization(
      parameterization_requested, prepared$standata$study_id,
      verbose = verbose
    )
    if (!identical(resolved, prepared$parameterization)) {
      .ctma_try_message(verbose,
        "Repatching the prepared model from '%s' to '%s' parameterization.",
        prepared$parameterization, resolved)
      prepared$stanmodeltext <- .ctma_patch_stan(
        prepared$original$stanmodeltext, prepared$covariance, resolved,
        prepared$covariance_blocks
      )
      prepared$parameterization <- resolved
      prepared$parameterization_requested <- parameterization_requested
    }
  }
  # Resolve automatic initialization only after any parameterization change;
  # the required Stan parameter is studyre_effect for centered models and
  # studyre_z for noncentered models.
  if (is.character(init) && identical(init, "auto"))
    init <- .ctma_init_function(prepared, init_jitter)
  if (.Platform$OS.type == "windows") stop(
    "This experimental implementation uses forked workers and currently runs ",
    "only on macOS/Linux. A Windows PSOCK implementation would need to load or ",
    "compile the Stan model independently in every worker."
  )
  cores <- as.integer(cores)
  n_starts <- as.integer(n_starts)
  if (cores < 1L || n_starts < 1L)
    stop("cores and n_starts must be positive integers.")
  cores <- min(cores, n_starts)
  if (!is.numeric(poll_seconds) || length(poll_seconds) != 1L ||
      !is.finite(poll_seconds) || poll_seconds <= 0)
    stop("poll_seconds must be one positive finite number.")

  forbidden <- intersect(names(optim_args), c(
    "object", "data", "seed", "init", "as_vector", "hessian", "verbose"
  ))
  if (length(forbidden)) stop(
    "These optim_args are managed by the wrapper and may not be supplied: ",
    paste(forbidden, collapse = ", ")
  )
  optim_args <- utils::modifyList(list(
    algorithm = algorithm, iter = as.integer(iter),
    refresh = as.integer(refresh)
  ), optim_args)

  .ctma_try_message(verbose,
    "[parallel ML 1/4] Compiling the modified Stan model once in the parent process.")
  compile_call <- utils::modifyList(list(
    model_code = prepared$stanmodeltext,
    model_name = paste0("cotima_study_re_", prepared$covariance, "_",
      prepared$parameterization, "_ml_parallel"),
    auto_write = TRUE
  ), compile_args)
  sm <- do.call(rstan::stan_model, compile_call)

  standata <- prepared$standata
  if (is.null(standata$priors)) stop(
    "The ctsem Stan data do not contain the expected 'priors' switch."
  )
  standata$priors <- 0L
  standata$studyre_use_hyperpriors <- 0L

  .ctma_try_message(verbose,
    "[parallel ML 2/4] Launching %d independent starts on up to %d cores.",
    n_starts, cores)
  worker <- function(i) .ctma_try_one_start(
    i = i, sm = sm, standata = standata, seed = seed,
    init = init, optim_args = optim_args, prepared = prepared
  )
  attempts <- .ctma_try_parallel_fork(
    n_starts = n_starts, cores = cores, worker = worker,
    poll_seconds = poll_seconds, verbose = verbose
  )

  ok <- vapply(attempts, function(x) isTRUE(x$ok), logical(1L))
  if (!any(ok)) stop(
    "All parallel ML starts failed. Inspect attr(error, 'attempts') is not ",
    "possible after a simple stop; rerun with fewer starts and refresh = 1 ",
    "to diagnose the optimizer output."
  )
  values <- vapply(attempts, `[[`, numeric(1L), "value")
  best_i <- which.max(values)
  best <- attempts[[best_i]]$result

  .ctma_try_message(verbose,
    "[parallel ML 3/4] Selected start %d with objective %.10g (%d/%d succeeded).",
    best_i, values[best_i], sum(ok), n_starts)

  # rstan::optimizing(as_vector = FALSE) returns constrained parameters as a
  # named list. Add the leading one-draw dimension expected by the existing
  # CoTiMA study-RE summary function.
  draws <- .ctma_add_draw_dimension(best$par)
  out <- structure(list(
    prepared = prepared,
    stanmodel = sm,
    estimation = "ml",
    optimizer_uncertainty = "none",
    method = "parallel_multistart_joint_ml",
    engine_fit = best,
    draws = draws,
    attempts = attempts,
    best_start = best_i,
    objective = values[best_i],
    parameterization = prepared$parameterization,
    ml_definition = paste(
      "Joint/conditional mode with ctsem priors and study hyperpriors disabled;",
      "study effects are not integrated out."
    )
  ), class = c("CoTiMAFit", "CoTiMAStudyREFit", "ctmaStudyREFit"))

  .ctma_try_message(verbose,
    "[parallel ML 4/4] Experimental fit object is ready. Use summary().")
  out
}

# Backward-compatible alias. New code should use ctmaReMLFit().
ctma_fit_study_re_ml_parallel <- ctmaReMLFit

attr(ctmaReMLFit, "script_version") <- "2026-08-11.5"

# Example
# -------
# source("/Users/cdormann/Documents/CoTiMA/CoTiMA_RE_hierarchical.R")
# source("/Users/cdormann/Documents/CoTiMA/CoTiMA_RE_hierarchical_parallel.R")
#
# prepared <- ctmaRePrep(fit, covariance = "independent")
# mlfit <- ctmaReMLFit(
#   prepared,
#   cores = 4,
#   n_starts = 8,       # two waves of four independent optimizations
#   parameterization = "centered",
#   iter = 2000,
#   refresh = 100
# )
# results <- summary(mlfit)
#
# Compare objective values and failures from all starts:
# vapply(mlfit$attempts, `[[`, numeric(1L), "value")
# lapply(mlfit$attempts, function(x) x[c("ok", "error", "warnings")])
