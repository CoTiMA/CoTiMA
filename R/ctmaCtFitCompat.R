# Internal compatibility helpers for ctsem fitting functions.

.ctmaSanitizeCtFitArgs <- function(args) {
  if (!is.list(args)) {
    stop("Internal ctsem fit arguments must be supplied as a list.", call. = FALSE)
  }
  if ("vb" %in% names(args)) {
    if (isTRUE(args$vb)) {
      stop(
        "Variational Bayes is no longer supported by ctsem. Use ",
        "optimize = TRUE for maximum likelihood/MAP estimation or ",
        "optimize = FALSE for Bayesian sampling.",
        call. = FALSE
      )
    }
    # Recent ctsem versions reject the argument even when it is FALSE.
    args$vb <- NULL
  }
  args
}

# ctsem 3.x uses Stan without a backend argument. Newer versions may
# default to Julia, so select Stan explicitly when the API permits it.
.ctmaStanCall <- function(fit, args) {
  args <- .ctmaSanitizeCtFitArgs(args)
  if ("backend" %in% names(formals(fit))) args$backend <- "stan"
  do.call(fit, .ctmaRetiredOptimControls(args))
}

.ctmaCtFitCompat <- function(...) {
  .ctmaStanCall(ctsem::ctFit, list(...))
}

.ctmaCtStanFitCompat <- function(...) {
  .ctmaStanCall(ctsem::ctStanFit, list(...))
}

# Keep historical defaults usable when newer ctsem removes inert controls.
# Explicit requests for retired behaviour are rejected rather than ignored.
.ctmaRetiredOptimControls <- function(args, checker = NULL) {
  controls <- args$optimcontrol
  if (!is.list(controls)) return(args)
  retired <- list(parstepsAutoModel = FALSE, groupFreeThreshold = 0.5)
  present <- intersect(names(retired), names(controls))
  if (!length(present)) return(args)
  if (is.null(checker)) {
    ns <- asNamespace("ctsem")
    if (!exists(".ctFitCheckControls", ns, inherits = FALSE)) return(args)
    checker <- get(".ctFitCheckControls", ns, inherits = FALSE)
  }
  backend <- if (is.null(args$backend)) "stan" else args$backend
  for (name in present) {
    probe <- tryCatch(checker(stats::setNames(list(controls[[name]]), name),
      backend = backend), error = function(e) e)
    removed <- inherits(probe, "error") &&
      grepl("Unrecognised optimcontrol name", conditionMessage(probe), fixed = TRUE)
    if (!removed) next
    if (!isTRUE(all.equal(controls[[name]], retired[[name]]))) stop(
      "This ctsem version no longer supports optimcontrol$", name,
      ". Remove the retired control; use supported explicit parsteps if needed. ",
      "A non-default request was not silently discarded.", call. = FALSE)
    controls[[name]] <- NULL
  }
  args$optimcontrol <- controls
  args
}
