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

.ctmaCtFitCompat <- function(...) {
  do.call(ctsem::ctFit, .ctmaSanitizeCtFitArgs(list(...)))
}

.ctmaCtStanFitCompat <- function(...) {
  do.call(ctsem::ctStanFit, .ctmaSanitizeCtFitArgs(list(...)))
}
