#' Summary method for CoTiMA fit objects
#'
#' @description
#' Provides the standard summary for conventional CoTiMA fit objects and
#' redirects moment-based and hierarchical study-level random-effects fits to
#' their corresponding summary methods.
#'
#' @param object A fitted object inheriting from class \code{"CoTiMAFit"}.
#' @param probs Numeric vector containing the lower, median, and upper
#'   probabilities used for random-effects posterior summaries.
#' @param ... Further arguments passed to the applicable summary method.
#'
#' @return
#' For a conventional CoTiMA fit, the stored \code{object$summary} is printed
#' and returned. For a hierarchical study-level random-effects fit, a list
#' containing population, study-specific, raw-scale and transformed-scale
#' heterogeneity, correlation, and Time-0 parameter summaries is returned.
#'
#' @method summary CoTiMAFit
#' @export
#'
summary.CoTiMAFit <- function(
    object,
    ...,
    probs = c(.025, .5, .975)) {

  if (!inherits(object, "CoTiMAFit")) {
    stop("Not a CoTiMAFit object!", call. = FALSE)
  }

  # Post-hoc empirical moment estimates returned by ctmaReMoments().
  if (inherits(object, "CoTiMAReMoments")) {
    return(.ctmaReMomentsSummary(
      object,
      probs = object$probs,
      print_note = TRUE
    ))
  }

  # Hierarchical study-level random-effects fits returned by ctmaReFit()
  # or ctmaReMLFit().
  if (inherits(object, "CoTiMAStudyREFit")) {

    if (!exists(".ctmaReSummary", mode = "function")) {
      stop(
        "The internal CoTiMA random-effects summary function ",
        "'.ctmaReSummary()' is not available.",
        call. = FALSE
      )
    }

    return(
      .ctmaReSummary(
        object,
        probs = probs
      )
    )
  }

  # Established behaviour for all conventional CoTiMAFit objects.
  print(object$summary)
}
