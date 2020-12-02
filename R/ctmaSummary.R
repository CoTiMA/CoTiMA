# Summary function for CoTiMA results
#' summary.CoTiMAFit
#'
#' @param object one CoTiMAFit object or more as ctmaFitList(object1, object2, ...)
#' @param ... further arguments to be passed through to summary()
#'
#' @method summary CoTiMAFit
#'
#' @export
#'
summary.CoTiMAFit <- function(object, ...) {
  return(print(object$summary))
}

