#' summary.CoTiMAFit
#'
#' @description defines summary for CoTiMAFit objects
#'
#' @param object one CoTiMAFit object or more as ctmaFitList(object1, object2, ...)
#' @param ... further arguments to be passed through to summary()
#'
#' @method summary CoTiMAFit
#'
#' @export
#'
summary.CoTiMAFit <- function(object, ...) {

  if(!'CoTiMAFit' %in% class(object)) stop('Not a CoTiMAFit object!')

  return(print(object$summary))
}

