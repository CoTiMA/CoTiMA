# Summary function for CoTiMA results
#' summary.CoTiMAFit
#'
#' @param x líst
#' @method summary CoTiMAFit
#' @export
#'
summary.CoTiMAFit <- function(x, ...) {
  return(print(x$summary))
}

