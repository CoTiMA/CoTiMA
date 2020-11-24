# Summary function for CoTiMA results
#' summary.CoTiMAFit
#'
#' @param x ""
#' @method summary CoTiMAFit
#' @export
#'
summary.CoTiMAFit <- function(x) {
  #UseMethod("print")
  print(x$summary)
  }
