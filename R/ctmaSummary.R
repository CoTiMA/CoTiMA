# Summary function for CoTiMA results
#' ctmaSummary
#'
#' @param x A CoTiMA Fit object
#' @method summary ctStanFit
#' @export
#'
summary.CoTiMAFit <- function(x)
  {
  #UseMethod("summary")
  print(x$summary)
  }
