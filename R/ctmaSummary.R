# Summary function for CoTiMA results
#' ctmaSummary
#'
#' @param x ""
#'
#'
summary.CoTiMAFit <- function(x)
  {
  #UseMethod("summary")
  print(x$summary)
  }
