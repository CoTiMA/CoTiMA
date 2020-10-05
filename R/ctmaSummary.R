# Summary function for CoTiMA results
#' ctmaSummary
#'
#' @param x ""
#'
#' @return
#'
summary.CoTiMAFit <- function(x)
  {
  UseMethod("summary")
  print(x$summary)
  }
