#' ctmaFitList
#'
#' @description Combines CoTiMAFit objects into a list with class CoTiMAFit to inform generic functions what to do
#'
#' @param ... any number of CoTiMAFit objects
#'
#' @examples
#' \dontrun{
#' CoTiMAInitFit_3$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAFullFit_3$activeDirectory <- "/Users/tmp/" # adapt!
#' plot(ctmaFitList(CoTiMAInitFit_3, CoTiMAFullFit_3),
#'      timeUnit="Months",
#'      timeRange=c(1, 144, 1) )
#'      }
#'
#' @export ctmaFitList
#'
#' @return a list that combines all objects supplied and is assigned the class 'CoTiMAFit'
#'
ctmaFitList <- function(...) {
  x <- list(...)
  class(x) <- "CoTiMAFit"
  return(x)
}

