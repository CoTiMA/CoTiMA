#' ctmaFitList
#'
#' @description Combines CoTiMAFit objects into a list with class CoTiMAFit to inform generic functions what to do
#'
#' @param ... any number of CoTiMAFit objects
#'
#' @examples
#' plot(ctmaFitList(CoTiMAInitFit_3, CoTiMAFullFit_3),
#'      timeUnit="Months",
#'      timeRange=c(1, 144, 1) )
#'
#' @export ctmaFitList
#'
ctmaFitList <- function(...) {
  x <- list(...)
  class(x) <- "CoTiMAFit"
  return(x)
}

