#' ctmaFitList
#'
#' @description Combines CoTiMAFit objects into a list with class CoTiMAFit to inform generic functions what to do
#'
#' @param ... any number of CoTiMAFit objects
#'
#' @export ctmaFitList
#'
ctmaFitList <- function(...) {
  x <- list(...)
  class(x) <- "CoTiMAFit"
  return(x)
}

