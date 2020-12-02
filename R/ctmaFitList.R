# Function to combine CoTiMAFits into a list with class CoTiMAFit (so that generic functions know what to do)
#' ctmaFitList
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

