# Plot function for CoTiMA results
#' plot.CoTiMAFit
#'
#' @param x líst
#' @param ... further arguments to be passed through to summary()
#'
#' @method plot CoTiMAFit
#'
#' @export
#'
plot.CoTiMAFit <- function(x, ...) {
  return(ctmaPlot(x, ...))
}

