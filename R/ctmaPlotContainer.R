#' plot.CoTiMAFit
#'
#' @description call ctmaPlot if a CoTiMAFit object is supplied to plot()
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

