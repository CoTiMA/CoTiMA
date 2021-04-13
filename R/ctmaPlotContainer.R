#' plot.CoTiMAFit
#'
#' @description call \code{\link{ctmaPlot}} if a CoTiMAFit object is supplied to plot()
#'
#' @param x líst
#' @param ... further arguments to be passed through to summary()
#'
#' @method plot CoTiMAFit
#'
#' @return returns a call to 'ctmaPlot', which is used to plot CoTiMA fit objects
#'
#' @export
#'
plot.CoTiMAFit <- function(x, ...) {
  return(ctmaPlot(x, ...))
}

