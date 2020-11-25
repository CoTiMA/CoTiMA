# Summary function for CoTiMA results
#' summary.CoTiMAFit
#'
#' @param x líst
#' @method summary CoTiMAFit
#' @export
#'
summary.CoTiMAFit <- function(x) {
  #UseMethod("print")
  #y <- x$summary
  #class(y) <- "character"
  #for (i in 1:length(y)) class(y[[i]]) <- "print"
  print(x$summary)
  #return(y)
}

