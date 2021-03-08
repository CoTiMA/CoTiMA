#' ctmaCorRel
#'
#' @description Disattenuates the entries in a correlation matrix using a vector of reliabilities.
#'
#' @param empcov Empirical correlation matrix
#' @param alphas Vector reliabilities
#'
#' @examples
#' empcov313new <- ctmaCorRel(empcov=empcov313, alphas=alphas313)
#'
#' @export ctmaCorRel
#'
#' @return A corrected correlation matrix (corEmpcov). Corrections leading to r > 1.0 are set to 1.0.
#'
ctmaCorRel <- function(empcov=NULL, alphas=NULL) {
  if (length(alphas) != dim(empcov)[1]) {
    ErrorMsg <- "\nNumber of alphas provided does not equal number of variables in correlation matrix. \nGood luck for the next try!"
    stop(ErrorMsg)
  }
  if (any(is.na(alphas))) {
    Msg <- "One or more alpha was NA. I will treat this at 1.0 (no disattenuation for these varaibles). \n"
    message(Msg)
  }
  alphasMat <- matrix((1/(alphas^.5 %x% alphas^.5)), length(alphas))
  corEmpcov <- empcov * alphasMat
  corEmpcov[corEmpcov > 1] <- 1
  return(corEmpcov)
}
