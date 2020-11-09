#' ctmaCorRel
#'
#' @param empcov "Empirical correlation matrix"
#' @param alpha "Vector reliabilities"
#'
#' @return "Empirical correlation matrix corrected for unreliability"
#' @export ctmaCorRel
#'
#'
ctmaCorRel <- function(empcov=NULL, alpha=NULL) {
  alphaMat <- matrix((1/(alpha^.5 %x% alpha^.5)), length(alpha))
  corEmpcov <- empcov * alphaMat
  corEmpcov[corEmpcov > 1] <- 1
  return(corEmpcov)
}
