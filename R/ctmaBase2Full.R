#' ctmaBase2Full
#'
#' @param baseMat ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaBase2Full <- function(baseMat=NULL) {
  result <- OpenMx::vec2diag(exp(OpenMx::diag2vec(baseMat))) + baseMat - OpenMx::vec2diag(OpenMx::diag2vec(baseMat))
  result <- result %*% t(result)
  return(result) # @ instead of $ introduced
}
