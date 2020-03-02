#' ctmaChangeSVLab
#'
#' @param startValues ?
#' @param ctmodelobj ?
#' @param fixedmodel ?
#' @param noOfStudies ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaChangeSVLab <- function(startValues = NULL, ctmodelobj = NULL, fixedmodel = NULL, noOfStudies = NULL) {
  targetNames1 <- targetNames2 <- c()
  labelPattern <- c(ctmodelobj$DRIFT, ctmodelobj$DIFFUSION, ctmodelobj$T0VAR, ctmodelobj$CINT)
  fixedPattern <- c(fixedmodel$DRIFT, fixedmodel$DIFFUSION, fixedmodel$T0VAR, fixedmodel$CINT)
  targetNames1 <- labelPattern[fixedPattern == "groupfixed"]
  labelPattern <- labelPattern[fixedPattern != "groupfixed"]
  labelPattern <- gsub("T0", "Tx", labelPattern)
  labelPattern <- labelPattern[-grep("0", labelPattern)]
  labelPattern <- gsub("Tx", "T0", labelPattern)
  for (k in 1:noOfStudies) targetNames2 <- c(targetNames2, paste0("Study_No_", k, "_", labelPattern))
  targetNames <- c(targetNames1, targetNames2)
  names(startValues) <- targetNames
  return(startValues)
}
