#' ctmaFitToPrep
#'
#' @description Extracts information from fitted CoTiMA objects to (re-)crearte list of primary studies originally created with \code{\link{ctmaPrep}}
#'
#' @param ctmaFitObject  ctmaFitObject
#' @param reUseEmprawData  whether data should be transferred (will be re-used in subsequent fit attempts)
#'
#' @examples
#' newStudyList <- ctmaFitToPrep(CoTiMAInitFit_3)
#'
#' @export ctmaFitToPrep
#'
#' @return list that could be used for fitting new CoTiMA models with \code{\link{ctmaInit}} or \code{\link{ctmaFit}}.
#'
ctmaFitToPrep <- function(ctmaFitObject=NULL, reUseEmprawData=FALSE)
{
  if (is.null(ctmaFitObject$originalStudyNo)) {
    if (!is.null(ctmaFitObject$studyList)) {
      tmp <- ctmaFitObject
      ctmaFitObject <- ctmaFitObject$studyList
    }
  }
  newStudyList <- list()
  newStudyList$studyNumbers <- lapply(ctmaFitObject , function(x) x$originalStudyNo)
  newStudyList$empcovs <- lapply(ctmaFitObject, function(x) x$empcov)
  newStudyList$deltas <- lapply(ctmaFitObject, function(x) x$delta_t)
  newStudyList$sampleSizes <- lapply(ctmaFitObject, function(x) x$sampleSize)
  newStudyList$moderators <- lapply(ctmaFitObject, function(x) x$moderators)
  newStudyList$source <- lapply(ctmaFitObject, function(x) x$source)
  newStudyList$pairwiseNs <- lapply(ctmaFitObject, function(x) x$pairwiseN)
  newStudyList$rawData <- lapply(ctmaFitObject, function(x) x$rawData)
  newStudyList$n.studies  <- length(newStudyList$rawData)
  newStudyList$summary <- newStudyList

  # CHD 13.6.2023
  if (reUseEmprawData == TRUE) newStudyList$emprawList <- lapply(tmp$emprawList, function(x) x[])

  class(newStudyList) <-  "CoTiMAFit"

  invisible(newStudyList)
}

