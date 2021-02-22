#' ctmaFitToPrep
#'
#'@description Extracts information from fitted CoTiMA object to (re-)crearte list of primary studies originally created with ctmaPrep
#'
#' @param ctmaFitObject  ctmaFitObject
#'
#' @examples
#' newStudyList <- ctmaFitToPrep(CoTiMAInitFit_3)
#'
#' @export ctmaFitToPrep
#'
ctmaFitToPrep <- function(ctmaFitObject=NULL)
{
  if (is.null(ctmaFitObject$originalStudyNo)) {
    if (!is.null(ctmaFitObject$studyList)) {
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

  class(newStudyList) <-  "CoTiMAFit"

  invisible(newStudyList)
}

