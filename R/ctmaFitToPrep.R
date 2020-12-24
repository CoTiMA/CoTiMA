ctmaFitToPrep <- function(ctmaFitObject=NULL)
{
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

  return(newStudyList)
}
