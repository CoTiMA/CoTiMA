##### CoTiMA fitting procedure - singleStudyModelFit

#' ctmaFittingSSMF
#'
#' @param coresToUse ?
#' @param empraw ?
#' @param currentModel ?
#' @param refits ?
#' @param retryattempts ?
#' @param singleModelStartValues ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaFittingSSMF <- function(coresToUse,
                              empraw,
                              currentModel,
                              refits,
                              retryattempts,
                              singleModelStartValues
                              )
{
  currentTime <- Sys.time()
  timeValue <- Sys.time() - currentTime; timeValue

  set.seed(as.numeric(timeValue))

  currentLabels <- c(currentModel$DRIFT, currentModel$DIFFUSION, currentModel$T0VAR); currentLabels
  currentLabels <- currentLabels[!(currentLabels %in% "0")]; currentLabels

 if (!(is.null(singleModelStartValues)))  {
    currentStartValues <- singleModelStartValues; currentStartValues
    names(currentStartValues) <- currentLabels
    results <- parallel::mclapply(seq(1, refits, by=1),
                        function(x) ctsem::ctFit(dat=empraw, currentModel, retryattempts=retryattempts,
                                          omxStartValues=currentStartValues),
                        mc.cores=coresToUse)
  } else {
    results <- parallel::mclapply(seq(1, refits, by=1),
                        function(x) ctsem::ctRefineTo(dat=empraw, currentModel, retryattempts=retryattempts),
                        mc.cores=coresToUse)
  }

  return(results)
}
