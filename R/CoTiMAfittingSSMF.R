##### CoTiMA fitting procedure - singleStudyModelFit

#' CoTiMAfittingSSMF
#'
#' @param coresToUse
#' @param empraw
#' @param currentModel
#' @param refits
#' @param retryattempts
#' @param singleModelStartValues
#'
#' @return
#' @export
#'
#' @examples
CoTiMAfittingSSMF <- function(coresToUse,
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
    results <- mclapply(seq(1, refits, by=1),
                        function(x) ctFit(dat=empraw, currentModel, retryattempts=retryattempts,
                                          omxStartValues=currentStartValues),
                        mc.cores=coresToUse)
  } else {
    results <- mclapply(seq(1, refits, by=1),
                        function(x) ctRefineTo(dat=empraw, currentModel, retryattempts=retryattempts),
                        mc.cores=coresToUse)
  }

  return(results)
}
