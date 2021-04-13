#' ctmaStanResample
#'
#' @description re-sample from a fitted stanct model to achieve desired number of finishsamples (could be useful to prevent exhausted memory)
#'
#' @param ctmaFittedModel a 'CoTiMA' fit object, usually with few 'finishsamples' to prevent memory exhaustion
#' @param nsamples sample size per run
#' @param overallSamples overall samples size to be achieved
#'
#' @importFrom ctsem ctAddSamples
#' @importFrom abind abind
#'
#' @return returns a CoTiMA fit object with an increased number of finish samples
#'
ctmaStanResample <- function(ctmaFittedModel=NULL,
                             nsamples=25,
                             overallSamples=500) {

  sm1 <- ctmaFittedModel
  sm1$stanfit$rawposterior <- sm1$stanfit$rawposterior[1:2,]
  runs <- round(overallSamples/nsamples+0.5); runs
  currentResults <- list()
  for (i in 1:runs) {
    cat(i*nsamples, "out of", overallSamples, "\n")
    currentResults[[i]] <- ctsem::ctAddSamples(sm1, nsamples=nsamples, cores=1) # ctmaAddSamples
  }
  ## combine resampled results
  # extract 2 rawposteriors
  tmp <- ctmaFittedModel$stanfit$rawposterior[ , ]; tmp
  for (i in 1:length(currentResults)) tmp <- rbind(tmp, currentResults[[i]]$stanfit$rawposterior[-c(1,2), ])
  ctmaFittedModel$stanfit$rawposterior <- tmp
  # 4 different object types to store transformedpars (vector, matrix, 2-dim array, 4-dim array)
  # determine elements, dims and names automatically to prevent failure after C. Driver changes list
  currenNames <- currentDims <- c()
  for (i in 1:length(ctmaFittedModel$stanfit$transformedpars)) {
    tmp <- (ctmaFittedModel$stanfit$transformedpars[i])
    currenNames[i] <- names(tmp); currenNames
    currentDims[i] <- length(dim(tmp[[1]])); currentDims
  }
  n.diff.elements <- max(currentDims) - min(currentDims) +  1; n.diff.elements
  targetElementXdim <- list()
  targetElementXdim[[n.diff.elements+1]] <- ""
  for (i in 1:length(ctmaFittedModel$stanfit$transformedpars)) {
    targetElementXdim[[currentDims[i]]] <- append(targetElementXdim[[currentDims[i]]],
                                                  currenNames[i])
  }
  # combine results of all runs
  for (h in 1:n.diff.elements) {
    for (j in 1:length((targetElementXdim[[h]]))) {
      tmp2 <- ctmaFittedModel$stanfit$transformedpars[[ targetElementXdim[[h]][j] ]]
      for (i in 1:length(currentResults)) {
        tmp3 <- currentResults[[i]]$stanfit$transformedpars[[ targetElementXdim[[h]][j] ]]#[-c(1,2)]
        tmp4 <- CoTiMAStanctArgs$optimcontrol$finishsamples; tmp4
        tmp5 <- dim(tmp3)
        tmp5[1] <- tmp5[1] - tmp4; tmp5
        if(length(dim(tmp3)) == 1) tmp3 <- tmp3[-c(1:tmp4)]
        if(length(dim(tmp3)) == 2) tmp3 <- tmp3[-c(1:tmp4), ]
        if(length(dim(tmp3)) == 3) tmp3 <- tmp3[-c(1:tmp4), , ]
        if(length(dim(tmp3)) == 4) tmp3 <- tmp3[-c(1:tmp4), , , ]
        # add empty dimension that could be lost (was emppty) in last step before
        if (length(dim(tmp2)) != length(dim(tmp3))) tmp3 <- array(tmp3, dim=tmp5)
        if (i > 1) tmp2 <- abind::abind(tmp2, tmp3, along=1) else tmp2 <- tmp3
      }
      ctmaFittedModel$stanfit$transformedpars[[ targetElementXdim[[h]][j] ]] <- tmp2
    }
  }
  return(ctmaFittedModel)

}

