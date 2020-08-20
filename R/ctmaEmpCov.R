
# ctmaTargetMat changes a full covariance matrix: select target variables, recode them, combine them (add), and
#    add rows/columns with NA if focal concepts are not available.

# depends on mvrnorm, ctmapRaw, psych

ctmaEmpCov <- function(targetVariables=c(), recodeVariables=c(), combineVariables=c(),
                             combineVariablesNames=c(), missingVariables=c(),
                             nlatents=NULL, Tpoints=NULL, sampleSize=NULL, pairwiseN=NULL, empcov=NULL) {
  # select subset of variables
  empcov <- empcov[targetVariables, targetVariables]; empcov
  # recode variables (correlations)
  for (i in recodeVariables) {
    empcov[i, ] <- -1 * empcov[i, ]
    empcov[ ,i] <- -1 * empcov[ ,i]
    #print(i)
  }
  # combine variables
  if (!(is.null(sampleSize))) tmpData <- as.data.frame(mvrnorm(n=sampleSize, mu=rep(0, dim(empcov)[1]), Sigma=empcov, empirical = TRUE))
  if (!(is.null(pairwiseN))) {
    results <- pRaw(empCovMat=empcov, empNMat=pairwiseN, empN=NULL, studyNumber=1,
                             empMeanVector=NULL, empVarVector=NULL, activateRPB=FALSE)
    tmpData <- as.data.frame(results$data)
    names(tmpData) <- unlist(dimnames(empcov)[1])
  }

  if(length(combineVariables) > 0) {
    for (i in 1:length(combineVariables)) {
      tmpData$tmpName <- apply(tmpData[combineVariables[[i]]], 1, mean, na.rm=TRUE)
      tmp <- names(tmpData); tmp
      tmp[length(tmp)] <- combineVariablesNames[i]; tmp
      names(tmpData) <- tmp; names(tmpData)
    }
    data <- tmpData[, -c(1:(dim(empcov)[1]))]
  } else {
    data <- tmpData
  }

  tmp <- corr.test(data, ci=FALSE)
  empcovNew <- tmp$r
  pairwiseNNew  <- tmp$n
  # the var of combined variables != 0 => cov to cor
  dc <- diag(empcovNew)
  So <- sqrt(dc)
  ds <- diag(1/So)
  Mo <- ds %*% empcovNew %*% ds
  Mo <- Mo - 0.5*(Mo-t(Mo)) # Correct any round-off error
  empcovNew <- tmpMat <- Mo
  # impute missingness pattern
  numberOfMissingVars <- length(missingVariables); numberOfMissingVars
  if (numberOfMissingVars > 0) {
    presentVariables <- which(!((1:(Tpoints*nlatents)) %in% missingVariables)); presentVariables
    tmpMat <- matrix(NA, Tpoints*nlatents, Tpoints*nlatents)
    counter <- 1
    for (i in 1:(Tpoints*nlatents)) {
      if (i %in% missingVariables) {
        tmpMat[i, ] <- NA
        tmpMat[ ,i] <- NA
      } else {
        tmpMat[presentVariables[counter], presentVariables] <- empcovNew[counter, ]
        tmpMat[presentVariables ,presentVariables[counter]] <- empcovNew[ ,counter]
        counter <- counter + 1
      }
    }
  }
  results <- list()
  results$r <- tmpMat
  results$pairwiseN <- pairwiseNNew
  return(results)
}
