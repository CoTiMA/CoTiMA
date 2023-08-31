#' ctmaEmpCov
#'
#' @description changes a full covariance matrix by selecting target variables, recoding them, combining them (compute the
#' mean of two or more variables), and by adding rows/columns with NA if focal variables are not available.
#'
#' @param targetVariables (col-/row-) number or names of the target variables
#' @param recodeVariables (col-/row-) number or names of the target variables require inverse coding
#' @param combineVariables list of vectors, which put together the targeted variables that should be used for composite variables
#' @param combineVariablesNames new names for combined variables - not really important
#' @param missingVariables missing variables
#' @param nlatents number of (latent) variables - actually it is the number of all variables
#' @param Tpoints number of time points.
#' @param sampleSize sample size
#' @param pairwiseN matrix of same dimensions as empcov containing possible pairwiseN.
#' @param empcov empirical correlation matrix
#'
#' @importFrom psych corr.test
#' @importFrom MASS mvrnorm
#' @importFrom crayon red blue
#'
#' @export ctmaEmpCov
#'
#' @examples
#' source17 <- c()
#' delta_t17 <- c(12)
#' sampleSize17 <- 440
#' empcov17 <- matrix(
#'   c( 1.00, -0.60, -0.36,  0.20,  0.62, -0.47, -0.18,  0.20,
#'     -0.60,  1.00,  0.55, -0.38, -0.43,  0.52,  0.27, -0.21,
#'     -0.36,  0.55,  1.00, -0.47, -0.26,  0.37,  0.51, -0.28,
#'      0.20, -0.38, -0.47,  1.00,  0.15, -0.28, -0.35,  0.56,
#'      0.62, -0.43, -0.26,  0.15,  1.00, -0.63, -0.30,  0.27,
#'     -0.47,  0.52,  0.37, -0.28, -0.63,  1.00,  0.55, -0.37,
#'     -0.18,  0.27,  0.51, -0.35, -0.30,  0.55,  1.00, -0.51,
#'      0.20, -0.21, -0.28,  0.56,  0.27, -0.37, -0.51,  1.00),
#'  nrow=8, ncol=8)
#' moderator17 <- c(3, 2)
#' rownames(empcov17) <- colnames(empcov17) <-
#'   c("Workload_1", "Exhaustion_1", "Cynicism_1", "Values_1",
#'     "Workload_2", "Exhaustion_2", "Cynicism_2", "Values_2")
#' targetVariables17 <-
#'   c("Workload_1", "Exhaustion_1", "Cynicism_1",
#'     "Workload_2", "Exhaustion_2", "Cynicism_2")
#' recodeVariables17 <- c("Workload_1", "Workload_2")
#' combineVariables17 <- list("Workload_1", c("Exhaustion_1", "Cynicism_1"),
#'                            "Workload_2", c("Exhaustion_2", "Cynicism_2"))
#' combineVariablesNames17 <- c("Demands_1",  "Burnout_1",
#'                              "Demands_2",  "Burnout_2")
#' missingVariables17 <- c();
#' results17 <- ctmaEmpCov(targetVariables = targetVariables17,
#'                         recodeVariables = recodeVariables17,
#'                         combineVariables = combineVariables17,
#'                         combineVariablesNames = combineVariablesNames17,
#'                         missingVariables = missingVariables17,
#'                         nlatents = 2, sampleSize = sampleSize17,
#'                         Tpoints = 2, empcov = empcov17)
#' empcov17 <- results17$r
#'
#' @return returns a list with two elements. The first element (results$r) contains the adapted correlation matrix, and
#' the second element (results$pairwiseNNew) an adapted version of a matrix of pairwise N if pariwiseN was provided for
#' the original correlation matrix supplied.
#'

ctmaEmpCov <- function(targetVariables=NULL, recodeVariables=c(), combineVariables=c(),
                             combineVariablesNames=c(), missingVariables=c(),
                             nlatents=NULL, Tpoints=NULL, sampleSize=NULL, pairwiseN=NULL, empcov=NULL) {

  # check
  if (is.null(sampleSize) & is.null(pairwiseN) ) {
    ErrorMsg <- "\nSample size (sampleSize) or pairwise N (pairwiseN) has to be provided! \nGood luck for the next try!"
    stop(ErrorMsg)
  }
  if (is.null(empcov)) {
    ErrorMsg <- "\nEmpirical correlation matrix has to be provided! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  # correct and then replace pairwise N (not required for data generation; only to provide corrected pairwise N after variables are combined)
  if (!(is.null(pairwiseN))) {
    tmp1 <- which(colnames(empcov) %in% targetVariables)
    pairwiseN <- pairwiseN[tmp1, tmp1]
    pairwiseNbackup <- pairwiseN
    if (is.null(sampleSize)) sampleSize <- max(pairwiseN) # max to avoid 0
    pairwiseN <- NULL
  } else {
    pairwiseNbackup <- matrix(sampleSize, dim(empcov)[1], dim(empcov)[2])
  }

  # select subset of variables
  if (!(is.null(targetVariables))) {
    tmp1 <- which(colnames(empcov) %in% targetVariables)
    empcov <- empcov[targetVariables, targetVariables]; empcov
  } else {
    tmp1 <- 1:dim(empcov)[1]
  }

  # recode variables (correlations)
  for (i in recodeVariables) {
    empcov[i, ] <- -1 * empcov[i, ]
    empcov[ ,i] <- -1 * empcov[ ,i]
  }

  # combine variables
  tmpData <- as.data.frame(MASS::mvrnorm(n=sampleSize, mu=rep(0, dim(empcov)[1]), Sigma=empcov, empirical = TRUE))
  newPairwiseN <- c()

  if(length(combineVariables) > 0) {
    for (i in 1:length(combineVariables)) {
      tmp1 <- which(targetVariables %in% combineVariables[[i]]); tmp1
      if (length(tmp1) > 1) {
        newPairwiseN[i] <- min(pairwiseNbackup[tmp1, tmp1]); newPairwiseN[i]
        tmp2 <- pairwiseNbackup; tmp2
        diag(tmp2) <- 99999999
        lowerPart <- 1:(min(tmp1)-1); lowerPart
        upperPart <- (max(tmp1)+1):(dim(pairwiseNbackup)[1]); upperPart
        if (min(lowerPart) >= 1) {
          tmp3lower <- min(tmp2[tmp1, lowerPart]); tmp3lower
          pairwiseNbackup[lowerPart, tmp1] <- tmp3lower
        pairwiseNbackup[tmp1, lowerPart] <- tmp3lower
        }
        if (max(upperPart) <= (dim(empcov)[1])) {
          tmp3upper <- min(tmp2[tmp1, upperPart]); tmp3upper
          pairwiseNbackup[upperPart, tmp1] <- tmp3upper
        pairwiseNbackup[tmp1, upperPart] <- tmp3upper
        }
        pairwiseNbackup[tmp1, tmp1] <- newPairwiseN[i]
      }
      tmpData$tmpName <- apply(tmpData[combineVariables[[i]]], 1, mean, na.rm=TRUE)
      tmp <- names(tmpData); tmp
      tmp[length(tmp)] <- combineVariablesNames[i]; tmp
      names(tmpData) <- tmp; names(tmpData)
    }
    data <- tmpData[, -c(1:(dim(empcov)[1]))]
    tmp1 <- unlist(lapply(combineVariables,function(extract) length(extract))); tmp1
    targetCols <- c()
    for (j in 1:(length(tmp1))) targetCols[j] <- sum(tmp1[1:j])
    pairwiseN <- pairwiseNbackup[targetCols , targetCols]
  } else {
    data <- tmpData
  }

  tmp <- psych::corr.test(data, ci=FALSE)
  empcovNew <- tmp$r
  pairwiseNNew  <- tmp$n
  if (exists("pairwiseNbackup")) pairwiseNNew <- pairwiseNbackup

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
    tmpMat <- tmpMatN <- matrix(NA, Tpoints*nlatents, Tpoints*nlatents)
    counter <- 1
    for (i in 1:(Tpoints*nlatents)) {
      if (i %in% missingVariables) {
        tmpMat[i, ] <- NA
        tmpMat[ ,i] <- NA
        if (exists("pairwiseN")) {
          tmpMatN[i, ] <- 0
          tmpMatN[, i] <- 0
        }
      } else {
        tmpMat[presentVariables[counter], presentVariables] <- empcovNew[counter, ]
        tmpMat[presentVariables ,presentVariables[counter]] <- empcovNew[ ,counter]
        if (exists("pairwiseN")) {
          tmpMatN[presentVariables[counter], presentVariables] <- pairwiseN[counter, ]
          tmpMatN[presentVariables ,presentVariables[counter]] <- pairwiseN[ ,counter]
          }
        counter <- counter + 1
      }
    }
    pairwiseNNew <- tmpMatN
  }
  results <- list()
  results$r <- tmpMat
  results$rNew <- tmpMat   # CHD added 31. Auf 2023 (to be consistent with pairwiseNNew)
  results$pairwiseNNew <- pairwiseNNew
  return(results)
}
