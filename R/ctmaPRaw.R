#' ctmaPRaw
#'
#' @description Converts empirical correlation matrices to pseudo raw data (i.e. random data, that perfectly reproduce the correlations)
#'
#' @param empCovMat empirical primary study covariance matrix
#' @param empNMat matrix of (possibly pairwise) N
#' @param empN N (in case of listwise N)
#' @param studyNumber internal number
#' @param empMeanVector vector of means for all variables, usually 0
#' @param empVarVector vector of variances for all variables, usually 1
#' @param activateRPB set TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param experimental set TRUE to try new pairwise N function
#'
#' @importFrom RPushbullet pbPost
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#' @importFrom psych corr.test
#'
#'
ctmaPRaw <- function(empCovMat=NULL, empNMat=matrix(0,0,0), empN=NULL, studyNumber=NULL,
                     empMeanVector=NULL, empVarVector=NULL, activateRPB=FALSE,
                     experimental=FALSE)
{  # begin function definition (until end of file)

  if (is.null(empCovMat)) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nNo empirical covariance matrix provided for pseudo raw data generation! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (!(isSymmetric(empCovMat))) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nThe empirical covariance matrix provided is not symmetric! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (!(is.null(empNMat))) {
    if (!(isSymmetric(empNMat))) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nThe pairwise N matrix provided is not symmetrix! \nGood luck for the next try!"
      stop(ErrorMsg)
    }
  }

  if ( (is.null(empNMat) & is.null(empN)) ) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nEITHER a matrix with pairwise N OR an overall N has to be provided pseudo raw data generation! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if ( (!(is.null(empMeanVector))) & (length(empMeanVector) != (dim(empCovMat)[1]) )  ){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nThe number of means provided does not match the number of variables in the empirical covariance matrix! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  rowNACounter <- colNACounter <- c()
  for (i in 1:(dim(empCovMat)[1])) {
    rowNACounter[i] <- length(which(is.na(empCovMat[i, ]) == TRUE))
    colNACounter[i] <- length(which(is.na(empCovMat[ ,i]) == TRUE))
  }
  if (any(rowNACounter != colNACounter)) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nCurrently missing correlations can only be handled if a variable is entirely missing. \nThe NA-pattern provided implies this is not the case. \nConsider setting all correlations involving a critical variable to NA. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  # Define objects
  if (dim(empNMat)[1] == 0) currentN <- matrix(empN, dim(empCovMat)[1], dim(empCovMat)[1]) else  currentN <- empNMat
  currentR <- empCovMat; currentR
  if (is.null(studyNumber)) studyNumber <- 1

  ###### Handling of missing correlations (NA) ######
  origN <- currentN; origN
  eigenValueProblems <- "Congratulations! This covariance matrix is positive definite (no negative eigenvalues)."
  newMat <- currentR; newMat
  k <- which(is.na(newMat), arr.ind=TRUE); k
  if (dim(k)[1] > 0) {
    sortedK <-  t(apply(k,1,sort)); sortedK
    uniqueK <- unique(sortedK); uniqueK
    reversedK <- cbind(uniqueK[,2], uniqueK[,1]); reversedK
    k <- rbind(uniqueK, reversedK); k
    randomCors <- stats::rnorm(length(uniqueK)/2, 0, .0); randomCors
    randomCors <- c(randomCors, randomCors); randomCors
    newMat[k] <- randomCors
    diag(newMat) <- 1; newMat
  }

  ### solve problem with non-positive definite matrices (https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/)
  #counter <- 0
  ##eigen(newMat)$values
  #while(any(eigen(newMat)$values < 0)) {
  #  counter <- counter + 1
  #  print(cat("Study", studyNumber, "has negative eigenvalues. Try to fix this problem at iteration =", counter, "!"))
  #  newEig <- eigen(newMat); newEig
  #  newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
  #  newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
  #  newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
  #  eigenValueProblems <- paste0("The correlation matrix of Study ", i, "had negative Eigenvalues. I tried to fix this, but better check the matrix.")
  #}

  toBeMadeNA <- makeNAdueToNACors <- makeNAdueToZeroN <- NULL
  # store rows and cols with correlations == NA
  if (dim(k)[1] > 0) makeNAdueToNACors <- which(table(c(k[,2])) > dim(newMat)[2]/2); makeNAdueToNACors # which columns in the final data set should be made NA
  makeNAdueToNACors <- as.numeric(names(makeNAdueToNACors)); makeNAdueToNACors ##### CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # store rows and cols with N == 0
  makeNAdueToZeroN <- which((apply(currentN, 1, sum) == 0) & (apply(currentN, 2, sum) == 0) == TRUE); makeNAdueToZeroN
  toBeMadeNA <- unique(c(makeNAdueToNACors, makeNAdueToZeroN)); toBeMadeNA

  currentR <- newMat; currentR
  currentN[currentN == 0] <- max(currentN[currentN > 0]); currentN #replace N=0 (missing correlations) by max N

  ###### Pairwise deletion ######
  if (experimental == FALSE) {
    numberOfMatrices <- 2^(dim(currentN)[1]); numberOfMatrices
    if (numberOfMatrices > 10^6) numberOfMatrices <- 10^6     # prevent memory overflow
    matrixL <- lapply(1:numberOfMatrices, list); matrixL
    matrixL[[1]] <- list(r=currentR, n=currentN); matrixL[[1]]
    dimnames(matrixL[[1]]$r) <- dimnames(matrixL[[1]]$n) <- list(1:(dim(matrixL[[1]]$r)[1]), 1:(dim(matrixL[[1]]$n)[1]))
    whileCounter <- matCounter <- 1
    ncounter <- 0
    while(length(matrixL[[whileCounter]]) > 1) {
      minN <- min(matrixL[[whileCounter]]$n); minN
      matDim <- dim(matrixL[[whileCounter]]$n)[1]; matDim
      if (is.null(matDim)) matDim <- 1; matDim
      if ( (minN == 1 & matDim == 1) | (minN == 0) ) {
        matrixL[[whileCounter]]$data <- (matrix(NA, 1, matDim))
      } else {
        if (matDim > 1) minNcoord <- which(min(matrixL[[whileCounter]]$n) == matrixL[[whileCounter]]$n, arr.ind=TRUE)[1, ] else minNcoord <- c(1,1)
        # correction for minN < dimMat by deleting rows cols with N < dimMat
        while (minN <= matDim) {
          # retain submatrix with largest sum of N
          tmp1 <- matrixL[[whileCounter]]$n[-minNcoord[1], -minNcoord[1]]; tmp1
          tmp2 <- matrixL[[whileCounter]]$n[-minNcoord[2], -minNcoord[2]]; tmp2
          tmp1 <- sum(tmp1[lower.tri(tmp1)]); tmp1
          tmp2 <- sum(tmp2[lower.tri(tmp2)]); tmp2
          if (tmp1 >= tmp2) {
            matrixL[[whileCounter]]$r <- matrixL[[whileCounter]]$r[-minNcoord[1], -minNcoord[1]]; matrixL[[whileCounter]]$r
            matrixL[[whileCounter]]$n <- matrixL[[whileCounter]]$n[-minNcoord[1], -minNcoord[1]]; matrixL[[whileCounter]]$n
          } else {
            matrixL[[whileCounter]]$r <- matrixL[[whileCounter]]$r[-minNcoord[2], -minNcoord[2]]; matrixL[[whileCounter]]$r
            matrixL[[whileCounter]]$n <- matrixL[[whileCounter]]$n[-minNcoord[2], -minNcoord[2]]; matrixL[[whileCounter]]$n
          }
          minN <- min(matrixL[[whileCounter]]$n); minN
          if (is.null(dim(matrixL[[whileCounter]]$n)[1])) matDim <- 1 else matDim <- dim(matrixL[[whileCounter]]$n)[1]
          if (matDim > 1) minNcoord <- which(min(matrixL[[whileCounter]]$n) == matrixL[[whileCounter]]$n, arr.ind=TRUE)[1, ] else minNcoord <- c(1,1)
          if (minN == 1 & matDim == 1) break #leave while loop
        }
        if (!(minN == 1 & matDim == 1)) {
          r <- matrixL[[whileCounter]]$r; r
          n <- matrixL[[whileCounter]]$n; n
          # generate data
          matrixL[[whileCounter]]$data <- MASS::mvrnorm(n=minN, mu=rep(0, matDim), Sigma=r, empirical = TRUE); matrixL[[whileCounter]]$data
          ncounter <- ncounter + minN
        }
        if (matDim > 1) {
          matrixL[[matCounter + 1]] <- matrixL[[matCounter + 2]] <- list()
          # reduce correlation matrices
          matrixL[[matCounter + 1]]$r <- r[-minNcoord[1], -minNcoord[1]]; matrixL[[matCounter + 1]]$r
          matrixL[[matCounter + 2]]$r <- r[-minNcoord[2], -minNcoord[2]]; matrixL[[matCounter + 2]]$r
          # preserve names and make matrix if dim = NUL (only 1 element)
          if (is.null(dim(matrixL[[matCounter + 1]]$r))) {
            matrixL[[matCounter + 1]]$r <- matrix(matrixL[[matCounter + 1]]$r, 1, 1,
                                                  dimnames=list(unlist(dimnames(r)[1])[-minNcoord[1]], unlist(dimnames(r)[1])[-minNcoord[1]]))
          }
          if (is.null(dim(matrixL[[matCounter + 2]]$r))) {
            matrixL[[matCounter + 2]]$r <- matrix(matrixL[[matCounter + 2]]$r, 1, 1,
                                                  dimnames=list(unlist(dimnames(r)[1])[-minNcoord[2]], unlist(dimnames(r)[1])[-minNcoord[2]]))
          }
          ## reduce matrices with pairwise sample size and reduce sample sizes by n already created
          matrixL[[matCounter + 1]]$n <- n[-minNcoord[1], -minNcoord[1]] - minN; matrixL[[matCounter + 1]]$n
          matrixL[[matCounter + 2]]$n <- n[-minNcoord[2], -minNcoord[2]] - minN; matrixL[[matCounter + 2]]$n
          # preserve names
          if (is.null(dim(matrixL[[matCounter + 1]]$n))) {
            matrixL[[matCounter + 1]]$n <- matrix(matrixL[[matCounter + 1]]$n, 1, 1,
                                                  dimnames=list(unlist(dimnames(r)[1])[-minNcoord[1]], unlist(dimnames(r)[1])[-minNcoord[1]]))
          }
          if (is.null(dim(matrixL[[matCounter + 2]]$n))) {
            matrixL[[matCounter + 2]]$n <- matrix(matrixL[[matCounter + 2]]$n, 1, 1,
                                                  dimnames=list(unlist(dimnames(r)[1])[-minNcoord[2]], unlist(dimnames(r)[1])[-minNcoord[2]]))
          }
          # variables included in both matrices
          targets <- (colnames(matrixL[[matCounter + 1]]$n)[colnames(matrixL[[matCounter + 1]]$n)
                                                            %in% colnames(matrixL[[matCounter + 2]]$n)]);targets
          # take only half of n in overlapping variables
          matrixL[[matCounter + 1]]$n[targets, targets] <- round(matrixL[[matCounter + 1]]$n[targets, targets]/2 +.01); matrixL[[matCounter + 1]]$n
          matrixL[[matCounter + 2]]$n[targets, targets] <- round(matrixL[[matCounter + 2]]$n[targets, targets]/2 -.01); matrixL[[matCounter + 2]]$n
          matCounter <- matCounter + 2; matCounter
        }
      }
      whileCounter <- whileCounter + 1; whileCounter
    }

    # combine data
    newData <- as.data.frame(matrix(NA, ncounter, ncol=(dim(currentR)[1]) ) ); newData
    colnames(newData) <- 1:(dim(currentR)[1]); colnames(newData)
    rowCounter <- 1
    for (i in 1:(whileCounter-1)){
      if (!(is.null(dim(matrixL[[i]]$data)[1]))) tmpN <- dim(matrixL[[i]]$data)[1] else tmpN <- 1
      currentRows <- rowCounter:(rowCounter+tmpN-1); currentRows
      currentNames <- colnames(matrixL[[i]]$data); currentNames
      newData[currentRows, currentNames] <- matrixL[[i]]$data
      rowCounter <- max(currentRows) + 1; rowCounter
    }

    ### too small N in the diagonal can be corrected
    nDiff <- psych::corr.test(newData, ci=FALSE)$n - origN; nDiff
    rDiff <- round(psych::corr.test(newData, ci=FALSE)$r - currentR, 6); rDiff
    currentData <- list()
    counter <- 0
    for (i in 1:(dim(nDiff)[1])) {
      if (nDiff[i,i] <= -1) {
        counter <- counter +1
        currentData[[counter]] <- MASS::mvrnorm(n=-nDiff[i,i], mu=rep(0), Sigma=matrix(currentR[i,i], 1, 1), empirical = TRUE)
        #colnames(currentData[[counter]]) <- i; currentData[[counter]]
      } else {
        counter <- counter +1
        currentData[[counter]] <- matrix(NA, 1, 1)
      }
    }
    # add data
    if (counter > 0) {
      rowCounter <- dim(newData)[1]+1; rowCounter
      for (i in 1:counter) {
        if (!(is.null(dim(currentData[[i]])[1]))) tmpN <- dim(currentData[[i]])[1] else tmpN <- 1
        currentRows <- rowCounter:(rowCounter+tmpN-1); currentRows
        #currentNames <- colnames(currentData[[i]]); currentNames
        tmpMat <- as.data.frame(matrix(NA, tmpN, dim(newData)[2])); tmpMat
        colnames(tmpMat) <- colnames(newData)
        newData <- rbind(newData, tmpMat)
        #newData[currentRows, currentNames] <- currentData[[i]]
        newData[currentRows, i] <- currentData[[i]] # inserted
        rowCounter <- max(currentRows) + 1; rowCounter
      }
    }
    if (all(diag(currentR) == 1)) newData <- scale(newData)

  } # END if (experimental == FALSE)

  if (experimental == TRUE)  {
    tmpRMat <- currentR; tmpRMat
    tmpNMat <- currentN
    tmpNMat[upper.tri(tmpNMat)] <- NA; tmpNMat
    colnames(tmpRMat) <- rownames(tmpRMat) <- seq(1, dim(tmpRMat)[1], 1)

    # d
    newData <- matrix(NA, ncol=dim(tmpRMat)[1], nrow=max(tmpNMat, na.rm=T)); newData  # all variables have 0
    colnames(newData) <- tmpColNames <- colnames(tmpRMat)

    #
    counter <- 0
    currentN2 <- 0          # currentN already in use in ctmaPRaw
    currentStartCol <- 1 # where start inserting praw data

    # try making positive definite by adding .01 to diag if necessary
    while (any(eigen(tmpRMat)$values < 0)) {
      print("Adding .01 to diagonal of correlation matrix to make it positive definite (ridge constant)")
      tmpRMat <- tmpRMat + diag(.01, dim(tmpRMat)[1])
    }

    allCollectors <- list()
    while (any(tmpNMat > 1, na.rm=T) & dim(tmpRMat)[1] >1 & dim(tmpNMat)[1] > 1) {
      counter <- counter + 1; counter
      collectorCounter <- 0
      min1 <- unique(which(tmpNMat == min(tmpNMat, na.rm=T), arr.ind = TRUE)); min1      # which matrix index is currenlty min(N)
      currentN2 <- tmpNMat[min1[1,1], min1[1,2]]; currentN2                                    # what is the current min(N)

      if (counter > 1) {
        for (i in length(allCollectors):1) {
          if (allCollectors[[i]][1] == length(allCollectors[[i]])-1) {
            tmpRMat <- tmpRMat[-allCollectors[[i]][1], -allCollectors[[i]][1]]
            tmpColNames <- tmpColNames[tmpColNames != allCollectors[[i]][1]]; tmpColNames
          } else {
            remove <- c()
            if (allCollectors[[i]][1] == allCollectors[[i]][length(allCollectors[[i]])]) remove <- allCollectors[[i]][1]
            tmpRMat <- tmpRMat[!(rownames(tmpRMat) %in% remove), !(colnames(tmpRMat) %in% remove)]
            tmpColNames <- tmpColNames[!(tmpColNames %in% remove)]; tmpColNames
          }
        }
      }
      if (is.null(dim(tmpRMat))) tmpRMat <- matrix(tmpRMat, 1, 1)
      if (currentN2 >= dim(tmpRMat)[1]) {                                                  # if not, some cases are lost
        data <- MASS::mvrnorm(n=currentN2, mu=rep(0, dim(tmpRMat)[1]),
                              Sigma=tmpRMat, empirical = TRUE)                            # create praw
        if (dim(newData)[1] < currentStartCol+dim(data)[1]) {
          tmpData <- newData[1:dim(data)[1] , ]
          tmpData[!(is.null(tmpData))] <- NA
          newData <- rbind(newData, tmpData)
        }
        newData[currentStartCol:(currentStartCol+dim(data)[1]-1) , tmpColNames] <-  data  # insert in in newData
      }

      currentStartCol <- currentStartCol + currentN2; currentStartCol

      # correction
      psych::corr.test(newData)
      tmpNMat[!(is.na(tmpNMat))] <- tmpNMat[!(is.na(tmpNMat))] - currentN2
      tmpNMat[tmpNMat == 0] <- NA

      ## collect fully filled rows and columns (per row)
      collector <- list()
      counter2 <- 0
      # collect
      for (r in unique(min1[ ,1])) {
        counter2 <- counter2 + 1
        collector[[counter2]] <- min1[min1[,1]==r , ]
      }
      # create vector in which the first value is a critical row and all subsequent ones are the forbidden columns
      allCollectors <- list() # collects lists of variable sets for which the remaining N is 0 after previous data computation
      for (r in 1:length(collector)) {
        collectorCounter <- collectorCounter + 1
        if (is.null(dim(collector[[r]]))) collector[[r]] <- matrix(collector[[r]], ncol=2, nrow=1)
        allCollectors[[collectorCounter]] <- c(collector[[r]][1,1], collector[[r]][ ,2])
      }
    }
    if (all(diag(currentR) == 1)) newData <- scale(newData)
  } # END if (experimental == TRUE)


  # replace values which cannot exist
  if (!(is.null("toBeMadeNA"))) newData[,toBeMadeNA] <- NA


  # add means if they are provided
  if (!(is.null(empMeanVector))) {
    for (i in 1:(dim(newData)[2])) {
      newData[,i] <- newData[,i] + empMeanVector[i]
    }
  }

  nDiff <- psych::corr.test(newData, ci=FALSE)$n - origN; nDiff
  rDiff <- round(psych::corr.test(newData, ci=FALSE)$r - currentR, 6); rDiff
  overallNDiff <- sum(nDiff[lower.tri(nDiff, diag=TRUE)]); overallNDiff
  relativeNDiff <- overallNDiff/sum(origN[lower.tri(origN, diag=TRUE)]); relativeNDiff
  overallRDiff <- sum(rDiff[lower.tri(rDiff, diag=TRUE)], na.rm=TRUE); overallRDiff
  results <- list(newData, eigenValueProblems, nDiff, overallNDiff, relativeNDiff)
  names(results) <- c("data", "problems", "lostN", "overallLostN", "relativeLostN")

  invisible(results)
}
