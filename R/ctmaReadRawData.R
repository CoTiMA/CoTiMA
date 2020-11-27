#debug
#listOfSpecs <- rawData11
# This function is created by adaptation of the respective code in ctmaInit
ctmaReadRawData <- function(listOfSpecs=list(fileName=NULL, studyNumbers=NULL, missingValues=c(-99),
                                             standardize=TRUE, header=FALSE, dec=".", sep=" ", n.latent=NULL)) {
  #listOfSpecs <- rawData11
  tmpData <- utils::read.table(file=listOfSpecs$fileName,
                               header=listOfSpecs$header,
                               dec=listOfSpecs$dec,
                               sep=listOfSpecs$sep)
  n.latent <- listOfSpecs$n.latent
  manifestNames <- paste0("V", 1:n.latent); manifestNames
  minInterval <- .00001

  # replace missing values
  tmpData <- as.matrix(tmpData) # important: line below will not work without having data as a matrix
  #tmpData[tmpData %in% studyList[[i]]$rawData$missingValues] <- NA
  tmpData[tmpData %in% listOfSpecs$missingValues] <- NA
  empraw <- as.data.frame(tmpData)
  ## START correction of current lags if entire time point is missing for a case
  # change variable names
  tmp1 <- dim(empraw)[2]; tmp1
  currentTpoints <- (tmp1 + 1)/(n.latent+1); currentTpoints
  colnames(empraw)[1:(currentTpoints * n.latent)] <- paste0(paste0("V", 1:n.latent), "_T", rep(0:(currentTpoints-1), each=n.latent))
  # wide to long
  emprawLongTmp <- ctsem::ctWideToLong(empraw, Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=manifestNames)
  emprawLongTmp <- suppressMessages(ctsem::ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))
  # eliminate rows where ALL latents are NA
  emprawLongTmp <- emprawLongTmp[apply(emprawLongTmp[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
  # eliminate rows where time is NA
  emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
  # make wide format
  emprawWide <- suppressMessages(ctsem::ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
  # intervalise
  emprawWide <- suppressMessages(ctsem::ctIntervalise(emprawWide,
                                                      Tpoints=currentTpoints,
                                                      n.manifest=n.latent,
                                                      manifestNames=manifestNames,
                                                      mininterval=minInterval))
  # restore
  empraw <- as.data.frame(emprawWide)
  # END correction
  #View(empraw)

  # Compute deltas if raw data are loaded
  #for (h in 1:(currentTpoints-1)) {
  #  #h <- 1
  #  # temporarily replace dT = minInterval (e.g., .0001) by NA
  #  tmp1 <- grep("dT", colnames(empraw)); tmp1
  #  colnamesTmp <- colnames(empraw)[tmp1]; colnamesTmp
  #  emprawTmp <- as.matrix(empraw[, tmp1])
  #  colnames(emprawTmp) <- colnamesTmp
  #  emprawTmp[emprawTmp == minInterval] <- NA
  #  #studyList[[i]]$delta_t[h] <- mean(emprawTmp[, paste0("dT", h)], na.rm=TRUE)
  #}
  #View(emprawTmp)

  # change sample size if entire cases were deleted
  #studyList[[i]]$sampleSize <- (dim(empraw))[1]
  #allSampleSizes[[i]] <- dim(empraw)[1]
  #currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
  #currentTpoints <- allTpoints[[i]]; currentTpoints
  currentVarnames <- c()
  for (j in 1:(currentTpoints)) {
    for (h in 1:n.latent) {
      currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
    }
  }
  #colnames(empraw)
  #colnames(empraw) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))

  # standardize (variables - not time lags) if option is chosen
  if (listOfSpecs$standardize == TRUE) empraw[, currentVarnames] <- scale(empraw[, currentVarnames])

  # replace missing values for time lags dTi by minInterval (has to be so because dTi are definition variables)
  tmpData <- empraw[, paste0("dT", seq(1:(currentTpoints-1)))]
  tmpData[is.na(tmpData)] <- minInterval
  empraw[, paste0("dT", seq(1:(currentTpoints-1)))] <- tmpData
  View(empraw)

  # add moderators to loaded raw data
  # Save raw data  on request
  #if ( i %in% saveRawData$studyNumbers ) {
  #  x1 <- paste0(saveRawData$fileName, i, ".dat"); x1
  #  utils::write.table(empraw, file=x1, row.names=saveRawData$row.names, col.names=saveRawData$col.names,
  #                     sep=saveRawData$sep, dec=saveRawData$dec)
  #}

  # augment pseudo raw data for stanct model
  {
    dataTmp <- empraw
    dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=currentTpoints, n.manifest=n.latent, #n.TIpred = (n.studies-1),
                                    manifestNames=manifestNames)
    dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
    #dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
    # eliminate rows where ALL latents are NA
    dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
    emprawLong <- dataTmp3

    listOfTimeLags <- dataTmp2[, "dT"]
    listOfTimeLags <- listOfTimeLags[listOfTimeLags != 0] # eliminate lag for Time 0
    table(listOfTimeLags)
  }

  View(dataTmp2)
} ### END for i ...
} ### END Read user provided data and create list with all study information ###
