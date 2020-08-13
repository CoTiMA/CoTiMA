#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

ctmaInit <- function(
  # Primary Study Information
  primaryStudies=NULL,                    #NEW: list of lists: list(deltas, sampleSizes, empcovs, moderators, startValues, studyNumbers)

  # Directory names and file names
  activeDirectory=NULL,
  resultsFilePrefix="ctmaInit",           # the prefix for the result file that is created
  saveFilePrefix="ctmaInit",
  loadFilePrefix="ctmaInit",


  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,                      #set to TRUE to receive push messages with CoTiMA notifications on your phone
  checkSingleStudyResults=TRUE,          # displays estimates from single study ctsem models and waits for user inoput to continue
  digits=4,

  # General Model Setup
  n.latent=NULL,

  # Load/save Raw Data
  loadRawData=list(),
  saveRawData=list(),

  # Fitting Parameters
  retryattempts=30,                       # number of iterations to be used by ctsem
  coresToUse=c(1),                        # could be set to -n (usually -1) which is subtracted from the overall number available (on non Windows machines)

  ## Save computed model fits or load previous fits (then provide (old) file prefix)
  silentOverwrite=FALSE,

  # ctsem models for all primary studies (always required for getting good starting values (and model fit for heterogeneity model without fitting it))
  saveSingleStudyModelFit=c(),            # save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
  loadSingleStudyModelFit=c(),            # load the fit of single study ctsem models

  CoTiMAStanctArgs=list(test=TRUE,
                        scaleTI=TRUE, scaleMod=TRUE, scaleLongData=FALSE,
                        scaleTime=1/1,
                        savesubjectmatrices=FALSE, verbose=1,
                        datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
                        iter=1000, intoverstates=TRUE,
                        binomial=FALSE, fit=TRUE,
                        intoverpop=FALSE, stationary=FALSE,
                        plot=FALSE, derrind="all",
                        optimize=TRUE, optimcontrol=list(is=F, stochastic=FALSE, finishsamples=1000),
                        nlcontrol=list(),
                        nopriors=TRUE,
                        chains=2,
                        cores=1,
                        inits=NULL, forcerecompile=FALSE,
                        savescores=FALSE, gendata=FALSE,
                        control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
                        verbose=0)

)

{  # begin function definition (until end of file)

  options(scipen = 1)


  #######################################################################################################################
  ########################################### Check Model Specification #################################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("########################## Check Model Specification ############################"))
    print(paste0("#################################################################################"))

    if (is.null(primaryStudies)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold(" List with lists of primary study information not specified!", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat(red$bold("Should I try to get primary study information (e.g., empcov1, delta_t22) from the global environment?", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat(red$bold("(This is the way how earlier version of CoTiMA were run but no longer recommended.)", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat(blue("Press 'q' to quit and specify or 'c' to give it a try. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and specify primaryStudies or 'c' to try reading from global environment Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') {
        stop("Good luck for the next try!")
      } else {
        if (!(is.numeric(tryCatch(get(paste("sampleSize", 1, sep = "")), error = function(e) e)))) {
          cat(red$bold("Getting primary study information from the global environment failed", sep="\n"))
          cat(red$bold(" ", " ", sep="\n"))
          cat(blue("To test, I searched for sampleSize1 and could not find it!", "\n"))
          cat(red$bold(" ", " ", sep="\n"))
          stop("Good luck for the next try!")
        } else {
          cat(blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n"))
          char <- as.numeric(readline(""))
          while ((is.na(char))) {
            cat((blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n")))
            char <- as.numeric(readline(""))
          }
          primaryStudies <- prep(selectedStudies=1:char)
        } # END else (catch failed)
      } # END else
    } #


    if (is.null(n.latent)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("Number of variables (n.latent) not specified!", sep="\n"))
      stop("Good luck for the next try!")
    }


    if (is.null(n.latent)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("Number of variables (n.latent) not specified!", sep="\n"))
      stop("Good luck for the next try!")
    }

    if (is.null(activeDirectory)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("No working directory has been specified!", sep="\n"))
      stop("Good luck for the next try!")
    }

    if (resultsFilePrefix=="ctmaInit") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default results file prefix (ctmaInit) has been chosen.", "\n"))
      cat(blue("Press 'q' to quit and change or'c'to continue. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }
    rsultsFileName <- paste0(resultsFilePrefix, ".txt")

    if (saveFilePrefix=="ctmaInit") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default save file prefix (ctmaInit) has been chosen.", "\n"))
      cat(blue("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards "))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change filename or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    resultsFileName <-paste0(resultsFilePrefix, ".txt")

    if (is.null(loadFilePrefix)) {ncharLoadFilePrefix <- 0} else {ncharLoadFilePrefix <- nchar(loadFilePrefix)} # check maximum path + filename length
    ncharResultsFileName <- nchar(resultsFileName)
    ncharFilePrefix <- max(nchar(saveFilePrefix), nchar(loadFilePrefix))
    ncharactiveDirectory <- nchar(activeDirectory)
    ncharModelFits <- max( (nchar(saveSingleStudyModelFit[1])+2),
                           nchar(loadSingleStudyModelFit),
                           #nchar(saveDRIFTAllModelFit),
                           #nchar(loadDRIFTAllModelFit),
                           #nchar(saveDRIFTSingleModelFit),
                           #nchar(loadDRIFTSingleModelFit),
                           #nchar(saveDRIFTCrossModelFit),
                           #nchar(loadDRIFTCrossModelFit),
                           #nchar(saveDRIFTAutoModelFit),
                           #nchar(loadDRIFTAutoModelFit),
                           #nchar(saveHeterogeneityModelFit),
                           #nchar(loadHeterogeneityModelFit),
                           1) # 1 added to prevent empty vector
    if (nchar("_CoTiMA results for moderated auto-regressive (all effects in the model were moderated) effects of V2.png") +
        max(ncharLoadFilePrefix, ncharResultsFileName, ncharFilePrefix, ncharModelFits) +
        ncharactiveDirectory > 259) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The overall length of the path name (working directory) + the file names to be created exceeds 260.", "\n"))
      cat(blue("Shorten the file names or prefixes, or better work in a different working directory.", sep="\n"))
      stop("Good luck for the next try!")
    }

    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  {
        if (.Platform$OS.type == "unix") {
          coresToUse <- parallel::detectCores() + coresToUse
          #coresToUse <- detectCores() + coresToUse
          cat("     #################################################################################\n")
          cat("     # You decided using multiple CPU cores for parallel fitting CoTiMA models. This #\n")
          cat("     # may work well in many cases, but we highly recommend not to run this script   #\n")
          cat("     # in RStudio. Instead, go to the console (e.g., console.app on MAC OS), then    #\n")
          cat("     # change to the direcrtory where your R-script is located (by typing, e.g.,     #\n")
          cat("     # \'cd \"/Users/trump/only/temp\"\'), and then, e.g., \'Rscript \"CoTiMA1.R\"\'.#\n")
          cat("     #################################################################################")
          Sys.sleep(1)
        } else {
          coresToUse <- 1
        }
      }
    } else {
      coresToUse <- 1
    }

    if ((-1) * coresToUse >= parallel::detectCores()) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- parallel::detectCores() - 1
      cat(red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
    }

    numOfThreads <- round((coresToUse), 0); numOfThreads
    if(numOfThreads < 1) numOfThreads <- 1
    tmp <- paste0("No. of Threads set to ", numOfThreads); tmp
    cat(red(tmp,"\n"))

  } ### END Check Model Specification ###


  #######################################################################################################################
  ##################### Read user provided data and create list with all study information ##############################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("###### Read user provided data and create list with all study information #######"))
    print(paste0("#################################################################################"))

    #setwd(activeDirectory)
    start.time <- Sys.time(); start.time

    studyList <- list()
    minInterval <- .0001
    maxLengthModeratorVector <- 0
    loadRawDataStudyNumbers <- unlist(lapply(primaryStudies$rawData,
                                             function(extract) extract$studyNumbers)); loadRawDataStudyNumbers

    # resorting primary study information
    if (!(exists("moderatorNumber"))) moderatorNumber <- 1; moderatorNumber
    for (i in 1:length(unlist(primaryStudies$studyNumbers))) {
      studyList[[i]] <- list(studyNumber=i, empcov=primaryStudies$empcovs[[i]], delta_t=primaryStudies$deltas[[i]],
                             sampleSize=primaryStudies$sampleSizes[[i]], originalStudyNo=primaryStudies$studyNumber[[i]],
                             timePoints=sum(length(primaryStudies$deltas[[i]]), 1), moderators=primaryStudies$moderators[[i]],
                             maxModerators=length(primaryStudies$moderators[[i]]), startValues=primaryStudies$startValues[[i]],
                             rawData=primaryStudies$rawData[[i]], pairwiseN=primaryStudies$pairwiseNs[[i]],
                             source=paste(primaryStudies$source[[i]], collapse=", "))
      if (length(primaryStudies$moderators[[i]]) > maxLengthModeratorVector) maxLengthModeratorVector <- length(primaryStudies$moderators[[i]])
      # check matrix symmetry if matrix is provided
      if (!(primaryStudies$studyNumbers[i] %in% loadRawDataStudyNumbers)) {
        if (isSymmetric(primaryStudies$empcovs[[i]]) == FALSE) {
          cat(red$bold("The correlation matrix of study no.", i, "is not symmetric. Check and re-start!", "\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    # determine number of studies (if list is created by ctmaInit.R it is 1 element too long)
    tmp <- length(unlist(primaryStudies$studyNumbers)); tmp
    if ( is.na(primaryStudies$deltas[tmp]) &
         is.na(primaryStudies$sampleSizes[tmp]) &
         is.na(primaryStudies$deltas[tmp]) &
         is.null(dim(primaryStudies$pairwiseNs[tmp])) &
         is.null(dim(primaryStudies$emcovs[tmp])) &
         all(is.na(unlist(primaryStudies$moderators[tmp]))) &
         is.null(primaryStudies$rawData$fileName[tmp]) ) {
      n.studies <- tmp-1
      primaryStudies$studyNumbers[[tmp]] <- NULL
    } else {
      n.studies <- tmp
    }
    n.studies

    ### create pseudo raw data for all studies or load raw data if available & specified
    empraw <- lags <- moderators <- emprawMod <- allSampleSizes <- lostN <- overallNDiff <- relativeNDiff <- list()
    emprawLong <- list()
    if (is.na(primaryStudies$deltas[length(primaryStudies$deltas)])) primaryStudies$deltas <- primaryStudies$deltas[-length(primaryStudies$deltas)]
    allTpoints <-unlist(lapply(primaryStudies$deltas, function(extract) length(extract)+1)); allTpoints # more complicated thasn expected
    manifestNames <- paste0("V", 1:n.latent); manifestNames
    for (i in 1:n.studies) {
      if (!(i %in% loadRawDataStudyNumbers)) {
        currentVarnames <- c()
        currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        currentEmpcov <- (lapply(studyList, function(extract) extract$empcov))[[i]]; currentEmpcov
        currentLags <- (lapply(studyList, function(extract) extract$delta_t))[[i]]; currentLags
        currentPairwiseN <- (lapply(studyList, function(extract) extract$pairwiseN))[[i]]; currentPairwiseN
        currentModerators <- (lapply(studyList, function(extract) extract$moderators))[[i]]; currentModerators
        currentVarnames <- c()
        for (j in 1:(currentTpoints)) {
          for (h in 1:n.latent) {
            currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
          }
        }
        tmp <- suppressWarnings(pRaw(empCovMat=currentEmpcov, empN=currentSampleSize, empNMat=currentPairwiseN))
        empraw[[i]] <- tmp$data
        lostN[[i]] <- tmp$lostN
        overallNDiff[[i]] <- tmp$overallLostN
        relativeNDiff[[i]] <- tmp$relativeLostN
        lags[[i]] <- matrix(currentLags, nrow=dim(empraw[[i]])[1], ncol=currentTpoints-1, byrow=TRUE)
        empraw[[i]] <- cbind(empraw[[i]], lags[[i]]);
        colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))
        empraw[[i]] <- as.data.frame(empraw[[i]])
        ## START correction of current lags if entire time point is missing for a case
        # wide to long
        emprawLongTmp <- ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=manifestNames)
        emprawLongTmp <- suppressMessages(ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        emprawLongTmp <- emprawLongTmp[, ][ apply(emprawLongTmp[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        # eliminate rows where time is NA
        emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
        # inrervalise
        emprawWide <- suppressMessages(ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=manifestNames))
        # restore
        empraw[[i]] <- as.data.frame(emprawWide)
        # END correction
        # add potential moderators
      }

      # load raw data on request
      if ( i %in% loadRawDataStudyNumbers ) {

        tmp1 <- studyList[[i]]$rawData$fileName; tmp1
        tmp2 <- studyList[[i]]$rawData$header; tmp2
        tmp3 <- studyList[[i]]$rawData$dec; tmp3
        tmp4 <- studyList[[i]]$rawData$sep; tmp4
        tmpData <- read.table(file=tmp1,
                              header=tmp2,
                              dec=tmp3,
                              sep=tmp4)
        # replace missing values
        tmpData <- as.matrix(tmpData) # important: line below will not work without having data as a matrix
        tmpData[tmpData %in% studyList[[i]]$rawData$missingValues] <- NA
        empraw[[i]] <- as.data.frame(tmpData)
        ## START correction of current lags if entire time point is missing for a case
        # change variable names
        tmp1 <- dim(empraw[[i]])[2]; tmp1
        currentTpoints <- (tmp1 + 1)/(n.latent+1); currentTpoints
        colnames(empraw[[i]])[1:(currentTpoints * n.latent)] <- paste0(paste0("V", 1:n.latent), "_T", rep(0:(currentTpoints-1), each=n.latent))
        # wide to long
        emprawLongTmp <- ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=manifestNames)
        emprawLongTmp <- suppressMessages(ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        emprawLongTmp <- emprawLongTmp[apply(emprawLongTmp[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        # eliminate rows where time is NA
        emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
        # intervalise
        emprawWide <- suppressMessages(ctIntervalise(emprawWide,
                                                     Tpoints=currentTpoints,
                                                     n.manifest=n.latent,
                                                     manifestNames=manifestNames,
                                                     mininterval=minInterval))
        # restore
        empraw[[i]] <- as.data.frame(emprawWide)
        # END correction

        # Change the NAs provided for deltas if raw data are loaded
        for (h in 1:(currentTpoints-1)) {
          # temporarily replace dT = .001 by NA
          tmp1 <- grep("dT", colnames(empraw[[i]])); tmp1
          colnamesTmp <- colnames(empraw[[i]])[tmp1]; colnamesTmp
          emprawTmp <- as.matrix(empraw[[i]][, tmp1])
          colnames(emprawTmp) <- colnamesTmp
          emprawTmp[emprawTmp == minInterval] <- NA
          studyList[[i]]$delta_t[h] <- mean(emprawTmp[, paste0("dT", h)], na.rm=TRUE)
        }
      }

      # change sample size if entire cases were deleted
      studyList[[i]]$sampleSize <- (dim(empraw[[i]]))[1]
      allSampleSizes[[i]] <- dim(empraw[[i]])[1]
      currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
      currentTpoints <- allTpoints[[i]]; currentTpoints
      currentVarnames <- c()
      for (j in 1:(currentTpoints)) {
        for (h in 1:n.latent) {
          currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
        }
      }
      colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))

      # standardize (variables - not time lags) if option is chosen
      if (studyList[[i]]$rawData$standardize == TRUE) empraw[[i]][, currentVarnames] <- scale(empraw[[i]][, currentVarnames])

      # replace missing values for time lags dTi by minInterval (has to be so because dTi are definition variables)
      tmpData <- empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))]
      tmpData[is.na(tmpData)] <- minInterval
      empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))] <- tmpData

      # add moderators to loaded raw data
      # Save raw data  on request
      if ( i %in% saveRawData$studyNumbers ) {
        x1 <- paste0(saveRawData$fileName, i, ".dat"); x1
        write.table(empraw[[i]], file=x1, row.names=saveRawData$row.names, col.names=saveRawData$col.names,
                    sep=saveRawData$sep, dec=saveRawData$dec)
      }

      # augment pseudo raw data for stanct model
      {
        dataTmp <- empraw[[i]]
        dataTmp2 <- ctWideToLong(dataTmp, Tpoints=currentTpoints, n.manifest=n.latent, #n.TIpred = (n.studies-1),
                                 manifestNames=manifestNames)
        dataTmp3 <- ctDeintervalise(dataTmp2)
        dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
        # eliminate rows where ALL latents are NA
        dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        emprawLong[[i]] <- dataTmp3
      }

    } ### END for i ...
  } ### END Read user provided data and create list with all study information ###

  #######################################################################################################################
  ################################################### Some Statistics ###################################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("################# Compute Summary Statistics of Primary Studies #################"))
    print(paste0("#################################################################################"))

    ### some stats
    # Sample size
    allSampleSizes <-unlist(lapply(studyList, function(extract) extract$sampleSize)); allSampleSizes
    if (is.na(allSampleSizes[length(allSampleSizes)])) allSampleSizes <- allSampleSizes[-length(allSampleSizes)]; allSampleSizes
    overallSampleSize <- sum(allSampleSizes, na.rm=TRUE); overallSampleSize
    meanSampleSize <- mean(allSampleSizes, na.rm=TRUE); meanSampleSize
    maxSampleSize <- max(allSampleSizes, na.rm=TRUE); maxSampleSize
    minSampleSize <- min(allSampleSizes, na.rm=TRUE); minSampleSize
    # lags
    # Will be computed again below if raw data are loaded #
    allDeltas <-unlist(lapply(studyList, function(extract) extract$delta_t)); allDeltas
    #if (length(loadRawDataStudyNumbers) < 1) allDeltas <- allDeltas[-length(allDeltas)]; allDeltas # already corrected if any raw data are used
    if (is.na(allDeltas[length(allDeltas)])) allDeltas <- allDeltas[-length(allDeltas)]; allDeltas
    meanDelta <- mean(allDeltas, na.rm=TRUE); meanDelta
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    minDelta <- min(allDeltas, na.rm=TRUE); minDelta
    # Time points
    # do not overwrite all Tpoints - has been solved before
    # allTpoints <- unlist(lapply(studyList, function(extract) extract$timePoints)); allTpoints; length(allTpoints)
    # if (is.na(allTpoints[length(allTpoints)])) allTpoints <- allTpoints[-length(allTpoints)]; allTpoints
    allTpoints; length(allTpoints)
    overallTpoints <- sum(allTpoints, na.rm=TRUE); overallTpoints
    meanTpoints  <- mean(allTpoints, na.rm=TRUE); meanTpoints
    maxTpoints  <- max(allTpoints, na.rm=TRUE); maxTpoints
    minTpoints  <- min(allTpoints, na.rm=TRUE); minTpoints
    # Moderators Information
    statisticsList <- list(originalStudyNumbers=unlist(primaryStudies$studyNumber),
                           allSampleSizes=allSampleSizes, overallSampleSize=overallSampleSize, meanSampleSize=meanSampleSize,
                           maxSampleSize=maxSampleSize, minSampleSize=minSampleSize,
                           allDeltas=allDeltas, meanDelta=meanDelta, maxDelta=maxDelta, minDelta=minDelta,
                           allTpoints=allTpoints, overallTpoints=overallTpoints, meanTpoints=meanTpoints, maxTpoints=maxTpoints, minTpoints=minTpoints)

  } ### END Some Statistics ###

  #######################################################################################################################
  ############################# Create ctsem model template to fit all primary studies ##################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("############ Create ctsem Model Template to fit all Primary Studies #############"))
    print(paste0("#################################################################################"))

    manifestNames <- c()
    for (i in 1:n.latent) manifestNames <- c(manifestNames, paste0("V",i))
    driftNames <- c()
    for (i in 1:(n.latent)) {
      for (j in 1:(n.latent)) {
        driftNames <- c(driftNames, paste0("V",i,"toV", j))
      }
    }
    cintNames <- c()
    for (i in 1:(n.latent)) {
      cintNames <- c(cintNames, paste0("cint",i))
    }
    # general ctsem model template
    ctsemModelTemplate <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=2, manifestNames=manifestNames,    # 2 waves in the template only
                                  DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent),
                                  LAMBDA=diag(n.latent),
                                  type='stanct',
                                  #CINT=matrix(cintNames, nrow=n.latent, ncol=1),
                                  CINT=matrix(0, nrow=n.latent, ncol=1),
                                  T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                  MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                  MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent))


    # ctsem models for each primary study with the correct number of time points
    ctsemModel <- list()
    counter <- 1
    for (i in unique(unlist(allTpoints))) {
      helper <- ctsemModelTemplate
      helper$Tpoints <- i
      ctsemModel[[counter]] <- helper
      counter <- counter +1
    }

  } ### END Create ctsem model template to fit all primary studies ###

  #######################################################################################################################
  ##################################### Check Specification of Primary Studies ##########################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("#################### Check Specification of Primary Studies #####################"))
    print(paste0("#################################################################################"))

    if (length(saveSingleStudyModelFit) == 1){
      if ((activateRPB==TRUE) &  (silentOverwrite==FALSE)) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      if (silentOverwrite==FALSE) {
        cat(red$bold("You have indicated that you want to save SingleStudyModelFits, but have not selected any study/studies to save the fit for!","\n"))
        cat(red("Would you like to save the SingleStudyModelFits for ALL studies??","\n"))
        cat(blue("Press 'y' to save ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(blue("select the study/studies you whish to save the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((blue("Please press 'y' to save ALL, 's' to specify the study/studies to save, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in unlist(primaryStudies$studyNumber)) saveSingleStudyModelFit <- c(saveSingleStudyModelFit, i)
        } else if (char == 's' | char == 'S') {
          cat(blue("Which SingleStudyModelFits would you like to save?", "\n"))
          cat(blue("Please enter the no. of study/studies separated by commas!", "\n"))
          chars <- readline(" ")
          chars <- gsub(" ", "", chars, fixed = TRUE)
          chars <- unlist(strsplit(chars, ","))
          chars
          saveSingleStudyModelFit <- c(saveSingleStudyModelFit, chars)
        } else if (char == 'q' | char == 'Q'){
          stop("Good luck for the next try!")
        }
      } else {
        saveSingleStudyModelFit <- c(saveSingleStudyModelFit, unlist(primaryStudies$studyNumber))
      }
    }

    if (length(loadSingleStudyModelFit) == 1){
      if ((activateRPB==TRUE) & (silentOverwrite==FALSE)) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      if (silentOverwrite==FALSE) {
        cat(red$bold("You have indicated that you want to load SingleStudyModelFits, but have not selected any study/studies to load the fit for!","\n"))
        cat(red("Would you like to load the SingleStudyModelFits for ALL studies??","\n"))
        cat(blue("Press 'y' to load ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(blue("select the study/studies you whish to load the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((blue("Please press 'y' to load ALL, 's' to specify the study/studies to load, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in unlist(primaryStudies$studyNumber)) loadSingleStudyModelFit <- c(loadSingleStudyModelFit, i)
        } else if (char == 's' | char == 'S') {
          cat(blue("Which SingleStudyModelFits would you like to load?", "\n"))
          cat(blue("Please enter the no. of study/studies separated by commas!", "\n"))
          chars <- readline(" ")
          chars <- gsub(" ", "", chars, fixed = TRUE)
          chars <- unlist(strsplit(chars, ","))
          loadSingleStudyModelFit <- c(loadSingleStudyModelFit, chars)
        } else if (char == 'q' | char == 'Q') {
          stop("Good luck for the next try!")
        }
      } else {
        loadSingleStudyModelFit <- c(loadSingleStudyModelFit, unlist(primaryStudies$studyNumber))
      }
    }
  } ### END check specification of primary studies ###

  #######################################################################################################################
  ##################################### Fit ctsem Model to each Primary Study ###########################################
  #######################################################################################################################
  {
    # loop through all primary studies
    studyFit <- studyFitCI <- studyFit_Minus2LogLikelihood <- studyFit_estimatedParameters <- list()
    model_Drift_Coef <- model_Drift_SE <- model_Drift_CI <- list()
    model_Diffusion_Coef <- model_Diffusion_SE <- model_Diffusion_CI <- list()
    model_T0var_Coef <- model_T0var_SE <- model_T0var_CI <- list()
    model_Cint_Coef <- model_Cint_SE <- model_Cint_CI <- list()
    resultsSummary <- list()
    for (i in 1:n.studies) {
      notLoadable <- TRUE
      if ( (length(loadSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) {
        tmp1 <- paste0(" LOADING SingleStudyFit ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ")
        tmp2 <- nchar(tmp1); tmp2
        tmp3 <- (81 - tmp2)/2; tmp3
        tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
        tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
        tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
        print(paste0("#################################################################################"))
        print(tmp6)
        print(paste0("#################################################################################"))
        x1 <- paste0(activeDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        file.exists(x1)
        if (file.exists(x1)) {
          notLoadable <- FALSE
          studyFit[[i]] <- readRDS(file=x1)
        } else {
          notLoadable <- TRUE
        }

      } # END if ( (length(loadSingleStudyModelFit) > 1) &  ...

      tmpLogic <- 0
      if ((studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) & (notLoadable == TRUE) ) tmpLogic <- 1
      if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) tmpLogic <- 1
      if (tmpLogic == 1) {

        print(paste0("#################################################################################"))
        print(paste0("################### Fitting SingleStudyModel ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ######################"))
        print(paste0("#################################################################################"))
        # select correct template
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        modelToSelect <- which(unique(allTpoints) == currentTpoints); modelToSelect
        currentModel <- ctsemModel[[modelToSelect]]; currentModel

        # FIT STANCT MODEL
        results <- ctStanFit(
          datalong = emprawLong[[i]],
          ctstanmodel = currentModel,
          savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
          stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
          iter=CoTiMAStanctArgs$iter,
          intoverstates=CoTiMAStanctArgs$intoverstates,
          binomial=CoTiMAStanctArgs$binomial,
          fit=CoTiMAStanctArgs$fit,
          intoverpop=CoTiMAStanctArgs$intoverpop,
          stationary=CoTiMAStanctArgs$stationary,
          plot=CoTiMAStanctArgs$plot,
          derrind=CoTiMAStanctArgs$derrind,
          optimize=CoTiMAStanctArgs$optimize,
          optimcontrol=CoTiMAStanctArgs$optimcontrol,
          nlcontrol=CoTiMAStanctArgs$nlcontrol,
          nopriors=CoTiMAStanctArgs$nopriors,
          chains=CoTiMAStanctArgs$chains,
          forcerecompile=CoTiMAStanctArgs$forcerecompile,
          savescores=CoTiMAStanctArgs$savescores,
          gendata=CoTiMAStanctArgs$gendata,
          control=CoTiMAStanctArgs$control,
          verbose=CoTiMAStanctArgs$verbose,
          warmup=CoTiMAStanctArgs$warmup,
          cores=coresToUse)

        studyFit[[i]] <- results
        studyFit[[i]]$resultsSummary <- summary(studyFit[[i]])

      } # END if (!(studyList[[i]]$originalStudyNo %in% ...

      # SAVE
      if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
        x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/"); x2
        CoTiMASaveFile(activateRPB, activeDirectory, studyFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
      } else {
        # SAVE 2 (previously unfitted models - not yet on request but by default)
        if ( (notLoadable == TRUE) ) {
          if (!(exists("saveSingleStudyModelFit[1]"))) saveSingleStudyModelFit[1] <- "safeSave"
          x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
          x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/")
          CoTiMASaveFile(activateRPB, activeDirectory, studyFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
        }
      }

      studyFit_Minus2LogLikelihood[[i]] <- studyFit[[i]]$stanfit$optimfit$f
      studyFit_estimatedParameters[[i]] <- length(studyFit[[i]]$stanfit$optimfit$par)
      resultsSummary <- summary(studyFit[[i]]); resultsSummary

      tmp <- grep("toV", rownames(resultsSummary$popmeans)); tmp
      model_Drift_Coef[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"], n.latent, byrow=TRUE)); model_Drift_Coef[[i]]
      names(model_Drift_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_Drift_Coef[[i]]

      model_Drift_SE[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"], n.latent, byrow=TRUE)); model_Drift_SE[[i]]
      names(model_Drift_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_Drift_SE[[i]]

      tmp1 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=TRUE)); tmp2
      model_Drift_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Drift_CI[[i]]
      tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                      paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
      names(model_Drift_CI[[i]]) <- tmp3; model_Drift_CI[[i]]

      tmp <- grep("diff", rownames(resultsSummary$popmeans)); tmp
      model_Diffusion_Coef[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"]); model_Diffusion_Coef[[i]]
      names(model_Diffusion_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_Diffusion_Coef[[i]]

      model_Diffusion_SE[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Sd"]); model_Diffusion_SE[[i]]
      names(model_Diffusion_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_Diffusion_SE[[i]]

      tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "2.5%"]; tmp1
      tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "97.5%"]; tmp2
      model_Diffusion_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Diffusion_CI[[i]]

      tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                      paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
      names(model_Diffusion_CI[[i]]) <- tmp3; model_Diffusion_CI[[i]]

      tmp <- grep("0var", rownames(resultsSummary$popmeans)); tmp
      model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Mean"]); model_T0var_Coef[[i]]
      names(model_T0var_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_T0var_Coef[[i]]

      model_T0var_SE[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Sd"]); model_T0var_SE[[i]]
      names(model_T0var_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_T0var_SE[[i]]

      tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]; tmp1
      tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "97.5%"]; tmp2
      model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]

      tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                      paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
      names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]

    } # END     for (i in 1:n.studies)

    # Combine summary information and fit statistics
    allStudies_Minus2LogLikelihood <- sum(unlist(studyFit_Minus2LogLikelihood)); allStudies_Minus2LogLikelihood
    allStudies_estimatedParameters <- sum(unlist(studyFit_estimatedParameters)); allStudies_estimatedParameters
    #allStudies_df <- ((n.latent * allTpoints) %*% ((n.latent * allTpoints) +1 )) / 2 -
    #  allStudies_estimatedParameters; allStudies_df
    n.par.first.lag <- ((2 * n.latent) * (2 * n.latent + 1)) / 2; n.par.first.lag
    n.par.later.lag <- ((2 * n.latent) * (2 * n.latent - 1)) / 2; n.par.later.lag
    n.later.lags <- allTpoints - n.latent; n.later.lags
    allStudies_df <- sum(n.later.lags * n.par.later.lag); allStudies_df

    allStudiesDRIFT_effects <- matrix(t(cbind(unlist(model_Drift_Coef), unlist(model_Drift_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    tmp1 <- names(model_Drift_Coef[[1]]); tmp1
    tmp2 <- rep("SE", length(tmp1)); tmp2
    colnames(allStudiesDRIFT_effects) <- c(rbind(tmp1, tmp2))

    source <- lapply(primaryStudies$source, function(extract) paste(extract, collapse=", ")); source
    source <- source[-length(source)]
    allStudiesDRIFT_effects_ext <- cbind(unlist(source), allStudiesDRIFT_effects)
    tmp <- allStudiesDRIFT_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesDRIFT_effects_ext <- tmp
    allStudiesDRIFT_effects_ext

    allStudiesDriftCI <- matrix(unlist(model_Drift_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesDriftCI) <- names(model_Drift_CI[[1]])
    allStudiesDiffusionCI <- matrix(unlist(model_Diffusion_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesDiffusionCI) <- names(model_Diffusion_CI[[1]])
    allStudiesT0varCI <- matrix(unlist(model_T0var_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesT0varCI) <- names(model_T0var_CI[[1]])
    allStudiesCI <- t(rbind(t(allStudiesDriftCI), t(allStudiesDiffusionCI), t(allStudiesT0varCI)))
    allStudiesCI <- cbind(allStudiesDRIFT_effects_ext, allStudiesCI)
    rownames(allStudiesCI) <- rownames(allStudiesDRIFT_effects_ext)

    # Label summary table
    rownames(allStudiesDRIFT_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDRIFT_effects_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
    newColNames <- c()
    for (j in 1:n.latent) {
      for (h in 1:n.latent) {
        newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)")
      }
    }
    colnames(allStudiesDRIFT_effects) <- newColNames;
    newColNames <- c("Source", newColNames)
    colnames(allStudiesDRIFT_effects_ext) <- newColNames;

    # check single study results
    if (checkSingleStudyResults == TRUE) {
      print(allStudiesDRIFT_effects_ext)
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(blue(" Press 'q' to quit or any other key to continue. Press ENTER afterwards."))
      char <- readline(" ")
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    DRIFTCoeff <- matrix(unlist(model_Drift_Coef), n.studies, n.latent^2, byrow=TRUE); DRIFTCoeff
    DRIFTSE <- matrix(unlist(model_Drift_SE), n.studies, n.latent^2, byrow=TRUE); DRIFTSE

    if (n.studies < 2) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(blue("Only a single primary study was handed over to CoTiMA. No further (meta-) analyses can be conducted."))
      cat(" ", sep="\n")
      cat(blue("I guess this stop is intended! You could ignore further warning messages such as \"sqrt(n-2): NaNs produced\""))
      cat(" ", sep="\n")
      #cat(blue("Hopefully the solution is reasonable. Set refits to a higher value if this is not the case."))
      cat(" ", sep="\n")
      cat(" ", sep="\n")
      stop("No further errors!")
      cat(" ", sep="\n")

    }

  } ### END fit ctsem model to each primary study

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  st <- paste0("Computation started at: ", start.time); st
  et <- paste0("Computation ended at: ", end.time); et
  tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt



  ###########################################################################################################
  ###################################### SAVE RESULTS AND FIGURES ###########################################
  ###########################################################################################################

  print(paste0("#################################################################################"))
  print(paste0("########################### Save Results and Figures ############################"))
  print(paste0("#################################################################################"))

  resultsFileName <- paste0(resultsFilePrefix, ".txt"); resultsFileName
  sink(file = resultsFileName, append = TRUE, type = c("output"), split = TRUE)

  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("##################### Continuous Time Meta Analysis (CoTiMA)  #########################", "", sep="\n"))
  cat(paste("####################  (cf. Dormann, Guthier, & Cortina, 2019) ##########################", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat(st, sep="\n")
  cat(et, sep="\n")
  cat(tt, sep="\n")
  cat(" ", "", sep="\n")
  cat(paste("########################### contact: cdormann@uni-mainz.de #############################", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("--------------------------------- CoTiMA Parameters ------------------------------------", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat("Number of latent variables (n.latent): ", n.latent,  sep="")
  cat(" ", "", sep="\n")
  cat("Allows parallel processing on unix-like computers (e.g., Mac;) using # cores: ", coresToUse,  sep="")
  cat(" ", "", sep="\n")
  cat("Selects from the vector of moderator values (the 1st is used if not otherwise specificed) (moderatorNumber): ", moderatorNumber,  sep="")
  cat(" ", "", sep="\n")
  cat("Rounding used in output (digits): ", digits,  sep="")
  cat(" ", "", sep="\n")
  cat("Displays estimates from single study ctsem models and waits for user inoput to continue (checkSingleStudyResults): ", checkSingleStudyResults,  sep="")
  cat(" ", "", sep="\n")
  if (is.null(saveSingleStudyModelFit)) x <- "no" else {
    x <- paste(rep(saveSingleStudyModelFit[1], length(saveSingleStudyModelFit)-1),
               " studyFit", saveSingleStudyModelFit[-1], ".rds", sep="") }
  cat("Save the fit of single study ctsem models (saveSingleStudyModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadSingleStudyModelFit)) x <- "no" else {
    x <- paste(rep(loadSingleStudyModelFit[1], length(loadSingleStudyModelFit)-1),
               " studyFit", loadSingleStudyModelFit[-1], ".rds", sep="") }
  cat("load the fit of single study ctsem models (loadSingleStudyModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  cat("The active directory (activeDirectory) is: ", activeDirectory,  sep="\n")
  cat(" ", "", sep="\n")
  cat("The prefix for the result file that is created (resultsFileName): ", resultsFileName,  sep="")
  cat(" ", "", sep="\n")
  cat("The prefix for the figure files and the model fit files that are created (saveFilePrefix): ", saveFilePrefix,  sep="")
  cat(" ", "", sep="\n")
  cat(" ", "", sep="\n")
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("------------------------------ Primary Study Statistics---------------------------------", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat(paste("-------------------------------- Sample Sizes ------------------------------------------", "", sep="\n"))
  cat("Sample Sizes (k):       ", unlist(allSampleSizes), "", sep="\n")
  cat(paste("Reported sample sizes should be treated with some caution. Primary studies with         ", "", sep="\n"))
  cat(paste("correlations matrices with pairwise deleted missing values or raw data with missings    ", "", sep="\n"))
  cat(paste("do not have a single number representing the sample size. However, the internal         ", "", sep="\n"))
  cat(paste("algorithm used to transform the correlation matrix with pairwise deleted missing values ", "", sep="\n"))
  cat(paste("into pseudo raw data estimates how large the sample size had to be at least in order to ", "", sep="\n"))
  cat(paste("reproduce the pattern of pairwise N provided for the respective primary studies.", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat("Overall Sample Size (N):", overallSampleSize)
  cat(" ", "", sep="\n")
  for (i in 1:n.studies) {
    if (dim(primaryStudies$pairwiseNs[[i]])[1] > 0) {
      cat(paste("Primary study No.", i, "had a correlation matrix with pairwise deleted cases. N was:    ", "", sep=" "))
      print(round(primaryStudies$pairwiseNs[[i]], digits))
      cat(" ", "", sep="\n")
      cat(paste("Some cases were lost during pseudo raw data generation. Lost N was                      ", "", sep="\n"))
      print(round(lostN[[i]], digits))
      cat(paste("Overall lost N was: ", overallNDiff[[i]], "", sep="\n"))
      cat(paste("Relative lost N was: ", round(relativeNDiff[[i]], digits), "", sep="\n"))
      cat(" ", "", sep="\n")
    }
  }
  cat(paste("---------------------------------- Time Lags -------------------------------------------", "", sep="\n"))
  cat("All Time Lags (deltas): ", sort(unlist(allDeltas)))
  cat(" ", "", sep="\n")
  cat("Mean Time Lag (delta):  ", meanDelta)
  cat(" ", "", sep="\n")
  cat(paste("--------------------------------- Time Points ------------------------------------------", "", sep="\n"))
  cat("All Time Points (Waves):", sort(unlist(allTpoints)))
  cat(" ", "", sep="\n")
  cat("Sum of all Time Points: ", overallTpoints)
  cat(" ", "", sep="\n")
  cat("Mean No. of Time Points:", round(meanTpoints, digits))
  cat(" ", "", sep="\n")
  cat(" ", "", sep="\n")
  cat(" ", "", sep="\n")

  if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n", st, "\n", et, "\nAnalysis successfully completed. \nThank you for using CoTiMA.\nHave a nice day!"))}

  sink()

  if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

  if (exists("allDriftCI")) tmp <- allDriftCI else tmp <- NULL

  results <- list(activeDirectory=activeDirectory, #sourceDirectory=NULL,
                  plot.type="drift", model.type="stanct",
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=studyList, studyFitList=studyFit,
                  emprawList=empraw, statisticsList=statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=list(DRIFT=names(model_Drift_Coef[[1]]), DIFFUSION=names(model_Diffusion_Coef[[1]]), T0VAR=names(model_T0var_Coef[[1]])),
                  summary=(list(model="all drift free (het. model)",
                                estimates=allStudiesDRIFT_effects_ext,
                                confidenceIntervals=allStudiesCI,
                                minus2ll= round(allStudies_Minus2LogLikelihood, digits),
                                n.parameters = round(allStudies_estimatedParameters, digits),
                                df= c(round(allStudies_df, digits)))))

  class(results) <- "CoTiMAFit"

  invisible(results)

}   # end function definition
