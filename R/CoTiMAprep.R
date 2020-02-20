### Version 1.0.0.3



#debug <- 0
#if (debug == 1) {
#  ##### ENTER DEBUG INFORMATION BELOW THE FOLLOWING #######
#  activeDirectory=NULL
#  sourceDirectory= NULL
#  resultsFilePrefix="CoTiMAprep"           # the prefix for the result file that is created
#  saveFilePrefix="CoTiMAprep"
#  loadFilePrefix="CoTiMAprep"
#  primaryStudies=NULL                    #NEW: list of lists: list(deltas, sampleSizes, empcovs, moderators, startValues, studyNumbers)
#  activateRPB=FALSE                      #set to TRUE to receive push messages with CoTiMA notifications on your phone
#  digits=4
#  checkSingleStudyResults=TRUE          # displays estimates from single study ctsem models and waits for user inoput to continue
#  n.latent=NULL
#  loadRawData=list()
#  saveRawData=list()
#  retryattempts=30                       # number of iterations to be used by ctsem
#  refits=5                               # number of re-fits used by this CoTiMA Function
#  extraRefits=c()                        # (vector of) study number(s) that should be re-fitted even more (difficult models)
#  factorExtraRefits=1                    # factor by which "refits" are multiplied for difficult models
#  NPSOL=FALSE                            # request NPSOL optimizer (sometimes better estimates, e.g., for confidence intervals)
#  coresToUse=c(1)                        # could be set to -n (usually -1) which is subtracted from the overall number available (on non Windows machines)
#  tryHard=TRUE
#  confidenceIntervals=FALSE              # estimate confidence intervals (for all requested models)
#  silentOverwrite=FALSE
#  saveSingleStudyModelFit=c()            # save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#  loadSingleStudyModelFit=c()            # load the fit of single study ctsem models

#  ##### ENTER DEBUG INFO HERE #######
#  activeDirectory="/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/ARTIKEL/COTIMA BURNOUT/NEU 2019/PB SUBMISSION (BUL-2019-1709)/REVISION/_TRYOUT/"
#  sourceDirectory="/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/_TRYOUT/"
#  resultsFilePrefix="CoTiMAprep"           # the prefix for the result file that is created
#  saveFilePrefix="CoTiMAprep"
#  loadFilePrefix="CoTiMAprep"
#
#  # Primary Study Information
#  primaryStudies=compiledListOfPrimaryStudies                    #NEW: list of lists: list(deltas, sampleSizes, empcovs, moderators, startValues, studyNumbers)
#
#  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
#  #activateRPB=FALSE,                      #set to TRUE to receive push messages with CoTiMA notifications on your phone
#  #checkSingleStudyResults=TRUE,          # displays estimates from single study ctsem models and waits for user inoput to continue
#  digits=4
#
#  # General Model Setup
#  n.latent=2
#
#  # Load/save Raw Data
#  #loadRawData=list(),
#  #saveRawData=list(),
#
#  # Fitting Parameters
#  #retryattempts=30,                       # number of iterations to be used by ctsem
#  #refits=5,                               # number of re-fits used by this CoTiMA Function
#  #extraRefits=c(),                        # (vector of) study number(s) that should be re-fitted even more (difficult models)
#  #factorExtraRefits=1,                    # factor by which "refits" are multiplied for difficult models
#  #NPSOL=FALSE,                            # request NPSOL optimizer (sometimes better estimates, e.g., for confidence intervals)
# coresToUse=c(-1)                        # could be set to -n (usually -1) which is subtracted from the overall number available (on non Windows machines)
#  #tryHard=TRUE,
#  # Further Analyses
#  #confidenceIntervals=FALSE,              # estimate confidence intervals (for all requested models)
#
#  ## Save computed model fits or load previous fits (then provide (old) file prefix)
#  #silentOverwrite=FALSE,
#
#  # ctsem models for all primasry studies (always required for getting good starting values (and model fit for heterogeneity model without fitting it))
#  #saveSingleStudyModelFit=c("prep1")            # save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#  loadSingleStudyModelFit=c("prep1")
#}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' CoTiMAprep
#'
#' @param activeDirectory
#' @param sourceDirectory
#' @param resultsFilePrefix
#' @param saveFilePrefix
#' @param loadFilePrefix
#' @param primaryStudies
#' @param activateRPB
#' @param checkSingleStudyResults
#' @param digits
#' @param n.latent
#' @param loadRawData
#' @param saveRawData
#' @param retryattempts
#' @param refits
#' @param extraRefits
#' @param factorExtraRefits
#' @param NPSOL
#' @param coresToUse
#' @param tryHard
#' @param confidenceIntervals
#' @param silentOverwrite
#' @param saveSingleStudyModelFit
#' @param loadSingleStudyModelFit
#'
#' @return
#' @export
#'
#' @examples
#'
CoTiMAprep <- function(
  # Directory names and file names
  activeDirectory=NULL,
  sourceDirectory= NULL,
  resultsFilePrefix="CoTiMAprep",           # the prefix for the result file that is created
  saveFilePrefix="CoTiMAprep",
  loadFilePrefix="CoTiMAprep",

  # Primary Study Information
  primaryStudies=NULL,                    #NEW: list of lists: list(deltas, sampleSizes, empcovs, moderators, startValues, studyNumbers)

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
  refits=5,                               # number of re-fits used by this CoTiMA Function
  extraRefits=c(),                        # (vector of) study number(s) that should be re-fitted even more (difficult models)
  factorExtraRefits=1,                    # factor by which "refits" are multiplied for difficult models
  NPSOL=FALSE,                            # request NPSOL optimizer (sometimes better estimates, e.g., for confidence intervals)
  coresToUse=c(1),                        # could be set to -n (usually -1) which is subtracted from the overall number available (on non Windows machines)
  tryHard=TRUE,

  # Further Analyses
  confidenceIntervals=FALSE,              # estimate confidence intervals (for all requested models)

  ## Save computed model fits or load previous fits (then provide (old) file prefix)
  silentOverwrite=FALSE,

  # ctsem models for all primasry studies (always required for getting good starting values (and model fit for heterogeneity model without fitting it))
  saveSingleStudyModelFit=c(),            # save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
  loadSingleStudyModelFit=c()            # load the fit of single study ctsem models
)

{  # begin function definition (until end of file)

  # #######################################################################################################################
  # ############################################ load required subfunctions ###############################################
  # #######################################################################################################################
  #
  # {
  #   print(paste0("#################################################################################"))
  #   print(paste0("########################## DEFINE Required Functions ############################"))
  #   print(paste0("#################################################################################"))
  #
  #   source(paste0(sourceDirectory, "CoTiMASaveFile.R"))
  #   source(paste0(sourceDirectory, "CoTiMAfittingSSMF.R"))
  #   source(paste0(sourceDirectory, "CoTiMApseudoRawData.R"))
  #   source(paste0(sourceDirectory, "CoTiMAcombinePseudoRawData.R"))
  #
  # } ### END DEFINE required subfunctions ###



  # #######################################################################################################################
  # ############################################### attach further packages ###############################################
  # #######################################################################################################################
  # {
  #   print(paste0("#################################################################################"))
  #   print(paste0("############################ Attach Further Packages ############################"))
  #   print(paste0("#################################################################################"))
  #
  #
  #   # If OpenMx and ctsem are already attached, detaching them is required to enable OpenMx to run in parallel mode.
  #   if ("OpenMx" %in% (.packages())) suppressWarnings(detach("package:OpenMx", force=TRUE)) #, unload=TRUE)
  #   if ("ctsem" %in% (.packages())) suppressWarnings(detach("package:ctsem", force=TRUE)) #, unload=TRUE))
  #   Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
  #   library(ctsem)
  #   mxOption(key='Number of Threads', value=parallel::detectCores()) #now
  #
  #   # Attach all required packages that were not already attach with "library(ctsem)" before
  #   if (!("MASS" %in% (.packages()))) library(MASS)
  #   if (!("MBESS" %in% (.packages()))) library(MBESS)
  #   if (!("rootSolve" %in% (.packages()))) library(rootSolve)
  #   if (!("doParallel" %in% (.packages()))) library(doParallel)  # this is for parallel processing loops (not for internal parallel processing of OpenMx)
  #   if (!("crayon" %in% (.packages()))) library(crayon)
  #   if (!("psych" %in% (.packages()))) library(psych)
  #
  #   if (activateRPB==TRUE) {
  #     if("RPushbullet" %in% rownames(installed.packages()) == FALSE) {install.packages("RPushbullet")}
  #     if (!("RPushbullet" %in% (.packages()))) library(RPushbullet)
  #
  #   }
  # } ### END install and attach further packages ###

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
          primaryStudies <- compileListOfPrimaryStudies(selectedStudies=1:char)
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

    # check statistical power

    # check time unit

    if (resultsFilePrefix=="CoTiMAprep") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default results file prefix (CoTiMAprep) has been chosen.", "\n"))
      cat(blue("Press 'q' to quit and change or'c'to continue. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }
    rsultsFileName <- paste0(resultsFilePrefix, ".txt")

    if (saveFilePrefix=="CoTiMAprep") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default save file prefix (CoTiMAprep) has been chosen.", "\n"))
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

    # select moderatorNumber if no number is specified by user

    #if (length(saveSingleStudyFit) > 1 | length(loadSingleStudyFit) > 1 | saveHeterogeneityStudyFit == TRUE | loadHeterogeneityStudyFit == TRUE
    #    | saveHomogeneityStudyFit == TRUE | loadHomogeneityStudyFit == TRUE | !(is.null(loadFilePrefix))) {
    #  if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    #  cat(red$bold("DEPRECATED parameters detected!","\n"))
    #  cat(red("You're using parameters of an older Version of CoTiMA (prior v1.6)","\n"))
    #  stop(red("Please adjust your file to the newest version or use an older function"))
    #}

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

    numOfThreads <- round((coresToUse/refits - .4), 0); numOfThreads
    if(numOfThreads < 1) numOfThreads <- 1
    tmp <- paste0("No of Threads set to ", numOfThreads); tmp
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
    maxLengthModeratorVector <- 0
    loadRawDataStudyNumbers <- unlist(lapply(primaryStudies$rawData,
                                             function(extract) extract$studyNumbers)); loadRawDataStudyNumbers

    # resorting primary study information
    if (!(exists("moderatorNumber"))) moderatorNumber <- 1; moderatorNumber
    for (i in 1:length(unlist(primaryStudies$studyNumbers))) {
      studyList[[i]] <- list(studyNumber=i, empcov=primaryStudies$empcovs[[i]], delta_t=primaryStudies$deltas[[i]],
                             sampleSize=primaryStudies$sampleSizes[[i]], originalStudyNo=primaryStudies$studyNumber[[i]],
                             timePoints=(length(primaryStudies$deltas[[i]]) + 1), moderators=primaryStudies$moderators[[i]],
                             maxModerators=length(primaryStudies$moderators[[i]]), startValues=primaryStudies$startValues[[i]],
                             rawData=primaryStudies$rawData[[i]], pairwiseN=primaryStudies$pairwiseNs[[i]])
      if (length(primaryStudies$moderators[[i]]) > maxLengthModeratorVector) maxLengthModeratorVector <- length(primaryStudies$moderators[[i]])
      # check matrix symmetry if matrix is provided
      if (!(primaryStudies$studyNumbers[i] %in% loadRawDataStudyNumbers)) {
        if (isSymmetric(primaryStudies$empcovs[[i]]) == FALSE) {
          cat(red$bold("The correlation matrix of study no.", i, "is not symmetric. Check and re-start!", "\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    # determine number of studies (if list is created by CoTiMAcompileListOfPrimaryStudies.R it is 1 element too long)
    tmp <- length(unlist(primaryStudies$studyNumbers)); tmp
    if ( is.na(primaryStudies$deltas[tmp]) & is.na(primaryStudies$sampleSizes[tmp]) & is.na(primaryStudies$deltas[tmp]) &
         is.null(dim(primaryStudies$pairwiseNs[tmp])) & is.null(dim(primaryStudies$emcovs[tmp])) &
         is.na(primaryStudies$moderators[tmp]) & is.null(primaryStudies$rawData$fileName[tmp]) ) {
      n.studies <- tmp-1
      primaryStudies$studyNumbers[[tmp]] <- NULL
    } else {
      n.studies <- tmp
    }
    n.studies

    # correction if no moderatorNumber is specified

    ### create pseudo raw data for all studies or load raw data if available & specified
    empraw <- lags <- moderators <- emprawMod <- allSampleSizes <- lostN <- overallNDiff <- relativeNDiff <- list()
    allTpoints <-lapply(studyList, function(extract) extract$timePoints); allTpoints
    for (i in 1:n.studies) {
      if (!(i %in% loadRawDataStudyNumbers)) {
        currentVarnames <- c()
        currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        currentEmpcov <- (lapply(studyList, function(extract) extract$empcov))[[i]]; currentEmpcov
        currentLags <- (lapply(studyList, function(extract) extract$delta_t))[[i]]; currentLags
        currentPairwiseN <- (lapply(studyList, function(extract) extract$pairwiseN))[[i]]; currentPairwiseN
        #if (testModeratorModel == TRUE) {
        currentModerators <- (lapply(studyList, function(extract) extract$moderators))[[i]]; currentModerators
        #}
        currentVarnames <- c()
        for (j in 1:(currentTpoints)) {
          for (h in 1:n.latent) {
            currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
          }
        }
        tmp <- suppressWarnings(pseudoRawData(empCovMat=currentEmpcov, empN=currentSampleSize, empNMat=currentPairwiseN))
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
        emprawLong <- ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=paste0("V", 1:n.latent))
        emprawLong <- suppressMessages(ctDeintervalise(datalong = emprawLong, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        emprawLong <- emprawLong[, ][ apply(emprawLong[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        # eliminate rows where time is NA
        emprawLong <- emprawLong[which(!(is.na(emprawLong[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctLongToWide(emprawLong, id='id', time='time', manifestNames=paste0("V", 1:n.latent)))
        # inrervalise
        emprawWide <- suppressMessages(ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=paste0("V", 1:n.latent)))
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
        # replace mising values
        tmpData <- as.matrix(tmpData) # important: line below will not work without having data as a matrix
        tmpData[tmpData %in% studyList[[i]]$rawData$missingValues] <- NA
        empraw[[i]] <- as.data.frame(tmpData)
        ## START correction of current lags if entire time point is missing for a case
        # change variable names
        tmp1 <- dim(empraw[[i]])[2]; tmp1
        currentTpoints <- (tmp1 + 1)/(n.latent+1); currentTpoints
        colnames(empraw[[i]])[1:(currentTpoints * n.latent)] <- paste0(paste0("V", 1:n.latent), "_T", rep(0:(currentTpoints-1), each=n.latent))
        # wide to long
        emprawLong <- ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=paste0("V", 1:n.latent))
        emprawLong <- suppressMessages(ctDeintervalise(datalong = emprawLong, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        emprawLong <- emprawLong[, ][ apply(emprawLong[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        # eliminate rows where time is NA
        emprawLong <- emprawLong[which(!(is.na(emprawLong[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctLongToWide(emprawLong, id='id', time='time', manifestNames=paste0("V", 1:n.latent)))
        # inrervalise
        emprawWide <- suppressMessages(ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=n.latent, manifestNames=paste0("V", 1:n.latent)))
        # restore
        empraw[[i]] <- as.data.frame(emprawWide)
        # END correction

        # Change the NAs provided for deltas if raw data are loaded
        for (h in 1:(currentTpoints-1)) studyList[[i]]$delta_t[h] <- mean(empraw[[i]][, paste0("dT", 1)], na.rm=TRUE)

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

      # replace missing values for time lags dTi by .000001 (has to be so because dTi are definition variables)
      tmpData <- empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))]
      tmpData[is.na(tmpData)] <- .001
      empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))] <- tmpData
      # add moderators to loaded raw data
      # Save raw data  on request
      if ( i %in% saveRawData$studyNumbers ) {
        x1 <- paste0(saveRawData$fileName, i, ".dat"); x1
        write.table(empraw[[i]], file=x1, row.names=saveRawData$row.names, col.names=saveRawData$col.names,
                    sep=saveRawData$sep, dec=saveRawData$dec)
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
    allSampleSizes <- allSampleSizes[-length(allSampleSizes)]; allSampleSizes
    overallSampleSize <- sum(allSampleSizes, na.rm=TRUE); overallSampleSize
    meanSampleSize <- mean(allSampleSizes, na.rm=TRUE); meanSampleSize
    maxSampleSize <- max(allSampleSizes, na.rm=TRUE); maxSampleSize
    minSampleSize <- min(allSampleSizes, na.rm=TRUE); minSampleSize
    # lags
    # Will be computed again below if raw data are loaded #
    allDeltas <-unlist(lapply(studyList, function(extract) extract$delta_t)); allDeltas
    allDeltas <- allDeltas[-length(allDeltas)]; allDeltas
    meanDelta <- mean(allDeltas, na.rm=TRUE); meanDelta
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    minDelta <- min(allDeltas, na.rm=TRUE); minDelta
    # Time points
    allTpoints <-unlist(lapply(studyList, function(extract) extract$timePoints)); allTpoints
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
    statisticsList
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

    #originalStudyNo <- unlist(lapply(studyList, function(extract) extract$originalStudyNo)); originalStudyNo

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
    #FitList=list()

    studyFit <- studyFitCI <- studyFit_Minus2LogLikelihood <- studyFit_estimatedParameters <- list()
    study_Drift_Coef <- study_Drift_SE <- list()
    study_Diffusion_Coef <- study_Diffusion_SE <- list()
    study_T0var_Coef <- study_T0var_SE <- list()
    study_Cint_Coef <- study_Cint_SE <- list()
    origRefits <- refits; origRefits
    #################################################################################
    for (i in 1:n.studies) {
      #i <- 1
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
        #loadSingleStudyModelFit
        x1 <- paste0(activeDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        if (file.exists(x1)) {
          notLoadable <- FALSE
          studyFit[[i]] <- readRDS(file=x1)
        } else {
          notLoadable <- TRUE
        }

        ## STORE
        #if (exists("studyFit[[i]]")) {
        #  FitList$tmpName <-studyFit[[i]]
        #  tmp <- names(FitList); tmp
        #  tmp[length(tmp)] <- paste0("SingleStudyFit", i); tmp
        #  names(FitList) <- tmp; names(FitList)
        #}
      }

      if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) | (notLoadable == TRUE) ) {

        print(paste0("#################################################################################"))
        print(paste0("################### Fitting SingleStudyModel ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ######################"))
        print(paste0("#################################################################################"))
        if (studyList[[i]]$originalStudyNo %in% extraRefits) refits <- factorExtraRefits * refits
        if (!(studyList[[i]]$originalStudyNo %in% extraRefits)) refits <- origRefits
        # select correct template
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        modelToSelect <- which(unique(allTpoints) == currentTpoints); modelToSelect
        currentModel <- ctsemModel[[modelToSelect]]; currentModel

        # Set Start Values
        currentStartValues <- (lapply(studyList, function(extract) extract$startValues))[[i]]; currentStartValues
        if (is.na(currentStartValues[1])) currentStartValues <- NULL

        # FIT (NOT LOAD)
        mxOption(NULL, 'Number of Threads', numOfThreads)
        results <- CoTiMAfittingSSMF(coresToUse=coresToUse, empraw=empraw[[i]], currentModel=currentModel,
                                     refits=refits, retryattempts=retryattempts,
                                     singleModelStartValues=currentStartValues)
        mxOption(key='Number of Threads', value=parallel::detectCores())

        # Extract estimates & statistics
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        studyFit[[i]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]

        # STORE
        #FitList$tmpName <-studyFit[[i]]
        #tmp <- names(FitList); tmp
        #tmp[length(tmp)] <- paste0("SingleStudyFit", studyList[[i]]$originalStudyNo); tmp
        #names(FitList) <- tmp; names(FitList)

        # SAVE
        if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
          x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
          x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/"); x2
          CoTiMASaveFile(activateRPB, activeDirectory, studyFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
        } else {
          # SAVE 2 (previously unfitted models - not yet on request but by defauls)
          if ( (notLoadable == TRUE) ) {
            if (!(exists("saveSingleStudyModelFit[1]"))) saveSingleStudyModelFit[1] <- "safeSave"
            x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
            x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/")
            CoTiMASaveFile(activateRPB, activeDirectory, studyFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
          }
        }
      }
      studyFit_Minus2LogLikelihood[[i]] <- studyFit[[i]]$mxobj$output$Minus2LogLikelihood
      studyFit_estimatedParameters[[i]] <- length(studyFit[[i]]$mxobj$output$estimate)
      study_Drift_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[1:(n.latent^2)]; study_Drift_Coef[[i]]
      study_Drift_SE[[i]] <- studyFit[[i]]$mxobj$output$standardErrors[1:(n.latent^2)]; study_Drift_SE[[i]]
      study_Diffusion_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("diffusion", names(studyFit[[1]]$mxobj$output$estimate))]
      study_Diffusion_SE[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("diffusion", rownames(studyFit[[1]]$mxobj$output$standardErrors))]
      study_T0var_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("T0var", names(studyFit[[1]]$mxobj$output$estimate))]
      study_T0var_SE[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("T0var", rownames(studyFit[[1]]$mxobj$output$standardErrors))]
      study_Cint_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("cint", names(studyFit[[1]]$mxobj$output$estimate))]
      study_Cint_SE[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("cint", rownames(studyFit[[1]]$mxobj$output$standardErrors))]
    }

    # Compute confidence intervals
    if (confidenceIntervals == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("######################## Computing Confidence Intervals #########################"))
      print(paste0("#################################################################################"))

      # LOAD
      if ( (length(loadSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) {
        notLoadable <- list()
        for (i in 1:n.studies) {
          print(paste0("#################################################################################"))
          tmp1 <- paste0(" LOADING SingleStudyFit with confidence intervals ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ")
          tmp2 <- nchar(tmp1); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
          print(paste0("#################################################################################"))
          print(tmp6)
          print(paste0("#################################################################################"))
          x1 <- paste0(activeDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFitCI", studyList[[i]]$originalStudyNo, ".rds"); x1
          if (file.exists(x1)) {
            notLoadable[[i]] <- FALSE
            studyFitCI[[i]] <- readRDS(file=x1)
          } else {
            notLoadable[[i]] <- TRUE
          }

          # STORE
          #if (notLoadable[[i]] == FALSE) {
          #  FitList$tmpName <-studyFitCI[[i]]
          #  tmp <- names(FitList); tmp
          #  tmp[length(tmp)] <- paste0("SingleStudyFitCI", i); tmp
          #  names(FitList) <- tmp; names(FitList)
          #}
        }
      }

      # FIT (NOT LOAD)
      if ( (length(loadSingleStudyModelFit) < 1) ) {
        mxOption(NULL, 'Number of Threads', numOfThreads)
        results <- mclapply(seq(1, n.studies, by=1),
                            function(studyNo) ctCI(ctfitobj=studyFit[[studyNo]], confidenceintervals = driftNames),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())

        # STORE
        #for (k in 1:n.studies) {
        #  # fits
        #  FitList$tmpName <- results[[k]]
        #  tmp <- names(FitList); tmp
        #  tmp[length(tmp)] <- paste0("SingleStudyFitCI", k); tmp
        #  names(FitList) <- tmp; names(FitList)
        #}

        # SAVE
        studyFitCI <- results; studyFitCI
        for (h in 1:n.studies) {
          if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[h]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
            x1 <- paste0(saveSingleStudyModelFit[1], " studyFitCI", studyList[[h]]$originalStudyNo, ".rds"); x1
            x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/")
            CoTiMASaveFile(activateRPB, activeDirectory, results[[h]], x1, x2, silentOverwrite=silentOverwrite)
          }
        }
      }
    } # END Computing Confidence Intervals

    # Combine summary information and fit statistics
    allStudies_Minus2LogLikelihood <- sum(unlist(studyFit_Minus2LogLikelihood)); allStudies_Minus2LogLikelihood
    allStudies_estimatedParameters <- sum(unlist(studyFit_estimatedParameters)); allStudies_estimatedParameters
    allStudies_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
      allStudies_estimatedParameters; allStudies_df
    if (confidenceIntervals == TRUE) allDriftCI <-lapply(studyFitCI, function(extract) extract$mxobj$output$confidenceIntervals)
    allStudiesDRIFT_effects <- matrix(t(cbind(unlist(study_Drift_Coef), unlist(study_Drift_SE)) ), n.studies, 2*n.latent^2, byrow=T)

    # Label summary table
    rownames(allStudiesDRIFT_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    newColNames <- c()
    for (j in 1:n.latent) {
      for (h in 1:n.latent) {
        newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)")
      }
    }
    colnames(allStudiesDRIFT_effects) <- newColNames; round(allStudiesDRIFT_effects, digits)
    #FitList$allStudiesDRIFT_effects <- allStudiesDRIFT_effects

    # check single study results
    if (checkSingleStudyResults == TRUE) {
      print(round(allStudiesDRIFT_effects, digits))
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(blue(" Press 'q' to quit or any other key to continue. Press ENTER afterwards."))
      char <- readline(" ")
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    # Extract parameter estimates and paramter names to be used as start values later
    est <- unlist(lapply(studyFit, function(extract) extract$mxobj$output$estimate)); est
    startValues <- c(t(est)); startValues
    noOfParams <- length(studyFit[[1]]$mxobj$output$estimate); noOfParams
    headNames <- c(matrix(paste0(rep(c(paste0("Study_No_", primaryStudies$studyNumbers, "_")), noOfParams)),
                          nrow=noOfParams, ncol=n.studies, byrow=T)); headNames
    names(startValues) <- paste0(headNames, names(est)); startValues

    DRIFTCoeff <- matrix(unlist(study_Drift_Coef), n.studies, n.latent^2, byrow=TRUE); DRIFTCoeff
    DRIFTSE <- matrix(unlist(study_Drift_SE), n.studies, n.latent^2, byrow=TRUE); DRIFTSE

    if (n.studies < 2) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(blue("Only a single primary study was handed over to CoTiMA. No further (meta-) analyses can be conducted."))
      cat(" ", sep="\n")
      cat(blue("I guess this stop is intended! You could ignore further warning messages such as \"sqrt(n-2): NaNs produced\""))
      cat(" ", sep="\n")
      cat(blue("Hopefully the solution is reasonable. Set refits to a higher value if this is not the case."))
      cat(" ", sep="\n")
      #return(FitList)
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
  cat("Number of iterations to be used by ctsem (retryattempts): ", retryattempts,  sep="")
  cat(" ", "", sep="\n")
  cat("Number of re-fits used by this CoTiMA Function (refits): ", refits,  sep="")
  cat(" ", "", sep="\n")
  if(is.null(extraRefits)) {
    x <- "empty"
  } else {
    x <- extraRefits
  }
  cat("(Vector of) study number(s) that should be re-fitted even more (difficult models) (extraRefits): ", x,  sep="")
  cat(" ", "", sep="\n")
  cat("Factor by which refits are multiplied for difficult models (factorExtraRefits): ", factorExtraRefits,  sep="")
  cat(" ", "", sep="\n")
  cat("Request NPSOL optimizer (NPSOL): ", NPSOL,  sep="")
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

  allResults <- list(activeDirectory=activeDirectory, sourceDirectory=sourceDirectory,
                     coresToUse=coresToUse, n.studies=n.studies,
                     studyList=studyList, studyFitList=studyFit,
                     emprawList=empraw, statisticsList=statisticsList,
                     studyResults=list(DRIFT=study_Drift_Coef, DIFFUSION=study_Diffusion_Coef, T0VAR=study_T0var_Coef, CINT=study_Cint_Coef))
  saveRDS(allResults, paste0(activeDirectory, "CoTiMAprep.RDS"))

  invisible(list(activeDirectory=activeDirectory, sourceDirectory=sourceDirectory,
                 coresToUse=coresToUse, n.studies=n.studies,
                 studyList=studyList, studyFitList=studyFit,
                 emprawList=empraw, statisticsList=statisticsList,
                 studyResults=list(DRIFT=study_Drift_Coef, DIFFUSION=study_Diffusion_Coef, T0VAR=study_T0var_Coef, CINT=study_Cint_Coef)))
}   # end function definition
