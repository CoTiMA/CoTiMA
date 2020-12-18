#' ctmaInit
#'
#' @description Fits a ctsem model to a list of primary studies prepared by ctmaPrep.
#'
#' @param primaryStudies list of primary study information created with ctmaPrep
#' @param activeDirectory defines another active directory than the one used in ctmaPrep
#' @param activateRPB set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param checkSingleStudyResults Displays estimates from single study ctsem models and waits for user input to continue. Useful to check estimates before they are saved.
#' @param digits Number of digits used for rounding (in outputs)
#' @param n.latent Number of latent variables of the model (hast to be specified)!
#' @param n.manifest Number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param manifestVars Define the error variances of the manifests with a single time point using R-type matrix with nrow=n.manifest & ncol=n.manifest.
#' @param drift Labels for drift effects. Have to be either of the type V1toV2 or 0 for effects to be excluded, which is usually not recommended)
#' @param indVarying Control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means
#' @param saveRawData Save (created pseudo) raw date. List: saveRawData$studyNumbers, $fileName, $row.names, col.names, $sep, $dec
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#' @param silentOverwrite Overwrite old files without asking
#' @param saveSingleStudyModelFit save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#' @param loadSingleStudyModelFit load the fit of single study ctsem models
#' @param scaleTI scale TI predictors
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param nopriors if TRUE, any priors are disabled – sometimes desirable for optimization
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param iter number of interation (defaul = 1000). Sometimes larger values could be reqiured fom Baysian estimation
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#'
#' @importFrom  RPushbullet pbPost
#' @importFrom  crayon red blue
#' @importFrom  parallel detectCores
#' @importFrom  ctsem ctDeintervalise ctLongToWide ctIntervalise ctWideToLong ctModel ctStanFit
#' @importFrom  utils read.table write.table
#' @importFrom openxlsx addWorksheet writeData createWorkbook openXL saveWorkbook
#'
#' @export ctmaInit
#'
#' @examples
#' # Fit a ctsem model to all primary studies summarized in
#' # studyList_Ex1 and save fitted models
#' \dontrun{
#' CoTiMAInitFit_Ex1 <- ctmaInit(
#' primaryStudies=studyList_Ex1,
#' n.latent=2,
#' saveSingleStudyModelFit=c(2, 4, 17))
#'
#' summary(CoTiMAInitFit_Ex1)
#' }
#'
ctmaInit <- function(
  primaryStudies=NULL,
  activeDirectory=NULL,
  activateRPB=FALSE,
  checkSingleStudyResults=TRUE,
  digits=4,
  n.latent=NULL,
  n.manifest=0,
  lambda=NULL,
  manifestVars=NULL,
  drift=NULL,
  indVarying=FALSE,
  saveRawData=list(),
  coresToUse=c(1),
  silentOverwrite=FALSE,
  saveSingleStudyModelFit=c(),
  loadSingleStudyModelFit=c(),
  scaleTI=NULL,
  scaleTime=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  finishsamples=NULL,
  chains=NULL,
  iter=NULL,
  verbose=NULL
)

{  # begin function definition (until end of file)

  start.time <- Sys.time()
  options(scipen = 999) # turn scientific notation off.


  #######################################################################################################################
  ########################################### Check Model Specification #################################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("########################## Check Model Specification ############################"))
    print(paste0("#################################################################################"))

    if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMAStanctArgs$verbose}

    if (is.null(primaryStudies)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold(" List with lists of primary study information not specified!", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold("Should I try to get primary study information (e.g., empcov1, delta_t22) from the global environment?", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold("(This is the way how earlier version of CoTiMA were run but no longer recommended.)", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::blue("Press 'q' to quit and specify or 'c' to give it a try. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and specify primaryStudies or 'c' to try reading from global environment Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') {
        stop("Good luck for the next try!")
      } else {
        if (!(is.numeric(tryCatch(get(paste("sampleSize", 1, sep = "")), error = function(e) e)))) {
          cat(crayon::red$bold("Getting primary study information from the global environment failed", sep="\n"))
          cat(crayon::red$bold(" ", " ", sep="\n"))
          cat(crayon::blue("To test, I searched for sampleSize1 and could not find it!", "\n"))
          cat(crayon::red$bold(" ", " ", sep="\n"))
          stop("Good luck for the next try!")
        } else {
          cat(crayon::blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n"))
          char <- as.numeric(readline(""))
          while ((is.na(char))) {
            cat((crayon::blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n")))
            char <- as.numeric(readline(""))
          }
          primaryStudies <- ctmaPrep(selectedStudies=1:char)
        } # END else (catch failed)
      } # END else
    } #


    if (is.null(n.latent)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("Number of variables (n.latent) not specified!", sep="\n"))
      stop("Good luck for the next try!")
    }


    if (is.null(activeDirectory)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("No working directory has been specified!", sep="\n"))
      stop("Good luck for the next try!")
    }

    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
    }

    if (coresToUse >= parallel::detectCores()) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- parallel::detectCores() - 1
      cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
    }

    if (n.manifest > n.latent ) {
      if (is.null(lambda)) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
        cat(crayon::red("Manifest variables specified, but not matrix of loadings (lambda) specified, which is required.","\n"))
        stop("Good luck for the next try!")
      } else {
        if ( (dim(lambda)[1] != n.manifest) | (dim(lambda)[2] != n.latent) ) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
          cat(crayon::red("Dimensions of loadings matrix (lambda) do not match n.latent and n.manifest.","\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    { # fitting params
      if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
      #if (!(is.null(scaleMod))) CoTiMAStanctArgs$scaleMod <- scaleMod
      if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
      if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
      if (!(is.null(nopriors))) CoTiMAStanctArgs$nopriors <- nopriors
      if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
      if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
      if (!(is.null(iter))) CoTiMAStanctArgs$iter <- iter
      if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
    }

  } ### END Check Model Specification ###


  #######################################################################################################################
  ##################### Read user provided data and create list with all study information ##############################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("###### Read user provided data and create list with all study information #######"))
    print(paste0("#################################################################################"))

    studyList <- list()
    minInterval <- .0001 # replace missing values for time lags (cannot be NA because it is a definition variable)
    maxLengthModeratorVector <- 0
    loadRawDataStudyNumbers <- unlist(lapply(primaryStudies$rawData,
                                             function(extract) extract$studyNumbers)); loadRawDataStudyNumbers

    # resorting primary study information
    if (!(exists("moderatorNumber"))) moderatorNumber <- 1; moderatorNumber
    # determine number of studies (if list is created by ctmaInit.R it is 1 element too long)
    # Option 1: Applies for older versions of ctmaPrep
    tmp <- length(unlist(primaryStudies$studyNumbers)); tmp
    if  ( is.na(primaryStudies$deltas[tmp]) &
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
    # Option 2: versions of ctmaPrep after 14. August 2020
    if (!(is.null(primaryStudies$n.studies))) n.studies <- primaryStudies$n.studies; n.studies
    primaryStudies
    # delete empty list entries
    tmp1 <- which(names(primaryStudies) == "n.studies"); tmp1
    for (i in 1:(tmp1-1)) primaryStudies[[i]][[n.studies+1]] <- NULL


    for (i in 1:n.studies) {
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
          cat(crayon::red$bold("The correlation matrix of study no.", i, "is not symmetric. Check and re-start!", "\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    ### create pseudo raw data for all studies or load raw data if available & specified
    empraw <- lags <- moderators <- emprawMod <- allSampleSizes <- lostN <- overallNDiff <- relativeNDiff <- list()
    emprawLong <- list()
    if (is.na(primaryStudies$deltas[length(primaryStudies$deltas)])) primaryStudies$deltas <- primaryStudies$deltas[-length(primaryStudies$deltas)]
    allTpoints <-unlist(lapply(primaryStudies$deltas, function(extract) length(extract)+1)); allTpoints # more complicated thasn expected

    # ADDED TO DEAL WITH OBSERVED INDICATORS
    if (n.manifest > n.latent) {
      manifestNames <- paste0("y", 1:n.manifest); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    } else {
      manifestNames <- paste0("V", 1:n.latent); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    }
    #manifestNames; latentNames

    for (i in 1:n.studies) {
      if (!(i %in% loadRawDataStudyNumbers)) {
        currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        currentEmpcov <- (lapply(studyList, function(extract) extract$empcov))[[i]]; currentEmpcov
        currentLags <- (lapply(studyList, function(extract) extract$delta_t))[[i]]; currentLags
        currentPairwiseN <- (lapply(studyList, function(extract) extract$pairwiseN))[[i]]; currentPairwiseN
        currentModerators <- (lapply(studyList, function(extract) extract$moderators))[[i]]; currentModerators

        currentVarnames <- c()
        for (j in 1:(currentTpoints)) {
          if (n.manifest == 0) {
            for (h in 1:n.latent) {
              currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
            }
          } else {
            for (h in 1:n.manifest) {
              currentVarnames <- c(currentVarnames, paste0("y",h,"_T", (j-1)))
            }
          }
        }
        #currentVarnames

        # CORRECT FUNCTION NAME USED
        tmp <- suppressWarnings(ctmaPRaw(empCovMat=currentEmpcov, empN=currentSampleSize, empNMat=currentPairwiseN))

        empraw[[i]] <- tmp$data
        lostN[[i]] <- tmp$lostN
        overallNDiff[[i]] <- tmp$overallLostN
        relativeNDiff[[i]] <- tmp$relativeLostN
        lags[[i]] <- matrix(currentLags, nrow=dim(empraw[[i]])[1], ncol=currentTpoints-1, byrow=TRUE)
        empraw[[i]] <- cbind(empraw[[i]], lags[[i]]);
        #utils::head(empraw[[i]])
        #c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))
        colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))
        empraw[[i]] <- as.data.frame(empraw[[i]])

        # ADDED TO DEAL WITH MANIFEST VARIABLES
        n.var <- max(c(n.manifest, n.latent)); n.var
        ## START correction of current lags if entire time point is missing for a case
        emprawLongTmp <- ctsem::ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames)
        emprawLongTmp <- suppressMessages(ctsem::ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))
        ## eliminate rows where ALL latents are NA
        if (n.manifest > n.latent) targetCols <- paste0("y", 1:n.manifest) else targetCols <- paste0("V", 1:n.latent)
        emprawLongTmp <- emprawLongTmp[, ][ apply(emprawLongTmp[, targetCols], 1, function(x) sum(is.na(x)) != n.var ), ]
        # eliminate rows where time is NA
        emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctsem::ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
        # intervalise
        emprawWide <- suppressMessages(ctsem::ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames))
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
        tmpData <- utils::read.table(file=tmp1,
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
        currentTpoints <- (tmp1 + 1)/(n.var+1); currentTpoints
        if (n.manifest > n.latent) {
          colnames(empraw[[i]])[1:(currentTpoints * n.manifest)] <- paste0(paste0("x", 1:n.manifest), "_T", rep(0:(currentTpoints-1), each=n.manifest))
        } else {
          colnames(empraw[[i]])[1:(currentTpoints * n.latent)] <- paste0(paste0("V", 1:n.latent), "_T", rep(0:(currentTpoints-1), each=n.latent))
        }
        # wide to long
        emprawLongTmp <- ctsem::ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames)
        emprawLongTmp <- suppressMessages(ctsem::ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        if (n.manifest > n.latent) {
          emprawLongTmp <- emprawLongTmp[apply(emprawLongTmp[, paste0("x", 1:n.manifest)], 1, function(x) sum(is.na(x)) != n.manifest ), ]
        } else {
          emprawLongTmp <- emprawLongTmp[apply(emprawLongTmp[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        }
        # eliminate rows where time is NA
        emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctsem::ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
        # intervalise
        emprawWide <- suppressMessages(ctsem::ctIntervalise(emprawWide,
                                                            Tpoints=currentTpoints,
                                                            n.manifest=n.var,
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
        utils::write.table(empraw[[i]], file=x1, row.names=saveRawData$row.names, col.names=saveRawData$col.names,
                           sep=saveRawData$sep, dec=saveRawData$dec)
      }

      # augment pseudo raw data for stanct model
      {
        dataTmp <- empraw[[i]]
        dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=currentTpoints, n.manifest=n.var, #n.TIpred = (n.studies-1),
                                        manifestNames=manifestNames)
        dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
        dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
        # eliminate rows where ALL latents are NA
        if (n.manifest > n.latent) {
          dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("y", 1:n.manifest)], 1, function(x) sum(is.na(x)) != n.manifest ), ]
        } else {
          dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
        }
        emprawLong[[i]] <- dataTmp3
      }

    } ### END for i ...
  } ### END Read user provided data and create list with all study information ###

  # Check if sample sizes specified in prep file deviate from cases provided in possible raw data files
  N1 <- sum(unlist((lapply(empraw, function(extract) dim(extract)[1]))) , na.rm=TRUE); N1
  N2 <- sum(unlist(primaryStudies$sampleSizes), na.rm=TRUE); N2
  if (!(N1 == N2)) {
    tmp1 <- unlist(lapply(empraw, function(extract) dim(extract)[1])); tmp1
    tmp2 <- unlist(primaryStudies$sampleSizes); tmp2
    if (!(any(is.na(tmp2)))) {   # check if mismatch is because >= 1 study used pairwise N
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold(" There is a possible mismatch between sample sizes specified in the primary study list
    (created with the PREP R-file) and the cases pwovided in raw data files.", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold("N based on raw data:   ", sep="\n"))
      print(tmp1)
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold("N as specified in list:", sep="\n"))
      print(tmp2)
    }
  }


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
    overallSampleSize <- N1; overallSampleSize
    meanSampleSize <- mean(allSampleSizes, na.rm=TRUE); meanSampleSize
    maxSampleSize <- max(allSampleSizes, na.rm=TRUE); maxSampleSize
    minSampleSize <- min(allSampleSizes, na.rm=TRUE); minSampleSize
    # lags (will be computed again below if raw data are loaded #
    allDeltas <-unlist(lapply(studyList, function(extract) extract$delta_t)); allDeltas
    if (is.na(allDeltas[length(allDeltas)])) allDeltas <- allDeltas[-length(allDeltas)]; allDeltas
    meanDelta <- mean(allDeltas, na.rm=TRUE); meanDelta
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    minDelta <- min(allDeltas, na.rm=TRUE); minDelta
    # Time points
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
                           allTpoints=allTpoints, overallTpoints=overallTpoints, meanTpoints=meanTpoints,
                           maxTpoints=maxTpoints, minTpoints=minTpoints,
                           n.studies=n.studies)

  } ### END Some Statistics ###

  #######################################################################################################################
  ############################# Create ctsem model template to fit all primary studies ##################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("############ Create ctsem Model Template to fit all Primary Studies #############"))
    print(paste0("#################################################################################"))

    # allow user-specified drift matrix
    driftNames <- c()
      for (i in 1:(n.latent)) {
        for (j in 1:(n.latent)) {
          driftNames <- c(driftNames, paste0("V",i,"toV", j))
        }
      }
    # backup full names for labelling output later
    fullDriftNames <- driftNames
    if (!(is.null(drift))) {
      # check validity of user-provided drift names
      tmp1 <- which(c(drift) %in% driftNames); tmp1
      tmp2 <- which(c(drift) == "0"); tmp2
      if ( (length(tmp1)+length(tmp2)) != length(driftNames) ) {
        cat(crayon::red$bold("Drift names provided by user do not match requirements.", sep="\n"))
        cat(crayon::red$bold(" ", " ", sep="\n"))
        cat(crayon::blue("They should be of the type V1toV2 or just 0.", "\n"))
        cat(crayon::red$bold(" ", " ", sep="\n"))
        stop("Good luck for the next try!")

      } else {
      driftNames <- t(drift)
      }
      driftNames
    }

    # Adaptations if latent variables are measured with multiple indicators
    # CHD
    # loadings
    if (n.manifest > n.latent) {
      LAMBDA <- lambda
    } else {
      LAMBDA=diag(n.latent)
    }
    #LAMBDA

    # error variances
    if(!(is.null(manifestVars))) manifestVarPattern <- manifestVars else manifestVarPattern <- 0

    # T0 variance
    T0VAR <- "auto"
    skip <- 0
    if (skip == 1) {
      tmp1 <- which(LAMBDA == "0")
      tmp2 <- which(LAMBDA == "1")
      if ( (length(tmp1) + length(tmp2)) < n.var * n.latent ) {
        tmp3 <- suppressWarnings(matrix(as.numeric(LAMBDA), nrow=nrow(LAMBDA))); tmp3
        tmp3 <- as.data.frame(tmp3)
        targetVar <- which(is.na(colSums(tmp3))); targetVar
        T0VAR <- matrix(0, n.latent, n.latent); T0VAR
        for (k in 1:n.latent) {
          for (m in k:n.latent) {
            T0VAR[k, m] <- paste0("T0VAR", k, m)
            if ( (k == targetVar) & (k == m) ) T0VAR[k, m] <- 1
          }
        }
        T0VAR <- t(T0VAR)
      } else {
        T0VAR <- "auto"
      }
      T0VAR
    }

    # manifest means (as interindividually varying params, which replaces error auto-correlations)
    MANIFESTMEANS <- 0
    skip <- 0
    if (skip == 1) {
      if ( (length(tmp1) + length(tmp2)) < n.var * n.latent ) {
        MANIFESTMEANS <- rep("0", n.manifest); MANIFESTMEANS
        targetVar <- which(is.na(rowSums(tmp3))); targetVar
        MANIFESTMEANS[targetVar] <- paste0("mean_", targetVar); MANIFESTMEANS
      } else {
        MANIFESTMEANS <- 0
      }
    }

    # general ctsem model template
    ctsemModelTemplate <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=2, manifestNames=manifestNames,    # 2 waves in the template only
                                         DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent),
                                         LAMBDA=LAMBDA,
                                         T0VAR=T0VAR,
                                         type='stanct',
                                         CINT=matrix(0, nrow=n.latent, ncol=1),
                                         T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                         MANIFESTMEANS = matrix(MANIFESTMEANS, nrow = n.var, ncol = 1),
                                         MANIFESTVAR=matrix(manifestVarPattern, nrow=n.var, ncol=n.var)
    )

    if (indVarying == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
      print(paste0("#################################################################################"))

      MANIFESTMEANS <- paste0("mean_", manifestNames); MANIFESTMEANS # if provided, indVarying is the default

      ctsemModelTemplate <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=2, manifestNames=manifestNames,    # 2 waves in the template only
                                           DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent),
                                           LAMBDA=LAMBDA,
                                           T0VAR=T0VAR,
                                           type='stanct',
                                           CINT=matrix(0, nrow=n.latent, ncol=1),
                                           T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                           MANIFESTMEANS = matrix(MANIFESTMEANS, nrow = n.var, ncol = 1),
                                           MANIFESTVAR=matrix(manifestVarPattern, nrow=n.var, ncol=n.var)
      )


    }

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
      if ((activateRPB==TRUE) &  (silentOverwrite==FALSE)) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      if (silentOverwrite==FALSE) {
        cat(crayon::red$bold("You have indicated that you want to save SingleStudyModelFits, but have not selected any study/studies to save the fit for!","\n"))
        cat(crayon::red("Would you like to save the SingleStudyModelFits for ALL studies??","\n"))
        cat(crayon::blue("Press 'y' to save ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(crayon::blue("select the study/studies you whish to save the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((crayon::blue("Please press 'y' to save ALL, 's' to specify the study/studies to save, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in unlist(primaryStudies$studyNumber)) saveSingleStudyModelFit <- c(saveSingleStudyModelFit, i)
        } else if (char == 's' | char == 'S') {
          cat(crayon::blue("Which SingleStudyModelFits would you like to save?", "\n"))
          cat(crayon::blue("Please enter the no. of study/studies separated by commas!", "\n"))
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
      if ((activateRPB==TRUE) & (silentOverwrite==FALSE)) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      if (silentOverwrite==FALSE) {
        cat(crayon::red$bold("You have indicated that you want to load SingleStudyModelFits, but have not selected any study/studies to load the fit for!","\n"))
        cat(crayon::red("Would you like to load the SingleStudyModelFits for ALL studies??","\n"))
        cat(crayon::blue("Press 'y' to load ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(crayon::blue("select the study/studies you whish to load the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((crayon::blue("Please press 'y' to load ALL, 's' to specify the study/studies to load, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in unlist(primaryStudies$studyNumber)) loadSingleStudyModelFit <- c(loadSingleStudyModelFit, i)
        } else if (char == 's' | char == 'S') {
          cat(crayon::blue("Which SingleStudyModelFits would you like to load?", "\n"))
          cat(crayon::blue("Please enter the no. of study/studies separated by commas!", "\n"))
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
    model_popsd <- list()
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
        results <- ctsem::ctStanFit(
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
          verbose=verbose,
          warmup=CoTiMAStanctArgs$warmup,
          cores=coresToUse)

        studyFit[[i]] <- results
        studyFit[[i]]$resultsSummary <- summary(studyFit[[i]])

        n.par.first.lag <- ((2 * n.latent) * (2 * n.latent + 1)) / 2; n.par.first.lag
        n.par.later.lag <- ((2 * n.latent) * (2 * n.latent - 1)) / 2; n.par.later.lag
        n.later.lags <- currentTpoints - 2; n.later.lags
        df <- n.par.first.lag + sum(n.later.lags * n.par.later.lag) - studyFit[[i]]$resultsSummary$npars; df
        studyFit[[i]]$resultsSummary$'df (CoTiMA)' <- df
      } # END if (!(studyList[[i]]$originalStudyNo %in% ...

      # SAVE
      if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
        x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/"); x2
        ctmaSaveFile(activateRPB, activeDirectory, studyFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
      }

      resultsSummary <- studyFit[[i]]$resultsSummary; resultsSummary

      studyFit_Minus2LogLikelihood[[i]] <- -2 * resultsSummary$loglik
      studyFit_estimatedParameters[[i]] <- resultsSummary$npars

      tmp <- grep("toV", rownames(resultsSummary$popmeans)); tmp
      model_Drift_Coef[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"], n.latent, byrow=FALSE)); model_Drift_Coef[[i]]
      names(model_Drift_Coef[[i]]) <- c(fullDriftNames); model_Drift_Coef[[i]]

      model_Drift_SE[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"], n.latent, byrow=FALSE)); model_Drift_SE[[i]]
      names(model_Drift_SE[[i]]) <- c(fullDriftNames); model_Drift_SE[[i]]

      tmp1 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=FALSE)); tmp1
      tmp2 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=FALSE)); tmp2
      model_Drift_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Drift_CI[[i]]
      tmp3 <- c(rbind(paste0(fullDriftNames, "LL"),
                      paste0(fullDriftNames, "UL"))); tmp3
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

      if (indVarying == TRUE) model_popsd[[i]] <- resultsSummary$popsd

    } # END     for (i in 1:n.studies)

    # Combine summary information and fit statistics
    allStudies_Minus2LogLikelihood <- sum(unlist(studyFit_Minus2LogLikelihood)); allStudies_Minus2LogLikelihood
    allStudies_estimatedParameters <- sum(unlist(studyFit_estimatedParameters)); allStudies_estimatedParameters
    allStudies_df <- sum(unlist(lapply(studyFit, function(extract) extract$resultsSummary$`df (CoTiMA)`)))
    allStudies_df <- NULL
    allStudiesDRIFT_effects <- matrix(t(cbind(unlist(model_Drift_Coef), unlist(model_Drift_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    #tmp1 <- names(model_Drift_Coef[[1]]); tmp1
    tmp1 <- fullDriftNames
    tmp2 <- rep("SE", length(tmp1)); tmp2
    colnames(allStudiesDRIFT_effects) <- c(rbind(tmp1, tmp2))

    source <- lapply(primaryStudies$source, function(extract) paste(extract, collapse=", ")); source
    #source <- source[-length(source)]
    for (l in 1:length(source)) if ( source[[l]] == "NA") source[[l]] <- "Reference not provided"
    allStudiesDRIFT_effects_ext <- cbind(unlist(source), allStudiesDRIFT_effects)
    tmp <- allStudiesDRIFT_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesDRIFT_effects_ext <- tmp

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

    # check single study results
    if (checkSingleStudyResults == TRUE) {
      print(allStudiesDRIFT_effects_ext)
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::blue(" Press 'q' to quit or any other key to continue. Press ENTER afterwards."))
      char <- readline(" ")
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    DRIFTCoeff <- matrix(unlist(model_Drift_Coef), n.studies, n.latent^2, byrow=TRUE); DRIFTCoeff
    DRIFTSE <- matrix(unlist(model_Drift_SE), n.studies, n.latent^2, byrow=TRUE); DRIFTSE

    if (n.studies < 2) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::blue("Only a single primary study was handed over to ctmaInitFit. No further (meta-) analyses can be conducted."))
      cat(" ", sep="\n")
      cat(crayon::blue("I guess this stop is intended! You could ignore further warning messages such as \"sqrt(n-2): NaNs produced\""))
      cat(" ", sep="\n")
      cat(" ", sep="\n")
      cat(" ", sep="\n")
      stop("No further errors!")
      cat(" ", sep="\n")

    }

  } ### END fit ctsem model to each primary study


  ##############################################################################################################
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

  results <- list(activeDirectory=activeDirectory,
                  time=list(start.time=start.time, end.time=end.time, time.taken=time.taken),
                  plot.type="drift", model.type="stanct",
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  n.manifest=n.manifest,
                  primaryStudyList=primaryStudies,
                  studyList=studyList, studyFitList=studyFit,
                  emprawList=empraw, statisticsList=statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=list(DRIFT=names(model_Drift_Coef[[1]]), DIFFUSION=names(model_Diffusion_Coef[[1]]), T0VAR=names(model_T0var_Coef[[1]])),
                  summary=(list(model="all drift free (het. model)",
                                estimates=allStudiesDRIFT_effects_ext,
                                randomEffects=model_popsd,
                                confidenceIntervals=allStudiesCI,
                                minus2ll= round(allStudies_Minus2LogLikelihood, digits),
                                n.parameters = round(allStudies_estimatedParameters, digits),
                                #df= c(round(allStudies_df, digits))),
                                df= NULL )))
  class(results) <- "CoTiMAFit"

  ### prepare Excel Workbook with several sheets ################################################################
  {
  wb <- openxlsx::createWorkbook()
  sheet1 <- openxlsx::addWorksheet(wb, sheetName="model")
  sheet2 <- openxlsx::addWorksheet(wb, sheetName="modelResults")
  sheet3 <- openxlsx::addWorksheet(wb, sheetName="estimates")
  sheet4 <- openxlsx::addWorksheet(wb, sheetName="confidenceIntervals")
  sheet5 <- openxlsx::addWorksheet(wb, sheetName="randomEffects")
  sheet6 <- openxlsx::addWorksheet(wb, sheetName="stats")
  openxlsx::writeData(wb, sheet1, results$summary$model)

  ### modelResults
  # DRIFT
  startCol <- 2; startCol
  startRow <- 1; startRow
  openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, matrix(driftNames, nrow=1), colNames = FALSE)
  startCol <- 2; startCol
  startRow <- 2; startRow
  openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                      matrix(unlist(results$modelResults$DRIFT),
                             nrow=n.studies, ncol=n.latent^2, byrow=TRUE))
  startCol <- 1; startCol
  startRow <- 2; startRow
  openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                      matrix(unlist(primaryStudies$studyNumbers), ncol=1))
  # DIFFUSION
  offset <-  n.studies + 1
  startCol <- 2; startCol
  startRow <- 2 + offset; startRow
  openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow,
                      matrix(names(results$modelResults$DIFFUSION[[1]]), nrow=1), colNames = FALSE)
  startCol <- 2; startCol
  startRow <- 2 + offset + 1# offset; startRow
  openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                      matrix(unlist(results$modelResults$DIFFUSION),
                             nrow=n.studies, byrow=TRUE))
  startCol <- 1; startCol
  startRow <- 2 + offset + 1; startRow
  openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                      matrix(unlist(primaryStudies$studyNumbers), ncol=1))
  # T0Var
    offset <-  offset + n.studies + 2
    startCol <- 2; startCol
    startRow <- 2 + offset; startRow
    openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow,
                        matrix(names(results$modelResults$T0VAR[[1]]), nrow=1), colNames = FALSE)
    startCol <- 2; startCol
    startRow <- 2 + offset + 1# offset; startRow
    openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                        matrix(unlist(results$modelResults$T0VAR),
                               nrow=n.studies, byrow=TRUE))
    startCol <- 1; startCol
    startRow <- 2 + offset + 1; startRow
    openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                        matrix(unlist(primaryStudies$studyNumbers), ncol=1))
  ### estimates
    startCol <- 2; startCol
    startRow <- 1; startRow
    openxlsx::writeData(wb, sheet3, startCol=startCol, startRow = startRow,
                        t(colnames(results$summary$estimates)), colNames = FALSE)
    openxlsx::writeData(wb, sheet3, startCol=startCol, startRow = startRow + 1, results$summary$estimates, colNames = FALSE)
  ### confidence Intervals
    startCol <- 2; startCol
    startRow <- 1; startRow
    openxlsx::writeData(wb, sheet4, startCol=startCol, startRow = startRow,
                        t(colnames(results$summary$confidenceIntervals)), colNames = FALSE)
    openxlsx::writeData(wb, sheet4, startCol=startCol, startRow = startRow + 1, results$summary$confidenceIntervals, colNames = FALSE)
  ### random Effects
    startCol <- 2; startCol
    startRow <- 1; startRow
    openxlsx::writeData(wb, sheet5, startCol=startCol, startRow = startRow, results$summary$randomEffects, colNames = FALSE)
  ### stats
    startCol <- 2; startCol
    startRow <- 1; startRow
    tmp <- cbind("-2ll = ", results$summary$minus2ll, "Number of Parameters = ", results$summary$n.parameters)
    openxlsx::writeData(wb, sheet6, startCol=startCol, startRow = startRow, t(tmp), colNames = FALSE)
}

  results$excelSheets <- wb

  invisible(results)

}   # end function definition
