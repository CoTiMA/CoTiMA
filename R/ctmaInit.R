#' ctmaInit
#'
#' @description Fits ctsem models to each primary study in the supplied list of primary studies prepared by \code{\link{ctmaPrep}}.
#'
#' @param activateRPB set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param activeDirectory defines another active directory than the one used in \code{\link{ctmaPrep}}
#' @param binaries which manifest is a binary. Still experimental
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param checkSingleStudyResults Displays estimates from single study ctsem models and waits for user input to continue. Useful to check estimates before they are saved.
#' @param cint default 'auto' (= 0). Are set free if random intercepts model with varying cints is requested (by indVarying='cint')
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function
#' @param customPar logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param diff labels for diffusion effects. Have to be either of the character strings of the type "diff_eta1" or "diff_eta2_eta1" (= freely estimated) or values (e.g., 0 for effects to be excluded, which is usually not recommended)
#' @param digits number of digits used for rounding (in outputs)
#' @param doPar parallel and multiple fitting if single studies. A value > 1 will fit each study doPar times in parallel mode during which no output is generated (screen remains silent). Useful to obtain best fit.
#' @param drift labels for drift effects. Have to be either of the character strings of the type V1toV2 (= freely estimated) or values (e.g., 0 for effects to be excluded, which is usually not recommended)
#' @param experimental used for debugging puposes (default = FALSE)
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param fit TRUE (default) fits the requested model. FALSE returns the \code{\link{ctsem}} code CoTiMA uses to set up the model, the ctsemmodelbase which can be modified to match users requirements, and the data set (in long format created). The model can then be fitted using \code{\link{ctStanFit}})
#' @param indVarying control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means
#' @param indVaryingT0 deprecated. Automatically set to NULL.
#' @param iter number of interation (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param loadSingleStudyModelFit load the fit of single study ctsem models
#' @param manifestMeans Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param manifestVars define the error variances of the manifests within a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param n.latent number of latent variables of the model (hast to be specified)!
#' @param n.manifest number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).
#' @param optimize if set to FALSE, Stan's Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}}
#' @param priors if FALSE, any priors are disabled – sometimes desirable for optimization
#' @param randomIntercepts (default = FALSE) Experimental. Overrides ctsem's default mode for modelling indVarying cints.
#' @param sameInitialTimes Only important for raw data. If TRUE (default=FALSE), T0MEANS occurs for every subject at the same time, rather than just at the earliest observation.
#' @param saveRawData save (created pseudo) raw date. List: saveRawData$studyNumbers, $fileName, $row.names, col.names, $sep, $dec
#' @param saveSingleStudyModelFit save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#' @param scaleTI scale TI predictors
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param silentOverwrite overwrite old files without asking
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param T0var (default = 'auto')
#' @param useSV if TRUE (default=FALSE) start values will be used if provided in the list of primary studies
#' @param verbose integer from 0 to 2. Higher values print more information during model fit - for debugging

#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red blue
#' @importFrom parallel detectCores
#' @importFrom ctsem ctDeintervalise ctLongToWide ctIntervalise ctWideToLong ctModel ctStanFit ctExtract ctCollapse
#' @importFrom utils read.table write.table
#' @importFrom openxlsx addWorksheet writeData createWorkbook openXL saveWorkbook
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom foreach %dopar%
#' @importFrom stats cov2cor
#' @importFrom OpenMx expm
#'
#' @export ctmaInit
#'
#' @examples
#' # Fit a ctsem model to all three primary studies summarized in
#' # CoTiMAstudyList_3 and save the three fitted models
#' \dontrun{
#' CoTiMAInitFit_3 <- ctmaInit(primaryStudies=CoTiMAstudyList_3,
#'                             n.latent=2,
#'                             checkSingleStudyResults=FALSE,
#'                             activeDirectory="/Users/tmp/") # adapt!
#' summary(CoTiMAInitFit_3)
#' }
#'
#' @return ctmaFit returns a list containing some arguments supplied, the fitted models, different elements summarizing the main results,
#' model type, and the type of plot that could be performed with the returned object. The arguments in the returned object are activeDirectory,
#' coresToUse, n.latent, n.manifest, and primaryStudyList. The study count is returned as n.studies, the created matrix of loadings of
#' manifest on latent factors is returned as lambda, and a re-organized list of primary studies with some information ommited is returned as
#' studyList. The fitted models for each primary study are found in studyFitList, which is a large list with many elements (e.g., the ctsem
#' model specified by CoTiMA, the rstan model created by ctsem, the fitted rstan model etc.). Further results returned are emprawList
#' (containing the pseudo raw data created), statisticsList (comprising baisc stats such as average sample size, no. of measurement points,
#' etc.), a list with modelResults (i.e., DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef,
#' CINT=model_Cint_Coef), and the paramter names internally used. The summary list,  which is printed if the summary function is applied to the
#' returned object, comprises "estimates" (the aggregated effects), possible randomIntercepts,confidenceIntervals, the
#' minus2ll value and its n.parameters, and possible warning messages (message). Plot type is plot.type=c("drift") and model.type="stanct"
#' ("omx" was deprecated).
#'
ctmaInit <- function(
    activateRPB=FALSE,
    activeDirectory=NULL,
    binaries=NULL,
    chains=NULL,
    checkSingleStudyResults=FALSE,
    cint=0,
    coresToUse=c(2),
    CoTiMAStanctArgs=NULL,
    customPar=FALSE,
    diff=NULL,
    digits=4,
    doPar=1,
    drift=NULL,
    experimental=FALSE,
    finishsamples=NULL,
    fit=TRUE,
    indVarying=FALSE,
    indVaryingT0=NULL,
    iter=NULL,
    lambda=NULL,
    loadSingleStudyModelFit=c(),
    manifestMeans=0,
    manifestVars=NULL,
    n.latent=NULL,
    n.manifest=0,
    optimize=TRUE,
    primaryStudies=NULL,
    priors=FALSE,
    randomIntercepts=FALSE,
    sameInitialTimes=FALSE,
    saveRawData=list(),
    saveSingleStudyModelFit=c(),
    scaleTI=NULL,
    scaleTime=NULL,
    silentOverwrite=FALSE,
    T0means=0,
    T0var='auto',
    useSV=FALSE,
    verbose=0
)

{  # begin function definition (until end of file)


  #start.time <- Sys.time()

  original.options <- options("scipen"); original.options
  options(scipen = 999); options("scipen") # turn scientific notation off.
  on.exit(options(scipen = original.options))  # scientific notation as user's original


  #######################################################################################################################
  ########################################### Check Model Specification #################################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n########################## Check Model Specification ############################ \n#################################################################################"
    message(Msg)


    if (!(is.null(indVaryingT0))) {
      Msg <- "The argument \"indVaryingT0\" was deprecated and set to NULL. Try \"indVarying = TRUE\" or \"randomIntercepts = TRUE\" instead.\n"
      message(Msg)
    }

    if ( (indVarying == "CINT") | (indVarying == "Cint") | (indVarying == "cint")) indVarying <- "CINT"
    if ( (indVarying == "MANIFEST") | (indVarying == "Manifest") | (indVarying == "manifest")) indVarying <- TRUE
    #
    randomInterceptsSettings <- randomIntercepts
    if (is.null(randomIntercepts)) randomIntercepts <- FALSE
    if ( (randomIntercepts == 'CINT') | (randomIntercepts == 'cint')  | (randomIntercepts == 'Cint')) randomIntercepts <- TRUE
    if ( (randomIntercepts == "Manifest") | (randomIntercepts == "manifest") | (randomIntercepts == "MANIFEST")) randomIntercepts <- "MANIFEST"
    if ( (randomIntercepts == "MANIFEST") | (randomIntercepts == TRUE) ) {
      indVarying <- FALSE
      indVaryingT0 <- NULL
    }

    #if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMA::CoTiMAStanctArgs$verbose}
    if (is.null(verbose)) tmp1 <- 1 else tmp1 <- 0
    if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0}
    if (tmp1 == 1) {verbose <- CoTiMA::CoTiMAStanctArgs$verbose}

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
        ErrorMsg <- "\nGood luck for the next try!"
        stop(ErrorMsg)
      } else {
        if (!(is.numeric(tryCatch(get(paste("sampleSize", 1, sep = "")), error = function(e) e)))) {
          ErrorMsg <- "\nGetting primary study information from the global environment failed \nTo test, I searched for sampleSize1 and could not find it! \nGood luck for the next try!"
          stop(ErrorMsg)
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
      ErrorMsg <- "\nNumber of variables (n.latent) not specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (is.null(binaries)) binaries <- rep(0, max(n.manifest, n.latent))
    if (length(binaries) != max(n.manifest, n.latent)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nThe number of binaries provided is incorrect! \nGood luck for the next try!"
      stop(ErrorMsg)
    }


    if (is.null(activeDirectory)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nNo active directory has been specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
    }

    if (coresToUse >= parallel::detectCores()) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- parallel::detectCores() - 1
      Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.\n"
      message(Msg)
    }

    # CHD added Aug 2023 because on github nopriors was replaced by priors argument
    #tmp1 <- formals(ctsem::ctStanFit)
    #if (is.na(tmp1$nopriors)) {
    #  nopriors <-NA
    #  CoTiMAStanctArgs$nopriors <- NA
    #}

    # CHD changed on 19 Sep 2022
    if (doPar > 1) {
      #doParallel::registerDoParallel(detectCores()-1)
      myCluster <- parallel::makeCluster(coresToUse)
      on.exit(parallel::stopCluster(myCluster))
      #doParallel::registerDoParallel(coresToUse)
      doParallel::registerDoParallel(myCluster)
      '%dopar%' <- foreach::'%dopar%'
    }

    if (n.manifest > n.latent ) {
      if (is.null(lambda)) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
        ErrorMsg <- "\nManifest variables specified, but not matrix of loadings (lambda) specified, which is required. \nGood luck for the next try!"
        stop(ErrorMsg)
      } else {
        if ( (dim(lambda)[1] != n.manifest) | (dim(lambda)[2] != n.latent) ) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
          ErrorMsg <- "\nDimensions of loadings matrix (lambda) do not match n.latent and n.manifest. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
    }

    if ((n.manifest <= n.latent ) &  (is.null(lambda)))  lambda=diag(n.latent)


    { # fitting params

      invariantDrift <- FALSE
      moderatedDrift <- NULL

      # Added 17. Aug 2022
      tmp1 <- names(CoTiMA::CoTiMAStanctArgs) %in% names(CoTiMAStanctArgs); tmp1
      tmp2 <- CoTiMA::CoTiMAStanctArgs
      if (!(is.null(CoTiMAStanctArgs))) tmp2[tmp1] <- CoTiMAStanctArgs
      CoTiMAStanctArgs <- tmp2

      if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
      if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
      if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
      #if ( (!(is.null(nopriors))) & (!(is.null(nopriors))) ) CoTiMAStanctArgs$nopriors <- nopriors # changed Aug 2023
      if (!(is.null(priors))) CoTiMAStanctArgs$priors <- priors # added Aug 2023
      if (!(is.null(priors))) CoTiMAStanctArgs$priors <- priors # added Aug 2023
      if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
      if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
      if (!(is.null(iter))) CoTiMAStanctArgs$iter <- iter
      if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
    }

    # correction of misspecified user input
    if (is.null(T0means)) T0means <- 0
    if (is.null(manifestMeans)) manifestMeans <- 0


  } ### END Check Model Specification ###


  #######################################################################################################################
  ##################### Read user provided data and create list with all study information ##############################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n###### Read user provided data and create list with all study information ####### \n#################################################################################"
    message(Msg)

    studyList <- list()
    minInterval <- .0001 # replace missing values for time lags (cannot be NA because it is a definition variable)
    maxLengthModeratorVector <- 0
    loadRawDataStudyNumbers <- unlist(lapply(primaryStudies$rawData,
                                             function(extract) extract$studyNumbers)); loadRawDataStudyNumbers

    # resorting primary study information
    if (!(exists("moderatorNumber"))) moderatorNumber <- 1; moderatorNumber
    # determine number of studies (if list is created by ctmaInit.R it is 1 element too long)
    # Option 1: Applies for older versions of ctmaPrep
    #tmp <- length(unlist(primaryStudies$studyNumbers)); tmp
    tmp <- unlist(primaryStudies$studyNumbers); tmp
    tmp2 <- which(is.na(tmp)); tmp2
    for (i in tmp2) primaryStudies$studyNumbers[[i]] <- NULL
    tmp <- length(tmp[!(is.na(tmp))]); tmp
    if  ( is.na(primaryStudies$deltas[tmp]) &
          is.na(primaryStudies$sampleSizes[tmp]) &
          is.na(primaryStudies$deltas[tmp]) &
          is.null(dim(primaryStudies$pairwiseNs[tmp])) &
          is.null(dim(primaryStudies$emcovs[tmp])) &
          all(is.na(unlist(primaryStudies$moderators[tmp]))) &
          is.null(primaryStudies$rawData$fileName[tmp])
          # 9. Aug. 2022:
          & (length(tmp) > 1) ) {
      if (tmp != 1) n.studies <- tmp-1 else n.studies <- 1
      primaryStudies$studyNumbers[[tmp]] <- NULL
    } else {
      n.studies <- tmp
    }
    # Option 2: versions of ctmaPrep after 14. August 2020
    if (!(is.null(primaryStudies$n.studies))) n.studies <- primaryStudies$n.studies; n.studies

    # delete empty list entries
    if (n.studies > 1) { # may not apply if a single study is fitted (e.g., ctmaOptimizeInit) # RECENT CHANGE
      tmp1 <- which(names(primaryStudies) == "n.studies"); tmp1
      for (i in 1:(tmp1-1)) primaryStudies[[i]][[n.studies+1]] <- NULL
    }

    for (i in 1:n.studies) {
      studyList[[i]] <- list(studyNumber=i, empcov=primaryStudies$empcovs[[i]], delta_t=primaryStudies$deltas[[i]],
                             sampleSize=primaryStudies$sampleSizes[[i]], originalStudyNo=primaryStudies$studyNumber[[i]],
                             timePoints=sum(length(primaryStudies$deltas[[i]]), 1), moderators=primaryStudies$moderators[[i]],
                             # CHD changed next line on 29 Sep 2022
                             #maxModerators=length(primaryStudies$moderators[[i]]), startValues=primaryStudies$inits[[i]],
                             maxModerators=length(primaryStudies$moderators[[i]]), startValues=primaryStudies$startValues[[i]],
                             rawData=primaryStudies$rawData[[i]], pairwiseN=primaryStudies$pairwiseNs[[i]],
                             source=paste(primaryStudies$source[[i]], collapse=", "))
      if (useSV == FALSE) studyList[[i]]$startValues <- NULL
      # CHD added next line on 30 Sep 2022, changed 7 Oct 2022 , see also line 893 & 939
      #if ((useSV == TRUE) & (is.null(studyList[[i]]$startValues)) ) studyList[[i]]$startValues <- NA
      #if ((useSV == TRUE) & (is.na(studyList[[i]]$startValues)) ) studyList[[i]]$startValues <- NULL DOES NOT WORK
      if (is.null(studyList[[i]]$startValues)) studyList[[i]]$startValues <- NA


      if (length(primaryStudies$moderators[[i]]) > maxLengthModeratorVector) maxLengthModeratorVector <- length(primaryStudies$moderators[[i]])
      # check matrix symmetry if matrix is provided
      if (!(primaryStudies$studyNumbers[i] %in% loadRawDataStudyNumbers)) {
        if (isSymmetric(primaryStudies$empcovs[[i]]) == FALSE) {
          ErrorMsg <- paste0("\nThe correlation matrix of study no.", i, " is not symmetric. Check and re-start! \nGood luck for the next try!")
          stop(ErrorMsg)
        }
      }
    }

    ### create pseudo raw data for all studies or load raw data if available & specified
    empraw <- lags <- moderators <- emprawMod <- allSampleSizes <- lostN <- overallNDiff <- relativeNDiff <- list()
    emprawLong <- list()
    empraw.ind.mod <- list() # CHD 19.6.2023

    if (n.studies > 1) {
      if (is.na(primaryStudies$deltas[length(primaryStudies$deltas)])) primaryStudies$deltas <- primaryStudies$deltas[-length(primaryStudies$deltas)]
    }
    if (is.na(primaryStudies$deltas[length(primaryStudies$deltas)])) primaryStudies$deltas[length(primaryStudies$deltas)] <- NULL
    allTpoints <-unlist(lapply(primaryStudies$deltas, function(extract) length(extract)+1)); allTpoints # more complicated than expected

    # manifests
    if (n.manifest > n.latent) {
      manifestNames <- paste0("y", 1:n.manifest); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    } else {
      manifestNames <- paste0("V", 1:n.latent); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    }

    for (i in 1:n.studies) {
      #i <- 1
      if (!(studyList[[i]]$originalStudyNo %in% loadRawDataStudyNumbers)) {
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

        # Create Pseudo Raw Data
        if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit)) {
          tmp1 <- paste0(" Create Pseudo Raw Data for Study No. ", i, ".    Could take long !!! ")
          tmp2 <- nchar(tmp1); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
          Msg <- paste0("################################################################################# \n", tmp6, "\n#################################################################################")
          message(Msg)

          #Msg <- paste0("################################################################################# \n###### Create Pseudo Raw Data for Study No. ", i, ".    Could take long !!! ####### \n#################################################################################")
          #message(Msg)
        }


        # CHD ADDED 7.9.2022 reduce computation time by creating Pseudo Raw Data with small sample size because the data are loaded anyway
        currentSampleSizeTmp <- currentSampleSize
        currentPairwiseNTmp <- currentPairwiseN
        currentEmpcovTmp <- currentEmpcov
        if (length(loadSingleStudyModelFit) > 0) {
          tmp <- as.numeric(loadSingleStudyModelFit[2:length(loadSingleStudyModelFit)]); tmp
          if (studyList[[i]]$originalStudyNo %in% tmp) {
            currentSampleSizeTmp <- (n.latent * currentTpoints)^2; currentSampleSizeTmp
            currentPairwiseNTmp <- matrix(0,0,0)
            currentEmpcovTmp <- matrix(0, n.latent * currentTpoints, n.latent * currentTpoints); currentEmpcovTmp
          }
        }

        # CHD: changed 7.9.2022 "Tmp"
        tmp <- suppressWarnings(ctmaPRaw(empCovMat=currentEmpcovTmp, empN=currentSampleSizeTmp, empNMat=currentPairwiseNTmp,
                                         experimental=FALSE))

        empraw[[i]] <- tmp$data
        lostN[[i]] <- tmp$lostN
        overallNDiff[[i]] <- tmp$overallLostN
        relativeNDiff[[i]] <- tmp$relativeLostN
        lags[[i]] <- matrix(currentLags, nrow=dim(empraw[[i]])[1], ncol=currentTpoints-1, byrow=TRUE)
        empraw[[i]] <- cbind(empraw[[i]], lags[[i]]);
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
      # CHD 13.6.2023 changed to allow raw data that are recovered by ctmaFitToList to be used (is.na was added)
      if (studyList[[i]]$originalStudyNo %in% loadRawDataStudyNumbers) {
        # CHD 13.6.2023 changed to allow raw data that are recovered by ctmaFitToList to be used (is.na was added)

        if ( (!(is.null(primaryStudies$emprawList[[i]]))) &
             (
               ((is.null(primaryStudies$empcovs[[i]]))
                |
                (length(primaryStudies$empcovs[[i]]) <= 1) # could be NA or 0x0 matrix
               )
             )
        ) {
          # if the function list of primary studies is already post-processed (ctmaSV) and called from ctmaOptimizeINit)
          empraw[[i]] <- primaryStudies$emprawList[[i]]
          if (!(exists("n.var"))) n.var <- max(c(n.latent, n.manifest))
          tmp1 <- dim(empraw[[i]])[2]; tmp1
          currentTpoints <- (tmp1 + 1)/(n.var+1); currentTpoints
          currentTpointsBackup <- currentTpoints # in case function is called from ctmaOptimizeInit
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
          tmp1 <- empraw[[i]]
          tmp2 <- tmp1[,grep("dT", colnames(tmp1))] == minInterval
          tmp1[ ,grep("dT", colnames(tmp1))][tmp2] <- NA
          for (h in 1:(currentTpoints-1)) studyList[[i]]$delta_t[h] <- mean(tmp1[, paste0("dT", h)], na.rm=TRUE)
          #studyList[[i]]$delta_t
        } else {
          if (!(exists("n.var"))) n.var <- max(c(n.latent, n.manifest))
          currentTpoints <- length((lapply(studyList, function(extract) extract$delta_t))[[i]])+1; currentTpoints

          tmp1 <- studyList[[i]]$rawData$fileName; tmp1
          tmp2 <- studyList[[i]]$rawData$header; tmp2
          tmp3 <- studyList[[i]]$rawData$dec; tmp3
          tmp4 <- studyList[[i]]$rawData$sep; tmp4
          if (!(is.null(primaryStudies$activeDirectory))) {
            if (activeDirectory != primaryStudies$activeDirectory) {
              tmp1 <- gsub(primaryStudies$activeDirectory, activeDirectory, tmp1, fixed=TRUE)
              print(tmp1)
            }
          } else {
            #tmp1 <- gsub(primaryStudies$activeDirectory, activeDirectory, tmp1, fixed=TRUE)
          }
          tmpData <- utils::read.table(file=tmp1,
                                       header=tmp2,
                                       dec=tmp3,
                                       sep=tmp4)
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

          # CHD 19.6.2023 extract possible ind level moderators 28.6.2023
          if ("n.ind.mod" %in% names(studyList[[i]]$rawData)) {
            if (studyList[[i]]$rawData$n.ind.mod != 0) {
              if (!(exists("n.var"))) n.var <- max(c(n.latent, n.manifest))
              targetCols <- (n.var * currentTpoints + currentTpoints -1 + 1): ncol(tmpData); targetCols
              empraw.ind.mod[[i]] <- list()
              empraw.ind.mod[[i]] <- tmpData[, targetCols]
              empraw.ind.mod[[i]][empraw.ind.mod[[i]] == studyList[[i]]$rawData$missingValues] <- NA
              tmpData <- tmpData[, -targetCols]
            } else {
              empraw.ind.mod[[i]] <- NA
            }
          }
          #}


          # replace missing values
          tmpData <- as.matrix(tmpData) # important: line below will not work without having data as a matrix
          tmpData[tmpData %in% studyList[[i]]$rawData$missingValues] <- NA
          empraw[[i]] <- as.data.frame(tmpData)

          ## START correction of current lags if entire time point is missing for a case
          # if called from ctmaOptimize
          if (!(exists("n.var"))) n.var <- max(n.latent, n.manifest)
          # change variable names
          tmp1 <- dim(empraw[[i]])[2]; tmp1
          currentTpoints <- (tmp1 + 1)/(n.var+1); currentTpoints

          # CHD 19.6.2023 this is an error check (n.ind.mod set to 0 but data set includes TI at the end)
          if (currentTpoints != round(currentTpoints, 0)) {
            if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
            ErrorMsg <- "\nI have problems with the raw data set. Possibly you specified n.ind.mod incorrectly before doing ctmaPrep. The n.ind.mod should by > 0 if the raw data set does not have the last time interval as last column.\nGood luck for the next try!"
            stop(ErrorMsg)
          }

          currentTpointsBackup <- currentTpoints # in case function is called from ctmaOptimizeInit
          if (n.manifest > n.latent) {
            colnames(empraw[[i]])[1:(currentTpoints * n.manifest)] <- paste0(paste0("x", 1:n.manifest), "_T", rep(0:(currentTpoints-1), each=n.manifest))
          } else {
            colnames(empraw[[i]])[1:(currentTpoints * n.latent)] <- paste0(paste0("V", 1:n.latent), "_T", rep(0:(currentTpoints-1), each=n.latent))
          }

          # wide to long
          emprawLongTmp <- ctsem::ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames)
          emprawLongTmp <- suppressMessages(ctsem::ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))

          # eliminate rows where ALL latents are NA
          ## skipped on 31. Aug. 2022
          #if (n.manifest > n.latent) {
          #  emprawLongTmp <- emprawLongTmp[apply(emprawLongTmp[, paste0("x", 1:n.manifest)], 1, function(x) sum(is.na(x)) != n.manifest ), ]
          #} else {
          #  emprawLongTmp <- emprawLongTmp[apply(emprawLongTmp[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
          #}

          # eliminate rows where time is NA
          emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
          # 9. Aug. 2022: handle possibly reduced Tpoints (if some rows were eliminated for all cases)
          currentTpoints <- max(table(emprawLongTmp[,1])); currentTpoints

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


          # change sample size if entire cases were deleted
          studyList[[i]]$sampleSize <- (dim(empraw[[i]]))[1]
          allSampleSizes[[i]] <- dim(empraw[[i]])[1]; allSampleSizes[[i]]
          currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
          currentTpoints <- allTpoints[[i]]; currentTpoints
          if (is.null(currentTpoints)) currentTpoints <- currentTpointsBackup

          #currentVarnames
          colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))

          # standardize (variables - not time lags) if option is chosen
          if (studyList[[i]]$rawData$standardize == TRUE) empraw[[i]][, currentVarnames] <- scale(empraw[[i]][, currentVarnames])
          # CHD undone 10.11.2023
          #if ( (studyList[[i]]$rawData$standardize == TRUE) &
          #  (experimental2 == FALSE) ) empraw[[i]][, currentVarnames] <- scale(empraw[[i]][, currentVarnames])


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
        }
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

        #if (experimental2 == TRUE) {
        #  print(paste0("Showing time points avaliable for Study ", i,"."))
        #  tmp1 <- table(emprawLong[[i]][, "time"])
        #  print(tmp1)
        #}
      }
    } ### END for i ...

    # Check if sample sizes specified in prep file deviate from cases provided in possible raw data files
    N1 <- (unlist((lapply(empraw, function(extract) dim(extract)[1])))); N1
    N2 <- (unlist(primaryStudies$sampleSizes)); N2
    N1 <- sum(unlist((lapply(empraw, function(extract) dim(extract)[1]))) , na.rm=TRUE); N1
    N2 <- sum(unlist(primaryStudies$sampleSizes), na.rm=TRUE); N2
    if (!(N1 == N2)) {
      tmp1 <- unlist(lapply(empraw, function(extract) dim(extract)[1])); tmp1
      tmp2 <- unlist(primaryStudies$sampleSizes); tmp2
      if (!(any(is.na(tmp2)))) {   # check if mismatch is because >= 1 study used pairwise N
        #Msg <- paste0("There is a possible mismatch between sample sizes specified in the primary study list
        #(created with the PREP R-file) and the cases provided in raw data files.\nN based on raw data: \n", tmp1)
        #message(Msg)
        #Msg <- paste0("\nN as specified in list: \n", tmp2)
        #message(Msg)
      }
    }
  } ### END Read user provided data and create list with all study information ###


  #######################################################################################################################
  ################################################### Some Statistics ###################################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n################# Compute Summary Statistics of Primary Studies ################# \n#################################################################################"
    message(Msg)

    ### some stats
    # Sample size
    allSampleSizes <-unlist(lapply(studyList, function(extract) extract$sampleSize)); allSampleSizes
    overallSampleSize <- N1; overallSampleSize
    meanSampleSize <- mean(allSampleSizes, na.rm=TRUE); meanSampleSize
    maxSampleSize <- suppressWarnings(max(allSampleSizes, na.rm=TRUE)); maxSampleSize
    minSampleSize <- suppressWarnings(min(allSampleSizes, na.rm=TRUE)); minSampleSize
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
    Msg <- "################################################################################# \n############# Set ctsem Model Parameters to fit all Primary Studies ############# \n#################################################################################"
    message(Msg)

    namesAndParams <- ctmaLabels(
      n.latent=n.latent,
      n.manifest=n.manifest,
      lambda=lambda,
      drift=drift,
      diff=diff, # CHD 13.3.24
      invariantDrift=invariantDrift,
      moderatedDrift=NULL,
      equalDrift=NULL,
      T0means=T0means,
      manifestMeans=manifestMeans,
      manifestVars=manifestVars
    )
    #
    driftNames <- namesAndParams$driftNames; driftNames
    driftFullNames <- namesAndParams$driftFullNames; driftFullNames
    driftParams <- namesAndParams$driftParams; driftParams
    diffNames <- namesAndParams$diffNames; diffNames
    diffParams <- namesAndParams$diffParams; diffParams
    diffFullNames <- namesAndParams$diffFullNames; diffFullNames
    invariantDriftNames <- namesAndParams$invariantDriftNames; invariantDriftNames
    invariantDriftParams <- namesAndParams$invariantDriftParams; invariantDriftParams
    lambdaParams <- namesAndParams$lambdaParams; lambdaParams
    T0VARParams <- namesAndParams$T0VARParams; T0VARParams
    manifestmeansParams <- namesAndParams$manifestMeansParams; manifestmeansParams
    manifestMeansParams <- namesAndParams$manifestMeansParams; manifestMeansParams
    T0meansParams=namesAndParams$T0meansParams; T0meansParams
    manifestVarsParams <- namesAndParams$manifestVarsParams; manifestVarsParams
    #CHD 13.6.2023
    T0VARParams <- T0var; T0VARParams
    CINTParams <- cint; CINTParams

  } ### END Create ctsem model template to fit all primary studies ###


  #######################################################################################################################
  ##################################### Check Specification of Primary Studies ##########################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n#################### Check Specification of Primary Studies ##################### \n#################################################################################"
    message(Msg)

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
          ErrorMsg <- "Good luck for the next try!"
          stop(ErrorMsg)
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
          ErrorMsg <- "Good luck for the next try!"
          stop(ErrorMsg)
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
    model_popsd <- list() # model_popcov <- model_popcor <- list()
    resultsSummary <- list()
    model_popcov_m <- model_popcov_sd <- model_popcov_T <- model_popcov_025 <- model_popcov_50 <- model_popcov_975 <- list()
    model_popcor_m <- model_popcor_sd <- model_popcor_T <- model_popcor_025 <- model_popcor_50 <- model_popcor_975 <- list()
    estProb <- list()

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
        Msg <- paste0("################################################################################# \n" , tmp6 ,"\n#################################################################################")
        message(Msg)
        x1 <- paste0(activeDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        #file.exists(x1)
        if (file.exists(x1)) {
          notLoadable <- FALSE
          studyFit[[i]] <- readRDS(file=x1)
          # CHD added 21 Sep 2022
          empraw[[i]] <- studyFit[[i]]$empraw
        } else {
          notLoadable <- TRUE
        }
      } # END if ( (length(loadSingleStudyModelFit) > 1) &  ...

      tmpLogic <- 0
      if ((studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) & (notLoadable == TRUE) ) tmpLogic <- 1
      if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) tmpLogic <- 1
      if (tmpLogic == 1) {
        tmp1 <- paste0(" Fitting SingleStudyModel ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ")
        tmp2 <- nchar(tmp1); tmp2
        tmp3 <- (81 - tmp2)/2; tmp3
        tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
        tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
        tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
        Msg <- paste0("################################################################################# \n", tmp6, "\n#################################################################################")
        message(Msg)

        if (!(optimize)) {
          customPar <- FALSE
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}

          tmp1a <- paste0(" Bayesian sampling was selected, which does require appropriate scaling of time. ")
          tmp2 <- nchar(tmp1a); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6a <- paste0(tmp4, tmp1a, tmp5); tmp6

          tmp1b <- paste0(" See the end of the summary output ")
          tmp2 <- nchar(tmp1b); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6b <- paste0(tmp4, tmp1b, tmp5); tmp6

          Msg <- paste0("################################################################################# \n", tmp6a, "\n", tmp6b, "\n#################################################################################")
          message(Msg)


          #Msg <- "Bayesian sampling was selected, which does require appropriate scaling of time. See the end of the summary output \n"
          #message(Msg)

        }

        # select correct template
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        driftParamsTmp <- driftParams
        diffParamsTmp <- diffParams
        longestLag <- max((lapply(studyList, function(extract) extract$delta_t))[[i]]); longestLag
        if ((longestLag > 6)  & (customPar)) {
          counter <- 0
          for (h in 1:(n.latent)) {
            for (j in 1:(n.latent)) {
              counter <- counter + 1
              if (h == j) {
                driftParamsTmp[counter] <- paste0(driftParams[counter], paste0("|-log1p_exp(-param *.1 -2)"))
                diffParamsTmp[counter] <- paste0(diffParams[counter], paste0("|log1p_exp(-param *.1 -2)"))
              }
            }
          }
        }

        # CHD 13.6.2023
        if ((indVarying == 'cint') | (indVarying == 'Cint')) indVarying <- 'CINT'
        if ((indVarying == TRUE) & (is.null(indVaryingT0))) indVaryingT0 <- TRUE
        if ((indVarying == 'CINT') & (is.null(indVaryingT0))) indVaryingT0 <- TRUE
        if (is.null(indVaryingT0)) indVaryingT0 <- FALSE

        # CHD 9.6.2023

        #if (randomIntercepts == FALSE) {
        if ( (randomIntercepts != TRUE) & (randomIntercepts != "MANIFEST") ) {
          if ( (indVarying == 'CINT') & (indVaryingT0 == TRUE) ) {
            print(paste0("#################################################################################"))
            print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
            print(paste0("#################################################################################"))

            T0meansParams <- 'auto'

            print(paste0("#################################################################################"))
            print(paste0("######################### CT intercepts are set free.  ##########################"))
            print(paste0("#################################################################################"))

            CINTParams <- c()
            for (c in 1:n.latent) {
              CINTParams <- c(CINTParams, paste0("cintV", c))
            }
          }
          #CINTParams

          if ( (indVarying == 'CINT') & (indVaryingT0 == FALSE) ) {
            print(paste0("#################################################################################"))
            print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
            print(paste0("#################################################################################"))

            T0meansParams <- 0

            print(paste0("#################################################################################"))
            print(paste0("####################### CT intercepts are set free.  ########################"))
            print(paste0("#################################################################################"))

            CINTParams <- c()
            for (c in 1:n.latent) {
              CINTParams <- c(CINTParams, paste0("cintV", c))
            }
          }

          if ( (indVarying == TRUE) & (indVaryingT0 == TRUE) ) {
            print(paste0("#################################################################################"))
            print(paste0("###### Just a note: Individually varying manifest means model requested.  #######"))
            print(paste0("#################################################################################"))

            T0meansParams <- 'auto'

            print(paste0("#################################################################################"))
            print(paste0("######### Manifest means (as replacement for intercepts) are set free.  #########"))
            print(paste0("#################################################################################"))

            manifestMeansParams <- 'auto'
          }

          if ( (indVarying == TRUE) & (indVaryingT0 == FALSE) ) {
            print(paste0("#################################################################################"))
            print(paste0("###### Just a note: Individually varying manifest means model requested.  #######"))
            print(paste0("#################################################################################"))

            T0meansParams <- 0

            print(paste0("#################################################################################"))
            print(paste0("######### Manifest means (as replacement for intercepts) are set free.  #########"))
            print(paste0("#################################################################################"))

            manifestMeansParams <- 'auto'
          }
        } # end  if ( (randomIntercepts != TRUE) & (randomIntercepts != "MANIFEST") )

        currentModel <- suppressMessages(
          ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=currentTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                         DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                         DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                         LAMBDA=lambdaParams,
                         T0VAR=T0VARParams,
                         type='stanct',
                         CINT=matrix(CINTParams, nrow=n.latent, ncol=1),
                         T0MEANS = matrix(c(T0meansParams), nrow = n.latent, ncol = 1),
                         MANIFESTMEANS = matrix(manifestMeansParams, nrow = n.var, ncol = 1),
                         MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var)
          )
        )

        if (indVarying == FALSE) currentModel$pars[, "indvarying"] <- FALSE
        #CHD 13.6.2023
        if (indVaryingT0 == TRUE) {
          currentModel$pars[currentModel$pars$matrix %in% 'T0MEANS','indvarying'] <- TRUE
        } else {
          currentModel$pars[currentModel$pars$matrix %in% 'T0MEANS','indvarying'] <- FALSE
        }
        # CHD 13.6.2023
        if (indVarying == 'CINT') {
          currentModel$pars[currentModel$pars$matrix %in% 'CINT','indvarying'] <- TRUE
        } else {
          currentModel$pars[currentModel$pars$matrix %in% 'CINT','indvarying'] <- FALSE
        }
        if (indVarying == TRUE) {
          currentModel$pars[currentModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- TRUE
        } else {
          currentModel$pars[currentModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- FALSE
        }
        #currentModel$pars

        currentModel$manifesttype <- binaries

        #if (randomIntercepts == TRUE) { # override ctsem's default setup for indvarying cints
        if ((randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") ) {
          print(paste0("#################################################################################"))
          print(paste0("# Note: Random intercepts modelled as processes, not as cint or manifest means. #"))
          print(paste0("#################################################################################"))

          nullMat <- matrix(0, n.latent, n.latent); nullMat
          DIFFUSIONtmp <- matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent); DIFFUSIONtmp
          DIFFUSIONtmp <- rbind(cbind(DIFFUSIONtmp, nullMat),
                                cbind(nullMat, nullMat)); DIFFUSIONtmp
          DRIFTtmp <- matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent); DRIFTtmp
          DRIFTtmp <- rbind(cbind(DRIFTtmp, diag(2)),
                            cbind(nullMat, nullMat)); DRIFTtmp
          if (randomIntercepts == "MANIFEST") {
            DRIFTtmp <- matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent); DRIFTtmp
            DRIFTtmp <- rbind(cbind(DRIFTtmp, nullMat),
                              cbind(nullMat, nullMat)); DRIFTtmp
          }
          T0VARtmp = gsub("diff", "T0cov", matrix(diffParamsTmp, n.latent, n.latent)); T0VARtmp
          cint_cint_cov <- gsub("eta", "cint", T0VARtmp); cint_cint_cov
          cint_cint_cov <- gsub("T0c", "C", cint_cint_cov); cint_cint_cov
          T0eta_cint_cov <- matrix(NA, n.latent, n.latent); T0eta_cint_cov
          for (ii in 1:n.latent) {
            for (ll in 1:n.latent) {
              T0eta_cint_cov[ii, ll] <- paste0("Cov_T0eta", ll, "_cint", ii)
            }
          }
          #T0eta_cint_cov
          T0VARtmp <- rbind(cbind(T0VARtmp, nullMat),
                            cbind(T0eta_cint_cov, cint_cint_cov)); T0VARtmp
          LAMBDAtmp <- cbind(lambdaParams, nullMat); LAMBDAtmp
          if (randomIntercepts == "MANIFEST") {
            LAMBDAtmp <- cbind(lambdaParams, lambdaParams); LAMBDAtmp
          }
          manifestNamesTmp <- manifestNames; manifestNamesTmp
          latentNamesTmp <- c(latentNames, paste0(latentNames, "_cint")); latentNamesTmp
          T0MEANStmp <- matrix(paste0("Mean", latentNamesTmp), n.latent*2, 1); T0MEANStmp
          currentModel <- (
            ctsem::ctModel(n.latent=n.latent*2,
                           n.manifest=n.var,
                           manifestNames=manifestNamesTmp,
                           latentNames = latentNamesTmp,
                           DIFFUSION=DIFFUSIONtmp,
                           DRIFT=DRIFTtmp,
                           LAMBDA=LAMBDAtmp,
                           CINT=matrix(0, nrow=n.latent*2, ncol=1),
                           T0MEANS = T0MEANStmp,
                           MANIFESTMEANS = matrix(manifestMeansParams, nrow=n.latent, ncol=1),
                           MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var),
                           T0VAR = T0VARtmp,
                           type = 'stanct')
          )
          currentModel$pars$indvarying <- FALSE

          if (!(optimize)) {
            customPar <- FALSE
            if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
            tmp1a <- paste0(" Bayesian sampling was selected, which does require appropriate scaling of time. ")
            tmp2 <- nchar(tmp1a); tmp2
            tmp3 <- (81 - tmp2)/2; tmp3
            tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
            tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
            tmp6a <- paste0(tmp4, tmp1a, tmp5); tmp6a

            tmp1b <- paste0(" See the end of the summary output ")
            tmp2 <- nchar(tmp1b); tmp2
            tmp3 <- (81 - tmp2)/2; tmp3
            tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
            tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
            tmp6b <- paste0(tmp4, tmp1b, tmp5); tmp6b

            Msg <- paste0("################################################################################# \n", tmp6a, "\n", tmp6b, "\n#################################################################################")
            message(Msg)
            #Msg <- "Bayesian sampling was selected, which does require appropriate scaling of time. See the end of the summary output\n"
            #message(Msg)
          }

        }

        if (fit == FALSE) {
          print(paste0("#################################################################################"))
          print(paste0("#############  No model is fitted, only data and code are generated. ############"))
          print(paste0("#################################################################################"))
        }

        # FIT STANCT MODEL
        if (fit == TRUE) {
          if (doPar < 2) {
            # CHD changed 7 Oct 2022
            if (any(is.na(studyList[[i]]$startValues))) inits <- NULL else inits <- studyList[[i]]$startValues
            #results <- suppressMessages(ctsem::ctStanFit(
            results <- (ctsem::ctStanFit(
              datalong = emprawLong[[i]],
              ctstanmodel = currentModel,
              fit=fit,
              sameInitialTimes=sameInitialTimes,
              #inits=studyList[[i]]$startValues,
              inits=inits,
              savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
              stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
              iter=CoTiMAStanctArgs$iter,
              intoverstates=CoTiMAStanctArgs$intoverstates,
              binomial=CoTiMAStanctArgs$binomial,
              #fit=CoTiMAStanctArgs$fit,
              intoverpop=CoTiMAStanctArgs$intoverpop,
              stationary=CoTiMAStanctArgs$stationary,
              plot=CoTiMAStanctArgs$plot,
              #derrind=CoTiMAStanctArgs$derrind, # deprecated, CHD deleted Aug 2023
              optimize=CoTiMAStanctArgs$optimize,
              optimcontrol=CoTiMAStanctArgs$optimcontrol,
              nlcontrol=CoTiMAStanctArgs$nlcontrol,
              #nopriors=CoTiMAStanctArgs$nopriors,
              priors=CoTiMAStanctArgs$priors,
              chains=CoTiMAStanctArgs$chains,
              forcerecompile=CoTiMAStanctArgs$forcerecompile,
              savescores=CoTiMAStanctArgs$savescores,
              gendata=CoTiMAStanctArgs$gendata,
              control=CoTiMAStanctArgs$control,
              verbose=verbose,
              warmup=CoTiMAStanctArgs$warmup,
              cores=coresToUse) )
          }
          if (doPar > 1) {
            # parallel re-fitting of problem study

            tmp1 <- paste0(" Parallel fit attmepts requested. Screen remains silent for a while. ")
            tmp2 <- nchar(tmp1); tmp2
            tmp3 <- (81 - tmp2)/2; tmp3
            tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
            tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
            tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6

            Msg <- paste0("################################################################################# \n", tmp6, "\n#################################################################################")
            message(Msg)

            #Msg <- "Parallel fit attmepts requested. Screen remains silent for a while.\n"
            #message(Msg)

            allfits <- foreach::foreach(p=1:doPar) %dopar% {
              # CHD changed 7 Oct 2022
              if (any(is.na(studyList[[i]]$startValues))) inits <- NULL else inits <- studyList[[i]]$startValues

              fits <- suppressMessages(ctsem::ctStanFit(
                datalong = emprawLong[[i]],
                ctstanmodel = currentModel,
                sameInitialTimes=sameInitialTimes,
                #inits=studyList[[i]]$startValues,
                inits=inits,
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
                #nopriors=CoTiMAStanctArgs$nopriors,
                priors=CoTiMAStanctArgs$priors,
                chains=CoTiMAStanctArgs$chains,
                forcerecompile=CoTiMAStanctArgs$forcerecompile,
                savescores=CoTiMAStanctArgs$savescores,
                gendata=CoTiMAStanctArgs$gendata,
                control=CoTiMAStanctArgs$control,
                verbose=verbose,
                warmup=CoTiMAStanctArgs$warmup,
                cores=1) )
              return(fits)
            }
            all_loglik <- unlist(lapply(allfits, function(x) x$stanfit$optimfit$value)); all_loglik
            # CHD added 27 SEP 2022 to prevent neg -2ll fits
            #if(posLL == FALSE) {
            #  if (all(all_loglik > 0)) {
            #    ErrorMsg <- "\n All loglik values > 0, but you provided the argument posLL=FALSE, so no fit confirmed your expectations and I had to stop!"
            #    stop(ErrorMsg)
            #  }
            #  all_loglik <- all_loglik[-(which(all_loglik > 0))]
            #}
            # CHD added 27 SEP 2022: changed min to max
            bestFit <- which(abs(all_loglik) == max(abs(all_loglik)))[1]; bestFit
            #
            results <- allfits[[bestFit]]
          }

          gc()
          #
          studyFit[[i]] <- results
          if (is.null(studyFit[[i]]$standata$priors)) studyFit[[i]]$standata$priors <- 0 # CHD added Sep 2023
          studyFit[[i]]$resultsSummary <- summary(studyFit[[i]])

          df <- "deprecated"
          studyFit[[i]]$resultsSummary$'df (CoTiMA)' <- df
          #} # END if (!(studyList[[i]]$originalStudyNo %in% ...

          # SAVE
          if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
            x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
            x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/"); x2
            # CHD changed 17 Feb 2024
            tmp <- studyFit[[i]]
            tmp$empraw <- empraw[[i]]
            tmp$primaryStudyList=primaryStudies
            tmp$studyList=studyList
            tmp$argumentListOfCtmaInit=list(primaryStudies=deparse(substitute(primaryStudies)),
                                            activateRPB=activateRPB,
                                            activeDirectory=activeDirectory,
                                            chains=chains,
                                            checkSingleStudyResults=checkSingleStudyResults,
                                            cint=cint,
                                            coresToUse=coresToUse,
                                            customPar=customPar,
                                            diff=diff,
                                            digits=digits,
                                            doPar=doPar,
                                            drift=drift,
                                            experimental=experimental,
                                            finishsamples=finishsamples,
                                            indVarying=indVarying,
                                            indVaryingT0=indVaryingT0,
                                            iter=iter,
                                            lambda=lambda,
                                            loadSingleStudyModelFit=loadSingleStudyModelFit,
                                            manifestMeans=manifestMeans,
                                            manifestVars=manifestVars,
                                            n.latent=n.latent,
                                            n.manifest=n.manifest,
                                            optimize=optimize,
                                            primaryStudies=primaryStudies,
                                            priors=priors,
                                            sameInitialTimes=sameInitialTimes,
                                            saveRawData=saveRawData,
                                            saveSingleStudyModelFit=saveSingleStudyModelFit,
                                            scaleTI=scaleTI,
                                            scaleTime=scaleTime,
                                            silentOverwrite=silentOverwrite,
                                            T0means=T0means,
                                            T0var=T0var,
                                            useSV=useSV,
                                            verbose=verbose,
                                            randomIntercepts=randomInterceptsSettings,
                                            CoTiMAStanctArgs=CoTiMAStanctArgs)
            #
            #ctmaSaveFile(activateRPB, activeDirectory, studyFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
            ctmaSaveFile(activateRPB, activeDirectory, tmp, x1, x2, silentOverwrite=silentOverwrite)
          }
        }
      } # end if (tmpLogic == 1) moved up to ~ 1257

      if (fit == TRUE) {
        resultsSummary <- studyFit[[i]]$resultsSummary; resultsSummary
        if (!(exists("currenModel"))) currentModel <- studyFit[[i]]$ctstanmodelbase

        studyFit_Minus2LogLikelihood[[i]] <- -2 * resultsSummary$loglik
        studyFit_estimatedParameters[[i]] <- resultsSummary$npars

        # Subsequent if ... else ... became necessary with ctsem 3.4.1 where the namens of the rows (rownames) were moved into a column named "matrix"
        tmp <- grep("toV", rownames(resultsSummary$popmeans)); tmp
        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"]) == 0)) {
          model_Drift_Coef[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"], n.latent)); model_Drift_Coef[[i]]
        } else { # new version of ctsem
          model_Drift_Coef[[i]] <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "Mean"], n.latent)); model_Drift_Coef[[i]]
          tmp8 <- which(!(is.na(currentModel$pars[currentModel$pars$matrix == "DRIFT",]$param))); tmp8
          model_Drift_Coef[[i]] <- model_Drift_Coef[[i]][tmp8]; model_Drift_Coef[[i]]
        }
        names(model_Drift_Coef[[i]]) <- driftFullNames; model_Drift_Coef[[i]]

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"]) == 0)) {
          model_Drift_SE[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"], n.latent, byrow=FALSE)); model_Drift_SE[[i]]
        } else { # new version of ctsem
          model_Drift_SE[[i]] <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "sd"], n.latent, byrow=FALSE)); model_Drift_SE[[i]]
          tmp8 <- which(!(is.na(currentModel$pars[currentModel$pars$matrix == "DRIFT",]$param))); tmp8
          model_Drift_SE[[i]] <- model_Drift_SE[[i]][tmp8]
        }
        names(model_Drift_SE[[i]]) <- c(driftFullNames); model_Drift_SE[[i]]

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"]) == 0 )) {
          tmp1 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=FALSE)); tmp1
          tmp2 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=FALSE)); tmp2
        } else {
          tmp1 <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "2.5%"], n.latent, byrow=FALSE)); tmp1
          tmp2 <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "97.5%"], n.latent, byrow=FALSE)); tmp2
          tmp8 <- which(!(is.na(currentModel$pars[currentModel$pars$matrix == "DRIFT",]$param))); tmp8
          tmp1 <- tmp1[tmp8]
          tmp2 <- tmp2[tmp8]
        }
        # CHD 19. Nov. 2023 Check if LL == UL or LL == 0  or UL == 0, indicating estimation problems (estProb)
        #estProb[[i]] <- list()
        tmp3a <- which(tmp1 - tmp2 == 0); tmp3a
        tmp3b <- which(tmp1 == 0); tmp3b
        tmp3c <- which(tmp2 == 0); tmp3c
        tmp4 <- unique(c(tmp3a, tmp3b, tmp3c)); tmp4
        tmp5 <- grep("toV", rownames(resultsSummary$popmeans)); tmp5
        tmp6 <- rownames(resultsSummary$popmeans)[tmp]; tmp6
        #if (any(tmp4 != 0)) estProb$DRIFT[[i]] <- paste0("Possible problems for Study ", i, " in estimating: ", paste0(tmp6[tmp4], collapse=" "))
        if (any(tmp4 != 0)) estProb[[length(estProb)+1]] <- paste0("Possible problems for Study ", i, " in estimating: ", paste0(tmp6[tmp4], collapse=" "))

        model_Drift_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Drift_CI[[i]]
        tmp3 <- c(rbind(paste0(driftFullNames, "LL"),
                        paste0(driftFullNames, "UL"))); tmp3
        names(model_Drift_CI[[i]]) <- tmp3; model_Drift_CI[[i]]

        tmp <- grep("diff", rownames(resultsSummary$popmeans)); tmp
        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"]) == 0)) {
          model_Diffusion_Coef[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"])
          names(model_Diffusion_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]
        } else {
          model_Diffusion_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "Mean"])
          tmp8 <- which(!(is.na(currentModel$pars[currentModel$pars$matrix == "DRIFT",]$param))); tmp8
          model_Diffusion_Coef[[i]] <- model_Diffusion_Coef[[i]][tmp8]
          names(model_Diffusion_Coef[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
        }

        if (!(is.null(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Sd"]))) {
          model_Diffusion_SE[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Sd"]) #; model_Diffusion_SE[[i]]
          names(model_Diffusion_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]
        } else {
          model_Diffusion_SE[[i]] <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "sd"] #; model_Diffusion_SE[[i]]
          tmp8 <- which(!(is.na(currentModel$pars[currentModel$pars$matrix == "DRIFT",]$param))); tmp8
          model_Diffusion_SE[[i]] <- model_Diffusion_SE[[i]][tmp8]
          names(model_Diffusion_SE[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
        }

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "2.5%"])) == 0) {
          tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "97.5%"]; tmp2
          model_Diffusion_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Diffusion_CI[[i]]
          tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                          paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
          names(model_Diffusion_CI[[i]]) <- tmp3; model_Diffusion_CI[[i]]
        } else {
          tmp1 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "97.5%"]; tmp2
          tmp8 <- which(!(is.na(currentModel$pars[currentModel$pars$matrix == "DRIFT",]$param))); tmp8
          tmp1 <- tmp1[tmp8]
          tmp2 <- tmp2[tmp8]
          model_Diffusion_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Diffusion_CI[[i]]
          tmp3 <- c(rbind(paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "LL"),
                          paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "UL"))); tmp3
          names(model_Diffusion_CI[[i]]) <- tmp3; model_Diffusion_CI[[i]]
        }
        # CHD 19. Nov. 2023 Check if LL == UL or LL == 0  or UL == 0, indicating estimation problems (estProb)
        #estProb[[i]] <- list()
        tmp3a <- which(tmp1 - tmp2 == 0); tmp3a
        tmp3b <- which(tmp1 == 0); tmp3b
        tmp3c <- which(tmp2 == 0); tmp3c
        tmp4 <- unique(c(tmp3a, tmp3b, tmp3c)); tmp4
        tmp5 <- grep("diff", rownames(resultsSummary$popmeans)); tmp5
        tmp6 <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp5])); tmp6
        #if (any(tmp4 != 0)) estProb$DIFF[[i]] <- paste0("Possibly problems for Study ", i, " in estimating: ", paste0(tmp6[tmp4], collapse=" "))
        #if (any(tmp4 != 0)) estProb[[length(estProb)+1]] <- paste0("Possibly problems for Study ", i, " in estimating: ", paste0(tmp6[tmp4], collapse=" "))
        if (any(tmp4 != 0)) estProb[[length(estProb)+0]] <- paste0(estProb[[length(estProb)+0]], " ", paste0(tmp6[tmp4], collapse=" "))

        #if (randomIntercepts == TRUE) target <- "T0cov_" else target <- "0var"
        if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )  target <- "T0cov_" else target <- "0var"
        #if (randomIntercepts == TRUE) target <- "ov_" else target <- "0var"
        #tmp <- grep("0var", rownames(resultsSummary$popmeans)); tmp
        tmp <- grep(target, rownames(resultsSummary$popmeans)); tmp
        # added 9. Aug. 2022. Next one became neccessary because ctsem labeling changed from "var" to "cov"
        if (length(tmp) == 0) {
          tmp <- grep("0cov", resultsSummary$parmatrices$matrix)
          tmp2 <- (resultsSummary$parmatrices[tmp, c("matrix", "row", "col")]); tmp2
          T0covNames <- c()
          for (m in 1:nrow(tmp2)) T0covNames[m] <- paste0(tmp2[m, 1], tmp2[m, 2], tmp2[m, 3])
        } else {
          T0covNames <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
        }
        #T0covNames

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Mean"])) == 0 ) {
          model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Mean"])
          names(model_T0var_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_T0var_Coef[[i]]
        }  else {
          model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "Mean"])
          # added 9. Aug. 2022
          if (length(model_T0var_Coef[[i]]) != n.latent^2) {
            names(model_T0var_Coef[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
          } else {
            names(model_T0var_Coef[[i]]) <- T0covNames
          }
          #if (randomIntercepts == TRUE) {
          if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )  {
            model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "Mean"])
            model_T0var_Coef[[i]] <- c(matrix(model_T0var_Coef[[i]], n.latent^2, n.latent^2)[1:n.latent, 1:n.latent])
            names(model_T0var_Coef[[i]]) <- T0covNames
          }
        }


        if (!(is.null(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Sd"]))) {
          model_T0var_SE[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Sd"]); model_T0var_SE[[i]]
          names(model_T0var_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_T0var_SE[[i]]
        } else {
          model_T0var_SE[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "sd"])
          if (length(model_T0var_SE[[i]]) != n.latent^2) {
            names(model_T0var_SE[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
          } else {
            names(model_T0var_SE[[i]]) <- T0covNames
          }
          #if (randomIntercepts == TRUE) {
          if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") ) {
            model_T0var_SE[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "sd"])
            model_T0var_SE[[i]] <- c(matrix(model_T0var_SE[[i]], n.latent^2, n.latent^2)[1:n.latent, 1:n.latent])
            names(model_T0var_SE[[i]]) <- T0covNames
          }
        }
        #model_T0var_SE

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]) == 0)) {
          tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "97.5%"]; tmp2
          model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]
          tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                          paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
          names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
        } else {
          tmp1 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "97.5%"]; tmp2
          model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]
          if (length(tmp1) != n.latent^2) {
            tmp3 <- c(rbind(paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "LL"),
                            paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "UL"))); tmp3
            names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
          } else {
            tmp3 <- c(rbind(paste0(T0covNames, "LL"),
                            paste0(T0covNames, "UL"))); tmp3
            names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
          }
          #if (randomIntercepts == TRUE) {
          if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )  {
            tmp1 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "2.5%"]; tmp1
            tmp1 <- c(matrix(tmp1, n.latent^2, n.latent^2)[1:n.latent, 1:n.latent])
            tmp2 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "97.5%"]; tmp2
            tmp2 <- c(matrix(tmp2, n.latent^2, n.latent^2)[1:n.latent, 1:n.latent])
            model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]
            tmp3 <- c(rbind(paste0(T0covNames, "LL"),
                            paste0(T0covNames, "UL"))); tmp3
            names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
          }
        }
        #model_T0var_CI

        # CHD 19. Nov. 2023 Check if LL == UL or LL == 0  or UL == 0, indicating estimation problems (estProb)
        #estProb[[i]] <- list()
        tmp3a <- which(tmp1 - tmp2 == 0); tmp3a
        tmp3b <- which(tmp1 == 0); tmp3b
        tmp3c <- which(tmp2 == 0); tmp3c
        tmp4 <- unique(c(tmp3a, tmp3b, tmp3c)); tmp4
        tmp5 <- grep("0var", rownames(resultsSummary$popmeans)); tmp5
        tmp6 <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp5])); tmp6
        #if (any(tmp4 != 0)) estProb$DIFF[[i]] <- paste0("Possibly problems for Study ", i, " in estimating: ", paste0(tmp6[tmp4], collapse=" "))
        #if (any(tmp4 != 0)) estProb[[length(estProb)+1]] <- paste0("Possibly problems for Study ", i, " in estimating: ", paste0(tmp6[tmp4], collapse=" "))
        if (any(tmp4 != 0)) estProb[[length(estProb)+0]] <- paste0(estProb[[length(estProb)+0]], " ", paste0(tmp6[tmp4], collapse=" "))


        # changed 17. Aug. 2022
        #if ( ((indVarying == TRUE) | (indVarying == "CINT") | (indVarying == 'cint')) & (randomIntercepts == FALSE) )  {
        if ( ( (indVarying == TRUE) | (indVarying == 'CINT') ) & ( (randomIntercepts != TRUE) | (randomIntercepts != "MANIFETS") ) ) {
          # CHD 12.6.2023
          e <- ctsem::ctExtract(studyFit[[i]])
          model_popsd_tmp <- resultsSummary$popsd
          #if (indVaryingT0 == TRUE) {
          #if ( (indVaryingT0 == TRUE) & (T0meansParams[1] != 0) ) {
          if (dim(model_popsd_tmp)[1] != n.latent) {
            model_popsd[[i]] <- resultsSummary$popsd
            model_popcov_m[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits)
            model_popcov_sd[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits)
            model_popcov_T[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean)/ctsem::ctCollapse(e$popcov, 1, stats::sd), digits)
            model_popcov_025[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025))
            model_popcov_50[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .50))
            model_popcov_975[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975))
            # convert to correlations and do the same (array to list then list to array)
            e$popcor <- lapply(seq(dim(e$popcov)[1]), function(x) e$popcov[x , ,])
            e$popcor <- lapply(e$popcor, stats::cov2cor)
            e$popcor <- array(unlist(e$popcor), dim=c(n.latent*2, n.latent*2, length(e$popcor)))
            model_popcor_m[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
            model_popcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
            model_popcor_T[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
            model_popcor_025[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .025))
            model_popcor_50[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .50))
            model_popcor_975[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .975))
            #model_popcor <- stats::cov2cor(model_popcov_m)
          } #else {
          #  model_popcov_m[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits); model_popcov_m
          #  model_popcov_sd[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits)
          #  model_popcov_T[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean)/ctsem::ctCollapse(e$popcov, 1, stats::sd), digits)
          #  model_popcov_025[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025))
          #  model_popcov_50[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .50))
          #  model_popcov_975[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975))
          #  # convert to correlations and do the same (array to list then list to array)
          #  e$popcor <- lapply(seq(dim(e$popcov)[1]), function(x) e$popcov[x , ,])
          #  e$popcor <- lapply(e$popcor, stats::cov2cor)
          #  e$popcor <- array(unlist(e$popcor), dim=c(n.latent   , n.latent  , length(e$popcor)))
          #  model_popcor_m[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
          #  model_popcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
          #  model_popcor_T[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
          #  model_popcor_025[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .025))
          #  model_popcor_50[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .50))
          #  model_popcor_975[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .975))
          #}
        }

        # CHD 22.1.2024
        #if ( ( !((indVarying == TRUE) | (indVarying == 'cint') | (indVarying == 'CINT')) & (randomIntercepts == FALSE)) ) {
        if ( (indVarying != TRUE) & (indVarying != 'CINT') & (randomIntercepts != TRUE) & (randomIntercepts != "MANIFEST") ) {
          model_popsd[[i]] <- "no random intercepts estimated"
          model_popcov_m[[i]] <- model_popcov_sd[[i]] <- model_popcov_T[[i]] <- model_popcov_025[[i]] <- model_popcov_50[[i]] <- model_popcov_975[[i]] <- "no random intercepts estimated"
          model_popcor_m[[i]] <- model_popcor_sd[[i]] <- model_popcor_T[[i]] <- model_popcor_025[[i]] <- model_popcor_50[[i]] <- model_popcor_975[[i]] <- "no random intercepts estimated"
        }

        #if (randomIntercepts == TRUE) {
        if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST")  ) {
          e <- ctsem::ctExtract(studyFit[[i]])
          ctsem::ctCollapse(e$pop_T0cov, 1, mean)
          model_popcov_m[[i]] <- round(ctsem::ctCollapse(e$pop_T0cov, 1, mean), digits = digits)
          model_popcov_sd[[i]] <- round(ctsem::ctCollapse(e$pop_T0cov, 1, stats::sd), digits = digits)
          model_popcov_T[[i]] <- round(ctsem::ctCollapse(e$pop_T0cov, 1, mean)/ctsem::ctCollapse(e$pop_T0cov, 1, stats::sd), digits)
          model_popcov_025[[i]] <- ctsem::ctCollapse(e$pop_T0cov, 1, function(x) stats::quantile(x, .025))
          model_popcov_50[[i]] <- ctsem::ctCollapse(e$pop_T0cov, 1, function(x) stats::quantile(x, .50))
          model_popcov_975[[i]] <- ctsem::ctCollapse(e$pop_T0cov, 1, function(x) stats::quantile(x, .975))
          # convert to correlations and do the same (array to list then list to array)
          e$popcor <- lapply(seq(dim(e$pop_T0cov)[1]), function(x) e$pop_T0cov[x , ,])
          e$popcor <- lapply(e$popcor, stats::cov2cor)
          e$popcor <- array(unlist(e$popcor), dim=c(n.latent*2, n.latent*2, length(e$popcor)))
          model_popcor_m[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
          model_popcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
          model_popcor_T[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
          model_popcor_025[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .025))
          model_popcor_50[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .50))
          model_popcor_975[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .975))
          #model_popcor <- stats::cov2cor(model_popcov_m)
        } # end     if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST")  )

      } # end if (fit == TRUE)
      #} # end if (tmpLogic == 1) moved up to ~ 1257

    } # END     for (i in 1:n.studies)

  }

  #######################################################################################################################
  ################################ Combine summary information and fit statistics #######################################
  #######################################################################################################################

  if (fit == TRUE) {
    allStudies_Minus2LogLikelihood <- sum(unlist(studyFit_Minus2LogLikelihood)); allStudies_Minus2LogLikelihood
    allStudies_estimatedParameters <- sum(unlist(studyFit_estimatedParameters)); allStudies_estimatedParameters
    allStudies_df <- "deprecated"
    #
    tmp1 <- driftFullNames; tmp1
    tmp2 <- rep("SE", length(tmp1)); tmp2
    targetNames1 <- c(rbind(tmp1, tmp2)); targetNames1
    #
    targetNames2 <- c()
    for (i in 1:n.latent) {
      for (j in 1:n.latent) {
        if (i != j) {
          targetNames2 <- c(targetNames2, paste0("_eta", i, "_eta", j))
        } else {
          targetNames2 <- c(targetNames2, paste0("_eta", i))
        }
      }
    }
    targetNames2  <- paste0("diff", targetNames2); targetNames2
    targetNames2  <- c(rbind(targetNames2, rep("SE", 4))); targetNames2
    #
    targetNames3 <- targetNames2
    targetNames3 <- gsub("eta", "V", targetNames3)
    targetNames3 <- gsub("diff", "T0var", targetNames3); targetNames3

    allStudiesDRIFT_effects <- matrix(t(cbind(unlist(model_Drift_Coef), unlist(model_Drift_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    colnames(allStudiesDRIFT_effects) <- targetNames1; allStudiesDRIFT_effects

    allStudiesDIFF_effects <- matrix(t(cbind(unlist(model_Diffusion_Coef), unlist(model_Diffusion_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    colnames(allStudiesDIFF_effects) <- targetNames2; allStudiesDIFF_effects

    allStudiesT0VAR_effects <- matrix(t(cbind(unlist(model_T0var_Coef), unlist(model_T0var_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    colnames(allStudiesT0VAR_effects) <- targetNames3; allStudiesT0VAR_effects

    #
    if (!(is.null(scaleTime))) {
      allStudiesDRIFT_effects_original_time_scale <- matrix(t(cbind(unlist(model_Drift_Coef) * scaleTime,
                                                                    unlist(model_Drift_SE) * scaleTime) ), n.studies, 2*n.latent^2, byrow=T)
      colnames(allStudiesDRIFT_effects_original_time_scale) <- targetNames1
      allStudiesDRIFT_effects_original_time_scale <- round(allStudiesDRIFT_effects_original_time_scale, digits)
      #
      allStudiesDIFF_effects_original_time_scale <- matrix(t(cbind(unlist(model_Diffusion_Coef) * scaleTime,
                                                                   unlist(model_Diffusion_SE) * scaleTime) ), n.studies, 2*n.latent^2, byrow=T)
      colnames(allStudiesDIFF_effects_original_time_scale) <- targetNames2
      allStudiesDIFF_effects_original_time_scale <- round(allStudiesDIFF_effects_original_time_scale, digits)
    } else {
      allStudiesDRIFT_effects_original_time_scale <- NULL
      allStudiesDIFFUSION_effects_original_time_scale <- NULL
    }

    source <- lapply(primaryStudies$source, function(extract) paste(extract, collapse=", ")); source
    if (length(source) > n.studies) source[n.studies +1] <- NULL # new 6. July< 2022
    for (l in 1:length(source)) if ( source[[l]] == "NA") source[[l]] <- "Reference not provided"
    #
    allStudiesDRIFT_effects_ext <- cbind(unlist(source), allStudiesDRIFT_effects)
    tmp <- allStudiesDRIFT_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesDRIFT_effects_ext <- tmp; allStudiesDRIFT_effects_ext
    #
    allStudiesDIFF_effects_ext <- cbind(unlist(source), allStudiesDIFF_effects)
    tmp <- allStudiesDIFF_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesDIFF_effects_ext <- tmp; allStudiesDIFF_effects_ext
    #
    allStudiesT0VAR_effects_ext <- cbind(unlist(source), allStudiesT0VAR_effects)
    tmp <- allStudiesT0VAR_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesT0VAR_effects_ext <- tmp; allStudiesT0VAR_effects_ext

    if (!(is.null(allStudiesDRIFT_effects_original_time_scale))) {
      #
      allStudiesDRIFT_effects_original_time_scale_ext <- cbind(unlist(source), allStudiesDRIFT_effects_original_time_scale)
      tmp <- allStudiesDRIFT_effects_original_time_scale_ext
      tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
      allStudiesDRIFT_effects_original_time_scale_ext <- tmp; allStudiesDRIFT_effects_original_time_scale_ext
      #allStudiesDRIFT_effects_ext
      #
      allStudiesDIFF_effects_original_time_scale_ext <- cbind(unlist(source), allStudiesDIFF_effects_original_time_scale)
      tmp <- allStudiesDIFF_effects_original_time_scale_ext
      tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
      allStudiesDIFF_effects_original_time_scale_ext <- tmp; allStudiesDIFF_effects_original_time_scale_ext
      #allStudiesDIFF_effects_ext
    } else {
      allStudiesDRIFT_effects_original_time_scale_ext <- allStudiesDRIFT_effects_ext # new 8.7.2022
      allStudiesDIFF_effects_original_time_scale_ext <- allStudiesDIFF_effects_ext # new 8.7.2022
    }

    # confidence intervals
    allStudiesDriftCI <- matrix(unlist(model_Drift_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesDriftCI) <- names(model_Drift_CI[[1]]); allStudiesDriftCI
    allStudiesDiffusionCI <- matrix(unlist(model_Diffusion_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesDiffusionCI) <- names(model_Diffusion_CI[[1]]); allStudiesDiffusionCI
    allStudiesT0varCI <- matrix(unlist(model_T0var_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesT0varCI) <- names(model_T0var_CI[[1]]); allStudiesT0varCI

    if (!(is.null(scaleTime))) {
      #
      allStudiesDriftCI_original_time_scale <- matrix(round(unlist(model_Drift_CI) * scaleTime, digits),
                                                      nrow=n.studies, byrow=TRUE)
      colnames(allStudiesDriftCI_original_time_scale) <- names(model_Drift_CI[[1]]); allStudiesDriftCI_original_time_scale
      #
      allStudiesDiffCI_original_time_scale <- matrix(round(unlist(model_Diffusion_CI) * scaleTime, digits),
                                                     nrow=n.studies, byrow=TRUE)
      colnames(allStudiesDiffCI_original_time_scale) <- names(model_Diffusion_CI[[1]]); allStudiesDiffCI_original_time_scale
    } else {
      allStudiesDriftCI_original_time_scale <- allStudiesDriftCI  # new 8.7.2022
      allStudiesDiffCI_original_time_scale <- allStudiesDiffusionCI  # new 8.7.2022
    }


    allStudiesCI <- cbind(allStudiesDriftCI, allStudiesDiffusionCI, allStudiesT0varCI); allStudiesCI
    allStudiesCI_ext <- cbind(allStudiesDRIFT_effects_ext[,1], allStudiesCI); allStudiesCI_ext

    allStudiesCI_original_time_scale <- allStudiesCI_ext
    tmp1 <- grep("T0", colnames(allStudiesCI_original_time_scale)); tmp1
    tmp2 <- allStudiesCI_original_time_scale[, -c(1, tmp1)]; tmp2
    if (!(is.null(scaleTime))) scaleTime2 <- scaleTime else scaleTime2 <- 1
    if (is.null(ncol(tmp2))) tmp2b <- length(tmp2) else tmp2b <- ncol(tmp2)
    tmp3b <- matrix(round(as.numeric(tmp2) * scaleTime2, digits), ncol=tmp2b); tmp3b

    if (!(is.matrix(allStudiesCI_original_time_scale[, 1]))) {
      tmp3a <- matrix(allStudiesCI_original_time_scale[, 1], nrow=length(allStudiesCI_original_time_scale[, 1]))
    } else {
      tmp3a <- allStudiesCI_original_time_scale[, 1]
    }
    if (!(is.matrix(allStudiesCI_original_time_scale[, tmp1]))) {
      tmp3c <- matrix(allStudiesCI_original_time_scale[, tmp1], nrow=length(allStudiesCI_original_time_scale[, tmp1]))
    } else {
      tmp3c <- allStudiesCI_original_time_scale[, tmp1]
    }

    if (dim(tmp3c)[2] == 1) tmp3c <- c(tmp3c) # new 6.7.2022
    if (is.null(dim(tmp3c))) tmp3 <- cbind(tmp3a, tmp3b, t(tmp3c)) else tmp3 <- cbind(tmp3a, tmp3b,   tmp3c)
    colnames(tmp3) <- colnames(allStudiesCI_original_time_scale); tmp3
    allStudiesCI_original_time_scale <- tmp3

    # Label summary table
    rownames(allStudiesDRIFT_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDRIFT_effects_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDIFF_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDIFF_effects_ext) <- paste0("Study No ", primaryStudies$studyNumbers)

    if (!(is.null(scaleTime))) {
      #
      rownames(allStudiesDRIFT_effects_original_time_scale) <- paste0("Study No ", primaryStudies$studyNumbers)
      rownames(allStudiesDRIFT_effects_original_time_scale_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
      #
      rownames(allStudiesDIFF_effects_original_time_scale) <- paste0("Study No ", primaryStudies$studyNumbers)
      rownames(allStudiesDIFF_effects_original_time_scale_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
    }

    # dt effects
    allStudiesDRIFT_effects_ext_dt <- allStudiesDRIFT_effects_ext
    tmp1 <- grep("toV", colnames(allStudiesDRIFT_effects_ext_dt))
    tmp2 <- (allStudiesDRIFT_effects_ext[, tmp1])
    if (!(is.null(dim(tmp2)))) {
      for (l in 1:dim(tmp2)[1]) {
        tmp3 <- matrix(as.numeric(tmp2[l, ]), n.latent, n.latent, byrow=TRUE)
        # changed 16 Sep 2022
        if (is.null(scaleTime)) tmp4 <- c(t(OpenMx::expm(tmp3))) else tmp4 <- c(t(OpenMx::expm(tmp3 * scaleTime)))
        allStudiesDRIFT_effects_ext_dt[l, tmp1] <- round(tmp4, digits)
      }
    } else {
      tmp3 <- matrix(as.numeric(tmp2), n.latent, n.latent, byrow=TRUE)
      # changed 16 Sep 2022
      if (is.null(scaleTime)) tmp4 <- c(t(OpenMx::expm(tmp3))) else tmp4 <- c(t(OpenMx::expm(tmp3 * scaleTime)))
      allStudiesDRIFT_effects_ext_dt[tmp1]
      allStudiesDRIFT_effects_ext_dt[tmp1] <- round(tmp4, digits)
    }
    tmp1 <- grep("SE", colnames(allStudiesDRIFT_effects_ext_dt))
    allStudiesDRIFT_effects_ext_dt <- allStudiesDRIFT_effects_ext_dt[, -tmp1]
    if (!(is.null(dim(allStudiesDRIFT_effects_ext_dt)))) {
      colnames(allStudiesDRIFT_effects_ext_dt) <- paste0(colnames(allStudiesDRIFT_effects_ext_dt), " discrete time")
    } else {
      names(allStudiesDRIFT_effects_ext_dt)  <- paste0(names(allStudiesDRIFT_effects_ext_dt), " discrete time")
    }

    # check single study results
    if (checkSingleStudyResults == TRUE) {
      print(allStudiesDRIFT_effects_ext)
      print(allStudiesDRIFT_effects_ext_dt)
      print(allStudies_Minus2LogLikelihood)
      #
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::blue(" Press 'q' to quit and prevent saving or any other key to continue. Press ENTER afterwards."))
      char <- readline(" ")
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    DRIFTCoeff <- matrix(unlist(model_Drift_Coef), n.studies, n.latent^2, byrow=TRUE); DRIFTCoeff
    DRIFTSE <- matrix(unlist(model_Drift_SE), n.studies, n.latent^2, byrow=TRUE); DRIFTSE
    DIFFCoeff <- matrix(unlist(model_Diffusion_Coef), n.studies, n.latent^2, byrow=TRUE); DIFFCoeff
    DIFFSE <- matrix(unlist(model_Diffusion_SE), n.studies, n.latent^2, byrow=TRUE); DIFFSE
    T0VARCoeff <- matrix(unlist(model_T0var_Coef), n.studies, n.latent^2, byrow=TRUE); T0VARCoeff
    T0VARSE <- matrix(unlist(model_T0var_SE), n.studies, n.latent^2, byrow=TRUE); T0VARSE

    if (n.studies < 2) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      Msg <- "Only a single primary study was handed over to ctmaInitFit. No further (meta-) analyses can be conducted. \nI guess this stop is intended!"
      message(Msg)
    }
  } #end if (fit == TRUE)


  ##############################################################################################################
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  if (fit == TRUE) {
    message1 <- c()

    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

    if (!(all(is.na(unlist(primaryStudies$deltas))))) {
      maxDeltas <- max(unlist(primaryStudies$deltas), na.rm=TRUE)

      if (!(is.null(scaleTime))) maxDeltas <- maxDeltas * scaleTime
      largeDelta <- which(unlist(primaryStudies$deltas) >= maxDeltas); largeDelta

      if (!(is.null(scaleTime))) {
        tmp1 <- unlist(primaryStudies$deltas)
        tmp1 <- tmp1[!(is.na(tmp1))]
        largeDelta <- which( (tmp1 * scaleTime) >= maxDeltas)
      }
      tmp1 <- table(unlist(primaryStudies$deltas)[largeDelta]); tmp1
      if (!(is.null(scaleTime))) {
        #table(unlist(primaryStudies$deltas * scaleTime)[largeDelta])
        table((tmp1 * scaleTime)[largeDelta])
      }
      tmp2 <- which(tmp1 == (max(tmp1))); tmp2
      suggestedScaleTime <- as.numeric(names(tmp1[tmp2])); suggestedScaleTime
      # maxDeltas hast already been time scaled above
      if (maxDeltas > 6) {
        #tmp2 <- paste0("Maximum time interval was ", maxDeltas, "."); tmp2 # CHD 19.11.2023
        tmp2 <- paste0("Maximum time interval was ", round(maxDeltas * scaleTime2, digits), "."); tmp2
        if (is.null(scaleTime)) scaleTime3 <- scaleTime2 else scaleTime3 <- scaleTime / scaleTime2
        tmp3 <- paste0("timeScale=1/", round(suggestedScaleTime * scaleTime3, digits)); tmp3 # CHD 19.11.2023
        if (suggestedScaleTime != scaleTime2) {
          tmp4 <- paste0("It is recommended to fit the model again using the arguments ", tmp3, ". "); tmp4
        } else {
          tmp4 <- paste0(""); tmp4
        }
        message1 <- paste(tmp2, tmp4, "If the model fit (-2ll) is better (lower), continue using, e.g.,", tmp3, "in all subsequent models.", collapse="\n"); message1
      }
    } else {
      maxDeltas <- NA
    }; maxDeltas

    if (is.null(scaleTime)) scaleTime2 <- 1 else scaleTime2 <- scaleTime

    if (!(is.null(scaleTime))) {
      model_Drift_Coef_original_time_scale <- lapply(model_Drift_Coef, function(x) x * scaleTime)
      model_Diffusion_Coef_original_time_scale <- lapply(model_Diffusion_Coef, function(x) x * scaleTime)
    } else {
      model_Drift_Coef_original_time_scale <- model_Drift_Coef
      model_Diffusion_Coef_original_time_scale <- model_Diffusion_Coef
    }

    if (randomIntercepts == TRUE) {
      randomIntercepts <- list(popsd=model_popsd,
                               popcov_mean=model_popcov_m,
                               popcov_sd=model_popcov_sd)
      # move T0Means of phantom latents to cints
      # needs still to be done
      skip <- 0
      if (skip == 1) {
        tmp1 <- grep("T0MEAN", rownames(invariantDrift_Coeff)); tmp1
        tmp1 <- tmp1[(n.latent+1):(n.latent^2)]; tmp1
        tmp1 <- invariantDrift_Coeff[tmp1, ]; tmp1
        tmp2 <- grep("CINT", rownames(invariantDrift_Coeff)); tmp2
        invariantDrift_Coeff[tmp2[1:n.latent],] <- tmp1
        # delete irrelevant rows
        toDelete <- c()
        for (i in (n.latent+1):(n.latent^2)) {
          toDelete <- c(toDelete, grep(i, rownames(invariantDrift_Coeff)))
        }
        invariantDrift_Coeff <- invariantDrift_Coeff[-toDelete, ]
        #invariantDrift_Coeff
        #estimates_original_time_scaleBackup <- estimates_original_time_scale
        toDelete <- c()
        for (i in (n.latent+1):(n.latent^2)) {
          toDelete <- c(toDelete, grep(i, rownames(estimates_original_time_scale)))
        }
        estimates_original_time_scale <- estimates_original_time_scale[-toDelete, ]
        #estimates_original_time_scale
      } # end skip
    } else {
      if ( (indVarying == 'CINT') | (indVarying == TRUE) ){
        randomIntercepts <- list(popsd=model_popsd,
                                 popcov_mean=model_popcov_m, model_popcov_sd=model_popcov_sd,
                                 model_popcov_T=model_popcov_T, model_popcov_025=model_popcov_025,
                                 model_popcov_50=model_popcov_50, model_popcov_975=model_popcov_975,
                                 popcor_mean=model_popcor_m, model_popcor_sd=model_popcor_sd,
                                 model_popcor_T=model_popcor_T, model_popcor_025=model_popcor_025,
                                 model_popcor_50=model_popcor_50, model_popcor_975=model_popcor_975)
      } else {
        randomIntercepts <- "Not estimated."
      }

    }

    results <- list(activeDirectory=activeDirectory,
                    plot.type="drift", model.type="stanct",
                    coresToUse=coresToUse, n.studies=n.studies,
                    n.latent=n.latent,
                    n.manifest=n.manifest,
                    lambda=lambdaParams,
                    primaryStudyList=primaryStudies,
                    studyList=studyList, studyFitList=studyFit,
                    emprawList=empraw,
                    statisticsList=statisticsList,
                    ind.mod.List=empraw.ind.mod,
                    argumentList=list(ctmaInitFit=deparse(substitute(primaryStudies)),
                                      activateRPB=activateRPB,
                                      activeDirectory=activeDirectory,
                                      chains=chains,
                                      checkSingleStudyResults=checkSingleStudyResults,
                                      cint=cint,
                                      coresToUse=coresToUse,
                                      customPar=customPar,
                                      diff=diff,
                                      digits=digits,
                                      doPar=doPar,
                                      drift=drift,
                                      experimental=experimental,
                                      finishsamples=finishsamples,
                                      indVarying=indVarying,
                                      indVaryingT0=indVaryingT0,
                                      iter=iter,
                                      lambda=lambda,
                                      loadSingleStudyModelFit=loadSingleStudyModelFit,
                                      manifestMeans=manifestMeans,
                                      manifestVars=manifestVars,
                                      n.latent=n.latent,
                                      n.manifest=n.manifest,
                                      #nopriors=nopriors,
                                      optimize=optimize,
                                      primaryStudies=primaryStudies,
                                      priors=priors,
                                      sameInitialTimes=sameInitialTimes,
                                      saveRawData=saveRawData,
                                      saveSingleStudyModelFit=saveSingleStudyModelFit,
                                      scaleTI=scaleTI,
                                      scaleTime=scaleTime,
                                      silentOverwrite=silentOverwrite,
                                      T0means=T0means,
                                      T0var=T0var,
                                      useSV=useSV,
                                      verbose=verbose,
                                      randomIntercepts=randomInterceptsSettings,
                                      CoTiMAStanctArgs=CoTiMAStanctArgs),
                    modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef,
                                      DRIFToriginal_time_scale=model_Drift_Coef_original_time_scale,
                                      DIFFUSIONoriginal_time_scale=model_Diffusion_Coef_original_time_scale),
                    ctModel = currentModel,
                    parameterNames=list(DRIFT=names(model_Drift_Coef[[1]]), DIFFUSION=names(model_Diffusion_Coef[[1]]), T0VAR=names(model_T0var_Coef[[1]])),
                    summary=(list(model="all drift free (het. model)",
                                  estimates=allStudiesDRIFT_effects_ext, #allStudiesDRIFT_effects_ext, = estimates that would be obtained without the scaleTime argument
                                  estimates_discrete=allStudiesDRIFT_effects_ext_dt,
                                  scaleTime=scaleTime2,
                                  # changed 17. Aug. 2022
                                  #randomIntercepts=list(popsd=model_popsd, popcov=model_popcov, popcor=model_popcor),
                                  randomIntercepts=randomIntercepts,
                                  confidenceIntervals=allStudiesCI,
                                  minus2ll= round(allStudies_Minus2LogLikelihood, digits),
                                  n.parameters = round(allStudies_estimatedParameters, digits),
                                  drift_estimates_original_time_scale =allStudiesDRIFT_effects_original_time_scale_ext,
                                  drift_CI_original_time_scale=allStudiesDriftCI_original_time_scale,
                                  diff_estimates_original_time_scale=allStudiesDIFF_effects_original_time_scale_ext,
                                  diff_CI_original_time_scale=allStudiesDiffCI_original_time_scale,
                                  estimationProblems=sapply(estProb, paste, sep="/n"),
                                  note=message1)))
    # excel workbook is added later
  } # end if (fit == TRUE)

  if (fit == FALSE) {
    results <- list(summary=c("No model was fitted, only data and code were generated. See $data & $ctModel section."),
                    data = empraw,
                    ctModel = currentModel)
  }
  class(results) <- "CoTiMAFit"

  ### prepare Excel Workbook with several sheets ################################################################
  if (fit == TRUE) {
    {
      wb <- openxlsx::createWorkbook()
      sheet1 <- openxlsx::addWorksheet(wb, sheetName="model")
      sheet2 <- openxlsx::addWorksheet(wb, sheetName="modelResults")
      sheet3 <- openxlsx::addWorksheet(wb, sheetName="estimates")
      sheet4 <- openxlsx::addWorksheet(wb, sheetName="confidenceIntervals")
      sheet5 <- openxlsx::addWorksheet(wb, sheetName="randomIntercepts")
      sheet6 <- openxlsx::addWorksheet(wb, sheetName="stats")
      openxlsx::writeData(wb, sheet1, results$summary$model)

      ### modelResults
      # DRIFT
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, matrix(c(t(matrix(driftNames, nrow=n.latent))), nrow=1), colNames = FALSE)
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
      ### random Intercepts
      startCol <- 2; startCol
      startRow <- 1; startRow
      # CHD 7. Sep 2022: Quickfix: do not report all random effect matrices
      # CHD 12. Oct 2023: Quickfix: do not report all random effect matrices
      #results$summary$randomIntercepts[[1]][1] <- "random intercept results are no longer reported in excel sheets"
      #openxlsx::writeData(wb, sheet5, startCol=startCol, startRow = startRow, results$summary$randomIntercepts[[1]][1], colNames = FALSE)
      openxlsx::writeData(wb, sheet5, startCol=startCol, startRow = startRow, "random intercept results are no longer reported in excel sheets", colNames = FALSE)
      ### stats
      startCol <- 2; startCol
      startRow <- 1; startRow
      tmp <- cbind("-2ll = ", results$summary$minus2ll, "Number of Parameters = ", results$summary$n.parameters)
      openxlsx::writeData(wb, sheet6, startCol=startCol, startRow = startRow, t(tmp), colNames = FALSE)
    }
    results$excelSheets <- wb
  }

  invisible(results)

}   # end function definition
