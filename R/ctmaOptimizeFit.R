#' ctmaOptimizeFit
#'
#' @description Replaces deprecated \code{\link{ctmaOptimizeInit}}, which was limited to initial fitting
#' (i.e., applies \code{\link{ctmaInit}}) of a primary study reFits times to capitalize on chance for obtaining
#' a hard-to-find optimal fit.
#' Now, optimizing a CoTiMA model generated with \code{\link{ctmaFit}} can also be done.
#' Using \code{\link{ctmaOptimizeFit}} could be helpful if a model yields out-of-range estimates, which could happen if the fitting
#' algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
#' \code{\link{ctmaOptimizeFit}} is like gambling, hoping that at least one set of starting values (the number it tries is specified in the reFits argument)
#' enables finding the global optimal fit.
#'
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param activeDirectory activeDirectory
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link[ctsem]{ctStanFit}} Function
#' @param ctmaFitFit a object fitted with \code{\link{ctmaFit}}
#' @param ctmaInitFit the ctmaInitFit object that was used to create the ctmaFitFit object with \code{\link{ctmaFit}}
#' @param ctStanFit a fit object created with ctStanFit (default=NULL)
#' @param ctStanFitArgs list of arguments passed forward to ctStanFit except datalong, ctstanmodel, cores & verbose  (default=NULL)
#' @param ctStanData data to be used for ctStanFit  (default=NULL)
#' @param customPar logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param iter number of iterations (default = 5000)
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}
#' @param problemStudy number (position in list) where the problem study in primaryStudies is found
#' @param randomIV logical (default = FALSE). randomly varies between "indVarying='MANIFEST'" and "indVarying='CINT'"
#' @param randomPar logical (default = FALSE). Overrides arguments used for customPar and randomly sets customPar either TRUE or FALSE
#' @param randomRI logical (default = FALSE). randomly varies between "randomIntercepts='MANIFEST'" and "randomIntercepts='CINT'"
#' @param randomScaleTime lower and upper limit (default = c(1,1)) of uniform distribution from which timeScale argument for ctmaInit is uniformly shuffled (integer)
#' @param randomScaleTI logical (default = FALSE). Overrides arguments used for scaleTI and randomly sets scaleTI either TRUE or FALSE
#' @param reFits how many reFits should be done
#' @param saveModelFits save the fit of each Fit attempt (default = FALSE).
#' @param scaleTI scale TI predictors - not recommended until version 0.5.3.1. Does not change aggregated results anyways, just interpretation of effects for dummies representing primary studies.
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param shuffleStudyList (default = FALSE) randomly re-arranges studies in primaryStudyList. We encountered a few cases where this mattered, even though it should not. Only works if ctmaFit is optimized.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#'
#' @importFrom foreach %dopar%
#' @importFrom RPushbullet pbPost
#' @importFrom stats runif
#' @importFrom methods is
#' @importFrom ctsem ctStanFit
#'
#' @examples
#' \dontrun{
#' optimFit313 <- ctmaOptimizeFit(primaryStudies=CoTiMAstudyList_3,
#'                                 activeDirectory="/Users/tmp/",  # adapt!
#'                                 problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
#'                                 reFits=10,
#'                                 n.latent=2)
#' summary(optimFit313)
#' }
#'
#' @export ctmaOptimizeFit
#'
#' @return returns a list with bestFit (= the best fit achieved), all_minus2ll (= all -2ll values for all fitted models), and summary, which
#' is printed if the summary function is applied to the returned object, and which shows the summary information of the ctsem model with the
#' best fit.
#'
ctmaOptimizeFit <- function(activateRPB=FALSE,
                            activeDirectory=NULL,
                            coresToUse=c(2),
                            CoTiMAStanctArgs=NULL,
                            ctmaFitFit=NULL,
                            ctmaInitFit=NULL,
                            ctStanFit=NULL,
                            ctStanData=NULL,
                            customPar=FALSE,
                            finishsamples=NULL,
                            iter=5000,
                            primaryStudies=NULL,
                            problemStudy=NULL,
                            randomIV=FALSE,
                            randomPar=FALSE,
                            randomRI=FALSE,
                            randomScaleTI=FALSE,
                            randomScaleTime=c(1,1),
                            saveModelFits=FALSE,
                            shuffleStudyList=FALSE,
                            reFits=NULL,
                            scaleTime=NULL,
                            scaleTI=NULL,
                            verbose=1,
                            ctStanFitArgs=list(stanmodeltext = NA,
                                               iter = 1000,
                                               intoverstates = TRUE,
                                               binomial = FALSE,
                                               fit = TRUE,
                                               intoverpop = "auto",
                                               sameInitialTimes = FALSE,
                                               stationary = FALSE,
                                               plot = FALSE,
                                               derrind = NA,
                                               optimize = TRUE,
                                               optimcontrol = list(),
                                               nlcontrol = list(),
                                               nopriors = NA,
                                               priors = FALSE,
                                               chains = 2,
                                               #cores = ifelse(optimize, getOption("mc.cores", 2L), "maxneeded"),
                                               inits = NULL,
                                               compileArgs = list(),
                                               forcerecompile = FALSE,
                                               saveCompile = TRUE,
                                               savescores = FALSE,
                                               savesubjectmatrices = FALSE,
                                               saveComplexPars = FALSE,
                                               gendata = FALSE,
                                               control = list(),
                                               #verbose = 0,
                                               vb = FALSE)
)
{

  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  {
    {
      if  (length(coresToUse) > 0) {
        if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
      }

      if (coresToUse >= parallel::detectCores()) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
        coresToUse <- parallel::detectCores() - 1
        Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash. \n"
        message(Msg)
      }
    }

    if (!(is.null(ctStanFit))) {
      if (!(is(ctStanFit, "ctStanFit"))) {
        ErrorMsg <- "\nThe ctStanFit object provided was not created with ctStanFit! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (is.null(ctStanData)) {
        ErrorMsg <- "\nA ctStanFit object was provided. A data frame for the argument ctStanData is required, too! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    # Dealing with CoTiMAStanctArgs
    if (is.null(ctStanFit)) {
      CoTiMAStanctArgsTmp <- CoTiMAStanctArgs
      if( (!(is.null(ctmaFitFit))) & (is.null(CoTiMAStanctArgs)) ) {
        CoTiMAStanctArgs <- ctmaFitFit$argumentList$CoTiMAStanctArgs
      }
      if( (is.null(ctmaFitFit)) & (is.null(CoTiMAStanctArgs)) & (!(is.null(ctmaInitFit))) ) {
        CoTiMAStanctArgs <- ctmaInitFit$argumentList$CoTiMAStanctArgs
      }
      if (!(is.null(CoTiMAStanctArgsTmp))) {
        tmp1 <- which(names(CoTiMA::CoTiMAStanctArgs) %in% names(CoTiMAStanctArgsTmp)); tmp1
        tmp2 <- CoTiMA::CoTiMAStanctArgs
        tmp2[tmp1] <- CoTiMAStanctArgsTmp
        CoTiMAStanctArgs <- tmp2
      }
      if (is.null(CoTiMAStanctArgsTmp)) CoTiMAStanctArgs <- CoTiMA::CoTiMAStanctArgs
      if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
    }

  }

  ########################################################################################################################

  #'%dopar%' <- foreach::'%dopar%' deprecated

  if (!(is.null(scaleTime))) {
    randomScaleTime[1] <- randomScaleTime[2] <- scaleTime
    Msg <- paste0("You provded the argumend scaleTime. This will override the randomScaleTime argument, and both values of the randomScaleTime argument will be set to, ", scaleTime, ".\n")
    message(Msg)
  }

  if (randomScaleTime[2] < randomScaleTime[1]) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    ErrorMsg <- "\nrandomScaleTime[1] has to be <= randomScaleTime[2]! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & (!(is.null(primaryStudies))) ) {
    ErrorMsg <- "Arguments for both ctmaFitFit and primaryStudies were provided. Only one out of the two can be chosen!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & ((is.null(ctmaInitFit))) ) {
    ErrorMsg <- "Argument for ctmaFitFit was provided but not for ctmaInitFit. Need a ctmaInitFit, too!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & (!(is.null(ctmaInitFit))) ) {
    if (ctmaFitFit$argumentList$ctmaInitFit != deparse(substitute(ctmaInitFit)))  {
      ErrorMsg <- paste0("The wrong ctmaInitFit object was provided. I need ",  ctmaFitFit$argumentList$ctmaInitFit, "!")
      stop(ErrorMsg)
    }
  }

  if((!is.null(ctStanFit)) & ( (!(is.null(ctmaFitFit))) | ((!is.null(ctmaInitFit))) ) ) {
    ErrorMsg <- "A ctStanFit was provided together with a ctmaFitFit or ctmaIniFit object. Make a decision!"
    stop(ErrorMsg)
  }

  if ( is.null(ctStanFit) & is.null(ctmaFitFit) & is.null(ctmaInitFit) )  {
    ErrorMsg <- "The ctmaFitFit or ctmaIniFit (or ctStanFit) argument has to be used."
    stop(ErrorMsg)
  }


  # Moderator Checks Moved to Sectioon where ctmaInit is optimized (not relevant if ctmaFit is optimized) # CHD Auf 2023


  # INIT Fit
  if ((is.null(ctmaFitFit)) & (is.null(ctStanFit))) {
    # CHD changed 21 SEP 2022
    ErrorMsg <- "argument primaryStudies is missing"
    if (is.null(primaryStudies))  stop(ErrorMsg)
    ErrorMsg <- "argument problemStudy is missing"
    if (is.null(problemStudy)) stop(ErrorMsg)
    ErrorMsg <- "argument reFits is missing"
    if (is.null(reFits)) stop(ErrorMsg)
    ErrorMsg <- "argument activeDirectory is missing"
    if (is.null(activeDirectory)) stop(ErrorMsg)

    # CHD 14.7.2025
    if (is.null(ctmaInitFit$argumentList$indVarying)) ctmaInitFit$argumentList$indVarying <- FALSE
    if (is.null(ctmaInitFit$argumentList$randomIntercepts)) ctmaInitFit$argumentList$randomIntercepts <- FALSE
    indVarying <- ctmaInitFit$argumentList$indVarying
    randomIntercepts <- ctmaInitFit$argumentList$randomIntercepts

    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == FALSE) & (randomRI == FALSE)) {
        #if (is.null(ctmaInitFit$argumentList$randomIntercepts)) ctmaInitFit$argumentList$randomIntercepts <- FALSE
        #if (is.null(ctmaInitFit$argumentList$indVarying)) ctmaInitFit$argumentList$indVarying <- FALSE
        Msg <- paste0("Both \"randomIV\" and \"randomRI\" are set FALSE (default). I use the same settings as used for ctmaInit(), which was
      \"randomIntercepts = ", ctmaInitFit$argumentList$randomIntercepts, "\", and
      \"indVarying = ", ctmaInitFit$argumentList$indVarying, ".\n")
        indVarying <- ctmaInitFit$argumentList$indVarying
        randomIntercepts <- ctmaInitFit$argumentList$randomIntercepts
        message(Msg)
      }
    }


    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == TRUE) & (randomRI == TRUE)) {
        Msg <- "Both \"randomIV\" and \"randomRI\" are set TRUE. Both are statistically equivalent, but since the latter is more difficult to fit I set \"randomRI = FALSE\". DO \"randomIV = FALSE\" to prevent this.\n"
        message(Msg)
        randomRI <- FALSE
      }
    }

    # CHD 14.7.2025
    if (is.null(randomIV)) {
      Msg <- "\"randomIV\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"indVarying = FALSE\". Take care.\n"
      randomIV <- FALSE
      indVarying <- FALSE
      message(Msg)
    }

    # CHD 14.7.2025
    if (is.null(randomRI)) {
      Msg <- "\"randomRI\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"randomIntercepts = FALSE\". Take care.\n"
      randomRI <- FALSE
      randomIntercepts <- FALSE
      message(Msg)
    }


    # create new study list with a single problem study only
    listElements <- names(primaryStudies); listElements
    newStudyList <- as.list(listElements)
    validElements <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues", "studyNumbers", "rawData", "empMeans", "empVars",
                       "source", "ageM", "malePercent", "occupation", "country", "alphas", "targetVariables", "recodeVariables", "combineVariables",
                       "combineVariablesNames", "missingVariables", "inits", "emprawList") #, "n.studies", "summary", "excelSheets", "plot.type")
    counter <- 0
    for (i in listElements) {
      counter <- counter + 1
      if (i %in% validElements) {
        if (i %in% c("pairwiseNs", "empcovs", "rawData", "deltas", "emprawList")) {
          newStudyList[[counter]] <- primaryStudies[[counter]][problemStudy]
        } else {
          newStudyList[[counter]] <- list(unlist(primaryStudies[[counter]][problemStudy], recursive=TRUE))
        }
      } else {
        newStudyList[[counter]] <- unlist(primaryStudies[[counter]])
      }
      if (is.logical(newStudyList[[counter]])) newStudyList[[counter]] <- NA
    }
    names(newStudyList) <- names(primaryStudies)
    newStudyList$n.studies <- 1

    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    all_scaleTime <- all_customPar <- c()
    all_scaleTI <- scaleTI
    all_usedStudyList <- primaryStudies
    all_randomIV <- all_randomRI <- c()
    warns <- errs <- list()


    for (i in 1:reFits) {
      #i <- 1
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar = c(TRUE, FALSE)[tmp1]
      }
      if (!(is.null(randomScaleTime))) {
        Msg <- paste0("Argument scaleTime is set to: ", scaleTime, ".")
        message(Msg)
      }
      if (randomPar == TRUE) {
        Msg <- paste0("Argument customPar is set to: ", customPar, ".")
        message(Msg)
      }
      all_scaleTime <- c(all_scaleTime, scaleTime)
      all_customPar <- c(all_customPar, customPar)
      all_randomIV <- c(all_randomIV)

      if (randomRI == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        randomIntercepts <- c("MANIFEST", "CINT")[tmp1]
        all_randomRI <- c(all_randomRI, randomIntercepts)
      } else {
        if (randomIntercepts != FALSE) randomIntercepts <- ctmaInitFit$argumentList$randomIntercepts
        all_randomRI <- c(all_randomRI, randomIntercepts)
      }

      if (randomIV == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        indVarying <- c("MANIFEST", "CINT")[tmp1]
        all_randomIV <- c(all_randomIV, indVarying)
      } else {
        if (indVarying != FALSE) indVarying <- ctmaInitFit$argumentList$indVarying
        all_randomIV <- c(all_randomIV, indVarying)
      }
      #

      if (is.null(finishsamples)) finishsamples <- ctmaInitFit$argumentList$finishsamples

      # CHD 12.4.24
      #if (is.null(indVarying)) indVarying <- ctmaFitFit$argumentList$indVarying

      problem <- FALSE
      #print(randomIntercepts)
      #print(indVarying)
      fit <- tryCatch(ctmaInit(primaryStudies=newStudyList,
                               coresToUse = coresToUse, # changed Aug 2023
                               scaleTime = scaleTime,
                               scaleTI=scaleTI,
                               customPar=customPar,
                               finishsamples=finishsamples,
                               iter=iter,
                               activeDirectory = activeDirectory,
                               CoTiMAStanctArgs=CoTiMAStanctArgs,
                               n.latent=ctmaInitFit$argumentList$n.latent,
                               n.manifest=ctmaInitFit$argumentList$n.manifest,
                               #indVarying = ctmaInitFit$argumentList$indVarying,
                               indVarying = indVarying,
                               checkSingleStudyResults=FALSE,
                               T0means=ctmaInitFit$argumentList$T0means,
                               manifestMeans=ctmaInitFit$argumentList$manifestMeans,
                               manifestVars=ctmaInitFit$argumentList$manifestVars,
                               chains=ctmaInitFit$argumentList$chains,
                               cint=ctmaInitFit$argumentList$cint,
                               diff=ctmaInitFit$argumentList$diff,
                               digits=ctmaInitFit$argumentList$digits,
                               drift=ctmaInitFit$argumentList$drift,
                               experimental=ctmaInitFit$argumentList$experimental,
                               indVaryingT0=ctmaInitFit$argumentList$indVaryingT0,
                               lambda=ctmaInitFit$argumentList$lambda,
                               #loadSingleStudyModelFit=loadSingleStudyModelFit,
                               #nopriors=nopriors,
                               optimize=ctmaInitFit$argumentList$optimize,
                               #primaryStudies=primaryStudies,
                               priors=ctmaInitFit$argumentList$priors,
                               sameInitialTimes=ctmaInitFit$argumentList$sameInitialTimes,
                               #saveRawData=saveRawData,
                               #saveSingleStudyModelFit=saveSingleStudyModelFit,
                               #silentOverwrite=silentOverwrite,
                               T0var=ctmaInitFit$argumentList$T0var,
                               useSV=ctmaInitFit$argumentList$useSV,
                               verbose=verbose,
                               #randomIntercepts=ctmaInitFit$argumentList$randomInterceptsSettings,
                               randomIntercepts=randomIntercepts),
                      error = function(e) problem <- TRUE
      )


      if ( (problem == FALSE) & is.list(fit) ) {
        all_minus2ll <- c(all_minus2ll, fit$summary$minus2ll)

        if (saveModelFits != FALSE) {
          saveRDS(fit, paste0(activeDirectory, saveModelFits, " ", i, " .rds"))
        }

        if (fit$summary$minus2ll < currentLL) {
          currentLL <- fit$summary$minus2ll
          bestFit <- fit
          usedStudyList <- ctmaInitFit$primaryStudyList
          usedTimeScale <- scaleTime
          usedScaleTI <- scaleTI
        }
      } else {
        all_minus2ll <- c(all_minus2ll, -999)
      }


    }
  }

  # ctmaFitFit
  if (!(is.null(ctmaFitFit))) {
    if (!(is(ctmaFitFit, "CoTiMAFit"))) {
      ErrorMsg <- "The ctmaFitFit object provided is not of class CoTiMAFit. Probably it was not created with ctmaFit."
      stop(ErrorMsg)
    }

    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == FALSE) & (randomRI == FALSE)) {
        Msg <- paste0("Both \"randomIV\" and \"randomRI\" are set FALSE (default). I use the same settings as used for ctmaFit(), which was
      \"randomIntercepts = ", ctmaFitFit$argumentList$randomIntercepts, "\", and
      \"indVarying = ", ctmaFitFit$argumentList$indVarying, ".\n")
        indVarying <- ctmaFitFit$argumentList$indVarying
        randomIntercepts <- ctmaFitFit$argumentList$randomIntercepts
        message(Msg)
      }
    }

    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == TRUE) & (randomRI == TRUE)) {
        Msg <- "Both \"randomIV\" and \"randomRI\" are set TRUE. Both are statistically equivalent, but since the latter is more difficult to fit I set \"randomRI = FALSE\". DO \"randomIV = FALSE\" to prevent this.\n"
        message(Msg)
        randomRI <- FALSE
      }
    }

    # CHD 14.7.2025
    if (is.null(randomIV)) {
      Msg <- "\"randomIV\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"indVarying = FALSE\". Take care.\n"
      randomIV <- FALSE
      indVarying <- FALSE
      message(Msg)
    }

    # CHD 14.7.2025
    if (is.null(randomRI)) {
      Msg <- "\"randomRI\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"randomIntercepts = FALSE\". Take care.\n"
      randomRI <- FALSE
      randomIntercepts <- FALSE
      message(Msg)
    }


    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    all_scaleTime <- all_customPar <- c()
    all_scaleTI <- c()
    all_usedStudyList <- c()
    all_randomIV <- all_randomRI <- c()
    warns <- errs <- list()

    for (i in 1:reFits) {
      cat(paste0("This is fit attempt #", i, " out of ", reFits, "re-fits."))
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      all_scaleTime <- c(all_scaleTime, scaleTime)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(customPar)) customPar <- ctmaFitFit$argumentList$customPar
      }
      #
      all_customPar <- c(all_customPar, customPar)
      if (randomScaleTI == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        scaleTI <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(scaleTI)) scaleTI <- ctmaFitFit$argumentList$scaleTI
      }
      #
      if (randomIV == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        indVarying <- c("MANIFEST", "CINT")[tmp1]
        all_randomIV <- c(all_randomIV, indVarying)
      } else {
        indVarying <- ctmaFitFit$argumentList$indVarying
        all_randomIV <- c(all_randomIV, indVarying)
      }
      #
      if (randomRI == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        randomIntercepts <- c("MANIFEST", "CINT")[tmp1]
        all_randomRI <- c(all_randomRI, randomIntercepts)
      } else {
        randomIntercepts <- ctmaFitFit$argumentList$randomIntercepts
        all_randomRI <- c(all_randomRI, randomIntercepts)
      }
      #
      if (shuffleStudyList == TRUE) {
        #
        tmpStudyList <- ctmaInitFit$studyList; length(tmpStudyList)
        studyNumbers <- unlist(lapply(tmpStudyList, function(x) x$originalStudyNo)); studyNumbers
        newStudyOrder <- sample(studyNumbers, length(studyNumbers), replace=FALSE); newStudyOrder
        newStudyList <- list()
        for (s in 1:length(tmpStudyList)) {
          newStudyList[[s]] <- tmpStudyList[[which(studyNumbers %in% newStudyOrder[s])]]
          # CHD 2.1.2025
          newStudyList[[s]]$studyNumber <- newStudyOrder[s]
        }
        ctmaInitFit$studyList <- newStudyList
        #
        tmpPrimaryStudyList <- ctmaInitFit$primaryStudyList
        newPrimaryStudyList <- list()
        for (s in 1:length(tmpPrimaryStudyList)) {
          if (names(tmpPrimaryStudyList[s]) %in% c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues", "studyNumbers",
                                                   "rawData", "source")) {
            newPrimaryStudyList[[s]] <- list()
            for (t in 1:length(tmpPrimaryStudyList[[s]])) {
              newPrimaryStudyList[[s]][[t]] <- tmpPrimaryStudyList[[s]][[which(studyNumbers %in% newStudyOrder[t])]]
            }
          } else {
            newPrimaryStudyList[[s]] <- tmpPrimaryStudyList[[s]]
          }
        }
        names(newPrimaryStudyList) <- names(tmpPrimaryStudyList)
        ctmaInitFit$primaryStudyList <- newPrimaryStudyList

        #
        tmpEmprawList <- ctmaInitFit$emprawList
        newEmprawList <- list()
        tmpStudyFitList <- ctmaInitFit$studyFitList
        newStudyFitList <- list()
        for (s in 1:length(tmpEmprawList)) {
          newEmprawList[[s]] <- tmpEmprawList[[which(studyNumbers %in% newStudyOrder[s])]]
          newStudyFitList[[s]] <- tmpStudyFitList[[which(studyNumbers %in% newStudyOrder[s])]]
        }
        ctmaInitFit$emprawList <- newEmprawList
        ctmaInitFit$studyFitList <- newStudyFitList

        all_usedStudyList <- c(all_usedStudyList, newStudyOrder)
      }

      if (!(is.null(randomScaleTime))) {
        Msg <- paste0("Argument scaleTime is set to: ", scaleTime, ".")
        message(Msg)
      }
      if (randomPar == TRUE) {
        Msg <- paste0("Argument customPar is set to: ", customPar, ".")
        message(Msg)
      }
      if (randomScaleTI == TRUE) {
        Msg <- paste0("Argument scaleTI is set to: ", scaleTI, ".")
        message(Msg)
      }
      if (shuffleStudyList == TRUE ) {
        tmp <- unlist(ctmaInitFit$primaryStudyList$studyNumbers)
        #tmp <- tmp[-length(tmp)]
        tmp <- paste(tmp, collapse=" ")
        Msg <- paste0("Order of studies in the shuffled study list is: ", tmp, ".")
        message(Msg)
      }

      if (is.null(finishsamples)) finishsamples <- ctmaFitFit$argumentList$finishsamples
      if (is.null(iter)) iter <- 5000

      fit <- tryCatch(withCallingHandlers(
        expr = ctmaFit(ctmaInitFit=ctmaInitFit,
                       primaryStudyList=ctmaInitFit$primaryStudyList,
                       cluster=ctmaFitFit$argumentList$cluster,
                       activeDirectory=activeDirectory,
                       activateRPB=ctmaFitFit$argumentList$activateRPB,
                       digits=ctmaFitFit$argumentList$digits,
                       drift=ctmaFitFit$argumentList$drift,
                       invariantDrift=ctmaFitFit$argumentList$invariantDrift,
                       moderatedDrift=ctmaFitFit$argumentList$moderatedDrift,
                       equalDrift=ctmaFitFit$argumentList$equalDrift,
                       mod.number=ctmaFitFit$argumentList$mod.number,
                       mod.type=ctmaFitFit$argumentList$mod.type,
                       mod.names=ctmaFitFit$argumentList$mod.names,
                       #indVarying=ctmaFitFit$argumentList$indVarying,
                       indVarying=indVarying,
                       coresToUse=coresToUse, # changed Aug 2023
                       sameInitialTimes=ctmaFitFit$argumentList$sameInitialTimes,
                       scaleTI=scaleTI,
                       scaleMod=ctmaFitFit$argumentList$scaleMod,
                       transfMod=ctmaFitFit$argumentList$transfMod,
                       scaleClus=ctmaFitFit$argumentList$scaleClus,
                       #scaleTime=ctmaFitFit$argumentList$scaleTime,
                       scaleTime=scaleTime,
                       optimize=ctmaFitFit$argumentList$optimize,
                       #nopriors=ctmaFitFit$argumentList$nopriors,
                       finishsamples=finishsamples,
                       iter=iter,
                       chains=ctmaFitFit$argumentList$chains,
                       verbose=verbose,
                       allInvModel=ctmaFitFit$argumentList$allInvModel,
                       customPar=customPar,
                       inits=ctmaFitFit$argumentList$inits,
                       modsToCompare=ctmaFitFit$argumentList$modsToCompare,
                       catsToCompare=ctmaFitFit$argumentList$catsToCompare,
                       driftsToCompare=ctmaFitFit$argumentList$driftsToCompare,
                       useSampleFraction=ctmaFitFit$argumentList$useSampleFraction,
                       T0means=ctmaFitFit$argumentList$T0means,
                       manifestMeans=ctmaFitFit$argumentList$manifestMeans,
                       CoTiMAStanctArgs=CoTiMAStanctArgs,
                       #randomIntercepts=ctmaFitFit$argumentList$randomIntercepts,
                       randomIntercepts=randomIntercepts,
                       manifestVars=ctmaFitFit$argumentList$manifestVars,
                       WEC=ctmaFitFit$argumentList$WEC,
                       priors=ctmaFitFit$argumentList$priors,
                       binaries=ctmaFitFit$argumentList$binaries,
                       T0var=ctmaFitFit$argumentList$T0var,
                       ind.mod.names=ctmaFitFit$argumentList$ind.mod.names,
                       ind.mod.number=ctmaFitFit$argumentList$ind.mod.number,
                       ind.mod.type=ctmaFitFit$argumentList$ind.mod.type,
                       cint=ctmaFitFit$argumentList$cint,
                       indVaryingT0=ctmaFitFit$argumentList$indVaryingT0,
                       fit=ctmaFitFit$argumentList$fit
        ),
        warning = function(w) {warns[[i]] <<- w}),
        error = function(e) {errs[[i]] <<- e}
      )

      if (is.list(fit)) { # do not do in case fit was not successfull CHD 14. Jan 2026
        all_minus2ll <- c(all_minus2ll, fit$summary$minus2ll)

        if (saveModelFits != FALSE) {
          saveRDS(fit, paste0(activeDirectory, saveModelFits, " ", i, " .rds"))
        }

        if (fit$summary$minus2ll < currentLL) {
          currentLL <- fit$summary$minus2ll
          bestFit <- fit
          usedStudyList <- ctmaInitFit$primaryStudyList
          usedTimeScale <- scaleTime
          usedScaleTI <- scaleTI
        }
      } # end if(is.list(fit))
    }
  }

  if (!(is.null(ctStanFit))) {
    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    all_scaleTime <- all_customPar <- c()
    all_scaleTI <- "Available only if CoTiMA models are optimized."
    all_usedStudyList <- "Available only if CoTiMA models are optimized."
    warns <- errs <- list()
    for (i in 1:reFits) {
      #i <- 1
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      all_scaleTime <- c(all_scaleTime, scaleTime)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(customPar)) customPar <- FALSE # ctmaFitFit$argumentList$customPar
      }
      all_customPar <- c(all_customPar, customPar)
      scaleTI <- FALSE
      #
      if (!(is.null(randomScaleTime))) {
        Msg <- paste0("Argument scaleTime is set to: ", scaleTime, ".")
        message(Msg)
      }
      if (randomPar == TRUE) {
        Msg <- paste0("Argument customPar is set to: ", customPar, ".")
        message(Msg)
      }
      #if (randomScaleTI == TRUE) {
      #  Msg <- paste0("Argument scaleTI is set to: ", scaleTI, ".")
      #  message(Msg)
      #}
      #if (shuffleStudyList == TRUE ) {
      #  tmp <- unlist(ctmaInitFit$primaryStudyList$studyNumbers)
      #  #tmp <- tmp[-length(tmp)]
      #  tmp <- paste(tmp, collapse=" ")
      #  Msg <- paste0("Order of studies in the shuffled study list is: ", tmp, ".")
      #  message(Msg)
      #}

      if (is.null(finishsamples)) finishsamples <- 1000
      if (is.null(iter)) iter <- 5000

      datalong <- cbind(ctStanFit$data$subject, matrix(ctStanFit$data$time, ncol=1),
                        ctStanFit$data$Y, ctStanFit$data$tdpreds)

      if (ctStanFit$standata$ntipred == 0) {
        colnames(datalong) <- colnames(ctStanFit$ctdatastruct)[1:ncol(datalong)]
      }

      if (ctStanFit$standata$ntipred > 0) {
        datalong <- cbind(datalong, matrix(ctStanFit$data$tipreds,
                                           nrow=length(ctStanFit$data$subject)))
        colnames(datalong) <- colnames(ctStanFit$ctdatastruct)
      }

      datalong[,2] <- datalong[,2] * scaleTime
      ctStanModel <- ctStanFit$ctstanmodelbase
      ctStanFitArgs$optimcontrol$finishsamples <- finishsamples

      fit <- tryCatch(withCallingHandlers(
        expr = ctsem::ctStanFit(datalong=datalong,
                                ctstanmodel=ctStanModel,
                                coresToUse=coresToUse,
                                finishsamples=finishsamples,
                                stanmodeltext = ctStanFitArgs$stanmodeltext,
                                iter=iter,
                                intoverstates = ctStanFitArgs$intoverstates,
                                binomial = ctStanFitArgs$binomial,
                                fit = ctStanFitArgs$fit,
                                intoverpop = ctStanFitArgs$intoverpop,
                                sameInitialTimes = ctStanFitArgs$sameInitialTimes,
                                stationary = ctStanFitArgs$stationary,
                                plot = ctStanFitArgs$plot,
                                derrind = ctStanFitArgs$derrind,
                                optimize = ctStanFitArgs$optimize,
                                optimcontrol = ctStanFitArgs$optimcontrol,
                                nlcontrol = ctStanFitArgs$nlcontrol,
                                nopriors = ctStanFitArgs$nopriors,
                                priors = ctStanFitArgs$priors,
                                chains = ctStanFitArgs$chains,
                                cores = coresToUse,
                                inits = ctStanFitArgs$inits,
                                compileArgs = ctStanFitArgs$compileArgs,
                                forcerecompile = ctStanFitArgs$forcerecompile,
                                saveCompile = ctStanFitArgs$saveCompile,
                                savescores = ctStanFitArgs$savescores,
                                savesubjectmatrices = ctStanFitArgs$savesubjectmatrices,
                                saveComplexPars = ctStanFitArgs$saveComplexPars,
                                gendata = ctStanFitArgs$gendata,
                                control = ctStanFitArgs$control,
                                verbose = verbose,
                                vb = ctStanFitArgs$vb
        ),
        warning = function(w) {warns[[i]] <<- w}),
        error = function(e) {errs[[i]] <<- e}
      )
      fit$summary <- summary(fit)
      fit$summary$minus2ll <- fit$summary$logposterior * -2

      all_minus2ll <- c(all_minus2ll, fit$summary$minus2ll)

      if (saveModelFits != FALSE) {
        saveRDS(fit, paste0(activeDirectory, saveModelFits, " ", i, " .rds"))
      }

      if (fit$summary$minus2ll < currentLL) {
        currentLL <- fit$summary$minus2ll
        bestFit <- fit
        usedStudyList <- "Available only if CoTiMA models are optimized."
        usedTimeScale <- scaleTime
        usedScaleTI <- "Available only if CoTiMA models are optimized."
      }
    }
    resultsSummary <- "Available only if CoTiMA models are optimized."
  } else {
    resultsSummary <- bestFit$studyFitList[[1]]$resultsSummary
  }

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  #usedStudyList=ctmaInitFit$primaryStudyList,
                  usedStudyList=ctmaFitToPrep(bestFit, reUseEmprawData=TRUE),
                  #usedTimeScale=usedTimeScale,
                  usedTimeScale=bestFit$argumentList$scaleTime,
                  #usedScaleTI=usedScaleTI,
                  usedScaleTI=bestFit$argumentList$scaleTI,
                  #resultsSummary=bestFit$studyFitList[[1]]$resultsSummary
                  randomIV = randomIV,
                  randomRI = randomRI,
                  resultsSummary=resultsSummary,
                  all_scaleTime = all_scaleTime,
                  all_customPar = all_customPar,
                  all_scaleTI = all_scaleTI,
                  all_randomIV = all_randomIV,
                  all_randomRI = all_randomRI,
                  all_usedStudyList = all_usedStudyList

  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}
