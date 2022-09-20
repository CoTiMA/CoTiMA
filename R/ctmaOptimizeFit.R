#' ctmaOptimizeFit
#'
#' @description Replaces deprecated \code{\link{ctmaOptimizeInit}}, which was limited to initial fitting
#' (i.e., applies \code{\link{ctmaInit}}) of a primary study reFits times to capitalize on chance for obtaining
#' a hard-to-find optimal fit.
#' Now, optimizing a CoTiMA model generated with \code{\link{ctmaFit}} can also be done.
#' Using \code{\link{ctmaOptimizeFit}} could be helpful if a model yields out-of-range estimates, which could happen if the fitting
#' algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
#' \code{\link{ctmaOptimizeFit}} is like gambling, hoping that at least one set of starting values (the number it tries is specified in the reFits argument)
#' enables finding the global optimal fit. On unix-like machines (e.g. MacOS), this could be done in parallel mode if coresToUse > 1.
#'
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}
#' @param activeDirectory activeDirectory
#' @param problemStudy number (position in list) where the problem study in primaryStudies is found
#' @param reFits how many reFits should be done
#' @param n.latent number of latent variables of the model (hast to be specified)!
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param indVarying control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means
#' @param randomScaleTime lower and upper limit of uniform distribution from which timeScale argument for ctmaInit is uniformly shuffled (integer)
#' @param customPar logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param manifestMeans Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function
#' @param checkSingleStudyResults displays estimates from single study 'ctsem' models and waits for user input to continue.
#' @param CoTiMAFit a object fitted with \code{\link{ctmaFit}}
#' @param CoTiMAInitFit the ctmaInitFit object that was used to create the CoTiMAFit object with \code{\link{ctmaFit}}
#' @param randomPar logical. Overrides arguments used fo customPar and randomly selects customPar either TRUE or FALSE
#' Useful to check estimates before they are saved.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom RPushbullet pbPost
#' @importFrom stats runif
#'
#' @note All but one of multiple cores are used on unix-type machines for parallel fitting
#' @note During fitting, not output is generated. Be patient.
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
ctmaOptimizeFit <- function(primaryStudies=NULL,
                            activeDirectory=NULL,
                            problemStudy=NULL,
                            reFits=NULL,
                            finishsamples=NULL,
                            n.latent=NULL,
                            coresToUse=c(1),
                            indVarying=FALSE,
                            randomScaleTime=c(1,1),
                            activateRPB=FALSE,
                            checkSingleStudyResults=FALSE,
                            customPar=FALSE,
                            T0means=0,
                            manifestMeans=0,
                            CoTiMAStanctArgs=NULL,
                            CoTiMAFit=NULL,
                            CoTiMAInitFit=NULL,
                            randomPar=FALSE)
{

  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################

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

  # CHD deleted 20. Sep 2022
  #if (.Platform$OS.type == "unix") {
  #  doParallel::registerDoParallel(coresToUse)
  #} else {
  #  doParallel::registerDoParallel(1)
  #}

  # CHD added 20 Sep 2022
  if  (length(coresToUse) > 0) {
    if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
  }
  #
  if (coresToUse >= parallel::detectCores()) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    coresToUse <- parallel::detectCores() - 1
    Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.\n"
    message(Msg)
  }
  #
  if (doPar > 1) {
    doParallel::registerDoParallel(coresToUse)
    '%dopar%' <- foreach::'%dopar%'
  }


  if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples

  # Added 17. Aug 2022
  tmp1 <- names(CoTiMA::CoTiMAStanctArgs) %in% names(CoTiMAStanctArgs); tmp1
  tmp2 <- CoTiMA::CoTiMAStanctArgs
  if (!(is.null(CoTiMAStanctArgs))) tmp2[tmp1] <- CoTiMAStanctArgs
  CoTiMAStanctArgs <- tmp2

  ########################################################################################################################

  '%dopar%' <- foreach::'%dopar%'

  if (randomScaleTime[2] < randomScaleTime[1]) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    ErrorMsg <- "\nrandomScaleTime[1] has to be <= randomScaleTime[2]! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if( (!(is.null(CoTiMAFit))) & (!(is.null(primaryStudies))) ) {
    ErrorMsg <- "Arguments for both CoTiMAFit and primaryStudies were provided. Only one out of the two can be chosen!"
    stop(ErrorMsg)
  }

  if( (!(is.null(CoTiMAFit))) & ((is.null(CoTiMAInitFit))) ) {
    ErrorMsg <- "Argument for CoTiMAFit was provided but not for CoTiMAInitFit. Need the latter, too!"
    stop(ErrorMsg)
  }

  if( (!(is.null(CoTiMAFit))) & (!(is.null(CoTiMAInitFit)))  & (CoTiMAFit$argumentList$ctmaInitFit != deparse(substitute(CoTiMAInitFit)))  ) {
    ErrorMsg <- paste0("The wrong CoTiMAInitFit object was provided. I need ",  CoTiMAFit$argumentList$ctmaInitFit, "!")
    stop(ErrorMsg)
  }

  # INIT Fit
  if(is.null(CoTiMAFit)) {
    ErrorMsg <- "arguments are missing"
    if (is.null(primaryStudies) | is.null(problemStudy) | is.null(reFits) | is.null(activeDirectory) | is.null(n.latent)  ) stop(ErrorMsg)

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

    #
    # adaptations for dealing with raw data
    # taken out 9. Aug. 2022
    #if (!(is.null(newStudyList$rawData[[1]]$studyNumbers))) {
    #  newStudyList$rawData[[1]]$studyNumbers <- 1
    #  newStudyList$studyNumbers <- 1 # otherwise raw data will not be loaded
    #  newStudyList$deltas <- unlist(newStudyList$deltas)
    #}
    #newStudyList

    # parallel re-fitting of problem study
    allfits <- foreach::foreach(i=1:reFits) %dopar% {
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar = c(TRUE, FALSE)[tmp1]
      }
      fits <- ctmaInit(newStudyList, coresToUse = 1, n.latent=n.latent,
                       indVarying = indVarying,
                       scaleTime = scaleTime,
                       activeDirectory = activeDirectory,
                       checkSingleStudyResults=checkSingleStudyResults,
                       customPar=customPar,
                       T0means=0,
                       manifestMeans=0)
      return(fits)
    }
  }

  # CoTiMAFit
  if(!(is.null(CoTiMAFit))) {
    if (class(CoTiMAFit) != "CoTiMAFit") {
      ErrorMsg <- "The CoTiMAfit object provided is not of class CoTiMAFit. Probably it was not created with ctmaFit."
      stop(ErrorMsg)
    }

    allfits <- foreach::foreach(i=1:reFits) %dopar% {
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar = c(TRUE, FALSE)[tmp1]
      }
      #ctmaInitFitBackup <- CoTiMAInitFit
      #primaryStudyListBackup <- CoTiMAInitFit$primaryStudyList
      ##
      #for (l in 1:length(CoTiMAFit$argumentList)) {
      #  assign(names(CoTiMAFit$argumentList)[[l]], CoTiMAFit$argumentList[[l]]) #, envir = as.environment(-1))
      #  }
      ##
      #ctmaInitFit <- ctmaInitFitBackup
      #primaryStudyList <- primaryStudyListBackup
      #
      fits <- ctmaFit(ctmaInitFit=CoTiMAInitFit,
                      primaryStudyList=CoTiMAInitFit$primaryStudyList,
                      cluster=CoTiMAFit$argumentList$cluster,
                      activeDirectory=activeDirectory,
                      activateRPB=CoTiMAFit$argumentList$activateRPB,
                      digits=CoTiMAFit$argumentList$digits,
                      drift=CoTiMAFit$argumentList$drift,
                      invariantDrift=CoTiMAFit$argumentList$invariantDrift,
                      moderatedDrift=CoTiMAFit$argumentList$moderatedDrift,
                      equalDrift=CoTiMAFit$argumentList$equalDrift,
                      mod.number=CoTiMAFit$argumentList$mod.number,
                      mod.type=CoTiMAFit$argumentList$mod.type,
                      mod.names=CoTiMAFit$argumentList$mod.names,
                      #n.manifest=0,
                      indVarying=CoTiMAFit$argumentList$indVarying,
                      coresToUse=CoTiMAFit$argumentList$coresToUse,
                      scaleTI=CoTiMAFit$argumentList$scaleTI,
                      scaleMod=CoTiMAFit$argumentList$scaleMod,
                      transfMod=CoTiMAFit$argumentList$transfMod,
                      scaleClus=CoTiMAFit$argumentList$scaleClus,
                      scaleTime=CoTiMAFit$argumentList$scaleTime,
                      optimize=CoTiMAFit$argumentList$optimize,
                      nopriors=CoTiMAFit$argumentList$nopriors,
                      finishsamples=CoTiMAFit$argumentList$finishsamples,
                      iter=CoTiMAFit$argumentList$iter,
                      chains=CoTiMAFit$argumentList$chains,
                      verbose=CoTiMAFit$argumentList$verbose,
                      allInvModel=CoTiMAFit$argumentList$allInvModel,
                      customPar=CoTiMAFit$argumentList$customPar,
                      inits=CoTiMAFit$argumentList$inits,
                      modsToCompare=CoTiMAFit$argumentList$modsToCompare,
                      catsToCompare=CoTiMAFit$argumentList$catsToCompare,
                      driftsToCompare=CoTiMAFit$argumentList$driftsToCompare,
                      useSampleFraction=CoTiMAFit$argumentList$useSampleFraction,
                      T0means=CoTiMAFit$argumentList$T0means,
                      manifestMeans=CoTiMAFit$argumentList$manifestMeans,
                      CoTiMAStanctArgs=CoTiMAFit$argumentList$CoTiMAStanctArgs
      )
      return(fits)
    }
  }


  all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
  bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
  bestFit <- allfits[[bestFit]]

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  resultsSummary=bestFit$studyFitList[[1]]$resultsSummary
  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}


