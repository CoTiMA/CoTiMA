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
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param activeDirectory activeDirectory
#' @param checkSingleStudyResults displays estimates from single study 'ctsem' models and waits for user input to continue.
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function
#' @param ctmaFitFit a object fitted with \code{\link{ctmaFit}}
#' @param ctmaInitFit the ctmaInitFit object that was used to create the ctmaFitFit object with \code{\link{ctmaFit}}
#' @param customPar logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param indVarying control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param manifestMeans Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param manifestVars define the error variances of the manifests within a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest. Useful to check estimates before they are saved.
#' @param n.latent number of latent variables of the model (hast to be specified)!
#' @param posLL logical. Allows (default = TRUE) of positive loglik (neg -2ll) values
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}
#' @param problemStudy number (position in list) where the problem study in primaryStudies is found
#' @param randomPar logical. Overrides arguments used fo customPar and randomly selects customPar either TRUE or FALSE
#' @param randomScaleTime lower and upper limit of uniform distribution from which timeScale argument for ctmaInit is uniformly shuffled (integer)
#' @param reFits how many reFits should be done
#' @param scaleMod scale moderator variables - TRUE (default) recommended for continuous and categorical moderators, to separate withing and betwen efeccts
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param transfMod more general option to change moderator values. A vector as long as number of moderators analyzed (e.g., c("mean(x)", "x - median(x)"))

#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom foreach %dopar%
#' @importFrom RPushbullet pbPost
#' @importFrom stats runif
#' @importFrom methods is
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
ctmaOptimizeFit <- function(activateRPB=FALSE,
                            activeDirectory=NULL,
                            checkSingleStudyResults=FALSE,
                            coresToUse=c(2),
                            CoTiMAStanctArgs=NULL,
                            ctmaFitFit=NULL,
                            ctmaInitFit=NULL,
                            customPar=FALSE,
                            finishsamples=NULL,
                            indVarying=FALSE,
                            lambda=NULL,
                            manifestMeans=0,
                            manifestVars=NULL,
                            n.latent=NULL,
                            posLL=TRUE,
                            primaryStudies=NULL,
                            problemStudy=NULL,
                            randomPar=FALSE,
                            randomScaleTime=c(1,1),
                            reFits=NULL,
                            scaleMod=NULL,
                            scaleTime=NULL,
                            T0means=0,
                            transfMod=NULL
                            )
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

  # CHD ADDED Aug. 2023
  if (!(is.null(reFits))) {
    parProces <- min(c(coresToUse, reFits))
    coresToUse <- coresToUse%/%reFits # integer division
    if (coresToUse < 1) coresToUse <- 1
  } else {
    parProces <- coresToUse
  }

  # CHD added 20 SEP 2020, changed Aug 2023
  #myCluster <- parallel::makeCluster(coresToUse)
  myCluster <- parallel::makeCluster(parProces)
  on.exit(parallel::stopCluster(myCluster))
  # CHD deleted 20. Sep 2022
  #if (.Platform$OS.type == "unix") {
    #doParallel::registerDoParallel(coresToUse)
    doParallel::registerDoParallel(myCluster)
    #} else {
  #  doParallel::registerDoParallel(1)
  #}


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

  if( (!(is.null(ctmaFitFit))) & (!(is.null(primaryStudies))) ) {
    ErrorMsg <- "Arguments for both ctmaFitFit and primaryStudies were provided. Only one out of the two can be chosen!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & ((is.null(ctmaInitFit))) ) {
    ErrorMsg <- "Argument for ctmaFitFit was provided but not for ctmaInitFit. Need the latter, too!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & (!(is.null(ctmaInitFit))) ) {
    if (ctmaFitFit$argumentList$ctmaInitFit != deparse(substitute(ctmaInitFit)))  {
    ErrorMsg <- paste0("The wrong ctmaInitFit object was provided. I need ",  ctmaFitFit$argumentList$ctmaInitFit, "!")
    stop(ErrorMsg)
  }
  }

  # Moderator Checks Moved to Sectioon where ctmaInit is optimized (not relevant if ctmaFit is optimized) # CHD Auf 2023


  # INIT Fit
  if(is.null(ctmaFitFit)) {
    # CHD changed 21 SEP 2022
    ErrorMsg <- "argument primaryStudies is missing"
    if (is.null(primaryStudies))  stop(ErrorMsg)
    ErrorMsg <- "argument problemStudy is missing"
    if (is.null(problemStudy)) stop(ErrorMsg)
    ErrorMsg <- "argument reFits is missing"
    if (is.null(reFits)) stop(ErrorMsg)
    ErrorMsg <- "argument activeDirectory is missing"
    if (is.null(activeDirectory)) stop(ErrorMsg)
    ErrorMsg <- "argument n.latent is missing"
    if (is.null(n.latent)) stop(ErrorMsg)

    # Moderator Checks Moved from above into ctmaInit-optimizinf section (not relevant if ctmaFit is optimized)
    # Neither relevant for ctmaInit nor ctmaFit, therefore de-activated
    #{ # check if scaleMod is not used in combination with transfMod
    #  if ( (!(is.null(scaleMod))) & (!(is.null(transfMod))) ) {
    #    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    #    ErrorMsg <- "The arguments scaleMod and transfMod cannot be used in combination. Set one of them NULL (leave out)."
    #    stop(ErrorMsg)
    #  }
    #}

    #n.ind.moderators <- 0 # shortcut until n.ind.moderators is enabled as argument
    #if (n.ind.moderators == 0) { # proceed if only moderators at the study level are used
    #  n.moderators <- length(mod.number); n.moderators
    #  { # check if transfMod is as long as n.moderators
    #    if ( (!(is.null(transfMod))) ) {
    #      if ( length(transfMod) != n.moderators ) {
    #        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    #        ErrorMsg <- "More transformations for moderators (transfMod) provided than moderators."
    #        stop(ErrorMsg)
    #      }
    #    }
    #  }
    #  if (n.moderators > 0) {
    #    currentModerators <- matrix(as.numeric(unlist(lapply(ctmaInitFit$studyList, function(extract) extract$moderators[mod.number]))), ncol=n.moderators); currentModerators
    #    #if (!(is.null(primaryStudyList))) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudyList$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE)
    #    #currentModerators
    #    #if (!is.null(primaryStudyList)) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudyList$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE); currentModerators
    #    if (is.na((currentModerators[length(currentModerators)])[[1]][1])) currentModerators <- currentModerators[-dim(currentModerators)[1],]; currentModerators
    #    if (is.null(dim(currentModerators)[1])) currentModerators <- matrix(currentModerators, ncol=1); currentModerators
    #    #table(currentModerators)

    #    if (any(is.na(currentModerators)) == TRUE) {
    #      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    #      ErrorMsg <- "\nAt least one of the primary studies does not have a valid value for the requested moderator. \nGood luck for the next try!"
    #      stop(ErrorMsg)
    #    }
    #    if (var(currentModerators) == 0) {
    #      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    #      ErrorMsg <- "\nModerator is constant across cases.\nGood luck for the next try!"
    #      stop(ErrorMsg)
    #    }
    #  }
    #}


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
      fits <- ctmaInit(primaryStudies=newStudyList,
                       #coresToUse = 1,
                       coresToUse = coresToUse, # changed Aug 2023
                       n.latent=n.latent,
                       indVarying = indVarying,
                       scaleTime = scaleTime,
                       activeDirectory = activeDirectory,
                       checkSingleStudyResults=checkSingleStudyResults,
                       customPar=customPar,
                       T0means=T0means,
                       manifestMeans=manifestMeans,
                       manifestVars=manifestVars)
      return(fits)
    }
  }

  # ctmaFitFit
  if(!(is.null(ctmaFitFit))) {
    if (!(is(ctmaFitFit, "CoTiMAFit"))) {
      ErrorMsg <- "The ctmaFitFit object provided is not of class CoTiMAFit. Probably it was not created with ctmaFit."
      stop(ErrorMsg)
    }

    allfits <- foreach::foreach(i=1:reFits) %dopar% {
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar = c(TRUE, FALSE)[tmp1]
      }

      fits <- ctmaFit(ctmaInitFit=ctmaInitFit,
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
                      #n.manifest=0,
                      indVarying=ctmaFitFit$argumentList$indVarying,
                      #coresToUse=ctmaFitFit$argumentList$coresToUse,
                      #coresToUse=1,
                      coresToUse=coresToUse, # changed Aug 2023
                      sameInitialTimes=ctmaFitFit$argumentList$sameInitialTimes,
                      scaleTI=ctmaFitFit$argumentList$scaleTI,
                      scaleMod=ctmaFitFit$argumentList$scaleMod,
                      transfMod=ctmaFitFit$argumentList$transfMod,
                      scaleClus=ctmaFitFit$argumentList$scaleClus,
                      scaleTime=ctmaFitFit$argumentList$scaleTime,
                      optimize=ctmaFitFit$argumentList$optimize,
                      nopriors=ctmaFitFit$argumentList$nopriors,
                      finishsamples=ctmaFitFit$argumentList$finishsamples,
                      iter=ctmaFitFit$argumentList$iter,
                      chains=ctmaFitFit$argumentList$chains,
                      verbose=ctmaFitFit$argumentList$verbose,
                      allInvModel=ctmaFitFit$argumentList$allInvModel,
                      customPar=ctmaFitFit$argumentList$customPar,
                      inits=ctmaFitFit$argumentList$inits,
                      modsToCompare=ctmaFitFit$argumentList$modsToCompare,
                      catsToCompare=ctmaFitFit$argumentList$catsToCompare,
                      driftsToCompare=ctmaFitFit$argumentList$driftsToCompare,
                      useSampleFraction=ctmaFitFit$argumentList$useSampleFraction,
                      T0means=ctmaFitFit$argumentList$T0means,
                      manifestMeans=ctmaFitFit$argumentList$manifestMeans,
                      CoTiMAStanctArgs=ctmaFitFit$argumentList$CoTiMAStanctArgs
      )
      return(fits)
    }
  }


  all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
  # CHD added 27 SEP 2022 to prevent neg -2ll fits
  if(posLL == FALSE) {
    if (all(all_minus2ll < 0)) {
      ErrorMsg <- "\n All loglik values > 0, but you provided the argument posLL=FALSE, so no fit confirmed your expectations and I had to stop!"
      stop(ErrorMsg)
    }
    all_minus2ll <- all_minus2ll[-(which(all_minus2ll < 0))]
    }
  #
  bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
  bestFit <- allfits[[bestFit]]

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  resultsSummary=bestFit$studyFitList[[1]]$resultsSummary
  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}


