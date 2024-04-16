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
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function
#' @param ctmaFitFit a object fitted with \code{\link{ctmaFit}}
#' @param ctmaInitFit the ctmaInitFit object that was used to create the ctmaFitFit object with \code{\link{ctmaFit}}
#' @param customPar logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param iter number of iterations (default = 5000)
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}
#' @param problemStudy number (position in list) where the problem study in primaryStudies is found
#' @param randomPar logical (default = FALSE). Overrides arguments used for customPar and randomly sets customPar either TRUE or FALSE
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
#'
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
                            #checkSingleStudyResults=FALSE,
                            coresToUse=c(2),
                            CoTiMAStanctArgs=NULL,
                            ctmaFitFit=NULL,
                            ctmaInitFit=NULL,
                            customPar=FALSE,
                            finishsamples=NULL,
                            iter=5000,
                            #indVarying=NULL,
                            #lambda=NULL,
                            #manifestMeans=0,
                            #manifestVars=NULL,
                            #n.latent=NULL,
                            primaryStudies=NULL,
                            problemStudy=NULL,
                            randomPar=FALSE,
                            randomScaleTI=FALSE,
                            randomScaleTime=c(1,1),
                            saveModelFits=FALSE,
                            shuffleStudyList=FALSE,
                            reFits=NULL,
                            scaleTime=NULL,
                            scaleTI=NULL,
                            #T0means=0,
                            #parallel=FALSE,
                            verbose=1
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

    # Dealing with CoTiMAStanctArgs
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
  if (is.null(ctmaFitFit)) {
    # CHD changed 21 SEP 2022
    ErrorMsg <- "argument primaryStudies is missing"
    if (is.null(primaryStudies))  stop(ErrorMsg)
    ErrorMsg <- "argument problemStudy is missing"
    if (is.null(problemStudy)) stop(ErrorMsg)
    ErrorMsg <- "argument reFits is missing"
    if (is.null(reFits)) stop(ErrorMsg)
    ErrorMsg <- "argument activeDirectory is missing"
    if (is.null(activeDirectory)) stop(ErrorMsg)
    #ErrorMsg <- "argument n.latent is missing"
    #if (is.null(n.latent)) stop(ErrorMsg)


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
    for (i in 1:reFits) {
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

      if (is.null(finishsamples)) finishsamples <- ctmaInitFit$argumentList$finishsamples
      #if (is.null(iter)) iter <- 5000

      # CHD 12.4.24
      #if (is.null(indVarying)) indVarying <- ctmaFitFit$argumentList$indVarying

      problem <- FALSE
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
                               indVarying = ctmaInitFit$argumentList$indVarying,
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
                               randomIntercepts=ctmaInitFit$argumentList$randomInterceptsSettings),
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

    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    for (i in 1:reFits) {
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(customPar)) customPar <- ctmaFitFit$argumentList$customPar
      }
      #
      if (randomScaleTI == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        scaleTI <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(scaleTI)) scaleTI <- ctmaFitFit$argumentList$scaleTI
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

      fit <- ctmaFit(ctmaInitFit=ctmaInitFit,
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
                     indVarying=ctmaFitFit$argumentList$indVarying,
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
                     randomIntercepts=ctmaFitFit$argumentList$randomIntercepts,
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
      )

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
    }
  }


  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  usedStudyList=ctmaInitFit$primaryStudyList, usedTimeScale=usedTimeScale, usedScaleTI=usedScaleTI,
                  resultsSummary=bestFit$studyFitList[[1]]$resultsSummary
  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}
