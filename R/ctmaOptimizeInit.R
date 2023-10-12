#' ctmaOptimizeInit
#'
#' @description Initial fitting (i.e., applies \code{\link{ctmaInit}}) to a primary study reFit times to capitalize on chance for obtaining
#' a hard-to-find optimal fit. This could be very helpful if a primary yields out-of-range estimates, which could happen if the fitting
#' algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
#' ctmaOptimizeInit is like gambling, hoping that at leas one set of starting values (the number is tries is specified in the reFits argument)
#' eneables finding the global optimal fit. On unix-like machines (e.g. MacOS), this could be done in parallel mode if coresToUse > 1.
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
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param manifestVars define the error variances of the manifests with a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param checkSingleStudyResults displays estimates from single study 'ctsem' models and waits for user input to continue.
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
#' optimFit313 <- ctmaOptimizeInit(primaryStudies=CoTiMAstudyList_3,
#'                                 activeDirectory="/Users/tmp/",  # adapt!
#'                                 problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
#'                                 reFits=10,
#'                                 n.latent=2)
#' summary(optimFit313)
#' }
#'
#' @export ctmaOptimizeInit
#'
#' @return returns a list with bestFit (= the best fit achieved), all_minus2ll (= all -2ll values for all fitted models), and summary, which
#' is printed if the summary function is applied to the returned object, and which shows the summary information of the ctsem model with the
#' best fit.
#'
ctmaOptimizeInit <- function(primaryStudies=NULL,
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
                             manifestVars=NULL,
                             CoTiMAStanctArgs=NULL,
                             scaleTime=NULL)
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

  if (.Platform$OS.type == "unix") {
    doParallel::registerDoParallel(coresToUse)
  } else {
    doParallel::registerDoParallel(1)
  }

  if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
  if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime

  # Added 17. Aug 2022
  tmp1 <- names(CoTiMA::CoTiMAStanctArgs) %in% names(CoTiMAStanctArgs); tmp1
  tmp2 <- CoTiMA::CoTiMAStanctArgs
  if (!(is.null(CoTiMAStanctArgs))) tmp2[tmp1] <- CoTiMAStanctArgs
  CoTiMAStanctArgs <- tmp2

  ########################################################################################################################

  '%dopar%' <- foreach::'%dopar%'
  ErrorMsg <- "arguments are missing"
  if (is.null(primaryStudies) | is.null(problemStudy) | is.null(reFits) | is.null(activeDirectory) | is.null(n.latent)  ) stop(ErrorMsg)

  if (randomScaleTime[2] < randomScaleTime[1]) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    ErrorMsg <- "\nrandomScaleTime[1] has to be <= randomScaleTime[2]! \nGood luck for the next try!"
    stop(ErrorMsg)
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
    fits <- ctmaInit(primaryStudies=newStudyList,
                     coresToUse = 1, n.latent=n.latent,
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
  all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
  bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
  bestFit <- allfits[[bestFit]]

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  resultsSummary=bestFit$studyFitList[[1]]$resultsSummary
  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}


