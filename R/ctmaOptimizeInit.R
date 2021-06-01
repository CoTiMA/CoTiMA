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
#' @param customPar logical. Leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (could help with ctsem > 3.4)
#' @param checkSingleStudyResults displays estimates from single study 'ctsem' models and waits for user input to continue.
#' Useful to check estimates before they are saved.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom RPushbullet pbPost
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
                             n.latent=NULL,
                             coresToUse=c(1),
                             activateRPB=FALSE,
                             checkSingleStudyResults=FALSE,
                             customPar=FALSE)
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

  ########################################################################################################################

  '%dopar%' <- foreach::'%dopar%'
  ErrorMsg <- "arguments are missing"
  if (is.null(primaryStudies) | is.null(problemStudy) | is.null(reFits) | is.null(activeDirectory) | is.null(n.latent)  ) stop(ErrorMsg)

  # create new study list with a single problem study only
  listElements <- names(primaryStudies); listElements
  newStudyList <- as.list(listElements)
  validElements <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues", "studyNumbers", "rawData", "empMeans", "empVars",
                     "source", "ageM", "malePercent", "occupation", "country", "alphas", "targetVariables", "recodeVariables", "combineVariables",
                     "combineVariablesNames", "missingVariables") #, "n.studies", "summary", "excelSheets", "plot.type")
  counter <- 0
  for (i in listElements) {
    counter <- counter + 1
    if (i %in% validElements) {
      if (i %in% c("pairwiseNs", "empcovs", "rawData", "deltas")) {
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


  # adaptations for dealing with raw data
  if (!(is.na(newStudyList$rawData[[1]]$studyNumbers))) {
    newStudyList$rawData[[1]]$studyNumbers <- 1
    newStudyList$studyNumbers <- 1 # otherwise raw data will not be loaded
    newStudyList$deltas <- unlist(newStudyList$deltas)
  }

  # parallel re-fitting of problem study
  allfits <- foreach::foreach(i=1:reFits) %dopar% {
    fits <- ctmaInit(newStudyList, coresToUse = 1, n.latent=n.latent,
                     activeDirectory = activeDirectory,
                     checkSingleStudyResults=checkSingleStudyResults,
                     customPar=customPar)
    return(fits)
  }

  all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
  bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
  bestFit <- allfits[[bestFit]]

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary)
  class(results) <- "CoTiMAFit"

  invisible(results)
}

