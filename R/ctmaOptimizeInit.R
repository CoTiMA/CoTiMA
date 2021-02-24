#' ctmaOptimzeInit
#'
#' @description Initial fitting (i.e., applies ctmaInit) to a primary study reFit times to capitalize on chance for obtaining a hard-to-find optimal fit.
#' This could be very helpful if a primary yields out-of-range estimates, which could happen if the fitting algorithm unfortunately used random start
#' values that resulted in a locally but not globally optimal fit. Essentially, ctmaOptimzeInit is like gambling, hoping that at leas one set of starting
#' values (the number is tries is specified in the reFits argument) eneables finding the global optimal fit. On unix-like machines (e.g. MacOS), this
#' could be done in parallel mode if coresToUse > 1.
#'
#' @param oldStudyList oldStudyList
#' @param activeDirectory activeDirectory
#' @param problemStudy problemStudy
#' @param reFits reFits
#' @param n.latent n.latent
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#' @param activateRPB  set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param checkSingleStudyResults Displays estimates from single study ctsem models and waits for user input to continue. Useful to check estimates before they are saved.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom RPushbullet pbPost
#'
#' @note All but one of multiple cores are used on unix-type machines for parallel fitting
#' @note During fitting, not output is generated. Be patient.
#'
#' @examples
#' \donttest{
#' optimFit313 <- ctmaOptimzeInit(oldStudyList=CoTiMAstudyList_3,
#'                                activeDirectory="/Users/cdormann/tmp/",
#'                                problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
#'                                reFits=10,
#'                                n.latent=2)
#' summary(optimFit313)
#' }
#'
#' @export ctmaOptimzeInit
#'
#' @return returns a list with bestFit (= the best fit achieved), all_minus2ll (= all -2ll values for all fitted models), and summary, which
#' is printed if the summary function is applied to the returned object, and which shows the summary information of the ctsem model with the
#' best fit.
#'
ctmaOptimzeInit <- function(oldStudyList=NULL,
                            activeDirectory=NULL,
                            problemStudy=NULL,
                            reFits=NULL,
                            n.latent=NULL,
                            coresToUse=c(1),
                            activateRPB=FALSE,
                            checkSingleStudyResults=FALSE)
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
      cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
    }
  }

  if (.Platform$OS.type == "unix") {
    doParallel::registerDoParallel(coresToUse)
  } else {
    doParallel::registerDoParallel(1)
  }

  ########################################################################################################################

  '%dopar%' <- foreach::'%dopar%'

  if (is.null(oldStudyList) | is.null(problemStudy) | is.null(reFits) | is.null(activeDirectory) | is.null(n.latent)  ) stop("arguments are missing")

  # create new study list with a single problem study only
  newStudyList <- list()
  for (j in 1:(length(names(oldStudyList)))) {
    tmp1 <- lapply(oldStudyList[[names(oldStudyList)[j]]], function(x) x)
    if (length(tmp1) == 1) {
      newStudyList[[j]] <- tmp1
    } else {
      newStudyList[[j]] <- list(tmp1[[problemStudy]])
    }
  }
  names(newStudyList) <- names(oldStudyList)
  newStudyList$n.studies <- 1

  # parallel re-fitting of problem study
  allfits <- foreach::foreach(i=1:reFits) %dopar% {
    fits <- ctmaInit(newStudyList, coresToUse = 1, n.latent=n.latent,
                     activeDirectory = activeDirectory,
                     checkSingleStudyResults=checkSingleStudyResults)
    return(fits)
  }

  all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
  bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
  bestFit <- allfits[[bestFit]]

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary)
  class(results) <- "CoTiMAFit"

  invisible(results)
}


