#######################################################################################################################
############################## CoTiMA Multiple Drift Effects Invariant Across primary Studies #########################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   ctmaMultipleFit=ctmaMultipleFitMX
#   #type="mx"
#   activeDirectory=NULL
#   equalDrift=NULL
#   resultsFilePrefix="CoTiMAEqualDrift"
#   saveFilePrefix="CoTiMAEqualDrift"
#   activateRPB=FALSE
#   #checkSingleStudyResults=TRUE
#   retryattempts=30
#   refits=1
#   NPSOL=FALSE
#   coresToUse=c(-1)
#   fullDriftStartValues=NULL
#   compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2]
#   useCTMultiGroupAlt=FALSE
#   confidenceIntervals=FALSE
#   saveEqualDriftModelFit=saveFilePrefix
#   ##### ENTER DEBUG INFO HERE #######
#   silentOverwrite=FALSE
#   #CoTiMAprepFit=CoTiMAprepFit1
#   #useCTMultiGroupAlt=TRUE
#   #confidenceIntervals=FALSE
#   digits=4
#   CoTiMAStanctArgs=list(test=TRUE, scaleTI=TRUE, scaleTime=1/1,
#                         savesubjectmatrices=FALSE, verbose=1,
#                         datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
#                         iter=1000, intoverstates=TRUE,
#                         binomial=FALSE, fit=TRUE,
#                         intoverpop=FALSE, stationary=FALSE,
#                         plot=FALSE, derrind="all",
#                         optimize=TRUE, optimcontrol=list(is=F, stochastic=FALSE),
#                         nlcontrol=list(),
#                         nopriors=TRUE,
#                         chains=2,
#                         cores=1,
#                         #numOfThreads=numOfThreads,
#                         inits=NULL, forcerecompile=FALSE,
#                         savescores=FALSE, gendata=FALSE,
#                         control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
#                         verbose=0)
# }
# debug <- 0



#' ctmaEqual
#'
#' @param activeDirectory ?
#' @param sourceDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param ctmaMultipleFit ?
#' @param activateRPB ?
#' @param checkSingleStudyResults ?
#' @param silentOverwrite ?
#' @param digits ?
#' @param equalDrift ?
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param confidenceIntervals ?
#' @param saveEqualDriftModelFit ?
#' @param CoTiMAStanctArgs ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaEqual <- function(
  # Directory names and file names
  activeDirectory=NULL,
  sourceDirectory= NULL,
  resultsFilePrefix="CoTiMAEqualDrift",
  saveFilePrefix="CoTiMAEqualDrift",

  # Primary Study Fits
  ctmaMultipleFit=NULL,                    #list of lists: could be more than one fit object

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  checkSingleStudyResults=FALSE,
  silentOverwrite=FALSE,

  digits=4,

  # General Model Setup
  #compareDrift=c(),

  # Fitting Parameters
  equalDrift=NULL,
  retryattempts=30,
  refits=1,
  NPSOL=FALSE,
  coresToUse=c(-1),
  #fullDriftStartValues=NULL,
  #compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  #useCTMultiGroupAlt=FALSE,
  confidenceIntervals=FALSE,
  saveEqualDriftModelFit=saveFilePrefix,
  CoTiMAStanctArgs=list(test=TRUE, scaleTI=TRUE, scaleTime=1/1,
                        savesubjectmatrices=FALSE, verbose=1,
                        datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
                        iter=1000, intoverstates=TRUE,
                        binomial=FALSE, fit=TRUE,
                        intoverpop=FALSE, stationary=FALSE,
                        plot=FALSE, derrind="all",
                        optimize=TRUE, optimcontrol=list(is=F, stochastic=FALSE),
                        nlcontrol=list(),
                        nopriors=TRUE,
                        chains=2,
                        cores=1,
                        inits=NULL, forcerecompile=FALSE,
                        savescores=FALSE, gendata=FALSE,
                        control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
                        verbose=0)

)


{  # begin function definition (until end of file)

  # check if mutipleDirftFit object is supplied
  if (! ((ctmaMultipleFit$model.type == "mx") || (ctmaMultipleFit$model.type == "stanct")) ) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object with more than a single invariant drift effect (fit of CoTiMAMultipleDrift) has to be supplied compare the effects. \n"))
    stop("Good luck for the next try!")
  }

  # check if fit object is specified
  if (is.null(ctmaMultipleFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object with more than a single invariant drift effect (fit of CoTiMAMultipleDrift) has to be supplied compare the effects. \n"))
    stop("Good luck for the next try!")
  }

  if (resultsFilePrefix=="CoTiMAEqualDrift") {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat("The default results file prefix (CoTiMAEqualDrift) has been chosen.", "\n")
    cat("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards ", "\n")
    char <- readline(" ")
    while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
      cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
      char <- readline(" ")
    }
    if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
  }

  #if (saveFilePrefix=="CoTiMAEqualDrift") {
  #  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
  #  cat("The default save file prefix (CoTiMAEqualDrift) has been chosen.", "\n")
  #  cat("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards ")
  #  char <- readline(" ")
  #  while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
  #    cat((crayon::blue("Please press 'q' to quit and change filename or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
  #    char <- readline(" ")
  #  }
  #  if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
  #}




  #######################################################################################################################
  ############################################### attach further packages ###############################################
  #######################################################################################################################
  #{
  #  print(paste0("#################################################################################"))
  #  print(paste0("############################ Attach Further Packages ############################"))
  #  print(paste0("#################################################################################"))
  #
  #  # If OpenMx and ctsem are already attached, detaching them is required to enable OpenMx to run in parallel mode.
  #  if ("OpenMx" %in% (.packages())) suppressWarnings(detach("package:OpenMx", force=TRUE)) #, unload=TRUE)
  #  if ("ctsem" %in% (.packages())) suppressWarnings(detach("package:ctsem", force=TRUE)) #, unload=TRUE))
  #  Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
  #  library(ctsem)
  #  mxOption(key='Number of Threads', value=parallel::detectCores()) #now
  #
  #  # Attach all required packages that were not already attach with "library(ctsem)" before
  #  if (!("MASS" %in% (.packages()))) library(MASS)
  #  if (!("MBESS" %in% (.packages()))) library(MBESS)
  #  if (!("rootSolve" %in% (.packages()))) library(rootSolve)
  #  if (!("doParallel" %in% (.packages()))) library(doParallel)  # this is for parallel processing loops (not for internal parallel processing of OpenMx)
  #  if (!("crayon" %in% (.packages()))) library(crayon)
  #  if (!("psych" %in% (.packages()))) library(psych)
  #
  #  if (activateRPB==TRUE) {
  #    if("RPushbullet" %in% rownames(installed.packages()) == FALSE) {install.packages("RPushbullet")}
  #    if (!("RPushbullet" %in% (.packages()))) library(RPushbullet)
  #  }
  #  }

  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  {
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
    coresToUse
    CoTiMAStanctArgs$cores <- coresToUse

    if ((-1) * coresToUse >= parallel::detectCores()) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- parallel::detectCores() - 1
      cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
    }

    numOfThreads <- round((coresToUse/refits - .4), 0); numOfThreads
    if(numOfThreads < 1) numOfThreads <- 1
    tmp <- paste0("No of Threads set to ", numOfThreads); tmp
    cat("","\n")
    cat(tmp,"\n")
  }


  #######################################################################################################################
  ######## Extracting parameters from fitted primary studies created with CoTiMAMultipleDrift Function  #################
  ############### and estimating modified model in which invsariant (fixed) effcts are set equal. #######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  if (is.null(activeDirectory)) activeDirectory <- ctmaMultipleFit$activeDirectory; activeDirectory
  if (is.null(sourceDirectory)) sourceDirectory <- ctmaMultipleFit$sourceDirectory; sourceDirectory
  n.studies <- unlist(ctmaMultipleFit$n.studies); n.studies
  #drift_Coef <- ctmaMultipleFit$studyResults$DRIFT; drift_Coef
  allTpoints <- ctmaMultipleFit$statisticsList$allTpoints; allTpoints
  maxTpoints <- max(allTpoints); maxTpoints
  allDeltas <- ctmaMultipleFit$statisticsList$allDeltas; allDeltas
  maxDelta <- max(allDeltas); maxDelta
  #manifestNames <- c(ctmaMultipleFit$studyFitList[[1]]$ctmodelobj$manifestNames); manifestNames

  if (ctmaMultipleFit$model.type == "mx")  {
    results <- ctmaEqualMx(
      #type="multipleDrift",
      equalDrift=equalDrift,
      #datawide_all=datawide_all,
      #groups=groups,
      #groupsNamed=groupsNamed,
      activeDirectory=activeDirectory,
      sourceDirectory=sourceDirectory,
      resultsFilePrefix="CoTiMAEqualDrift",
      saveFilePrefix="CoTiMAEqualDrift",
      ctmaMultipleFit=ctmaMultipleFit,
      activateRPB=activateRPB,
      #checkSingleStudyResults=checkSingleStudyResults,
      digits=4,
      retryattempts=retryattempts,
      refits=refits,
      NPSOL=NPSOL,
      coresToUse=coresToUse,
      #fullDriftStartValues=fullDriftStartValues,
      #compSVMethod=compSVMethod,
      #useCTMultiGroupAlt=useCTMultiGroupAlt,
      confidenceIntervals=confidenceIntervals,
      saveEqualDriftModelFit=saveEqualDriftModelFit
    )
  }

  if (ctmaMultipleFit$model.type == "stanct") {
    results <- ctmaEqualStanct(
      #type="multipleDrift",
      equalDrift=equalDrift,
      ctmaMultipleFit=ctmaMultipleFit,
      #datalong_all=datalong_all,
      activeDirectory=activeDirectory,
      sourceDirectory=sourceDirectory,
      resultsFilePrefix="CoTiMAEqualDrift",
      saveFilePrefix="CoTiMAEqualDrift",
      activateRPB=activateRPB,
      silentOverwrite=silentOverwrite,
      digits=digits,
      #n.latent=n.latent,
      #n.studies=n.studies,
      retryattempts=retryattempts,
      refits=refits,
      NPSOL=NPSOL,
      coresToUse=coresToUse,
      #numOfThreads=numOfThreads,
      #fullDriftStartValues=fullDriftStartValues,
      saveEqualDriftModelFit=saveEqualDriftModelFit,
      CoTiMAStanctArgs=CoTiMAStanctArgs
    )
  }

  #class(results) <- "CoTiMAFit"
  invisible(results)

} ### END function definition
