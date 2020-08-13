#######################################################################################################################
############################## CoTiMA Multiple Drift Effects Invariant Across primary Studies #########################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   #type="mx"
#   activeDirectory=NULL
#   equalDrift=NULL
#   resultsFilePrefix="CoTiMAEqualDrift"
#   saveFilePrefix="CoTiMAEqualDrift"
#   ctmaMultipleFit=CoTiMAMultipleMxFit
#   activateRPB=FALSE
#   #checkSingleStudyResults=TRUE
#   retryattempts=30
#   refits=1
#   NPSOL=FALSE
#   coresToUse=6
#   fullDriftStartValues=NULL
#   compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2]
#   useCTMultiGroupAlt=FALSE
#   confidenceIntervals=FALSE
#   saveEqualDriftModelFit=saveFilePrefix
#   ##### ENTER DEBUG INFO HERE #######
#   silentOverwrite=FALSE
#   #ctmaInitFit=ctmaInitFit1
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

#' ctmaEqualMx
#'
#' @param ctmaMultipleFit ?
#' @param activeDirectory ?
#' @param sourceDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
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
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaEqualMx <- function(
  # Primary Study Fits
  ctmaMultipleFit=NULL,                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,
  sourceDirectory= NULL,
  resultsFilePrefix="CoTiMAEqualDrift",
  saveFilePrefix="CoTiMAEqualDrift",


  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  checkSingleStudyResults=FALSE,
  silentOverwrite=FALSE,

  digits=4,

  # General Model Setup
  #compareDrift=c(),

  # Fitting Parameters
  #type="mx",
  equalDrift=NULL,
  retryattempts=30,
  refits=1,
  NPSOL=FALSE,
  coresToUse=c(-1),
  #fullDriftStartValues=NULL,
  #compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  #useCTMultiGroupAlt=FALSE,
  confidenceIntervals=FALSE,
  saveEqualDriftModelFit=saveFilePrefix
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

  #if (resultsFilePrefix=="CoTiMAEqualDrift") {
  #  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
  #  cat("The default results file prefix (CoTiMAEqualDrift) has been chosen.", "\n")
  #  cat("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards ", "\n")
  #  char <- readline(" ")
  #  while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
  #    cat((blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
  #    char <- readline(" ")
  #  }
  #  if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
  #}


  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  #{
  #  if  (length(coresToUse) > 0) {
  #    if (coresToUse < 1)  {
  #      if (.Platform$OS.type == "unix") {
  #        coresToUse <- parallel::detectCores() + coresToUse
  #        #coresToUse <- detectCores() + coresToUse
  #        cat("     #################################################################################\n")
  #        cat("     # You decided using multiple CPU cores for parallel fitting CoTiMA models. This #\n")
  #        cat("     # may work well in many cases, but we highly recommend not to run this script   #\n")
  #        cat("     # in RStudio. Instead, go to the console (e.g., console.app on MAC OS), then    #\n")
  #        cat("     # change to the direcrtory where your R-script is located (by typing, e.g.,     #\n")
  #        cat("     # \'cd \"/Users/trump/only/temp\"\'), and then, e.g., \'Rscript \"CoTiMA1.R\"\'.#\n")
  #        cat("     #################################################################################")
  #        Sys.sleep(1)
  #      } else {
  #        coresToUse <- 1
  #      }
  #    }
  #  } else {
  #    coresToUse <- 1
  #  }
  #  coresToUse
  #  CoTiMAStanctArgs$cores <- coresToUse
  #
  #  if ((-1) * coresToUse >= parallel::detectCores()) {
  #    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
  #    coresToUse <- parallel::detectCores() - 1
  #    cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
  #  }
  #
  #  numOfThreads <- round((coresToUse/refits - .4), 0); numOfThreads
  #  if(numOfThreads < 1) numOfThreads <- 1
  #  tmp <- paste0("No of Threads set to ", numOfThreads); tmp
  #  cat("","\n")
  #  cat(tmp,"\n")
  #}


  #######################################################################################################################
  ######## Extracting parameters from fitted primary studies created with CoTiMAMultipleDrift Function  #################
  ############### and estimating modified model in which invsariant (fixed) effcts are set equal. #######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  extraTries <- retryattempts
  if (is.null(activeDirectory)) activeDirectory <- ctmaMultipleFit$activeDirectory; activeDirectory
  if (is.null(sourceDirectory)) sourceDirectory <- ctmaMultipleFit$sourceDirectory; sourceDirectory
  n.studies <- unlist(ctmaMultipleFit$n.studies); n.studies
  #drift_Coef <- ctmaInitFit$studyResults$DRIFT; drift_Coef
  allTpoints <- ctmaMultipleFit$statisticsList$allTpoints; allTpoints
  maxTpoints <- max(allTpoints); maxTpoints
  allDeltas <- ctmaMultipleFit$statisticsList$allDeltas; allDeltas
  maxDelta <- max(allDeltas); maxDelta
  #manifestNames <- c(ctmaMultipleFit$studyFitList[[1]]$ctmodelobj$manifestNames); manifestNames
  parameterNames <- ctmaMultipleFit$parameterNames; parameterNames
  driftNames <- parameterNames$DRIFT; driftNames

  prevMxModel <- ctmaMultipleFit$studyFitList
  if (prevMxModel$name == "ctsem multigroup") {
    n.studies <- length(prevMxModel$submodels); n.studies
    # make new drift labels
    driftLabels <- c(prevMxModel$Study_No_1$DRIFT$labels); driftLabels
    tmp1 <- grep("Study", driftLabels); tmp1
    targetLabels <- driftLabels[-tmp1]; targetLabels
    targetNames <- targetLabels
    newLabels <- paste(targetLabels, collapse = "_eq_"); newLabels
    driftLabels[-tmp1] <- newLabels; driftLabels
    # change submodels externally
    currentSubmodel <- list()
    for (i in 1:n.studies) {
      currentSubmodel[[i]] <- prevMxModel$submodels[[i]]; currentSubmodel[[i]]
      currentSubmodel[[i]]$DRIFT$labels[-tmp1] <- newLabels
      currentSubmodel[[i]]$DRIFT$values[-tmp1] <- mean(currentSubmodel[[i]]$DRIFT$values[-tmp1])
    }
    # combine changed submodels in new mx model
    equalDriftModel <- prevMxModel
    for (i in 1:n.studies) equalDriftModel <- OpenMx::mxModel(model=equalDriftModel, OpenMx::mxModel(currentSubmodel[[i]]))
    #equalDriftModel
    OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
    results <- parallel::mclapply(seq(1, refits, by=1), function(refits) OpenMx::mxTryHard(model=equalDriftModel, extraTries = extraTries),
                        mc.cores=coresToUse)
    OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
    # Select model with best fit
    allMinus2LogLikelihood <- lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
    equalDriftFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
  }

  if (prevMxModel$name == "ctsem") {
    n.studies <- length(grep("T0VAR_", names(prevMxModel))); n.studies
    targetMatrixNames <- paste0("DRIFT_", 1:n.studies); targetMatrixNames
    allMatrices <- prevMxModel$matrices
    targetMatrices <- allMatrices[names(allMatrices) %in% targetMatrixNames]
    # make new drift labels
    driftLabels <- targetMatrices[[1]]$labels; driftLabels
    tmp1 <- grep("_G", driftLabels); tmp1
    targetLabels <- driftLabels[-tmp1]; targetLabels
    newLabels <- paste(targetLabels, collapse = "_eq_"); newLabels
    driftLabels[-tmp1] <- newLabels; driftLabels
    # change submodels externally
    for (i in 1:n.studies) {
      #i <- 1
      targetMatrices[[i]]$labels[-tmp1] <- newLabels
      targetMatrices[[i]]$values[-tmp1] <- mean(targetMatrices[[i]]$values[-tmp1])
    }
    # combine changed submodels in new mx model
    equalDriftModel <- prevMxModel
    for (i in 1:n.studies) equalDriftModel <- OpenMx::mxModel(model=equalDriftModel, targetMatrices[[i]])
    OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
    results <- parallel::mclapply(seq(1, refits, by=1), function(refits) OpenMx::mxTryHard(model=equalDriftModel, extraTries = extraTries),
                        mc.cores=coresToUse)
    OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
    # Select model with best fit
    allMinus2LogLikelihood <- lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
    equalDriftFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
  }

  # SAVE
  if (length(saveEqualDriftModelFit) > 0) {
    x1 <- paste0(saveEqualDriftModelFit[1], " MX ", (Sys.time()), ".rds"); x1
    x2 <- c()
    ctmaSaveFile(activateRPB, activeDirectory, equalDriftFit, x1, x2, silentOverwrite=silentOverwrite)
  }

  # extract results
  tmp1 <- equalDriftFit$output$estimate; tmp1
  tmp2 <- grep("_eq_", names(tmp1)); tmp2
  equalDrift_Coef <- tmp1[tmp2]; equalDrift_Coef
  tmp1 <- equalDriftFit$output$standardErrors; tmp1
  tmp2 <- grep("_eq_", rownames(tmp1)); tmp2
  equalDrift_SE <- tmp1[tmp2]; equalDrift_SE
  Tvalue <- (equalDrift_Coef/equalDrift_SE); Tvalue
  equalDriftMinus2LogLikelihood  <- equalDriftFit$output$Minus2LogLikelihood; equalDriftMinus2LogLikelihood
  equalDriftestimatedParameters  <- length(equalDriftFit$output$estimate); equalDriftestimatedParameters
  equalDriftdf <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
    equalDriftestimatedParameters; equalDriftdf
  equalDrift_effects <- c(equalDrift_Coef, equalDrift_SE, Tvalue); equalDrift_effects
  names(equalDrift_effects) <- c(names(equalDrift_Coef), "SE", "Tvalues"); equalDrift_effects

  # Compute confidence intervals
  if (confidenceIntervals == TRUE ) {
    print(paste0("#################################################################################"))
    print(paste0("############### Computing Confidence Intervals for CoTiMA Model  ################"))
    print(paste0("#################################################################################"))

    equalDriftCI <- equalDriftFit
    ci <- OpenMx::mxCI(names(equalDrift_Coef)) # make mxCI object for every drift coefficients
    equalDriftCIModel <- OpenMx::mxModel(equalDriftFit, ci) # make OpenMx Models for all drift coefficients
    OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
    equalDriftCIFit <- OpenMx::mxRun(equalDriftCIModel, intervals=TRUE)
    equalDriftCI <- equalDriftCIFit$output$confidenceIntervals; equalDriftCI

    equalDrift_effects <- c(equalDrift_effects, equalDriftCI[c(1,3)]); equalDrift_effects
    names(equalDrift_effects) <- c(names(equalDrift_Coef), "SE", "Tvalues", "lbound", "ubound"); equalDrift_effects
  }

  # extract results
  model_Drift_Coef <- model_Diffusion_Coef <- model_T0var_Coef <- list()
  if (prevMxModel$name == "ctsem") {
    tmp1 <- equalDriftFit$output$estimate; tmp1
    tmp2 <- grep("toV", names(tmp1)); tmp2
    all_model_Drift_Coef <- tmp1[tmp2]; all_model_Drift_Coef
    tmp1 <- grep("_G", names(all_model_Drift_Coef)); tmp1
    group_model_Drift_Coef <- all_model_Drift_Coef[tmp1]; group_model_Drift_Coef
    fixed_model_Drift_Coef <- all_model_Drift_Coef[!(all_model_Drift_Coef %in% group_model_Drift_Coef)]; fixed_model_Drift_Coef
    for (i in 1:n.studies) {
      tmp1 <- grep(paste0("_G", i), names(group_model_Drift_Coef)); tmp1
      current_group_model_Drift_Coef <- group_model_Drift_Coef[tmp1]; current_group_model_Drift_Coef
      tmp2 <- rep(fixed_model_Drift_Coef, length(targetNames)); tmp2
      names(tmp2) <- targetNames; tmp2
      current_group_model_Drift_Coef <- c(current_group_model_Drift_Coef, tmp2); current_group_model_Drift_Coef
      current_group_model_Drift_Coef <- current_group_model_Drift_Coef[sort(names(current_group_model_Drift_Coef))]; current_group_model_Drift_Coef
      tmp3 <- which(names(current_group_model_Drift_Coef) %in% targetNames); tmp3
      names(current_group_model_Drift_Coef)[tmp3] <- rep(names(fixed_model_Drift_Coef), length(targetNames)); current_group_model_Drift_Coef
      current_group_model_Drift_Coef <- current_group_model_Drift_Coef[sort(names(current_group_model_Drift_Coef))]; current_group_model_Drift_Coef
      model_Drift_Coef[[i]] <- matrix(current_group_model_Drift_Coef, n.latent, n.latent); model_Drift_Coef[[i]]

      tmp1 <- grep(paste0("_G", i), names(equalDriftFit$output$estimate)); tmp1
      tmp2 <- equalDriftFit$output$estimate[tmp1]; tmp2
      tmp3 <- grep("diffusion", names(tmp2)); tmp3
      model_Diffusion_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_Diffusion_Coef[[i]]
      #model_Diffusion_Coef[[i]] <- base2Full(model_Diffusion_Coef[[i]]); model_Diffusion_Coef[[i]]

      tmp3 <- grep("T0var", names(tmp2)); tmp3
      model_T0var_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_T0var_Coef[[i]]
      #model_T0var_Coef[[i]] <- base2Full(model_T0var_Coef[[i]]); model_T0var_Coef[[i]]
    }
  } else {
    tmp1 <- equalDriftFit$output$estimate; tmp1
    tmp2 <- grep("toV", names(tmp1)); tmp2
    all_model_Drift_Coef <- tmp1[tmp2]; all_model_Drift_Coef
    tmp1 <- grep("Study_No", names(all_model_Drift_Coef)); tmp1
    group_model_Drift_Coef <- all_model_Drift_Coef[tmp1]; group_model_Drift_Coef
    fixed_model_Drift_Coef <- all_model_Drift_Coef[!(all_model_Drift_Coef %in% group_model_Drift_Coef)]; fixed_model_Drift_Coef
    for (i in 1:n.studies) {
      #i<- 1
      # to get the right order, relabeling and replacements are required
      tmp1 <- grep(paste0("Study_No_", i), names(group_model_Drift_Coef)); tmp1
      current_group_model_Drift_Coef <- group_model_Drift_Coef[tmp1]; current_group_model_Drift_Coef
      names(current_group_model_Drift_Coef) <- gsub(paste0("Study_No_", i, "_"), "", names(current_group_model_Drift_Coef)); current_group_model_Drift_Coef
      tmp2 <- rep(fixed_model_Drift_Coef, length(targetNames)); tmp2
      names(tmp2) <- targetNames; tmp2
      current_group_model_Drift_Coef <- c(current_group_model_Drift_Coef, tmp2); current_group_model_Drift_Coef
      current_group_model_Drift_Coef <- current_group_model_Drift_Coef[sort(names(current_group_model_Drift_Coef))]; current_group_model_Drift_Coef
      tmp3 <- which(names(current_group_model_Drift_Coef) %in% targetNames); tmp3
      names(current_group_model_Drift_Coef)[tmp3] <- rep(names(fixed_model_Drift_Coef), length(targetNames)); current_group_model_Drift_Coef
      current_group_model_Drift_Coef <- current_group_model_Drift_Coef[sort(names(current_group_model_Drift_Coef))]; current_group_model_Drift_Coef
      model_Drift_Coef[[i]] <- matrix(current_group_model_Drift_Coef, n.latent, n.latent); model_Drift_Coef[[i]]

      tmp1 <- grep(paste0("Study_No_", i), names(equalDriftFit$output$estimate)); tmp1
      tmp2 <- equalDriftFit$output$estimate[tmp1]; tmp2
      tmp3 <- grep("diffusion", names(tmp2)); tmp3
      model_Diffusion_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_Diffusion_Coef[[i]]
      #model_Diffusion_Coef[[i]] <- base2Full(model_Diffusion_Coef[[i]]); model_Diffusion_Coef[[i]]

      tmp3 <- grep("T0var", names(tmp2)); tmp3
      model_T0var_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_T0var_Coef[[i]]
      #model_T0var_Coef[[i]] <- base2Full(model_T0var_Coef[[i]]); model_T0var_Coef[[i]]
    }
  }

  results <- list(activeDirectory=activeDirectory, sourceDirectory=sourceDirectory,
                  plot.type=NULL, model.type="mx",
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=equalDriftFit, # , homDRIFTMultipleFitCI),
                  data=NULL, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=NULL),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model=paste(targetNames, "fixed", collapse=" equal to "),
                               estimates=round(equalDrift_effects, digits), #[]
                               minus2ll= round(equalDriftMinus2LogLikelihood, digits),
                               n.parameters = round(equalDriftestimatedParameters, digits),
                               df= c(round(equalDriftdf, digits))))
  class(results) <- "CoTiMAFit"

  tmp <- paste(targetNames, collapse="_&_"); tmp
  saveRDS(results, paste0(activeDirectory, "CoTiMAEqualDriftMX ", tmp, " allResults", " ", Sys.time(), ".rds"))

  invisible(results)



  #class(results) <- "CoTiMAFit"
  invisible(results)

} ### END function definition
