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
#   ctmaMultipleFit=CoTiMAFit
#   activateRPB=FALSE
#   #checkSingleStudyResults=TRUE
#   retryattempts=30
#   refits=1
#   NPSOL=FALSE
#   coresToUse=7
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

#' ctmaEqualStanct
#'
#' @param ctmaMultipleFit ?
#' @param equalDrift ?
#' @param activeDirectory ?
#' @param sourceDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param activateRPB ?
#' @param checkSingleStudyResults ?
#' @param silentOverwrite ?
#' @param digits ?
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
ctmaEqualStanct <- function(
  # Primary Study Fits
  ctmaMultipleFit=NULL,                    #list of lists: could be more than one fit object
  equalDrift=NULL,

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
  manifestNames <- ctmaMultipleFit$studyFitList$ctstanmodel$manifestNames; manifestNames
  parameterNames <- ctmaMultipleFit$parameterNames; parameterNames
  driftNames <- parameterNames$DRIFT; driftNames


  # copy previous model
  prevStanctModel <- ctmaMultipleFit$studyFitList$ctstanmodelbase

  # identify Drift coefficents that were fixed (across all TI, which is just a check)
  tmpRow <- which(prevStanctModel$pars$matrix == "DRIFT"); tmpRow
  tmp1 <- prevStanctModel$pars[tmpRow, paste0(prevStanctModel$TIpredNames,'_effect')]; tmp1
  tmp2 <- apply(tmp1, 1, unique); tmp2
  tmp3 <- unlist(lapply(tmp2, length)); tmp3
  targetDriftRow <- tmpRow[tmp3==1 & tmp2==FALSE]; targetDriftRow

  # new model
  stanctModel <- prevStanctModel
  newDriftLabel <- paste(stanctModel$pars[targetDriftRow, "param"], collapse = "_eq_"); newDriftLabel
  stanctModel$pars[targetDriftRow, "param"] <- newDriftLabel

  prevEst <- ctmaMultipleFit$studyFitList$stanfit$rawest; prevEst
  # correct for having appropriate inits
  tmp1 <- which(tmpRow %in% targetDriftRow); tmp1
  newInits <- mean(prevEst[tmp1]); newInits
  prevEst[tmp1] <- newInits; prevEst
  prevEst <- prevEst[-tmp1[-1]]; prevEst

  prevData <- ctmaMultipleFit$data

  fitStanctModel2 <- ctsem::ctStanFit(
    inits=prevEst,
    datalong = prevData,
    ctstanmodel = stanctModel,
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
    nopriors=CoTiMAStanctArgs$nopriors,
    chains=CoTiMAStanctArgs$chains,
    forcerecompile=CoTiMAStanctArgs$forcerecompile,
    savescores=CoTiMAStanctArgs$savescores,
    gendata=CoTiMAStanctArgs$gendata,
    control=CoTiMAStanctArgs$control,
    verbose=CoTiMAStanctArgs$verbose,
    cores=coresToUse)

  equalDriftStanctFit <- summary(fitStanctModel2, digits=digits)



  # SAVE
  if (length(saveEqualDriftModelFit) > 0) {
    x1 <- paste0(saveEqualDriftModelFit[1], " STANCT ", (Sys.time()), ".rds"); x1
    x2 <- c()
    ctmaSaveFile(activateRPB, activeDirectory, equalDriftStanctFit, x1, x2, silentOverwrite=silentOverwrite)
  }

  #for (i in 1:length(targetNames)){
  Tvalues <- equalDriftStanctFit$popmeans[,1]/equalDriftStanctFit$popmeans[,2]; Tvalues
  equalDrift_Coeff <- round(cbind(equalDriftStanctFit$popmeans, Tvalues), digits); equalDrift_Coeff
  # re-label
  tmp1 <- which(rownames(equalDrift_Coeff) %in% newDriftLabel); tmp1
  tmp2 <- paste0(rownames(equalDrift_Coeff)[tmp1] , " (fixed)"); tmp2
  rownames(equalDrift_Coeff)[tmp1] <- tmp2; equalDrift_Coeff
  #
  equalDrift_Minus2LogLikelihood  <- -2*equalDriftStanctFit$logprob; equalDrift_Minus2LogLikelihood
  equalDrift_estimatedParameters  <- equalDriftStanctFit$npars; equalDrift_estimatedParameters
  #equalDrift_estimatedParameters  <- length(fitStanctModel$stanfit$optimfit$par); equalDrift_estimatedParameters
  equalDrift_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
    equalDrift_estimatedParameters; equalDrift_df
  model_Drift_Coef <- equalDrift_Coeff[grep("toV", rownames(equalDrift_Coeff)), 1]; model_Drift_Coef
  #names(model_Drift_Coef) <- driftNames; model_Drift_Coef
  names(model_Drift_Coef) <- rownames(equalDrift_Coeff)[grep("toV", rownames(equalDrift_Coeff))]; model_Drift_Coef
  model_Diffusion_Coef <- equalDrift_Coeff[grep("diffusion", rownames(equalDrift_Coeff)), 1]; model_Diffusion_Coef
  names(model_Diffusion_Coef) <- rownames(equalDrift_Coeff)[grep("diffusion", rownames(equalDrift_Coeff))]; model_Diffusion_Coef
  model_T0var_Coef <- equalDrift_Coeff[grep("T0var", rownames(equalDrift_Coeff)), 1]; model_T0var_Coef
  names(model_T0var_Coef) <- rownames(equalDrift_Coeff)[grep("T0var", rownames(equalDrift_Coeff))]; model_T0var_Coef

  results <- list(activeDirectory=activeDirectory, sourceDirectory=sourceDirectory,
                  plot.type=NULL, model.type="stanct",
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=ctmaMultipleFit$studyList, studyFitList=fitStanctModel2, # , homDRIFTEqualFitCI),
                  data=NULL, statisticsList=ctmaMultipleFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=NULL),
                  parameterNames=ctmaMultipleFit$parameterNames,
                  summary=list(model=paste(targetNames, "fixed", collapse=" equal to "),
                               estimates=round(equalDrift_effects, digits), #[]
                               minus2ll= round(equalDriftMinus2LogLikelihood, digits),
                               n.parameters = round(equalDriftestimatedParameters, digits),
                               df= c(round(equalDriftdf, digits))))
  class(results) <- "CoTiMAFit"

  tmp <- paste(targetNames, collapse="_eq_"); tmp
  saveRDS(results, paste0(activeDirectory, "CoTiMAEqualDriftSTANCT ", tmp, " allResults", " ", Sys.time(), ".rds"))

  invisible(results)



  #class(results) <- "CoTiMAFit"
  invisible(results)

} ### END function definition
