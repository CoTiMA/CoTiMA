#######################################################################################################################
############################## CoTiMA Multiple Drift Effects Invariant Across primary Studies #########################
#######################################################################################################################

# takes a fitted CoTiMA as argumnt, and all parameters that were set invariant in these model are set equal

#' ctmaEqual
#'
#' @param ctmaInvariantFit ""
#' @param activeDirectory ""
#' @param activateRPB ""
#' @param digits ""
#' @param equalDrift ""
#' @param coresToUse ""
#' @param CoTiMAStanctArgs ""
#'
ctmaEqual <- function(
  # Primary Study Fits
  ctmaInvariantFit=NULL,                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,
  #sourceDirectory= NULL,
  #resultsFilePrefix="CoTiMAEqualDrift",
  #saveFilePrefix="CoTiMAEqualDrift",

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  #checkSingleStudyResults=FALSE,
  #silentOverwrite=FALSE,

  digits=4,

  # General Model Setup

  # Fitting Parameters
  equalDrift=NULL,
  coresToUse=1,
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

  # check if mutipleDriftFit object is supplied
  if (! ((ctmaInvariantFit$model.type == "mx") || (ctmaInvariantFit$model.type == "stanct")) ) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object with more than a single invariant drift effect (fit of CoTiMAMultipleDrift) has to be supplied compare the effects. \n"))
    stop("Good luck for the next try!")
  }

  # check if fit object is specified
  if (is.null(ctmaInvariantFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object with more than a single invariant drift effect (fit of CoTiMAMultipleDrift) has to be supplied compare the effects. \n"))
    stop("Good luck for the next try!")
  }

  if ( length(grep("invariant", names(ctmaInvariantFit$modelResults$DRIFT))) < 2) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object was supplied, but is has to have more than a single invariant drift effect to compare the effects. \n"))
    stop("Good luck for the next try!")
  }



  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  if  (length(coresToUse) > 0) {
    if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
  }

  if (coresToUse >= parallel::detectCores()) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    coresToUse <- parallel::detectCores() - 1
    cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
  }



  #######################################################################################################################
  ######## Extracting parameters from fitted primary studies created with CoTiMAMultipleDrift Function  #################
  ############### and estimating modified model in which invsariant (fixed) effcts are set equal. #######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time
{
  n.latent <- length(ctmaInvariantFit$modelResults$DRIFT)^.5; n.latent
  if (is.null(activeDirectory)) activeDirectory <- ctmaInvariantFit$activeDirectory; activeDirectory
  n.studies <- unlist(ctmaInvariantFit$n.studies); n.studies
  allTpoints <- ctmaInvariantFit$statisticsList$allTpoints; allTpoints
  maxTpoints <- max(allTpoints); maxTpoints
  allDeltas <- ctmaInvariantFit$statisticsList$allDeltas; allDeltas
  maxDelta <- max(allDeltas); maxDelta
  manifestNames <- ctmaInvariantFit$studyFitList$ctstanmodel$manifestNames; manifestNames
  parameterNames <- ctmaInvariantFit$parameterNames; parameterNames
  driftNames <- ctmaInvariantFit$parameterNames$DRIFT; driftNames
  #driftNames <- names(ctmaInvariantFit$modelResults$DRIFT); driftNames
  targetNames <- names(ctmaInvariantFit$modelResults$DRIFT[grep("invariant", names(ctmaInvariantFit$modelResults$DRIFT))]); targetNames
}

  # copy previous model
  prevStanctModel <- ctmaInvariantFit$studyFitList[[1]]$ctstanmodelbase
  prevStanctModelFit <- summary(ctmaInvariantFit$studyFitList[[1]])

  # identify Drift coefficents that were fixed (across all TI, which is just a check)
  tmpRow <- which(prevStanctModel$pars$matrix == "DRIFT"); tmpRow
  equalDriftPos <- grep("invariant", names(ctmaInvariantFit$modelResults$DRIFT)); equalDriftPos
  tmpRow <- tmpRow[equalDriftPos]; tmpRow
  tmp1 <- prevStanctModel$pars[tmpRow, paste0(prevStanctModel$TIpredNames,'_effect')]; tmp1
  tmp2 <- apply(tmp1, 1, unique); tmp2
  tmp3 <- unlist(lapply(tmp2, length)); tmp3
  targetDriftRow <- tmpRow[tmp3==1 & tmp2==FALSE]; targetDriftRow

  # new model
  stanctModel <- prevStanctModel
  newDriftLabel <- paste(stanctModel$pars[targetDriftRow, "param"], collapse = "_eq_"); newDriftLabel
  stanctModel$pars[targetDriftRow, "param"] <- newDriftLabel

  prevEst <- ctmaInvariantFit$studyFitList[[1]]$stanfit$rawest; prevEst
  # correct for having appropriate inits
  tmp1 <- which(tmpRow %in% targetDriftRow); tmp1
  newInits <- mean(prevEst[tmp1]); newInits
  prevEst[tmp1] <- newInits; prevEst
  prevEst <- prevEst[-tmp1[-1]]; prevEst

  prevData <- ctmaInvariantFit$data


  fitStanctModel <- ctsem::ctStanFit(
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

  if (!(is.null(CoTiMAStanctArgs$resample))) {
    fitStanctModel <- ctmaStanResample(fitStanctModel=fitStanctModel, CoTiMAStanctArgs=CoTiMAStanctArgs)
  }

  equalDriftStanctFit <- summary(fitStanctModel, digits=digits)


  Tvalues <- equalDriftStanctFit$parmatrices[,3]/equalDriftStanctFit$parmatrices[,4]; Tvalues
  equalDrift_Coeff <- round(cbind(equalDriftStanctFit$parmatrices, Tvalues), digits); equalDrift_Coeff
  # re-label
  tmp1 <- which(rownames(equalDrift_Coeff) == "DRIFT"); tmp1
  rownames(equalDrift_Coeff)[tmp1] <- "DRIFT " # space at the end for later extraction
  tmp2 <- tmp1[equalDriftPos]; tmp2
  tmp3 <- paste0(rownames(equalDrift_Coeff)[tmp2] , "(invariant & equal)"); tmp3
  #tmp1 <- which(rownames(equalDrift_Coeff) %in% newDriftLabel); tmp1
  #tmp2 <- paste0(rownames(equalDrift_Coeff)[tmp1] , " (invariant)"); tmp2
  #rownames(equalDrift_Coeff)[tmp1] <- tmp2; equalDrift_Coeff
  rownames(equalDrift_Coeff)[tmp2] <- tmp3; equalDrift_Coeff
  #
  equalDrift_Minus2LogLikelihood  <- -2*equalDriftStanctFit$loglik; equalDrift_Minus2LogLikelihood
  equalDrift_estimatedParameters  <- equalDriftStanctFit$npars; equalDrift_estimatedParameters
  #equalDrift_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
  #  equalDrift_estimatedParameters; equalDrift_df
  #n.par.first.lag <- ((2 * n.latent) * (2 * n.latent + 1)) / 2; n.par.first.lag
  #n.par.later.lag <- ((2 * n.latent) * (2 * n.latent - 1)) / 2; n.par.later.lag
  #n.later.lags <- allTpoints - n.latent; n.later.lags
  #equalDrift_df <- sum(n.later.lags * n.par.later.lag); equalDrift_df
  equalDrift_df <- NULL

  model_Drift_Coef <- equalDrift_Coeff[grep("DRIFT ", rownames(equalDrift_Coeff)), 3]; model_Drift_Coef
  names(model_Drift_Coef) <- driftNames
  names(model_Drift_Coef)[equalDriftPos] <- newDriftLabel; model_Drift_Coef

  model_Diffusion_Coef <- equalDrift_Coeff[(rownames(equalDrift_Coeff) == "DIFFUSIONcov"), 3]; model_Diffusion_Coef
  model_Diffusion_Coef <- c(OpenMx::vech2full(model_Diffusion_Coef)); model_Diffusion_Coef
  names(model_Diffusion_Coef) <- paste0("diff_", driftNames); model_Diffusion_Coef

  model_T0var_Coef <- equalDrift_Coeff[(rownames(equalDrift_Coeff) == "T0VAR"), 3]; model_T0var_Coef
  model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
  names(model_T0var_Coef) <- paste0("T0VAR_", driftNames); model_T0var_Coef

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  st <- paste0("Computation started at: ", start.time); st
  et <- paste0("Computation ended at: ", end.time); et
  tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt

  # results is built again after test for equality
  results <- list(activeDirectory=activeDirectory,
                  plot.type=NULL, model.type="stanct",
                  time=list(start.time=start.time, end.time=end.time, time.taken=time.taken),
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=ctmaInvariantFit$studyList, studyFitList=fitStanctModel,
                  data=NULL, statisticsList=ctmaInvariantFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=NULL),
                  parameterNames=ctmaInvariantFit$parameterNames,
                  CoTiMAStanctArgs=CoTiMAStanctArgs,
                  equalDrift=newDriftLabel,
                  summary=list(#model=paste(targetNames, collapse=" equal to "),
                               model=newDriftLabel,
                               estimates=round(equalDrift_Coeff, digits), #[]
                               minus2ll= round(equalDrift_Minus2LogLikelihood, digits),
                               n.parameters = round(equalDrift_estimatedParameters, digits),
                               df= NULL))
  class(results) <- "CoTiMAFit"

  # model comparison
  compFit <- ctmaCompFit(ctmaInvariantFit, results)[[1]][-7]
  results$summary$ll_diff_test <- compFit

  invisible(results)

} ### END function definition