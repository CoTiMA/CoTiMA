#' ctmaEqual
#'
#' @description test if the two or more invariant drift parameters in the CoTiMAFit object supplied are equal. The supplied CoTiMA fit-object
#' (ctmaInvariantFit) has to be a model fitted with \code{\link{ctmaFit}} where at least two parameters were set invariant across primary studies (e.g., 2 cross
#' effects). All parameters that are set invariant in the supplied model are then constrained to be equal by ctmaEqual (no user action
#' required), the model is fitted, and a log-liklihood ratio test is performed informing about the probability that equality applies.
#'
#' @param ctmaInvariantFit  object to which a CoTiMA fit has been assigned to (i.e., what has been returned by \code{\link{ctmaFit}}).
#' In most cases probably a model in which (only) two effects were specified with invariantDrift.
#' @param activeDirectory defines another active directory than the one used in ctmaInvariantFit
#' @param activateRPB  set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param digits Number of digits used for rounding (in outputs)
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#'
#' @importFrom  RPushbullet pbPost
#' @importFrom  parallel detectCores
#' @importFrom  ctsem ctStanFit
#' @importFrom  OpenMx vech2full
#'
#' @export ctmaEqual
#'
#' @examples
#' # Fit a CoTiMA with a set of parameters set equal that were set
#' # invariant in a previous model (of which the fit object is
#' # supplied in argument ctmaInvariantFit)
#' \dontrun{
#' CoTiMAFullInv23Fit_6$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAFullInvEq23Fit_6 <- ctmaEqual(ctmaInvariantFit=CoTiMAFullInv23Fit_6)
#' }
#'
#' @return returns a model where two or more parameters were set equal across primary studies and a log-likelihood difference test
#' informing about the probability that the equality assumption is correct.
#'
ctmaEqual <- function(
  ctmaInvariantFit=NULL,
  activeDirectory=NULL,
  activateRPB=FALSE,
  digits=4,
  coresToUse=2
  )

{  # begin function definition (until end of file)

  # use same fitting params as used to fit ctmaInvariantFit
  CoTiMAStanctArgs <- ctmaInvariantFit$CoTiMAStanctArgs

  # check if mutipleDriftFit object is supplied
  if (! ((ctmaInvariantFit$model.type == "mx") || (ctmaInvariantFit$model.type == "stanct")) ) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object with more than a single invariant drift effect (fit of ctmaFit) has to be supplied compare the effects. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  # check if fit object is specified
  if (is.null(ctmaInvariantFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object with more than a single invariant drift effect (fit of ctmaFit) has to be supplied compare the effects. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if ( length(grep("invariant", names(ctmaInvariantFit$modelResults$DRIFT))) < 2) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object was supplied, but is has to have more than a single invariant drift effect to compare the effects. \nGood luck for the next try!"
    stop(ErrorMsg)
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
    Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash."
    message(Msg)
  }



  #######################################################################################################################
  ############## Extracting parameters from fitted primary studies created with ctmaFit Function  #######################
  ############### and estimating modified model in which invariant (fixed) effects are set equal. #######################
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
  targetNames <- names(ctmaInvariantFit$modelResults$DRIFT[grep("invariant", names(ctmaInvariantFit$modelResults$DRIFT))]); targetNames
}

  # copy previous model
  prevStanctModel <- ctmaInvariantFit$studyFitList[[1]]$ctstanmodelbase
  if (is.null(prevStanctModel)) prevStanctModel <- ctmaInvariantFit$studyFitList$ctstanmodelbase
  prevStanctModelFit <- summary(ctmaInvariantFit$studyFitList[[1]])
  if (!("npars" %in% names(prevStanctModelFit))) prevStanctModelFit <- summary(ctmaInvariantFit$studyFitList)

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
  if (is.null(prevEst)) prevEst <- ctmaInvariantFit$studyFitList$stanfit$rawest
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
    #derrind=CoTiMAStanctArgs$derrind, # CHD deprecated, deleted Aug 2023
    optimize=CoTiMAStanctArgs$optimize,
    optimcontrol=CoTiMAStanctArgs$optimcontrol,
    nlcontrol=CoTiMAStanctArgs$nlcontrol,
    nopriors=CoTiMAStanctArgs$nopriors,
    priors=CoTiMAStanctArgs$priors, # added Aug 2023
    chains=CoTiMAStanctArgs$chains,
    forcerecompile=CoTiMAStanctArgs$forcerecompile,
    savescores=CoTiMAStanctArgs$savescores,
    gendata=CoTiMAStanctArgs$gendata,
    control=CoTiMAStanctArgs$control,
    verbose=CoTiMAStanctArgs$verbose,
    cores=coresToUse)

  if (!(is.null(CoTiMAStanctArgs$resample))) {
    fitStanctModel <- ctmaStanResample(ctmaFittedModel=fitStanctModel) #, CoTiMAStanctArgs=CoTiMAStanctArgs)
  }

  equalDriftStanctFit <- summary(fitStanctModel, digits=digits)

  # Extract estimates & statistics
  # account for changes in ctsem 3.4.1
  if ("matrix" %in% colnames(equalDriftStanctFit$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
  tmpMean <- grep("ean", colnames(equalDriftStanctFit$parmatrices)); tmpMean
  tmpSd <- tmpMean+1; tmpSd
  Tvalues <- equalDriftStanctFit$parmatrices[,tmpMean]/equalDriftStanctFit$parmatrices[,tmpSd]; Tvalues
  equalDrift_Coeff <- cbind(equalDriftStanctFit$parmatrices, Tvalues); equalDrift_Coeff
  equalDrift_Coeff[, tmpMean:(dim(equalDrift_Coeff)[2])] <- round(equalDrift_Coeff[, tmpMean:(dim(equalDrift_Coeff)[2])], digits); equalDrift_Coeff

  # re-label
  if (ctsem341) {
    tmp1 <- which(equalDrift_Coeff[, "matrix"] == "DRIFT")
    driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); driftNamesTmp
    rownames(equalDrift_Coeff) <- paste0(equalDrift_Coeff[, c("matrix")], "_",
                                         equalDrift_Coeff[, c("row")], "_",
                                         equalDrift_Coeff[, c("col")])
  } else {
    tmp1 <- which(rownames(equalDrift_Coeff) == "DRIFT")
    driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); driftNamesTmp
  }

  driftNamesTmp[equalDriftPos] <- paste0(driftNamesTmp[equalDriftPos], " (invariant & equal)"); driftNamesTmp
  rownames(equalDrift_Coeff)[tmp1] <- driftNamesTmp; equalDrift_Coeff
  tmp2 <- grep("toV", rownames(equalDrift_Coeff)); tmp2
  #tmp2 <- grep("invariant", rownames(equalDrift_Coeff)); tmp2
  #tmp3 <- paste0("DRIFT ", rownames(equalDrift_Coeff)[tmp2] , " (invariant & equal)"); tmp3
  #rownames(equalDrift_Coeff)[tmp2] <- tmp3; equalDrift_Coeff
  #tmp4 <- tmp1[which(!(tmp1 %in% tmp2))]; tmp4 # change to "DRIFT " for later extraction
  #rownames(invariantDrift_Coeff)[tmp4] <- paste0("DRIFT ", driftNames[which(!(tmp1 %in% tmp2))]); invariantDrift_Coeff
  rownames(equalDrift_Coeff)[tmp2] <- paste0("DRIFT ", rownames(equalDrift_Coeff)[tmp2]); equalDrift_Coeff


  #tmp1 <- which(rownames(equalDrift_Coeff) == "DRIFT"); tmp1
  #rownames(equalDrift_Coeff)[tmp1] <- "DRIFT " # space at the end for later extraction
  #tmp2 <- tmp1[equalDriftPos]; tmp2
  #tmp3 <- paste0(rownames(equalDrift_Coeff)[tmp2] , "(invariant & equal)"); tmp3
  ##tmp1 <- which(rownames(equalDrift_Coeff) %in% newDriftLabel); tmp1
  ##tmp2 <- paste0(rownames(equalDrift_Coeff)[tmp1] , " (invariant)"); tmp2
  ##rownames(equalDrift_Coeff)[tmp1] <- tmp2; equalDrift_Coeff
  #rownames(equalDrift_Coeff)[tmp2] <- tmp3; equalDrift_Coeff
  #
  equalDrift_Minus2LogLikelihood  <- -2*equalDriftStanctFit$loglik; equalDrift_Minus2LogLikelihood
  equalDrift_estimatedParameters  <- equalDriftStanctFit$npars; equalDrift_estimatedParameters
  equalDrift_df <- "deprecated"

  model_Drift_Coef <- equalDrift_Coeff[grep("DRIFT ", rownames(equalDrift_Coeff)), tmpMean]; model_Drift_Coef
  names(model_Drift_Coef) <- rownames(equalDrift_Coeff)[grep("DRIFT ", rownames(equalDrift_Coeff))]
  names(model_Drift_Coef)[equalDriftPos] <- newDriftLabel; model_Drift_Coef
  names(model_Drift_Coef) <- gsub("DRIFT ", "", names(model_Drift_Coef)); model_Drift_Coef

  model_Diffusion_Coef <- equalDrift_Coeff[grep("DIFFUSIONcov", substr(rownames(equalDrift_Coeff), 1, 12)), tmpMean] ; model_Diffusion_Coef
  if (!(ctsem341)) model_Diffusion_Coef <- c(OpenMx::vech2full(model_Diffusion_Coef)); model_Diffusion_Coef
  names(model_Diffusion_Coef) <- paste0("diff_", driftNames); model_Diffusion_Coef

  if (ctsem341) {
    model_T0var_Coef <- equalDrift_Coeff[grep("T0cov", substr(rownames(equalDrift_Coeff), 1, 5)), tmpMean]; model_T0var_Coef
  } else {
    model_T0var_Coef <- equalDrift_Coeff[grep("T0VAR", substr(rownames(equalDrift_Coeff), 1, 5)), tmpMean]; model_T0var_Coef
    model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
  }
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
                               estimates=equalDrift_Coeff,
                               minus2ll=equalDrift_Minus2LogLikelihood,
                               n.parameters = round(equalDrift_estimatedParameters, digits),
                               scaleTime=ctmaInvariantFit$summary$scaleTime))
  class(results) <- "CoTiMAFit"

  # model comparison
  compFit <- ctmaCompFit(ctmaInvariantFit, results)[[1]][-7]
  results$summary$ll_diff_test <- compFit

  invisible(results)

} ### END function definition
