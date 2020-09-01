#######################################################################################################################
############################################## CoTiMA InvariantDrift ##################################################
#######################################################################################################################

#' ctmaFull
#'
#' @param ctmaInitFit ""
#' @param primaryStudyList  ""
#' @param activeDirectory  ""
#' @param activateRPB  ""
#' @param digits  ""
#' @param type  ""
#' @param coresToUse  ""
#' @param CoTiMAStanctArgs  ""
#'
#' @return
#' @export
#'
ctmaFit <- function(
  # Primary Study Fits
  ctmaInitFit=NULL,                    #list of lists: could be more than one fit object
  primaryStudyList=NULL,               # created by the PREP file for moderator analyses
  cluster=NULL,                        # vecto with clust variables (e.g., countries)

  # Directory names and file names
  activeDirectory=NULL,

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,

  digits=4,

  # General Model Setup
  invariantDrift=NULL,

  # Fitting Parameters
  type="stanct",
  coresToUse=c(-1),
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
                        verbose=0,
                        warmup=500)

)


{  # begin function definition (until end of file)

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
    stop("Good luck for the next try!")
  }

  #######################################################################################################################
  ####### Copy/Change INIT File based on information delivered by different PREP files (e.g., moderator studies ) #######
  #######################################################################################################################

  if (!(is.null(primaryStudyList))) {
    ctmaTempFit <- ctmaInitFit
    targetStudyNumbers <- unlist(primaryStudyList$studyNumbers); targetStudyNumbers; length(targetStudyNumbers)
    for (i in (length(ctmaTempFit$studyFitList)):1) {
      if (!(ctmaTempFit$studyList[[i]]$originalStudyNo %in% targetStudyNumbers)) {
        ctmaTempFit$studyList[[i]] <- NULL
        ctmaTempFit$studyFitList[[i]] <- NULL
        ctmaTempFit$emprawList[[i]] <- NULL
        ctmaTempFit$statisticsList$originalStudyNumbers[i] <- NA
        ctmaTempFit$statisticsList$allSampleSizes[i+1] <- NA
        ctmaTempFit$statisticsList$allTpoints[i] <- NA
        ctmaTempFit$modelResults[[1]][[i]] <- NULL
        ctmaTempFit$modelResults[[2]][[i]] <- NULL
        ctmaTempFit$modelResults[[3]][[i]] <- NULL
      }
    }
    ctmaTempFit$n.studies <- length(targetStudyNumbers); ctmaTempFit$n.studies
    ctmaTempFit$statisticsList$allDeltas <- unlist(lapply(ctmaTempFit$studyList, function(extract) extract$delta_t))
    ctmaTempFit$statisticsList$allDeltas
    ctmaTempFit$statisticsList$minDelta <- min(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$minDelta
    ctmaTempFit$statisticsList$maxDelta <- max(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$maxDelta
    ctmaTempFit$statisticsList$meanDelta <- mean(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$meanDelta
    ctmaTempFit$statisticsList$overallSampleSize <- sum(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$overallSampleSize
    ctmaTempFit$statisticsList$meanSampleSize <- mean(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$meanSampleSize
    ctmaTempFit$statisticsList$maxSampleSize <- max(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$maxSampleSize
    ctmaTempFit$statisticsList$minSampleSize <- min(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$minSampleSize
    ctmaTempFit$statisticsList$overallTpoints <- sum(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$overallTpoints
    ctmaTempFit$statisticsList$meanTpoints <- mean(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$meanTpoints
    ctmaTempFit$statisticsList$maxTpoints <- max(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$maxTpoints
    ctmaTempFit$statisticsList$minTpoints <- min(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$minTpoints
    ctmaTempFit$summary$model <- "Moderator Model (for details see model summary)"
    tmpStudyNumber <- as.numeric(gsub("Study No ", "", rownames(ctmaTempFit$summary$estimates))); tmpStudyNumber
    targetRows <- which(tmpStudyNumber %in% targetStudyNumbers); targetRows; length(targetRows)
    ctmaTempFit$summary$estimates <- ctmaTempFit$summary$estimates[targetRows, ]
    ctmaTempFit$summary$estimates
    ctmaTempFit$summary$confidenceIntervals <- ctmaTempFit$summary$confidenceIntervals[targetRows, ]
    ctmaTempFit$summary$confidenceIntervals
    ctmaTempFit$summary$n.parameters <- ctmaTempFit$studyFitList[[1]]$resultsSummary$npars * length(targetRows)
    ctmaTempFit$summary$n.parameters
    ctmaTempFit$statisticsList$originalStudyNumbers <-
      ctmaTempFit$statisticsList$originalStudyNumbers[which(!(is.na(ctmaTempFit$statisticsList$originalStudyNumbers)))]
    ctmaTempFit$statisticsList$originalStudyNumbers
    ctmaTempFit$statisticsList$allSampleSizes <-
      ctmaTempFit$statisticsList$allSampleSizes[which(!(is.na(ctmaTempFit$statisticsList$allSampleSizes)))]
    ctmaTempFit$statisticsList$allSampleSizes
    ctmaTempFit$statisticsList$allTpoints <-
      ctmaTempFit$statisticsList$allTpoints[which(!(is.na(ctmaTempFit$statisticsList$allTpoints)))]
    ctmaTempFit$statisticsList$allTpoints

    ctmaInitFit <- ctmaTempFit
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
  ############# Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {
    n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    driftNames <- ctmaInitFit$parameterNames$DRIFT; driftNames
    if (is.null(invariantDrift)) invariantDrift <- driftNames
    usedTimeRange <- seq(0, 1.5*maxDelta, 1)

  }


  #######################################################################################################################
  ################################################# data preparation ####################################################
  #######################################################################################################################
  {
    # combine pseudo raw data for mx model
    tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit)
    datawide_all <- tmp$alldata
    groups <- tmp$groups
    names(groups) <- c("Study_No_"); groups
    groupsNamed <- (paste0("Study_No_", groups)); groupsNamed

    # auggment pseudo raw data for stanct model
    dataTmp <- cbind(datawide_all, groups)
    for (i in 1:(n.studies-1)) {
      tmp <- matrix(0, nrow=nrow(dataTmp)); tmp
      colnames(tmp) <- paste0("TI", i); tmp
      dataTmp <- cbind(dataTmp, tmp); dim(dataTmp)
      tmp <- which(dataTmp[,"groups"] == i); tmp
      dataTmp[tmp, ncol(dataTmp)] <- 1
      if (CoTiMAStanctArgs$scaleTI == TRUE) dataTmp[ , ncol(dataTmp)] <- scale(dataTmp[ , ncol(dataTmp)])
    }
    targetCols <- which(colnames(dataTmp) == "groups"); targetCols
    dataTmp <- dataTmp[ ,-targetCols]

    # add clusters as dummy moderators
    if (!(is.null(cluster))) {
      # determine number of required dummies
      targetCluster <- which(table(cluster) > 1); targetCluster  # no cluster if only one study is included
      targetCluster <- names(targetCluster); targetCluster
      clusCounter <- length(targetCluster); clusCounter
      # create dummies
      tmpTI <- matrix(0, dim(dataTmp)[1], clusCounter)
      for (i in 1:clusCounter) {
        targetGroups <- which(cluster == targetCluster[i]); targetGroups
        tmp2 <- which(groups %in% targetGroups); length(tmp2)
        tmpTI[tmp2, i] <- 1
      }
      if (CoTiMAStanctArgs$scaleClus == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)])
      currentStartNumber <- n.studies; currentStartNumber
      currentEndNumber <- currentStartNumber + clusCounter -1; currentEndNumber
      colnames(tmpTI) <- paste0("TI", currentStartNumber:currentEndNumber); colnames(tmpTI)
      dataTmp <- cbind(dataTmp, tmpTI); dim(dataTmp)
    } else {
      clusCounter <- 0
    }


    dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.latent, n.TIpred = (n.studies-1+clusCounter),
                                    manifestNames=manifestNames)
    dataTmp3 <- ctsem::ctDeintervalise(dataTmp2)
    dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
    # eliminate rows where ALL latents are NA
    dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
    datalong_all <- dataTmp3
  }


  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################


  # Make model with most time points
  {
    stanctModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                                  DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent),
                                  LAMBDA=diag(n.latent),
                                  CINT=matrix(0, nrow=n.latent, ncol=1),
                                  T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                  MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                  MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent),
                                  MANIFESTTRAITVAR = 'auto',
                                  type = 'stanct',
                                  n.TIpred = (n.studies-1+clusCounter),
                                  TIpredNames = paste0("TI", 1:(n.studies-1+clusCounter)),
                                  TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+clusCounter)))


    #stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'LAMBDA',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'CINT',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTMEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTVAR',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    if (!(is.null(cluster))) {
      stanctModel$pars[stanctModel$pars$matrix %in% 'DIFFUSION',paste0(stanctModel$TIpredNames[n.studies:(n.studies+clusCounter-1)],'_effect')] <- FALSE
      stanctModel$pars[stanctModel$pars$matrix %in% 'T0VAR',paste0(stanctModel$TIpredNames[n.studies:(n.studies+clusCounter-1)],'_effect')] <- FALSE
    }
  }

  # the target effects
  tmp1 <- which(stanctModel$pars$matrix == "DRIFT"); tmp1
  tmp2 <- which(stanctModel$pars[tmp1, "param"] %in% invariantDrift); tmp2
  stanctModel$pars[tmp1[tmp2], paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
  stanctModel$pars


  fitStanctModel <- ctsem::ctStanFit(
    datalong = datalong_all,
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
    warmup=CoTiMAStanctArgs$warmup,
    cores=coresToUse)

  ### resample in parcels to avoid memory crash and speed up
  if (!(is.null(CoTiMAStanctArgs$resample))) {
    fitStanctModel <- ctmaStanResample(fitStanctModel=fitStanctModel, CoTiMAStanctArgs=CoTiMAStanctArgs)
  }

  invariantDriftStanctFit <- summary(fitStanctModel, digits=2*digits, parmatrices=TRUE, residualcov=FALSE)

  # Extract estimates & statistics
  Tvalues <- invariantDriftStanctFit$parmatrices[,3]/invariantDriftStanctFit$parmatrices[,4]; Tvalues
  invariantDrift_Coeff <- round(cbind(invariantDriftStanctFit$parmatrices, Tvalues), digits); invariantDrift_Coeff
  # re-label
  tmp1 <- which(rownames(invariantDrift_Coeff) == "DRIFT"); tmp1
  rownames(invariantDrift_Coeff)[tmp1] <- driftNames; invariantDrift_Coeff
  tmp2 <- which(rownames(invariantDrift_Coeff) %in% invariantDrift); tmp2
  tmp3 <- paste0("DRIFT ", rownames(invariantDrift_Coeff)[tmp2] , " (invariant)"); tmp3
  rownames(invariantDrift_Coeff)[tmp2] <- tmp3; invariantDrift_Coeff
  tmp4 <- tmp1[which(!(tmp1 %in% tmp2))]; tmp4 # change to "DRIFT " for later extraction
  rownames(invariantDrift_Coeff)[tmp4] <- paste0("DRIFT ", driftNames[which(!(tmp1 %in% tmp2))]); invariantDrift_Coeff

  invariantDrift_Minus2LogLikelihood  <- -2*invariantDriftStanctFit$loglik; invariantDrift_Minus2LogLikelihood
  invariantDrift_estimatedParameters  <- invariantDriftStanctFit$npars; invariantDrift_estimatedParameters
  #n.par.first.lag <- ((2 * n.latent) * (2 * n.latent + 1)) / 2; n.par.first.lag
  #n.par.later.lag <- ((2 * n.latent) * (2 * n.latent - 1)) / 2; n.par.later.lag
  #n.later.lags <- allTpoints - n.latent; n.later.lags
  #invariantDrift_df <- sum(n.later.lags * n.par.later.lag); invariantDrift_df
  #invariantDrift_df <- invariantDrift_df + (n.studies-1) * n.latent^2; invariantDrift_df
  invariantDrift_df <- NULL

  model_Drift_Coef <- invariantDrift_Coeff[(grep("DRIFT ", rownames(invariantDrift_Coeff))), 3]; model_Drift_Coef

  model_Diffusion_Coef <- invariantDrift_Coeff[(rownames(invariantDrift_Coeff) == "DIFFUSIONcov"), 3]; model_Diffusion_Coef
  model_Diffusion_Coef <- c(OpenMx::vech2full(model_Diffusion_Coef)); model_Diffusion_Coef
  names(model_Diffusion_Coef) <- driftNames; model_Diffusion_Coef

  model_T0var_Coef <- invariantDrift_Coeff[(rownames(invariantDrift_Coeff) == "T0VAR"), 3]; model_T0var_Coef
  model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
  names(model_T0var_Coef) <- driftNames; model_Diffusion_Coef

  ## cluster effects
  if (!(is.null(cluster))) {
    tmp1 <- c()
    for (i in (n.studies):(n.studies+clusCounter-1)) tmp1 <- c(tmp1, (grep(i, rownames(invariantDriftStanctFit$tipreds))))
    Tvalues <- invariantDriftStanctFit$tipreds[tmp1, ][,6]; Tvalues
    clusTI_Coeff <- round(cbind(invariantDriftStanctFit$tipreds[tmp1, ], Tvalues), digits); clusTI_Coeff
    # re-label
    for (i in 1:clusCounter) {
      targetNamePart <- paste0("tip_TI", n.studies+i-1); targetNamePart
      newNamePart <- paste0(targetCluster[i], "_on_"); newNamePart
      rownames(clusTI_Coeff) <- sub(targetNamePart, paste0(targetCluster[i], "_on_"), rownames(clusTI_Coeff))
    }
  } else {
    clusTI_Coeff <- NULL
  }


  ### Numerically compute Optimal Time lag sensu Dormann & Griffin (2015)
  driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=T); driftMatrix # byrow set because order is different compared to mx model
  #driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=F); driftMatrix # byrow set because order is different compared to mx model
  OTL <- function(timeRange) {
    OpenMx::expm(driftMatrix * timeRange)[targetRow, targetCol]}
  # loop through all cross effects
  optimalCrossLag <- matrix(NA, n.latent, n.latent)
  maxCrossEffect <- matrix(NA, n.latent, n.latent)
  for (j in 1:n.latent) {
    for (h in 1:n.latent) {
      if (j != h) {
        targetRow <- j
        targetCol <- h
        targetParameters <- sapply(usedTimeRange, OTL)
        maxCrossEffect[j,h] <- max(abs(targetParameters))
        optimalCrossLag[j,h] <- which(abs(targetParameters)==maxCrossEffect[j,h])*1+0
      }
    }
  }
  #} ## END  fit stanct model

  #######################################################################################################################

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  st <- paste0("Computation started at: ", start.time); st
  et <- paste0("Computation ended at: ", end.time); et
  tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt



  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}


  tmp1 <- grep("CINT", rownames(invariantDriftStanctFit$parmatrices)); tmp1
  tmp2 <- grep("asym", rownames(invariantDriftStanctFit$parmatrices)); tmp2
  tmp3 <- grep("dt", rownames(invariantDriftStanctFit$parmatrices)); tmp3
  tmp4 <- tmp1[(tmp1 %in% c(tmp2, tmp3)) == FALSE]; tmp4
  model_Cint_Coef <- invariantDriftStanctFit$parmatrices[tmp4, 3]; model_Cint_Coef

  results <- list(activeDirectory=activeDirectory,
                  time=list(start.time=start.time, end.time=end.time, time.taken=time.taken),
                  plot.type="drift",  model.type="stanct",
                  coresToUse=coresToUse, n.studies=1,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=list(fitStanctModel),
                  data=datalong_all, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=ctmaInitFit$parameterNames,
                  CoTiMAStanctArgs=CoTiMAStanctArgs,
                  invariantDrift=invariantDrift,
                  summary=list(model="invariant drift (hom. model)",
                               estimates=round(invariantDrift_Coeff, digits),
                               minus2ll= round(invariantDrift_Minus2LogLikelihood, digits),
                               n.parameters = round(invariantDrift_estimatedParameters, digits),
                               df= invariantDrift_df,
                               opt.lag = optimalCrossLag,
                               max.effects = round(maxCrossEffect, digits),
                               clus.effects=clusTI_Coeff))
  class(results) <- "CoTiMAFit"

  invisible(results)

} ### END function definition
