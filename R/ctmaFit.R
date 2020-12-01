#######################################################################################################################
############################################## CoTiMA InvariantDrift ##################################################
#######################################################################################################################

#' ctmaFit
#'
#' @param ctmaInitFit object to which all single ctsem fits of primary studies has been assigned to (i.e., what has been returned by ctmaInit)
#' @param primaryStudyList  could be a list of primary studies compiled with ctmaPrep that defines the subset of studies in ctmaInitFit that should actually be used
#' @param cluster  vector with cluster variables (e.g., countries). Has to be set up carfully. Will be included in ctmaPrep later.
#' @param activeDirectory  defines another active directory than the one used in ctmaInitFit
#' @param activateRPB  set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param digits Number of digits used for rounding (in outputs)
#' @param invariantDrift  drift Labels for drift effects that are set invariant across primary studies (default = all drift effects).
#' @param drift Labels for drift effects. Have to be either of the type V1toV2 or 0 for effects to be excluded, which is usually not recommended)
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#' @param indVarying Allows ct intercepts to vary at the individual level (random effects model, accounts for unobserved heteregeneity)
#' @param scaleTI scale TI predictors - not recommended if TI are dummies representing primary studies as probably in most instances
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param nopriors if TRUE, any priors are disabled – sometimes desirable for optimization
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#'
#' @importFrom  RPushbullet pbPost
#' @importFrom  crayon red
#' @importFrom  parallel detectCores
#' @importFrom  ctsem ctWideToLong ctDeintervalise ctModel ctStanFit
#' @importFrom  OpenMx vech2full expm
#'
#' @export ctmaFit
#'
#' @examples Example 1. Fit a CoTiMA to all primary studies previously fitted one by one with the fits assigned to CoTiMAInitFit_Ex1
#' CoTiMAFullFit_Ex1 <- ctmaFit(ctmaInitFit=CoTiMAInitFit_Ex1)
#'
#' saveRDS(CoTiMAFullFit_Ex1, file=paste0(activeDirectory, "CoTiMAFullFit_Ex1.rds"))
#' summary(CoTiMAFullFit_Ex1)
#'
#' @examples Example 2. Fit a CoTiMA with only 2 cross effects invariant (not the auto effects) to all primary studies previously fitted one by one with the fits assigned to CoTiMAInitFit_Ex1
#' CoTiMA12lFit_Ex2 <- ctmaFit(ctmaInitFit=CoTiMAInitFit_Ex1,
#' invariantDrift=c("V1toV2", "V2toV1"))
#'
#' saveRDS(CoTiMA12Fit_Ex2, file=paste0(activeDirectory, "CoTiMA12Fit_Ex2.rds"))
#' summary(CoTiMA12lFit_Ex2)
#'
ctmaFit <- function(
  ctmaInitFit=NULL,
  primaryStudyList=NULL,
  cluster=NULL,
  activeDirectory=NULL,
  activateRPB=FALSE,
  digits=4,
  invariantDrift=NULL,
  drift=NULL,
  indVarying=FALSE,
  coresToUse=c(1),
  scaleTI=NULL,
  scaleTime=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  finishsamples=NULL,
  chains=NULL,
  verbose=NULL
)


{  # begin function definition (until end of file)

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
    stop("Good luck for the next try!")
  }

  { # fitting params
    if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
    if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
    if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
    if (!(is.null(nopriors))) CoTiMAStanctArgs$nopriors <- nopriors
    if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
    if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
    if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
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
    ctmaTempFit$statisticsList$minDelta <- min(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$maxDelta <- max(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$meanDelta <- mean(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$overallSampleSize <- sum(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$meanSampleSize <- mean(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$maxSampleSize <- max(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$minSampleSize <- min(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$overallTpoints <- sum(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$meanTpoints <- mean(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$maxTpoints <- max(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$minTpoints <- min(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$summary$model <- "Moderator Model (for details see model summary)"
    tmpStudyNumber <- as.numeric(gsub("Study No ", "", rownames(ctmaTempFit$summary$estimates))); tmpStudyNumber
    targetRows <- which(tmpStudyNumber %in% targetStudyNumbers); targetRows; length(targetRows)
    ctmaTempFit$summary$estimates <- ctmaTempFit$summary$estimates[targetRows, ]
    ctmaTempFit$summary$confidenceIntervals <- ctmaTempFit$summary$confidenceIntervals[targetRows, ]
    ctmaTempFit$summary$n.parameters <- ctmaTempFit$studyFitList[[1]]$resultsSummary$npars * length(targetRows)
    ctmaTempFit$statisticsList$originalStudyNumbers <-
      ctmaTempFit$statisticsList$originalStudyNumbers[which(!(is.na(ctmaTempFit$statisticsList$originalStudyNumbers)))]
    ctmaTempFit$statisticsList$allSampleSizes <-
      ctmaTempFit$statisticsList$allSampleSizes[which(!(is.na(ctmaTempFit$statisticsList$allSampleSizes)))]
    ctmaTempFit$statisticsList$allTpoints <-
      ctmaTempFit$statisticsList$allTpoints[which(!(is.na(ctmaTempFit$statisticsList$allTpoints)))]

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
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    if (is.null(manifestNames)) n.manifest <- 0 else n.manifest <- length(manifestNames)
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

    # augment pseudo raw data for stanct model
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
    n.var <- max(n.manifest, n.latent); n.var
    #dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.latent, n.TIpred = (n.studies-1+clusCounter),
    #                                manifestNames=manifestNames)
    dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.var, n.TIpred = (n.studies-1+clusCounter),
                                    manifestNames=manifestNames)
    dataTmp3 <- ctsem::ctDeintervalise(dataTmp2)
    dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
    utils::head(dataTmp3)
    # eliminate rows where ALL latents are NA
    if (n.manifest > n.latent) namePart <- "y" else namePart <- "V"
    #dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
    dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0(namePart, 1:n.var)], 1, function(x) sum(is.na(x)) != n.var ), ]
    datalong_all <- dataTmp3
  }
  datalong_all <- as.data.frame(datalong_all)

  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################

  tmp1 <- ctmaInitFit$studyFitList[[1]]$ctstanmodelbase$pars

  tmp2 <- tmp1$matrix=="DRIFT"; tmp2
  DRIFT <- tmp1[tmp2,]$param; DRIFT
  DRIFT[is.na(DRIFT)] <- "0"; DRIFT

  tmp2 <- tmp1$matrix=="LAMBDA"; tmp2
  LAMBDA <- tmp1[tmp2,]$value; LAMBDA
  LAMBDA <- matrix(LAMBDA, nrow= n.var, byrow=TRUE); LAMBDA

  tmp2 <- tmp1$matrix=="MANIFESTVAR"; tmp2
  MANIFESTVAR <- tmp1[tmp2,]$value; MANIFESTVAR
  MANIFESTVAR <- matrix(MANIFESTVAR, nrow= n.var); MANIFESTVAR

  # Make model with max time points
  {
    stanctModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=maxTpoints, manifestNames=manifestNames,
                                  DRIFT=matrix(DRIFT, nrow=n.latent, ncol=n.latent, byrow=TRUE),
                                  LAMBDA=LAMBDA,
                                  CINT=matrix(0, nrow=n.latent, ncol=1),
                                  T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                  MANIFESTMEANS = matrix(c(0), nrow = n.var, ncol = 1),
                                  MANIFESTVAR=MANIFESTVAR,
                                  type = 'stanct',
                                  n.TIpred = (n.studies-1+clusCounter),
                                  TIpredNames = paste0("TI", 1:(n.studies-1+clusCounter)),
                                  TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+clusCounter)))

    if (indVarying == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
      print(paste0("#################################################################################"))

      manifestNames <- paste0("y", 1:n.manifest); manifestNames
      latentNames <- paste0("V", 1:2); latentNames
      if (n.manifest == n.latent) manifestNames <- latentNames
      MANIFESTMEANS <- manifestNames; MANIFESTMEANS

      stanctModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=maxTpoints, manifestNames=manifestNames,
                                    DRIFT=matrix(DRIFT, nrow=n.latent, ncol=n.latent, byrow=TRUE),
                                    LAMBDA=LAMBDA,
                                    CINT=matrix(0, nrow=n.latent, ncol=1),
                                    T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                    MANIFESTMEANS = matrix(MANIFESTMEANS, nrow = n.manifest, ncol = 1),
                                    MANIFESTVAR=MANIFESTVAR,
                                    type = 'stanct',
                                    n.TIpred = (n.studies-1+clusCounter),
                                    TIpredNames = paste0("TI", 1:(n.studies-1+clusCounter)),
                                    TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+clusCounter)))
    }

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
  #stanctModel$pars[tmp1[tmp2], paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
  stanctModel$pars[tmp1[tmp2], paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- FALSE

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
    #control=CoTiMAStanctArgs$control,
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
  #driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); driftNamesTmp
  driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); driftNamesTmp
  #rownames(invariantDrift_Coeff)[tmp1] <- driftNames; invariantDrift_Coeff
  rownames(invariantDrift_Coeff)[tmp1] <- driftNamesTmp; invariantDrift_Coeff
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
  names(model_T0var_Coef) <- driftNames; model_T0var_Coef

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
  driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=FALSE); driftMatrix # byrow set because order is different compared to mx model
  #driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=TRUE); driftMatrix # byrow set because order is different compared to mx model
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
        if (driftMatrix[j, h] != 0) { # an effect that is zero has no optimal lag
          targetParameters <- sapply(usedTimeRange, OTL)
          maxCrossEffect[j,h] <- max(abs(targetParameters))
          optimalCrossLag[j,h] <- which(abs(targetParameters)==maxCrossEffect[j,h])*1+0
        } else {
          optimalCrossLag[j,h] <- NA
        }
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
                  summary=list(model=paste(invariantDrift, "unequal but invariant across samples", collapse=" "),
                               estimates=round(invariantDrift_Coeff, digits),
                               randomEffects=invariantDriftStanctFit$popsd,
                               minus2ll= round(invariantDrift_Minus2LogLikelihood, digits),
                               n.parameters = round(invariantDrift_estimatedParameters, digits),
                               df= invariantDrift_df,
                               opt.lag = optimalCrossLag,
                               max.effects = round(maxCrossEffect, digits),
                               clus.effects=clusTI_Coeff))
  class(results) <- "CoTiMAFit"

  invisible(results)

} ### END function definition
