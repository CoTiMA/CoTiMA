#' ctmaLCS
#'
#' @description Transforms estimates obtained with \code{\link{ctmaFit}} into LCS (latent change score) terminology.
#' LCS models can be estimated with CT CLPM, but results have to be transformed. When time intervals vary much
#' between and within persons, LCS models are virtually impossible to fit. However, CT CLPM models can be
#' fitted, and the results - after transformation - show what LCS estimates would have been (cf Voelke & Oud,
#' 2015; their terminology to label LCS effects is used in the output created by ctmaLCS)
#'
#' @param CoTiMAFiT Fitted CoTiMA object.
#' @param undoTimeScaling Whether (TRUE) or not (FALSE) LCS results should be provided ignoring the scaleTime argument used in ctmaFit.
#' @param digits Number of digits used for rounding (in outputs)
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#'
#' @importFrom OpenMx expm
#' @importFrom RPushbullet pbPost
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#' LCSresults <- ctmaLCS(CoTiMAFullFit_6)
#' }
#'
#' @export ctmaLCS
#'
#' @return Returns LCS effects derived from CT CoTiMA CLPM estimates.
#'
#'
ctmaLCS <- function(CoTiMAFiT=NULL, undoTimeScaling=TRUE, digits=4, activateRPB=FALSE) {
  if (is.null("CoTiMAFiT$studyFitList$stanfit")) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object has to be supplied to compute LCS results \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  e <- ctExtract(CoTiMAFiT$studyFitList)
  if (undoTimeScaling == TRUE) {
    e$pop_DRIFT <- e$pop_DRIFT * CoTiMAFiT$argumentList$scaleTime
    #e$pop_DIFFUSION <- e$pop_DIFFUSION * CoTiMAFiT$argumentList$scaleTime # done later after GG'
  }

  # if init fit or fit before 6.6.2023 is used as input
  if (is.null(CoTiMAFiT$summary$randomEffects$model_popcov_025)) oldFit <- TRUE else oldFit <- FALSE

  driftNames <- CoTiMAFiT$parameterNames$DRIFT; driftNames
  diffNames <- CoTiMAFiT$parameterNames$DIFFUSION
  T0varNames <- CoTiMAFiT$parameterNames$T0VAR

  n.latent <- CoTiMAFiT$studyFitList$ctstanmodelbase$n.latent; n.latent
  n.manifest <- CoTiMAFiT$studyFitList$ctstanmodelbase$n.manifest; n.manifest

  driftNamesMatrix <- matrix(driftNames, n.latent, n.latent, byrow=T)

  tmp1 <- colnames(CoTiMAFiT$studyFitList$ctdatastruct); tmp1
  tmp1 <- tmp1[-(c(1:2))]; tmp1 # eliminate id & time name
  latentNames <- tmp1[(1:n.latent)]; latentNames

  # get results
  est <- CoTiMAFiT$summary$estimates; est
  colNames <- colnames(est); colNames
  rowNames <- rownames(est); rowNames
  targetColEst <- grep("ean", colNames); targetColEst
  targetColSD <- grep("sd", colNames); targetColSD
  targetColLL <- grep("2.5%", colNames); targetColLL
  targetColUL <- grep("97.5%", colNames); targetColUL
  #
  targetRow <- grep("T0MEAN", rowNames); targetRow
  initialMeans <- est[targetRow, targetColEst]; initialMeans
  initialMeansSD <- est[targetRow, targetColSD]; initialMeansSD
  initialMeansLL <- est[targetRow, targetColLL]; initialMeansLL
  initialMeansUL <- est[targetRow, targetColUL]; initialMeansUL
  names(initialMeans) <- paste0("InitialMean_", latentNames); initialMeans
  #
  targetRow <- grep("MANIFESTMEANS", rowNames); targetRow
  manifestMeans <- est[targetRow, targetColEst]; manifestMeans
  manifestMeansSD <- est[targetRow, targetColSD]; manifestMeansSD
  manifestMeansLL <- est[targetRow, targetColLL]; manifestMeansLL
  manifestMeansUL <- est[targetRow, targetColUL]; manifestMeansUL
  names(manifestMeans) <- paste0("ManifestMean_", latentNames); manifestMeans
  #
  tmp1 <- grep("CINT", rowNames); tmp1
  tmp2 <- grep("asym", rowNames); tmp2
  targetRow <- tmp1[!(tmp1 %in% tmp2)]; targetRow
  cints <- est[targetRow, targetColEst]; cints
  cintsSD <- est[targetRow, targetColSD]; cintsSD
  cintsLL <- est[targetRow, targetColLL]; cintsLL
  cintsUL <- est[targetRow, targetColUL]; cintsUL
  #
  tmp1 <- which(CoTiMAFiT$studyFitList$ctstanmodelbase$pars$matrix == "MANIFESTMEANS"); tmp1
  tmp2 <- CoTiMAFiT$studyFitList$ctstanmodelbase$pars$param[tmp1]; tmp2
  tmp1 <- which(CoTiMAFiT$studyFitList$ctstanmodelbase$pars$matrix == "CINT"); tmp1
  tmp3 <- CoTiMAFiT$studyFitList$ctstanmodelbase$pars$param[tmp1]; tmp3
  if ( (n.latent == n.manifest) & (!(any(is.na(manifestMeans)))) & (all(is.na(cints))) ) {
    slopeMeans_Cint <- manifestMeans
    slopeMeans_CintSD <- manifestMeansSD
    slopeMeans_CintLL <- manifestMeansLL
    slopeMeans_CintUL <- manifestMeansUL
  } else {
    slopeMeans_Cint <- cints
    slopeMeans_CintSD <- cintsSD
    slopeMeans_CintLL <- cintsLL
    slopeMeans_CintUL <- cintsUL
  }
  names(slopeMeans_Cint) <- paste0("SlopeMean_Cint", latentNames); slopeMeans_Cint
  #
  if (undoTimeScaling == TRUE) {
    drift <- matrix(CoTiMAFiT$modelResults$DRIFToriginal_time_scale, n.latent, n.latent, byrow=TRUE); drift
  } else {
    drift <- matrix(CoTiMAFiT$modelResults$DRIFT, n.latent, n.latent, byrow=TRUE); drift
  }
  ## get SDs
  tmp1 <- grep("DRIFT", rowNames); tmp1
  tmp2 <- grep("dt", rowNames); tmp2
  driftSD <- est[tmp1[!(tmp1 %in% tmp2)], targetColSD]; driftSD
  driftLL <- est[tmp1[!(tmp1 %in% tmp2)], targetColLL]; driftLL
  driftUL <- est[tmp1[!(tmp1 %in% tmp2)], targetColUL]; driftUL
  # rescale
  if (undoTimeScaling == TRUE) {
    driftSD <- matrix(driftSD * (CoTiMAFiT$argumentList$scaleTime), n.latent, n.latent); driftSD
    driftLL <- matrix(driftLL * (CoTiMAFiT$argumentList$scaleTime), n.latent, n.latent); driftLL
    driftUL <- matrix(driftUL * (CoTiMAFiT$argumentList$scaleTime), n.latent, n.latent); driftUL
  }  else {
    driftSD <- matrix(driftSD, n.latent, n.latent); driftSD
    driftLL <- matrix(driftLL, n.latent, n.latent); driftLL
    driftUL <- matrix(driftUL, n.latent, n.latent); driftUL
  }
  #
  autos <- diag(drift); autos
  autosSD <- diag(driftSD); autosSD
  autosLL <- diag(driftLL); autosLL
  autosUL <- diag(driftUL); autosUL
  names(autos) <- diag(driftNamesMatrix); autos
  autoNames <- names(autos) ; autoNames
  #
  crosss <- drift[!(drift %in% diag(drift))]; crosss
  crosssSD <- driftSD[!(driftSD %in% diag(driftSD))]; crosssSD
  crosssLL <- driftLL[!(driftLL %in% diag(driftLL))]; crosssLL
  crosssUL <- driftUL[!(driftUL %in% diag(driftUL))]; crosssUL
  names(crosss) <- sort(driftNames[!(driftNames %in% names(autos))]); crosss
  crossNames <- names(crosss) ; crossNames
  names(crosss) <- paste0("CTcross_", names(crosss) ); crosss
  names(autos) <- paste0("CTauto_", names(autos)); autos
  #
  expm_drift_minus1_samples <- expm_drift_samples <- drift_samples <- list()
  proportions <- autoregressions <- proportionsSD <- autoregressionsSD <- c()
  proportionsLL <- autoregressionsLL <- proportionsUL <- autoregressionsUL <- c()
  for (j in 1:dim(e$pop_DRIFT)[1]) {
    drift_samples[[j]] <- e$pop_DRIFT[j,,]
    expm_drift_samples[[j]] <- expm(e$pop_DRIFT[j,,])
    expm_drift_minus1_samples[[j]] <- expm_drift_samples[[j]]  - diag(n.latent)
  }
  for (i in 1:n.latent) {
    proportions[i] <- mean(unlist(lapply(expm_drift_minus1_samples, function(x) x[i,i])))
    proportionsSD[i] <- sd(unlist(lapply(expm_drift_minus1_samples, function(x) x[i,i])))
    proportionsLL[i] <- stats::quantile(unlist(lapply(expm_drift_minus1_samples, function(x) x[i,i])), prob=.025)
    proportionsUL[i] <- stats::quantile(unlist(lapply(expm_drift_minus1_samples, function(x) x[i,i])), prob=.975)
    autoregressions[i] <- mean(unlist(lapply(expm_drift_samples, function(x) x[i,i])))
    autoregressionsSD[i] <- sd(unlist(lapply(expm_drift_samples, function(x) x[i,i])))
    autoregressionsLL[i] <- stats::quantile(unlist(lapply(expm_drift_samples, function(x) x[i,i])), prob=.025)
    autoregressionsUL[i] <- stats::quantile(unlist(lapply(expm_drift_samples, function(x) x[i,i])), prob=.975)
  }
  names(proportions) <- paste0("Proportion_", latentNames); proportions
  names(autoregressions) <- paste0("Autoregressions_", latentNames, " (Delta=1)"); autoregressions
  #
  couplings <- crosslaggeds <- couplingsSD <- crosslaggedsSD <- c()
  couplingsLL <- crosslaggedsLL <- couplingsUL <- crosslaggedsUL <- c()
  counter <- 0
  for (i in 1:n.latent) {
    for (j in 1:n.latent) {
      if (i != j) {
        counter <- counter + 1
        couplings[counter] <- mean(unlist(lapply(drift_samples, function(x) x[j,i])))
        couplingsSD[counter] <- sd(unlist(lapply(drift_samples, function(x) x[j,i])))
        couplingsLL[counter] <- stats::quantile(unlist(lapply(drift_samples, function(x) x[j,i])), prob=.025)
        couplingsUL[counter] <- stats::quantile(unlist(lapply(drift_samples, function(x) x[j,i])), prob=.975)
        crosslaggeds[counter] <- mean(unlist(lapply(expm_drift_samples, function(x) x[j,i])))
        crosslaggedsSD[counter] <- sd(unlist(lapply(expm_drift_samples, function(x) x[j,i])))
        crosslaggedsLL[counter] <- stats::quantile(unlist(lapply(expm_drift_samples, function(x) x[j,i])), prob=.025)
        crosslaggedsUL[counter] <- stats::quantile(unlist(lapply(expm_drift_samples, function(x) x[j,i])), prob=.975)
      }
    }
  }
  names(couplings) <- paste0("Coupling_", crossNames); couplings
  names(crosslaggeds) <- paste0("CrossLagged_", crossNames, " (Delta=1)"); crosslaggeds
  #
  ## this block is more complex and involves Kroneker products and matrix inverse
  diffusion_samples <- list()
  for (j in 1:dim(e$pop_DIFFUSIONcov)[1]) {
    if (undoTimeScaling == TRUE) {
      diffusion_samples[[j]] <- e$pop_DIFFUSIONcov[j,,] * CoTiMAFiT$CoTiMAStanctArgs$scaleTime # GG'
    } else {
      diffusion_samples[[j]] <- e$pop_DIFFUSIONcov[j,,]
    }
  }
  DynErrVarsMatrix_samples <- list()
  for (j in 1:dim(e$pop_DRIFT)[1]) {
    kronek <- e$pop_DRIFT[j,,] %x% diag(n.latent) + diag(n.latent) %x% e$pop_DRIFT[j,,]
    DynErrVarsMatrix_samples[[j]] <- matrix(solve(kronek) %*% (OpenMx::expm(kronek) - diag(nrow(kronek))) %*% c(diffusion_samples[[j]]), n.latent, n.latent)
  }
  #
  DynErrVars <- DynErrVarsSD <- diffusionVars <- diffusionVarsSD <- c()
  DynErrVarsLL <- DynErrVarsUL <- diffusionVarsLL <- diffusionVarsUL <- c()
  for (i in 1:n.latent) {
    DynErrVars[i] <- mean(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[i,i])))
    DynErrVarsSD[i] <- sd(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[i,i])))
    DynErrVarsLL[i] <- stats::quantile(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[i,i])), prob=.025)
    DynErrVarsUL[i] <- stats::quantile(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[i,i])), prob=.975)
    diffusionVars[i] <- mean(unlist(lapply(diffusion_samples, function(x) x[i,i])))
    diffusionVarsSD[i] <- sd(unlist(lapply(diffusion_samples, function(x) x[i,i])))
    diffusionVarsLL[i] <- stats::quantile(unlist(lapply(diffusion_samples, function(x) x[i,i])), prob=.025)
    diffusionVarsUL[i] <- stats::quantile(unlist(lapply(diffusion_samples, function(x) x[i,i])), prob=.975)
  }
  names(diffusionVars) <- paste0("Diffusion_", latentNames); diffusionVars
  names(DynErrVars) <- paste0("DynErrVar_", latentNames, " (Delta=1)"); DynErrVars

  #
  DynErrCovs <- DynErrCovsSD <- diffusionCovs <- diffusionCovsSD <- c()
  DynErrCovsLL <- DynErrCovsUL <- diffusionCovsLL <- diffusionCovsUL <- c()
  counter <- 0
  for (i in 1:n.latent) {
    for (j in 1:n.latent) {
      if (i != j) {
        counter <- counter + 1
        DynErrCovs[counter] <- mean(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[j,i])))
        DynErrCovsSD[counter] <- sd(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[j,i])))
        DynErrCovsLL[counter] <- stats::quantile(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[j,i])), prob=.025)
        DynErrCovsUL[counter] <- stats::quantile(unlist(lapply(DynErrVarsMatrix_samples, function(x) x[j,i])), prob=.975)
        diffusionCovs[counter] <- mean(unlist(lapply(diffusion_samples, function(x) x[j,i])))
        diffusionCovsSD[counter] <- sd(unlist(lapply(diffusion_samples, function(x) x[j,i])))
        diffusionCovsLL[counter] <- stats::quantile(unlist(lapply(diffusion_samples, function(x) x[j,i])), prob=.025)
        diffusionCovsUL[counter] <- stats::quantile(unlist(lapply(diffusion_samples, function(x) x[j,i])), prob=.975)
      }
    }
  }
  names(diffusionCovs) <- names(couplings); diffusionCovs
  names(diffusionCovs) <- gsub("Coupling_", "Diffusion_", names(diffusionCovs)); diffusionCovs
  names(diffusionCovs) <- gsub("to", "with", names(diffusionCovs)); diffusionCovs
  names(DynErrCovs) <- names(diffusionCovs); DynErrVars
  names(DynErrCovs) <- gsub("Diffusion", "DynErrCov", names(DynErrCovs)); DynErrCovs
  names(DynErrCovs) <- paste0(names(DynErrCovs), " (Delta=1)" ); DynErrCovs
  #
  targetRow <- grep("MANIFESTcov", rowNames); targetRow
  measurementErrorVarsMatrix <- matrix(est[targetRow, targetColEst], n.latent, n.latent); measurementErrorVarsMatrix
  measurementErrorVarsSDMatrix <- matrix(est[targetRow, targetColSD], n.latent, n.latent); measurementErrorVarsSDMatrix
  measurementErrorVarsLLMatrix <- matrix(est[targetRow, targetColLL], n.latent, n.latent); measurementErrorVarsLLMatrix
  measurementErrorVarsULMatrix <- matrix(est[targetRow, targetColUL], n.latent, n.latent); measurementErrorVarsULMatrix
  measurementErrorVars <- measurementErrorVarsSD <- c()
  measurementErrorVarsLL <- measurementErrorVarsUL <- c()
  for (i in 1:n.latent) {
    measurementErrorVars[i] <- measurementErrorVarsMatrix[i,i]
    measurementErrorVarsSD[i] <- measurementErrorVarsSDMatrix[i,i]
    measurementErrorVarsLL[i] <- measurementErrorVarsLLMatrix[i,i]
    measurementErrorVarsUL[i] <- measurementErrorVarsULMatrix[i,i]
  }
  names(measurementErrorVars) <- paste0("measurementErrorVar_", latentNames); measurementErrorVars
  #
  measurementErrorCovs <- measurementErrorCovsSD <- c()
  measurementErrorCovsLL <- measurementErrorCovsUL <- c()
  counter <- 0
  for (i in 1:n.latent) {
    for (j in 1:n.latent) {
      if (i != j) {
        counter <- counter + 1
        measurementErrorCovs[i] <- measurementErrorVarsMatrix[j,i]
        measurementErrorCovsSD[i] <- measurementErrorVarsSDMatrix[j,i]
        measurementErrorCovsLL[i] <- measurementErrorVarsLLMatrix[j,i]
        measurementErrorCovsUL[i] <- measurementErrorVarsULMatrix[j,i]
      }
    }
  }
  names(measurementErrorCovs) <- names(diffusionCovs); measurementErrorCovs
  names(measurementErrorCovs) <- gsub("Diffusion", "MeasurementErrorCov", names(measurementErrorCovs)); measurementErrorCovs
  #
  targetRow <- grep("T0cov", rowNames); targetRow
  T0covMatrix <- matrix(est[targetRow, targetColEst], n.latent, n.latent); T0covMatrix
  T0covSDMatrix <- matrix(est[targetRow, targetColSD], n.latent, n.latent); T0covSDMatrix
  T0covLLMatrix <- matrix(est[targetRow, targetColLL], n.latent, n.latent); T0covLLMatrix
  T0covULMatrix <- matrix(est[targetRow, targetColUL], n.latent, n.latent); T0covULMatrix
  T0Vars <- T0VarsSD <- T0VarsLL <- T0VarsUL <- c()
  for (i in 1:n.latent) {
    T0Vars[i] <- T0covMatrix[i,i]
    T0VarsSD[i] <- T0covSDMatrix[i,i]
    T0VarsLL[i] <- T0covLLMatrix[i,i]
    T0VarsUL[i] <- T0covULMatrix[i,i]
  }
  names(T0Vars) <- paste0("T0Var_", latentNames); T0Vars
  #
  T0Covs <- T0CovsSD <- T0CovsLL <- T0CovsUL <- c()
  counter <- 0
  for (i in 1:n.latent) {
    for (j in 1:n.latent) {
      if (i != j) {
        counter <- counter + 1
        T0Covs[i] <- T0covMatrix[j,i]
        T0CovsSD[i] <- T0covSDMatrix[j,i]
        T0CovsLL[i] <- T0covLLMatrix[j,i]
        T0CovsUL[i] <- T0covULMatrix[j,i]
      }
    }
  }
  names(T0Covs) <- names(diffusionCovs); T0Covs
  names(T0Covs) <- gsub("Diffusion", "T0Cov", names(T0Covs)); T0Covs
  #
  # overwrite previous in case random intercepts were modelled
  InitialVars <- InitialVarsSD <- InitialVarsLL <- InitialVarsUL <- c()
  #if (CoTiMAFiT$summary$randomEffects$popsd[1] != "no random effects estimated") {
  if (length(length(grep("random", CoTiMAFiT$summary$randomEffects$popsd))) > 0) {
    for (i in 1:n.latent) {
      InitialVars[i] <- CoTiMAFiT$summary$randomEffects$popcov_mean[i,i]
      InitialVarsSD[i] <- CoTiMAFiT$summary$randomEffects$model_popcov_sd[i,i]
      if (oldFit) InitialVarsLL[i] <- T0VarsLL[i] else InitialVarsLL[i] <- CoTiMAFiT$summary$randomEffects$model_popcov_025[i,i]
      if (oldFit) InitialVarsUL[i] <- T0VarsUL[i] else InitialVarsUL[i] <- CoTiMAFiT$summary$randomEffects$model_popcov_975[i,i]
    }
    names(InitialVars) <- paste0("InitialVar_", latentNames); InitialVars
    InitialCovs <- InitialCovsSD <- InitialCovsLL <- InitialCovsUL <- c()
    counter <- 0
    for (i in 1:n.latent) {
      for (j in 1:n.latent) {
        if (i != j) {
          counter <- counter + 1
          T0CovsLL
          InitialCovs[counter] <- CoTiMAFiT$summary$randomEffects$popcov_mean[j,i]
          InitialCovsSD[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_sd[i,i]
          if (oldFit) InitialCovsLL[counter] <- T0CovsLL[counter] else InitialCovsLL[i] <- CoTiMAFiT$summary$randomEffects$model_popcov_025[i,i]
          if (oldFit) InitialCovsUL[counter] <- T0CovsUL[counter] else InitialCovsUL[i] <- CoTiMAFiT$summary$randomEffects$model_popcov_975[i,i]
        }
      }
    }
    names(InitialCovs) <- names(diffusionCovs); InitialCovs
    names(InitialCovs) <- gsub("Diffusion", "InitialCov", names(InitialCovs)); InitialCovs
    #
    SlopeVariance_TraitVariance <- SlopeVariance_TraitVarianceSD <- c()
    SlopeVariance_TraitVarianceLL <- SlopeVariance_TraitVarianceUL <- c()
    counter <- 0
    for (i in (n.latent+1):(2*n.latent)) {
      counter <- counter + 1
      SlopeVariance_TraitVariance[counter] <- CoTiMAFiT$summary$randomEffects$popcov_mean[i,i]
      SlopeVariance_TraitVarianceSD[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_sd[i,i]
      if (oldFit) SlopeVariance_TraitVarianceLL[counter] <- NA else SlopeVariance_TraitVarianceLL[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_025[i,i]
      if (oldFit) SlopeVariance_TraitVarianceUL[counter] <- NA else SlopeVariance_TraitVarianceUL[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_975[i,i]
    }
    names(SlopeVariance_TraitVariance) <-paste0("SlopeVariance_TraitVariance_", latentNames); SlopeVariance_TraitVariance
    #
    SlopeCov_TraitCov <- SlopeCov_TraitCovSD <- SlopeCor_TraitCor <- SlopeCor_TraitCorSD <- c()
    SlopeCov_TraitCovLL <- SlopeCov_TraitCovUL <- c()
    SlopeCor_TraitCorLL <- SlopeCor_TraitCorUL <- c()
    counter <- 0
    for (i in (n.latent+1):(2*n.latent)) {
      for (j in (n.latent+1):(2*n.latent)) {
        if (i != j) {
          counter <- counter + 1
          SlopeCov_TraitCov[counter] <- CoTiMAFiT$summary$randomEffects$popcov_mean[j,i]
          SlopeCov_TraitCovSD[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_sd[j,i]
          if (oldFit) SlopeCov_TraitCovLL[counter] <- NA else SlopeCov_TraitCovLL[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_025[j,i]
          if (oldFit) SlopeCov_TraitCovUL[counter] <- NA else SlopeCov_TraitCovUL[counter] <- CoTiMAFiT$summary$randomEffects$model_popcov_975[j,i]
          SlopeCor_TraitCor[counter] <- CoTiMAFiT$summary$randomEffects$popcor_mean[j,i]
          SlopeCor_TraitCorSD[counter]  <- NA
          SlopeCor_TraitCorLL[counter]  <- NA
          SlopeCor_TraitCorUL[counter]  <- NA
        }
      }
    }
    SlopeCov_TraitCovLL
    names(SlopeCov_TraitCov) <- names(diffusionCovs); SlopeCov_TraitCov
    names(SlopeCov_TraitCov) <- gsub("Diffusion", "SlopeCov_TraitCov", names(SlopeCov_TraitCov)); SlopeCov_TraitCov
    names(SlopeCor_TraitCor) <- names(diffusionCovs); SlopeCov_TraitCov
    names(SlopeCor_TraitCor) <- gsub("Diffusion", "SlopeCor_TraitCor", names(SlopeCor_TraitCor)); SlopeCor_TraitCor
    SlopeCor_TraitCorSD <- rep(NA, (n.latent*n.latent-n.latent))
  } else {
    InitialVar <- T0Vars; InitialVar
    InitialVarLL <- T0VarsLL; InitialVarLL
    InitialVarUL <- T0VarsUL; InitialVarUL
    InitialCovs <- T0Covs; InitialCovs
    InitialCovsLL <- T0CovsLL; InitialCovsLL
    InitialCovsUL <- T0CovsUL; InitialCovsUL
    SlopeVariance_TraitVariance <- NA
    SlopeVariance_TraitVarianceLL <- NA
    SlopeVariance_TraitVarianceUL <- NA
    SlopeCov_TraitCov <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCov_TraitCovLL <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCov_TraitCovUL <- rep(NA, (n.latent*n.latent-n.latent))
    InitialVarSD <- T0VarsSD; InitialVarSD
    InitialVarLL <- T0VarsLL; InitialVarLL
    InitialVarUL <- T0VarsUL; InitialVarUL
    InitialCovsSD <- T0CovsSD; InitialCovsSD
    InitialCovsLL <- T0CovsLL; InitialCovsLL
    InitialCovsUL <- T0CovsUL; InitialCovsUL
    SlopeVariance_TraitVarianceSD <- NA
    SlopeVariance_TraitVarianceLL <- NA
    SlopeVariance_TraitVarianceUL <- NA
    SlopeCov_TraitCovSD <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCov_TraitCovLL <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCov_TraitCovUL <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCor_TraitCor <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCor_TraitCorSD <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCor_TraitCorLL <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCor_TraitCorUL <- rep(NA, (n.latent*n.latent-n.latent))
  }
  Minus2LL <- CoTiMAFiT$summary$minus2ll; Minus2LL
  names(Minus2LL) <- "Minus2LL"
  NumEstParams <- CoTiMAFiT$summary$n.parameters; NumEstParams
  names(NumEstParams) <- "NumEstParams"

  #
  tmp1 <- round(c(initialMeans,
                  manifestMeans,
                  slopeMeans_Cint,
                  autos,
                  crosss,
                  proportions,
                  autoregressions,
                  couplings,
                  crosslaggeds,
                  diffusionVars,
                  diffusionCovs,
                  DynErrVars,
                  DynErrCovs,
                  measurementErrorVars,
                  measurementErrorCovs,
                  InitialVars,
                  InitialCovs,
                  SlopeVariance_TraitVariance,
                  SlopeCov_TraitCov,
                  SlopeCor_TraitCor,
                  Minus2LL,
                  NumEstParams), digits); length(tmp1)
  tmp2 <- round(c(initialMeansSD,
                  manifestMeansSD,
                  slopeMeans_CintSD,
                  autosSD,
                  crosssSD,
                  proportionsSD,
                  autoregressionsSD,
                  couplingsSD,
                  crosslaggedsSD,
                  diffusionVarsSD,
                  diffusionCovsSD,
                  DynErrVarsSD,
                  DynErrCovsSD,
                  measurementErrorVarsSD,
                  measurementErrorCovsSD,
                  InitialVarsSD,
                  InitialCovsSD,
                  SlopeVariance_TraitVarianceSD,
                  SlopeCov_TraitCovSD,
                  SlopeCor_TraitCorSD,
                  NA,
                  NA), digits); length(tmp2)
  tmp3 <- round(c(initialMeansLL,
                  manifestMeansLL,
                  slopeMeans_CintLL,
                  autosLL,
                  crosssLL,
                  proportionsLL,
                  autoregressionsLL,
                  couplingsLL,
                  crosslaggedsLL,
                  diffusionVarsLL,
                  diffusionCovsLL,
                  DynErrVarsLL,
                  DynErrCovsLL,
                  measurementErrorVarsLL,
                  measurementErrorCovsLL,
                  InitialVarsLL,
                  InitialCovsLL,
                  SlopeVariance_TraitVarianceLL,
                  SlopeCov_TraitCovLL,
                  SlopeCor_TraitCorLL,
                  NA,
                  NA), digits); length(tmp3)
  tmp4 <- round(c(initialMeansUL,
                  manifestMeansUL,
                  slopeMeans_CintUL,
                  autosUL,
                  crosssUL,
                  proportionsUL,
                  autoregressionsUL,
                  couplingsUL,
                  crosslaggedsUL,
                  diffusionVarsUL,
                  diffusionCovsUL,
                  DynErrVarsUL,
                  DynErrCovsUL,
                  measurementErrorVarsUL,
                  measurementErrorCovsUL,
                  InitialVarsUL,
                  InitialCovsUL,
                  SlopeVariance_TraitVarianceUL,
                  SlopeCov_TraitCovUL,
                  SlopeCor_TraitCorUL,
                  NA,
                  NA), digits); length(tmp4)

  #length(tmp1)
  #tmp2
  resultsTable <- as.matrix(cbind(tmp1, tmp2, tmp3, tmp4))
  resultsTable

  rownames(resultsTable) <- names(tmp1)
  colnames(resultsTable) <- c("est", "SD", "LL", "UL")

  return(resultsTable)
}
