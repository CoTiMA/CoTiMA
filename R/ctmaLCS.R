#' ctmaLCS
#'
#' @description Transforms estimates obtained with \code{\link{ctmaFit}} into LCS (latent change score) terminology.
#' LCS models can be estimated with CT CLPM, but results have to be transformed. When time intervals vary much
#' between and within persons, LCS models are virtually impossible to fit. However, CT CLPM models can be
#' fitted, and the results - after transformation - show what LCS estimates would have been (cf Voelke & Oud,
#' 2015; their terminology to label LCS effects is used in the output created by ctmaLCS)
#'
#' @param CoTiMAFit Fitted CoTiMA object.
#' @param undoTimeScaling Whether (TRUE) or not (FALSE) LCS results should be provided ignoring the scaleTime argument used in ctmaFit.
#' @param digits Number of digits used for rounding (in outputs)
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#'
#' @importFrom OpenMx expm
#' @importFrom RPushbullet pbPost
#' @importFrom stats quantile
#' @importFrom ctsem ctExtract
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
ctmaLCS <- function(CoTiMAFit=NULL, undoTimeScaling=TRUE, digits=4, activateRPB=FALSE) {

  if (is.null(CoTiMAFit$studyFitList$stanfit) & is.null(CoTiMAFit$stanfit)) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object has to be supplied to compute LCS results \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (!(is.null((is.null(CoTiMAFit$studyFitList$stanfit))))) {
    fit <- CoTiMAFit$studyFitList
    CoTiMA <- TRUE
    } else  {
      fit <- CoTiMAFit
      CoTiMA <- FALSE
    }

  if ((undoTimeScaling == TRUE) & CoTiMA) scaleTime <- CoTiMAFit$argumentList$scaleTime else scaleTime <- 1

  e <- ctExtract(fit)
  e$pop_DRIFT <- e$pop_DRIFT * scaleTime

  n.latent <- fit$ctstanmodelbase$n.latent; n.latent
  n.manifest <- fit$ctstanmodelbase$n.manifest; n.manifest

  # get random intercept stats
  if (!(is.null(fit$stanfit$transformedpars$popsd))) {
    model_popsd <- apply(fit$stanfit$transformedpars$popsd, 2, mean); model_popsd
    #e <- ctsem::ctExtract(fitStanctModel)
    model_popcov_m <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits)
    model_popcov_sd <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits)
    model_popcov_T <- round(ctsem::ctCollapse(e$popcov, 1, mean)/ctsem::ctCollapse(e$popcov, 1, stats::sd), digits)
    model_popcov_025 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025))
    model_popcov_50 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .50))
    model_popcov_975 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975))
    # convert to correlations and do the same (array to list then list to array)
    e$popcor <- lapply(seq(dim(e$popcov)[1]), function(x) e$popcov[x , ,])
    e$popcor <- lapply(e$popcor, stats::cov2cor)
    e$popcor <- array(unlist(e$popcor), dim=c(n.latent*2, n.latent*2, length(e$popcor)))
    model_popcor_m <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
    model_popcor_sd <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
    model_popcor_T <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
    model_popcor_025 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .025))
    model_popcor_50 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .50))
    model_popcor_975 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .975))
    #model_popcor <- stats::cov2cor(model_popcov_m)
  } else {
    model_popsd <- "no random effects estimated"
    model_popcov_m <- model_popcov_sd <- model_popcov_T <- model_popcov_025 <- model_popcov_50 <- model_popcov_975 <- "no random effects estimated"
    model_popcor_m <- model_popcor_sd <- model_popcor_T <- model_popcor_025 <- model_popcor_50 <- model_popcor_975 <- "no random effects estimated"
  }

  latentNames <- fit$ctstanmodelbase$latentNames; latentNames
  manifestNames <- fit$ctstanmodelbase$manifestNames; manifestNames

  if (CoTiMA) {
    driftNames <- CoTiMAFit$parameterNames$DRIFT; driftNames
    diffNames <- CoTiMAFit$parameterNames$DIFFUSION; diffNames
    T0varNames <- CoTiMAFit$parameterNames$T0VAR; T0varNames
  } else {
    driftNames <- diffNames <- T0varNames <- c()
    tmp1 <- fit$ctstanmodelbase$latentNames
    for (i in 1:length(tmp1)) {
      for (j in 1:length(tmp1)) {
        driftNames <- c(driftNames, paste0(tmp1[j], "_to_", tmp1[i]))
        diffNames <- c(diffNames, paste0("diff_", tmp1[j], "_", tmp1[i]))
        T0varNames <- c(T0varNames, paste0("T0cov",i,j))
      }
    }
  }
  #driftNames; diffNames; T0varNames
  manifestVarNames <- c()
  for (i in 1:(length(manifestNames))) {
    manifestVarNames <- c(manifestVarNames, paste0(manifestNames[i], "_with_", manifestNames[i]))
  }
  #manifestVarNames
  manifestCovNames <- c()
  for (i in 1:(length(manifestNames)-1)) {
    for (j in (i+1):length(manifestNames)) {
      manifestCovNames <- c(manifestCovNames, paste0(manifestNames[j], "_with_", manifestNames[i]))
    }
  }
  #manifestCovNames

  lowerTri <- lower.tri(matrix(NA, n.latent, n.latent)); lowerTri
  driftNamesMatrix <- matrix(driftNames, n.latent, n.latent, byrow=T); driftNamesMatrix
  #manifestCovNamesMatrix <- matrix(manifestCovNames, length(which(lowerTri==TRUE)), length(which(lowerTri==TRUE)), byrow=T); manifestCovNamesMatrix
  #manifestVarNames

  # get results
  tmp <- summary(fit)
  fit$summary <- tmp

  est <- fit$summary$parmatrices; est
  colNames <- colnames(est); colNames
  rowNames <- est$matrix; rowNames
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
  tmp1 <- which(fit$ctstanmodelbase$pars$matrix == "MANIFESTMEANS"); tmp1
  tmp2 <- fit$ctstanmodelbase$pars$param[tmp1]; tmp2
  tmp1 <- which(fit$ctstanmodelbase$pars$matrix == "CINT"); tmp1
  tmp3 <- fit$ctstanmodelbase$pars$param[tmp1]; tmp3
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
  drift <- matrix(fit$summary$popmeans[grep("to", rownames(fit$summary$popmeans)), 1], n.latent, n.latent, byrow=TRUE) * scaleTime; drift

  ## get SDs
  tmp1 <- grep("DRIFT", rowNames); tmp1
  tmp2 <- grep("dt", rowNames); tmp2
  driftSD <- matrix(est[tmp1[!(tmp1 %in% tmp2)], targetColSD] * scaleTime, n.latent, n.latent, byrow=TRUE); driftSD
  driftLL <- matrix(est[tmp1[!(tmp1 %in% tmp2)], targetColLL] * scaleTime, n.latent, n.latent, byrow=TRUE); driftLL
  driftUL <- matrix(est[tmp1[!(tmp1 %in% tmp2)], targetColUL] * scaleTime, n.latent, n.latent, byrow=TRUE); driftUL
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
  for (j in 1:dim(e$pop_DIFFUSIONcov)[1]) diffusion_samples[[j]] <- e$pop_DIFFUSIONcov[j,,] * scaleTime
  #
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
  for (i in 1:(n.latent-1)) {
    for (j in (i+1):n.latent) {
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

  names(diffusionCovs) <- matrix(diffNames, n.latent, n.latent)[lowerTri]
  names(diffusionCovs) <- gsub("diff_", "Diffusion_", names(diffusionCovs)); diffusionCovs
  names(DynErrCovs) <- names(diffusionCovs); DynErrCovs
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
  for (i in 1:(n.latent-1)) {
    for (j in (i+1):n.latent) {
      if (i != j) {
        counter <- counter + 1
        measurementErrorCovs[i] <- measurementErrorVarsMatrix[j,i]
        measurementErrorCovsSD[i] <- measurementErrorVarsSDMatrix[j,i]
        measurementErrorCovsLL[i] <- measurementErrorVarsLLMatrix[j,i]
        measurementErrorCovsUL[i] <- measurementErrorVarsULMatrix[j,i]
      }
    }
  }
  names(measurementErrorCovs) <- paste0("measurementError_", manifestCovNames); measurementErrorCovs

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
  for (i in 1:(n.latent-1)) {
    for (j in (i+1):n.latent) {
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
  # overwrite previous results in case random intercepts were modelled
  InitialVars <- InitialVarsSD <- InitialVarsLL <- InitialVarsUL <- c()
  if(model_popsd[1] != "no random effects estimated") {
    for (i in 1:n.latent) {
      InitialVars[i] <- model_popcov_m[i,i]
      InitialVarsSD[i] <- model_popcov_sd[i,i]
      InitialVarsLL[i] <- model_popcov_025[i,i]
      InitialVarsUL[i] <- model_popcov_975[i,i]
    }
    names(InitialVars) <- paste0("InitialVar_", latentNames); InitialVars
    InitialCovs <- InitialCovsSD <- InitialCovsLL <- InitialCovsUL <- c()
    counter <- 0
    for (i in 1:(n.latent-1)) {
      for (j in (i+1):n.latent) {
        if (i != j) {
          counter <- counter + 1
          T0CovsLL
          InitialCovs[counter] <- model_popcov_m[j,i]
          InitialCovsSD[counter] <- model_popcov_sd[i,i]
          InitialCovsLL[i] <- model_popcov_025[i,i]
          InitialCovsUL[i] <- model_popcov_975[i,i]
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
      SlopeVariance_TraitVariance[counter] <- model_popcov_m[i,i]
      SlopeVariance_TraitVarianceSD[counter] <- model_popcov_sd[i,i]
      SlopeVariance_TraitVarianceLL[counter] <- model_popcov_025[i,i]
      SlopeVariance_TraitVarianceUL[counter] <- model_popcov_975[i,i]
    }
    names(SlopeVariance_TraitVariance) <-paste0("SlopeVariance_TraitVariance_", latentNames); SlopeVariance_TraitVariance
    #
    SlopeCov_TraitCov <- SlopeCov_TraitCovSD <- SlopeCor_TraitCor <- SlopeCor_TraitCorSD <- c()
    SlopeCov_TraitCovLL <- SlopeCov_TraitCovUL <- c()
    SlopeCor_TraitCorLL <- SlopeCor_TraitCorUL <- c()
    counter <- 0
    for (i in (n.latent+1):(2*n.latent-1)) {
      for (j in (i+1):(2*n.latent)) {
        if (i != j) {
          counter <- counter + 1
          SlopeCov_TraitCov[counter] <- model_popcov_m[j,i]
          SlopeCov_TraitCovSD[counter] <- model_popcov_sd[j,i]
          SlopeCov_TraitCovLL[counter] <- model_popcov_025[j,i]
          SlopeCov_TraitCovUL[counter] <- model_popcov_975[j,i]
          SlopeCor_TraitCor[counter] <- model_popcor_m[j,i]
          SlopeCor_TraitCorSD[counter]  <- NA
          SlopeCor_TraitCorLL[counter]  <- NA
          SlopeCor_TraitCorUL[counter]  <- NA
        }
      }
    }
    names(SlopeCov_TraitCov) <- names(diffusionCovs); SlopeCov_TraitCov
    names(SlopeCov_TraitCov) <- gsub("Diffusion", "SlopeCov_TraitCov", names(SlopeCov_TraitCov)); SlopeCov_TraitCov
    names(SlopeCor_TraitCor) <- names(diffusionCovs); SlopeCov_TraitCov
    names(SlopeCor_TraitCor) <- gsub("Diffusion", "SlopeCor_TraitCor", names(SlopeCor_TraitCor)); SlopeCor_TraitCor
    #SlopeCor_TraitCorSD <- rep(NA, (n.latent*n.latent-n.latent))
    SlopeCor_TraitCorSD <- rep(NA, length(SlopeCor_TraitCor)); SlopeCor_TraitCorSD
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
    SlopeCov_TraitCov <- rep(NA, length(InitialCovs))
    SlopeCov_TraitCovLL <- rep(NA, length(InitialCovs))
    SlopeCov_TraitCovUL <- rep(NA, length(InitialCovs))
    InitialVarSD <- T0VarsSD; InitialVarSD
    InitialVarLL <- T0VarsLL; InitialVarLL
    InitialVarUL <- T0VarsUL; InitialVarUL
    InitialCovsSD <- T0CovsSD; InitialCovsSD
    InitialCovsLL <- T0CovsLL; InitialCovsLL
    InitialCovsUL <- T0CovsUL; InitialCovsUL
    SlopeVariance_TraitVarianceSD <- NA
    SlopeVariance_TraitVarianceLL <- NA
    SlopeVariance_TraitVarianceUL <- NA
    SlopeCov_TraitCovSD <- rep(NA, length(InitialCovs))
    SlopeCov_TraitCovLL <- rep(NA, length(InitialCovs))
    SlopeCov_TraitCovUL <- rep(NA, length(InitialCovs))
    SlopeCor_TraitCor <- rep(NA, length(InitialCovs))
    SlopeCor_TraitCorSD <- rep(NA, length(InitialCovs))
    SlopeCor_TraitCorLL <- rep(NA, length(InitialCovs))
    SlopeCor_TraitCorUL <- rep(NA, length(InitialCovs))
  }
  Minus2LL <- fit$summary$loglik *-2; Minus2LL
  names(Minus2LL) <- "Minus2LL"
  NumEstParams <- fit$summary$npars; NumEstParams
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
                  SlopeCov_TraitCovSD,
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

  resultsTable <- as.matrix(cbind(tmp1, tmp2, tmp3, tmp4))
  resultsTable

  rownames(resultsTable) <- names(tmp1)
  colnames(resultsTable) <- c("est", "SD", "LL", "UL")

  return(resultsTable)
}
