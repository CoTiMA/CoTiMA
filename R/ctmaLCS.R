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

  if (!(is.null(CoTiMAFit$studyFitList$stanfit))) {
    fit <- CoTiMAFit$studyFitList
    CoTiMA <- TRUE
  } else  {
    fit <- CoTiMAFit
    CoTiMA <- FALSE
  }

  if (is.null(fit$stanfit$transformedpars$popsd)) {
    ErrorMsg <- "\nThe fitted CoTiMA object or ctsem object did not include T0means and ct intercepts. LCS results cannot be computed. Fit an appropriate CoTiMA/ctsem model. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  n.latent <- fit$ctstanmodelbase$n.latent; n.latent
  n.manifest <- fit$ctstanmodelbase$n.manifest; n.manifest

  if (dim(fit$stanfit$transformedpars$popsd)[2] != (n.latent*2)) {
    ErrorMsg <- "\nAThe fitted CoTiMA object or ctsem object did not included both(!) T0means and ct intercepts. LCS results cannot be computed. Fit an appropriate CoTiMA/ctsem model. nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if ((undoTimeScaling == TRUE) & (CoTiMA == TRUE)) scaleTime <- CoTiMAFit$argumentList$scaleTime else scaleTime <- 1
  if (is.numeric(undoTimeScaling)) scaleTime <- undoTimeScaling
  if (is.null(scaleTime)) scaleTime <- 1

  #e <- ctExtract(CoTiMAFullFit_Ov4_0531_CHD_mmri$studyFitList)
  # CHD 28. Aug 2023 (not really necessary, but for documentation that nopriors argument in ctsem changed to priors argument)
  if (is.null(fit$standata$priors)) fit$standata$priors <- 0
  #
  e <- ctExtract(fit)
  e$pop_DRIFT <- e$pop_DRIFT * scaleTime

  # get random intercept stats
  tmp1 <- ctsem::ctCollapse(e$pop_CINT, 1, mean); tmp1
  tmp2 <- ctsem::ctCollapse(e$pop_CINT, 1, sd); tmp2
  tmp3 <- ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, mean); tmp3
  tmp4 <- ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, sd); tmp4
  if ( (n.latent == n.manifest) & (!(any(c(tmp3, tmp4) == 0))) & (all(c(tmp1, tmp2) == 0)) ) mmRI <- TRUE else mmRI <- FALSE


  if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
    #### IDEA: Transformation matrix describing popcov_MM into popciv_cint transformations and then popcov_cint <- trans %*% popcov_mm %*% t(trans) (see: #https://stats.stackexchange.com/questions/113700/covariance-of-a-random-vector-after-a-linear-transformation)
    print("Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) are calculated based on a model with individually varying manifest means instead of Cints.")
    e$popcov_est <- e$popcov
    e$popcov_est[!(is.na(e$popcov_est))] <- NA
    UL <- UR <- diag(1, n.latent, n.latent); UL
    LL <- matrix(0, n.latent, n.latent); LL
    for (k in 1:(dim(e$popcov_est)[1])) {
      LR <- -e$pop_DRIFT[k,,]; LR
      trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
      e$popcov_est[k, , ] <- trans %*% e$popcov[k, , ] %*% t(trans)
    }
    e$popcov <- e$popcov_est
    message <- "Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) were inferred from a model with individually varying manifest means instead of Cints."
  }
  ctCollapse(e$popcov, 1, mean)

  model_popcov_m <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits); model_popcov_m
  model_popcov_sd <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits); model_popcov_sd
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

  # get results
  tmp <- summary(fit)
  fit$summary <- tmp

  if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
    e$pop_T0MEANS_est <- e$pop_T0MEANS[ ,1:n.latent, ]        # 1.n.latent => eliminates last diminsion, which is not required
    e$pop_T0MEANS_est[!(is.na(e$pop_T0MEANS_est))] <- NA
    e$pop_T0MEANS <- e$pop_T0MEANS[ ,1:n.latent, ]
    e$pop_MANIFESTMEANS <- e$pop_MANIFESTMEANS[ ,1:n.latent,]
    for (i in 1:(dim(e$pop_T0MEANS_est)[1])) e$pop_T0MEANS_est[i,] <- e$pop_T0MEANS[i,] + e$pop_MANIFESTMEANS[i,]
    e$pop_T0MEANS <- e$pop_T0MEANS_est
    e$pop_MANIFESTMEANS_backup <- e$pop_MANIFESTMEANS
    e$pop_MANIFESTMEANS[e$pop_MANIFESTMEANS != 0] <- 0
  }
  #
  if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
    initialMeans <- round(ctsem::ctCollapse(e$pop_T0MEANS, 1, mean), digits = digits); initialMeans
    initialMeansSD <- round(ctsem::ctCollapse(e$pop_T0MEANS, 1, stats::sd), digits = digits); initialMeansSD
    initialMeansLL <- ctsem::ctCollapse(e$pop_T0MEANS, 1, function(x) stats::quantile(x, .025)); initialMeansLL
    initialMeansUL <- ctsem::ctCollapse(e$pop_T0MEANS, 1, function(x) stats::quantile(x, .975)); initialMeansUL
  } else {
    initialMeans <- round(ctsem::ctCollapse(e$pop_T0MEANS, 1, mean), digits = digits)[1:n.latent]; initialMeans
    initialMeansSD <- round(ctsem::ctCollapse(e$pop_T0MEANS, 1, stats::sd), digits = digits)[1:n.latent]; initialMeansSD
    initialMeansLL <- ctsem::ctCollapse(e$pop_T0MEANS, 1, function(x) stats::quantile(x, .025))[1:n.latent]; initialMeansLL
    initialMeansUL <- ctsem::ctCollapse(e$pop_T0MEANS, 1, function(x) stats::quantile(x, .975))[1:n.latent]; initialMeansUL
  }
  names(initialMeans) <- paste0("InitialMean_", latentNames); initialMeans

  #
  manifestMeans <- round(ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, mean), digits = digits); manifestMeans
  manifestMeansSD <- round(ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, stats::sd), digits = digits); manifestMeansSD
  manifestMeansLL <- ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, function(x) stats::quantile(x, .025)); manifestMeansLL
  manifestMeansUL <- ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, function(x) stats::quantile(x, .975)); manifestMeansUL
  names(manifestMeans) <- paste0("ManifestMean_", latentNames); manifestMeans

  #
  if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
    for (k in 1:dim(e$pop_MANIFESTMEANS_backup)[1]) {
      e$pop_CINT[k,,1] <- -e$pop_DRIFT[k,,] %*% e$pop_MANIFESTMEANS_backup[k,]
    }
  }
  slopeMeans_Cint <- round(ctsem::ctCollapse(e$pop_CINT, 1, mean), digits = digits); slopeMeans_Cint
  slopeMeans_CintSD <- round(ctsem::ctCollapse(e$pop_CINT, 1, stats::sd), digits = digits); slopeMeans_CintSD
  slopeMeans_CintLL <- ctsem::ctCollapse(e$pop_CINT, 1, function(x) stats::quantile(x, .025)); slopeMeans_CintLL
  slopeMeans_CintUL <- ctsem::ctCollapse(e$pop_CINT, 1, function(x) stats::quantile(x, .975)); slopeMeans_CintUL
  names(slopeMeans_Cint) <- paste0("SlopeMean_Cint_", latentNames); slopeMeans_Cint

  #
  #drift <- matrix(fit$summary$popmeans[grep("to", rownames(fit$summary$popmeans)), 1], n.latent, n.latent, byrow=TRUE) * scaleTime; drift
  drift <- ctCollapse(e$pop_DRIFT, 1, mean); drift
  driftSD <- ctCollapse(e$pop_DRIFT, 1, sd); driftSD
  driftLL <- ctsem::ctCollapse(e$pop_DRIFT, 1, function(x) stats::quantile(x, .025)); driftLL
  driftUL <- ctsem::ctCollapse(e$pop_DRIFT, 1, function(x) stats::quantile(x, .975)); driftUL
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
        crosslaggeds[counter] <- mean(unlist(lapply(expm_drift_samples, function(x) x[j,i])))
        crosslaggedsSD[counter] <- sd(unlist(lapply(expm_drift_samples, function(x) x[j,i])))
        crosslaggedsLL[counter] <- stats::quantile(unlist(lapply(expm_drift_samples, function(x) x[j,i])), prob=.025)
        crosslaggedsUL[counter] <- stats::quantile(unlist(lapply(expm_drift_samples, function(x) x[j,i])), prob=.975)
        couplings[counter] <- crosslaggeds[counter]
        couplingsSD[counter] <- crosslaggedsSD[counter]
        couplingsLL[counter] <- crosslaggedsLL[counter]
        couplingsUL[counter] <- crosslaggedsUL[counter]
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
  measurementErrorVarsMatrix <- round(ctsem::ctCollapse(e$pop_MANIFESTcov, 1, mean), digits = digits); measurementErrorVarsMatrix
  measurementErrorVarsSDMatrix <- round(ctsem::ctCollapse(e$pop_MANIFESTcov, 1, sd), digits = digits); measurementErrorVarsSDMatrix
  measurementErrorVarsLLMatrix <- ctsem::ctCollapse(e$pop_MANIFESTcov, 1, function(x) stats::quantile(x, .025)); measurementErrorVarsLLMatrix
  measurementErrorVarsULMatrix <- ctsem::ctCollapse(e$pop_MANIFESTcov, 1, function(x) stats::quantile(x, .975)); measurementErrorVarsULMatrix
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
  InitialVars <- InitialVarsSD <- InitialVarsLL <- InitialVarsUL <- c()
  for (i in 1:n.latent) {
    InitialVars[i] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits)[i,i];
    InitialVarsSD[i] <- round(ctsem::ctCollapse(e$popcov, 1, sd), digits = digits)[i,i]
    InitialVarsLL[i] <- round(ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025)), digits = digits)[i,i]
    InitialVarsUL[i] <- round(ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975)), digits = digits)[i,i]
  }
  names(InitialVars) <- paste0("InitialVar_", latentNames); InitialVars
  #
  InitialCovs <- InitialCovsSD <- InitialCovsLL <- InitialCovsUL <- c()
  counter <- 0
  for (i in 1:(n.latent-1)) {
    for (j in (i+1):n.latent) {
      if (i != j) {
        counter <- counter + 1
        InitialCovs[counter] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits)[j,i];
        InitialCovsSD[counter] <- round(ctsem::ctCollapse(e$popcov, 1, sd), digits = digits)[j,i]
        InitialCovsLL[counter] <- round(ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025)), digits = digits)[j,i]
        InitialCovsUL[counter] <- round(ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975)), digits = digits)[j,i]
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
        SlopeCor_TraitCorSD[counter]  <- model_popcor_sd[j,i]
        SlopeCor_TraitCorLL[counter]  <- model_popcor_025[j,i]
        SlopeCor_TraitCorUL[counter]  <- model_popcor_975[j,i]
      }
    }
  }
  names(SlopeCov_TraitCov) <- names(diffusionCovs); SlopeCov_TraitCov
  names(SlopeCov_TraitCov) <- gsub("Diffusion", "SlopeCov_TraitCov", names(SlopeCov_TraitCov)); SlopeCov_TraitCov
  names(SlopeCor_TraitCor) <- names(diffusionCovs); SlopeCov_TraitCov
  names(SlopeCor_TraitCor) <- gsub("Diffusion", "SlopeCor_TraitCor", names(SlopeCor_TraitCor)); SlopeCor_TraitCor
  names(SlopeCor_TraitCorSD) <- names(SlopeCov_TraitCov)
  names(SlopeCor_TraitCorSD) <- gsub("Cov", "Cor", names(SlopeCor_TraitCorSD) )


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
  tmp1
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

  resultsTable <- as.matrix(cbind(tmp1, tmp2, tmp3, tmp4))
  #resultsTable

  rownames(resultsTable) <- names(tmp1)
  colnames(resultsTable) <- c("est", "SD", "LL", "UL")

  return(resultsTable)
}
