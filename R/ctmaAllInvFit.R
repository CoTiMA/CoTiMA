#' ctmaAllInvFit
#'
#' #' @description Fit a CoTiMA model with all params (drift, T0var, diffusion) invariant across primary studies
#'
#' @param ctmaInitFit ctmaInitFit
#' @param activeDirectory activeDirectory
#' @param activateRPB activateRPB
#' @param digits digits
#' @param coresToUse coresToUse
#' @param drift Labels for drift effects. Have to be either of the type V1toV2 or 0 for effects to be excluded, which is usually not recommended)
#' @param n.manifest Number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).
#' @param indVarying Allows ct intercepts to vary at the individual level (random effects model, accounts for unobserved heterogeneity)
#' @param indVaryingT0 Forces T0MEANS (T0 scores) to vary interindividually, which undos the nesting of T0(co-)variances in primary studies (default = TRUE). Was standard until Aug. 2022. Could provide better/worse estimates if set to FALSE.
#' @param scaleTime scaleTime
#' @param optimize optimize
#' @param nopriors nopriors (TRUE, but deprecated)
#' @param priors priors (FALSE)
#' @param finishsamples finishsamples
#' @param iter iter
#' @param chains chains
#' @param verbose verbose
#' @param loadAllInvFit loadAllInvFit
#' @param saveAllInvFit saveAllInvFit
#' @param silentOverwrite silentOverwrite
#' @param customPar logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param manifestMeans Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param manifestVars define the error variances of the manifests with a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param indVaryingT0 Allows ct intercepts to vary at the individual level (random effects model, accounts for unobserved heterogeneity)
#'
#' @return returns a fitted CoTiMA object, in which all drift parameters, Time 0 variances and covariances, and diffusion parameters were set invariant across primary studies
#'
ctmaAllInvFit <- function(
  ctmaInitFit=NULL,
  activeDirectory=NULL,
  activateRPB=FALSE,
  digits=4,
  drift=drift,
  coresToUse=c(1),
  n.manifest=0,
  indVarying=FALSE,
  scaleTime=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  priors=FALSE,
  finishsamples=NULL,
  iter=NULL,
  chains=NULL,
  verbose=NULL,
  loadAllInvFit=c(),
  saveAllInvFit=c(),
  silentOverwrite=FALSE,
  customPar=FALSE,
  T0means=0,
  manifestMeans=0,
  CoTiMAStanctArgs=NULL,
  lambda=NULL,
  manifestVars=NULL,
  indVaryingT0=TRUE
)
{

  if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMAStanctArgs$verbose}

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object has to be supplied to plot something. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  { # fitting params
    # Added 10. Oct 2022 (17. Aug 2022 in Init fit similar)
    tmp1 <- names(CoTiMA::CoTiMAStanctArgs) %in% names(CoTiMAStanctArgs); tmp1
    tmp2 <- CoTiMA::CoTiMAStanctArgs
    if (!(is.null(CoTiMAStanctArgs))) tmp2[tmp1] <- CoTiMAStanctArgs
    CoTiMAStanctArgs <- tmp2

    if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
    if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
    if (!(is.null(nopriors))) CoTiMAStanctArgs$nopriors <- nopriors
    if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
    if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
    if (!(is.null(iter))) CoTiMAStanctArgs$iter <- iter
    if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
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
    Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.\n"
    message(Msg)
  }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  #####################
  #######################################################################################################################

  {
    n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    lambda <- ctmaInitFit$statisticsList$lambda; lambda
  }

  #######################################################################################################################
  ################################################# data preparation ####################################################
  #######################################################################################################################
  {
    # manifests & latent names
    if (n.manifest > n.latent) {
      manifestNames <- paste0("y", 1:n.manifest); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    } else {
      manifestNames <- paste0("V", 1:n.latent); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    }

    # combine pseudo raw data
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
      #if (CoTiMAStanctArgs$scaleTI == TRUE) dataTmp[ , ncol(dataTmp)] <- scale(dataTmp[ , ncol(dataTmp)]) # CHD 10 Oct 2022 not necessary in AllInv Model
    }
    targetCols <- which(colnames(dataTmp) == "groups"); targetCols
    dataTmp <- dataTmp[ ,-targetCols]

    # add clusters as dummy moderators
    clusCounter <- 0
    cluster.weights <- c()
    cluster.sizes <- c()
    cluster.note <- c()

    n.var <- max(n.manifest, n.latent); n.var

    dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.var, n.TIpred = (n.studies-1+clusCounter),
                                    manifestNames=manifestNames)
    dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
    dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime

    # eliminate rows where ALL latents are NA
    if (n.manifest > n.latent) namePart <- "y" else namePart <- "V"
    dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0(namePart, 1:n.var)], 1, function(x) sum(is.na(x)) != n.var ), ]
    datalong_all <- dataTmp3

  }

  datalong_all <- as.data.frame(datalong_all)

  # get names and params
  namesAndParams <- ctmaLabels(
    n.latent=n.latent,
    n.manifest=n.manifest,
    lambda=lambda,
    T0means=T0means,
    manifestMeans=manifestMeans,
    manifestVars=manifestVars)#,
    #drift=drift)
  driftNames <- namesAndParams$driftNames; driftNames
  driftNames <- gsub(" \\(invariant\\)" , "", driftNames); driftNames
  driftFullNames <- namesAndParams$driftFullNames; driftFullNames
  driftParams <- namesAndParams$driftParams; driftParams
  diffNames <- namesAndParams$diffNames; diffNames
  diffParams <- namesAndParams$diffParams; diffParams
  diffFullNames <- namesAndParams$diffFullNames; diffFullNames
  lambdaParams <- namesAndParams$lambdaParams; lambdaParams
  T0VARParams <- namesAndParams$T0VARParams; T0VARParams
  manifestMeansParams <- namesAndParams$manifestMeansParams; manifestMeansParams
  T0meansParams=namesAndParams$T0meansParams; T0meansParams
  manifestVarsParams <- namesAndParams$manifestVarsParams; manifestVarsParams


  driftParamsTmp <- driftParams; driftParamsTmp
  diffParamsTmp  <- diffParams
  meanLag <- mean(allDeltas, na.rm=TRUE); meanLag
  if ((meanLag > 6) & (customPar)) {
    counter <- 0
    for (h in 1:(n.latent)) {
      for (j in 1:(n.latent)) {
        counter <- counter + 1
        if (h == j) {
          driftParamsTmp[counter] <- paste0(driftParamsTmp[counter], paste0("|-log1p_exp(-param *.1 -2)"))
          diffParamsTmp[counter] <- paste0(diffParamsTmp[counter], paste0("|log1p_exp(param *.1 +2)"))
        }
      }
    }
  }

  #tmp0 <- matrix(diffNamesTmp, n.latent); tmp0
  #tmp0[upper.tri(tmp0, diag=FALSE)] <- 0; tmp0
  #diffNamesTmp <- tmp0; diffNamesTmp

  # taken out and replaced 28 Sep 2022
  #if (indVarying == TRUE) {
  #  T0MEANS <- "auto"
  #  #T0MEANS <- matrix(0, nrow=n.latent, ncol=1)
  #  MANIFESTMEANS <- "auto"
  #  #MANIFESTVAR <- "auto"
  #  MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent)
  #} else {
  #  T0MEANS <- matrix(c(0), nrow = n.latent, ncol = 1)
  #  MANIFESTMEANS <- matrix(c(0), nrow = n.latent, ncol = 1)
  #  MANIFESTVAR <- matrix(0, nrow=n.latent, ncol=n.latent)
  #}
  T0MEANS <- T0meansParams
  MANIFESTMEANS <- manifestMeansParams
  MANIFESTVAR <- manifestVarsParams

  # CHD 13.6.2023
  # CHD 9.6.2023
  if ((indVarying == 'cint') | (indVarying == 'Cint')) indVarying <- 'CINT'

    # CHD 9.6.2023
    if ( (indVarying == 'CINT') & (indVaryingT0 == TRUE) ) {
      print(paste0("#################################################################################"))
      print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
      print(paste0("#################################################################################"))

      print(paste0("#################################################################################"))
      print(paste0("# T0means are set to \'auto\'. T0(co-)variances not modelled nested in primaries.#"))
      print(paste0("#################################################################################"))
      T0meansParams <- 'auto'

      print(paste0("#################################################################################"))
      print(paste0("####################### CT intercepts are set free.  ########################"))
      print(paste0("#################################################################################"))

      CINTParams <- c()
      for (c in 1:n.latent) {
        CINTParams <- c(CINTParams, paste0("cintV", c))
      }
    }
    #CINTParams

    if ( (indVarying == 'CINT') & (indVaryingT0 == FALSE) ) {
      print(paste0("#################################################################################"))
      print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
      print(paste0("#################################################################################"))

      print(paste0("#################################################################################"))
      print(paste0("### T0means are set to 0. T0(co-)variances are modelled nested in primaries. ####"))
      print(paste0("#################################################################################"))
      T0meansParams <- 0

      print(paste0("#################################################################################"))
      print(paste0("####################### CT intercepts are set free.  ########################"))
      print(paste0("#################################################################################"))

      CINTParams <- c()
      for (c in 1:n.latent) {
        CINTParams <- c(CINTParams, paste0("cintV", c))
      }
    }

    if ( (indVarying == TRUE) & (indVaryingT0 == TRUE) ) {
      print(paste0("#################################################################################"))
      print(paste0("###### Just a note: Individually varying manifest means model requested.  #######"))
      print(paste0("#################################################################################"))

      print(paste0("#################################################################################"))
      print(paste0("## T0means set to \'auto\'. T0(co-)variances not modelled nested in primaries. ##"))
      print(paste0("#################### Consider setting \'indVaryingT0 = FALSE\' ####################"))
      print(paste0("#################################################################################"))
      T0meansParams <- 'auto'

      print(paste0("#################################################################################"))
      print(paste0("######### Manifest means (as replacement for intercepts) are set free.  #########"))
      print(paste0("#################################################################################"))

      manifestMeansParams <- 'auto'
    }


    if ( (indVarying == TRUE) & (indVaryingT0 == FALSE) ) {
      print(paste0("#################################################################################"))
      print(paste0("###### Just a note: Individually varying manifest means model requested.  #######"))
      print(paste0("#################################################################################"))

      print(paste0("#################################################################################"))
      print(paste0("### T0means are set to 0. T0(co-)variances are modelled nested in primaries. ####"))
      print(paste0("#################################################################################"))
      T0meansParams <- 0

      print(paste0("#################################################################################"))
      print(paste0("######### Manifest means (as replacement for intercepts) are set free.  #########"))
      print(paste0("#################################################################################"))

      manifestMeansParams <- 'auto'
    }


  # all fixed model is a model with no TI predictors (identical to ctsemModel)
  # CHD added on 24.6.2023 to prevent warning messages
  if (T0MEANS == 0) T0MEANS <- matrix(0, nrow=n.latent, ncol=1)
  if (MANIFESTMEANS == 0) MANIFESTMEANS <- matrix(0, nrow=n.latent, ncol=1)
  if (MANIFESTVAR == 0) MANIFESTVAR <- matrix(0, nrow=n.latent, ncol=n.latent)
  # CHD allFixedModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpointsModel, manifestNames=manifestNames,    # 2 waves in the template only
  allFixedModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                           DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                           DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                           LAMBDA=lambdaParams,
                           type='stanct',
                           CINT=matrix(0, nrow=n.latent, ncol=1),
                           T0MEANS=T0MEANS,
                           MANIFESTMEANS=MANIFESTMEANS,
                           MANIFESTVAR=MANIFESTVAR)

  if (indVarying == FALSE) allFixedModel$pars[, "indvarying"] <- FALSE
  #CHD 13.6.2023
  if (indVaryingT0 == TRUE) {
    allFixedModel$pars[allFixedModel$pars$matrix %in% 'T0MEANS','indvarying'] <- TRUE
  } else {
    allFixedModel$pars[allFixedModel$pars$matrix %in% 'T0MEANS','indvarying'] <- FALSE
  }
  # CHD 13.6.2023
  if (indVarying == 'CINT') {
    allFixedModel$pars[allFixedModel$pars$matrix %in% 'CINT','indvarying'] <- TRUE
  } else {
    allFixedModel$pars[allFixedModel$pars$matrix %in% 'CINT','indvarying'] <- FALSE
  }
  if (indVarying == TRUE) {
    allFixedModel$pars[allFixedModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- TRUE
  } else {
    allFixedModel$pars[allFixedModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- FALSE
  }



  # LOAD or Fit
  if (length(loadAllInvFit) > 0) {
    x1 <- paste0(activeDirectory, loadAllInvFit[1], ".rds"); x1
    results <- readRDS(file=x1)
  } else {
    allFixedModelFit <- suppressMessages(ctsem::ctStanFit(
      datalong = datalong_all,
      ctstanmodel = allFixedModel,
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
      cores=coresToUse))

    Msg <- "\nComputing results summary of all invariant model.\n"
    message(Msg)
    allFixedModelFitSummary <- summary(allFixedModelFit, digits=digits)
  }

  # SAVE
  if (length(saveAllInvFit) > 0)  {
    x1 <- paste0(saveAllInvFit[1], ".rds"); x1
    x2 <- paste0(activeDirectory); x2
    ctmaSaveFile(activateRPB, "", allFixedModelFit, x1, x2, silentOverwrite=silentOverwrite)
  }

  ### Extract estimates & statistics
  # account for changes in ctsem 3.4.1
  {
    if ("matrix" %in% colnames(allFixedModelFitSummary$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
    #if (ctsem341) driftNamesTmp <-c(matrix(driftNames, n.latent, byrow=TRUE)) else driftNamesTmp <- driftNames
    tmpMean <- grep("ean", colnames(allFixedModelFitSummary$parmatrices)); tmpMean
    tmpSd <- tmpMean+1; tmpSd
    Tvalues <- allFixedModelFitSummary$parmatrices[,tmpMean]/allFixedModelFitSummary$parmatrices[,tmpSd]; Tvalues
    allFixedDrift_Coeff <- cbind(allFixedModelFitSummary$parmatrices, Tvalues); allFixedDrift_Coeff
    allFixedDrift_Coeff[, tmpMean:(dim(allFixedDrift_Coeff)[2])] <- round(allFixedDrift_Coeff[, tmpMean:(dim(allFixedDrift_Coeff)[2])], digits); allFixedDrift_Coeff
    # re-label
    if (ctsem341) {
      tmp1 <- which(allFixedDrift_Coeff[, "matrix"] == "DRIFT")
      driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); driftNamesTmp
      rownames(allFixedDrift_Coeff) <- paste0(allFixedDrift_Coeff[, c("matrix")], "_",
                                              allFixedDrift_Coeff[, c("row")], "_",
                                              allFixedDrift_Coeff[, c("col")])
    } else {
      tmp1 <- which(rownames(allFixedDrift_Coeff) == "DRIFT")
      driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); driftNamesTmp
    }
    rownames(allFixedDrift_Coeff)[tmp1] <- paste0("DRIFT ", driftNamesTmp); allFixedDrift_Coeff

    #tmp <- grep("toV", rownames(allFixedModelFitSummary$popmeans)); tmp
    if (ctsem341) {
      homAll_Drift_Coef <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DRIFT", tmpMean]; homAll_Drift_Coef
      homAll_Drift_SE <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DRIFT", tmpSd]; homAll_Drift_SE
      tmp1 <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DRIFT", "2.5%"]; tmp1
      tmp2 <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DRIFT", "97.5%"]; tmp2
    } else {
      homAll_Drift_Coef <- allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DRIFT", "Mean"]; homAll_Drift_Coef
      homAll_Drift_SE <-allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DRIFT", "Sd"]; homAll_Drift_SE
      tmp1 <- c(matrix(allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=TRUE)); tmp2
    }
    names(homAll_Drift_Coef) <- driftNamesTmp; homAll_Drift_Coef
    names(homAll_Drift_SE) <- driftNamesTmp; homAll_Drift_SE
    homAll_Drift_CI <- c(rbind(tmp1, tmp2)); homAll_Drift_CI
    tmp3 <- c(rbind(paste0(driftNamesTmp, "LL"),
                    paste0(driftNamesTmp, "UL"))); tmp3
    names(homAll_Drift_CI) <- tmp3; homAll_Drift_CI
    homAll_Drift_Tvalue <- homAll_Drift_Coef/homAll_Drift_SE; homAll_Drift_Tvalue


    if (ctsem341) {
      homAll_Diffusion_Coef <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", tmpMean]; homAll_Diffusion_Coef
      homAll_Diffusion_SE <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", tmpSd]; homAll_Diffusion_SE
      tmp1 <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "2.5%"]; tmp1
      tmp2 <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "97.5%"]; tmp2
    } else {
      homAll_Diffusion_Coef <- allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DIFFUSIONcov", "Mean"]; homAll_Diffusion_Coef
      homAll_Diffusion_SE <-allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DIFFUSIONcov", "Sd"]; homAll_Diffusion_SE
      tmp1 <- c(matrix(allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DIFFUSIONcov", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "DIFFUSIONcov", "97.5%"], n.latent, byrow=TRUE)); tmp2
      homAll_Diffusion_Coef <- OpenMx::vech2full(homAll_Diffusion_Coef)
      homAll_Diffusion_SE <- OpenMx::vech2full(homAll_Diffusion_SE)
      tmp1 <- OpenMx::vech2full(tmp1)
      tmp2 <- OpenMx::vech2full(tmp2)
    }
    names(homAll_Diffusion_Coef) <- driftNamesTmp; homAll_Diffusion_Coef
    names(homAll_Diffusion_SE) <- driftNamesTmp; homAll_Diffusion_SE
    homAll_Diffusion_CI <- c(rbind(tmp1, tmp2)); homAll_Diffusion_CI
    tmp3 <- c(rbind(paste0(driftNamesTmp, "LL"),
                    paste0(driftNamesTmp, "UL"))); tmp3
    names(homAll_Diffusion_CI) <- tmp3; homAll_Diffusion_CI
    homAll_Diffusion_Tvalue <- homAll_Diffusion_Coef/homAll_Diffusion_SE; homAll_Diffusion_Tvalue


    if (ctsem341) {
      homAll_T0VAR_Coef <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "T0cov", tmpMean]; homAll_T0VAR_Coef
      homAll_T0VAR_SE <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "T0cov", tmpSd]; homAll_T0VAR_SE
      tmp1 <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "T0cov", "2.5%"]; tmp1
      tmp2 <- allFixedModelFitSummary$parmatrices[allFixedModelFitSummary$parmatrices[, "matrix"] == "T0cov", "97.5%"]; tmp2
    } else {
      homAll_T0VAR_Coef <- allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "T0VAR", "Mean"]; homAll_T0VAR_Coef
      homAll_T0VAR_SE <-allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "T0VAR", "Sd"]; homAll_T0VAR_SE
      tmp1 <- c(matrix(allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "T0VAR", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(allFixedModelFitSummary$parmatrices[rownames(allFixedModelFitSummary$parmatrices) == "T0VAR", "97.5%"], n.latent, byrow=TRUE)); tmp2
      homAll_T0VAR_Coef <- OpenMx::vech2full(homAll_T0VAR_Coef)
      homAll_T0VAR_SE <- OpenMx::vech2full(homAll_T0VAR_SE)
      tmp1 <- OpenMx::vech2full(tmp1)
      tmp2 <- OpenMx::vech2full(tmp2)
    }
    names(homAll_T0VAR_Coef) <- driftNamesTmp; homAll_T0VAR_Coef
    names(homAll_T0VAR_SE) <- driftNamesTmp; homAll_T0VAR_SE
    homAll_T0VAR_CI <- c(rbind(tmp1, tmp2)); homAll_T0VAR_CI
    tmp3 <- c(rbind(paste0(driftNamesTmp, "LL"),
                    paste0(driftNamesTmp, "UL"))); tmp3
    names(homAll_Diffusion_CI) <- tmp3; homAll_Diffusion_CI
    homAll_T0VAR_Tvalue <- homAll_T0VAR_Coef/homAll_T0VAR_SE; homAll_T0VAR_Tvalue

    ## Extract Model Fit
    homAll_Minus2LogLikelihood <- -2 * allFixedModelFitSummary$loglik; homAll_Minus2LogLikelihood
    homAll_estimatedParameters <- allFixedModelFitSummary$npars; homAll_estimatedParameters
    homAll_df <- NULL

    homAll_effects <- matrix(t(cbind((homAll_Drift_Coef), (homAll_Drift_SE),
                                     (homAll_Drift_Tvalue))), 1, 3*length(driftNames), byrow=T); homAll_effects

    homAll_effects <- rbind(homAll_effects,
                            matrix(t(cbind(homAll_Diffusion_Coef, homAll_Diffusion_SE,
                                           homAll_Diffusion_Tvalue)), 1, 3*length(driftNames), byrow=T)); homAll_effects

    homAll_effects <- rbind(homAll_effects,
                            matrix(t(cbind(homAll_T0VAR_Coef, homAll_T0VAR_SE,
                                           homAll_T0VAR_Tvalue)), 1, 3*length(driftNames), byrow=T)); homAll_effects

    # Label summary table
    rownames(homAll_effects) <- c("Fixed Effects Drift", "Fixed Effects Diffusion", "Fixed Effects T0Var")
    newColNames <- c()
    for (j in 1:n.latent) {
      for (h in 1:n.latent) {
        newColNames <- c(newColNames, paste0("V",h,"toV", j), "(SE)", "Tvalue")
      }
    }
    colnames(homAll_effects) <- newColNames; homAll_effects
  }

  DRIFT <- matrix(homAll_Drift_Coef, n.latent, n.latent, byrow=TRUE); DRIFT
  DIFFUSION <- matrix(homAll_Diffusion_Coef, n.latent); DIFFUSION
  T0VAR <- matrix(homAll_T0VAR_Coef, n.latent); T0VAR

  results <- list(plot.type="drift",  model.type="stanct",
                  coresToUse=coresToUse, n.studies=1,
                  n.latent=n.latent,
                  #studyList=ctmaInitFit$studyList,
                  studyFitList=list(allFixedModelFit),
                  data=datalong_all,
                  #statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=DRIFT, DIFFUSION=DIFFUSION, T0VAR=T0VAR, CINT=NULL),
                  #parameterNames=ctmaInitFit$parameterNames,
                  CoTiMAStanctArgs=CoTiMAStanctArgs,
                  #invariantDrift=invariantDrift,
                  summary=list(model="All effect (drift, diffusion, T0var) invariant across primary studies.",
                               parmatrices=allFixedModelFitSummary$parmatrices,
                               #randomEffects=invariantDriftStanctFit$popsd,
                               minus2ll= homAll_Minus2LogLikelihood,
                               n.parameters = homAll_estimatedParameters))
  class(results) <- "CoTiMAFit"

  invisible(results)

}
