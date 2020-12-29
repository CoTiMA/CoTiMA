#' ctmaAllInvFit
#'
#' @param ctmaInitFit ctmaInitFit
#' @param activeDirectory activeDirectory
#' @param activateRPB activateRPB
#' @param digits digits
#' @param coresToUse coresToUse
#' @param scaleTime scaleTime
#' @param optimize optimize
#' @param nopriors nopriors
#' @param finishsamples finishsamples
#' @param iter iter
#' @param chains chains
#' @param verbose verbose
#' @param loadAllInvFit loadAllInvFit
#' @param saveAllInvFit saveAllInvFit
#' @param silentOverwrite silentOverwrite
#' @param customPar logical. Leverages the first pass using priors and ensure that the drift diagonal cannott easily go too negative (could help with ctsem > 3.4)
#'
#' @export ctmaAllInvFit
#'
ctmaAllInvFit <- function(
  ctmaInitFit=NULL,
  activeDirectory=NULL,
  activateRPB=FALSE,
  digits=4,
  coresToUse=c(1),
  scaleTime=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  finishsamples=NULL,
  iter=NULL,
  chains=NULL,
  verbose=NULL,
  loadAllInvFit=c(),
  saveAllInvFit=c(),
  silentOverwrite=FALSE,
  customPar=TRUE
)
{

  if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMAStanctArgs$verbose}

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
    stop("Good luck for the next try!")
  }

  { # fitting params
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
    cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
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
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    if (is.null(manifestNames)) n.manifest <- 0 else n.manifest <- length(manifestNames)
    driftNames <- ctmaInitFit$parameterNames$DRIFT; driftNames
    #if (!(is.null(drift))) driftNames <- c(t(matrix(drift, n.latent, n.latent)))
    tmp1 <- names(ctmaInitFit$modelResults$DIFFUSION[[1]]); tmp1
    if (length(tmp1) != n.latent^2) {
      diffNames <- OpenMx::vech2full(tmp1); diffNames
    } else {
      diffNames <- matrix(names(ctmaInitFit$modelResults$DIFFUSION[[1]]), n.latent); diffNames
    }

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

  # scale Drift to cover changes in ctsem 3.4.1 (this would be for ctmaFit/ctmaModFit, but for Init individual study modification is done later)
  driftNamesTmp <- driftNames
  diffNamesTmp  <- diffNames
  meanLag <- mean(allDeltas, na.rm=TRUE); meanLag
  if ((meanLag > 6) & (customPar)) {
    counter <- 0
    for (h in 1:(n.latent)) {
      for (j in 1:(n.latent)) {
        counter <- counter + 1
        #if (h == j) driftNamesTmp[counter] <- paste0(driftNamesTmp[counter], paste0("|-log1p_exp(-param *.1 -2)"))
        if (h == j) {
          driftNamesTmp[counter] <- paste0(driftNamesTmp[counter], paste0("|-log1p_exp(-param *.1 -2)"))
          diffNamesTmp[counter] <- paste0(diffNamesTmp[counter], paste0("|log1p_exp(param *.1 +2)"))
        }
      }
    }
  }

  tmp0 <- matrix(diffNamesTmp, n.latent); tmp0
  tmp0[upper.tri(tmp0, diag=FALSE)] <- 0; tmp0
  diffNamesTmp <- tmp0; diffNamesTmp

  # all fixed model is a model with no TI predictors (identical to ctsemModel)
  # CHD allFixedModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpointsModel, manifestNames=manifestNames,    # 2 waves in the template only
  allFixedModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                           DIFFUSION=matrix(diffNamesTmp, nrow=n.latent, ncol=n.latent),
                           DRIFT=matrix(driftNamesTmp, nrow=n.latent, ncol=n.latent, byrow=FALSE),
                           LAMBDA=diag(n.latent),
                           type='stanct',
                           #CINT=matrix(cintNames, nrow=n.latent, ncol=1),
                           CINT=matrix(0, nrow=n.latent, ncol=1),
                           T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                           MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                           MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent))

  allFixedModel$pars[, "indvarying"] <- FALSE


  # LOAD or Fit
  if (length(loadAllInvFit) > 0) {
    x1 <- paste0(activeDirectory, loadAllInvFit[1], ".rds"); x1
    results <- readRDS(file=x1)
  } else {
    allFixedModelFit <- ctStanFit(
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
      cores=coresToUse)

    cat( "\n", "Computing results summary of all invariant model.", "\n")
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
    #allFixedModelFitSummary
    homAll_Minus2LogLikelihood <- -2 * allFixedModelFitSummary$loglik; homAll_Minus2LogLikelihood
    homAll_estimatedParameters <- allFixedModelFitSummary$npars; homAll_estimatedParameters
    #homAll_df <- ctmaInitFit$summary$df+(ctmaInitFit$summary$n.parameters-homAll_estimatedParameters); homAll_df
    homAll_df <- NULL

    homAll_effects <- matrix(t(cbind((homAll_Drift_Coef), (homAll_Drift_SE),
                                     (homAll_Drift_Tvalue))), 1, 3*length(driftNames), byrow=T); homAll_effects

    #homAll_effects <- rbind(homAll_effects,
    #                        matrix(t(cbind((c(OpenMx::vech2full(homAll_Diffusion_Coef))),
    #                                       c(OpenMx::vech2full((homAll_Diffusion_SE))),
    #                                       c(OpenMx::vech2full((homAll_Diffusion_Tvalue))) )), 1, 3*length(driftNames), byrow=T)); homAll_effects
    homAll_effects <- rbind(homAll_effects,
                            matrix(t(cbind(homAll_Diffusion_Coef, homAll_Diffusion_SE,
                                           homAll_Diffusion_Tvalue)), 1, 3*length(driftNames), byrow=T)); homAll_effects

    #homAll_effects <- rbind(homAll_effects,
    #                        matrix(t(cbind(c(OpenMx::vech2full((homAll_T0Var_Coef))),
    #                                       c(OpenMx::vech2full((homAll_T0Var_SE))),
    #                                       c(OpenMx::vech2full((homAll_T0Var_Tvalue))) )), 1, 3*length(driftNames), byrow=T)); homAll_effects
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
