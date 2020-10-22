
# now requires package lavaan

#######################################################################################################################
############################################ CoTiMA Statistical Power #################################################
#######################################################################################################################

ctmaPower <- function(
  # Primary Study Fits
  ctmaInitFit=NULL,                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,
  #resultsFilePrefix="CoTiMPower",
  #saveFilePrefix="ctmaPower",

  statisticalPower=c(),
  timeRange=NULL,                     # e.g., seq(0, 120, 1)
  useMBESS=FALSE,

  # Fitting Parameters
  #type="stanc",
  coresToUse=1,

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  silentOverwrite=FALSE,

  digits=8,

  saveStatPower="CoTiMAPower_AllResults",
  loadStatPower=NULL,
  loadAllInvFit=c(),
  saveAllInvFit=NULL,
  loadAllInvWOSingFit=c(),
  saveAllInvWOSingFit=NULL,

  skipScaling=TRUE,

  useSampleFraction=NULL, # select subsample (percent of the full sample, e.g. 20 selects every 5th case)

  CoTiMAStanctArgs=list(test=TRUE, scaleTI=TRUE, scaleTime=1/1, scaleMod=TRUE, scaleLongData=FALSE,
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

  { ### CHECKS
    if (is.null(activeDirectory)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("No working directory has been specified!", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat("Press 'q' to quit and change or 'c' to continue with the active directory of the ctmaInitFit file. Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
      activeDirectory <- ctmaInitFit$activeDirectory
    }

    # check if fit object is specified
    if (is.null(ctmaInitFit)){
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("A fitted CoTiMA object (\"ctmaInitFit\") has to be supplied to analyse something. \n"))
      stop("Good luck for the next try!")
    }

    if (length(statisticalPower) < 1) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold(" ", sep="\n"))
      cat(red$bold("Levels (betas) for statistical power (\"statisticalPower\") have not been suppiled \n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat("Press 'q' to quit and change or 'c' to continue with the statisticalPower=c(.90, .80). Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
      statisticalPower <- c(.90, .80)
    }

    if (length(statisticalPower) > 0) {
      for (i in 1:length(statisticalPower)) {
        if (statisticalPower[i] < 0 | statisticalPower[i] > 1){
          if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          cat(red$bold("Parameter for statistical poweranalysis outside the allowed interval!\n"))
          cat(red("Values have to be between 0 and 1\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    #if (resultsFilePrefix=="ctmaPower") {
    #  if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    #  cat("The default results file prefix (ctmaPower) has been chosen.", "\n")
    #  cat("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards ", "\n")
    #  char <- readline(" ")
    #  while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
    #    cat((blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
    #    char <- readline(" ")
    #  }
    #  if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    #}
  } # END Checks


  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
{
  if  (length(coresToUse) > 0) {
    if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
  }

  if (coresToUse >= parallel::detectCores()) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
    coresToUse <- parallel::detectCores() - 1
    cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
  }
}

  #######################################################################################################################
  ############## Extracting Parameters from Fitted Primary Studies created with ctmaINIT Function  ######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {
    n.latent <- ctmaInitFit$n.latent; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- ctmaInitFit$n.studies; n.studies
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    driftNames <- driftNamesBackup <- ctmaInitFit$parameterNames$DRIFT; driftNames
    diffusionNames <- diffusionNamesBackup <- ctmaInitFit$parameterNames$DIFFUSION; diffusionNames
    T0varNames <- T0varNamesBackup <- ctmaInitFit$parameterNames$T0VAR; T0varNames

    # all parameter estimates
    tmp0 <- ctmaInitFit$studyFitList[[1]]$resultsSummary$parmatrices; tmp0
    driftRows <- which(rownames(tmp0) == "DRIFT"); driftRows
    diffusionRows <- which(rownames(tmp0)  == "DIFFUSIONcov"); diffusionRows
    T0varRows <- which(rownames(tmp0) =="T0VAR"); T0varRows
    {
      tmp1 <- lapply(ctmaInitFit$studyFitList, function(extract) c(extract$resultsSummary$parmatrices[driftRows, 3]) ); tmp1
      tmp1 <- matrix(unlist(tmp1), nrow=n.studies, byrow=TRUE); tmp1
      tmp2 <- lapply(ctmaInitFit$studyFitList, function(extract) c(extract$resultsSummary$parmatrices[diffusionRows, 3]) ); tmp2
      tmp2 <- matrix(unlist(tmp2), nrow=n.studies, byrow=TRUE); tmp2
      tmp3 <- lapply(ctmaInitFit$studyFitList, function(extract) c(extract$resultsSummary$parmatrices[T0varRows, 3]) ); tmp3
      tmp3 <- matrix(unlist(tmp3), nrow=n.studies, byrow=TRUE); tmp3

      all_Coeff <- cbind(tmp1, tmp2, tmp3); all_Coeff
      colnames(all_Coeff) <- c(names(ctmaInitFit$modelResults$DRIFT[[1]]),
                               names(ctmaInitFit$modelResults$DIFFUSION[[1]]),
                               names(ctmaInitFit$modelResults$T0VAR[[1]]))
    }
    # all standard errors
    {
      tmp1 <- lapply(ctmaInitFit$studyFitList, function(extract) c(extract$resultsSummary$parmatrices[driftRows, 4]) ); tmp1
      tmp1 <- matrix(unlist(tmp1), nrow=n.studies, byrow=TRUE); tmp1
      tmp2 <- lapply(ctmaInitFit$studyFitList, function(extract) c(extract$resultsSummary$parmatrices[diffusionRows, 4]) ); tmp2
      tmp2 <- matrix(unlist(tmp2), nrow=n.studies, byrow=TRUE); tmp2
      tmp3 <- lapply(ctmaInitFit$studyFitList, function(extract) c(extract$resultsSummary$parmatrices[T0varRows, 4]) ); tmp3
      tmp3 <- matrix(unlist(tmp3), nrow=n.studies, byrow=TRUE); tmp3

      all_SE <- cbind(tmp1, tmp2, tmp3); all_SE
      colnames(all_SE) <- c(names(ctmaInitFit$modelResults$DRIFT[[1]]),
                            names(ctmaInitFit$modelResults$DIFFUSION[[1]]),
                            names(ctmaInitFit$modelResults$T0VAR[[1]]))
      #all_SE
    }

    allSampleSizes <- ctmaInitFit$statisticsList$allSampleSizes; allSampleSizes

    # CHD maxTpointsModel <- which(ctmaInitFit$statisticsList$allTpoints == max(ctmaInitFit$statisticsList$allTpoints)); maxTpointsModel
    maxTpoints <- max(allTpoints); maxTpoints # replacement

    # CHD ctsemModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpointsModel, manifestNames=manifestNames,    # 2 waves in the template only
    ctsemModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                          DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent, byrow=TRUE), # byrow because names are in stanct order
                          LAMBDA=diag(n.latent),
                          type='stanct',
                          #CINT=matrix(cintNames, nrow=n.latent, ncol=1),
                          CINT=matrix(0, nrow=n.latent, ncol=1),
                          T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                          MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                          MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent))


    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas); maxDelta
    if (is.null(timeRange)) usedTimeRange <- seq(0, 1.5*maxDelta, 1) else usedTimeRange <- timeRange
    # augment by all existing time lags
    usedTimeRange <- sort(unique(c(usedTimeRange, allDeltas))); usedTimeRange
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints

    stepWidth <- max(usedTimeRange)/(length(usedTimeRange)-1)
  }


  #######################################################################################################################
  ################################################# data preparation ####################################################
  #######################################################################################################################
  {
    # combine pseudo raw data for model
    tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit)

    if (skipScaling == FALSE) {
      tmp1 <- grep("_T", colnames(tmp$alldata)); tmp1
      tmp$alldata[, tmp1] <- scale(tmp$alldata[, tmp1])
      }

    datawide_all <- tmp$alldata
    groups <- tmp$groups

    # possible subsample selection
    if (!(is.null(useSampleFraction))) {
      N <- dim(datawide_all)[1]; N
      stepwidth <- 100/useSampleFraction
      targetCases <- round(seq(1, N, stepwidth), 0); targetCases
      datawide_all <- datawide_all[targetCases, ]
      groups <- groups[targetCases]
    }

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
    dataTmp2 <- ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.latent, n.TIpred = (n.studies-1),
                             manifestNames=manifestNames)
    dataTmp3 <- ctDeintervalise(dataTmp2)
    dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
    # eliminate rows where ALL latents are NA
    dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
    datalong_all <- dataTmp3
  }


  #######################################################################################################################
  ##################################### Analyses of Statistical Power ###################################################
  #######################################################################################################################

  print(paste0("#################################################################################"))
  print(paste0("########################## Statistical Power Analysis ###########################"))
  print(paste0("#################################################################################"))

  print(paste0("#################################################################################"))
  print(paste0("######## Fitting all fixed CoTiMA - ALL parameters equal across groups ##########"))
  print(paste0("#################################################################################"))

  datalong_all <- datalong_all[, -grep("TI", colnames(datalong_all))]

  # all fixed model is a model with no TI predictors (identical to ctsemModel)
  # CHD allFixedModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpointsModel, manifestNames=manifestNames,    # 2 waves in the template only
  allFixedModel <- ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                           DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent, byrow=TRUE), # byrow because names are in stanct order
                           LAMBDA=diag(n.latent),
                           type='stanct',
                           #CINT=matrix(cintNames, nrow=n.latent, ncol=1),
                           CINT=matrix(0, nrow=n.latent, ncol=1),
                           T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                           MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                           MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent))
  # LOAD or Fit
  if (length(loadAllInvFit) > 0) {
    x1 <- paste0(activeDirectory, loadAllInvFit[1], ".rds"); x1
    results <- readRDS(file=x1)
  } else {
    results <- ctStanFit(
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
    resultsSummary <- summary(results, digits=8)

  }

  # SAVE
  if (length(saveAllInvFit) > 0)  {
    x1 <- paste0(saveAllInvFit[1], ".rds"); x1
    x2 <- paste0(activeDirectory); x2
    CoTiMASaveFile(activateRPB, "", results, x1, x2, silentOverwrite=silentOverwrite)
  }

  ### Extract estimates & statistics
  tmp <- grep("toV", rownames(resultsSummary$popmeans)); tmp
  homAll_Drift_Coef <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"], n.latent, byrow=TRUE)); homAll_Drift_Coef
  names(homAll_Drift_Coef) <- rownames(resultsSummary$popmeans)[tmp]; homAll_Drift_Coef
  homAll_Drift_SE <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"], n.latent, byrow=TRUE)); homAll_Drift_SE
  names(homAll_Drift_SE) <- rownames(resultsSummary$popmeans)[tmp]; homAll_Drift_SE
  tmp1 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=TRUE)); tmp1
  tmp2 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=TRUE)); tmp2
  homAll_Drift_CI <- c(rbind(tmp1, tmp2)); homAll_Drift_CI
  tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                  paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
  names(homAll_Drift_CI) <- tmp3; homAll_Drift_CI
  homAll_Drift_Tvalue <- homAll_Drift_Coef/homAll_Drift_SE; homAll_Drift_Tvalue

  tmp <- grep("diff", rownames(resultsSummary$popmeans)); tmp
  homAll_Diffusion_Coef <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"]); homAll_Diffusion_Coef
  names(homAll_Diffusion_Coef) <- rownames(resultsSummary$popmeans)[tmp]; homAll_Diffusion_Coef
  homAll_Diffusion_SE <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Sd"]); homAll_Diffusion_SE
  names(homAll_Diffusion_SE) <- rownames(resultsSummary$popmeans)[tmp]; homAll_Diffusion_SE
  tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "2.5%"]; tmp1
  tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "97.5%"]; tmp2
  homAll_Diffusion_CI <- c(rbind(tmp1, tmp2)); homAll_Diffusion_CI
  tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                  paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
  names(homAll_Diffusion_CI) <- tmp3; homAll_Diffusion_CI
  homAll_Diffusion_Tvalue <- homAll_Diffusion_Coef/homAll_Diffusion_SE; homAll_Diffusion_Tvalue

  tmp <- grep("T0var", rownames(resultsSummary$popmeans)); tmp
  homAll_T0Var_Coef <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Mean"]); homAll_T0Var_Coef
  #homAll_T0Var_Coef <- cov2cor(vech2full(homAll_T0Var_Coef)); homAll_T0Var_Coef
  #homAll_T0Var_Coef <- homAll_T0Var_Coef[lower.tri(homAll_T0Var_Coef, diag=TRUE)]; homAll_T0Var_Coef
  names(homAll_T0Var_Coef) <- rownames(resultsSummary$popmeans)[tmp]; homAll_T0Var_Coef
  homAll_T0Var_SE <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Sd"]); homAll_T0Var_SE
  names(homAll_T0Var_SE) <- rownames(resultsSummary$popmeans)[tmp]; homAll_T0Var_SE
  tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]; tmp1
  tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "97.5%"]; tmp2
  homAll_T0Var_CI <- c(rbind(tmp1, tmp2)); homAll_T0Var_CI
  tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                  paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
  names(homAll_T0Var_CI) <- tmp3; homAll_T0Var_CI
  homAll_T0Var_Tvalue <- homAll_T0Var_Coef/homAll_T0Var_SE; homAll_T0Var_Tvalue


  ## Extract Model Fit
  homAll_Minus2LogLikelihood <- 2* results$stanfit$optimfit$f; homAll_Minus2LogLikelihood
  homAll_estimatedParameters <-length(results$stanfit$optimfit$par); homAll_estimatedParameters
  #homAll_df <- ctmaInitFit$summary$df+(ctmaInitFit$summary$n.parameters-homAll_estimatedParameters); homAll_df
  homAll_df <- NULL

  # Combine summary information
  homAll_effects <- matrix(t(cbind((homAll_Drift_Coef), (homAll_Drift_SE),
                                   (homAll_Drift_Tvalue))), 1, 3*length(driftNames), byrow=T); homAll_effects
  homAll_effects <- rbind(homAll_effects,
                          matrix(t(cbind((c(vech2full(homAll_Diffusion_Coef))),
                                         c(vech2full((homAll_Diffusion_SE))),
                                         c(vech2full((homAll_Diffusion_Tvalue))) )), 1, 3*length(driftNames), byrow=T)); homAll_effects
  homAll_effects <- rbind(homAll_effects,
                          matrix(t(cbind(c(vech2full((homAll_T0Var_Coef))),
                                         c(vech2full((homAll_T0Var_SE))),
                                         c(vech2full((homAll_T0Var_Tvalue))) )), 1, 3*length(driftNames), byrow=T)); homAll_effects
  # Label summary table
  rownames(homAll_effects) <- c("Fixed Effects Drift", "Fixed Effects Diffusion", "Fixed Effects T0Var")
  newColNames <- c()
  for (j in 1:n.latent) {
    for (h in 1:n.latent) {
      newColNames <- c(newColNames, paste0("V",h,"toV", j), "(SE)", "Tvalue")
    }
  }
  colnames(homAll_effects) <- newColNames; homAll_effects

  DRIFT <- matrix(homAll_Drift_Coef, n.latent, n.latent, byrow=TRUE); DRIFT
  DIFFUSION <- vech2full(homAll_Diffusion_Coef); DIFFUSION
  T0VAR <- vech2full(homAll_T0Var_Coef); T0VAR

  print(paste0("#################################################################################"))
  print(paste0("######## Computing implied covariance matrices for different time lags ##########"))
  print(paste0("#################################################################################"))

  {
    # functions to compute dt-coefficients
    discreteDriftFunction <- function(driftMatrix, timeScale, number) {
      discreteDriftValue <- expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }
    discreteDiffusionFunction <- function(diffusionMatrix, driftMatrix, timeScale, number) {
      driftHatch <- driftMatrix %x% diag(dim(diffusionMatrix)[1]) + diag(dim(diffusionMatrix)[1]) %x% driftMatrix
      discreteDiffusionValue <- solve(driftHatch) %*% (expm(timeScale %x% driftHatch) - diag(dim(driftHatch)[1])) %*% c(diffusionMatrix)
      discreteDiffusionValue[number] }

    implCov <- list()

    for (t in 1:length(usedTimeRange)) {
      # compute the covariance matrix implied by the dt-coefficients across requested time intervals
      # for equations see Neudecker, H. & Satorra, A. (1991). Linear structural relations: Gradient and ...
      # ... Hessian of the fitting function. Statistics & Probability Letters 11 (1991) 57-61. North-Holland
      beta <- discreteDriftFunction(DRIFT, t); beta
      tmpMat <- matrix(0, n.latent, n.latent); tmpMat
      B <- diag(1, 2*n.latent) - cbind(rbind(tmpMat, beta), rbind(tmpMat, tmpMat)); B
      phi <- T0VAR; phi
      psi <- matrix(discreteDiffusionFunction(DIFFUSION, DRIFT, t, 1:4), n.latent); psi
      phi <- cbind(rbind(phi, tmpMat), rbind(tmpMat, tmpMat)); #phi
      psi <- cbind(rbind(tmpMat, tmpMat), rbind(tmpMat, psi)); #psi
      tmp <- solve(B) %*% phi %*% t(solve(B)) + psi; #tmp
      implCov[[t]] <- cov2cor(tmp); implCov[[t]]
      rownames(implCov[[t]]) <- varNames
      colnames(implCov[[t]]) <- varNames
    }

  }

  print(paste0("#################################################################################"))
  print(paste0("################# Set up required discrete time lavaan models ###################"))
  print(paste0("#################################################################################"))

  {

    # full lavaan model setup
    {
      modelText <- c()
      counter <- 0
      for (i in 1:n.latent) {
        counter <- counter + 1
        modelText[counter] <- ""
        for (j in 1:n.latent) {
          if (j == 1) modelText[counter] <- paste0(modelText[counter], paste0("V", i, "T1 ~ V", j, "T0"))
          if (j != 1) modelText[counter] <- paste0(modelText[counter], paste0(" + V", j, "T0"))
        }
      }
      counter <- n.latent; counter
      for (i in 1:n.latent) {
        for (j in i:n.latent) {
          if (i != j) {
            counter <- counter + 1; counter
            modelText[counter] <- paste0("V", i, "T0 ~~ V", j, "T0")
            counter <- counter + 1; counter
            modelText[counter] <- paste0("V", i, "T1 ~~ V", j, "T1")
          }
        }
      }
    }
    tmp1 <- paste(modelText, sep="", collapse="\n "); tmp1
    model.full <- paste0("\n ", tmp1, "\n"); tmp1

    # lavaan model setup without single cross effects
    model.wo <- list()
    counter <- 0
    for (i in 1:n.latent) {
      # split text of full model
      tmp1 <- unlist(strsplit(model.full, "\n")); tmp1
      tmp2 <- tmp1[tmp1 != ""]; tmp2
      for (j in 1:n.latent) {
        #j <- 2
        if (i != j) {
          counter <- counter + 1
          if (i == 1) {
            toDelete <- paste0(" \\+", " V", j, "T0"); toDelete
          } else {
            toDelete <- paste0("V", j, "T0 ", "\\+"); toDelete
          }
          tmp2[i] <- gsub(toDelete, "", tmp2[i]); tmp2[i]
          tmp3 <- paste(tmp2, sep="", collapse="\n "); tmp3
          model.wo[[counter]] <- paste0("\n ", tmp3, "\n"); model.wo[[counter]]
        }
      }
    }
    #model.wo

  }

  print(paste0("#################################################################################"))
  print(paste0("########### Compute required sample size to achieve requested power #############"))
  print(paste0("#################################################################################"))


  # Fast function to calculate required sample sizes later (as optional replacement for ss.power.reg.coef)
  nestedProbFunT <- function (fvalue, alpha=.05, power=.80, p=2, x) (1-
                                                                       pt(
                                                                         qt((1 - alpha/2), df = (x)-p-1,
                                                                            lower.tail = TRUE, log.p = FALSE),
                                                                         df = (x)-p-1, ncp = sqrt(x) * abs(fvalue),
                                                                         lower.tail = TRUE, log.p = FALSE)) - power

  # Create table: sampleSizes x deltas (of primary studies) for post hoc power calculations
  tableNxDeltas <- matrix(NA, nrow=n.studies, ncol=maxTpoints); tableNxDeltas
  tableNxDeltas[ ,1]  <- unlist(allSampleSizes); tableNxDeltas
  counter <- 0
  for (j in 1:n.studies) {
    for (h in 1:(allTpoints[j]-1)) {
      counter <- counter + 1; counter
      tableNxDeltas[j , (1+h)] <- allDeltas[[counter]]
    }
  }
  tableNxDeltas[is.na(tableNxDeltas)] <- -99; tableNxDeltas
  tableNxPowerAlpha05 <- tableNxDeltas; tableNxPowerAlpha05
  tableNxPowerAlpha05[ , 2:maxTpoints] <- NA; tableNxPowerAlpha05
  tableNxPowerAlpha01 <- tableNxPowerAlpha05; tableNxPowerAlpha01
  listPowerAlpha05 <- list()
  listPowerAlpha01 <- list()

  effectSizes <- matrix(NA, nrow=length(usedTimeRange), ncol=(n.latent^2-n.latent))

  # Loop through a range of lags to determine sample sizes (same parameters as for plotting the effects furter below)
  plotPairs <- array(dim=c(n.latent^2, length(statisticalPower), length(usedTimeRange), 2))  # all drift effects, all powers, time range, timePoint+SampleSize

  # Plot required sample size for cross effects.
  for (h in 1:length(statisticalPower)) {
    counter <- 0
    for (j1 in 1:(n.latent)) {
      for (j2 in 1:(n.latent)) {
        #counter <- counter + 1; counter
        if (j1 != j2) {
          counter <- counter + 1; counter
          #j <- counter; j
          for (k in 1:(length(usedTimeRange)-1)) {
            delta_t <- usedTimeRange[k+1]; delta_t
            #plotPairs[j, h, k, 1] <- usedTimeRange[k+1] # time point
            plotPairs[counter, h, k, 1] <- usedTimeRange[k+1] # time point

            # R2 in terms of Kelley & Maxwell 2008
            #A <- expm(DRIFT %x% delta_t) %*% T0VAR %*% t(expm(DRIFT %x% delta_t)); A         # implied variance at later Tpoint
            #S <- matrix(discreteDiffusionFunction(DIFFUSION, DRIFT, delta_t, 1:(n.latent^2)),
            #            n.latent, n.latent); S   # residual variance at later Tpoint
            #R2 <- A[j1,j1]/( (A + S)[j1,j1] ); R2                                            # explained variance at later Tpoint
            # betas & psis for model with all effects included
            sample.nobs <- 1000 # large enough to prevent shrinkage
            model.full.fit <- lavaan::sem(model.full, sample.cov = implCov[[k]], sample.nobs = sample.nobs)
            tmp <- lavaan::inspect(model.full.fit, "est"); tmp
            psi <- diag(tmp$psi)[grep("T1", names(diag(tmp$psi)))]; psi
            beta <- tmp$beta[1:n.latent ,(n.latent+1):(2*n.latent)]; beta
            R2 <- 1 - psi[j1]; R2
            #j1; j2
            effectSizes[k, counter] <- beta[j2, j1];

            # R2 without j (cross effect) in terms of Kelley & Maxwell 2008
            #A.j <- expm(DRIFT.without.j[[counter]] %x% delta_t) %*%
            #  T0VAR.without.j[[j]] %*%
            #  t(expm(DRIFT.without.j[[j]] %x% delta_t)); A.j           # implied variance without j (counter) as predictor at later Tpoint
            #S.j <- matrix(
            #  discreteDiffusionFunction(DIFFUSION.without.j[[j]],      # implied residual variance without j (counter) as predictor at later Tpoint
            #                            DRIFT.without.j[[j]],
            #                            delta_t, 1:(n.latent^2)), n.latent,n.latent); S.j
            #
            #R2.j <- A.j[j1,j1]/( (A.j + S.j)[j1,j1] ); R2.j           # explained variance without j (counter) as predictor at later Tpoint
            model.wo.fit <- lavaan::sem(unlist(model.wo[[counter]]),
                                sample.cov = implCov[[k]],
                                sample.nobs = sample.nobs)
            #tmp <- inspect(model.wo.fit[[counter]], "est"); tmp
            tmp <- lavaan::inspect(model.wo.fit, "est"); tmp
            R2.j <- 1 - tmp$psi[j1,j1]; R2.j

            # Rounding to prevent endless computations
            #powerRounding <- 8
            #R2 <- round(R2, powerRounding); R2
            #R2.j <- round(R2.j, powerRounding); R2.j

            # Skip in case that more variance is explained after predictor is removed
            if (R2 < R2.j) {
              plotPairs[counter, h, k, 2] <- 1000000
            } else {
              # The following is Kelley's function, which is replaced below by our own, which is > 50 times faster
              # Kelley, K. (2019). The MBESS R Package. R package version 4.6.0. Retrieved from:
              # https://cran.r-project.org/web/packages/MBESS/MBESS.pdf
              if (useMBESS == TRUE) {
                plotPairs[counter, h, k, 2] <- MBESS::ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                                           p = n.latent, desired.power = statisticalPower[h],
                                                           alpha.level = 0.05)[[1]] #
              } else {
                # The following uses our own function
                signalToNoiseRatios <- sqrt((R2-R2.j)/(1-R2)); signalToNoiseRatios
                helper <- round(rootSolve::uniroot.all(nestedProbFunT, c(n.latent+2,999999999),
                                            fvalue=signalToNoiseRatios, alpha=.05,
                                            power=statisticalPower[h], p=n.latent) + .49999, 0)
                if (length(helper) < 1) helper <- NA
                plotPairs[counter, h, k, 2] <- helper
              }

              # Post hoc power computations
              if ( (delta_t %in% tableNxDeltas[ ,-1]) & (h == 1) ){  # do only once (not for all a priori powers)
                M <- tableNxDeltas[ ,-1] == delta_t; M # temp matrix used below
                empiricalN <- matrix(tableNxDeltas[apply(M, 1, any), ], ncol=maxTpoints)[,1]; empiricalN
                empiricalN <- na.omit(empiricalN); empiricalN
                for (l in empiricalN) {
                  p05 <- MBESS::ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                           p = n.latent, Specified.N = l, alpha.level = 0.05)[2]
                  p01 <- MBESS::ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                           p = n.latent, Specified.N = l, alpha.level = 0.01)[2]
                  for (m in 1:n.studies) { # poke power into matrices
                    if (tableNxDeltas[m,1] == l) { # if study has current empirical N ...
                      for (n in 2:(maxTpoints)) {  # ... loop through al lags
                        if (tableNxDeltas[m,n] == delta_t) { # if lag corresponds to current lag ...
                          tableNxPowerAlpha05[m,n] <- as.numeric(p05)  # ... poke
                          tableNxPowerAlpha01[m,n] <- as.numeric(p01)  # ... poke
                        }
                      } # (poke power into matrices) end maxTpoints loops
                    } # (poke power into matrices) end if tableNxDeltas[m,1] == l (delta)
                  } # (poke power into matrices) end n.studies loop
                } # (compute power if necessary) end empiricalN loop
              } # (compute power if necessary) end delta_t %in% tableNxDeltas[-1,]
            } # (compute power if possible) ende else
          } # end usedTimeRange loop
          if (h == 1) {
            listPowerAlpha05[[counter]] <- tableNxPowerAlpha05
            listPowerAlpha01[[counter]] <- tableNxPowerAlpha01
          }
        } # end j1!=j2 condition
      } # end j2 loop
    } # end j1 loop
    print(paste0("#################################################################################"))
    print(paste0("###################### power calculation for ", statisticalPower[h], " completed ######################"))
    print(paste0("#################################################################################"))
  } # end h loop (length(statisticalPower))

  # Table of required sample sizes for range of different time lags (a priori power)
  requiredSampleSizes <- list()
  currentDriftNames <- c()
  counter1 <- counter2 <- 0
  for (j1 in 1:(n.latent)) {
    for (j2 in 1:(n.latent)) {
      counter1 <- counter1 + 1 # does not seem to be necessary
      if (j1 != j2 ) {
        counter2 <- counter2 + 1
        #requiredSampleSizes[[counter2]] <- plotPairs[counter1, , , 2]
        requiredSampleSizes[[counter2]] <- plotPairs[counter2, , , 2]
        currentDriftNames <- c(currentDriftNames, driftNamesBackup[counter1])
        #rowNames  <- plotPairs[counter1, 1, , 1]
        rowNames  <- plotPairs[counter2, 1, , 1]
      }
    }
  }

  # re-structure into a single table and replace 100000 by NA
  numberOfEffects <- n.latent^2 - n.latent; numberOfEffects
  tmp <- requiredSampleSizes[[1]]; tmp
  for (h in 2:(numberOfEffects)) tmp <- rbind(tmp, requiredSampleSizes[[h]])
  requiredSampleSizes <- t(tmp); requiredSampleSizes
  requiredSampleSizes[requiredSampleSizes==1000000] <- NA

  # Label columns
  columnNames <- c()
  for (h in 1:numberOfEffects) {
    for (j in 1:length(statisticalPower)) {
      columnNames <- c(columnNames, paste0(currentDriftNames[h], " Power=", statisticalPower[j]))
    }
  }
  colnames(requiredSampleSizes) <- columnNames
  rownames(requiredSampleSizes) <- round(rowNames, digits)
  #requiredSampleSizes

  # add (not really standardized) effect sizes based on matrix exponentiation
  tmp1 <- as.numeric(rownames(requiredSampleSizes))
  timeLags <- tmp1[!(is.na(tmp1))]; timeLags
  effectSizes2 <- matrix(NA, nrow=length(usedTimeRange), ncol=(n.latent*(n.latent-1))); effectSizes2
  effectCounter <- 0
  for (i in 1:n.latent) {
    for (j in 1:n.latent) {
      if (i != j) {
        effectCounter <- effectCounter + 1
        rowCounter <- 0
        for (k in 1:(length(usedTimeRange)-1)) {
          rowCounter <- rowCounter + 1; rowCounter
          A <- expm(DRIFT %x% usedTimeRange[k+1])
          effectSizes2[rowCounter, effectCounter] <- round(A[i,j], digits)
        }
      }
    }
  }

  diffDim <- dim(requiredSampleSizes)[1] - dim(effectSizes)[1]; diffDim
  helperMat <- matrix(NA, nrow=diffDim, ncol=dim(effectSizes)[2]); helperMat
  effectSizes <- rbind(effectSizes, helperMat)
  tmp1 <- matrix(driftNames, n.latent, n.latent)
  diag(tmp1) <- NA
  tmp1 <- tmp1[!(is.na(tmp1))]; tmp1
  colnames(effectSizes) <- tmp1
  requiredSampleSizes <- cbind(requiredSampleSizes, effectSizes)

  # Determine optimal time lag in terms of min sample size required
  rowNames <- c(rownames(requiredSampleSizes), "Min N", "Opt. Lag"); rowNames
  minN <- (apply(requiredSampleSizes, 2, min, na.rm=TRUE)); minN
  optimalCrossLagN <- c()
  for (h in 1:(dim(requiredSampleSizes)[2])) optimalCrossLagN[h] <- mean(which(requiredSampleSizes[ ,h ] == minN[h]))
  optimalCrossLagN <- optimalCrossLagN*stepWidth; optimalCrossLagN
  requiredSampleSizes <- rbind(requiredSampleSizes, minN, optimalCrossLagN)
  rownames(requiredSampleSizes) <- rowNames

  # Formatting of post hoc results
  postHocPowerList <- list()
  tableNxDeltas[tableNxDeltas == -99] <- NA; tableNxDeltas
  columnNames <- c("N", rep(c("Time Lag", "Power (\u03b1=.05)", "Power (\u03b1=.01)"), (maxTpoints-1))); columnNames

  # Remove empty elements form list and combine
  listPowerAlpha05 <- listPowerAlpha05[!sapply(listPowerAlpha05, is.null)]
  listPowerAlpha01 <- listPowerAlpha01[!sapply(listPowerAlpha01, is.null)]
  for (j in 1:length(currentDriftNames)) {
    tmp05 <- round(listPowerAlpha05[[j]], digits)
    tmp01 <- round(listPowerAlpha01[[j]], digits)
    postHocPower <- tableNxDeltas[, 1]; postHocPower
    for (k in 1:(maxTpoints-1)) { postHocPower <- cbind(postHocPower,
                                                        tableNxDeltas[,k+1],
                                                        tmp05[,k+1],
                                                        tmp01[,k+1])
    }
    colnames(postHocPower) <- columnNames; postHocPower

    # Compute mean & median power (for first lag only) with and without temporarily replacing NA with 0
    tmp <- postHocPower
    tmp[is.na(tmp)] <- 0; tmp
    targetCols <- grep("Power", columnNames); targetCols
    meanPower0 <- apply(tmp[, targetCols[1:2]], 2, mean, na.rm=TRUE); meanPower0
    meanPowerNA <- apply(postHocPower[, targetCols[1:2]], 2, mean, na.rm=TRUE); meanPowerNA
    medianPower0 <- apply(tmp[, targetCols[1:2]], 2, median, na.rm=TRUE); medianPower0
    medianPowerNA <- apply(postHocPower[, targetCols[1:2]], 2, median, na.rm=TRUE); medianPowerNA
    postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
    postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
    postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
    postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower # four times is correct
    postHocPower[dim(postHocPower)[1]-3, targetCols[1:2]] <- round(meanPower0, digits)
    postHocPower[dim(postHocPower)[1]-2, targetCols[1:2]] <- round(meanPowerNA, digits)
    postHocPower[dim(postHocPower)[1]-1, targetCols[1:2]] <- round(medianPower0, digits)
    postHocPower[dim(postHocPower)[1], targetCols[1:2]] <- round(medianPowerNA, digits)
    newNames <- c(paste0("Study_No_", 1:n.studies),
                  "Mean (NA = 0 Power)", "Mean (NA = NA)",
                  "Median (NA = 0 Power)", "Median (NA = NA)")
    rownames(postHocPower) <- newNames
    postHocPowerList[[j]] <- postHocPower
    names(postHocPowerList)[[j]] <- currentDriftNames[j]
  }

  results <- list(activeDirectory=activeDirectory,
                  plot.type=c("power"), model.type="stanct", #model.type="mx",
                  coresToUse=NULL, n.studies=1,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, #studyFitList=list(homAllFit), #fullWOSingleFit)
                  emprawList=NULL,
                  statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=DRIFT, DIFFUSION=DIFFUSION, T0VAR=T0VAR, CINT=NULL,
                                    DRIFTwoJ=DRIFT.without.j, DIFFUSIONwoJ=DIFFUSION.without.j,
                                    T0VARwoJ=T0VAR.without.j),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="Analysis of Statistical Power and Required Sample Sizes",
                               estimates=list("Estimates of Model with all Effects Invariant"=round(homAll_effects, digits),
                                              "Requested Statistical Power"=statisticalPower,
                                              "Power (post hoc) for Drift Effects"=postHocPowerList,
                                              "Required Sample Sizes"=round(requiredSampleSizes, digits))))
  class(results) <- "CoTiMAFit"

  invisible(results)

} ### END function definition

