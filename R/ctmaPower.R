#' ctmaPower
#'
#' @description Fits a full invariant model to a list of primary studies and performs analyses of expected (post hoc) power and required sample sizes.
#'
#' @param ctmaInitFit object to which all single 'ctsem' fits of primary studies has been assigned to (i.e., what has been returned by \code{\link{ctmaInit}})
#' @param activeDirectory defines another active directory than the one used in \code{\link{ctmaInit}}
#' @param statisticalPower vector of requested statistical power values
#' @param failSafeN  sample size used to determine across which time intervals effects become non-significant
#' @param failSafeP  p-value used to determine across which time intervals effects become non-significant
#' @param timeRange vector describing the time range for x-axis as sequence from/to/stepSize (e.g., c(1, 144, 1))
#' @param useMBESS use 'MBESS' package to calculate statistical power (slower)
#' @param coresToUse if negative, the value is subtracted from available cores, else value = cores to use
#' @param indVarying Allows continuous time intercepts to vary at the individual level (random effects model, accounts for unobserved heterogeneity)
#' @param activateRPB set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param silentOverwrite overwrite old files without asking
#' @param digits number of digits used for rounding (in outputs)
#' @param loadAllInvFit load the fit of fully constrained 'CoTiMA' model
#' @param saveAllInvFit save the fit of fully constrained 'CoTiMA' model
#' @param loadAllInvWOSingFit load series of fits of fully constrained 'CoTiMA' model with single cross effects excluded, respectively
#' @param saveAllInvWOSingFit save series of fits of fully constrained 'CoTiMA' model with single cross effects excluded, respectively
#' @param skipScaling does not (re-)scale raw data (re-scaling of imported pseudo raw data achieves correlations = 1)
#' @param useSampleFraction to speed up debugging. Provided as fraction (e.g., 1/10)
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param nopriors if TRUE, any priors are disabled – sometimes desirable for optimization
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param iter number of iterations (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#' @param customPar logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting

#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red blue
#' @importFrom parallel detectCores
#' @importFrom ctsem ctModel ctWideToLong ctDeintervalise
#' @importFrom OpenMx vech2full expm
#' @importFrom stats cov2cor pt qt na.omit median
#' @importFrom lavaan sem inspect
#' @importFrom MBESS ss.power.reg.coef
#' @importFrom rootSolve uniroot.all
#'
#' @examples \dontrun{
#' CoTiMAInitFit_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAPower_D_BO <- ctmaPower(ctmaInitFit=CoTiMAInitFit_D_BO,
#'                               statisticalPower = c(.50, .80, .95),
#'                               finishsamples = 10000)
#' summary(CoTiMAPower_D_BO)
#' }
#'
#' @export ctmaPower
#'
#' @return ctmaPower returns a list containing some arguments supplied, a fitted model with all (!) parameters invariant across primary
#' studies, different elements summarizing the main results, model type, and the type of plot that could be performed with the returned object.
#' The arguments in the returned object are activeDirectory, coresToUse, n.latent, n.manifest, and primaryStudyList. A further result
#' returned is n.studies = 1 (required for proper plotting). Further arguments, which are just copied from the init-fit object supplied, are,
#' n.latent, studyList, and the statisticsList. The fitted model is found in studyFitList, which is a large list with many elements (e.g., the
#' ctsem model specified by CoTiMA, the rstan model created by ctsem, the fitted rstan model etc.). Further results returned are  a list with
#' modelResults (i.e., DRIFT=DRIFT, DIFFUSION=DIFFUSION, T0VAR=T0VAR, CINT=NULL) and the paramter names internally used. The summary list,
#' which is printed if the summary function is applied to the returned object, contains "estimates", which is itself a list comprising
#' "Estimates of Model with all Effects Invariant", "Requested Statistical Power" (which just returns the argument statisticalPower),
#' "Power (post hoc) for Drift Effects",  "Required Sample Sizes" "Effect Sizes (based on discrete-time calcs; used for power calcs.)", and
#' "Range of significant effects" (across which intervals effects were significant). Plot type is plot.type=c("power") and model.type="stanct"
#' ("omx" was deprecated).

#'
ctmaPower <- function(
  ctmaInitFit=NULL,
  activeDirectory=NULL,
  statisticalPower=c(),
  failSafeN =NULL,
  failSafeP=NULL,
  timeRange=NULL,
  useMBESS=FALSE,
  coresToUse=1,
  digits=4,
  indVarying=FALSE,
  activateRPB=FALSE,
  silentOverwrite=FALSE,
  loadAllInvFit=c(),
  saveAllInvFit=c(),
  loadAllInvWOSingFit=c(),
  saveAllInvWOSingFit=c(),
  skipScaling=TRUE,
  useSampleFraction=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  finishsamples=NULL,
  iter=NULL,
  chains=NULL,
  verbose=NULL,
  customPar=FALSE,
  scaleTime=NULL
)

{  # begin function definition (until end of file)

  { ### CHECKS

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

    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory

    # check if fit object is specified
    if (is.null(ctmaInitFit)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA object (\"ctmaInitFit\") has to be supplied to analyse something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (length(statisticalPower) < 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold(" ", sep="\n"))
      cat(crayon::red$bold("Levels (betas) for statistical power (\"statisticalPower\") have not been suppiled \n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat("Press 'q' to quit and change or 'c' to continue with the statisticalPower=c(.90, .80). Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      ErrorMsg <- "Good luck for the next try!"
      if (char == 'q' | char == 'Q') stop(ErrorMsg)
      statisticalPower <- c(.90, .80)
    }
    if (length(statisticalPower) > 0) {
      for (i in 1:length(statisticalPower)) {
        if (statisticalPower[i] < 0 | statisticalPower[i] > 1){
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          ErrorMsg <-"\nParameter for statistical poweranalysis outside the allowed interval! \nValues have to be between 0 and 1 \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
    }

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
      Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash. \n"
      message(Msg)
    }
  }

  #######################################################################################################################
  ############## Extracting Parameters from Fitted Primary Studies created with ctmaINIT Function  ######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {
    n.latent <- ctmaInitFit$n.latent; n.latent
    n.manifest <- ctmaInitFit$n.manifest; n.manifest
    if (is.null(n.manifest)) n.manifest <- 0
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- ctmaInitFit$n.studies; n.studies
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    driftNames <- driftNamesBackup <- ctmaInitFit$parameterNames$DRIFT; driftNames
    # check if drift effects were fixed to 0
    tmp1 <- rownames(ctmaInitFit$studyFitList[[1]]$resultsSummary$popmeans); tmp1
    tmp2 <- which(!(driftNames %in% tmp1)); tmp2
    if (length(tmp2) != 0) driftNames[tmp2] <- "0"

    driftNames <- c(matrix(driftNames, n.latent, byrow=FALSE));     driftNames

    diffusionNames <- diffusionNamesBackup <- ctmaInitFit$parameterNames$DIFFUSION; diffusionNames
    T0varNames <- T0varNamesBackup <- ctmaInitFit$parameterNames$T0VAR; T0varNames

    # all parameter estimates
    tmp0 <- ctmaInitFit$studyFitList[[1]]$resultsSummary$parmatrices; tmp0
    driftRows <- which(rownames(tmp0) == "DRIFT"); driftRows
    if (length(driftRows) == 0) driftRows <- which(tmp0[,1] =="DRIFT"); driftRows

    diffusionRows <- which(rownames(tmp0)  == "DIFFUSIONcov"); diffusionRows
    if (length(diffusionRows) == 0) diffusionRows <- which(tmp0[,1] =="DIFFUSIONcov"); diffusionRows

    T0varRows <- which(rownames(tmp0) =="T0VAR"); T0varRows
    if (length(T0varRows) == 0) T0varRows <- which(rownames(tmp0) =="T0cov"); T0varRows
    if (length(T0varRows) == 0) T0varRows <- which(tmp0[,1] =="T0cov"); T0varRows

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
    #failSafeN
    failSafeNhelper <- ""
    if ( (!(is.null(failSafeN)) & ((is.null(failSafeP))))
         | ((is.null(failSafeN)) & (!(is.null(failSafeP)))) )  {
      if (activateRPB == TRUE) {
        RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
      }
      cat(crayon::red$bold("Only one argument (failSafeN OR failSafeP) has been provided, but both are required!", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      #round(mean(allSampleSizes+.5),0)
      # Changed 26 Nov 2022
      #failSafeNhelper <- ""
      failSafeNhelper <-  round(mean(allSampleSizes+.5),0)
      if (is.null(failSafeN)) {
        failSafeN <- round(mean(allSampleSizes+.5),0)
        failSafeNhelper <- "( = avg. N)"
      }
      if (is.null(failSafeP)) failSafeP <- .01
      cat(crayon::red$bold("I will use failSafeN =", failSafeN, " and failSafeP = ", failSafeP, "!", sep=" "))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat("Press 'q' to quit and change or 'c' to continue with the suggested values. Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      ErrorMsg <- "Good luck for the next try!"
      if (char == 'q' | char == 'Q') stop(ErrorMsg)
    }
    if (is.null(failSafeN)) {
      failSafeN <- round(mean(allSampleSizes+.5),0)
      failSafeNhelper <- "( = avg. N)"
    }
    if (is.null(failSafeP)) failSafeP <- .01

    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints # replacement
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas); maxDelta
    #if (is.null(timeRange)) usedTimeRange <- seq(0, 1.5*maxDelta, 1) else usedTimeRange <- timeRange
    if (is.null(timeRange)) usedTimeRange <- seq(0, 1.5*maxDelta, 1) else usedTimeRange <- seq(timeRange[1], timeRange[2], timeRange[3])
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
    # CHD adeed 26 Nov 2022
    if (length(loadAllInvFit) < 1) {
      # combine pseudo raw data for model
      tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit)
      #var(tmp$alldata[,1:8], na.rm=TRUE)

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
      #var(datawide_all[,1:8], na.rm=TRUE)


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
      dataTmp2 <- suppressMessages(ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.latent, n.TIpred = (n.studies-1),
                                                manifestNames=manifestNames))
      dataTmp3 <- suppressMessages(ctDeintervalise(dataTmp2))
      dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
      # eliminate rows where ALL latents are NA
      dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
      datalong_all <- dataTmp3
      datalong_all <- datalong_all[, -grep("TI", colnames(datalong_all))]
    }
  }


  #######################################################################################################################
  ##################################### Analyses of Statistical Power ###################################################
  #######################################################################################################################

  Msg <- "################################################################################# \n########################## Statistical Power Analysis ########################### \n#################################################################################"
  message(Msg)

  Msg <- "################################################################################# \n######## Fitting all fixed CoTiMA - ALL parameters equal across groups ########## \n#################################################################################"
  message(Msg)

  # CHD moved up 26 Nov 2022
  #datalong_all <- datalong_all[, -grep("TI", colnames(datalong_all))]


  # LOAD or Fit
  if (length(loadAllInvFit) > 0) {
    x1 <- paste0(activeDirectory, loadAllInvFit[1], ".rds"); x1
    # CHD changed Nov 2022
    #results <- readRDS(file=x1)
    allInvModelFit <- readRDS(file=x1)
    allInvModelFitSummary <- summary(allInvModelFit, digits=digits)
  } else {
    allInvModelFit <- ctmaAllInvFit(ctmaInitFit=ctmaInitFit,
                                    activeDirectory=activeDirectory,
                                    activateRPB=activateRPB,
                                    digits=digits,
                                    coresToUse=coresToUse,
                                    indVarying=indVarying,
                                    scaleTime=CoTiMAStanctArgs$scaleTime,
                                    optimize=optimize,
                                    nopriors=nopriors,
                                    finishsamples=finishsamples,
                                    iter=iter,
                                    chains=chains,
                                    verbose=verbose,
                                    customPar=customPar)
    Msg <- "\nComputing results summary of all invariant model. \n"
    message(Msg)
    allInvModelFitSummary <- summary(allInvModelFit, digits=digits)
  }

  # SAVE
  if (length(saveAllInvFit) > 0)  {
    x1 <- paste0(saveAllInvFit[1], ".rds"); x1
    x2 <- paste0(activeDirectory); x2
    ctmaSaveFile(activateRPB, "", allInvModelFit, x1, x2, silentOverwrite=silentOverwrite)
  }

  ### Extract estimates & statistics
  # account for changes in ctsem 3.4.1
  {
    if ("matrix" %in% colnames(allInvModelFitSummary$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
    tmpMean <- grep("ean", colnames(allInvModelFitSummary$parmatrices)); tmpMean
    tmpSd <- tmpMean+1; tmpSd
    Tvalues <- allInvModelFitSummary$parmatrices[,tmpMean]/allInvModelFitSummary$parmatrices[,tmpSd]; Tvalues
    allInvDrift_Coeff <- cbind(allInvModelFitSummary$parmatrices, Tvalues); allInvDrift_Coeff
    allInvDrift_Coeff[, tmpMean:(dim(allInvDrift_Coeff)[2])] <- round(allInvDrift_Coeff[, tmpMean:(dim(allInvDrift_Coeff)[2])], digits); allInvDrift_Coeff
    # re-label
    if (ctsem341) {
      tmp1 <- which(allInvDrift_Coeff[, "matrix"] == "DRIFT")
      driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); driftNamesTmp
      rownames(allInvDrift_Coeff) <- paste0(allInvDrift_Coeff[, c("matrix")], "_",
                                            allInvDrift_Coeff[, c("row")], "_",
                                            allInvDrift_Coeff[, c("col")])
    } else {
      tmp1 <- which(rownames(allInvDrift_Coeff) == "DRIFT")
      driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); driftNamesTmp
    }
    rownames(allInvDrift_Coeff)[tmp1] <- paste0("DRIFT ", driftNamesTmp); allInvDrift_Coeff

    #tmp <- grep("toV", rownames(allInvModelFitSummary$popmeans)); tmp
    if (ctsem341) {
      homAll_Drift_Coef <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DRIFT", tmpMean]; homAll_Drift_Coef
      homAll_Drift_SE <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DRIFT", tmpSd]; homAll_Drift_SE
      tmp1 <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DRIFT", "2.5%"]; tmp1
      tmp2 <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DRIFT", "97.5%"]; tmp2
    } else {
      homAll_Drift_Coef <- allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DRIFT", "Mean"]; homAll_Drift_Coef
      homAll_Drift_SE <-allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DRIFT", "Sd"]; homAll_Drift_SE
      tmp1 <- c(matrix(allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=TRUE)); tmp2
    }
    names(homAll_Drift_Coef) <- driftNamesTmp; homAll_Drift_Coef
    names(homAll_Drift_SE) <- driftNamesTmp; homAll_Drift_SE
    homAll_Drift_CI <- c(rbind(tmp1, tmp2)); homAll_Drift_CI
    tmp3 <- c(rbind(paste0(driftNamesTmp, "LL"),
                    paste0(driftNamesTmp, "UL"))); tmp3
    names(homAll_Drift_CI) <- tmp3; homAll_Drift_CI
    homAll_Drift_Tvalue <- homAll_Drift_Coef/homAll_Drift_SE; homAll_Drift_Tvalue


    if (ctsem341) {
      homAll_Diffusion_Coef <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", tmpMean]; homAll_Diffusion_Coef
      homAll_Diffusion_SE <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", tmpSd]; homAll_Diffusion_SE
      tmp1 <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "2.5%"]; tmp1
      tmp2 <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "97.5%"]; tmp2
    } else {
      homAll_Diffusion_Coef <- allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DIFFUSIONcov", "Mean"]; homAll_Diffusion_Coef
      homAll_Diffusion_SE <-allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DIFFUSIONcov", "Sd"]; homAll_Diffusion_SE
      tmp1 <- c(matrix(allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DIFFUSIONcov", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "DIFFUSIONcov", "97.5%"], n.latent, byrow=TRUE)); tmp2
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
      homAll_T0VAR_Coef <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "T0cov", tmpMean]; homAll_T0VAR_Coef
      homAll_T0VAR_SE <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "T0cov", tmpSd]; homAll_T0VAR_SE
      tmp1 <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "T0cov", "2.5%"]; tmp1
      tmp2 <- allInvModelFitSummary$parmatrices[allInvModelFitSummary$parmatrices[, "matrix"] == "T0cov", "97.5%"]; tmp2
    } else {
      homAll_T0VAR_Coef <- allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "T0VAR", "Mean"]; homAll_T0VAR_Coef
      homAll_T0VAR_SE <-allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "T0VAR", "Sd"]; homAll_T0VAR_SE
      tmp1 <- c(matrix(allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "T0VAR", "2.5%"], n.latent, byrow=TRUE)); tmp1
      tmp2 <- c(matrix(allInvModelFitSummary$parmatrices[rownames(allInvModelFitSummary$parmatrices) == "T0VAR", "97.5%"], n.latent, byrow=TRUE)); tmp2
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
    #allInvModelFitSummary
    homAll_Minus2LogLikelihood <- -2 * allInvModelFitSummary$loglik; homAll_Minus2LogLikelihood
    homAll_estimatedParameters <- allInvModelFitSummary$npars; homAll_estimatedParameters
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

  DRIFT <- allInvModelFit$modelResults$DRIFT; DRIFT
  DIFFUSION <- allInvModelFit$modelResults$DIFFUSION; DIFFUSION
  T0VAR <- allInvModelFit$modelResults$T0VAR; T0VAR

  Msg <- "################################################################################# \n# Use estimates from all fixed model to compute corr-matrices for all time lags # \n################################################################################ \n################# Set up required discrete time lavaan models ################### \n#################################################################################"
  message(Msg)
  {
    # full lavaan model setup
    {
      modelText <- c()
      # check for drift effects fixed to 0 (reverse order here)
      driftNamesTmp <- driftNamesBackup; driftNamesTmp
      tmp1 <- rownames(ctmaInitFit$studyFitList[[1]]$resultsSummary$popmeans); tmp1
      tmp2 <- which(!(driftNamesTmp %in% tmp1)); tmp2
      if (length(tmp2) != 0) driftNamesTmp[tmp2] <- "0"
      driftNamesTmp <- c(matrix(driftNamesTmp, n.latent, byrow=FALSE)); driftNamesTmp

      counter <- 0
      for (i in 1:n.latent) {
        counter <- counter + 1
        modelText[counter] <- ""
        for (j in 1:n.latent) {
          tmp1 <- paste0("V", j, "toV", i); tmp1
          if(tmp1 %in% driftNamesTmp) {  # if drift effect is not fixed to 0
            if (j == 1) modelText[counter] <- paste0(modelText[counter], paste0("V", i, "T1 ~ V", j, "T0"))
            if (j != 1) modelText[counter] <- paste0(modelText[counter], paste0(" + V", j, "T0"))
          }
          modelText
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
        if (i != j) {
          counter <- counter + 1
          if (i == 1) {
            toDelete <- paste0("\\+", " V", j, "T0"); toDelete
          } else {
            toDelete <- paste0("V", j, "T0 ", "\\+"); toDelete
          }
          tmp2[i] <- gsub(toDelete, "", tmp2[i]); tmp2[i]
          tmp3 <- paste(tmp2, sep="", collapse="\n "); tmp3
          model.wo[[counter]] <- paste0("\n ", tmp3, "\n"); model.wo[[counter]]
          model.wo[[counter]] <- gsub("\\ \\ ", " ", model.wo[[counter]])
        }
      }
    }
  }

  Msg <- "################################################################################# \n###### Now computing implied correlation matrices for different time lags ####### \n#################################################################################"
  message(Msg)
  {
    # functions to compute dt-coefficients
    discreteDriftFunction <- function(driftMatrix, timeScale, number) {
      discreteDriftValue <- OpenMx::expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }
    discreteDiffusionFunction <- function(diffusionMatrix, driftMatrix, timeScale, number) {
      driftHatch <- driftMatrix %x% diag(dim(diffusionMatrix)[1]) + diag(dim(diffusionMatrix)[1]) %x% driftMatrix
      discreteDiffusionValue <- solve(driftHatch) %*% (OpenMx::expm(timeScale %x% driftHatch) - diag(dim(driftHatch)[1])) %*% c(diffusionMatrix)
      discreteDiffusionValue[number] }

    implCov <- list()
    tmp1 <- paste0("V", 1:n.latent); tmp1
    varNames <- paste0(rep(tmp1, 2), c(rep("T0", n.latent), rep("T1", n.latent))); varNames
    pValues <- matrix(NA, nrow=length(usedTimeRange), ncol=(1+n.latent^2-n.latent)); pValues[1,]
    #tmp1 <- diag(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); tmp1
    tmp1 <- diag(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); tmp1
    #tmp2 <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); tmp2
    tmp2 <- c(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); tmp2
    #colnames(pValues) <- c("Interval", paste0("p(", tmp2[!(tmp2 %in% tmp1)], ")") ); colnames(pValues)
    colnames(pValues) <- c("Interval", paste0("(", tmp2[!(tmp2 %in% tmp1)], ")") ); colnames(pValues)

    #for (t in 2:(length(usedTimeRange)-1)) {
    for (t in 2:(length(usedTimeRange))) {
      #t <- 13
      # compute the covariance matrix implied by the dt-coefficients across requested time intervals
      # for equations see Neudecker, H. & Satorra, A. (1991). Linear structural relations: Gradient and ...
      # ... Hessian of the fitting function. Statistics & Probability Letters 11 (1991) 57-61. North-Holland
      #
      beta <- discreteDriftFunction(DRIFT, usedTimeRange[t]); beta
      tmpMat <- matrix(0, n.latent, n.latent); tmpMat
      B <- diag(1, 2*n.latent) - cbind(rbind(tmpMat, beta), rbind(tmpMat, tmpMat)); B
      phi <- T0VAR; phi
      psi <- matrix(discreteDiffusionFunction(DIFFUSION, DRIFT, usedTimeRange[t], 1:4), n.latent); psi
      phi <- cbind(rbind(phi, tmpMat), rbind(tmpMat, tmpMat)); #phi
      psi <- cbind(rbind(tmpMat, tmpMat), rbind(tmpMat, psi)); #psi
      tmp <- solve(B) %*% phi %*% t(solve(B)) + psi; #tmp
      implCov[[t]] <- stats::cov2cor(tmp); #implCov[[t]]
      rownames(implCov[[t]]) <- varNames
      colnames(implCov[[t]]) <- varNames

      # Compute p-value (used in next section for min max intervals for which effects are sign.)
      if ( (is.null(failSafeN)) & ((is.null(failSafeP)))) {
        failSafeP <- .01
        failSafeN <- mean(allSampleSizes)
      }
      model.full.fit2 <- lavaan::sem(model.full, sample.cov = implCov[[t]], sample.nobs = failSafeN)

      # The following lines just extract p-values from lavaans results, but the result is delivered in 'strange' format.
      # Strange format means the R can easily handle the fit objects, but NOT within a package.
      # Therefore, we developed some weird code that finally turned out to work.
      # str(model.full.fit2)
      tmp3 <- which(model.full.fit2@ParTable$op == '~'); tmp3
      tmp4 <- model.full.fit2@ParTable$lhs; tmp4
      tmp4a <- gsub("T1", "", tmp4[tmp3]); tmp4a
      tmp4 <- model.full.fit2@ParTable$rhs; tmp4
      tmp4b <- gsub("T0", "", tmp4[tmp3]); tmp4b
      tmp5 <- which(tmp4a != tmp4b); tmp5
      est <- model.full.fit2@Fit@est[tmp5]; est
      se <- model.full.fit2@Fit@se[tmp5]; se
      z <- est/se; z
      p <- 2*pnorm(z, lower.tail=FALSE); p
      pValues[t, ] <- c(usedTimeRange[t], p); pValues[t, ]
      #
    }
  }


  Msg <- "################################################################################# \n# Compute min and max discrete time intervals for which effects are significant # \n#################################################################################"
  message(Msg)
  {
    targetNames <- colnames(pValues)[-1]; targetNames

    # eliminate drift effects that were fixed to 0
    tmp1 <- which(targetNames == "(0)"); tmp1
    if (length(tmp1) != 0) targetNames <- targetNames[-tmp1]; targetNames

    significanceRange <- c()
    for (i in 1:(length(targetNames))) {
      #i <- 1
      #pValues
      #targetNames[i]
      #pValues[,targetNames[i]]
      tmp1 <- suppressWarnings(usedTimeRange[min(which(pValues[,targetNames[i]] < failSafeP))]); tmp1
      tmp2 <- suppressWarnings(usedTimeRange[max(which(pValues[,targetNames[i]] < failSafeP))]); tmp2
      if (!(is.na(tmp1))) {
        tmp3 <- paste0("The shortest interval across which the effect ", targetNames[i], " is significant "); tmp3
        tmp4 <- paste0("with p < ", failSafeP, " assuming N = ", round(failSafeN, 0), " ", failSafeNhelper, " is ", tmp1, ". "); tmp4
      } else {
        tmp3 <- paste0("There is no shortest interval across which the effect ", targetNames[i], " is significant "); tmp3
        tmp4 <- paste0("with p < ", failSafeP, " assuming N = ", round(failSafeN, 0), " ", failSafeNhelper, ". "); tmp4
      }
      if (!(is.na(tmp2))) {
        tmp5 <- paste0("The longest interval across which the effect ", targetNames[i], " is significant "); tmp5
        tmp6 <- paste0("with p < ", failSafeP, " assuming N = ", round(failSafeN, 0), " ", failSafeNhelper, " is ", tmp2, ". "); tmp6
      } else {
        tmp5 <- paste0("There is no longest interval across which the effect ", targetNames[i], " is significant "); tmp5
        tmp6 <- paste0("with p < ", failSafeP, " assuming N = ", round(failSafeN, 0), ". "); tmp6
      }
      tmp7 <- NULL
      if (is.null(timeRange)) {
        tmp7 <- paste0("Note that you have not provided an explicit time range for analysis of statistical power. "); tmp7
        tmp7 <- paste0(tmp7, "The time intervals used ranged from 1 to 1.5 times the longest interval used "); tmp7
        tmp7 <- paste0(tmp7, "in the primary studies, using integer steps of 1.0. These intervals were then "); tmp7
        tmp7 <- paste0(tmp7, "augmented by time intervals found in primary studies that were non-integers."); tmp7
      }
      significanceRange[i] <- paste(tmp3, tmp4, tmp5, tmp6, tmp7); significanceRange[i]
    }
  }

  Msg <- "################################################################################# \n########### Compute required sample sizes to achieve requested power ############ \n#################################################################################"
  message(Msg)

  # Fast function to calculate required sample sizes later (as optional replacement for ss.power.reg.coef)
  nestedProbFunT <- function (fvalue, alpha=.05, power=.80, p=2, x) (1-
                                                                       stats::pt(
                                                                         stats::qt((1 - alpha/2), df = (x)-p-1,
                                                                                   lower.tail = TRUE, log.p = FALSE),
                                                                         df = (x)-p-1, ncp = sqrt(x) * abs(fvalue),
                                                                         lower.tail = TRUE, log.p = FALSE)) - power

  # Create table: sampleSizes x deltas (of primary studies) for post hoc power calculations
  {
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
  }

  # Loop through a range of lags to determine sample sizes (same parameters as for plotting the effects furter below)
  plotPairs <- array(dim=c(n.latent^2, length(statisticalPower), length(usedTimeRange), 2))  # all drift effects, all powers, time range, timePoint+SampleSize

  # Plot required sample size for cross effects.
  for (h in 1:length(statisticalPower)) {
    counter <- 0
    #h <- 1
    for (j1 in 1:(n.latent)) {
      #j1 <- 1
      for (j2 in 1:(n.latent)) {
        #j2 <- 2
        if (j1 != j2) {
          counter <- counter + 1; counter
          for (k in 1:(length(usedTimeRange)-1)) {
            #k <- 1
            delta_t <- usedTimeRange[k+1]; delta_t
            plotPairs[counter, h, k, 1] <- usedTimeRange[k+1]; plotPairs[counter, h, k, 1] # time point
            #plotPairs[1, , ,]
            # R2 in terms of Kelley & Maxwell 2008
            # betas & psis for model with all effects included
            sample.nobs <- 1000 # large enough to prevent shrinkage
            model.full.fit <- lavaan::sem(model.full, sample.cov = implCov[[k+1]], sample.nobs = sample.nobs)
            tmp <- lavaan::inspect(model.full.fit, "est"); tmp
            psi <- diag(tmp$psi)[grep("T1", names(diag(tmp$psi)))]; psi
            beta <- tmp$beta[1:n.latent ,(n.latent+1):(2*n.latent)]; beta
            R2 <- 1 - psi[j1]; R2

            effectSizes[k, (n.latent^2-n.latent+1-counter)] <- beta[j2, j1]; effectSizes[k, (n.latent^2-n.latent+1-counter)]

            # R2 without j (cross effect) in terms of Kelley & Maxwell 2008
            model.wo.fit <- lavaan::sem(unlist(model.wo[[counter]]),
                                        sample.cov = implCov[[k+1]],
                                        sample.nobs = sample.nobs)
            tmp <- lavaan::inspect(model.wo.fit, "est"); tmp
            R2.j <- 1 - tmp$psi[j1,j1]; R2.j

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
                # The following uses our own function (much faster)
                signalToNoiseRatios <- sqrt((R2-R2.j)/(1-R2)); signalToNoiseRatios
                helper <- round(rootSolve::uniroot.all(nestedProbFunT, c(n.latent+2,999999999),
                                                       fvalue=signalToNoiseRatios, alpha=.05,
                                                       power=statisticalPower[h], p=n.latent) + .49999, 0)
                if (length(helper) < 1) helper <- NA
                plotPairs[counter, h, k, 2] <- helper; plotPairs[counter, h, k, 2]
              }

              # Post hoc power computations
              if ( (delta_t %in% tableNxDeltas[ ,-1]) & (h == 1) ){  # do only once (not for all a priori powers)
                M <- tableNxDeltas[ ,-1] == delta_t; M # temp matrix used below
                empiricalN <- matrix(tableNxDeltas[apply(M, 1, any), ], ncol=maxTpoints)[,1]; empiricalN
                empiricalN <- stats::na.omit(empiricalN); empiricalN
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
    Msg <- paste0("################################################################################# \n###################### Power calculation for ", statisticalPower[h], " completed ###################### \n#################################################################################")
    message(Msg)
  } # end h loop (length(statisticalPower))

  # shortcut: eliminate effects that were fixed to zero (if props = .25 or .005 = by chance)
  # deactivated on 2.6.2023
  skip <- 1
  if (skip != 1) {
  for (l in length(listPowerAlpha05):1) {
    tmp1 <- apply(listPowerAlpha05[[l]], 2, mean, na.rm=TRUE); tmp1
    (round(tmp1[2], 4) == .0250)
    if (round(tmp1[2], 4) == .0250) listPowerAlpha05[[l]] <- NULL
    tmp1 <- apply(listPowerAlpha01[[l]], 2, mean, na.rm=TRUE); tmp1
    if (round(tmp1[2], 4) == .0050) listPowerAlpha01[[l]] <- NULL
  }
  }

  # Table of required sample sizes for range of different time lags (a priori power)
  requiredSampleSizes <- list()
  currentDriftNames <- c()
  counter1 <- counter2 <- 0
  for (j1 in 1:(n.latent)) {
    for (j2 in 1:(n.latent)) {
      counter1 <- counter1 + 1
      if (j1 != j2 ) {
        counter2 <- counter2 + 1
        if ( paste0("V", j2, "toV", j1) %in% driftNames) {
          requiredSampleSizes[[counter2]] <- plotPairs[counter2, , , 2]
          currentDriftNames <- c(currentDriftNames, driftNames[counter1])
          rowNames  <- plotPairs[counter2, 1, , 1]
        }
      }
    }
  }
  # eliminate RSS for drift effects that were fixed to zero
  for (l in length(requiredSampleSizes):1) if(is.null(requiredSampleSizes[[l]])) requiredSampleSizes[[l]] <- NULL

  # re-structure into a single table and replace 100000 by NA
  tmp1 <- n.latent^2-n.latent; tmp1
  tmp2 <- length(which(driftNames == "0")); tmp2
  numberOfEffects <- tmp1 - tmp2; numberOfEffects

  if (!(is.null(dim(requiredSampleSizes[[1]])[1]))){
    nrows <- dim(requiredSampleSizes[[1]])[1]
  } else {
    nrows <- 1
  }

  tmp <- matrix(requiredSampleSizes[[1]], nrow=nrows); tmp
  if (numberOfEffects > 1) for (h in 2:(numberOfEffects)) tmp <- rbind(tmp, requiredSampleSizes[[h]])

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

  # add (not really standardized) effect sizes based on matrix exponentiation
  tmp1 <- as.numeric(rownames(requiredSampleSizes)); tmp1
  timeLags <- tmp1[!(is.na(tmp1))]; timeLags
  effectSizes2 <- matrix(NA, nrow=length(usedTimeRange), ncol=numberOfEffects); effectSizes2
  effectCounter <- 0
  for (i in 1:n.latent) {
    for (j in 1:n.latent) {
      if (i != j) {
        if ( paste0("V", j, "toV", i) %in% driftNames) {
          effectCounter <- effectCounter + 1; effectCounter
          rowCounter <- 0
          for (k in 1:(length(usedTimeRange)-1)) {
            rowCounter <- rowCounter + 1; rowCounter
            A <- OpenMx::expm(DRIFT %x% usedTimeRange[k+1])
            effectSizes2[rowCounter, effectCounter] <- round(A[i,j], digits)
          }
        }
      }
    }
  }

  # shortcut: eliminate effects sizes fixed to 0
  tmp1 <- apply(effectSizes, 2, mean, na.rm=TRUE); tmp1
  tmp2 <- which(tmp1 == 0); tmp2

  if (length(tmp2) != 0) effectSizes <- effectSizes[, -tmp2]
  if (is.null(dim(effectSizes)[1])) effectSizes <- matrix(effectSizes, ncol=1)

  diffDim <- dim(requiredSampleSizes)[1] - dim(effectSizes)[1]; diffDim
  if (diffDim != 0) {
    helperMat <- matrix(NA, nrow=diffDim, ncol=dim(effectSizes)[2]); helperMat
    effectSizes <- rbind(effectSizes, helperMat); effectSizes
  }
  tmp1 <- matrix(driftNames, n.latent, n.latent); tmp1
  diag(tmp1) <- NA
  tmp1 <- tmp1[!(is.na(tmp1))]; tmp1
  tmp1 <- tmp1[tmp1 != "0"]; tmp1
  colnames(effectSizes) <- tmp1
  requiredSampleSizes <- cbind(requiredSampleSizes, effectSizes)

  # Determine optimal time lag in terms of min sample size required
  rowNames <- c(rownames(requiredSampleSizes), "Min N", "Opt. Lag"); rowNames
  minN <- (apply(requiredSampleSizes, 2, min, na.rm=TRUE)); minN
  optimalCrossLagN <- c()
  for (h in 1:(dim(requiredSampleSizes)[2])) optimalCrossLagN[h] <- mean(which(requiredSampleSizes[ ,h ] == minN[h]))
  # CHD next line replace by line below 10. OCT 2022
  #optimalCrossLagN <- optimalCrossLagN*stepWidth; optimalCrossLagN
  optimalCrossLagN <- as.numeric(rownames(requiredSampleSizes))[round(optimalCrossLagN-.4, 0)]
  requiredSampleSizes <- rbind(requiredSampleSizes, minN, optimalCrossLagN)
  rownames(requiredSampleSizes) <- rowNames
  tmp1 <- grep("Power", colnames(requiredSampleSizes)); tmp1
  tmp2 <- 1:ncol(requiredSampleSizes); tmp2
  tmp3 <- tmp2[!(tmp2 %in% tmp1)]; tmp3
  requiredSampleSizes[c("Min N", "Opt. Lag"), tmp3] <-NA


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
    meanPower <- apply(tmp[, targetCols[1:2]], 2, mean, na.rm=TRUE); meanPower
    medianPower <- apply(tmp[, targetCols[1:2]], 2, stats::median, na.rm=TRUE); medianPower
    postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
    postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
    postHocPower[dim(postHocPower)[1]-1, targetCols[1:2]] <- round(meanPower, digits)
    postHocPower[dim(postHocPower)[1], targetCols[1:2]] <- round(medianPower, digits)
    newNames <- c(paste0("Study_No_", 1:n.studies),
                  "Mean", "Median")
    rownames(postHocPower) <- newNames
    postHocPowerList[[j]] <- postHocPower
    names(postHocPowerList)[[j]] <- currentDriftNames[j]
  }

  allInvModelFit$resultsSummary <- allInvModelFitSummary

  results <- list(activeDirectory=activeDirectory,
                  plot.type=c("power"), model.type="stanct", #model.type="mx",
                  coresToUse=coresToUse, n.studies=1,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=allInvModelFit,
                  emprawList=NULL,
                  statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=DRIFT, DIFFUSION=DIFFUSION, T0VAR=T0VAR, CINT=NULL),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="Analysis of Statistical Power and Required Sample Sizes",
                               estimates=list("Estimates of Model with all Effects Invariant"=round(homAll_effects, digits),
                                              "Requested Statistical Power"=statisticalPower,
                                              "Power (post hoc) for Drift Effects"=postHocPowerList,
                                              "Required Sample Sizes"=round(requiredSampleSizes, digits),
                                              "Effect Sizes (based on discrete-time calcs; used for power calcs.)"=round(effectSizes, digits),
                                              "Range of significant effects"=significanceRange)))
  class(results) <- "CoTiMAFit"

  invisible(results)

} ### END function definition
