#######################################################################################################################
############################################ CoTiMA Publication Bias ##################################################
#######################################################################################################################
#' ctmaBias
#'
#' @param ctmaInitFit ""
#' @param activeDirectory ""
#' @param PETPEESEalpha ""
#' @param activateRPB ""
#' @param digits ""
#'
#' @return
#' @export
#'
ctmaBias <- function(
  # Primary Study Fits
  ctmaInitFit=NULL,                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,

  # probability level (condition) below which to switch from PET to PEESE (Stanley, 2017, SPPS,p. 582, below Eq. 2)
  PETPEESEalpha=.10,

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  digits=4
)


{  # begin function definition (until end of file)

  ### check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA (\"ctmaInitFit\") object has to be supplied to analyse something. \n"))
    stop("Good luck for the next try!")
  }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  { # start extracting
    n.latent <- ctmaInitFit$n.latent; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- ctmaInitFit$n.studies; n.studies
    tmp1 <- matrix(unlist(ctmaInitFit$modelResults$DRIFT), ncol=n.latent^2, byrow=TRUE); tmp1
    tmp2 <- matrix(unlist(ctmaInitFit$modelResults$DIFFUSION), ncol=(n.latent * (n.latent+1) / 2), byrow=TRUE); tmp2
    tmp3 <- matrix(unlist(ctmaInitFit$modelResults$T0VAR), ncol=(n.latent * (n.latent+1) / 2), byrow=TRUE); tmp3
    all_Coeff <- cbind(tmp1, tmp2, tmp3); all_Coeff
    colnames(all_Coeff) <- c(names(ctmaInitFit$modelResults$DRIFT[[1]]),
                             names(ctmaInitFit$modelResults$DIFFUSION[[1]]),
                             names(ctmaInitFit$modelResults$T0VAR[[1]]))

    all_SE <- matrix(NA, ncol=dim(all_Coeff)[2], nrow=dim(all_Coeff)[1]); all_SE
    # identify lower.tri elements
    tmp1 <- matrix(99, ncol=n.latent, nrow=n.latent); tmp1
    toSelect <- which(lower.tri(tmp1, diag=TRUE)); toSelect
    for (i in 1:n.studies) {
      tmp1 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DRIFT[, , 1],
                    ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DRIFT[, , 2])
      tmp1 <- tmp1[, c(1,3,2,4)]
      tmp2 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DIFFUSION[, , 1],
                    ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DIFFUSION[, , 2])
      tmp2 <- tmp2[, toSelect]
      tmp3 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_T0VAR[, , 1],
                    ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_T0VAR[, , 2])
      tmp3 <- tmp3[, toSelect]
      tmp4 <- cbind(tmp1, tmp2, tmp3)
      apply(tmp4, 2, mean)
      tmp4 <- apply(tmp4, 2, stats::var)^.5
      all_SE[i,] <- t(as.matrix(tmp4))
    }
    colnames(all_SE) <- c(names(ctmaInitFit$modelResults$DRIFT[[1]]),
                          names(ctmaInitFit$modelResults$DIFFUSION[[1]]),
                          names(ctmaInitFit$modelResults$T0VAR[[1]]))

    all_SE
    allSampleSizes <- ctmaInitFit$statisticsList$allSampleSizes; allSampleSizes
  } # end extracting


  #######################################################################################################################
  ##################################### Analyses of Publication Bias ####################################################
  #######################################################################################################################

  #{ # Start Analysis of Publication Bias
    DRIFTCoeff <- DIFFUSIONCoeff <- T0VARCoeff <- DRIFTSE <- DIFFUSIONSE <- T0VARSE <- matrix(NA, n.studies, n.latent^2)

    DRIFTCoeff <- all_Coeff[, grep("toV", colnames(all_Coeff))]; DRIFTCoeff
    DRIFTSE <- all_SE[, grep("toV", colnames(all_SE))]; DRIFTSE
    DIFFUSIONCoeff <- all_Coeff[, grep("diff", colnames(all_Coeff))]; DIFFUSIONCoeff
    DIFFUSIONSE <- all_SE[, grep("diff", colnames(all_SE))]; DIFFUSIONSE
    T0VARCoeff <- all_Coeff[, grep("T0var", colnames(all_Coeff))]; T0VARCoeff
    T0VARSE <- all_SE[, grep("T0var", colnames(all_SE))]; T0VARSE

    DRIFTCoeffSND <- DRIFTCoeff / DRIFTSE; DRIFTCoeffSND
    DRIFTPrecision <- c(rep(1, n.latent^2))/(DRIFTSE); DRIFTPrecision
    colnames(DRIFTPrecision) <- colnames(DRIFTCoeffSND); DRIFTPrecision

    message1 <- "The pos. & sign. intercept indicates that SMALLER studies produced more positive (or less negative) effects"
    message2 <- "The neg. & sign. intercept indicates that LARGER studies produced more positive (or less negative) effects"
    tmp <- c()

    eggerDrift <- list()
    for (j in 1:(n.latent^2)) {
      eggerDrift[[j]] <- stats::lm(DRIFTCoeffSND[,j]~DRIFTPrecision[,j]) # This is identical to a weighted regression of drift on se ...
      eggerDrift[[j]]$message <- "No sign. evidence for publication bias."
      if (summary(eggerDrift[[j]])$coefficients[1,1] > 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message1
      }
      if (summary(eggerDrift[[j]])$coefficients[1,1] < 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message2
      }
    }

    FREAResults <- list()

    FREAResults[[1]] <- "############# Eggers Test for DRIFT Parameter Estimates  ###############################"
    FREACounter <- 1
    for (j in 1:(n.latent^2)) {
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- paste0("-------------------------------- Eggers Test for ",
                                           colnames(DRIFTCoeff)[j], "--------------------------------")
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- eggerDrift[[j]]$message
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- summary(eggerDrift[[j]])
    }


    ################################### Fixed & Random Effects Analyses ###################################################

    # FIXED EFFECTS ANALYSIS ###############################################################################
    DriftMeans <- colMeans(DRIFTCoeff); DriftMeans
    # Sum of within weights  and weight * effect size
    T_DriftWeights <- colSums(DRIFTPrecision^2); T_DriftWeights
    T_DriftMeans <- colSums(DRIFTCoeff * DRIFTPrecision^2); T_DriftMeans
    names(T_DriftMeans) <- names(T_DriftWeights); T_DriftMeans
    # Fixed effects results
    FixedEffect_Drift <- T_DriftMeans/T_DriftWeights; FixedEffect_Drift
    FixedEffect_DriftVariance <- 1/T_DriftWeights; FixedEffect_DriftVariance
    FixedEffect_DriftSE <- FixedEffect_DriftVariance^.5; FixedEffect_DriftSE
    FixedEffect_DriftUpperLimit <- FixedEffect_Drift + 1.96*FixedEffect_DriftSE; FixedEffect_DriftUpperLimit
    FixedEffect_DriftLowerLimit <- FixedEffect_Drift - 1.96*FixedEffect_DriftSE; FixedEffect_DriftLowerLimit
    FixedEffect_DriftZ <- FixedEffect_Drift/FixedEffect_DriftSE; FixedEffect_DriftZ
    FixedEffect_DriftProb <- round(1-stats::pnorm(abs(FixedEffect_DriftZ),
                                           mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log=F), digits=digits); FixedEffect_DriftProb
    Q_Drift <- colSums(DRIFTPrecision^2 * DRIFTCoeff^2)- (colSums(DRIFTPrecision^2 * DRIFTCoeff))^2 / colSums(DRIFTPrecision^2); Q_Drift
    H2_Drift <- Q_Drift/(n.studies-1); H2_Drift
    I2_Drift <- (H2_Drift-1)/H2_Drift*100; I2_Drift
    # Tau square
    #T2_DriftWeights <- colSums(DRIFTPrecision^2); T2_DriftWeights
    T2_DriftWeights <- colSums(DRIFTPrecision^2^2); T2_DriftWeights # Borenstein et al., 2007, p. 98
    cDrift <- T_DriftWeights-T2_DriftWeights/T_DriftWeights; cDrift
    tau2Drift <- (Q_Drift-(n.studies-1))/cDrift; tau2Drift
    SElnHDrift <- c()
    SElnHDrift[] <- 0
    for (j in 1:(n.latent^2)) {
      if (Q_Drift[j] > n.studies) SElnHDrift[j] <- 1/2*(log(Q_Drift[j])-log(n.studies-1))/((2*Q_Drift[j])^.5-(2*(n.studies-1)-1)^.5)
      if (Q_Drift[j] <= n.studies) SElnHDrift[j] <-  (1/(2*(n.studies-2)) * (1 - 1/(3*(n.studies-2)^.5)) )^.5
    }

    H2DriftUpperLimit <- exp(log(H2_Drift) + 1.96*SElnHDrift); H2DriftUpperLimit
    H2DriftLowerLimit <- exp(log(H2_Drift) - 1.96*SElnHDrift); H2DriftLowerLimit
    L <- exp(0.5*log(Q_Drift/(n.studies-1))-1.96*SElnHDrift)
    U <- exp(0.5*log(Q_Drift/(n.studies-1))+1.96*SElnHDrift)
    I2DriftUpperLimit <- (U^2-1)/U^2 * 100; I2DriftUpperLimit
    I2DriftLowerLimit <- (L^2-1)/L^2 * 100; I2DriftLowerLimit

    MeanOfDriftValues <- DriftMeans
    fixedEffectDriftResults <- rbind(MeanOfDriftValues, FixedEffect_Drift, FixedEffect_DriftVariance, FixedEffect_DriftSE,
                                     FixedEffect_DriftUpperLimit, FixedEffect_DriftLowerLimit,
                                     FixedEffect_DriftZ, FixedEffect_DriftProb, tau2Drift, Q_Drift, H2_Drift,
                                     H2DriftUpperLimit, H2DriftLowerLimit, I2_Drift,
                                     I2DriftUpperLimit, I2DriftLowerLimit)

    # RANDOM EFFECTS ANALYSIS ###############################################################################
    # Total variance weighting
    Ttot_DriftWeights <- 0
    Ttot_DriftMeans <- 0
    tau2DriftExtended <- do.call(rbind, replicate(n.studies, tau2Drift, simplify=FALSE))
    Ttot_DriftWeights <-colSums(1/ (DRIFTSE^2 + tau2DriftExtended)); Ttot_DriftWeights
    Ttot_DriftMeans <- colSums(DRIFTCoeff * 1/ (DRIFTSE^2 + tau2DriftExtended)); Ttot_DriftMeans
    # Random effects results
    RandomEffecttot_Drift <- Ttot_DriftMeans/Ttot_DriftWeights; RandomEffecttot_Drift
    RandomEffecttot_DriftVariance <- 1/Ttot_DriftWeights; RandomEffecttot_DriftVariance
    RandomEffecttot_DriftSE <- RandomEffecttot_DriftVariance^.5; RandomEffecttot_DriftSE
    RandomEffecttot_DriftUpperLimit <- RandomEffecttot_Drift + 1.96*RandomEffecttot_DriftSE; RandomEffecttot_DriftUpperLimit
    RandomEffecttot_DriftLowerLimit <- RandomEffecttot_Drift - 1.96*RandomEffecttot_DriftSE; RandomEffecttot_DriftLowerLimit
    RandomEffecttot_DriftZ <- RandomEffecttot_Drift/RandomEffecttot_DriftSE; RandomEffecttot_DriftZ
    RandomEffecttot_DriftProb <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ),
                                               mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log=F), digits=digits); RandomEffecttot_DriftProb
    RandomEffectDriftResults <- rbind(RandomEffecttot_Drift, RandomEffecttot_DriftVariance, RandomEffecttot_DriftSE,
                                      RandomEffecttot_DriftUpperLimit, RandomEffecttot_DriftLowerLimit,
                                      RandomEffecttot_DriftZ, RandomEffecttot_DriftProb)
    RandomEffecttot_DriftUpperLimitPI <- RandomEffecttot_Drift + 1.96*(tau2Drift^.5); RandomEffecttot_DriftUpperLimitPI
    RandomEffecttot_DriftLowerLimitPI <- RandomEffecttot_Drift - 1.96*(tau2Drift^.5); RandomEffecttot_DriftLowerLimitPI
    RandomEffectDriftResults <- rbind(RandomEffecttot_Drift, RandomEffecttot_DriftVariance, RandomEffecttot_DriftSE,
                                      RandomEffecttot_DriftUpperLimit, RandomEffecttot_DriftLowerLimit,
                                      RandomEffecttot_DriftZ, RandomEffecttot_DriftProb,
                                      RandomEffecttot_DriftUpperLimitPI, RandomEffecttot_DriftLowerLimitPI)


    ### PET, PEESE & WLS approaches to correct for bias
    PETDrift_fit <- list()
    PEESEDrift_fit <- list()
    WLSDrift_fit <- list()
    PET_PEESEDrift_fit <- list()
    Egger2Drift_fit <- list()

    sampleSizes <- unlist(allSampleSizes); sampleSizes
    for (ii in 1:(n.latent^2)) {
      driftCoeff <- DRIFTCoeff[ , ii]; driftCoeff
      driftSE <- DRIFTSE[, ii]; driftSE

      # PET
      IV <- driftSE; IV
      DV <- driftCoeff; DV
      currentWeigths <- (1/(driftSE^2)); currentWeigths
      PETDrift_fit[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PETDrift_fit[[ii]]
      summary(PETDrift_fit[[ii]])

      # Egger's Test (alternative but algebraically identical model)
      Egger2Drift_fit[[ii]] <- t(c(summary(PETDrift_fit[[ii]])$coefficients[2,1:4]))

      # PEESE
      IV <- driftSE^2; IV
      DV <- driftCoeff
      currentWeigths <- (1/(driftSE^2)); currentWeigths

      PEESEDrift_fit[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PEESEDrift_fit[[ii]]
      summary(PEESEDrift_fit[[ii]])

      # PET-PEESE
      if ( (summary(PETDrift_fit[[ii]]))$coefficients[1,4] > PETPEESEalpha) {
        PET_PEESEDrift_fit[[ii]] <- PETDrift_fit[[ii]]
      }
      if ( (summary(PETDrift_fit[[ii]]))$coefficients[1,4] <= PETPEESEalpha) {
        PET_PEESEDrift_fit[[ii]] <- PEESEDrift_fit[[ii]]
      }

      # WLS
      DV <- driftCoeff/driftSE
      IV <- 1/driftSE
      WLSDrift_fit[[ii]] <- stats::lm(DV ~ IV + 0); WLSDrift_fit[[ii]]; FixedEffect_Drift[[ii]] # should be identical to FixedEffect_Drift
      WLSDriftSE_fit <- summary(WLSDrift_fit[[ii]])$coefficients[2]; WLSDriftSE_fit; FixedEffect_DriftSE[[ii]] # should outperform FixedEffect_DriftSE
    }


    # Combine results
    {
    PET_Drift <-unlist(lapply(PETDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift
    PET_SE <- c()
    for (k in 1:(n.latent^2)) PET_SE <- c(PET_SE, summary(PETDrift_fit[[k]])$coefficients[1,2])

    PEESE_Drift <-unlist(lapply(PEESEDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift
    PEESE_SE <- c()
    for (k in 1:(n.latent^2)) PEESE_SE <- c(PEESE_SE, summary(PEESEDrift_fit[[k]])$coefficients[1,2])

    PET_PEESE_Drift <-unlist(lapply(PET_PEESEDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift
    PET_PEESE_SE <- c()
    for (k in 1:(n.latent^2)) PET_PEESE_SE <- c(PET_PEESE_SE, summary(PET_PEESEDrift_fit[[k]])$coefficients[1,2])

    WLS_Drift <- unlist(lapply(WLSDrift_fit, function(extract) extract$coefficients)); WLS_Drift
    WLS_SE <- c()
    for (k in 1:(n.latent^2)) WLS_SE <- c(WLS_SE, summary(WLSDrift_fit[[k]])$coefficients[1,2])

    Egger2Drift_results <- matrix(unlist(Egger2Drift_fit), ncol=n.latent^2, nrow=4); Egger2Drift_results

    PET_PEESE_DRIFTresults <- rbind(PET_Drift, PET_SE,
                                    PEESE_Drift, PEESE_SE,
                                    PET_PEESE_Drift, PET_PEESE_SE,
                                    WLS_Drift, WLS_SE,
                                    Egger2Drift_results)
    colnames(PET_PEESE_DRIFTresults) <- colnames(DRIFTCoeff)
    rownames(PET_PEESE_DRIFTresults) <- c(rownames(PET_PEESE_DRIFTresults)[1:8], "Egger's b0", "SE(b0)", "T", "p")
    PET_PEESE_DRIFTresults

    ### some corrections for the output
    heterogeneity <- fixedEffectDriftResults[9:16,]; heterogeneity
    fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults
    eggerTest <- PET_PEESE_DRIFTresults[9:12,]; eggerTest
    PET_PEESE_DRIFTresults <- PET_PEESE_DRIFTresults[1:8,]; PET_PEESE_DRIFTresults
    }

   # } # End Analysis of Publication Bias

  results <- list(activeDirectory=activeDirectory,
                  plot.type=c("funnel", "forest"), model.type="bias",
                  coresToUse=NULL, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=NULL, # , homDRIFTallFitCI),
                  emprawList=NULL,
                  statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
                                    DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="Analysis of Publication Bias",
                               estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
                                              "Heterogeneity"=round(heterogeneity, digits),
                                              "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
                                              "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
                                              "Egger's tests"=round(eggerTest, digits),
                                              "Egger's tests Alt. Version"= FREAResults)))
  class(results) <- "CoTiMAFit"

  #saveRDS(results, paste0(activeDirectory, "CoTiMAbias_allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
