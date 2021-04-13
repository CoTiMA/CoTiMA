#' ctmaBiGOMX
#'
#' @description Analysis of publication bias and fixed and ranom effects analysis of single drift coefficients if OLD OpenMx fit files are supplied
#'
#' @param ctmaInitFit fit object created with ctmaInti containing the fitted ctsem model of each primary study
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param PETPEESEalpha   # probability level (condition) below which to switch from PET to PEESE (Stanley, 2017, SPPS,p. 582, below Eq. 2; (default p = .10)
#' @param activateRPB if TRUE, messages (warning, finishs) could be send to smart phone (default = FALSE)
#' @param digits rounding (default = 4)
#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red
#' @importFrom stats lm pnorm
#' @importFrom OpenMx vec2diag diag2vec
#' @importFrom utils capture.output
#' @importFrom zcurve zcurve
#'
#' @return returns a CoTiMA fit object with results of publication bias analysis, fixed and random effect analysis, Egger's tests, PET-PEESE corrections.
#'
ctmaBiGOMX <- function(
  ctmaInitFit=NULL,
  activeDirectory=NULL,
  PETPEESEalpha=.10,
  activateRPB=FALSE,
  digits=4
)

{  # begin function definition (until end of file)

  { ### CHECKS
    # check if fit object is specified
    if (is.null(ctmaInitFit)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA (\"ctmaInitFit\") object has to be supplied to analyse something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

  }

  base2Full <- function(baseMat=NULL) {
    result <- OpenMx::vec2diag(exp(OpenMx::diag2vec(baseMat))) + baseMat - OpenMx::vec2diag(OpenMx::diag2vec(baseMat))
    result <- result %*% t(result)
    return(result) # @ instead of $ introduced
  }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {
    n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    all_Coeff <- (lapply(ctmaInitFit$studyFitList, function(extract) extract$mxobj$output$estimate)); all_Coeff
    all_SE <- (lapply(ctmaInitFit$studyFitList, function(extract) extract$mxobj$output$standardErrors)); all_SE
    allSampleSizes <- unlist(lapply(ctmaInitFit$studyList, function(extract) extract$sampleSize)); allSampleSizes
    allSampleSizes <- allSampleSizes[-length(allSampleSizes)]; allSampleSizes
  }


  #######################################################################################################################
  ##################################### Analyses of Publication Bias ####################################################
  #######################################################################################################################

  {
    DRIFTCoeff <- DIFFUSIONCoeff <- T0VARCoeff <- DRIFTSE <- DIFFUSIONSE <- T0VARSE <- matrix(NA, n.studies, n.latent^2)
    for (i in 1:n.studies) {
      #i <- 1
      tmp1 <- all_Coeff[[i]]; tmp1
      DRIFTCoeff[i, ] <- tmp1[grep("toV", names(tmp1))]; DRIFTCoeff
      tmp2  <- tmp1[grep("diff", names(tmp1))]; tmp2
      DIFFUSIONCoeff[i, ] <- base2Full(OpenMx::vech2full(tmp2)); DIFFUSIONCoeff
      tmp2  <- tmp1[grep("T0var", names(tmp1))]; tmp2
      T0VARCoeff[i, ] <- base2Full(OpenMx::vech2full(tmp2)); T0VARCoeff

      tmp1 <- all_SE[[i]]; tmp1
      DRIFTSE[i, ] <- tmp1[grep("toV", rownames(tmp1))]; DRIFTSE
      tmp2  <- tmp1[grep("diff", rownames(tmp1))]; tmp2
      DIFFUSIONSE[i, ] <- OpenMx::vech2full(tmp2); DIFFUSIONSE
      tmp2  <- tmp1[grep("T0var", rownames(tmp1))]; tmp2
      T0VARSE[i, ] <- OpenMx::vech2full(tmp2); T0VARSE
    }

    tmp1 <- names(all_Coeff[[1]][grep("toV", names(all_Coeff[[1]]))]); tmp1
    colnames(DRIFTCoeff) <- colnames(DRIFTSE) <- tmp1; DRIFTCoeff
    tmp2 <- c(OpenMx::vech2full(names(all_Coeff[[1]][grep("diff", names(all_Coeff[[1]]))]))); tmp2
    colnames(DIFFUSIONCoeff) <- colnames(DIFFUSIONSE) <- tmp2; DIFFUSIONCoeff
    tmp3 <- c(OpenMx::vech2full(names(all_Coeff[[1]][grep("T0var", names(all_Coeff[[1]]))]))); tmp3
    colnames(T0VARCoeff) <- colnames(T0VARSE) <- tmp3; T0VARCoeff


    DRIFTCoeffSND <- DRIFTCoeff / DRIFTSE; DRIFTCoeffSND
    DRIFTPrecision <- c(rep(1, n.latent^2))/(DRIFTSE); DRIFTPrecision
    colnames(DRIFTPrecision) <- colnames(DRIFTCoeffSND); DRIFTPrecision

    message1 <- "The pos. & sign. intercept indicates that SMALLER studies produced more positive (or less negative) effects"
    message2 <- "The neg. & sign. intercept indicates that LARGER studies produced more positive (or less negative) effects"
    tmp <- c()

    eggerDrift <- list()
    for (j in 1:(n.latent^2)) {
      eggerDrift[[j]] <- stats::lm(DRIFTCoeffSND[,j]~DRIFTPrecision[,j]) # This is identical to a weighted regression of drift on se ...
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


    #######################################################################################################################
    ################################### Fixed & Random Effects Analyses ###################################################
    #######################################################################################################################

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

    ############################################## zcurve Analysis ###################################################
    zFit <- list()
    for (i in 1: dim(DRIFTCoeffSND)[2]) {
      tmp1 <- abs(DRIFTCoeffSND[, i]); tmp1
      zFit[[i]] <- summary(zcurve::zcurve(z=tmp1))
    }
    names(zFit) <- paste0("Z-Curve 2.0 analysis of ", colnames(DRIFTCoeffSND)); zFit
    ## format results
    #z.CurveResults <- matrix(NA, ncol=length(zFit), nrow=19)
    #for (h in 1:length(zFit)) {
    #  stopRow = 0
    #  for (j in 2:length(zFit[[i]])) {
    #    startRow  <- stopRow + 1; startRow
    #    stopRow <- startRow + length(unlist((zFit[[h]][j]))) - 1; stopRow
    #    unlist((zFit[[h]][j]))
    #    if (j == 2) z.CurveResults[startRow:stopRow, h] <- round(c((matrix(unlist((zFit[[h]][j])), ncol=2, byrow=TRUE))), digits)
    #    if (j != 2) z.CurveResults[startRow:stopRow, h] <- unlist((zFit[[h]][j]))
    #  }
    #}
    #rownamesPart <- c("ERR", "ERR lower CI", "ERR upper CI", "EDR", "EDR lower CI", "EDR upper CI" )
    #rownamesPart <- c(rownamesPart, names(unlist((zFit[[1]][3]))), names(unlist((zFit[[1]][4]))), names(unlist((zFit[[1]][5]))) )
    #rownames(z.CurveResults) <- rownamesPart


    # Combine results
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


    ### some corrections for the output
    heterogeneity <- fixedEffectDriftResults[9:16,]; heterogeneity
    fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults
    eggerTest <- PET_PEESE_DRIFTresults[9:12,]; eggerTest
    PET_PEESE_DRIFTresults <- PET_PEESE_DRIFTresults[1:8,]; PET_PEESE_DRIFTresults

  } ### END Fixed & Random Effects Analyses ###

  results <- list(activeDirectory=activeDirectory, #sourceDirectory=sourceDirectory,
                  plot.type=c("funnel", "forest"), model.type="BiG",
                  coresToUse=NULL, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=NULL, # , homDRIFTallFitCI),
                  emprawList=NULL,
                  statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
                                    DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="Analysis of Publication Bias & Generalizability",
                               estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
                                              "Heterogeneity"=round(heterogeneity, digits),
                                              "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
                                              "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
                                              "Egger's tests"=round(eggerTest, digits),
                                              "Z-Curve 2.0 Results:"=zFit)))

  class(results) <- "CoTiMAFit"

  invisible(results)

} ### END function definition
