#' ctmaRedHet
#'
#' @description Computes the Reduction in Heterogeneity in drift effects after introducing study-level moderators
#'
#' @param activateRPB if TRUE, messages (warning, finished) could be send to smart phone (default = FALSE)
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param ctmaFitObject ctmaFit Object WITHOUT Moderators (obtained from \code{\link{ctmaFit}} with the arguments WEC=\'TRUE\' and scaleTI=FALSE)
#' @param ctmaFitObjectMod ctmaFit Object WITH Moderators (obtained from \code{\link{ctmaFit}} with the arguments WEC=\'TRUE\' and scaleTI=FALSE)
#' @param digits rounding (default = 4)
#' @param dt A scalar indicating a time interval across which discrete time effects should be estimated and then used for ctmaBiG.
#' @param undoTimeScaling if TRUE, the original time scale is used (timeScale argument possibly used in \code{\link{ctmaInit}} is undone )
#'
#' @importFrom RPushbullet pbPost
#' @importFrom ctsem ctExtract
#' @importFrom stats var lm pnorm
#'
#' @export ctmaRedHet
#'
ctmaRedHet <- function(activateRPB=FALSE,
                       activeDirectory=NULL,
                       ctmaFitObject=NULL,
                       ctmaFitObjectMod=NULL,
                       digits=4,
                       dt=NULL,
                       undoTimeScaling=TRUE)
{  # begin function definition (until end of file)

  {
    # check if correct fit objects are specified ######################
    if (is.null(ctmaFitObject)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA (\"ctmaFit\") object has to be supplied to analyse something.
    \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (is.null(ctmaFitObjectMod)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA (\"ctmaFit\") object with moderators has to be supplied to analyse something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (ctmaFitObject$argumentList$WEC == FALSE){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA (\"ctmaFit\") object estimated using the argument \'WEC=TRUE\' has to be supplied.
    \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (ctmaFitObjectMod$argumentList$WEC == FALSE){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA (\"ctmaFit\") object estimated using the argument \'WEC=TRUE\' has to be supplied.
    \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (( length(ctmaFitObjectMod$argumentList$ind.mod.number) == 0) & ( length(ctmaFitObjectMod$argumentList$mod.number) == 0)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA The ctmaFitObjectMod object does neither involve study-level nor individual-level moderators.
    \nGood luck for the next try!"
      stop(ErrorMsg)
    }
  }

  # Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################

  start.time <- Sys.time(); start.time
  {
    driftNames1 <- ctmaFitObject$parameterNames$DRIFT; driftNames1
    n.latent1 <- ctmaFitObject$n.latent; n.latent1
    n.studies1 <- length(ctmaFitObject$studyList); n.studies1
    finishsamples1 = dim(ctmaFitObject$studyFitList$stanfit$rawposterior)[1]; finishsamples1
    #
    DRIFTCoeff1 <- matrix(unlist(ctmaFitObject$modelResults$DRIFTCoeffViaWECMean), nrow=n.studies1, byrow=T); DRIFTCoeff1
    DRIFTSE1 <- matrix(unlist(ctmaFitObject$modelResults$DRIFTCoeffViaWECSD), nrow=n.studies1, byrow=T); DRIFTSE1
    DRIFTCoeffSND1 <- DRIFTCoeff2 / DRIFTSE2; DRIFTCoeffSND1
    DRIFTPrecision1 <- c(rep(1, n.latent1^2))/(DRIFTSE1); DRIFTPrecision1
    colnames(DRIFTPrecision1) <- driftNames1; DRIFTPrecision1
    #
    DRIFTCoeff2 <- matrix(unlist(ctmaFitObjectMod$modelResults$DRIFTCoeffViaWECMean), nrow=n.studies1, byrow=T); DRIFTCoeff2
    DRIFTSE2 <- matrix(unlist(ctmaFitObjectMod$modelResults$DRIFTCoeffViaWECSD), nrow=n.studies1, byrow=T); DRIFTSE2
    DRIFTCoeffSND2 <- DRIFTCoeff2 / DRIFTSE2; DRIFTCoeffSND2
    DRIFTPrecision2 <- c(rep(1, n.latent1^2))/(DRIFTSE2); DRIFTPrecision2
    colnames(DRIFTPrecision2) <- driftNames2; DRIFTPrecision2
  }

  ## FIXED EFFECTS ANALYSIS ###############################################################################
  {
    DriftMeans1 <- colMeans(DRIFTCoeff1); DriftMeans1
    # Sum of within weights  and weight * effect size
    T_DriftWeights1 <- colSums(DRIFTPrecision1^2); T_DriftWeights1
    # DRIFTPrecision
    T_DriftMeans1 <- colSums(DRIFTCoeff1 * DRIFTPrecision1^2); T_DriftMeans1
    names(T_DriftMeans1) <- names(T_DriftWeights1); T_DriftMeans1
    #
    DriftMeans2 <- colMeans(DRIFTCoeff2); DriftMeans2
    # Sum of within weights  and weight * effect size
    T_DriftWeights2 <- colSums(DRIFTPrecision2^2); T_DriftWeights2
    # DRIFTPrecision
    T_DriftMeans2 <- colSums(DRIFTCoeff2 * DRIFTPrecision2^2); T_DriftMeans2
    names(T_DriftMeans2) <- names(T_DriftWeights2); T_DriftMeans2
  }

  ### Fixed effects results for ctmaFitObjectMod without moderators ######
  {
    FixedEffect_Drift1 <- T_DriftMeans1/T_DriftWeights1; FixedEffect_Drift1
    FixedEffect_DriftVariance1 <- 1/T_DriftWeights1; FixedEffect_DriftVariance1
    FixedEffect_DriftSE1 <- FixedEffect_DriftVariance1^.5; FixedEffect_DriftSE1
    FixedEffect_DriftUpperLimit1 <- FixedEffect_Drift1 + 1.96*FixedEffect_DriftSE1; FixedEffect_DriftUpperLimit1
    FixedEffect_DriftLowerLimit1 <- FixedEffect_Drift1 - 1.96*FixedEffect_DriftSE1; FixedEffect_DriftLowerLimit1
    FixedEffect_DriftZ1 <- FixedEffect_Drift1/FixedEffect_DriftSE1; FixedEffect_DriftZ1
    FixedEffect_DriftProb1 <- round(1-stats::pnorm(abs(FixedEffect_DriftZ1),
                                                   mean=c(rep(0, (n.latent1^2))), sd=c(rep(1, (n.latent1^2))), log.p=F), digits=digits); FixedEffect_DriftProb1
    Q_Drift1 <- colSums(DRIFTPrecision1^2 * DRIFTCoeff1^2)- (colSums(DRIFTPrecision1^2 * DRIFTCoeff1))^2 / colSums(DRIFTPrecision1^2); Q_Drift1
    H2_Drift1 <- Q_Drift1/(n.studies1-1); H2_Drift1
    I2_Drift1 <- (H2_Drift1-1)/H2_Drift1*100; I2_Drift1
    # Tau square
    T2_DriftWeights1 <- colSums(DRIFTPrecision1^2^2); T2_DriftWeights1 # Borenstein et al., 2007, p. 98
    cDrift1 <- T_DriftWeights1 - T2_DriftWeights1/T_DriftWeights1; cDrift1
    tau2Drift1 <- (Q_Drift1 - (n.studies1-1))/cDrift1; tau2Drift1
    SElnHDrift1 <- c()
    SElnHDrift1[] <- 0
    for (j in 1:(n.latent1^2)) {
      if (Q_Drift1[j] > n.studies1) SElnHDrift1[j] <- 1/2*(log(Q_Drift1[j])-log(n.studies1-1))/((2*Q_Drift1[j])^.5-(2*(n.studies1-1)-1)^.5)
      if (Q_Drift1[j] <= n.studies1) SElnHDrift1[j] <-  (1/(2*(n.studies1-2)) * (1 - 1/(3*(n.studies1-2)^.5)) )^.5
    }
    H2DriftUpperLimit1 <- exp(log(H2_Drift1) + 1.96*SElnHDrift1); H2DriftUpperLimit1
    H2DriftLowerLimit1 <- exp(log(H2_Drift1) - 1.96*SElnHDrift1); H2DriftLowerLimit1
    L1 <- exp(0.5*log(Q_Drift1/(n.studies1-1))-1.96*SElnHDrift1)
    U1 <- exp(0.5*log(Q_Drift1/(n.studies1-1))+1.96*SElnHDrift1)
    I2DriftUpperLimit1 <- (U1^2-1)/U1^2 * 100; I2DriftUpperLimit1
    I2DriftLowerLimit1 <- (L1^2-1)/L1^2 * 100; I2DriftLowerLimit1
    # same for dt
  }

  ### Fixed effects results for ctmaFitObjectMod with moderators ######
  {FixedEffect_Drift2 <- T_DriftMeans2/T_DriftWeights2; FixedEffect_Drift2
  FixedEffect_DriftVariance2 <- 1/T_DriftWeights2; FixedEffect_DriftVariance2
  FixedEffect_DriftSE2 <- FixedEffect_DriftVariance2^.5; FixedEffect_DriftSE2
  FixedEffect_DriftUpperLimit2 <- FixedEffect_Drift2 + 1.96*FixedEffect_DriftSE2; FixedEffect_DriftUpperLimit2
  FixedEffect_DriftLowerLimit2 <- FixedEffect_Drift2 - 1.96*FixedEffect_DriftSE2; FixedEffect_DriftLowerLimit2
  FixedEffect_DriftZ2 <- FixedEffect_Drift2/FixedEffect_DriftSE2; FixedEffect_DriftZ2
  FixedEffect_DriftProb2 <- round(1-stats::pnorm(abs(FixedEffect_DriftZ2),
                                                 mean=c(rep(0, (n.latent2^2))), sd=c(rep(1, (n.latent2^2))), log.p=F), digits=digits); FixedEffect_DriftProb2
  Q_Drift2 <- colSums(DRIFTPrecision2^2 * DRIFTCoeff2^2)- (colSums(DRIFTPrecision2^2 * DRIFTCoeff2))^2 / colSums(DRIFTPrecision2^2); Q_Drift2
  H2_Drift2 <- Q_Drift2/(n.studies2-1); H2_Drift2
  I2_Drift2 <- (H2_Drift2-1)/H2_Drift2*100; I2_Drift2
  # same for dt
  # Tau square
  T2_DriftWeights2 <- colSums(DRIFTPrecision2^2^2); T2_DriftWeights2 # Borenstein et al., 2007, p. 98
  cDrift2 <- T_DriftWeights2 - T2_DriftWeights2/T_DriftWeights2; cDrift2
  tau2Drift2 <- (Q_Drift2 - (n.studies2-1))/cDrift2; tau2Drift2
  SElnHDrift2 <- c()
  SElnHDrift2[] <- 0
  for (j in 1:(n.latent2^2)) {
    if (Q_Drift2[j] > n.studies2) SElnHDrift2[j] <- 1/2*(log(Q_Drift2[j])-log(n.studies2-1))/((2*Q_Drift2[j])^.5-(2*(n.studies2-1)-1)^.5)
    if (Q_Drift2[j] <= n.studies2) SElnHDrift2[j] <-  (1/(2*(n.studies2-2)) * (1 - 1/(2*(n.studies2-2)^.5)) )^.5
  }
  H2DriftUpperLimit2 <- exp(log(H2_Drift2) + 1.96*SElnHDrift2); H2DriftUpperLimit2
  H2DriftLowerLimit2 <- exp(log(H2_Drift2) - 1.96*SElnHDrift2); H2DriftLowerLimit2
  L2 <- exp(0.5*log(Q_Drift2/(n.studies2-1))-1.96*SElnHDrift2)
  U2 <- exp(0.5*log(Q_Drift2/(n.studies2-1))+1.96*SElnHDrift2)
  I2DriftUpperLimit2 <- (U2^2-1)/U2^2 * 100; I2DriftUpperLimit2
  I2DriftLowerLimit2 <- (L2^2-1)/L2^2 * 100; I2DriftLowerLimit2
  # same for dt
  }

  fixedEffectDriftResults <- rbind(cbind(t(MeanOfDriftValues1), t(MeanOfDriftValues2)), cbind(t(FixedEffect_Drift1), t(FixedEffect_Drift2)),
                                   cbind(t(FixedEffect_DriftVariance1), t(FixedEffect_DriftVariance2)), cbind(t(FixedEffect_DriftSE1), t(FixedEffect_DriftSE2)),
                                   cbind(t(FixedEffect_DriftUpperLimit1), t(FixedEffect_DriftUpperLimit2)), cbind(t(FixedEffect_DriftLowerLimit1), t(FixedEffect_DriftLowerLimit2)),
                                   cbind(t(FixedEffect_DriftZ1), t(FixedEffect_DriftZ2)), cbind(t(FixedEffect_DriftProb1), t(FixedEffect_DriftProb2)),
                                   cbind(t(tau2Drift1), t(tau2Drift2)), cbind(t(Q_Drift1), t(Q_Drift2)), cbind(t(H2_Drift1), t(H2_Drift2)),
                                   cbind(t(H2DriftUpperLimit1), t(H2DriftUpperLimit2)), cbind(t(H2DriftLowerLimit1), t(H2DriftLowerLimit2)),
                                   cbind(t(I2_Drift1), t(I2_Drift2)), cbind(t(I2DriftUpperLimit1), t(I2DriftUpperLimit2)), cbind(t(I2DriftLowerLimit1), t(I2DriftLowerLimit2)))
  rownames(fixedEffectDriftResults) <- c("MeanOfDriftValues", "FixedEffect_Drift",
                                         "FixedEffect_DriftVariance", "FixedEffect_DriftSE", "FixedEffect_DriftUpperLimit",
                                         "FixedEffect_DriftLowerLimit", "FixedEffect_DriftZ", "FixedEffect_DriftProb",
                                         "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
                                         "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")

  fixedEffectDriftResults1 <- rbind(MeanOfDriftValues1, FixedEffect_Drift1, FixedEffect_DriftVariance1, FixedEffect_DriftSE1,
                                    FixedEffect_DriftUpperLimit1, FixedEffect_DriftLowerLimit1,
                                    FixedEffect_DriftZ1, FixedEffect_DriftProb1, tau2Drift1, Q_Drift1, H2_Drift1,
                                    H2DriftUpperLimit1, H2DriftLowerLimit1, I2_Drift1,
                                    I2DriftUpperLimit1, I2DriftLowerLimit1)
  fixedEffectDriftResults2 <- rbind(MeanOfDriftValues2, FixedEffect_Drift2, FixedEffect_DriftVariance2, FixedEffect_DriftSE2,
                                    FixedEffect_DriftUpperLimit2, FixedEffect_DriftLowerLimit2,
                                    FixedEffect_DriftZ2, FixedEffect_DriftProb2, tau2Drift2, Q_Drift2, H2_Drift2,
                                    H2DriftUpperLimit2, H2DriftLowerLimit2, I2_Drift2,
                                    I2DriftUpperLimit2, I2DriftLowerLimit2)

  fixedEffectDriftMessage <- c()
  if ( (any(I2_Drift1 < 0)) | (any(I2_Drift2 < 0)) ) fixedEffectDriftMessage <- "Negative I2 values can be set to 0.0."
  tau2DriftMessage <- c()
  if ( (any(tau2Drift1 < 0)) | (any(tau2Drift1 < 0)) ) tau2DriftMessage <- "Some tau-squared are negative. Random effects cannot be computed. Possibly a small-k-problem."


  ## RANDOM EFFECTS ANALYSIS ##############################################################################

  # Total variance weighting
  tau2DriftExtended1 <- do.call(rbind, replicate(n.studies1, tau2Drift1, simplify=FALSE))
  Ttot_DriftWeights1 <-colSums(1/ (DRIFTSE1^2 + tau2DriftExtended1)); Ttot_DriftWeights1
  Ttot_DriftMeans1 <- colSums(DRIFTCoeff1 * 1/ (DRIFTSE1^2 + tau2DriftExtended1)); Ttot_DriftMeans1
  # same for dt # tbd

  tau2DriftExtended2 <- do.call(rbind, replicate(n.studies2, tau2Drift2, simplify=FALSE))
  Ttot_DriftWeights2 <-colSums(1/ (DRIFTSE2^2 + tau2DriftExtended2)); Ttot_DriftWeights2
  Ttot_DriftMeans2 <- colSums(DRIFTCoeff2 * 1/ (DRIFTSE2^2 + tau2DriftExtended2)); Ttot_DriftMeans2
  # same for dt # tbd

  ### Random effects results for ctmaFitObject without moderators ####
  {
    RandomEffecttot_Drift1 <- Ttot_DriftMeans1/Ttot_DriftWeights1; RandomEffecttot_Drift1
    RandomEffecttot_DriftVariance1 <- 1/Ttot_DriftWeights1; RandomEffecttot_DriftVariance1
    RandomEffecttot_DriftSE1 <- RandomEffecttot_DriftVariance1^.5; RandomEffecttot_DriftSE1
    RandomEffecttot_DriftUpperLimit1 <- RandomEffecttot_Drift1 + 1.96*RandomEffecttot_DriftSE1; RandomEffecttot_DriftUpperLimit1
    RandomEffecttot_DriftLowerLimit1 <- RandomEffecttot_Drift1 - 1.96*RandomEffecttot_DriftSE1; RandomEffecttot_DriftLowerLimit1
    RandomEffecttot_DriftZ1 <- RandomEffecttot_Drift1/RandomEffecttot_DriftSE1; RandomEffecttot_DriftZ1
    RandomEffecttot_DriftProb1 <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ1),
                                                       mean=c(rep(0, (n.latent1^2))), sd=c(rep(1, (n.latent1^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb1
    RandomEffecttot_DriftUpperLimitPI1 <- RandomEffecttot_Drift1 + 1.96*(tau2Drift1^.5); RandomEffecttot_DriftUpperLimitPI1
    RandomEffecttot_DriftLowerLimitPI1 <- RandomEffecttot_Drift1 - 1.96*(tau2Drift1^.5); RandomEffecttot_DriftLowerLimitPI1
    RandomEffectDriftResults1 <- rbind(RandomEffecttot_Drift1, RandomEffecttot_DriftVariance1, RandomEffecttot_DriftSE1,
                                       RandomEffecttot_DriftUpperLimit1, RandomEffecttot_DriftLowerLimit1,
                                       RandomEffecttot_DriftZ1, RandomEffecttot_DriftProb1,
                                       RandomEffecttot_DriftUpperLimitPI1, RandomEffecttot_DriftLowerLimitPI1)
    # same for dt

    ### some corrections for the output
    heterogeneity1 <- fixedEffectDriftResults1[9:16,]; heterogeneity1

    # same for dt # tbd
  }

  ### Random effects results for ctmaFitObjectMod with moderators #######
  {
    RandomEffecttot_Drift2 <- Ttot_DriftMeans2/Ttot_DriftWeights2; RandomEffecttot_Drift2
    RandomEffecttot_DriftVariance2 <- 1/Ttot_DriftWeights2; RandomEffecttot_DriftVariance2
    RandomEffecttot_DriftSE2 <- RandomEffecttot_DriftVariance2^.5; RandomEffecttot_DriftSE2
    RandomEffecttot_DriftUpperLimit2 <- RandomEffecttot_Drift2 + 1.96*RandomEffecttot_DriftSE2; RandomEffecttot_DriftUpperLimit2
    RandomEffecttot_DriftLowerLimit2 <- RandomEffecttot_Drift2 - 1.96*RandomEffecttot_DriftSE2; RandomEffecttot_DriftLowerLimit2
    RandomEffecttot_DriftZ2 <- RandomEffecttot_Drift2/RandomEffecttot_DriftSE2; RandomEffecttot_DriftZ2
    RandomEffecttot_DriftProb2 <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ2),
                                                       mean=c(rep(0, (n.latent2^2))), sd=c(rep(1, (n.latent2^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb2
    RandomEffecttot_DriftUpperLimitPI2 <- RandomEffecttot_Drift2 + 1.96*(tau2Drift2^.5); RandomEffecttot_DriftUpperLimitPI2
    RandomEffecttot_DriftLowerLimitPI2 <- RandomEffecttot_Drift2 - 1.96*(tau2Drift2^.5); RandomEffecttot_DriftLowerLimitPI2
    RandomEffectDriftResults2 <- rbind(RandomEffecttot_Drift2, RandomEffecttot_DriftVariance2, RandomEffecttot_DriftSE2,
                                       RandomEffecttot_DriftUpperLimit2, RandomEffecttot_DriftLowerLimit2,
                                       RandomEffecttot_DriftZ2, RandomEffecttot_DriftProb2,
                                       RandomEffecttot_DriftUpperLimitPI2, RandomEffecttot_DriftLowerLimitPI2)
    # same for dt

    ### some corrections for the output
    heterogeneity2 <- fixedEffectDriftResults2[9:16,]; heterogeneity2

    # same for dt # tbd
  }

  randomEffectDriftResults <- rbind(cbind(t(MeanOfDriftValues1), t(MeanOfDriftValues2)), cbind(t(RandomEffecttot_Drift1), t(RandomEffecttot_Drift2)),
                                    cbind(t(RandomEffecttot_DriftVariance1), t(RandomEffecttot_DriftVariance2)), cbind(t(RandomEffecttot_DriftSE1), t(RandomEffecttot_DriftSE2)),
                                    cbind(t(RandomEffecttot_DriftUpperLimit1), t(RandomEffecttot_DriftUpperLimit2)), cbind(t(RandomEffecttot_DriftLowerLimit1), t(RandomEffecttot_DriftLowerLimit2)),
                                    cbind(t(RandomEffecttot_DriftZ1), t(RandomEffecttot_DriftZ2)), cbind(t(RandomEffecttot_DriftProb1), t(RandomEffecttot_DriftProb2)),
                                    cbind(t(tau2Drift1), t(tau2Drift2)), cbind(t(Q_Drift1), t(Q_Drift2)), cbind(t(H2_Drift1), t(H2_Drift2)),
                                    cbind(t(H2DriftUpperLimit1), t(H2DriftUpperLimit2)), cbind(t(H2DriftLowerLimit1), t(H2DriftLowerLimit2)),
                                    cbind(t(I2_Drift1), t(I2_Drift2)), cbind(t(I2DriftUpperLimit1), t(I2DriftUpperLimit2)), cbind(t(I2DriftLowerLimit1), t(I2DriftLowerLimit2)))
  rownames(randomEffectDriftResults) <- c("MeanOfDriftValues", "RandomEffecttot_Drift",
                                          "RandomEffect_DriftVariance", "RandomEffect_DriftSE", "RandomEffect_DriftUpperLimit",
                                          "RandomEffect_DriftLowerLimit", "RandomEffect_DriftZ", "RandomEffect_DriftProb",
                                          "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
                                          "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")
  round(randomEffectDriftResults, 4)

  # Collect Results for both fixed effect analysis ######
  modelResultsList1 <- list(DRIFT_ctmaFitObject = DRIFTCoeff1,
                            DRIFTSE1_ctmaFitObject = DRIFTSE1); modelResultsList1
  modelResultsList2 <- list(DRIFT_ctmaFitObjectMod = DRIFTCoeff2,
                            DRIFTSE1_ctmaFitObjectMod = DRIFTSE2); modelResultsList2
  summaryList1 <- list(model="Analysis of Heterogeneity for ctmaFitObject",
                       estimates1=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults1, digits),
                                       "Heterogeneity"=round(heterogeneity1, digits),
                                       "I2 message" = fixedEffectDriftMessage,
                                       "Tau2 message" = tau2DriftMessage,
                                       "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults1, digits))); summaryList1
  summaryList2 <- list(model="Analysis of Heterogeneity for ctmaFitObjectMod",
                       estimates1=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults2, digits),
                                       "Heterogeneity"=round(heterogeneity2, digits),
                                       "I2 message" = fixedEffectDriftMessage,
                                       "Tau2 message" = tau2DriftMessage,
                                       "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults2, digits))); summaryList2

  #fixedEffectDriftResults
  FE <- round(fixedEffectDriftResults[1:8,], digits); FE
  RE <- round(randomEffectDriftResults[1:8,], digits); RE

  Het <- round(fixedEffectDriftResults[9:16,], digits); Het

  names11 <- names(ctmaFitObjectMod$modelResults$DRIFT); names11

  fitObjDriftAndSE <- round(cbind(modelResultsList1[[1]], modelResultsList1[[2]]), digits); fitObjDriftAndSE
  colnames(fitObjDriftAndSE) <- paste0("w/o Mod ", c(names11, names11)); fitObjDriftAndSE
  colnames(fitObjDriftAndSE) <- gsub("DRIFT ", "", colnames(fitObjDriftAndSE)); fitObjDriftAndSE
  colnames(fitObjDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)] <- paste0(colnames(fitObjDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)], " SE"); fitObjDriftAndSE
  #fitObjDriftAndSE

  fitObjModDriftAndSE <- round(cbind(modelResultsList2[[1]], modelResultsList2[[2]]), digits); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE) <- paste0("with Mod ", c(names11, names11)); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)] <- paste0(colnames(fitObjModDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)], " SE"); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE) <- gsub("DRIFT ", "", colnames(fitObjModDriftAndSE)); fitObjModDriftAndSE
  #fitObjModDriftAndSE

  # Analysis of Reduction in Heterogeneity by means of moderators #####
  heterogeneity1[heterogeneity1 < 0] <- .00001
  heterogeneity2[heterogeneity2 < 0] <- .00001
  #redHet <- round(heterogeneity1/heterogeneity2, digits)[c(2,3,6),]*100; redHet
  redHet <- round((heterogeneity1-heterogeneity2)/(heterogeneity1), digits)[c(2,3,6),]*100; redHet
  colnames(redHet) <- names11; redHet
  colnames(redHet) <- gsub("DRIFT ", "", colnames(redHet)); redHet
  rownames(redHet) <- gsub("_Drift2", "", rownames(redHet)); redHet
  rownames(redHet) <- paste0(rownames(redHet), ": % reduction (neg. values = % increase)"); redHet

  results <- list(activeDirectory=activeDirectory,
                  plot.type=NULL, model.type="BiG",
                  n.studies=n.studies2,
                  n.latent=n.latent2,
                  studyList=list(ctmaFitObject, ctmaFitObjectMod),
                  summary=list('fixedEffects (w/o Mod | with Mod)'=FE, 'randomEffects (w/o Mod | with Mod)'=RE,
                               fitObj_DriftAndSE=fitObjDriftAndSE, fitObjMod_DriftAndSE=fitObjModDriftAndSE,
                               'heterogeneity (w/o Mod | with Mod)'=Het,
                               HeterogeneityReduction=redHet,
                               Note="Negative I2 values were set to .00001 for computation of reduction in heterogeneity."))

  class(results) <- "CoTiMAFit"

  invisible(results)
} ### END function definition
