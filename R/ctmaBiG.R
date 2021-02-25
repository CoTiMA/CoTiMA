#' ctmaBiG
#'
#' @description Analysis of publication bias and generalizability. The function takes a CoTiMA fit object (created with ctmaInit)
#' and estimates fixed and random effects of single drift coefficients, heterogeneity  (Q, I square, H square, tau square),
#' PET-PEESE corrections, Egger's tests, and z-curve analysis yielding expected  replication and detection rates (ERR, EDR).
#'
#'
#' @param ctmaInitFit fit object created with ctmaInit containing the fitted ctsem model of each primary study
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param PETPEESEalpha probability level (condition) below which to switch from PET to PEESE (cf. Stanley, 2017, p. 582, below Eq. 2; default p = .10)
#' @param activateRPB if TRUE, messages (warning, finished) could be send to smart phone (default = FALSE)
#' @param digits rounding (default = 4)
#' @param zcurve performs z-curve analysis. Could fail if too few studies (e.g. around 10) are supplied. default=FALSE
#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red
#' @importFrom stats var lm pnorm
#' @importFrom zcurve zcurve
#'
#' @export ctmaBiG
#'
#' @examples
#' \dontrun{
#' # perform analyses of publication bias and generalizability
#' CoTiMAInitFit_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMABiG_D_BO <- ctmaBiG(ctmaInitFit=CoTiMAInitFit_D_BO, zcurve=FALSE)
#' }
#'
#' @examples
#' # display results
#' summary(CoTiMABiG_D_BO)
#'
#' @examples
#' \dontrun{
#' # get funnel & forest plots
#' CoTiMABiG_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' plot(CoTiMABiG_D_BO)
#' }
#' @return ctmaBiG returns a list containing some arguments supplied, the results of analyses of publication bias and generalizability,
#' model type, and the type of plot that could be performed with the returned object. The arguments in the returned object are activeDirectory,
#' and coresToUse. Further arguments, which are just copied from the init-fit object supplied, are, n.studies, n.latent, studyList,
#' statisticsList, modelResults (all parameter estimates and their standard error), and parameter names. All new results are returned
#' as the list element "summary", which is printed if the summary function is applied to the returned object. The summary list element
#' comprises a title (model='Analysis of Publication Bias & Generalizability') and "estimates", which is another list comprising
#' "Fixed Effects of Drift Coefficients", "Heterogeneity", "Random Effects of Drift Coefficients", "PET-PEESE corrections",
#' "Egger's tests" (constant of the WLS regression of drift coefficients on their standard errors (SE) with 1/SE^2 as weights),
#' "Egger's tests Alt. Version" (constant of the OLS regression of the standard normal deviates of the drift coefficients on their
#' precision), and "Z-Curve 2.0 Results". Plot type is plot.type=c("funnel", "forest") and model.type="BiG".
#'
ctmaBiG <- function(
  ctmaInitFit=NULL,
  activeDirectory=NULL,
  PETPEESEalpha=.10,
  activateRPB=FALSE,
  digits=4,
  zcurve=FALSE
)


{  # begin function definition (until end of file)

  ### check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA (\"ctmaInitFit\") object has to be supplied to analyse something. \n"))
    stop("Good luck for the next try!")
  }

  if (ctmaInitFit$model.type == "mx") {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Your attention is required."))}
    cat(crayon::red$bold(" I found old OpenMx-fitted Init fits. Will try ctmaBiGOMX!", sep="\n"))
    cat(crayon::red$bold(" ", " ", sep="\n"))
    result <- ctmaBiGOMX(ctmaInitFit=ctmaInitFit,
                         activeDirectory=activeDirectory,
                         PETPEESEalpha=PETPEESEalpha,
                         activateRPB=activateRPB,
                         digits=digits)

    class(results) <- "CoTiMAFit"

    return(result)
    #stop("Please check results carefully")
  }

  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  if (ctmaInitFit$model.type != "mx") {

    { # start extracting
      n.latent <- ctmaInitFit$n.latent; n.latent
      if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
      n.studies <- ctmaInitFit$n.studies; n.studies
      names1 <- names(ctmaInitFit$modelResults$DRIFT[[1]]); names1
      names2 <- names(ctmaInitFit$modelResults$DIFFUSION[[1]]); names2
      names3 <- names(ctmaInitFit$modelResults$T0VAR[[1]]); names3

      if (length(names1) != length(names2)) { # old vs. new ctsem version
        names2 <- c(OpenMx::vech2full(names2))
        names3 <- c(OpenMx::vech2full(names3))
      }

      all_Coeff <- matrix(NA, ncol=(3*(n.latent^2)), nrow=n.studies); all_Coeff
      all_SE <- matrix(NA, ncol=(3*(n.latent^2)), nrow=n.studies); all_SE

      if ("pop_T0cov" %in% names(ctmaInitFit$studyFitList[[1]]$stanfit$transformedpars)) ctsem341 <- TRUE else ctsem341 <- FALSE
      for (i in 1:n.studies) {
        if ("transformedpars" %in% names(ctmaInitFit$studyFitList[[i]]$stanfit)) {        # if maximum likelihood was used
          tmp1 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DRIFT[, , 1],
                        ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DRIFT[, , 2])
          #if (ctsem341 == TRUE) tmp1 <- tmp1[, c(1,3,2,4)]
          tmp1 <- tmp1[, c(1,3,2,4)]
          tmp2 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DIFFUSIONcov[, , 1],
                        ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_DIFFUSIONcov[, , 2])
          if (ctsem341 == TRUE) {
          tmp3 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_T0cov[, , 1],
                        ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_T0cov[, , 2])
          } else {
            tmp3 <- cbind(ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_T0VAR[, , 1],
                          ctmaInitFit$studyFitList[[i]]$stanfit$transformedpars$pop_T0VAR[, , 2])
          }
          tmp4 <- cbind(tmp1, tmp2, tmp3); tmp4
          all_Coeff[i,] <- apply(tmp4, 2, mean); all_Coeff[i,]
          all_SE[i,] <- apply(tmp4, 2, stats::sd); all_SE[i,]
        } else if ("ll" %in% names(ctmaInitFit$studyFitList[[i]]$stanfit)) { # if NUTS sampler was used
          tmp1 <- matrix(99, ncol=n.latent, nrow=n.latent); tmp1
          toSelect <- which(lower.tri(tmp1, diag=TRUE)); toSelect
          stanSummary <- summary(ctmaInitFit$studyFitList[[i]]$stanfit)
          stanSummary <- stanSummary$summary
          tmp <- grep("pop_DRIFT", rownames(stanSummary)); tmp
          tmp1 <- stanSummary[tmp, 3]; tmp1
          tmp <- grep("pop_DIFFUSIONcov", rownames(stanSummary)); tmp
          tmp2 <- stanSummary[tmp, 3][toSelect]; tmp2
          if (ctsem341 == FALSE) {
            tmp <- grep("pop_T0VAR", rownames(stanSummary))
          } else {
            tmp <- grep("pop_T0cov", rownames(stanSummary))
          }
          tmp3 <- stanSummary[tmp, 3][toSelect]; tmp3
          tmp4 <- c(tmp1, tmp2, tmp3); tmp4
          all_SE[i,] <- t(as.matrix(tmp4))
        }
      }
      colnames(all_SE) <- colnames(all_Coeff) <- c(names1, names2, names3)
      allSampleSizes <- ctmaInitFit$statisticsList$allSampleSizes; allSampleSizes
    } # end extracting

    #######################################################################################################################
    ##################################### Analyses of Publication Bias ####################################################
    #######################################################################################################################


    print(paste0("#################################################################################"))
    print(paste0("################################# Egger's test ##################################"))
    print(paste0("#################################################################################"))

    DRIFTCoeff <- DIFFUSIONCoeff <- T0VARCoeff <- DRIFTSE <- DIFFUSIONSE <- T0VARSE <- matrix(NA, n.studies, n.latent^2)

    #all_Coeff

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

    print(paste0("#################################################################################"))
    print(paste0("########## Fixed effect analysis of each drift coefficient separately  ##########"))
    print(paste0("#################################################################################"))

    # FIXED EFFECTS ANALYSIS ###############################################################################
    DriftMeans <- colMeans(DRIFTCoeff); DriftMeans
    # Sum of within weights  and weight * effect size
    T_DriftWeights <- colSums(DRIFTPrecision^2); T_DriftWeights
    #DRIFTPrecision
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
                                                  mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); FixedEffect_DriftProb
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

    print(paste0("#################################################################################"))
    print(paste0("########## Random effect analysis of each drift coefficient separately  #########"))
    print(paste0("#################################################################################"))

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

    RandomEffectDriftResults
    ### PET, PEESE & WLS approaches to correct for bias

    print(paste0("#################################################################################"))
    print(paste0("################ PET, PEESE & WLS approaches to correct for bias  ###############"))
    print(paste0("#################################################################################"))

    PETDrift_fit <- list()
    PEESEDrift_fit <- list()
    WLSDrift_fit <- list()
    PET_PEESEDrift_fit <- list()
    Egger2Drift_fit <- list()

    sampleSizes <- unlist(allSampleSizes); sampleSizes
    for (ii in 1:(n.latent^2)) {
      #ii <- 1
      driftCoeff <- DRIFTCoeff[ , ii]; driftCoeff
      driftSE <- DRIFTSE[, ii]; driftSE

      # PET; Egger's test = constant of the weigthed least squares regression of drift coefficients on their standard errors (SE) with 1/SE^2 as weights
      IV <- driftSE; IV
      DV <- driftCoeff; DV
      currentWeigths <- (1/(driftSE^2)); currentWeigths
      PETDrift_fit[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PETDrift_fit[[ii]]

      # Egger's Test (alternative but algebraically identical model)
      Egger2Drift_fit[[ii]] <- t(c(summary(PETDrift_fit[[ii]])$coefficients[2,1:4]))

      # PEESE
      IV <- driftSE^2; IV
      DV <- driftCoeff
      currentWeigths <- (1/(driftSE^2)); currentWeigths

      PEESEDrift_fit[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PEESEDrift_fit[[ii]]

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
    zFit <- "Set zcurve=TRUE for z-curve analysis results"
    if (zcurve == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("############################## Z-curve 2.0 analysis #############################"))
      print(paste0("#################################################################################"))

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
    }


    ############################################## Combine Results ###################################################
    {
      PETDrift_fit[[1]]$coefficients
      PETDrift_fit[[2]]$coefficients
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
                                                "Egger's tests Alt. Version"= FREAResults,
                                                "Z-Curve 2.0 Results:"=zFit)))
    class(results) <- "CoTiMAFit"

    invisible(results)
  }
} ### END function definition
