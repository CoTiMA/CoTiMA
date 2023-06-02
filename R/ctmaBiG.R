#' ctmaBiG
#'
#' @description Analysis of publication bias and generalizability. The function takes a CoTiMA fit object (created with \code{\link{ctmaInit}})
#' and estimates fixed and random effects of single drift coefficients, heterogeneity  (Q, I square, H square, tau square),
#' PET-PEESE corrections, Egger's tests, and z-curve analysis yielding expected replication and detection rates (ERR, EDR).
#'
#'
#' @param ctmaInitFit fit object created with \code{\link{ctmaInit}} containing the fitted ctsem model of each primary study
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param PETPEESEalpha probability level (condition) below which to switch from PET to PEESE (cf. Stanley, 2017, p. 582, below Eq. 2; default p = .10)
#' @param activateRPB if TRUE, messages (warning, finished) could be send to smart phone (default = FALSE)
#' @param digits rounding (default = 4)
#' @param zcurve performs z-curve analysis. Could fail if too few studies (e.g. around 10) are supplied. default=FALSE
#' @param undoTimeScaling if TRUE, the original time scale is used (timeScale argument possibly used in \code{\link{ctmaInit}} is undone )
#' @param dt A scalar indicating a time interval across which discrete time effects should be estimated and then used for ctmaBiG.
#'
#' @importFrom RPushbullet pbPost
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
  zcurve=FALSE,
  undoTimeScaling=TRUE,
  dt=NULL   # try BiG with dt effects. Specify the Time for which dT effects should analyzed. Always does it with original time scale
)


{  # begin function definition (until end of file)

  ### check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA (\"ctmaInitFit\") object has to be supplied to analyse something. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (ctmaInitFit$model.type == "mx") {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Your attention is required."))}
    Msg <- "I found old OpenMx-fitted Init fits. Will try ctmaBiGOMX! \n"
    message(Msg)
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
          if (dim(tmp2)[2] != dim(tmp3)[2]) { # if random intercepts were included and T0var includes manifestvar
            tmp3 <- tmp3[, c(1:n.latent, (1+n.latent):(n.latent+n.latent))]
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


    # undo time scaling
    all_Coeff_timeScaled <- all_Coeff
    all_SE_timeScaled <- all_SE
    if (undoTimeScaling) {
        if(!(is.null(ctmaInitFit$summary$scaleTime))) {
          all_Coeff <- all_Coeff * ctmaInitFit$summary$scaleTime
          all_SE <- all_SE * ctmaInitFit$summary$scaleTime
        }
    }

    # CHD 24.2.2023
    if (any(all_SE == 0)) {
      tmp <- which(all_SE == 0, arr.ind = TRUE); tmp
      colnames(tmp) <- c("Study", "Drift Effect"); tmp
      ErrorMsg <- paste0("\nAt least one SE was zero. Analysis of heterogeneity and bias cannot be performed.
                         Possible problem in fitting the single studies listed above. \nGood luck for the next try!")
      print(tmp)
      stop(ErrorMsg)
    }

    # CHD 24.2.2023
    if (!(is.null(dt))) {
      nsamples <- 1000
      drift_Coeff_dt <- matrix(NA, ncol=(1*(n.latent^2)), nrow=n.studies); drift_Coeff_dt
      drift_SE_dt <- matrix(NA, ncol=(1*(n.latent^2)), nrow=n.studies); drift_SE_dt
      #
      tmpTimeScale <- ctmaInitFit$summary$scaleTime; tmpTimeScale
      tmpTimeScale <- 1/tmpTimeScale * dt
      for (i in 1:n.studies) {
        tmpFit <- ctmaInitFit$studyFitList[[i]]
        tmpDrift_dt <- ctsem::ctStanDiscretePars(tmpFit,
                           subjects = "popmean",
                           times = tmpTimeScale,
                           nsamples = nsamples,
                           plot=FALSE,
                           indices="ALL",
                           cores='maxneeded')
        tmp <- cbind(tmpDrift_dt[ , 1, 1, , 1], tmpDrift_dt[ , 1, 1, , 2]) # columnwise
        dimnames(tmp)[[2]] <- c(matrix(names1, 2, 2, byrow=TRUE))
        drift_Coeff_dt[i, ] <- apply(tmp, 2, mean)
        drift_SE_dt[i, ] <- apply(tmp, 2, sd)
      }
    }

    #######################################################################################################################
    ##################################### Analyses of Publication Bias ####################################################
    #######################################################################################################################


    print(paste0("#################################################################################"))
    print(paste0("################################# Egger's test ##################################"))
    print(paste0("#################################################################################"))

    DRIFTCoeff <- DIFFUSIONCoeff <- T0VARCoeff <- DRIFTSE <- DIFFUSIONSE <- T0VARSE <- matrix(NA, n.studies, n.latent^2)

    DRIFTCoeff <- all_Coeff[, grep("toV", colnames(all_Coeff))]; DRIFTCoeff
    DRIFTSE <- all_SE[, grep("toV", colnames(all_SE))]; DRIFTSE
    DIFFUSIONCoeff <- all_Coeff[, grep("diff", colnames(all_Coeff))]; DIFFUSIONCoeff
    DIFFUSIONSE <- all_SE[, grep("diff", colnames(all_SE))]; DIFFUSIONSE
    tmp <- grep("T0var", colnames(all_Coeff)); tmp
    if (length(tmp) == 0) tmp <- grep("T0cov", colnames(all_Coeff)); tmp
    T0VARCoeff <- all_Coeff[, tmp]; T0VARCoeff
    tmp <- grep("T0var", colnames(all_SE)); tmp
    if (length(tmp) == 0) tmp <- grep("T0cov", colnames(all_Coeff)); tmp
    T0VARSE <- all_SE[, tmp]; T0VARSE

    # just added in case time scaled funnel and forest plots should be done with ctmaPlot
    DRIFTCoeff_timeScaled <- all_Coeff_timeScaled[, grep("toV", colnames(all_Coeff))]; DRIFTCoeff_timeScaled
    DRIFTSE_timeScaled <- all_SE_timeScaled[, grep("toV", colnames(all_SE))]; DRIFTSE_timeScaled
    DIFFUSIONCoeff_timeScaled <- all_Coeff_timeScaled[, grep("diff", colnames(all_Coeff))]; DIFFUSIONCoeff_timeScaled
    DIFFUSIONSE_timeScaled <- all_SE_timeScaled[, grep("diff", colnames(all_SE))]; DIFFUSIONSE_timeScaled


    DRIFTCoeffSND <- DRIFTCoeff / DRIFTSE; DRIFTCoeffSND
    DRIFTPrecision <- c(rep(1, n.latent^2))/(DRIFTSE); DRIFTPrecision
    colnames(DRIFTPrecision) <- colnames(DRIFTCoeffSND); DRIFTPrecision
    if (!(is.null(dt))) {
      DRIFTPrecision_dt <- c(rep(1, n.latent^2))/(drift_SE_dt); DRIFTPrecision_dt
      colnames(DRIFTPrecision_dt) <- colnames(DRIFTPrecision_dt); DRIFTPrecision_dt
    }


    message1 <- "The pos. & sign. intercept indicates that SMALLER studies produced more positive (or less negative) effects"
    message2 <- "The neg. & sign. intercept indicates that LARGER studies produced more positive (or less negative) effects"
    tmp <- c()

    eggerDrift <- list()
    for (j in 1:(n.latent^2)) {
      #j <- 1
      tmp1 <- stats::lm(DRIFTCoeffSND[,j]~DRIFTPrecision[,j]) # This is identical to a weighted regression of drift on se ...
      tmp2 <- summary(tmp1)
      eggerDrift[[j]] <- list()
      eggerDrift[[j]]$message <- "No sign. evidence for publication bias."
      eggerDrift[[j]]$summary <- tmp2[[4]]
      #if (summary(eggerDrift[[j]])$coefficients[1,1] > 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
      if (tmp2$coefficients[1,1] > 0 & tmp2$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message1
      }
      #if (summary(eggerDrift[[j]])$coefficients[1,1] < 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
      if (tmp2$coefficients[1,1] < 0 & tmp2$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message2
      }
    }

    FREAResults <- list()
    FREAResults[[1]] <- "############# Eggers Test for DRIFT Parameter Estimates  ###############################"
    FREACounter <- 1
    for (j in 1:(n.latent^2)) {
      #j <- 1
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- paste0("-------------------------------- Eggers Test for ",
                                           colnames(DRIFTCoeff)[j], "--------------------------------")
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- eggerDrift[[j]]$message
      FREACounter <- FREACounter + 1
      #FREAResults[[FREACounter]] <- summary(eggerDrift[[j]])
      FREAResults[[FREACounter]] <- eggerDrift[[j]]$summary
    }
    #str(FREAResults

    ################################### Fixed & Random Effects Analyses ###################################################

    print(paste0("#################################################################################"))
    print(paste0("########## Fixed effect analysis of each drift coefficient separately  ##########"))
    print(paste0("#################################################################################"))

    # FIXED EFFECTS ANALYSIS ###############################################################################
    DriftMeans <- colMeans(DRIFTCoeff); DriftMeans
    if (!(is.null(dt))) DriftMeans_dt <- colMeans(drift_Coeff_dt)
    # Sum of within weights  and weight * effect size
    T_DriftWeights <- colSums(DRIFTPrecision^2); T_DriftWeights
    if (!(is.null(dt))) T_DriftWeights_dt <- colSums(DRIFTPrecision_dt^2)
    #DRIFTPrecision
    T_DriftMeans <- colSums(DRIFTCoeff * DRIFTPrecision^2); T_DriftMeans
    names(T_DriftMeans) <- names(T_DriftWeights); T_DriftMeans
    if (!(is.null(dt))) {
      T_DriftMeans_dt <- colSums(drift_Coeff_dt * DRIFTPrecision_dt^2); T_DriftMeans_dt
      names(T_DriftMeans_dt) <- names(T_DriftWeights); T_DriftMeans_dt
    }
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
    # same for dt
    if (!(is.null(dt))) {
      FixedEffect_Drift_dt <- T_DriftMeans_dt/T_DriftWeights_dt; FixedEffect_Drift_dt
      FixedEffect_DriftVariance_dt <- 1/T_DriftWeights_dt; FixedEffect_DriftVariance_dt
      FixedEffect_DriftSE_dt <- FixedEffect_DriftVariance_dt^.5; FixedEffect_DriftSE_dt
      FixedEffect_DriftUpperLimit_dt <- FixedEffect_Drift_dt + 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftUpperLimit_dt
      FixedEffect_DriftLowerLimit_dt <- FixedEffect_Drift_dt - 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftLowerLimit_dt
      FixedEffect_DriftZ_dt <- FixedEffect_Drift_dt/FixedEffect_DriftSE_dt; FixedEffect_DriftZ_dt
      FixedEffect_DriftProb_dt <- round(1-stats::pnorm(abs(FixedEffect_DriftZ_dt),
                                                       mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); FixedEffect_DriftProb_dt
      Q_Drift_dt <- colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt^2)- (colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt))^2 / colSums(DRIFTPrecision_dt^2); Q_Drift_dt
      H2_Drift_dt <- Q_Drift_dt/(n.studies-1); H2_Drift_dt
      I2_Drift_dt <- (H2_Drift_dt-1)/H2_Drift_dt*100; I2_Drift_dt
    }

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
    # same for dt
    if (!(is.null(dt))) {
      T2_DriftWeights_dt <- colSums(DRIFTPrecision_dt^2^2); T2_DriftWeights_dt # Borenstein et al., 2007, p. 98
      cDrift_dt <- T_DriftWeights_dt-T2_DriftWeights_dt/T_DriftWeights_dt; cDrift_dt
      tau2Drift_dt <- (Q_Drift_dt-(n.studies-1))/cDrift_dt; tau2Drift_dt
      SElnHDrift_dt <- c()
      SElnHDrift_dt[] <- 0
      for (j in 1:(n.latent^2)) {
        if (Q_Drift_dt[j] > n.studies) SElnHDrift_dt[j] <- 1/2*(log(Q_Drift_dt[j])-log(n.studies-1))/((2*Q_Drift_dt[j])^.5-(2*(n.studies-1)-1)^.5)
        if (Q_Drift_dt[j] <= n.studies) SElnHDrift_dt[j] <-  (1/(2*(n.studies-2)) * (1 - 1/(3*(n.studies-2)^.5)) )^.5
      }
      H2DriftUpperLimit_dt <- exp(log(H2_Drift_dt) + 1.96*SElnHDrift_dt); H2DriftUpperLimit_dt
      H2DriftLowerLimit_dt <- exp(log(H2_Drift_dt) - 1.96*SElnHDrift_dt); H2DriftLowerLimit_dt
      L_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))-1.96*SElnHDrift_dt)
      U_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))+1.96*SElnHDrift_dt)
      I2DriftUpperLimit_dt <- (U_dt^2-1)/U_dt^2 * 100; I2DriftUpperLimit_dt
      I2DriftLowerLimit_dt <- (L_dt^2-1)/L_dt^2 * 100; I2DriftLowerLimit_dt
    }

    MeanOfDriftValues <- DriftMeans
    fixedEffectDriftResults <- rbind(MeanOfDriftValues, FixedEffect_Drift, FixedEffect_DriftVariance, FixedEffect_DriftSE,
                                     FixedEffect_DriftUpperLimit, FixedEffect_DriftLowerLimit,
                                     FixedEffect_DriftZ, FixedEffect_DriftProb, tau2Drift, Q_Drift, H2_Drift,
                                     H2DriftUpperLimit, H2DriftLowerLimit, I2_Drift,
                                     I2DriftUpperLimit, I2DriftLowerLimit)
    if (!(is.null(dt))) {
      fixedEffectDriftResults <- rbind(fixedEffectDriftResults,
                                       FixedEffect_Drift_dt, FixedEffect_DriftVariance_dt, FixedEffect_DriftSE_dt,
                                       FixedEffect_DriftUpperLimit_dt, FixedEffect_DriftLowerLimit_dt,
                                       FixedEffect_DriftZ_dt, FixedEffect_DriftProb_dt, tau2Drift_dt, Q_Drift_dt, H2_Drift_dt,
                                       H2DriftUpperLimit_dt, H2DriftLowerLimit_dt, I2_Drift_dt,
                                       I2DriftUpperLimit_dt, I2DriftLowerLimit_dt)
    }


    # CHD 24.2.2023
    fixedEffectDriftMessage <- c()
    if (any(I2_Drift < 0)) fixedEffectDriftMessage <- "Negative I2 values can be set to 0.0."
    tau2DriftMessage <- c()
    if (any(tau2Drift < 0)) tau2DriftMessage <- "Some tau-squared are negative. Random effects cannot be computed. Possibly a small-k-problem."

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
    # same for dt
    if (!(is.null(dt))) {
      Ttot_DriftWeights_dt <- 0
      Ttot_DriftMeans_dt <- 0
      tau2DriftExtended_dt <- do.call(rbind, replicate(n.studies, tau2Drift_dt, simplify=FALSE))
      Ttot_DriftWeights_dt <-colSums(1/ (drift_SE_dt^2 + tau2DriftExtended_dt)); Ttot_DriftWeights_dt
      Ttot_DriftMeans_dt <- colSums(drift_Coeff_dt * 1/ (drift_SE_dt^2 + tau2DriftExtended_dt)); Ttot_DriftMeans_dt
    }
    # Random effects results
    RandomEffecttot_Drift <- Ttot_DriftMeans/Ttot_DriftWeights; RandomEffecttot_Drift
    RandomEffecttot_DriftVariance <- 1/Ttot_DriftWeights; RandomEffecttot_DriftVariance
    RandomEffecttot_DriftSE <- RandomEffecttot_DriftVariance^.5; RandomEffecttot_DriftSE
    RandomEffecttot_DriftUpperLimit <- RandomEffecttot_Drift + 1.96*RandomEffecttot_DriftSE; RandomEffecttot_DriftUpperLimit
    RandomEffecttot_DriftLowerLimit <- RandomEffecttot_Drift - 1.96*RandomEffecttot_DriftSE; RandomEffecttot_DriftLowerLimit
    RandomEffecttot_DriftZ <- RandomEffecttot_Drift/RandomEffecttot_DriftSE; RandomEffecttot_DriftZ
    RandomEffecttot_DriftProb <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ),
                                                      mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb
    RandomEffecttot_DriftUpperLimitPI <- RandomEffecttot_Drift + 1.96*(tau2Drift^.5); RandomEffecttot_DriftUpperLimitPI
    RandomEffecttot_DriftLowerLimitPI <- RandomEffecttot_Drift - 1.96*(tau2Drift^.5); RandomEffecttot_DriftLowerLimitPI
    RandomEffectDriftResults <- rbind(RandomEffecttot_Drift, RandomEffecttot_DriftVariance, RandomEffecttot_DriftSE,
                                      RandomEffecttot_DriftUpperLimit, RandomEffecttot_DriftLowerLimit,
                                      RandomEffecttot_DriftZ, RandomEffecttot_DriftProb,
                                      RandomEffecttot_DriftUpperLimitPI, RandomEffecttot_DriftLowerLimitPI)
    # same for dt
    if (!(is.null(dt))) {
      RandomEffecttot_Drift_dt <- Ttot_DriftMeans_dt/Ttot_DriftWeights_dt; RandomEffecttot_Drift_dt
      RandomEffecttot_DriftVariance_dt <- 1/Ttot_DriftWeights_dt; RandomEffecttot_DriftVariance_dt
      RandomEffecttot_DriftSE_dt <- RandomEffecttot_DriftVariance_dt^.5; RandomEffecttot_DriftSE_dt
      RandomEffecttot_DriftUpperLimit_dt <- RandomEffecttot_Drift_dt + 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftUpperLimit_dt
      RandomEffecttot_DriftLowerLimit_dt <- RandomEffecttot_Drift_dt - 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftLowerLimit_dt
      RandomEffecttot_DriftZ_dt <- RandomEffecttot_Drift_dt/RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftZ_dt
      RandomEffecttot_DriftProb_dt <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ_dt),
                                                        mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb_dt
      RandomEffecttot_DriftUpperLimitPI_dt <- RandomEffecttot_Drift_dt + 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftUpperLimitPI_dt
      RandomEffecttot_DriftLowerLimitPI_dt <- RandomEffecttot_Drift_dt - 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftLowerLimitPI_dt
      RandomEffectDriftResults_dt <- rbind(RandomEffecttot_Drift_dt, RandomEffecttot_DriftVariance_dt, RandomEffecttot_DriftSE_dt,
        RandomEffecttot_DriftUpperLimit_dt, RandomEffecttot_DriftLowerLimit_dt,
        RandomEffecttot_DriftZ_dt, RandomEffecttot_DriftProb_dt,
        RandomEffecttot_DriftUpperLimitPI_dt, RandomEffecttot_DriftLowerLimitPI_dt)
    }

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
      #WLSDriftSE_fit <- summary(WLSDrift_fit[[ii]])$coefficients[2]; WLSDriftSE_fit; FixedEffect_DriftSE[[ii]] # should outperform FixedEffect_DriftSE
    }
    # same for dt
    if (!(is.null(dt))) {
      PETDrift_fit_dt <- list()
      PEESEDrift_fit_dt <- list()
      WLSDrift_fit_dt <- list()
      PET_PEESEDrift_fit_dt <- list()
      Egger2Drift_fit_dt <- list()

      #sampleSizes <- unlist(allSampleSizes); sampleSizes
      for (ii in 1:(n.latent^2)) {
        #ii <- 1
        driftCoeff <- drift_Coeff_dt[ , ii]; driftCoeff
        driftSE <- drift_SE_dt[, ii]; driftSE

        # PET; Egger's test = constant of the weigthed least squares regression of drift coefficients on their standard errors (SE) with 1/SE^2 as weights
        IV <- driftSE; IV
        DV <- driftCoeff; DV
        currentWeigths <- (1/(driftSE^2)); currentWeigths
        PETDrift_fit_dt[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PETDrift_fit_dt[[ii]]

        # Egger's Test (alternative but algebraically identical model)
        Egger2Drift_fit_dt[[ii]] <- t(c(summary(PETDrift_fit_dt[[ii]])$coefficients[2,1:4]))

        # PEESE
        IV <- driftSE^2; IV
        DV <- driftCoeff
        currentWeigths <- (1/(driftSE^2)); currentWeigths

        PEESEDrift_fit_dt[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PEESEDrift_fit_dt[[ii]]

        # PET-PEESE
        if ( (summary(PETDrift_fit_dt[[ii]]))$coefficients[1,4] > PETPEESEalpha) {
          PET_PEESEDrift_fit_dt[[ii]] <- PETDrift_fit_dt[[ii]]
        }
        if ( (summary(PETDrift_fit_dt[[ii]]))$coefficients[1,4] <= PETPEESEalpha) {
          PET_PEESEDrift_fit_dt[[ii]] <- PEESEDrift_fit_dt[[ii]]
        }

        # WLS
        DV <- driftCoeff/driftSE
        IV <- 1/driftSE
        WLSDrift_fit_dt[[ii]] <- stats::lm(DV ~ IV + 0); WLSDrift_fit_dt[[ii]]; FixedEffect_Drift_dt[[ii]] # should be identical to FixedEffect_Drift
        #WLSDriftSE_fit_dt <- summary(WLSDrift_fit_dt[[ii]])$coefficients[2]; WLSDriftSE_fit_dt; FixedEffect_DriftSE_dt[[ii]] # should outperform FixedEffect_DriftSE
      }
    }


    ############################################## zcurve Analysis ###################################################
    zFit <- zFit_dt <- "Set zcurve=TRUE for z-curve analysis results"
    if (zcurve == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("############################## Z-curve 2.0 analysis #############################"))
      print(paste0("#################################################################################"))
      print(paste0("                                                                                 "))
      print(paste0("Note: The fitting range is from |1.96| - |5|, as z-statistics > |5| have almost "))
      print(paste0("guaranteed > 99% power (the model treats them separately and includes them in   "))
      print(paste0("the EDR/ERR calculations, but does not fit the mixture of them, so too few      "))
      print(paste0("usable significant effects could perhaps be present in your data. In that case, "))
      print(paste0("you can be quite confident that the EDR and ERR will be close 1, but you have   "))
      print(paste0("to set the argument zcurve=FALSE!                                               "))

      zFit <- list()
      for (i in 1: dim(DRIFTCoeffSND)[2]) {
        tmp1 <- abs(DRIFTCoeffSND[, i]); tmp1
        tmp1
        zcurve::zcurve(z=tmp1)
        zFit[[i]] <- summary(zcurve::zcurve(z=tmp1))
      }
      names(zFit) <- paste0("Z-Curve 2.0 analysis of ", colnames(DRIFTCoeffSND)); zFit

      # same for dt
      if (!(is.null(dt))) {
        DRIFTCoeffSND_dt <- drift_Coeff_dt / drift_SE_dt; DRIFTCoeffSND_dt
        zFit_dt <- list()
        for (i in 1: dim(DRIFTCoeffSND)[2]) {
          tmp1 <- abs(DRIFTCoeffSND[, i]); tmp1
          zFit_dt[[i]] <- summary(zcurve::zcurve(z=tmp1))
        }
        names(zFit_dt) <- paste0("Z-Curve 2.0 analysis of ", colnames(DRIFTCoeffSND)); zFit_dt
      }

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
      PET_Drift <- unlist(lapply(PETDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift
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
      #PET_PEESE_DRIFTresults

      ### some corrections for the output
      heterogeneity <- fixedEffectDriftResults[9:16,]; heterogeneity
      #fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults # moved below because still needed
      eggerTest <- PET_PEESE_DRIFTresults[9:12,]; eggerTest
      PET_PEESE_DRIFTresults <- PET_PEESE_DRIFTresults[1:8,]; PET_PEESE_DRIFTresults
    }
    # same for dt
    if (!(is.null(dt))) {
      PET_Drift_dt <- unlist(lapply(PETDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift_dt
      PET_SE_dt <- c()
      for (k in 1:(n.latent^2)) PET_SE_dt <- c(PET_SE_dt, summary(PETDrift_fit_dt[[k]])$coefficients[1,2])

      PEESE_Drift_dt <-unlist(lapply(PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift_dt
      PEESE_SE_dt <- c()
      for (k in 1:(n.latent^2)) PEESE_SE_dt <- c(PEESE_SE_dt, summary(PEESEDrift_fit_dt[[k]])$coefficients[1,2])

      PET_PEESE_Drift_dt <-unlist(lapply(PET_PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift_dt
      PET_PEESE_SE_dt <- c()
      for (k in 1:(n.latent^2)) PET_PEESE_SE_dt <- c(PET_PEESE_SE_dt, summary(PET_PEESEDrift_fit_dt[[k]])$coefficients[1,2])

      WLS_Drift_dt <- unlist(lapply(WLSDrift_fit_dt, function(extract) extract$coefficients)); WLS_Drift_dt
      WLS_SE_dt <- c()
      for (k in 1:(n.latent^2)) WLS_SE_dt <- c(WLS_SE_dt, summary(WLSDrift_fit_dt[[k]])$coefficients[1,2])

      Egger2Drift_results_dt <- matrix(unlist(Egger2Drift_fit_dt), ncol=n.latent^2, nrow=4); Egger2Drift_results_dt

      PET_PEESE_DRIFTresults_dt <- rbind(PET_Drift_dt, PET_SE_dt,
                                      PEESE_Drift_dt, PEESE_SE_dt,
                                      PET_PEESE_Drift_dt, PET_PEESE_SE_dt,
                                      WLS_Drift_dt, WLS_SE_dt,
                                      Egger2Drift_results_dt)
      colnames(PET_PEESE_DRIFTresults_dt) <- colnames(DRIFTCoeff)
      rownames(PET_PEESE_DRIFTresults_dt) <- c(rownames(PET_PEESE_DRIFTresults_dt)[1:8], "Egger's b0", "SE(b0)", "T", "p")
      ### some corrections for the output
      heterogeneity_dt <- fixedEffectDriftResults[24:31,]; heterogeneity_dt # correct! fixedEffectDriftResults_dt does not yet exists
      fixedEffectDriftResults_dt <- fixedEffectDriftResults[17:23,]; fixedEffectDriftResults_dt
      eggerTest_dt <- PET_PEESE_DRIFTresults_dt[9:12,]; eggerTest_dt
      PET_PEESE_DRIFTresults_dt <- PET_PEESE_DRIFTresults_dt[1:8,]; PET_PEESE_DRIFTresults_dt
      #
      fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults
    }

    # } # End Analysis of Publication Bias
    if (is.null(dt)) {
      modelResultsList <- list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
                      DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE,
                      DRIFT_timeScaled=DRIFTCoeff_timeScaled, DIFFUSION_timeScaled=DIFFUSIONCoeff_timeScaled,
                      DRIFTSE_timeScaled=DRIFTSE_timeScaled, DIFFUSIONSE_timeScaled=DIFFUSIONSE_timeScaled)
      summaryList <- list(model="Analysis of Publication Bias & Generalizability",
                          estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
                                         "Heterogeneity"=round(heterogeneity, digits),
                                         "I2 message" = fixedEffectDriftMessage,
                                         "Tau2 message" = tau2DriftMessage,
                                         "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
                                         "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
                                         "Egger's tests"=round(eggerTest, digits),
                                         #"Egger's tests Alt. Version"= FREAResults,
                                         "Z-Curve 2.0 Results:"=zFit))
    }
    if (!(is.null(dt))) {
      modelResultsList <- list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
                           DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE,
                           DRIFT_timeScaled=DRIFTCoeff_timeScaled, DIFFUSION_timeScaled=DIFFUSIONCoeff_timeScaled,
                           DRIFTSE_timeScaled=DRIFTSE_timeScaled, DIFFUSIONSE_timeScaled=DIFFUSIONSE_timeScaled,
                           DRIFT_dt=drift_Coeff_dt,
                           DRIFT_dt_SE=drift_SE_dt)
      summaryList <- list(model="Analysis of Publication Bias & Generalizability",
                          estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
                                         "Heterogeneity"=round(heterogeneity, digits),
                                         "I2 message" = fixedEffectDriftMessage,
                                         "Tau2 message" = tau2DriftMessage,
                                         "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
                                         "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
                                         "Egger's tests"=round(eggerTest, digits),
                                         #"Egger's tests Alt. Version"= FREAResults,
                                         "Z-Curve 2.0 Results:"=zFit,
                                         "discreteTimeTimeScale" = dt,
                                         "Fixed Effects of DISCRETE TIME Drift Coefficients"=round(fixedEffectDriftResults_dt, digits),
                                         "Heterogeneity in DISCRETE TIME "=round(heterogeneity_dt, digits),
                                         "Random Effects of DISCRETE TIME Drift Coefficients"=round(RandomEffectDriftResults_dt, digits),
                                         "PET-PEESE corrections IN DISCRETE TIME"=round(PET_PEESE_DRIFTresults_dt, digits),
                                         "Egger's tests"=round(eggerTest_dt, digits),
                                         "Z-Curve 2.0 Results in DISCRTE TIME:"=zFit_dt))
      }

    results <- list(activeDirectory=activeDirectory,
                    plot.type=c("funnel", "forest"), model.type="BiG",
                    coresToUse=NULL, n.studies=n.studies,
                    n.latent=n.latent,
                    studyList=ctmaInitFit$studyList, studyFitList=NULL, # , homDRIFTallFitCI),
                    emprawList=NULL,
                    statisticsList=ctmaInitFit$statisticsList,
                    modelResults=modelResultsList,
                    parameterNames=ctmaInitFit$parameterNames,
                    summary=summaryList)

    class(results) <- "CoTiMAFit"

    invisible(results)
  }
} ### END function definition
