#' ctmaRedHet
#'
#' @description Computes the Reduction in Heterogeneity in drift effects after introducing study-level moderators
#'
#' @param activateRPB if TRUE, messages (warning, finished) could be send to smart phone (default = FALSE)
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param ctmaFitObject ctmaFit Object WITHOUT Moderators (obtained from \code{\link{ctmaFit}} with the arguments WEC=\'TRUE\' and scaleTI=FALSE)
#' @param ctmaFitObjectMod ctmaFit Object WITH Moderators (obtained from \code{\link{ctmaFit}} with the arguments WEC=\'TRUE\' and scaleTI=FALSE)
#' @param digits rounding (default = 4)
#' @param dt A vector of scalars indicating a time interval across which discrete time effects should be estimated and then used for ctmaBiG.
#' @param undoTimeScaling if TRUE, the original time scale is used (timeScale argument possibly used in \code{\link{ctmaInit}} is undone )
#'
#' @importFrom RPushbullet pbPost
#' @importFrom ctsem ctExtract
#' @importFrom OpenMx expm
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

    if ( (!(is.null(dt))) & (any(dt <= 0)) ) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA The argument \"dt\" was provided. All values must by > 0. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    scaleTime1 <- scaleTime2 <- 1
    if (ctmaFitObjectMod$argumentList$scaleTime != ctmaFitObject$argumentList$scaleTime) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      Msg <- "\nA The two fit objects provided were estimated with different scaleTime-arguments. Undoing time scaling for both fit objects is enforced!"
      message(Msg)
      scaleTime1 <- ctmaFitObjectMod$argumentList$scaleTime; scaleTime1
      scaleTime2 <- ctmaFitObject$argumentList$scaleTime; scaleTime2
    } else {
      if (undoTimeScaling == TRUE) scaleTime1 <- scaleTime2 <- ctmaFitObjectMod$argumentList$scaleTime
    }
  }

  # Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################

  start.time <- Sys.time(); start.time
  {
    driftNames1 <- driftNames2 <- ctmaFitObject$parameterNames$DRIFT; driftNames1
    n.latent1 <- n.latent2 <- ctmaFitObject$n.latent; n.latent1
    n.studies1 <- n.studies2 <- length(ctmaFitObject$studyList); n.studies1
    finishsamples1 = dim(ctmaFitObject$studyFitList$stanfit$rawposterior)[1]; finishsamples1
    #
    DRIFTCoeff1 <- matrix(unlist(ctmaFitObject$modelResults$DRIFTCoeffViaWECMean), nrow=n.studies1, byrow=T); DRIFTCoeff1
    DRIFTSE1 <- matrix(unlist(ctmaFitObject$modelResults$DRIFTCoeffViaWECSD), nrow=n.studies1, byrow=T); DRIFTSE1
    DRIFTCoeffSND1 <- DRIFTCoeff1 / DRIFTSE1; DRIFTCoeffSND1
    DRIFTPrecision1 <- c(rep(1, n.latent1^2))/(DRIFTSE1); DRIFTPrecision1
    colnames(DRIFTPrecision1) <- driftNames1; DRIFTPrecision1
    #
    DRIFTCoeff2 <- matrix(unlist(ctmaFitObjectMod$modelResults$DRIFTCoeffViaWECMean), nrow=n.studies1, byrow=T); DRIFTCoeff2
    DRIFTSE2 <- matrix(unlist(ctmaFitObjectMod$modelResults$DRIFTCoeffViaWECSD), nrow=n.studies1, byrow=T); DRIFTSE2
    DRIFTCoeffSND2 <- DRIFTCoeff2 / DRIFTSE2; DRIFTCoeffSND2
    DRIFTPrecision2 <- c(rep(1, n.latent2^2))/(DRIFTSE2); DRIFTPrecision2
    colnames(DRIFTPrecision2) <- driftNames2; DRIFTPrecision2
  }

  if (!(is.null(dt))) { # this is complex because we need to re-estimate dt effects based on raw data to obtain their SD

    ### model without moderator
    # part of the following code is taken from ctmaFit
    fitStanctModel <- ctmaFitObject$studyFitList
    n.mod.values.to.plot <- length(colnames(fitStanctModel$data$tipredsdata)[1:(n.studies1-1)])+1; n.mod.values.to.plot
    modPos <- 1:(n.studies1-1); modPos
    TIpred.values <- matrix(fitStanctModel$data$tipredsdata[, modPos], ncol=length(modPos)); head(TIpred.values)
    effectCodingWeights <- unique(TIpred.values); effectCodingWeights
    #
    tmp1 <- which(!(is.na(fitStanctModel$ctstanmodelbase$pars$param))); tmp1; length(tmp1)
    tmpPars <- fitStanctModel$ctstanmodelbase$pars[tmp1,]; tmpPars
    driftPos <- which(tmpPars$matrix == "DRIFT"); driftPos
    rawDrift <- fitStanctModel$stanfit$rawposterior[ , driftPos]; dim(rawDrift)
    #
    TIpredCols <- grep("_effect", colnames(tmpPars))[modPos]; TIpredCols
    tmp2 <- (length(tmp1)+(min(driftPos)-1)*length(TIpredCols)+1); tmp2
    TIpredPos <- tmp2:(tmp2+length(driftPos)*length(TIpredCols)-1); TIpredPos
    # resort (needed later)
    TIpredPos <- c(matrix(TIpredPos, ncol=length(TIpredCols), byrow=T)); TIpredPos
    #
    tmpNames <- paste0("Drift for Study No ", unlist(lapply(ctmaFitObject$studyList, function(x) x$originalStudyNo)), "."); tmpNames
    #
    TIpredEffTmp <- fitStanctModel$stanfit$rawposterior[,TIpredPos]; dim(TIpredEffTmp)
    studyDRIFTCoeffRaw <- list()
    for (d in 1:nrow(effectCodingWeights)) {
      multiplier <- c(matrix(effectCodingWeights[d, ], nrow=(n.latent1^2), ncol=length(effectCodingWeights[d, ]), byrow=T)); multiplier
      tmp3 <- t((multiplier) * t(TIpredEffTmp)); tmp3; dim(tmp3)
      tmp1 <- rawDrift; dim(tmp1)
      startCol <- 1
      endCol <- n.latent1^2; endCol
      for (i in 1:(n.studies1-1)) {
        tmp1 <- tmp1 + tmp3[ , startCol:endCol]
        startCol <- startCol + n.latent1^2; startCol
        endCol <- endCol + n.latent1^2; endCol
      }
      studyDRIFTCoeffRaw[[tmpNames[d]]] <- tmp1
    }
    # tform
    studyDRIFTCoeff <- studyDRIFTCoeffRaw
    pars <- fitStanctModel$ctstanmodelbase$pars
    pars <- pars[which(pars$matrix == "DRIFT"), ]
    tmp1a <- pars[, "transform"]; tmp1a
    tmp1b <- pars[, "param"]; tmp1b
    tmp1c <- which(!(is.na(tmp1b))); tmp1c
    transforms <- tmp1a[tmp1c]; transforms
    # compute tformed studyDRIFTCoeff & scale time
    for (k in 1:(length(studyDRIFTCoeffRaw))) {
      counter <- 0
      for (l in 1:(n.latent1)) {
        for (m in 1:n.latent1) {
          counter <- counter + 1
          paramTmp <- studyDRIFTCoeffRaw[[k]][, counter]; paramTmp
          for (p in 1:length(paramTmp)) {
            param <- paramTmp[p]
            studyDRIFTCoeff[[k]][p , counter] <- eval(parse(text=transforms[counter])); studyDRIFTCoeff[[k]][p, counter]
            studyDRIFTCoeff[[k]][p , counter] <- studyDRIFTCoeff[[k]][p , counter] * scaleTime1
          }
        }
      }
    }
    # compute dt estimates
    studyDRIFTCoeffDT_Mean <- studyDRIFTCoeffDT_SD <- list()
    for (k in 1:(length(studyDRIFTCoeff))) {
      studyDRIFTCoeffDT_Mean[[k]] <- studyDRIFTCoeffDT_SD[[k]] <- list()
      for (d in 1:length(dt)) {
        tmp1 <- dt[d] * studyDRIFTCoeff[[k]]
        tmp1 <- t(apply(tmp1, 1, function(x) OpenMx::expm(matrix(x, nrow=n.latent1))))
        studyDRIFTCoeffDT_Mean[[k]][[d]] <- apply(tmp1, 2, mean)
        studyDRIFTCoeffDT_SD[[k]][[d]] <- apply(tmp1, 2, sd)
      }
    }
    #studyDRIFTCoeffDT_Mean

    ### model with moderator
    fitStanctModel <- ctmaFitObjectMod$studyFitList
    n.mod.values.to.plot <- length(colnames(fitStanctModel$data$tipredsdata)[1:(n.studies2-1)])+1; n.mod.values.to.plot
    modPos <- 1:(n.studies2-1); modPos
    TIpred.values <- matrix(fitStanctModel$data$tipredsdata[, modPos], ncol=length(modPos)); head(TIpred.values)
    effectCodingWeights <- unique(TIpred.values); effectCodingWeights
    #
    tmp1 <- which(!(is.na(fitStanctModel$ctstanmodelbase$pars$param))); tmp1; length(tmp1)
    tmpPars <- fitStanctModel$ctstanmodelbase$pars[tmp1,]; tmpPars
    driftPos <- which(tmpPars$matrix == "DRIFT"); driftPos
    rawDrift <- fitStanctModel$stanfit$rawposterior[ , driftPos]; dim(rawDrift)
    # below is different from before because TIpredEffects are numbered by row and substantive moderators have to be excluded
    TIpredCols <- grep("_effect", colnames(tmpPars)); TIpredCols
    # determine positions of moderated drift
    tmpMat <- tmpPars[ , TIpredCols]; tmpMat
    counter <- 0
    toTake <- c()
    toDrop <- c()
    for (i in 1:max(driftPos)) {
      for (j in 1:ncol(tmpMat)) {
        if ( (tmpMat[i,j] == TRUE) & (i < min(driftPos)) ) counter <- c(counter, max(counter) + 1)
        if ( (tmpMat[i,j] == TRUE) & (i %in% driftPos) ) {
          counter <- c(counter, max(counter) + 1)
          toTake <- c(toTake, counter[length(counter)])
        }
        if ( (tmpMat[i,j] == TRUE) & (i %in% driftPos) & (j > max(modPos) )) {
        toDrop <- c(toDrop, counter[length(counter)])
        }
      }
    }
    counter <- counter[-1]; counter
    tmp2 <- counter[toTake[!toTake %in% toDrop]]; tmp2
    TIpredPos <- tmp2+length(tmp1); TIpredPos
    tmpNames <- paste0("Drift for Mod Study No ", unlist(lapply(ctmaFitObjectMod$studyList, function(x) x$originalStudyNo)), "."); tmpNames
    #
    TIpredEffTmp <- fitStanctModel$stanfit$rawposterior[,TIpredPos]; dim(TIpredEffTmp)
    modStudyDRIFTCoeffRaw <- list()
    for (d in 1:nrow(effectCodingWeights)) {
      multiplier <- c(matrix(effectCodingWeights[d, ], nrow=(n.latent2^2), ncol=length(effectCodingWeights[d, ]), byrow=T)); multiplier
      tmp3 <- t((multiplier) * t(TIpredEffTmp)); tmp3; dim(tmp3)
      tmp1 <- rawDrift; dim(tmp1)
      startCol <- 1
      endCol <- n.latent2^2; endCol
      for (i in 1:(n.studies2-1)) {
        tmp1 <- tmp1 + tmp3[ , startCol:endCol]
        startCol <- startCol + n.latent2^2; startCol
        endCol <- endCol + n.latent2^2; endCol
      }
      modStudyDRIFTCoeffRaw[[tmpNames[d]]] <- tmp1
    }
    # tform
    modStudyDRIFTCoeff <- modStudyDRIFTCoeffRaw
    pars <- fitStanctModel$ctstanmodelbase$pars
    pars <- pars[which(pars$matrix == "DRIFT"), ]
    tmp1a <- pars[, "transform"]; tmp1a
    tmp1b <- pars[, "param"]; tmp1b
    tmp1c <- which(!(is.na(tmp1b))); tmp1c
    transforms <- tmp1a[tmp1c]; transforms
    # compute tformed studyDRIFTCoeff and scale time
    for (k in 1:(length(modStudyDRIFTCoeffRaw))) {
      counter <- 0
      for (l in 1:(n.latent2)) {
        for (m in 1:n.latent2) {
          counter <- counter + 1
          paramTmp <- modStudyDRIFTCoeffRaw[[k]][, counter]; paramTmp
          for (p in 1:length(paramTmp)) {
            param <- paramTmp[p]
            modStudyDRIFTCoeff[[k]][p , counter] <- eval(parse(text=transforms[counter])); modStudyDRIFTCoeff[[k]][p, counter]
            modStudyDRIFTCoeff[[k]][p , counter] <- modStudyDRIFTCoeff[[k]][p , counter] * scaleTime1
          }
        }
      }
    }
    # compute dt estimates
    modStudyDRIFTCoeffDT_Mean <- modStudyDRIFTCoeffDT_SD <- list()
    for (k in 1:(length(modStudyDRIFTCoeff))) {
      modStudyDRIFTCoeffDT_Mean[[k]] <- modStudyDRIFTCoeffDT_SD[[k]] <- list()
      for (d in 1:length(dt)) {
        tmp1 <- dt[d] * modStudyDRIFTCoeff[[k]]
        tmp1 <- t(apply(tmp1, 1, function(x) OpenMx::expm(matrix(x, nrow=n.latent1))))
        modStudyDRIFTCoeffDT_Mean[[k]][[d]] <- apply(tmp1, 2, mean)
        modStudyDRIFTCoeffDT_SD[[k]][[d]] <- apply(tmp1, 2, sd)
      }
    }
    # compute estimates required for computing fixed & random effects and heterogeneity
    DRIFTCoeff1DT <- DRIFTSE1DT <- DRIFTCoeffSND1DT <- DRIFTPrecision1DT <- list()
    DRIFTCoeff2DT <- DRIFTSE2DT <- DRIFTCoeffSND2DT <- DRIFTPrecision2DT <- list()
    for (i in 1:length(dt)) {
      DRIFTCoeff1DT[[i]] <- matrix(unlist(lapply(studyDRIFTCoeffDT_Mean, function(x) x[[i]])), nrow=n.studies1, byrow=T)
      DRIFTSE1DT[[i]] <- matrix(unlist(lapply(studyDRIFTCoeffDT_SD, function(x) x[[i]])), nrow=n.studies1, byrow=T)
      DRIFTCoeffSND1DT[[i]] <- DRIFTCoeff1DT[[i]] / DRIFTSE1DT[[i]]; DRIFTCoeffSND1DT[[i]]
      DRIFTPrecision1DT[[i]] <- c(rep(1, n.latent1^2))/(DRIFTSE1DT[[i]]); DRIFTPrecision1DT[[i]]
      colnames(DRIFTPrecision1DT[[i]]) <- driftNames1; DRIFTPrecision1DT[[i]]
      #
      DRIFTCoeff2DT[[i]] <- matrix(unlist(lapply(modStudyDRIFTCoeffDT_Mean, function(x) x[[i]])), nrow=n.studies2, byrow=T)
      DRIFTSE2DT[[i]] <- matrix(unlist(lapply(modStudyDRIFTCoeffDT_SD, function(x) x[[i]])), nrow=n.studies2, byrow=T)
      DRIFTCoeffSND2DT[[i]] <- DRIFTCoeff2DT[[i]] / DRIFTSE2DT[[i]]; DRIFTCoeffSND2DT[[i]]
      DRIFTPrecision2DT[[i]] <- c(rep(1, n.latent2^2))/(DRIFTSE2DT[[i]]); DRIFTPrecision2DT[[i]]
      colnames(DRIFTPrecision2DT[[i]]) <- driftNames2; DRIFTPrecision2DT[[i]]
    }
    names(DRIFTCoeff1DT) <- names(DRIFTSE1DT) <- paste0("Time Interval = ", dt)
    names(DRIFTCoeff2DT) <- names(DRIFTSE2DT) <- paste0("Time Interval = ", dt)
  } # end if (!(is.null(dt)))

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

    ### Fixed effects results for ctmaFitObjectMod without moderators ######

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
  }

    # same for dt
    if (!(is.null(dt))) {
      DriftMeans1DT <- T_DriftWeights1DT<- T_DriftMeans1DT <- list()
      DriftMeans2DT <- T_DriftWeights2DT <- T_DriftMeans2DT <- list()

      FixedEffect_Drift1DT <- FixedEffect_DriftVariance1DT <- FixedEffect_DriftSE1DT <- list()
      FixedEffect_DriftUpperLimit1DT <- FixedEffect_DriftLowerLimit1DT <- FixedEffect_DriftZ1DT <- list()
      FixedEffect_DriftProb1DT <- Q_Drift1DT <- H2_Drift1DT <- I2_Drift1DT <- list()
      T2_DriftWeights1DT <- tau2Drift1DT <- list()
      H2DriftUpperLimit1DT <- H2DriftLowerLimit1DT <- I2DriftUpperLimit1DT <- I2DriftLowerLimit1DT <- list()
      for (i in 1:length(dt)) {
        ##
        DriftMeans1DT[[i]] <- colMeans(DRIFTCoeff1DT[[i]]); DriftMeans1DT[[i]]
        # Sum of within weights  and weight * effect size
        T_DriftWeights1DT[[i]] <- colSums(DRIFTPrecision1DT[[i]]^2); T_DriftWeights1DT[[i]]
        T_DriftMeans1DT[[i]] <- colSums(DRIFTCoeff1DT[[i]] * DRIFTPrecision1DT[[i]]^2); T_DriftMeans1DT[[i]]
        names(T_DriftMeans1DT[[i]]) <- names(T_DriftWeights1DT[[i]]); T_DriftMeans1DT[[i]]
        #
        DriftMeans2DT[[i]] <- colMeans(DRIFTCoeff2DT[[i]]); DriftMeans2DT[[i]]
        # Sum of within weights  and weight * effect size
        T_DriftWeights2DT[[i]] <- colSums(DRIFTPrecision2DT[[i]]^2); T_DriftWeights2DT[[i]]
        T_DriftMeans2DT[[i]] <- colSums(DRIFTCoeff2DT[[i]] * DRIFTPrecision2DT[[i]]^2); T_DriftMeans2DT[[i]]
        names(T_DriftMeans2DT[[i]]) <- names(T_DriftWeights2DT[[i]]); T_DriftMeans2DT[[i]]
        #
        ### Fixed effects results for ctmaFitObjectMod without moderators ######
        #
        FixedEffect_Drift1DT[[i]] <- T_DriftMeans1DT[[i]]/T_DriftWeights1DT[[i]]; FixedEffect_Drift1DT[[i]]
        FixedEffect_DriftVariance1DT[[i]] <- 1/T_DriftWeights1DT[[i]]; FixedEffect_DriftVariance1DT[[i]]
        FixedEffect_DriftSE1DT[[i]] <- FixedEffect_DriftVariance1DT[[i]]^.5; FixedEffect_DriftSE1DT[[i]]
        FixedEffect_DriftUpperLimit1DT[[i]] <- FixedEffect_Drift1DT[[i]] + 1.96*FixedEffect_DriftSE1DT[[i]]; FixedEffect_DriftUpperLimit1DT[[i]]
        FixedEffect_DriftLowerLimit1DT[[i]] <- FixedEffect_Drift1DT[[i]] - 1.96*FixedEffect_DriftSE1DT[[i]]; FixedEffect_DriftLowerLimit1DT[[i]]
        FixedEffect_DriftZ1DT[[i]] <- FixedEffect_Drift1DT[[i]]/FixedEffect_DriftSE1DT[[i]]; FixedEffect_DriftZ1DT[[i]]
        FixedEffect_DriftProb1DT[[i]] <- round(1-stats::pnorm(abs(FixedEffect_DriftZ1DT[[i]]),
                                                       mean=c(rep(0, (n.latent1^2))), sd=c(rep(1, (n.latent1^2))), log.p=F), digits=digits); FixedEffect_DriftProb1DT[[i]]
        Q_Drift1DT[[i]] <- colSums(DRIFTPrecision1DT[[i]]^2 * DRIFTCoeff1DT[[i]]^2)-
          (colSums(DRIFTPrecision1DT[[i]]^2 * DRIFTCoeff1DT[[i]]))^2 / colSums(DRIFTPrecision1DT[[i]]^2); Q_Drift1DT[[i]]
        H2_Drift1DT[[i]] <- Q_Drift1DT[[i]]/(n.studies1-1); H2_Drift1DT[[i]]
        I2_Drift1DT[[i]] <- (H2_Drift1DT[[i]]-1)/H2_Drift1DT[[i]]*100; I2_Drift1DT[[i]]
        # Tau square
        T2_DriftWeights1DT[[i]] <- colSums(DRIFTPrecision1DT[[i]]^2^2); T2_DriftWeights1DT[[i]] # Borenstein et al., 2007, p. 98
        cDrift1 <- T_DriftWeights1DT[[i]] - T2_DriftWeights1DT[[i]]/T_DriftWeights1DT[[i]]; cDrift1
        tau2Drift1DT[[i]] <- (Q_Drift1DT[[i]]  - (n.studies1-1))/cDrift1; tau2Drift1DT[[i]]
        SElnHDrift1 <- c()
        SElnHDrift1[] <- 0
        for (j in 1:(n.latent1^2)) {
          if (Q_Drift1DT[[i]][j] > n.studies1) SElnHDrift1[j] <- 1/2*(log(Q_Drift1DT[[i]][j])-log(n.studies1-1))/((2*Q_Drift1DT[[i]][j])^.5-(2*(n.studies1-1)-1)^.5)
          if (Q_Drift1DT[[i]][j] <= n.studies1) SElnHDrift1[j] <-  (1/(2*(n.studies1-2)) * (1 - 1/(3*(n.studies1-2)^.5)) )^.5
        }
        H2DriftUpperLimit1DT[[i]] <- exp(log(H2_Drift1DT[[i]]) + 1.96*SElnHDrift1); H2DriftUpperLimit1DT[[i]]
        H2DriftLowerLimit1DT[[i]] <- exp(log(H2_Drift1DT[[i]]) - 1.96*SElnHDrift1); H2DriftLowerLimit1DT[[i]]
        L1 <- exp(0.5*log(Q_Drift1DT[[i]]/(n.studies1-1))-1.96*SElnHDrift1)
        U1 <- exp(0.5*log(Q_Drift1DT[[i]]/(n.studies1-1))+1.96*SElnHDrift1)
        I2DriftUpperLimit1DT[[i]] <- (U1^2-1)/U1^2 * 100; I2DriftUpperLimit1DT[[i]]
        I2DriftLowerLimit1DT[[i]] <- (L1^2-1)/L1^2 * 100; I2DriftLowerLimit1DT[[i]]
      } # end for (i in 1:length(dt))
    } # end if (!(is.null(dt)))


  ## now all with moderators ################################################################


  ### Fixed effects results for ctmaFitObjectMod with moderators ######
  {
    FixedEffect_Drift2 <- T_DriftMeans2/T_DriftWeights2; FixedEffect_Drift2
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
  }

  # same for dt
  if (!(is.null(dt))) {
    FixedEffect_Drift2DT <- FixedEffect_DriftVariance2DT <- FixedEffect_DriftSE2DT <- list()
    FixedEffect_DriftUpperLimit2DT <- FixedEffect_DriftLowerLimit2DT <- FixedEffect_DriftZ2DT <- list()
    FixedEffect_DriftProb2DT <- Q_Drift2DT <- H2_Drift2DT <- I2_Drift2DT <- list()
    T2_DriftWeights2DT <- tau2Drift2DT <- list()
    H2DriftUpperLimit2DT <- H2DriftLowerLimit2DT <- I2DriftUpperLimit2DT <- I2DriftLowerLimit2DT <- list()
    for (i in 1:length(dt)) {
      FixedEffect_Drift2DT[[i]] <- T_DriftMeans2DT[[i]]/T_DriftWeights2DT[[i]]; FixedEffect_Drift2DT[[i]]
      FixedEffect_DriftVariance2DT[[i]] <- 1/T_DriftWeights2DT[[i]]; FixedEffect_DriftVariance2DT[[i]]
      FixedEffect_DriftSE2DT[[i]] <- FixedEffect_DriftVariance2DT[[i]]^.5; FixedEffect_DriftSE2DT[[i]]
      FixedEffect_DriftUpperLimit2DT[[i]] <- FixedEffect_Drift2DT[[i]] + 1.96*FixedEffect_DriftSE2DT[[i]]; FixedEffect_DriftUpperLimit2DT[[i]]
      FixedEffect_DriftLowerLimit2DT[[i]] <- FixedEffect_Drift2DT[[i]] - 1.96*FixedEffect_DriftSE2DT[[i]]; FixedEffect_DriftLowerLimit2DT[[i]]
      FixedEffect_DriftZ2DT[[i]] <- FixedEffect_Drift2DT[[i]]/FixedEffect_DriftSE2DT[[i]]; FixedEffect_DriftZ2DT[[i]]
      FixedEffect_DriftProb2DT[[i]] <- round(1-stats::pnorm(abs(FixedEffect_DriftZ2DT[[i]]),
                                                            mean=c(rep(0, (n.latent1^2))), sd=c(rep(1, (n.latent1^2))), log.p=F), digits=digits); FixedEffect_DriftProb2DT[[i]]
      Q_Drift2DT[[i]] <- colSums(DRIFTPrecision2DT[[i]]^2 * DRIFTCoeff2DT[[i]]^2)-
        (colSums(DRIFTPrecision2DT[[i]]^2 * DRIFTCoeff2DT[[i]]))^2 / colSums(DRIFTPrecision2DT[[i]]^2); Q_Drift2DT[[i]]
      H2_Drift2DT[[i]] <- Q_Drift2DT[[i]]/(n.studies1-1); H2_Drift2DT[[i]]
      I2_Drift2DT[[i]] <- (H2_Drift2DT[[i]]-1)/H2_Drift2DT[[i]]*100; I2_Drift2DT[[i]]
      # Tau square
      T2_DriftWeights2DT[[i]] <- colSums(DRIFTPrecision2DT[[i]]^2^2); T2_DriftWeights2DT[[i]] # Borenstein et al., 2007, p. 98
      cDrift1 <- T_DriftWeights2DT[[i]] - T2_DriftWeights2DT[[i]]/T_DriftWeights2DT[[i]]; cDrift1
      tau2Drift2DT[[i]] <- (Q_Drift2DT[[i]]  - (n.studies1-1))/cDrift1; tau2Drift2DT[[i]]
      SElnHDrift1 <- c()
      SElnHDrift1[] <- 0
      for (j in 1:(n.latent1^2)) {
        if (Q_Drift2DT[[i]][j] > n.studies1) SElnHDrift1[j] <- 1/2*(log(Q_Drift2DT[[i]][j])-log(n.studies1-1))/((2*Q_Drift2DT[[i]][j])^.5-(2*(n.studies1-1)-1)^.5)
        if (Q_Drift2DT[[i]][j] <= n.studies1) SElnHDrift1[j] <-  (1/(2*(n.studies1-2)) * (1 - 1/(3*(n.studies1-2)^.5)) )^.5
      }
      H2DriftUpperLimit2DT[[i]] <- exp(log(H2_Drift2DT[[i]]) + 1.96*SElnHDrift1); H2DriftUpperLimit2DT[[i]]
      H2DriftLowerLimit2DT[[i]] <- exp(log(H2_Drift2DT[[i]]) - 1.96*SElnHDrift1); H2DriftLowerLimit2DT[[i]]
      L1 <- exp(0.5*log(Q_Drift2DT[[i]]/(n.studies1-1))-1.96*SElnHDrift1)
      U1 <- exp(0.5*log(Q_Drift2DT[[i]]/(n.studies1-1))+1.96*SElnHDrift1)
      I2DriftUpperLimit2DT[[i]] <- (U1^2-1)/U1^2 * 100; I2DriftUpperLimit2DT[[i]]
      I2DriftLowerLimit2DT[[i]] <- (L1^2-1)/L1^2 * 100; I2DriftLowerLimit2DT[[i]]
    } # end for (i in 1:length(dt))
  } # end if (!(is.null(dt)))


  fixedEffectDriftResults <- rbind(#cbind(t(MeanOfDriftValues1), t(MeanOfDriftValues2)),
    cbind(t(FixedEffect_Drift1), t(FixedEffect_Drift2)),
    cbind(t(FixedEffect_DriftVariance1), t(FixedEffect_DriftVariance2)), cbind(t(FixedEffect_DriftSE1), t(FixedEffect_DriftSE2)),
    cbind(t(FixedEffect_DriftUpperLimit1), t(FixedEffect_DriftUpperLimit2)), cbind(t(FixedEffect_DriftLowerLimit1), t(FixedEffect_DriftLowerLimit2)),
    cbind(t(FixedEffect_DriftZ1), t(FixedEffect_DriftZ2)), cbind(t(FixedEffect_DriftProb1), t(FixedEffect_DriftProb2)),
    cbind(t(tau2Drift1), t(tau2Drift2)), cbind(t(Q_Drift1), t(Q_Drift2)), cbind(t(H2_Drift1), t(H2_Drift2)),
    cbind(t(H2DriftUpperLimit1), t(H2DriftUpperLimit2)), cbind(t(H2DriftLowerLimit1), t(H2DriftLowerLimit2)),
    cbind(t(I2_Drift1), t(I2_Drift2)), cbind(t(I2DriftUpperLimit1), t(I2DriftUpperLimit2)), cbind(t(I2DriftLowerLimit1), t(I2DriftLowerLimit2)))
  rownames(fixedEffectDriftResults) <- c(#"MeanOfDriftValues",
    "FixedEffect_Drift",
    "FixedEffect_DriftVariance", "FixedEffect_DriftSE", "FixedEffect_DriftUpperLimit",
    "FixedEffect_DriftLowerLimit", "FixedEffect_DriftZ", "FixedEffect_DriftProb",
    "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
    "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")

  fixedEffectDriftResults1 <- rbind(#MeanOfDriftValues1,
    FixedEffect_Drift1, FixedEffect_DriftVariance1, FixedEffect_DriftSE1,
    FixedEffect_DriftUpperLimit1, FixedEffect_DriftLowerLimit1,
    FixedEffect_DriftZ1, FixedEffect_DriftProb1, tau2Drift1, Q_Drift1, H2_Drift1,
    H2DriftUpperLimit1, H2DriftLowerLimit1, I2_Drift1,
    I2DriftUpperLimit1, I2DriftLowerLimit1)
  fixedEffectDriftResults2 <- rbind(#MeanOfDriftValues2,
    FixedEffect_Drift2, FixedEffect_DriftVariance2, FixedEffect_DriftSE2,
    FixedEffect_DriftUpperLimit2, FixedEffect_DriftLowerLimit2,
    FixedEffect_DriftZ2, FixedEffect_DriftProb2, tau2Drift2, Q_Drift2, H2_Drift2,
    H2DriftUpperLimit2, H2DriftLowerLimit2, I2_Drift2,
    I2DriftUpperLimit2, I2DriftLowerLimit2)

  fixedEffectDriftMessage <- c()
  if ( (any(I2_Drift1 < 0)) | (any(I2_Drift2 < 0)) ) fixedEffectDriftMessage <- "Negative I2 values can be set to 0.0."
  tau2DriftMessage <- c()
  if ( (any(tau2Drift1 < 0)) | (any(tau2Drift1 < 0)) ) tau2DriftMessage <- "Some tau-squared are negative. Random effects cannot be computed. Possibly a small-k-problem."

  if (!(is.null(dt))) {
    fixedEffectDriftResultsDT <- fixedEffectDriftResults1DT <- fixedEffectDriftResults2DT <- list()
    for (i in 1:length(dt)) {
      #i <- 1
      fixedEffectDriftResultsDT[[i]] <- rbind(#cbind(t(MeanOfDriftValues1), t(MeanOfDriftValues2)),
        cbind(t(FixedEffect_Drift1DT[[i]]), t(FixedEffect_Drift2DT[[i]])),
        cbind(t(FixedEffect_DriftVariance1DT[[i]]), t(FixedEffect_DriftVariance2DT[[i]])), cbind(t(FixedEffect_DriftSE1DT[[i]]), t(FixedEffect_DriftSE2DT[[i]])),
        cbind(t(FixedEffect_DriftUpperLimit1DT[[i]]), t(FixedEffect_DriftUpperLimit2DT[[i]])), cbind(t(FixedEffect_DriftLowerLimit1DT[[i]]), t(FixedEffect_DriftLowerLimit2DT[[i]])),
        cbind(t(FixedEffect_DriftZ1DT[[i]]), t(FixedEffect_DriftZ2DT[[i]])), cbind(t(FixedEffect_DriftProb1DT[[i]]), t(FixedEffect_DriftProb2DT[[i]])),
        cbind(t(tau2Drift1DT[[i]]), t(tau2Drift2DT[[i]])), cbind(t(Q_Drift1DT[[i]]), t(Q_Drift2DT[[i]])), cbind(t(H2_Drift1DT[[i]]), t(H2_Drift2DT[[i]])),
        cbind(t(H2DriftUpperLimit1DT[[i]]), t(H2DriftUpperLimit2DT[[i]])), cbind(t(H2DriftLowerLimit1DT[[i]]), t(H2DriftLowerLimit2DT[[i]])),
        cbind(t(I2_Drift1DT[[i]]), t(I2_Drift2DT[[i]])), cbind(t(I2DriftUpperLimit1DT[[i]]), t(I2DriftUpperLimit2DT[[i]])),
        cbind(t(I2DriftLowerLimit1DT[[i]]), t(I2DriftLowerLimit2DT[[i]])))
      rownames(fixedEffectDriftResultsDT[[i]]) <- c(#"MeanOfDriftValues",
        "FixedEffect_Drift",
        "FixedEffect_DriftVariance", "FixedEffect_DriftSE", "FixedEffect_DriftUpperLimit",
        "FixedEffect_DriftLowerLimit", "FixedEffect_DriftZ", "FixedEffect_DriftProb",
        "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
        "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")

      fixedEffectDriftResults1DT[[i]] <- rbind(#MeanOfDriftValues1,
        FixedEffect_Drift1DT[[i]],  FixedEffect_DriftVariance1DT[[i]],  FixedEffect_DriftSE1DT[[i]],
        FixedEffect_DriftUpperLimit1DT[[i]],  FixedEffect_DriftLowerLimit1DT[[i]],
        FixedEffect_DriftZ1DT[[i]],  FixedEffect_DriftProb1DT[[i]],  tau2Drift1DT[[i]],  Q_Drift1DT[[i]],  H2_Drift1DT[[i]],
        H2DriftUpperLimit1DT[[i]],  H2DriftLowerLimit1DT[[i]],  I2_Drift1DT[[i]],
        I2DriftUpperLimit1DT[[i]],  I2DriftLowerLimit1)
      rownames(fixedEffectDriftResults1DT[[i]]) <- c(#"MeanOfDriftValues",
        "FixedEffect_Drift",
        "FixedEffect_DriftVariance", "FixedEffect_DriftSE", "FixedEffect_DriftUpperLimit",
        "FixedEffect_DriftLowerLimit", "FixedEffect_DriftZ", "FixedEffect_DriftProb",
        "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
        "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")

      fixedEffectDriftResults2DT[[i]] <- rbind(#MeanOfDriftValues2DT[[i]],
        FixedEffect_Drift2DT[[i]],  FixedEffect_DriftVariance2DT[[i]],  FixedEffect_DriftSE2DT[[i]],
        FixedEffect_DriftUpperLimit2DT[[i]],  FixedEffect_DriftLowerLimit2DT[[i]],
        FixedEffect_DriftZ2DT[[i]],  FixedEffect_DriftProb2DT[[i]],  tau2Drift2DT[[i]],  Q_Drift2DT[[i]],  H2_Drift2DT[[i]],
        H2DriftUpperLimit2DT[[i]],  H2DriftLowerLimit2DT[[i]],  I2_Drift2DT[[i]],
        I2DriftUpperLimit2DT[[i]],  I2DriftLowerLimit2)
      rownames(fixedEffectDriftResults2DT[[i]]) <- c(#"MeanOfDriftValues",
        "FixedEffect_Drift",
        "FixedEffect_DriftVariance", "FixedEffect_DriftSE", "FixedEffect_DriftUpperLimit",
        "FixedEffect_DriftLowerLimit", "FixedEffect_DriftZ", "FixedEffect_DriftProb",
        "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
        "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")

    } # end for (i in 1:length(dt))
    #lapply(fixedEffectDriftResultsDT, function(x) round(x, digits))
    names(fixedEffectDriftResultsDT) <-  paste0("Time Interval = ", dt)
    names(fixedEffectDriftResults1DT) <-  paste0("Time Interval = ", dt)
    names(fixedEffectDriftResults2DT) <-  paste0("Time Interval = ", dt)
  } # end  if (!(is.null(dt)))


  ## RANDOM EFFECTS ANALYSIS ##############################################################################

  # Total variance weighting
  tau2DriftExtended1 <- do.call(rbind, replicate(n.studies1, tau2Drift1, simplify=FALSE))
  Ttot_DriftWeights1 <-colSums(1/ (DRIFTSE1^2 + tau2DriftExtended1)); Ttot_DriftWeights1
  Ttot_DriftMeans1 <- colSums(DRIFTCoeff1 * 1/ (DRIFTSE1^2 + tau2DriftExtended1)); Ttot_DriftMeans1

  tau2DriftExtended2 <- do.call(rbind, replicate(n.studies2, tau2Drift2, simplify=FALSE))
  Ttot_DriftWeights2 <-colSums(1/ (DRIFTSE2^2 + tau2DriftExtended2)); Ttot_DriftWeights2
  Ttot_DriftMeans2 <- colSums(DRIFTCoeff2 * 1/ (DRIFTSE2^2 + tau2DriftExtended2)); Ttot_DriftMeans2

  # same for dt
  if (!(is.null(dt))) {
    tau2DriftExtended1DT <- Ttot_DriftWeights1DT <- Ttot_DriftMeans1DT <- list()
    tau2DriftExtended2DT <- Ttot_DriftWeights2DT <- Ttot_DriftMeans2DT <- list()
    for (i in 1:length(dt)) {
      tau2DriftExtended1DT[[i]] <- do.call(rbind, replicate(n.studies1, tau2Drift1DT[[i]], simplify=FALSE))
      Ttot_DriftWeights1DT[[i]] <-colSums(1/ (DRIFTSE1DT[[i]]^2 + tau2DriftExtended1DT[[i]])); Ttot_DriftWeights1DT[[i]]
      Ttot_DriftMeans1DT[[i]] <- colSums(DRIFTCoeff1DT[[i]] * 1/ (DRIFTSE1DT[[i]]^2 + tau2DriftExtended1DT[[i]])); Ttot_DriftMeans1DT[[i]]

      tau2DriftExtended2DT[[i]] <- do.call(rbind, replicate(n.studies1, tau2Drift2DT[[i]], simplify=FALSE))
      Ttot_DriftWeights2DT[[i]] <-colSums(1/ (DRIFTSE2DT[[i]]^2 + tau2DriftExtended2DT[[i]])); Ttot_DriftWeights2DT[[i]]
      Ttot_DriftMeans2DT[[i]] <- colSums(DRIFTCoeff2DT[[i]] * 1/ (DRIFTSE2DT[[i]]^2 + tau2DriftExtended2DT[[i]])); Ttot_DriftMeans2DT[[i]]
    }
  }

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
    rownames(RandomEffectDriftResults1) <- c(#"MeanOfDriftValues",
      "RandomEffect_Drift", "RandomEffect_DriftVariance", "RandomEffect_DriftSE",
      "RandomEffect_DriftUpperLimit", "RandomEffect_DriftLowerLimit",
      "RandomEffect_DriftZ", "RandomEffect_DriftProb",
      "H2DriftUpperLimit", "H2DriftLowerLimit")

     ### some corrections for the output
    heterogeneity1 <- fixedEffectDriftResults1[8:15,]; heterogeneity1
  }

    # same for dt
    if (!(is.null(dt))) {
      RandomEffecttot_Drift1DT <- RandomEffecttot_DriftVariance1DT <- RandomEffecttot_DriftSE1DT <- list()
      RandomEffecttot_DriftUpperLimit1DT <- RandomEffecttot_DriftLowerLimit1DT <- RandomEffecttot_DriftProb12DT <- list()
      RandomEffecttot_DriftUpperLimitPI1DT <- RandomEffecttot_DriftLowerLimitPI1DT <- RandomEffectDriftResults1DT <- list()
      RandomEffecttot_DriftZ1DT <- RandomEffecttot_DriftProb1DT <- list()
      RandomEffecttot_DriftUpperLimitPI1DT <- RandomEffecttot_DriftLowerLimitPI1DT <- RandomEffectDriftResults1DT <- list()
      heterogeneity1DT <- list()
      for (i in 1:length(dt)) {
        Ttot_DriftMeans1[[i]]
        Ttot_DriftWeights1[[i]]
        RandomEffecttot_Drift1DT[[i]] <- Ttot_DriftMeans1DT[[i]]/Ttot_DriftWeights1DT[[i]]; RandomEffecttot_Drift1DT[[i]]
        RandomEffecttot_DriftVariance1DT[[i]] <- 1/Ttot_DriftWeights1DT[[i]]; RandomEffecttot_DriftVariance1DT[[i]]
        RandomEffecttot_DriftSE1DT[[i]] <- RandomEffecttot_DriftVariance1DT[[i]]^.5; RandomEffecttot_DriftSE1DT[[i]]
        RandomEffecttot_DriftUpperLimit1DT[[i]] <- RandomEffecttot_Drift1DT[[i]] + 1.96*RandomEffecttot_DriftSE1DT[[i]]; RandomEffecttot_DriftUpperLimit1DT[[i]]
        RandomEffecttot_DriftLowerLimit1DT[[i]] <- RandomEffecttot_Drift1DT[[i]] - 1.96*RandomEffecttot_DriftSE1DT[[i]]; RandomEffecttot_DriftLowerLimit1DT[[i]]
        RandomEffecttot_DriftZ1DT[[i]] <- RandomEffecttot_Drift1DT[[i]]/RandomEffecttot_DriftSE1DT[[i]]; RandomEffecttot_DriftZ1DT[[i]]
        RandomEffecttot_DriftProb1DT[[i]] <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ1DT[[i]]),
                                                                  mean=c(rep(0, (n.latent1^2))), sd=c(rep(1, (n.latent1^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb1DT[[i]]
        RandomEffecttot_DriftUpperLimitPI1DT[[i]] <- RandomEffecttot_Drift1DT[[i]] + 1.96*(tau2Drift1DT[[i]]^.5); RandomEffecttot_DriftUpperLimitPI1DT[[i]]
        RandomEffecttot_DriftLowerLimitPI1DT[[i]] <- RandomEffecttot_Drift1DT[[i]] - 1.96*(tau2Drift1DT[[i]]^.5); RandomEffecttot_DriftLowerLimitPI1DT[[i]]
        RandomEffectDriftResults1DT[[i]] <- rbind(RandomEffecttot_Drift1DT[[i]], RandomEffecttot_DriftVariance1DT[[i]], RandomEffecttot_DriftSE1DT[[i]],
                                                  RandomEffecttot_DriftUpperLimit1DT[[i]], RandomEffecttot_DriftLowerLimit1DT[[i]],
                                                  RandomEffecttot_DriftZ1DT[[i]], RandomEffecttot_DriftProb1DT[[i]],
                                                  RandomEffecttot_DriftUpperLimitPI1DT[[i]], RandomEffecttot_DriftLowerLimitPI1DT[[i]])
        rownames(RandomEffectDriftResults1DT[[i]]) <- c(#"MeanOfDriftValues",
          "RandomEffect_Drift", "RandomEffect_DriftVariance", "RandomEffect_DriftSE",
          "RandomEffect_DriftUpperLimit", "RandomEffect_DriftLowerLimit",
          "RandomEffect_DriftZ", "RandomEffect_DriftProb",
          "H2DriftUpperLimit", "H2DriftLowerLimit")

        heterogeneity1DT[[i]] <- fixedEffectDriftResults1DT[[i]][8:15,]; heterogeneity1DT[[i]]
      } # end for (i in 1:length(dt))
      names(RandomEffectDriftResults1DT) <-  paste0("Time Interval = ", dt)
      names(heterogeneity1DT) <-  paste0("Time Interval = ", dt)
    } # end if (!(is.null(dt)))

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
    rownames(RandomEffectDriftResults2) <- c(#"MeanOfDriftValues",
      "RandomEffect_Drift", "RandomEffect_DriftVariance", "RandomEffect_DriftSE",
      "RandomEffect_DriftUpperLimit", "RandomEffect_DriftLowerLimit",
      "RandomEffect_DriftZ", "RandomEffect_DriftProb",
      "H2DriftUpperLimit", "H2DriftLowerLimit")

    ### some corrections for the output
    heterogeneity2 <- fixedEffectDriftResults2[8:15,]; heterogeneity2
  }

    # same for dt
    if (!(is.null(dt))) {
      RandomEffecttot_Drift2DT <- RandomEffecttot_DriftVariance2DT <- RandomEffecttot_DriftSE2DT <- list()
      RandomEffecttot_DriftUpperLimit2DT <- RandomEffecttot_DriftLowerLimit2DT <- RandomEffecttot_DriftProb12DT <- list()
      RandomEffecttot_DriftUpperLimitPI2DT <- RandomEffecttot_DriftLowerLimitPI2DT <- RandomEffectDriftResults2DT <- list()
      RandomEffecttot_DriftZ2DT <- RandomEffecttot_DriftProb2DT <- list()
      RandomEffecttot_DriftUpperLimitPI2DT <- RandomEffecttot_DriftLowerLimitPI2DT <- RandomEffectDriftResults2DT <- list()
      heterogeneity2DT <- list()
      for (i in 1:length(dt)) {
        RandomEffecttot_Drift2DT[[i]] <- Ttot_DriftMeans2DT[[i]]/Ttot_DriftWeights2DT[[i]]; RandomEffecttot_Drift2DT[[i]]
        RandomEffecttot_DriftVariance2DT[[i]] <- 1/Ttot_DriftWeights2DT[[i]]; RandomEffecttot_DriftVariance2DT[[i]]
        RandomEffecttot_DriftSE2DT[[i]] <- RandomEffecttot_DriftVariance2DT[[i]]^.5; RandomEffecttot_DriftSE2DT[[i]]
        RandomEffecttot_DriftUpperLimit2DT[[i]] <- RandomEffecttot_Drift2DT[[i]] + 1.96*RandomEffecttot_DriftSE2DT[[i]]; RandomEffecttot_DriftUpperLimit2DT[[i]]
        RandomEffecttot_DriftLowerLimit2DT[[i]] <- RandomEffecttot_Drift2DT[[i]] - 1.96*RandomEffecttot_DriftSE2DT[[i]]; RandomEffecttot_DriftLowerLimit2DT[[i]]
        RandomEffecttot_DriftZ2DT[[i]] <- RandomEffecttot_Drift2DT[[i]]/RandomEffecttot_DriftSE2DT[[i]]; RandomEffecttot_DriftZ2DT[[i]]
        RandomEffecttot_DriftProb2DT[[i]] <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ2DT[[i]]),
                                                                  mean=c(rep(0, (n.latent1^2))), sd=c(rep(1, (n.latent1^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb2DT[[i]]
        RandomEffecttot_DriftUpperLimitPI2DT[[i]] <- RandomEffecttot_Drift2DT[[i]] + 1.96*(tau2Drift2DT[[i]]^.5); RandomEffecttot_DriftUpperLimitPI2DT[[i]]
        RandomEffecttot_DriftLowerLimitPI2DT[[i]] <- RandomEffecttot_Drift2DT[[i]] - 1.96*(tau2Drift2DT[[i]]^.5); RandomEffecttot_DriftLowerLimitPI2DT[[i]]
        RandomEffectDriftResults2DT[[i]] <- rbind(RandomEffecttot_Drift2DT[[i]], RandomEffecttot_DriftVariance2DT[[i]], RandomEffecttot_DriftSE2DT[[i]],
                                                  RandomEffecttot_DriftUpperLimit2DT[[i]], RandomEffecttot_DriftLowerLimit2DT[[i]],
                                                  RandomEffecttot_DriftZ2DT[[i]], RandomEffecttot_DriftProb2DT[[i]],
                                                  RandomEffecttot_DriftUpperLimitPI2DT[[i]], RandomEffecttot_DriftLowerLimitPI2DT[[i]])
        rownames(RandomEffectDriftResults2DT[[i]]) <- c(#"MeanOfDriftValues",
          "RandomEffect_Drift", "RandomEffect_DriftVariance", "RandomEffect_DriftSE",
          "RandomEffect_DriftUpperLimit", "RandomEffect_DriftLowerLimit",
          "RandomEffect_DriftZ", "RandomEffect_DriftProb",
          "H2DriftUpperLimit", "H2DriftLowerLimit")

        heterogeneity2DT[[i]] <- fixedEffectDriftResults2DT[[i]][8:15,]; heterogeneity2DT[[i]]
      } # end for (i in 1:length(dt))
      names(heterogeneity2DT) <-  paste0("Time Interval = ", dt)
    } # end if (!(is.null(dt)))

  randomEffectDriftResults <- rbind(#cbind(t(MeanOfDriftValues1), t(MeanOfDriftValues2)),
    cbind(t(RandomEffecttot_Drift1), t(RandomEffecttot_Drift2)),
    cbind(t(RandomEffecttot_DriftVariance1), t(RandomEffecttot_DriftVariance2)), cbind(t(RandomEffecttot_DriftSE1), t(RandomEffecttot_DriftSE2)),
    cbind(t(RandomEffecttot_DriftUpperLimit1), t(RandomEffecttot_DriftUpperLimit2)), cbind(t(RandomEffecttot_DriftLowerLimit1), t(RandomEffecttot_DriftLowerLimit2)),
    cbind(t(RandomEffecttot_DriftZ1), t(RandomEffecttot_DriftZ2)), cbind(t(RandomEffecttot_DriftProb1), t(RandomEffecttot_DriftProb2)),
    cbind(t(tau2Drift1), t(tau2Drift2)), cbind(t(Q_Drift1), t(Q_Drift2)), cbind(t(H2_Drift1), t(H2_Drift2)),
    cbind(t(H2DriftUpperLimit1), t(H2DriftUpperLimit2)), cbind(t(H2DriftLowerLimit1), t(H2DriftLowerLimit2)),
    cbind(t(I2_Drift1), t(I2_Drift2)), cbind(t(I2DriftUpperLimit1), t(I2DriftUpperLimit2)), cbind(t(I2DriftLowerLimit1), t(I2DriftLowerLimit2)))
  rownames(randomEffectDriftResults) <- c(#"MeanOfDriftValues",
    "RandomEffecttot_Drift",
    "RandomEffect_DriftVariance", "RandomEffect_DriftSE", "RandomEffect_DriftUpperLimit",
    "RandomEffect_DriftLowerLimit", "RandomEffect_DriftZ", "RandomEffect_DriftProb",
    "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
    "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")
  #round(randomEffectDriftResults, 4)

  # same for dt
  if (!(is.null(dt))) {
    randomEffectDriftResultsDT <- list()
    for (i in 1:length(dt)) {
      randomEffectDriftResultsDT[[i]] <- rbind(#cbind(t(MeanOfDriftValues1), t(MeanOfDriftValues2)),
        cbind(t(RandomEffecttot_Drift1DT[[i]]), t(RandomEffecttot_Drift2DT[[i]])),
        cbind(t(RandomEffecttot_DriftVariance1DT[[i]]), t(RandomEffecttot_DriftVariance2DT[[i]])), cbind(t(RandomEffecttot_DriftSE1DT[[i]]), t(RandomEffecttot_DriftSE2DT[[i]])),
        cbind(t(RandomEffecttot_DriftUpperLimit1DT[[i]]), t(RandomEffecttot_DriftUpperLimit2DT[[i]])), cbind(t(RandomEffecttot_DriftLowerLimit1DT[[i]]), t(RandomEffecttot_DriftLowerLimit2DT[[i]])),
        cbind(t(RandomEffecttot_DriftZ1DT[[i]]), t(RandomEffecttot_DriftZ2DT[[i]])), cbind(t(RandomEffecttot_DriftProb1DT[[i]]), t(RandomEffecttot_DriftProb2DT[[i]])),
        cbind(t(tau2Drift1DT[[i]]), t(tau2Drift2DT[[i]])), cbind(t(Q_Drift1DT[[i]]), t(Q_Drift2DT[[i]])), cbind(t(H2_Drift1DT[[i]]), t(H2_Drift2DT[[i]])),
        cbind(t(H2DriftUpperLimit1DT[[i]]), t(H2DriftUpperLimit2DT[[i]])), cbind(t(H2DriftLowerLimit1DT[[i]]), t(H2DriftLowerLimit2DT[[i]])),
        cbind(t(I2_Drift1DT[[i]]), t(I2_Drift2DT[[i]])), cbind(t(I2DriftUpperLimit1DT[[i]]), t(I2DriftUpperLimit2DT[[i]])), cbind(t(I2DriftLowerLimit1DT[[i]]), t(I2DriftLowerLimit2DT[[i]])))
      rownames(randomEffectDriftResultsDT[[i]]) <- c(#"MeanOfDriftValues",
        "RandomEffecttot_Drift",
        "RandomEffect_DriftVariance", "RandomEffect_DriftSE", "RandomEffect_DriftUpperLimit",
        "RandomEffect_DriftLowerLimit", "RandomEffect_DriftZ", "RandomEffect_DriftProb",
        "tau2Drift", "Q_Drift", "H2_Drift", "H2DriftUpperLimit",
        "H2DriftLowerLimit", "I2_Drift", "I2DriftUpperLimit", "I2DriftLowerLimit")
    } # end for (i in 1:length(dt))
    names(randomEffectDriftResultsDT) <-  paste0("Time Interval = ", dt)
  } # end  if (!(is.null(dt)))

  # Collect Results for both fixed effect analysis ######
  modelResultsList1 <- list(DRIFT_ctmaFitObject = DRIFTCoeff1,
                            DRIFTSE1_ctmaFitObject = DRIFTSE1); modelResultsList1
  modelResultsList2 <- list(DRIFT_ctmaFitObjectMod = DRIFTCoeff2,
                            DRIFTSE1_ctmaFitObjectMod = DRIFTSE2); modelResultsList2
  # same for dt
  if (!(is.null(dt))) {
  modelResultsList1DT <- list(DRIFT_ctmaFitObject = DRIFTCoeff1DT,
                            DRIFTSE1_ctmaFitObject = DRIFTSE1DT); modelResultsList1DT
  modelResultsList2DT <- list(DRIFT_ctmaFitObjectMod = DRIFTCoeff2DT,
                            DRIFTSE1_ctmaFitObjectMod = DRIFTSE2DT); modelResultsList2DT
  }

  summaryList1 <- list(model="Analysis of Heterogeneity for ctmaFitObject, i.e., model without moderator effects.",
                       estimates1=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults1, digits),
                                       "Heterogeneity"=round(heterogeneity1, digits),
                                       "I2 message" = fixedEffectDriftMessage,
                                       "Tau2 message" = tau2DriftMessage,
                                       "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults1, digits))); summaryList1
  summaryList2 <- list(model="Analysis of Heterogeneity for ctmaFitObjectMod, i.e., model with moderator effects.",
                       estimates1=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults2, digits),
                                       "Heterogeneity"=round(heterogeneity2, digits),
                                       "I2 message" = fixedEffectDriftMessage,
                                       "Tau2 message" = tau2DriftMessage,
                                       "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults2, digits))); summaryList2
  # same for dt
  if (!(is.null(dt))) {
    summaryList1DT <- list(model="Analysis of Heterogeneity of Discrete Time Effecst for ctmaFitObject, i.e., model without moderator effects.",
                       estimates1=list("Fixed Effects of Drift Coefficients"=lapply(fixedEffectDriftResults1DT, function(x) round(x,digits)),
                                       "Heterogeneity"=lapply(heterogeneity1DT, function(x) round(x,digits)),
                                       "I2 message" = fixedEffectDriftMessage,
                                       "Tau2 message" = tau2DriftMessage,
                                       "Random Effects of Drift Coefficients"=lapply(RandomEffectDriftResults1DT, function(x) round(x, digits)))); summaryList1DT
  summaryList2DT <- list(model="Analysis of Heterogeneity of Discrete Time Effecst for ctmaFitObjectMod, i.e., model with moderator effects.",
                         estimates1=list("Fixed Effects of Drift Coefficients"=lapply(fixedEffectDriftResults2DT, function(x) round(x,digits)),
                                         "Heterogeneity"=lapply(heterogeneity2DT, function(x) round(x,digits)),
                                         "I2 message" = fixedEffectDriftMessage,
                                         "Tau2 message" = tau2DriftMessage,
                                         "Random Effects of Drift Coefficients"=lapply(RandomEffectDriftResults2DT, function(x) round(x, digits)))); summaryList2DT
  }
  #fixedEffectDriftResults
  FE <- round(fixedEffectDriftResults[1:8,], digits); FE
  RE <- round(randomEffectDriftResults[1:8,], digits); RE
  if (!(is.null(dt))) {
    FEDT <- lapply(fixedEffectDriftResultsDT, function(x) round(x[1:8,], digits)); FEDT
    REDT <- lapply(randomEffectDriftResultsDT, function(x) round(x[1:8,], digits)); REDT
  }

  #Het <- round(fixedEffectDriftResults[9:16,], digits); Het
  Het <- round(fixedEffectDriftResults[8:15,], digits); Het
  if (!(is.null(dt))) HetDT <- lapply(fixedEffectDriftResultsDT, function(x) round(x[8:15,], digits))

  names11 <- names(ctmaFitObjectMod$modelResults$DRIFT); names11

  fitObjDriftAndSE <- round(cbind(modelResultsList1[[1]], modelResultsList1[[2]]), digits); fitObjDriftAndSE
  colnames(fitObjDriftAndSE) <- paste0("w/o Mod ", c(names11, names11)); fitObjDriftAndSE
  colnames(fitObjDriftAndSE) <- gsub("DRIFT ", "", colnames(fitObjDriftAndSE)); fitObjDriftAndSE
  colnames(fitObjDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)] <- paste0(colnames(fitObjDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)], " SE"); fitObjDriftAndSE
  #fitObjDriftAndSE

  # same for dt
  if (!(is.null(dt))) {
    fitObjDriftAndSEDT <- list()
    for (i in 1:length(dt)) {
      #fitObjDriftAndSEDT[[i]] <- cbind(modelResultsList1DT[[i]][[1]], modelResultsList1DT[[i]][[2]]); fitObjDriftAndSEDT[[i]]
      fitObjDriftAndSEDT[[i]] <- cbind(modelResultsList1DT[[1]][[i]], modelResultsList1DT[[2]][[i]]); fitObjDriftAndSEDT[[i]]
      colnames(fitObjDriftAndSEDT[[i]]) <- paste0("w/o Mod ", c(names11, names11)); fitObjDriftAndSEDT[[i]]
      colnames(fitObjDriftAndSEDT[[i]]) <- gsub("DRIFT ", "", colnames(fitObjDriftAndSEDT[[i]])); fitObjDriftAndSEDT[[i]]
      colnames(fitObjDriftAndSEDT[[i]])[(n.latent1^2+1):(2*n.latent1^2)] <- paste0(colnames(fitObjDriftAndSEDT[[i]])[(n.latent1^2+1):(2*n.latent1^2)], " SE"); fitObjDriftAndSEDT[[i]]
    }
    names(fitObjDriftAndSEDT) <-  paste0("Time Interval = ", dt)
    #fitObjDriftAndSEDT
  }

  fitObjModDriftAndSE <- round(cbind(modelResultsList2[[1]], modelResultsList2[[2]]), digits); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE) <- paste0("with Mod ", c(names11, names11)); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)] <- paste0(colnames(fitObjModDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)], " SE"); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE) <- gsub("DRIFT ", "", colnames(fitObjModDriftAndSE)); fitObjModDriftAndSE
  #fitObjModDriftAndSE

  # same for dt
  if (!(is.null(dt))) {
    fitObjModDriftAndSEDT <- list()
    for (i in 1:length(dt)) {
      #i <- 1
      #fitObjModDriftAndSEDT[[i]] <- cbind(modelResultsList2DT[[i]][[1]], modelResultsList2DT[[i]][[2]]); fitObjModDriftAndSEDT[[i]]
      fitObjModDriftAndSEDT[[i]] <- cbind(modelResultsList2DT[[1]][[i]], modelResultsList2DT[[2]][[i]]); fitObjModDriftAndSEDT[[i]]
      colnames(fitObjModDriftAndSEDT[[i]]) <- paste0("with Mod ", c(names11, names11)); fitObjModDriftAndSEDT[[i]]
      colnames(fitObjModDriftAndSEDT[[i]]) <- gsub("DRIFT ", "", colnames(fitObjModDriftAndSEDT[[i]])); fitObjModDriftAndSEDT[[i]]
      colnames(fitObjModDriftAndSEDT[[i]])[(n.latent1^2+1):(2*n.latent1^2)] <- paste0(colnames(fitObjModDriftAndSEDT[[i]])[(n.latent1^2+1):(2*n.latent1^2)], " SE"); fitObjModDriftAndSEDT[[i]]
    }
    names(fitObjModDriftAndSEDT) <-  paste0("Time Interval = ", dt)
    #fitObjModDriftAndSEDT
  }


  # Analysis of Reduction in Heterogeneity by means of moderators #####
  heterogeneity1[heterogeneity1 < 0] <- .00001
  heterogeneity2[heterogeneity2 < 0] <- .00001
  #redHet <- round(heterogeneity1/heterogeneity2, digits)[c(2,3,6),]*100; redHet
  redHet <- round((heterogeneity1-heterogeneity2)/(heterogeneity1), digits)[c(2,3,6),]*100; redHet
  colnames(redHet) <- names11; redHet
  colnames(redHet) <- gsub("DRIFT ", "", colnames(redHet)); redHet
  rownames(redHet) <- gsub("_Drift2", "", rownames(redHet)); redHet
  rownames(redHet) <- paste0(rownames(redHet), ": % reduction (neg. values = % increase)"); redHet

  # same for dt
  if (!(is.null(dt))) {
    redHetDT <- list()
    for (i in 1:length(dt)) {
      heterogeneity1DT[[i]][heterogeneity1DT[[i]] < 0] <- .00001
      heterogeneity2DT[[i]][heterogeneity2DT[[i]] < 0] <- .00001
      redHetDT[[i]] <- round((heterogeneity1DT[[i]]-heterogeneity2DT[[i]])/(heterogeneity1DT[[i]]), digits)[c(2,3,6),]*100
      colnames(redHetDT[[i]]) <- names11; redHetDT[[i]]
      colnames(redHetDT[[i]]) <- gsub("DRIFT ", "", colnames(redHetDT[[i]])); redHetDT[[i]]
      rownames(redHetDT[[i]]) <- gsub("_Drift2", "", rownames(redHetDT[[i]])); redHetDT[[i]]
      rownames(redHetDT[[i]]) <- paste0(rownames(redHetDT[[i]]), ": % reduction (neg. values = % increase)"); redHetDT[[i]]
    }
    names(redHetDT) <-  paste0("Time Interval = ", dt)
    #redHetDT
  }

  if (is.null(dt)) {
    discreteTime <- "No dt effects estimated."
  } else {
    discreteTime <- list('fixedEffects (w/o Mod | with Mod)'=FEDT, 'randomEffects (w/o Mod | with Mod)'=REDT,
                         fitObj_DriftAndSE=fitObjDriftAndSEDT, fitObjMod_DriftAndSE=fitObjModDriftAndSEDT,
                         'heterogeneity (w/o Mod | with Mod)'=HetDT,
                         HeterogeneityReduction=redHetDT,
                         Note="Negative I2 values were set to .00001 for computation of reduction in heterogeneity.")
    }

  results <- list(activeDirectory=activeDirectory,
                  plot.type=NULL, model.type="BiG",
                  n.studies=n.studies2,
                  n.latent=n.latent2,
                  studyList=list(ctmaFitObject, ctmaFitObjectMod),
                  summary=list("continuousTime" = list('fixedEffects (w/o Mod | with Mod)'=FE, 'randomEffects (w/o Mod | with Mod)'=RE,
                                                       fitObj_DriftAndSE=fitObjDriftAndSE, fitObjMod_DriftAndSE=fitObjModDriftAndSE,
                                                       'heterogeneity (w/o Mod | with Mod)'=Het,
                                                       HeterogeneityReduction=redHet,
                                                       Note="Negative I2 values were set to .00001 for computation of reduction in heterogeneity."),
                               "discreteTime" = discreteTime))

  class(results) <- "CoTiMAFit"

  invisible(results)
} ### END function definition
