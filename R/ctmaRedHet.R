#' ctmaRedHet
#'
#' @description Computes the Reduction in Heterogeneity in drift effects after introducing study-level moderators
#'
#' @param activateRPB if TRUE, messages (warning, finished) could be send to smart phone (default = FALSE)
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param ctmaFitObject ctmaFit Object WITHOUT Moderators (obtained from \code{\link{ctmaFit}} with the arguments invariantDrift=\'none\' and scaleTI=FALSE)
#' @param ctmaFitObjectMod ctmaFit Object WITH Moderators (obtained from \code{\link{ctmaFit}} with the arguments invariantDrift=\'none\' and scaleTI=FALSE)
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

  # check if correct fit objects are specified ######################
  if (is.null(ctmaFitObject)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA (\"ctmaFit\") object has to be supplied to analyse something.
    \nYou could try setting \"ctmaFitObject=\"experimental\"\" to get estimates of heterogeneity reduction anyway.
    \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (is.null(ctmaFitObjectMod)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA (\"ctmaFit\") object with Moderators has to be supplied to analyse something. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  # Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################

  start.time <- Sys.time(); start.time

  ## extract model params ctmaFit #######################################################################
  {
    if (!(ctmaFitObject[[1]] == "experimental")) {
      driftNames1 <- ctmaFitObject$parameterNames$DRIFT; driftNames1
      n.latent1 <- ctmaFitObject$n.latent; n.latent1
      n.studies1 <- length(ctmaFitObject$studyList); n.studies1
      finishsamples1 = dim(ctmaFitObject$studyFitList$stanfit$rawposterior)[1]; finishsamples1


      if ((is.null(ctmaFitObject$argumentList$invariantDrift)) |
          (length(ctmaFitObject$argumentList$invariantDrift) == n.latent1^2) ){
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
        ErrorMsg <- "The ctmaFitObject provided war created with not a single drift parameter varying across primary study.
            Reduction in heterogeneity of drift effects cannot be computed because there is no heterogeneity."
        stop(ErrorMsg)
      }
    }

    if (!(ctmaFitObject[[1]] == "experimental")) {
      ### main effects estimated ##############################################################################
      tmp1 <- which(!(is.na(ctmaFitObject$studyFitList$ctstanmodelbase$pars$param))); tmp1
      pars1 <- ctmaFitObject$studyFitList$ctstanmodelbase$pars[tmp1, ]; pars1

      #### determine the positions where the study_dummy moderators of the drift effects are located ####
      targetRows1a <- which(pars1$matrix == "DRIFT"); targetRows1a
      #targetRows1b <- which(pars1$TI1_effect == "TRUE"); targetRows1b
      targetRows1b <- 1:length(pars1$TI1_effect); targetRows1b
      targetRows1 <- targetRows1a[targetRows1a %in% targetRows1b]; targetRows1
      #### positions in the rawdata for drift in the reference group and the moderator effects (latter == TIpredeffects) ######
      #mainEffectsPerStudy1 <- length(tmp1); mainEffectsPerStudy1
      #referenceTargetCols1 <- 1:(n.latent1^2); referenceTargetCols1
      referenceTargetCols1 <- targetRows1
      #### extract sampled effects #######
      drift1 <- ctmaFitObject$studyFitList$stanfit$rawposterior[,referenceTargetCols1]; drift1[1,] # will be tformed last
      TIpredEffectsList1 <- list()
      e1 <- ctsem::ctExtract(ctmaFitObject$studyFitList)
      for (s in 2:n.studies1) TIpredEffectsList1[[s-1]] <- e1$TIPREDEFFECT[,referenceTargetCols1, (s-1)]
    }

    ## extract model params ctmaFitMod #####################################################################
    driftNames2 <- ctmaFitObjectMod$parameterNames$DRIFT; driftNames2
    n.latent2 <- ctmaFitObjectMod$n.latent; n.latent2
    n.studies2 <- length(ctmaFitObjectMod$studyList); n.studies2
    finishsamples2 = dim(ctmaFitObjectMod$studyFitList$stanfit$rawposterior)[1]; finishsamples2

    if ((is.null(ctmaFitObjectMod$argumentList$invariantDrift)) |
        (length(ctmaFitObjectMod$argumentList$invariantDrift) == n.latent2^2) ){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      Msg <- "The ctmaFitObjectMod provided war created with not a single drift parameter varying across primary study.
            Reduction in heterogeneity of drift effects cannot be computed because there is no heterogeneity."
      message(Msg)
    }

    ### main effects estimated ####
    tmp2 <- which(!(is.na(ctmaFitObjectMod$studyFitList$ctstanmodelbase$pars$param))); tmp2
    pars2 <- ctmaFitObjectMod$studyFitList$ctstanmodelbase$pars[tmp2, ]; pars2
    #### determine the positions where the study_dummy moderators of the drift effects are located #####
    targetRows2a <- which(pars2$matrix == "DRIFT"); targetRows2a
    #targetRows2b <- which(pars2$TI1_effect == "TRUE"); targetRows2b
    targetRows2b <- 1:length(pars2$TI1_effect); targetRows2b
    targetRows2 <- targetRows2a[targetRows2a %in% targetRows2b]; targetRows2
    referenceTargetCols2 <- targetRows2
    #### positions in the rawdata for drift in the reference group and the moderator effects (latter == TIpredeffects) #####
    #mainEffectsPerStudy2 <- length(tmp2); mainEffectsPerStudy2
    #referenceTargetCols2 <- 1:(n.latent2^2); referenceTargetCols2
    #### extract sampled effects #######
    drift2 <- ctmaFitObjectMod$studyFitList$stanfit$rawposterior[,referenceTargetCols2]; drift2[1,] # will be tformed

    TIpredEffectsList2 <- list()
    e2 <- ctsem::ctExtract(ctmaFitObjectMod$studyFitList)
    for (s in 2:n.studies2) TIpredEffectsList2[[s-1]] <- e2$TIPREDEFFECT[,referenceTargetCols2, (s-1)]
    #apply(TIpredEffectsList2[[1]], 2, mean)
    #apply(drift2, 2, mean)

    tipred.values2 <- list()
    for (i in 1:(n.studies2-1)) {
      tipred.values2[[i]] <- matrix(ctmaFitObjectMod$studyFitList$data$tipredsdata[,i])
      tipred.values2[[i]] <- mean(tipred.values2[[i]][tipred.values2[[i]] != 0])
    }
    #tipred.values2 # should be 1 for all in case TIpred was not scaled

    if (ctmaFitObject[[1]] == "experimental") {
      modEffectsList2 <- list()
      # how many additional moderator effects?
      n.mod.effects <- length(grep("_effect", colnames(pars2))) - n.studies2 + 1; n.mod.effects

      # raw mod effects
      for (s in n.studies2:(n.studies2+n.mod.effects-1)) {
        modEffectsList2[[s-n.studies2+1]] <- e2$TIPREDEFFECT[,referenceTargetCols2, s]
      }

      # transformed mod effects
      #n.mod.params <- dim(e2$TIPREDEFFECT)[2]; n.mod.params
      #startCol <- (n.studies2-1) * n.mod.params + 1; startCol
      #for (s in 1:n.mod.effects) {
      #  endCol <- startCol + length(targetRows2a) - 1; endCol
      #  modEffectsList2[[s]] <- ctmaFitObjectMod$studyFitList$stanfit$transformedpars$tipredeffectparams[, startCol:endCol]
      #  if (s > 1) startCol <- startCol + n.mod.params
      #}

      # test
      #modEffectsList2[[1]] <- e2$linearTIPREDEFFECT[, 1:4, 8]
      #apply(modEffectsList2[[1]], 2, mean)

      # cont
      mod.values2 <-list()
      for (i in 1:n.studies2) {
        mod.values2[[i]] <- ctmaFitObjectMod$studyList[[i]]$moderators[ctmaFitObjectMod$argumentList$mod.number]
      }

      # cat
      if (ctmaFitObjectMod$argumentList$mod.type == 'cat') {
        tmp1 <- grep("TI", colnames(ctmaFitObjectMod$data)); tmp1
        mod.pos <- tmp1[n.studies2:length(tmp1)]; mod.pos
        allSampleSizes <- ctmaFitObjectMod$statisticsList$allSampleSizes; allSampleSizes
        cumSampleSizes <- cumsum(allSampleSizes) - allSampleSizes[1] + 1; cumSampleSizes
        cumSampleSizes <- cumSampleSizes * ctmaFitObjectMod$statisticsList$maxTpoints - ctmaFitObjectMod$statisticsList$maxTpoints + 1; cumSampleSizes
        mod.values2 <-list()
        for (i in 1:n.studies2) {
          #i <- 1
          mod.values2[[i]] <- NA
          for (k in 1:n.mod.effects) {
            #k <- 2
            mod.values2[[i]][k] <- ctmaFitObjectMod$data[cumSampleSizes[i], mod.pos[k]]
          }
        }
        #mod.values2
      }


      if (ctmaFitObjectMod$argumentList$scaleMod == TRUE) {
        for (i in 1:length(mod.values2[[1]])) {
          length(mod.values2[[1]])
          #i <- 1
          tmp <- scale(unlist(lapply(mod.values2, function(x) x[i]))); tmp
          for (j in 1:length(tmp)) {
            #j <- 1
            mod.values2[[j]][i] <- tmp[j]
          }
        }
      }
    }
  }

  # do some checks here that fit are compatible
  # tbd
  # finishsamples may differ but a warning shall be issued


  print(paste0("#################################################################################"))
  print(paste0("## Extracting raw estimates for main and moderating effects: Will take a while! #"))
  print(paste0("#################################################################################"))

  # get effects from ctmaFit without moderators #########################################################
  if (!(ctmaFitObject[[1]] == "experimental")) {
    ## get tforms ######
    tmp1a <- ctmaFitObject$studyFitList$ctstanmodelbase$pars[, "transform"]; tmp1a
    tmp1b <- which(ctmaFitObject$studyFitList$ctstanmodelbase$pars$matrix == "DRIFT"); tmp1b
    transforms1 <- tmp1a[tmp1b]; transforms1

    ## apply tforms to reference group in ctmaFit #####
    modDriftMean1 <- modDriftSD1 <- list()
    modDriftMean1[[n.studies1]] <- matrix(NA, n.latent1, n.latent1)
    modDriftSD1[[n.studies1]] <- matrix(NA, n.latent1, n.latent1)
    counter <- 0
    for (l in 1:(n.latent1)) {
      for (m in 1:(n.latent1)) {
        counter <- counter + 1; counter
        transformedParam <- c()
        for (f in 1:finishsamples1) {
          param <- matrix(drift1[f, ], n.latent1, n.latent1, byrow=T)
          param <- param[l,m]; param
          transformedParam <- c(transformedParam, eval(parse(text=transforms1[counter]))); transformedParam
        }
        #drift1
        modDriftMean1[[n.studies1]][l,m] <- mean(transformedParam)
        modDriftSD1[[n.studies1]][l,m] <- sd(transformedParam)
        # undoTimScaling after tform
        if (undoTimeScaling & (!(is.null(ctmaFitObject$summary$scaledTime)) ) ) {
          modDriftMean1[[n.studies1]][l,m] <- modDriftMean1[[n.studies1]][l,m] * ctmaFitObject$summary$scaledTime
          modDriftSD1[[n.studies1]][l,m] <- modDriftSD1[[n.studies1]][l,m] * ctmaFitObject$summary$scaledTime
        }
      }
    }
    ## apply tforms to other groups
    for (k in 1:(n.studies1-1)) {
      counter <- 0
      modDriftMean1[[k]] <- matrix(NA, n.latent1, n.latent1)
      modDriftSD1[[k]] <- matrix(NA, n.latent1, n.latent1)
      for (l in 1:(n.latent1)) {
        for (m in 1:(n.latent1)) {
          counter <- counter + 1; counter
          transformedParam <- c()
          for (f in 1:finishsamples1) {
            #f <- 1
            tmp1 <- matrix(drift1[f, ], n.latent1, n.latent1, byrow=T); tmp1
            tmp2 <- matrix(TIpredEffectsList1[[k]][f, ], n.latent1, n.latent1, byrow=T); tmp2
            param <- tmp1[l,m] +tmp2[l,m]; param
            transformedParam <- c(transformedParam, eval(parse(text=transforms1[counter]))); transformedParam
          }
          modDriftMean1[[k]][l,m] <- mean(transformedParam)
          modDriftSD1[[k]][l,m] <- sd(transformedParam)
          # undoTimScaling after tform
          if (undoTimeScaling & (!(is.null(ctmaFitObject$summary$scaledTime)) ) ) {
            modDriftMean1[[k]][l,m] <- modDriftMean1[[k]][l,m] * ctmaFitObject$summary$scaledTime
            modDriftSD1[[k]][l,m] <- modDriftSD1[[k]][l,m] * ctmaFitObject$summary$scaledTime
          }
        }
      }
    }
    allDriftsMean1 <- matrix(unlist(modDriftMean1), ncol=n.latent1^2, byrow=T); allDriftsMean1
    allDriftsSD1 <- matrix(unlist(modDriftSD1), ncol=n.latent1^2, byrow=T); allDriftsSD1
  }

  # get effects from ctmaFitMod with moderators #########################################################
  {
    ## get tforms ####
    tmp1a <- ctmaFitObjectMod$studyFitList$ctstanmodelbase$pars[, "transform"]; tmp1a
    tmp1b <- which(ctmaFitObjectMod$studyFitList$ctstanmodelbase$pars$matrix == "DRIFT"); tmp1b
    transforms2 <- tmp1a[tmp1b]; transforms2

    ### apply tforms to reference group in ctmaFit ######
    modDriftMean2 <- modDriftSD2 <- list()
    modDriftMean2[[n.studies2]] <- matrix(NA, n.latent2, n.latent2)
    modDriftSD2[[n.studies2]] <- matrix(NA, n.latent2, n.latent2)
    if (ctmaFitObject[[1]] == "experimental") {
      modDriftMean3 <- modDriftSD3 <- list()
      modDriftMean3[[n.studies2]] <- matrix(NA, n.latent2, n.latent2)
      modDriftSD3[[n.studies2]] <- matrix(NA, n.latent2, n.latent2)
    }
    counter <- 0
    for (l in 1:(n.latent2)) {
      #l <- 1
      for (m in 1:(n.latent2)) {
        #m <- 2
        counter <- counter + 1; counter
        transformedParam <- c()
        if (ctmaFitObject[[1]] == "experimental") transformedParam2 <- c()
        for (f in 1:finishsamples2) {
          #f <- 1
          param <- matrix(drift2[f, ], n.latent2, n.latent2, byrow=T)
          param <- param[l,m]; param
          transformedParam <- c(transformedParam, eval(parse(text=transforms2[counter]))); transformedParam
          if (ctmaFitObject[[1]] == "experimental") {
            for (n in 1:n.mod.effects) {
              #n <- 1
              tmp3 <- matrix(modEffectsList2[[n]][f, ], n.latent2, n.latent2, byrow=T); tmp3
              param <- param + tmp3[l,m] * mod.values2[[n.studies2]][n] # reference study is last
              #param
            }
            transformedParam2 <- c(transformedParam2, eval(parse(text=transforms2[counter]))); transformedParam2
          }
        }
        modDriftMean2[[n.studies2]][l,m] <- mean(transformedParam)
        modDriftSD2[[n.studies2]][l,m] <- sd(transformedParam)
        modDriftMean3[[n.studies2]][l,m] <- mean(transformedParam2)
        modDriftSD3[[n.studies2]][l,m] <- sd(transformedParam2)

        # undoTimScaling after tform
        if (!(ctmaFitObject[[1]] == "experimental")) {
          if (undoTimeScaling & (!(is.null(ctmaFitObject$summary$scaledTime)) ) ) {
            modDriftMean2[[n.studies2]][l,m] <- modDriftMean2[[n.studies2]][l,m] * ctmaFitObjectMod$summary$scaledTime
            modDriftSD2[[n.studies2]][l,m] <- modDriftSD2[[n.studies2]][l,m] * ctmaFitObjectMod$summary$scaledTime
          }
        }
        if (ctmaFitObject[[1]] == "experimental") {
          if (undoTimeScaling & (!(is.null(ctmaFitObjectMod$summary$scaledTime)) ) ) {
            modDriftMean3[[n.studies2]][l,m] <- modDriftMean3[[n.studies2]][l,m] * ctmaFitObjectMod$summary$scaledTime
            modDriftSD3[[n.studies2]][l,m] <- modDriftSD3[[n.studies2]][l,m] * ctmaFitObjectMod$summary$scaledTime
          }
        }
      }
    }
    modDriftMean2 # this is the reported overall single drift matrix from ctmaFit with invariant=none + 1 moderator (i.e., reference study incl. its moderator effect)
    modDriftMean3 # this is the so-far unknown effect for last study with mod-effect partialled
    # there are deviations from ctmaInit because moderator cannot explain all variance
    modDriftSD2 # this is not the (reduced) variance!!
    modDriftSD3

    ## apply tforms to other groups ######
    for (k in 1:(n.studies2-1)) {
      #k <- 1
      counter <- 0
      modDriftMean2[[k]] <- matrix(NA, n.latent2, n.latent2)
      modDriftSD2[[k]] <- matrix(NA, n.latent2, n.latent2)
      if (ctmaFitObject[[1]] == "experimental") {
        modDriftMean3[[k]] <- matrix(NA, n.latent2, n.latent2)
        modDriftSD3[[k]] <- matrix(NA, n.latent2, n.latent2)
      }
      for (l in 1:(n.latent2)) {
        #l <- 1
        for (m in 1:(n.latent2)) {
          #m <- 1
          counter <- counter + 1; counter
          transformedParam <- c()
          if (ctmaFitObject[[1]] == "experimental") transformedParam2 <- c()
          for (f in 1:finishsamples2) {
            #f <- 1
            tmp1 <- matrix(drift2[f, ], n.latent2, n.latent2, byrow=T); tmp1
            tmp2 <- matrix(TIpredEffectsList2[[k]][f, ], n.latent2, n.latent2, byrow=T); tmp2
            param <- tmp1[l,m] + tmp2[l,m] * tipred.values2[[k]]; param # effects wo substantive moderator (only study dummies)
            transformedParam <- c(transformedParam, eval(parse(text=transforms2[counter]))); transformedParam
            if (ctmaFitObject[[1]] == "experimental") {
              for (n in 1:n.mod.effects) {
                tmp3 <- matrix(modEffectsList2[[n]][f, ], n.latent2, n.latent2, byrow=T); tmp3
                param <- param + tmp3[l,m] * mod.values2[[k]][n]
              }
              transformedParam2 <- c(transformedParam2, eval(parse(text=transforms2[counter]))); transformedParam2
            }
          }

          modDriftMean2[[k]][l,m] <- mean(transformedParam)
          modDriftSD2[[k]][l,m] <- sd(transformedParam)
          if (ctmaFitObject[[1]] == "experimental") {
            modDriftMean3[[k]][l,m] <- mean(transformedParam2)
            modDriftSD3[[k]][l,m] <- sd(transformedParam2)
          }
          # undoTimScaling after tform
          if (!(ctmaFitObject[[1]] == "experimental")) {
            if (undoTimeScaling & (!(is.null(ctmaFitObject$summary$scaledTime)) ) ) {
              modDriftMean2[[k]][l,m] <- modDriftMean2[[k]][l,m] * ctmaFitObjectMod$summary$scaledTime
              modDriftSD2[[k]][l,m] <- modDriftSD2[[k]][l,m] * ctmaFitObjectMod$summary$scaledTime
            }
          }
          if (ctmaFitObject[[1]] == "experimental") {
            if (undoTimeScaling & (!(is.null(ctmaFitObjectMod$summary$scaledTime)) ) ) {
              modDriftMean3[[k]][l,m] <- modDriftMean3[[k]][l,m] * ctmaFitObjectMod$summary$scaledTime
              modDriftSD3[[k]][l,m] <- modDriftSD3[[k]][l,m] * ctmaFitObjectMod$summary$scaledTime
            }
          }
        }
      }
    }
    # here was a }
    allDriftsMean2<- matrix(unlist(modDriftMean2), ncol=n.latent2^2, byrow=T); allDriftsMean2 # ~ ctmaInit result, but mod effect is included
    allDriftsSD2 <- matrix(unlist(modDriftSD2), ncol=n.latent2^2, byrow=T); allDriftsSD2
    if (ctmaFitObject[[1]] == "experimental") {
      allDriftsMean3<- matrix(unlist(modDriftMean3), ncol=n.latent2^2, byrow=T); allDriftsMean3
      allDriftsSD3 <- matrix(unlist(modDriftSD3), ncol=n.latent2^2, byrow=T); allDriftsSD3
    }
  }
  #allDriftsMean2
  #allDriftsMean3

  # 'translate' objects into objects as used in ctmaBiG #####
  {
    names11 <- names(ctmaFitObjectMod$modelResults$DRIFT); names11
    names12 <- names(ctmaFitObjectMod$modelResults$DIFFUSION); names12
    names13 <- names(ctmaFitObjectMod$modelResults$T0VAR); names13

    if (!(ctmaFitObject[[1]] == "experimental")) {
      all_Coeff1 <- allDriftsMean1; all_Coeff1
      all_SE1 <- allDriftsSD1; all_SE1
      colnames(all_SE1) <- colnames(all_Coeff1) <-c(names11)
    }
    #
    all_Coeff2 <- allDriftsMean2; all_Coeff2
    all_SE2 <- allDriftsSD2; all_SE2
    colnames(all_SE2) <- colnames(all_Coeff2) <- c(names11)
    #
    if (ctmaFitObject[[1]] == "experimental") {
      all_Coeff3 <- allDriftsMean3; all_Coeff3
      all_SE3 <- allDriftsSD3; all_SE3
      colnames(all_SE3) <- colnames(all_Coeff3) <- c(names11)
    }

    # check if all SE are < 0.0 #####
    if (!(ctmaFitObject[[1]] == "experimental")) {
      if (any(all_SE1 == 0)) {
        tmp <- which(all_SE1 == 0, arr.ind = TRUE); tmp
        colnames(tmp) <- c("Study", "Drift Effect"); tmp
        ErrorMsg <- paste0("\nAt least one SE of drift effects in ctmaFitObject was zero. Analysis of heterogeneity and bias cannot be performed.
                         Possible problem in fitting the single studies listed above. \nGood luck for the next try!")
        print(tmp)
        stop(ErrorMsg)
      }
    }
    if (any(all_SE2 == 0)) {
      tmp <- which(all_SE2 == 0, arr.ind = TRUE); tmp
      colnames(tmp) <- c("Study", "Drift Effect"); tmp
      ErrorMsg <- paste0("\nAt least one SE of drift effects in ctmaFitObjectMod was zero. Analysis of heterogeneity and bias cannot be performed.
                         Possible problem in fitting the single studies listed above. \nGood luck for the next try!")
      print(tmp)
      stop(ErrorMsg)
    }
    if (ctmaFitObject[[1]] == "experimental") {
      if (any(all_SE3 == 0)) {
        tmp <- which(all_SE3 == 0, arr.ind = TRUE); tmp
        colnames(tmp) <- c("Study", "Drift Effect"); tmp
        ErrorMsg <- paste0("\nAt least one SE of the primary study (sic) drift effects in ctmaFitModObject was zero. Analysis of heterogeneity and bias cannot be performed.
                         Possible problem in fitting the single studies listed above. \nGood luck for the next try!")
        print(tmp)
        stop(ErrorMsg)
      }
    }


    # dt effects
    # tbd
  }

  # Fixed & Random Effects Analyses #####################################################################

  #DRIFTCoeff <- DIFFUSIONCoeff <- T0VARCoeff <- DRIFTSE <- DIFFUSIONSE <- T0VARSE <- matrix(NA, n.studies1, n.latent1^2)

  if (!(ctmaFitObject[[1]] == "experimental")) {
    DRIFTCoeff1 <- all_Coeff1; DRIFTCoeff1
    DRIFTSE1 <- all_SE1; DRIFTSE1
  }
  {
    DRIFTCoeff2 <- all_Coeff2; DRIFTCoeff2
    DRIFTSE2 <- all_SE2; DRIFTSE2
  }
  if (ctmaFitObject[[1]] == "experimental") {
    DRIFTCoeff3 <- all_Coeff3; DRIFTCoeff3
    DRIFTSE3 <- all_SE3; DRIFTSE3
  }

  if (!(ctmaFitObject[[1]] == "experimental")) {
    DRIFTCoeffSND1 <- DRIFTCoeff1 / DRIFTSE1; DRIFTCoeffSND1
    DRIFTPrecision1 <- c(rep(1, n.latent1^2))/(DRIFTSE1); DRIFTPrecision1
    colnames(DRIFTPrecision1) <- colnames(DRIFTCoeffSND1); DRIFTPrecision1
    if (!(is.null(dt))) {
      # tbd
    }
  }
  {
    DRIFTCoeffSND2 <- DRIFTCoeff2 / DRIFTSE2; DRIFTCoeffSND2
    DRIFTPrecision2 <- c(rep(1, n.latent2^2))/(DRIFTSE2); DRIFTPrecision2
    colnames(DRIFTPrecision2) <- colnames(DRIFTCoeffSND2); DRIFTPrecision2
    if (!(is.null(dt))) {
      # tbd
    }
  }
  if (ctmaFitObject[[1]] == "experimental") {
    DRIFTCoeffSND3 <- DRIFTCoeff3 / DRIFTSE3; DRIFTCoeffSND3
    DRIFTPrecision3 <- c(rep(1, n.latent2^2))/(DRIFTSE3); DRIFTPrecision3
    colnames(DRIFTPrecision3) <- colnames(DRIFTCoeffSND3); DRIFTPrecision3
    if (!(is.null(dt))) {
      # tbd
    }
  }

  ## FIXED EFFECTS ANALYSIS ###############################################################################
  {
    if (!(ctmaFitObject[[1]] == "experimental")) {
      DriftMeans1 <- colMeans(DRIFTCoeff1); DriftMeans1
      T_DriftWeights1 <- colSums(DRIFTPrecision1^2); T_DriftWeights1
      T_DriftMeans1 <- colSums(DRIFTCoeff1 * DRIFTPrecision1^2); T_DriftMeans1
      names(T_DriftMeans1) <- names(T_DriftWeights1); T_DriftMeans1
    }
    DriftMeans2 <- colMeans(DRIFTCoeff2); DriftMeans2
    if (!(is.null(dt))) {} # tbd DriftMeans_dt <- colMeans(drift_Coeff_dt)
    # Sum of within weights  and weight * effect size
    T_DriftWeights2 <- colSums(DRIFTPrecision2^2); T_DriftWeights2
    if (!(is.null(dt))) {} # tbd T_DriftWeights_dt <- colSums(DRIFTPrecision_dt^2)
    # DRIFTPrecision
    T_DriftMeans2 <- colSums(DRIFTCoeff2 * DRIFTPrecision2^2); T_DriftMeans2
    names(T_DriftMeans2) <- names(T_DriftWeights2); T_DriftMeans2
    if (!(is.null(dt))) {} # tbd
    #T_DriftMeans_dt <- colSums(drift_Coeff_dt * DRIFTPrecision_dt^2); T_DriftMeans_dt
    #names(T_DriftMeans_dt) <- names(T_DriftWeights); T_DriftMeans_dt
    if (ctmaFitObject[[1]] == "experimental") {
      DriftMeans3 <- colMeans(DRIFTCoeff3); DriftMeans3
      T_DriftWeights3 <- colSums(DRIFTPrecision3^2); T_DriftWeights3
      T_DriftMeans3 <- colSums(DRIFTCoeff3 * DRIFTPrecision3^2); T_DriftMeans3
      names(T_DriftMeans3) <- names(T_DriftWeights3); T_DriftMeans3
    }
    # here was a }

    ### Fixed effects results for ctmaFitObject without moderators ######
    if (!(ctmaFitObject[[1]] == "experimental")) {
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
      # same for dt
      if (!(is.null(dt))) {
        #FixedEffect_Drift_dt <- T_DriftMeans_dt/T_DriftWeights_dt; FixedEffect_Drift_dt
        #FixedEffect_DriftVariance_dt <- 1/T_DriftWeights_dt; FixedEffect_DriftVariance_dt
        #FixedEffect_DriftSE_dt <- FixedEffect_DriftVariance_dt^.5; FixedEffect_DriftSE_dt
        #FixedEffect_DriftUpperLimit_dt <- FixedEffect_Drift_dt + 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftUpperLimit_dt
        #FixedEffect_DriftLowerLimit_dt <- FixedEffect_Drift_dt - 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftLowerLimit_dt
        #FixedEffect_DriftZ_dt <- FixedEffect_Drift_dt/FixedEffect_DriftSE_dt; FixedEffect_DriftZ_dt
        #FixedEffect_DriftProb_dt <- round(1-stats::pnorm(abs(FixedEffect_DriftZ_dt),
        #                                                 mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); FixedEffect_DriftProb_dt
        #Q_Drift_dt <- colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt^2)- (colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt))^2 / colSums(DRIFTPrecision_dt^2); Q_Drift_dt
        #H2_Drift_dt <- Q_Drift_dt/(n.studies-1); H2_Drift_dt
        #I2_Drift_dt <- (H2_Drift_dt-1)/H2_Drift_dt*100; I2_Drift_dt
      }
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
      if (!(is.null(dt))) { # tbd
        #T2_DriftWeights_dt <- colSums(DRIFTPrecision_dt^2^2); T2_DriftWeights_dt # Borenstein et al., 2007, p. 98
        #cDrift_dt <- T_DriftWeights_dt-T2_DriftWeights_dt/T_DriftWeights_dt; cDrift_dt
        #tau2Drift_dt <- (Q_Drift_dt-(n.studies-1))/cDrift_dt; tau2Drift_dt
        #SElnHDrift_dt <- c()
        #SElnHDrift_dt[] <- 0
        #for (j in 1:(n.latent^2)) {
        #  if (Q_Drift_dt[j] > n.studies) SElnHDrift_dt[j] <- 1/2*(log(Q_Drift_dt[j])-log(n.studies-1))/((2*Q_Drift_dt[j])^.5-(2*(n.studies-1)-1)^.5)
        #  if (Q_Drift_dt[j] <= n.studies) SElnHDrift_dt[j] <-  (1/(2*(n.studies-2)) * (1 - 1/(3*(n.studies-2)^.5)) )^.5
        #}
        #H2DriftUpperLimit_dt <- exp(log(H2_Drift_dt) + 1.96*SElnHDrift_dt); H2DriftUpperLimit_dt
        #H2DriftLowerLimit_dt <- exp(log(H2_Drift_dt) - 1.96*SElnHDrift_dt); H2DriftLowerLimit_dt
        #L_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))-1.96*SElnHDrift_dt)
        #U_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))+1.96*SElnHDrift_dt)
        #I2DriftUpperLimit_dt <- (U_dt^2-1)/U_dt^2 * 100; I2DriftUpperLimit_dt
        #I2DriftLowerLimit_dt <- (L_dt^2-1)/L_dt^2 * 100; I2DriftLowerLimit_dt
      }
    }


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
      # same for dt
      if (!(is.null(dt))) {
        #FixedEffect_Drift_dt <- T_DriftMeans_dt/T_DriftWeights_dt; FixedEffect_Drift_dt
        #FixedEffect_DriftVariance_dt <- 1/T_DriftWeights_dt; FixedEffect_DriftVariance_dt
        #FixedEffect_DriftSE_dt <- FixedEffect_DriftVariance_dt^.5; FixedEffect_DriftSE_dt
        #FixedEffect_DriftUpperLimit_dt <- FixedEffect_Drift_dt + 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftUpperLimit_dt
        #FixedEffect_DriftLowerLimit_dt <- FixedEffect_Drift_dt - 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftLowerLimit_dt
        #FixedEffect_DriftZ_dt <- FixedEffect_Drift_dt/FixedEffect_DriftSE_dt; FixedEffect_DriftZ_dt
        #FixedEffect_DriftProb_dt <- round(1-stats::pnorm(abs(FixedEffect_DriftZ_dt),
        #                                                 mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); FixedEffect_DriftProb_dt
        #Q_Drift_dt <- colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt^2)- (colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt))^2 / colSums(DRIFTPrecision_dt^2); Q_Drift_dt
        #H2_Drift_dt <- Q_Drift_dt/(n.studies-1); H2_Drift_dt
        #I2_Drift_dt <- (H2_Drift_dt-1)/H2_Drift_dt*100; I2_Drift_dt
      }
      # Tau square
      T2_DriftWeights2 <- colSums(DRIFTPrecision2^2^2); T2_DriftWeights2 # Borenstein et al., 2007, p. 98
      cDrift2 <- T_DriftWeights2 - T2_DriftWeights2/T_DriftWeights2; cDrift2
      tau2Drift2 <- (Q_Drift2 - (n.studies2-1))/cDrift2; tau2Drift2
      SElnHDrift2 <- c()
      SElnHDrift2[] <- 0
      for (j in 1:(n.latent2^2)) {
        if (Q_Drift2[j] > n.studies2) SElnHDrift2[j] <- 1/2*(log(Q_Drift2[j])-log(n.studies2-1))/((2*Q_Drift2[j])^.5-(2*(n.studies2-1)-1)^.5)
        if (Q_Drift2[j] <= n.studies2) SElnHDrift2[j] <-  (1/(2*(n.studies2-2)) * (1 - 1/(3*(n.studies2-2)^.5)) )^.5
      }
      H2DriftUpperLimit2 <- exp(log(H2_Drift2) + 1.96*SElnHDrift2); H2DriftUpperLimit2
      H2DriftLowerLimit2 <- exp(log(H2_Drift2) - 1.96*SElnHDrift2); H2DriftLowerLimit2
      L2 <- exp(0.5*log(Q_Drift2/(n.studies2-1))-1.96*SElnHDrift2)
      U2 <- exp(0.5*log(Q_Drift2/(n.studies2-1))+1.96*SElnHDrift2)
      I2DriftUpperLimit2 <- (U2^2-1)/U2^2 * 100; I2DriftUpperLimit2
      I2DriftLowerLimit2 <- (L2^2-1)/L2^2 * 100; I2DriftLowerLimit2
      # same for dt
      if (!(is.null(dt))) { # tbd
        #T2_DriftWeights_dt <- colSums(DRIFTPrecision_dt^2^2); T2_DriftWeights_dt # Borenstein et al., 2007, p. 98
        #cDrift_dt <- T_DriftWeights_dt-T2_DriftWeights_dt/T_DriftWeights_dt; cDrift_dt
        #tau2Drift_dt <- (Q_Drift_dt-(n.studies-1))/cDrift_dt; tau2Drift_dt
        #SElnHDrift_dt <- c()
        #SElnHDrift_dt[] <- 0
        #for (j in 1:(n.latent^2)) {
        #  if (Q_Drift_dt[j] > n.studies) SElnHDrift_dt[j] <- 1/2*(log(Q_Drift_dt[j])-log(n.studies-1))/((2*Q_Drift_dt[j])^.5-(2*(n.studies-1)-1)^.5)
        #  if (Q_Drift_dt[j] <= n.studies) SElnHDrift_dt[j] <-  (1/(2*(n.studies-2)) * (1 - 1/(3*(n.studies-2)^.5)) )^.5
        #}
        #H2DriftUpperLimit_dt <- exp(log(H2_Drift_dt) + 1.96*SElnHDrift_dt); H2DriftUpperLimit_dt
        #H2DriftLowerLimit_dt <- exp(log(H2_Drift_dt) - 1.96*SElnHDrift_dt); H2DriftLowerLimit_dt
        #L_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))-1.96*SElnHDrift_dt)
        #U_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))+1.96*SElnHDrift_dt)
        #I2DriftUpperLimit_dt <- (U_dt^2-1)/U_dt^2 * 100; I2DriftUpperLimit_dt
        #I2DriftLowerLimit_dt <- (L_dt^2-1)/L_dt^2 * 100; I2DriftLowerLimit_dt
      }
    }


    ### EXPERIMENTAL Fixed effects results for ctmaFitObjectMod with moderators ######
    if (ctmaFitObject[[1]] == "experimental") {
      FixedEffect_Drift3 <- T_DriftMeans3/T_DriftWeights3; FixedEffect_Drift3
      FixedEffect_DriftVariance3 <- 1/T_DriftWeights3; FixedEffect_DriftVariance3
      FixedEffect_DriftSE3 <- FixedEffect_DriftVariance3^.5; FixedEffect_DriftSE3
      FixedEffect_DriftUpperLimit3 <- FixedEffect_Drift3 + 1.96*FixedEffect_DriftSE3; FixedEffect_DriftUpperLimit3
      FixedEffect_DriftLowerLimit3 <- FixedEffect_Drift3 - 1.96*FixedEffect_DriftSE3; FixedEffect_DriftLowerLimit3
      FixedEffect_DriftZ3 <- FixedEffect_Drift3/FixedEffect_DriftSE3; FixedEffect_DriftZ3
      FixedEffect_DriftProb3 <- round(1-stats::pnorm(abs(FixedEffect_DriftZ3),
                                                     mean=c(rep(0, (n.latent2^2))), sd=c(rep(1, (n.latent2^2))), log.p=F), digits=digits); FixedEffect_DriftProb3
      Q_Drift3 <- colSums(DRIFTPrecision3^2 * DRIFTCoeff3^2)- (colSums(DRIFTPrecision3^2 * DRIFTCoeff3))^2 / colSums(DRIFTPrecision3^2); Q_Drift3
      H2_Drift3 <- Q_Drift3/(n.studies2-1); H2_Drift3
      I2_Drift3 <- (H2_Drift3-1)/H2_Drift3*100; I2_Drift3
      # same for dt
      if (!(is.null(dt))) {
        #FixedEffect_Drift_dt <- T_DriftMeans_dt/T_DriftWeights_dt; FixedEffect_Drift_dt
        #FixedEffect_DriftVariance_dt <- 1/T_DriftWeights_dt; FixedEffect_DriftVariance_dt
        #FixedEffect_DriftSE_dt <- FixedEffect_DriftVariance_dt^.5; FixedEffect_DriftSE_dt
        #FixedEffect_DriftUpperLimit_dt <- FixedEffect_Drift_dt + 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftUpperLimit_dt
        #FixedEffect_DriftLowerLimit_dt <- FixedEffect_Drift_dt - 1.96*FixedEffect_DriftSE_dt; FixedEffect_DriftLowerLimit_dt
        #FixedEffect_DriftZ_dt <- FixedEffect_Drift_dt/FixedEffect_DriftSE_dt; FixedEffect_DriftZ_dt
        #FixedEffect_DriftProb_dt <- round(1-stats::pnorm(abs(FixedEffect_DriftZ_dt),
        #                                                 mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); FixedEffect_DriftProb_dt
        #Q_Drift_dt <- colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt^2)- (colSums(DRIFTPrecision_dt^2 * drift_Coeff_dt))^2 / colSums(DRIFTPrecision_dt^2); Q_Drift_dt
        #H2_Drift_dt <- Q_Drift_dt/(n.studies-1); H2_Drift_dt
        #I2_Drift_dt <- (H2_Drift_dt-1)/H2_Drift_dt*100; I2_Drift_dt
      }
      # Tau square
      T2_DriftWeights3 <- colSums(DRIFTPrecision3^2^2); T2_DriftWeights3 # Borenstein et al., 2007, p. 98
      cDrift3 <- T_DriftWeights3 - T2_DriftWeights3/T_DriftWeights3; cDrift3
      tau2Drift3 <- (Q_Drift3 - (n.studies2-1))/cDrift3; tau2Drift3
      SElnHDrift3 <- c()
      SElnHDrift3[] <- 0
      for (j in 1:(n.latent2^2)) {
        if (Q_Drift3[j] > n.studies2) SElnHDrift3[j] <- 1/2*(log(Q_Drift3[j])-log(n.studies2-1))/((2*Q_Drift3[j])^.5-(2*(n.studies2-1)-1)^.5)
        if (Q_Drift3[j] <= n.studies2) SElnHDrift3[j] <-  (1/(2*(n.studies2-2)) * (1 - 1/(3*(n.studies2-2)^.5)) )^.5
      }
      H2DriftUpperLimit3 <- exp(log(H2_Drift3) + 1.96*SElnHDrift3); H2DriftUpperLimit3
      H2DriftLowerLimit3 <- exp(log(H2_Drift3) - 1.96*SElnHDrift3); H2DriftLowerLimit3
      L3 <- exp(0.5*log(Q_Drift3/(n.studies2-1))-1.96*SElnHDrift3)
      U3 <- exp(0.5*log(Q_Drift3/(n.studies2-1))+1.96*SElnHDrift3)
      I2DriftUpperLimit3 <- (U3^2-1)/U3^2 * 100; I2DriftUpperLimit3
      I2DriftLowerLimit3 <- (L3^2-1)/L3^2 * 100; I2DriftLowerLimit3
      # same for dt
      if (!(is.null(dt))) { # tbd
        #T2_DriftWeights_dt <- colSums(DRIFTPrecision_dt^2^2); T2_DriftWeights_dt # Borenstein et al., 2007, p. 98
        #cDrift_dt <- T_DriftWeights_dt-T2_DriftWeights_dt/T_DriftWeights_dt; cDrift_dt
        #taU3Drift_dt <- (Q_Drift_dt-(n.studies-1))/cDrift_dt; taU3Drift_dt
        #SElnHDrift_dt <- c()
        #SElnHDrift_dt[] <- 0
        #for (j in 1:(n.latent^2)) {
        #  if (Q_Drift_dt[j] > n.studies) SElnHDrift_dt[j] <- 1/2*(log(Q_Drift_dt[j])-log(n.studies-1))/((2*Q_Drift_dt[j])^.5-(2*(n.studies-1)-1)^.5)
        #  if (Q_Drift_dt[j] <= n.studies) SElnHDrift_dt[j] <-  (1/(2*(n.studies-2)) * (1 - 1/(3*(n.studies-2)^.5)) )^.5
        #}
        #H3DRiftUpperLimit_dt <- exp(log(H2_Drift_dt) + 1.96*SElnHDrift_dt); H3DRiftUpperLimit_dt
        #H3DRiftLowerLimit_dt <- exp(log(H2_Drift_dt) - 1.96*SElnHDrift_dt); H3DRiftLowerLimit_dt
        #L_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))-1.96*SElnHDrift_dt)
        #U_dt <- exp(0.5*log(Q_Drift_dt/(n.studies-1))+1.96*SElnHDrift_dt)
        #I3DRiftUpperLimit_dt <- (U_dt^2-1)/U_dt^2 * 100; I3DRiftUpperLimit_dt
        #I3DRiftLowerLimit_dt <- (L_dt^2-1)/L_dt^2 * 100; I3DRiftLowerLimit_dt
      }
    }


    ### Collect Results from fixed effect analysis #####
    if (!(ctmaFitObject[[1]] == "experimental")) MeanOfDriftValues1 <- DriftMeans1
    MeanOfDriftValues2 <- DriftMeans2
    if (ctmaFitObject[[1]] == "experimental") MeanOfDriftValues3 <- DriftMeans3

    newNames <- names(FixedEffect_Drift2) ; newNames
    newNames1 <- gsub("DRIFT ", "fitObj ", newNames); newNames1
    newNames2 <- gsub("DRIFT ", "fitObjMod ", newNames); newNames2

    if (!(ctmaFitObject[[1]] == "experimental")) {
      names(MeanOfDriftValues1) <- newNames1
      names(FixedEffect_Drift1) <- newNames1
      names(FixedEffect_DriftVariance1) <- newNames1
      names(FixedEffect_DriftSE1) <- newNames1
      names(FixedEffect_DriftUpperLimit1) <- newNames1
      names(FixedEffect_DriftLowerLimit1) <- newNames1
      names(FixedEffect_DriftZ1) <- newNames1
      names(FixedEffect_DriftProb1) <- newNames1
      names(tau2Drift1) <- newNames1
      names(Q_Drift1) <- newNames1
      names(H2_Drift1) <- newNames1
      names(H2DriftUpperLimit1) <- newNames1
      names(H2DriftLowerLimit1) <- newNames1
      names(I2_Drift1) <- newNames1
      names(I2DriftUpperLimit1) <- newNames1
      names(I2DriftLowerLimit1) <- newNames1
    }

    names(MeanOfDriftValues2) <- newNames2
    names(FixedEffect_Drift2) <- newNames2
    names(FixedEffect_DriftVariance2) <- newNames2
    names(FixedEffect_DriftSE2) <- newNames2
    names(FixedEffect_DriftUpperLimit2) <- newNames2
    names(FixedEffect_DriftLowerLimit2) <- newNames2
    names(FixedEffect_DriftZ2) <- newNames2
    names(FixedEffect_DriftProb2) <- newNames2
    names(tau2Drift2) <- newNames2
    names(Q_Drift2) <- newNames2
    names(H2_Drift2) <- newNames2
    names(H2DriftUpperLimit2) <- newNames2
    names(H2DriftLowerLimit2) <- newNames2
    names(I2_Drift2) <- newNames2
    names(I2DriftUpperLimit2) <- newNames2
    names(I2DriftLowerLimit2) <- newNames2

    if (ctmaFitObject[[1]] == "experimental") {
      names(MeanOfDriftValues3) <- newNames1
      names(FixedEffect_Drift3) <- newNames1
      names(FixedEffect_DriftVariance3) <- newNames1
      names(FixedEffect_DriftSE3) <- newNames1
      names(FixedEffect_DriftUpperLimit3) <- newNames1
      names(FixedEffect_DriftLowerLimit3) <- newNames1
      names(FixedEffect_DriftZ3) <- newNames1
      names(FixedEffect_DriftProb3) <- newNames1
      names(tau2Drift3) <- newNames1
      names(Q_Drift3) <- newNames1
      names(H2_Drift3) <- newNames1
      names(H2DriftUpperLimit3) <- newNames1
      names(H2DriftLowerLimit3) <- newNames1
      names(I2_Drift3) <- newNames1
      names(I2DriftUpperLimit3) <- newNames1
      names(I2DriftLowerLimit3) <- newNames1
      # replace for easier output formattiing
      MeanOfDriftValues1 <- MeanOfDriftValues3
      FixedEffect_Drift1 <- FixedEffect_Drift3
      FixedEffect_DriftVariance1 <- FixedEffect_DriftVariance3
      FixedEffect_DriftSE1 <- FixedEffect_DriftSE3
      FixedEffect_DriftUpperLimit1 <- FixedEffect_DriftUpperLimit3
      FixedEffect_DriftLowerLimit1 <- FixedEffect_DriftLowerLimit3
      FixedEffect_DriftZ1 <- FixedEffect_DriftZ3
      FixedEffect_DriftProb1 <- FixedEffect_DriftProb3
      tau2Drift1 <- tau2Drift3
      Q_Drift1 <- Q_Drift3
      H2_Drift1 <- H2_Drift3
      H2DriftUpperLimit1 <- H2DriftUpperLimit3
      H2DriftLowerLimit1 <- H2DriftLowerLimit3
      I2_Drift1 <- I2_Drift3
      I2DriftUpperLimit1 <- I2DriftUpperLimit3
      I2DriftLowerLimit1 <- I2DriftLowerLimit3
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

    if (!(is.null(dt))) { #tbd
      #fixedEffectDriftResults <- rbind(fixedEffectDriftResults,
      #                                 FixedEffect_Drift_dt, FixedEffect_DriftVariance_dt, FixedEffect_DriftSE_dt,
      #                                 FixedEffect_DriftUpperLimit_dt, FixedEffect_DriftLowerLimit_dt,
      #                                 FixedEffect_DriftZ_dt, FixedEffect_DriftProb_dt, tau2Drift_dt, Q_Drift_dt, H2_Drift_dt,
      #                                 H2DriftUpperLimit_dt, H2DriftLowerLimit_dt, I2_Drift_dt,
      #                                 I2DriftUpperLimit_dt, I2DriftLowerLimit_dt)
    }

    fixedEffectDriftMessage <- c()
    if ( (any(I2_Drift1 < 0)) | (any(I2_Drift2 < 0)) ) fixedEffectDriftMessage <- "Negative I2 values can be set to 0.0."
    tau2DriftMessage <- c()
    if ( (any(tau2Drift1 < 0)) | (any(tau2Drift1 < 0)) ) tau2DriftMessage <- "Some tau-squared are negative. Random effects cannot be computed. Possibly a small-k-problem."

  }
  #round(fixedEffectDriftResults1, 4)

  ## RANDOM EFFECTS ANALYSIS ##############################################################################

  # Total variance weighting
  {
    if (!(ctmaFitObject[[1]] == "experimental")) {
      Ttot_DriftWeights1 <- 0
      Ttot_DriftMeans1 <- 0
      tau2DriftExtended1 <- do.call(rbind, replicate(n.studies1, tau2Drift1, simplify=FALSE))
      Ttot_DriftWeights1 <-colSums(1/ (DRIFTSE1^2 + tau2DriftExtended1)); Ttot_DriftWeights1
      Ttot_DriftMeans1 <- colSums(DRIFTCoeff1 * 1/ (DRIFTSE1^2 + tau2DriftExtended1)); Ttot_DriftMeans1
      # same for dt # tbd
      if (!(is.null(dt))) {
      }
    }

    Ttot_DriftWeights2 <- 0
    Ttot_DriftMeans2 <- 0
    tau2DriftExtended2 <- do.call(rbind, replicate(n.studies2, tau2Drift2, simplify=FALSE))
    Ttot_DriftWeights2 <-colSums(1/ (DRIFTSE2^2 + tau2DriftExtended2)); Ttot_DriftWeights2
    Ttot_DriftMeans2 <- colSums(DRIFTCoeff2 * 1/ (DRIFTSE2^2 + tau2DriftExtended2)); Ttot_DriftMeans2
    # same for dt # tbd
    if (!(is.null(dt))) {
      #  Ttot_DriftWeights_dt <- 0
      #  Ttot_DriftMeans_dt <- 0
      #  tau2DriftExtended_dt <- do.call(rbind, replicate(n.studies, tau2Drift_dt, simplify=FALSE))
      #  Ttot_DriftWeights_dt <-colSums(1/ (drift_SE_dt^2 + tau2DriftExtended_dt)); Ttot_DriftWeights_dt
      #  Ttot_DriftMeans_dt <- colSums(drift_Coeff_dt * 1/ (drift_SE_dt^2 + tau2DriftExtended_dt)); Ttot_DriftMeans_dt
    }

    if (ctmaFitObject[[1]] == "experimental") {
      Ttot_DriftWeights3 <- 0
      Ttot_DriftMeans3 <- 0
      tau2DriftExtended3 <- do.call(rbind, replicate(n.studies2, tau2Drift3, simplify=FALSE))
      Ttot_DriftWeights3 <-colSums(1/ (DRIFTSE3^2 + tau2DriftExtended3)); Ttot_DriftWeights3
      Ttot_DriftMeans3 <- colSums(DRIFTCoeff3 * 1/ (DRIFTSE3^2 + tau2DriftExtended3)); Ttot_DriftMeans3
      # same for dt # tbd
      if (!(is.null(dt))) {
      }
    }


    ### Random effects results for ctmaFitObject without moderators ####
    if (!(ctmaFitObject[[1]] == "experimental")) {
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
      if (!(is.null(dt))) { # tbd
        #RandomEffecttot_Drift_dt <- Ttot_DriftMeans_dt/Ttot_DriftWeights_dt; RandomEffecttot_Drift_dt
        #RandomEffecttot_DriftVariance_dt <- 1/Ttot_DriftWeights_dt; RandomEffecttot_DriftVariance_dt
        #RandomEffecttot_DriftSE_dt <- RandomEffecttot_DriftVariance_dt^.5; RandomEffecttot_DriftSE_dt
        #RandomEffecttot_DriftUpperLimit_dt <- RandomEffecttot_Drift_dt + 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftUpperLimit_dt
        #RandomEffecttot_DriftLowerLimit_dt <- RandomEffecttot_Drift_dt - 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftLowerLimit_dt
        #RandomEffecttot_DriftZ_dt <- RandomEffecttot_Drift_dt/RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftZ_dt
        #RandomEffecttot_DriftProb_dt <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ_dt),
        #                                                     mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb_dt
        #RandomEffecttot_DriftUpperLimitPI_dt <- RandomEffecttot_Drift_dt + 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftUpperLimitPI_dt
        #RandomEffecttot_DriftLowerLimitPI_dt <- RandomEffecttot_Drift_dt - 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftLowerLimitPI_dt
        #RandomEffectDriftResults_dt <- rbind(RandomEffecttot_Drift_dt, RandomEffecttot_DriftVariance_dt, RandomEffecttot_DriftSE_dt,
        #                                     RandomEffecttot_DriftUpperLimit_dt, RandomEffecttot_DriftLowerLimit_dt,
        #                                     RandomEffecttot_DriftZ_dt, RandomEffecttot_DriftProb_dt,
        #                                     RandomEffecttot_DriftUpperLimitPI_dt, RandomEffecttot_DriftLowerLimitPI_dt)
      }

      ### some corrections for the output
      heterogeneity1 <- fixedEffectDriftResults1[9:16,]; heterogeneity1

      # same for dt # tbd
      if (!(is.null(dt))) {
        #  PET_Drift_dt <- unlist(lapply(PETDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift_dt
        #  PET_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PET_SE_dt <- c(PET_SE_dt, summary(PETDrift_fit_dt[[k]])$coefficients[1,2])
        #  PEESE_Drift_dt <-unlist(lapply(PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift_dt
        #  PEESE_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PEESE_SE_dt <- c(PEESE_SE_dt, summary(PEESEDrift_fit_dt[[k]])$coefficients[1,2])
        #  PET_PEESE_Drift_dt <-unlist(lapply(PET_PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift_dt
        #  PET_PEESE_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PET_PEESE_SE_dt <- c(PET_PEESE_SE_dt, summary(PET_PEESEDrift_fit_dt[[k]])$coefficients[1,2])
        #  WLS_Drift_dt <- unlist(lapply(WLSDrift_fit_dt, function(extract) extract$coefficients)); WLS_Drift_dt
        #  WLS_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) WLS_SE_dt <- c(WLS_SE_dt, summary(WLSDrift_fit_dt[[k]])$coefficients[1,2])
        #  Egger2Drift_results_dt <- matrix(unlist(Egger2Drift_fit_dt), ncol=n.latent^2, nrow=4); Egger2Drift_results_dt
        #  PET_PEESE_DRIFTresults_dt <- rbind(PET_Drift_dt, PET_SE_dt,
        #                                     PEESE_Drift_dt, PEESE_SE_dt,
        #                                     PET_PEESE_Drift_dt, PET_PEESE_SE_dt,
        #                                     WLS_Drift_dt, WLS_SE_dt,
        #                                     Egger2Drift_results_dt)
        #  colnames(PET_PEESE_DRIFTresults_dt) <- colnames(DRIFTCoeff)
        #  rownames(PET_PEESE_DRIFTresults_dt) <- c(rownames(PET_PEESE_DRIFTresults_dt)[1:8], "Egger's b0", "SE(b0)", "T", "p")
        #  ### some corrections for the output
        #  heterogeneity_dt <- fixedEffectDriftResults[24:31,]; heterogeneity_dt # correct! fixedEffectDriftResults_dt does not yet exists
        #  fixedEffectDriftResults_dt <- fixedEffectDriftResults[17:23,]; fixedEffectDriftResults_dt
        #  eggerTest_dt <- PET_PEESE_DRIFTresults_dt[9:12,]; eggerTest_dt
        #  PET_PEESE_DRIFTresults_dt <- PET_PEESE_DRIFTresults_dt[1:8,]; PET_PEESE_DRIFTresults_dt
        #  #
        #  fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults
      }
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
      if (!(is.null(dt))) { # tbd
        #RandomEffecttot_Drift_dt <- Ttot_DriftMeans_dt/Ttot_DriftWeights_dt; RandomEffecttot_Drift_dt
        #RandomEffecttot_DriftVariance_dt <- 1/Ttot_DriftWeights_dt; RandomEffecttot_DriftVariance_dt
        #RandomEffecttot_DriftSE_dt <- RandomEffecttot_DriftVariance_dt^.5; RandomEffecttot_DriftSE_dt
        #RandomEffecttot_DriftUpperLimit_dt <- RandomEffecttot_Drift_dt + 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftUpperLimit_dt
        #RandomEffecttot_DriftLowerLimit_dt <- RandomEffecttot_Drift_dt - 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftLowerLimit_dt
        #RandomEffecttot_DriftZ_dt <- RandomEffecttot_Drift_dt/RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftZ_dt
        #RandomEffecttot_DriftProb_dt <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ_dt),
        #                                                     mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb_dt
        #RandomEffecttot_DriftUpperLimitPI_dt <- RandomEffecttot_Drift_dt + 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftUpperLimitPI_dt
        #RandomEffecttot_DriftLowerLimitPI_dt <- RandomEffecttot_Drift_dt - 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftLowerLimitPI_dt
        #RandomEffectDriftResults_dt <- rbind(RandomEffecttot_Drift_dt, RandomEffecttot_DriftVariance_dt, RandomEffecttot_DriftSE_dt,
        #                                     RandomEffecttot_DriftUpperLimit_dt, RandomEffecttot_DriftLowerLimit_dt,
        #                                     RandomEffecttot_DriftZ_dt, RandomEffecttot_DriftProb_dt,
        #                                     RandomEffecttot_DriftUpperLimitPI_dt, RandomEffecttot_DriftLowerLimitPI_dt)
      }

      ### some corrections for the output
      heterogeneity2 <- fixedEffectDriftResults2[9:16,]; heterogeneity2


      # same for dt # tbd
      if (!(is.null(dt))) {
        #  PET_Drift_dt <- unlist(lapply(PETDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift_dt
        #  PET_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PET_SE_dt <- c(PET_SE_dt, summary(PETDrift_fit_dt[[k]])$coefficients[1,2])
        #  PEESE_Drift_dt <-unlist(lapply(PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift_dt
        #  PEESE_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PEESE_SE_dt <- c(PEESE_SE_dt, summary(PEESEDrift_fit_dt[[k]])$coefficients[1,2])
        #  PET_PEESE_Drift_dt <-unlist(lapply(PET_PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift_dt
        #  PET_PEESE_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PET_PEESE_SE_dt <- c(PET_PEESE_SE_dt, summary(PET_PEESEDrift_fit_dt[[k]])$coefficients[1,2])
        #  WLS_Drift_dt <- unlist(lapply(WLSDrift_fit_dt, function(extract) extract$coefficients)); WLS_Drift_dt
        #  WLS_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) WLS_SE_dt <- c(WLS_SE_dt, summary(WLSDrift_fit_dt[[k]])$coefficients[1,2])
        #  Egger2Drift_results_dt <- matrix(unlist(Egger2Drift_fit_dt), ncol=n.latent^2, nrow=4); Egger2Drift_results_dt
        #  PET_PEESE_DRIFTresults_dt <- rbind(PET_Drift_dt, PET_SE_dt,
        #                                     PEESE_Drift_dt, PEESE_SE_dt,
        #                                     PET_PEESE_Drift_dt, PET_PEESE_SE_dt,
        #                                     WLS_Drift_dt, WLS_SE_dt,
        #                                     Egger2Drift_results_dt)
        #  colnames(PET_PEESE_DRIFTresults_dt) <- colnames(DRIFTCoeff)
        #  rownames(PET_PEESE_DRIFTresults_dt) <- c(rownames(PET_PEESE_DRIFTresults_dt)[1:8], "Egger's b0", "SE(b0)", "T", "p")
        #  ### some corrections for the output
        #  heterogeneity_dt <- fixedEffectDriftResults[24:31,]; heterogeneity_dt # correct! fixedEffectDriftResults_dt does not yet exists
        #  fixedEffectDriftResults_dt <- fixedEffectDriftResults[17:23,]; fixedEffectDriftResults_dt
        #  eggerTest_dt <- PET_PEESE_DRIFTresults_dt[9:12,]; eggerTest_dt
        #  PET_PEESE_DRIFTresults_dt <- PET_PEESE_DRIFTresults_dt[1:8,]; PET_PEESE_DRIFTresults_dt
        #  #
        #  fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults
      }

    }

    if (ctmaFitObject[[1]] == "experimental") {
      RandomEffecttot_Drift3 <- Ttot_DriftMeans2/Ttot_DriftWeights3; RandomEffecttot_Drift3
      RandomEffecttot_DriftVariance3 <- 1/Ttot_DriftWeights3; RandomEffecttot_DriftVariance3
      RandomEffecttot_DriftSE3 <- RandomEffecttot_DriftVariance3^.5; RandomEffecttot_DriftSE3
      RandomEffecttot_DriftUpperLimit3 <- RandomEffecttot_Drift3 + 1.96*RandomEffecttot_DriftSE3; RandomEffecttot_DriftUpperLimit3
      RandomEffecttot_DriftLowerLimit3 <- RandomEffecttot_Drift3 - 1.96*RandomEffecttot_DriftSE3; RandomEffecttot_DriftLowerLimit3
      RandomEffecttot_DriftZ3 <- RandomEffecttot_Drift3/RandomEffecttot_DriftSE3; RandomEffecttot_DriftZ3
      RandomEffecttot_DriftProb3 <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ3),
                                                         mean=c(rep(0, (n.latent2^2))), sd=c(rep(1, (n.latent2^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb3
      RandomEffecttot_DriftUpperLimitPI3 <- RandomEffecttot_Drift3 + 1.96*(tau2Drift3^.5); RandomEffecttot_DriftUpperLimitPI3
      RandomEffecttot_DriftLowerLimitPI3 <- RandomEffecttot_Drift3 - 1.96*(tau2Drift3^.5); RandomEffecttot_DriftLowerLimitPI3
      RandomEffectDriftResults3 <- rbind(RandomEffecttot_Drift3, RandomEffecttot_DriftVariance3, RandomEffecttot_DriftSE3,
                                         RandomEffecttot_DriftUpperLimit3, RandomEffecttot_DriftLowerLimit3,
                                         RandomEffecttot_DriftZ3, RandomEffecttot_DriftProb3,
                                         RandomEffecttot_DriftUpperLimitPI3, RandomEffecttot_DriftLowerLimitPI3)
      # same for dt
      if (!(is.null(dt))) { # tbd
        #RandomEffecttot_Drift_dt <- Ttot_DriftMeans_dt/Ttot_DriftWeights_dt; RandomEffecttot_Drift_dt
        #RandomEffecttot_DriftVariance_dt <- 1/Ttot_DriftWeights_dt; RandomEffecttot_DriftVariance_dt
        #RandomEffecttot_DriftSE_dt <- RandomEffecttot_DriftVariance_dt^.5; RandomEffecttot_DriftSE_dt
        #RandomEffecttot_DriftUpperLimit_dt <- RandomEffecttot_Drift_dt + 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftUpperLimit_dt
        #RandomEffecttot_DriftLowerLimit_dt <- RandomEffecttot_Drift_dt - 1.96*RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftLowerLimit_dt
        #RandomEffecttot_DriftZ_dt <- RandomEffecttot_Drift_dt/RandomEffecttot_DriftSE_dt; RandomEffecttot_DriftZ_dt
        #RandomEffecttot_DriftProb_dt <- round(1-stats::pnorm(abs(RandomEffecttot_DriftZ_dt),
        #                                                     mean=c(rep(0, (n.latent^2))), sd=c(rep(1, (n.latent^2))), log.p=F), digits=digits); RandomEffecttot_DriftProb_dt
        #RandomEffecttot_DriftUpperLimitPI_dt <- RandomEffecttot_Drift_dt + 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftUpperLimitPI_dt
        #RandomEffecttot_DriftLowerLimitPI_dt <- RandomEffecttot_Drift_dt - 1.96*(tau2Drift_dt^.5); RandomEffecttot_DriftLowerLimitPI_dt
        #RandomEffectDriftResults_dt <- rbind(RandomEffecttot_Drift_dt, RandomEffecttot_DriftVariance_dt, RandomEffecttot_DriftSE_dt,
        #                                     RandomEffecttot_DriftUpperLimit_dt, RandomEffecttot_DriftLowerLimit_dt,
        #                                     RandomEffecttot_DriftZ_dt, RandomEffecttot_DriftProb_dt,
        #                                     RandomEffecttot_DriftUpperLimitPI_dt, RandomEffecttot_DriftLowerLimitPI_dt)
      }

      ### some corrections for the output# now re-labeld as "1"
      heterogeneity1 <- fixedEffectDriftResults1[9:16,]; heterogeneity1


      # same for dt # tbd
      if (!(is.null(dt))) {
        #  PET_Drift_dt <- unlist(lapply(PETDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift_dt
        #  PET_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PET_SE_dt <- c(PET_SE_dt, summary(PETDrift_fit_dt[[k]])$coefficients[1,2])
        #  PEESE_Drift_dt <-unlist(lapply(PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift_dt
        #  PEESE_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PEESE_SE_dt <- c(PEESE_SE_dt, summary(PEESEDrift_fit_dt[[k]])$coefficients[1,2])
        #  PET_PEESE_Drift_dt <-unlist(lapply(PET_PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift_dt
        #  PET_PEESE_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) PET_PEESE_SE_dt <- c(PET_PEESE_SE_dt, summary(PET_PEESEDrift_fit_dt[[k]])$coefficients[1,2])
        #  WLS_Drift_dt <- unlist(lapply(WLSDrift_fit_dt, function(extract) extract$coefficients)); WLS_Drift_dt
        #  WLS_SE_dt <- c()
        #  for (k in 1:(n.latent^2)) WLS_SE_dt <- c(WLS_SE_dt, summary(WLSDrift_fit_dt[[k]])$coefficients[1,2])
        #  Egger2Drift_results_dt <- matrix(unlist(Egger2Drift_fit_dt), ncol=n.latent^2, nrow=4); Egger2Drift_results_dt
        #  PET_PEESE_DRIFTresults_dt <- rbind(PET_Drift_dt, PET_SE_dt,
        #                                     PEESE_Drift_dt, PEESE_SE_dt,
        #                                     PET_PEESE_Drift_dt, PET_PEESE_SE_dt,
        #                                     WLS_Drift_dt, WLS_SE_dt,
        #                                     Egger2Drift_results_dt)
        #  colnames(PET_PEESE_DRIFTresults_dt) <- colnames(DRIFTCoeff)
        #  rownames(PET_PEESE_DRIFTresults_dt) <- c(rownames(PET_PEESE_DRIFTresults_dt)[1:8], "Egger's b0", "SE(b0)", "T", "p")
        #  ### some corrections for the output
        #  heterogeneity_dt <- fixedEffectDriftResults[24:31,]; heterogeneity_dt # correct! fixedEffectDriftResults_dt does not yet exists
        #  fixedEffectDriftResults_dt <- fixedEffectDriftResults[17:23,]; fixedEffectDriftResults_dt
        #  eggerTest_dt <- PET_PEESE_DRIFTresults_dt[9:12,]; eggerTest_dt
        #  PET_PEESE_DRIFTresults_dt <- PET_PEESE_DRIFTresults_dt[1:8,]; PET_PEESE_DRIFTresults_dt
        #  #
        #  fixedEffectDriftResults <- fixedEffectDriftResults[1:8,]; fixedEffectDriftResults
      }

    }
    if (ctmaFitObject[[1]] == "experimental") {
      MeanOfDriftValues1 <- MeanOfDriftValues3
      RandomEffecttot_Drift1 <- RandomEffecttot_Drift3
      RandomEffecttot_DriftVariance1 <- RandomEffecttot_DriftVariance3
      RandomEffecttot_DriftSE1 <- RandomEffecttot_DriftSE3
      RandomEffecttot_DriftUpperLimit1 <- RandomEffecttot_DriftUpperLimit3
      RandomEffecttot_DriftLowerLimit1 <- RandomEffecttot_DriftLowerLimit3
      RandomEffecttot_DriftZ1 <- RandomEffecttot_DriftZ3
      RandomEffecttot_DriftProb1 <- RandomEffecttot_DriftProb3
      tau2Drift1 <- tau2Drift3
      Q_Drift1 <- Q_Drift3
      H2_Drift1 <- H2_Drift3
      H2DriftUpperLimit1 <- H2DriftUpperLimit3
      H2DriftLowerLimit1 <- H2DriftLowerLimit3
      I2_Drift1 <- I2_Drift3
      I2DriftUpperLimit1 <- I2DriftUpperLimit3
      I2DriftLowerLimit1 <- I2DriftLowerLimit3
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

  }

  # Collect Results for both fixed effect analysis ######
  {
    #if (is.null(dt)) {

    if (ctmaFitObject[[1]] == "experimental") {
      DRIFTCoeff1 <- DRIFTCoeff3
      DRIFTSE1 <- DRIFTSE3
      RandomEffectDriftResults1 <- RandomEffectDriftResults3
    }
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
    #}
    if (!(is.null(dt))) { # tbd
      #modelResultsList <- list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
      #                         DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE,
      #                         DRIFT_timeScaled=DRIFTCoeff_timeScaled, DIFFUSION_timeScaled=DIFFUSIONCoeff_timeScaled,
      #                         DRIFTSE_timeScaled=DRIFTSE_timeScaled, DIFFUSIONSE_timeScaled=DIFFUSIONSE_timeScaled,
      #                         DRIFT_dt=drift_Coeff_dt,
      #                         DRIFT_dt_SE=drift_SE_dt)
      #summaryList <- list(model="Analysis of Publication Bias & Generalizability",
      #                    estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
      #                                   "Heterogeneity"=round(heterogeneity, digits),
      #                                   "I2 message" = fixedEffectDriftMessage,
      #                                   "Tau2 message" = tau2DriftMessage,
      #                                   "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
      #                                   "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
      #                                   "Egger's tests"=round(eggerTest, digits),
      #                                   #"Egger's tests Alt. Version"= FREAResults,
      #                                   "Z-Curve 2.0 Results:"=zFit,
      #                                   "discreteTimeTimeScale" = dt,
      #                                   "Fixed Effects of DISCRETE TIME Drift Coefficients"=round(fixedEffectDriftResults_dt, digits),
      #                                   "Heterogeneity in DISCRETE TIME "=round(heterogeneity_dt, digits),
      #                                   "Random Effects of DISCRETE TIME Drift Coefficients"=round(RandomEffectDriftResults_dt, digits),
      #                                   "PET-PEESE corrections IN DISCRETE TIME"=round(PET_PEESE_DRIFTresults_dt, digits),
      #                                   "Egger's tests"=round(eggerTest_dt, digits),
      #                                   "Z-Curve 2.0 Results in DISCRTE TIME:"=zFit_dt))
    }
  }




  #fixedEffectDriftResults
  FE <- round(fixedEffectDriftResults[1:8,], digits); FE
  RE <- round(randomEffectDriftResults[1:8,], digits); RE

  Het <- round(fixedEffectDriftResults[9:16,], digits); Het

  fitObjDriftAndSE <- round(cbind(modelResultsList1[[1]], modelResultsList1[[2]]), digits); fitObjDriftAndSE
  colnames(fitObjDriftAndSE) <- paste0("fitObject ", c(names11, names11)); fitObjDriftAndSE
  colnames(fitObjDriftAndSE) <- gsub("DRIFT ", "", colnames(fitObjDriftAndSE)); fitObjDriftAndSE
  colnames(fitObjDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)] <- paste0(colnames(fitObjDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)], " SE"); fitObjDriftAndSE
  #fitObjDriftAndSE

  fitObjModDriftAndSE <- round(cbind(modelResultsList2[[1]], modelResultsList2[[2]]), digits); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE) <- paste0("fitObjectMod ", c(names11, names11)); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)] <- paste0(colnames(fitObjModDriftAndSE)[(n.latent2^2+1):(2*n.latent2^2)], " SE"); fitObjModDriftAndSE
  colnames(fitObjModDriftAndSE) <- gsub("DRIFT ", "", colnames(fitObjModDriftAndSE)); fitObjModDriftAndSE
  #fitObjModDriftAndSE

  # Analysis of Reduction in Heterogeneity by means of moderators #####
  #heterogeneity2[heterogeneity2 < 0] <- NA
  #heterogeneity1[heterogeneity1 < 0] <- NA
  heterogeneity2[heterogeneity2 < 0] <- .00001 # CHD 1.9. 2023
  heterogeneity1[heterogeneity1 < 0] <- .00001 # CHD 1.9. 2023
  redHet <- round(heterogeneity2/heterogeneity1, digits)[c(2,3,6),]*100; redHet
  colnames(redHet) <- names11; redHet
  colnames(redHet) <- gsub("DRIFT ", "", colnames(redHet)); redHet
  rownames(redHet) <- gsub("_Drift2", "", rownames(redHet)); redHet
  rownames(redHet) <- paste0(rownames(redHet), " pct. of initial het."); redHet

  results <- list(activeDirectory=activeDirectory,
                  plot.type=NULL, model.type="BiG",
                  n.studies=n.studies2,
                  n.latent=n.latent2,
                  studyList=list(ctmaFitObject, ctmaFitObjectMod),
                  summary=list(fixedEffects=FE, randomEffects=RE,
                               fitObj_DriftAndSE=fitObjDriftAndSE, fitObjMod_DriftAndSE=fitObjModDriftAndSE,
                               heterogeneity=Het,
                               HeterogeneityReduction=redHet,
                               Note="Negative I2 values were set to .00001 for computation of reduction in heterogeneity."))

  class(results) <- "CoTiMAFit"

  invisible(results)
} ### END function definition
