#' ctmaMMtoCINT
#'
#' @description Compute covariance of CINT-based random intercepts obtained with a MANIFESTMEANS model specification and vice versa.
#'
#' @param ctmaFitObject fit object created with ctmaInit or ctmaFit (or with ctStanFit)
#' @param digits digits used for rounding
#' @param undoTimeScaling if FALSE (default), results will correspond to the $randomIntercepts part of the summary of the equivalent model with indVarying="CINT"
#' @param cintCov Covariance matrix of continuous time intercepts and T0MEANS (default NULL). Requires drift matrix of same dimensions to be supplied, too.
#' @param mmCov Covariance matrix of manifest means and T0MEANS (default NULL). Requires drift matrix of same dimensions to be supplied, too.
#' @param drift matrix of drift parameters (default NULL). When drift together with cintCov or mmCov is supplied, any ctmaFitObject supplied is ignored.
#'
#' @importFrom ctsem ctCollapse ctExtract
#' @importFrom stats quantile cov2cor
#'
#' @examples
#' \dontrun{
#' RI_cov <- ctmaMMtoCINT(ctmaFitObject=CoTiMAFullFit_3)
#' print(RI_cov)
#' }
#'
#' @export ctmaMMtoCINT
#'
#' @return returns covariance of CINT-based random intercepts.
#'
ctmaMMtoCINT <- function(ctmaFitObject=NULL, undoTimeScaling=FALSE, digits=4,
                         cintCov=NULL, mmCov=NULL, drift=NULL) {
  #cintCov=NULL
  #mmCov=NULL
  #drift=NULL

  #{
  if (any(unlist(lapply(list(cintCov, mmCov, drift), function(x) !is.null(x))) == TRUE)) {
    if (is.null(drift)) {
      ErrorMsg <- "You provided a covariance matrix of cints or manifestMeans with T0 variables. I need a drift matrix, too!"
      stop(ErrorMsg)
    }
    if ( (!is.null(cintCov)) &(!is.null(mmCov)) ) {
      ErrorMsg <- "You provided a covariance matrix of cints and of manifestMeans with T0 variables. Please chose only one of them!"
      stop(ErrorMsg)
    }
    if (!is.null(cintCov))  {
      if ( ( dim(cintCov)[1] != dim(drift)[1]*2)  | ( dim(cintCov)[2] != dim(drift)[2]*2) ) {
        ErrorMsg <- "The covariance matrix of cints with T0 variables has to have the twice as many rows and columns as the drift matrix!"
        stop(ErrorMsg)
      }
    }
    if (!is.null(mmCov))  {
      if ( ( dim(mmCov)[1] != dim(drift)[1]*2)  | ( dim(mmCov)[2] != dim(drift)[2]*2) ) {
        ErrorMsg <- "The covariance matrix of cints with T0 variables has to have the twice as many rows and columns as the drift matrix!!"
        stop(ErrorMsg)
      }
    }

    n.latent <- nrow(drift); n.latent
    UL <- UR <- diag(1, n.latent, n.latent); UL
    LL <- matrix(0, n.latent, n.latent); LL
    LR <- -drift; LR

    if (!is.null(mmCov)) {
      trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
      cintCov <- trans %*% mmCov %*% t(trans)
    }
    if (!is.null(cintCov)) {
      trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
      mmCov <- solve(trans) %*% cintCov %*% t(solve(trans))
    }

    return(list(CINTpopcov_mean=round(cintCov, digits),
                MMpopcov_mean=round(mmCov, digits),
                message1=c("MM refers to latent and manifest mean modeling (traits), CINT refers to modeling intercepts.")))
  } else {

    # if ctStanFit instead of CoTiMA fit object is provided
    if (is(ctmaFitObject)[1] == "ctStanFit") ctmaFitObject$studyFitList <- ctmaFitObject
    # if CoTiMA fit object contains one or more singleStudyFits
    if (is(ctmaFitObject$studyFitList[[1]])[1] == "ctStanFit") {
      n.studies <- length(ctmaFitObject$studyFitList)
    } else {
      n.studies <- 1
    }
    #n.studies
  }

  {
    # if CoTiMA fit object is provided
    if (is(ctmaFitObject)[1] == "CoTiMAFit") {
      arguments <- ctmaFitObject$argumentList
      n.latent <- ctmaFitObject$ctModel$n.latent; n.latent
      arguments$randomIntercepts
      if ( (arguments$randomIntercepts == "MANIFEST") |
           (arguments$randomIntercepts == "CINT") ) n.latent <- n.latent/2
      n.manifest <- ctmaFitObject$ctModel$n.manifest; n.manifest
      if (!is.null(ctmaFitObject$argumentList$checkSingleStudyResults)) {
        ErrorMsg <- "The fit object was probably created by ctmaInit. Please provide an object created with ctmaFit!"
        stop(ErrorMsg)
      }
    }
  }

  {
    if (is(ctmaFitObject)[1] == "ctStanFit") {
      arguments <- ctmaFitObject$ctstanmodelbase
      arguments$randomIntercepts <- "dummy"
      arguments$scaleTime <- 1
      n.latent <- arguments$n.latent; n.latent
      n.manifest <- arguments$n.manifest; n.manifest
      #
      pars <- arguments$pars
      if (all(pars[pars$matrix=="MANIFESTMEANS", "indvarying"] == TRUE)) arguments$indVarying <- "MANIFEST"
      if (all(pars[pars$matrix=="CINT", "indvarying"] == TRUE)) arguments$indVarying <- "CINT"
      if (is.null(arguments$indVarying)) {
        Msg1 <- "The fit object provided used neither CINT nor MANIFESTMEANS to model random intercepts."
        # are twice as many latents as processes among latents?
        tmp1 <- (pars[pars$matrix=="DRIFT", "param"]); tmp1
        tmp1 <- length(which(is.na(tmp1)))/length(tmp1); tmp1
        tmp2 <- (pars[pars$matrix=="LAMBDA", "value"]); tmp2
        tmp2 <- length(which(tmp2 != 0))/length(tmp2); tmp2
        if ((tmp1 == .75) & (tmp2 == .25)) {
          arguments$randomIntercepts <- "CINT"
          arguments$indVarying <- FALSE
          message(Msg1)
          Msg2 <- "The fit object seems to have modelled CINTs as separate processes. I continue!"
          message(Msg2)
        }
        if ((tmp1 == .75) & (tmp2 == .50)) {
          arguments$randomIntercepts <- "MANIFEST"
          arguments$indVarying <- FALSE
          message(Msg1)
          Msg2 <- "The fit object seems to have modelled CINTs as MANIFESTSMEANS AND as separate processes. I continue!"
          message(Msg2)
        }
        if (tmp1 != .75)  {
          ErrorMsg <- Msg1
          stop(ErrorMsg)
        }
        #else {
        #if (tmp != .75) arguments$randomIntercepts <- FALSE # arguments$randomIntercepts is used for alternative model specification in CoTiMA only
      }
    }
    #arguments$indVarying
    #arguments$randomIntercepts
  }

  if (!is.null(arguments$scaleTime)) scaleTime <- arguments$scaleTime else scaleTime <- 1

  {
    if ( (arguments$indVarying == TRUE) | (arguments$indVarying == "MANIFEST") |
         (arguments$randomIntercepts == "MANIFEST") ) {
      mmRI <- TRUE
    } else {
      mmRI <- FALSE
    }
    if ( (arguments$indVarying == "CINT") | (arguments$randomIntercepts == "CINT") ) {
      cintRI <- TRUE
    } else {
      cintRI <- FALSE
    }
    #mmRI; cintRI
  }

  {
    CINTpopcov_mean <- CINTpopcov_sd <- CINTpopcov_2.5 <- CINTpopcov_97.5 <- CINTpopcov_T <- list()
    CINTpopcor_mean <- CINTpopcor_sd <- CINTpopcor_2.5 <- CINTpopcor_97.5 <- CINTpopcor_T <- list()
    MMpopcov_mean <- MMpopcov_sd <- MMpopcov_2.5 <- MMpopcov_97.5 <- MMpopcov_T <- list()
    MMpopcor_mean <- MMpopcor_sd <- MMpopcor_2.5 <- MMpopcor_97.5 <- MMpopcor_T <- list()
  }

  for (i in 1:n.studies) {
    #i <- 1; n.studies
    if (is(ctmaFitObject$studyFitList[[i]])[1] == "ctStanFit") {
      fit <- ctmaFitObject$studyFitList[[i]]
    } else {
      fit <- ctmaFitObject$studyFitList
    }

    e <- ctsem::ctExtract(fit)

    if (undoTimeScaling == TRUE) e$pop_DRIFT <- e$pop_DRIFT * scaleTime

    ( mmRI & (arguments$randomIntercepts != "MANIFEST") )
    if ( mmRI & (arguments$randomIntercepts != "MANIFEST") ) { # if random intercepts are modelled as manifest means instead cint
      #### IDEA: Transformation matrix describing popcov_MM into popcov_cint transformations and then popcov_cint <- trans %*% popcov_mm %*% t(trans) (see: #https://stats.stackexchange.com/questions/113700/covariance-of-a-random-vector-after-a-linear-transformation)
      #print("Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) are calculated based on a model with individually varying manifest means instead of Cints.")
      popcov <- e$popcov
      e$popcov_est <- popcov
      e$popcov_est[!(is.na(e$popcov_est))] <- NA
      e$popcor <- e$popcov_est
      # save MM estimates
      for (j in 1:(dim(e$popcor)[1])) {
        e$popcor[j,,] <- stats::cov2cor(matrix(e$popcov[j,,], n.latent^2, n.latent^2))
      }

      MMpopcov_mean[[i]] <- round(ctsem::ctCollapse(popcov, 1, mean), digits)
      MMpopcov_sd[[i]] <- round(ctsem::ctCollapse(popcov, 1, sd), digits)
      MMpopcov_T[[i]] <- round(MMpopcov_mean[[i]]/MMpopcov_sd[[i]], digits)
      MMpopcov_2.5[[i]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.025), digits)
      MMpopcov_97.5[[i]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.975), digits)
      #
      MMpopcor_mean[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, mean), digits)
      MMpopcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, sd), digits)
      MMpopcor_T[[i]] <- round(MMpopcor_mean[[i]]/MMpopcor_sd[[i]], digits)
      MMpopcor_2.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.025), digits)
      MMpopcor_97.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.975), digits)
      #
      UL <- UR <- diag(1, n.latent, n.latent); UL
      LL <- matrix(0, n.latent, n.latent); LL
      for (k in 1:(dim(e$popcov_est)[1])) {
        LR <- -e$pop_DRIFT[k,,]; LR
        trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
        e$popcov_est[k, , ] <- trans %*% e$popcov[k, , ] %*% t(trans)
      }
      e$popcov <- e$popcor <- e$popcov_est
      for (j in 1:(dim(e$popcor)[1])) {
        e$popcor[j,,] <- stats::cov2cor(matrix(e$popcor[j,,], n.latent^2, n.latent^2))
      }
      message <- "Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) were inferred from a model with individually varying manifest means instead of Cints."
      CINTpopcov_mean[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits)
      CINTpopcov_sd[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, sd), digits)
      CINTpopcov_T[[i]] <- round(CINTpopcov_mean[[i]]/CINTpopcov_sd[[i]], digits)
      CINTpopcov_2.5[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.025), digits)
      CINTpopcov_97.5[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.975), digits)
      #
      CINTpopcor_mean[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, mean), digits)
      CINTpopcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, sd), digits)
      CINTpopcor_T[[i]] <- round(CINTpopcor_mean[[i]]/CINTpopcor_sd[[i]], digits)
      CINTpopcor_2.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.025), digits)
      CINTpopcor_97.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.975), digits)
    }

    ( cintRI & (arguments$randomIntercepts != "CINT") )
    if ( cintRI & (arguments$randomIntercepts != "CINT") ) { # if random intercepts are modelled as manifest means instead cint
      #### REVERSED IDEA: Transformation matrix describing popcov_MM into popcov_cint transformations and then popcov_cint <- trans %*% popcov_mm %*% t(trans) (see: #https://stats.stackexchange.com/questions/113700/covariance-of-a-random-vector-after-a-linear-transformation)
      #popcov <- e$pop_T0cov  # THIS YIELDS WRONG RESULTS
      popcov <- e$popcov
      e$popcov_est <- popcov
      e$popcov_est[!(is.na(e$popcov_est))] <- NA
      # save CINT estimates
      e$popcor <- e$popcov_est
      for (j in 1:(dim(e$popcor)[1])) {
        #e$popcor[j,,] <- stats::cov2cor(matrix(e$pop_T0cov[j,,], n.latent^2, n.latent^2))
        e$popcor[j,,] <- stats::cov2cor(matrix(popcov[j,,], n.latent^2, n.latent^2))
      }
      CINTpopcov_mean[[i]] <- round(ctsem::ctCollapse(popcov, 1, mean), digits)
      CINTpopcov_sd[[i]] <- round(ctsem::ctCollapse(popcov, 1, sd), digits)
      CINTpopcov_T[[i]] <- round(CINTpopcov_mean[[i]]/CINTpopcov_sd[[i]], digits)
      CINTpopcov_2.5[[i]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.025), digits)
      CINTpopcov_97.5[[i]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.975), digits)
      #
      CINTpopcor_mean[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, mean), digits)
      CINTpopcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, sd), digits)
      CINTpopcor_T[[i]] <- round(CINTpopcor_mean[[i]]/CINTpopcor_sd[[i]], digits)
      CINTpopcor_2.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.025), digits)
      CINTpopcor_97.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.975), digits)
      #
      UL <- UR <- diag(1, n.latent, n.latent); UL
      LL <- matrix(0, n.latent, n.latent); LL
      for (k in 1:(dim(e$popcov_est)[1])) {
        LR <- -e$pop_DRIFT[k,,]; LR
        trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
        #e$popcov_est[k, , ] <- solve(trans) %*% e$popcov[k, , ] %*% t(solve(trans)) # changed
        e$popcov_est[k, , ] <- solve(trans) %*% popcov[k, , ] %*% t(solve(trans)) # changed
      }
      e$popcov <- e$popcor <- e$popcov_est
      for (j in 1:(dim(e$popcor)[1])) {
        e$popcor[j,,] <- stats::cov2cor(matrix(e$popcor[j,,], n.latent^2, n.latent^2))
      }
      message <- "Manifest Means (Trait means), T0means (initial latent means), and T0covs (initial (co-)vars) were inferred from a model with individually varying Cints means instead of Manifest Means."
      MMpopcov_mean[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits)
      MMpopcov_sd[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, sd), digits)
      MMpopcov_T[[i]] <- round(MMpopcov_mean[[i]]/MMpopcov_sd[[i]], digits)
      MMpopcov_2.5[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.025), digits)
      MMpopcov_97.5[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.975), digits)
      #
      MMpopcor_mean[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, mean), digits)
      MMpopcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, sd), digits)
      MMpopcor_T[[i]] <- round(MMpopcor_mean[[i]]/MMpopcor_sd[[i]], digits)
      MMpopcor_2.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.025), digits)
      MMpopcor_97.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.975), digits)
    }


    ########################################################################################################
    ###################### GENERAL COMPUTATIONS FOR REAL CoTiMA randomIntercept MODELS #####################
    ########################################################################################################

    ( (arguments$randomIntercepts == "MANIFEST") | (arguments$randomIntercepts == "CINT") )
    if ( (arguments$randomIntercepts == "MANIFEST") | (arguments$randomIntercepts == "CINT") ) { # if random intercepts are modelled as separate process
      # parNames <- ctsem:::getparnames(fitStanctModel); parNames
      # since getparnames is not exported, I took part fo the function and replicated it here
      ms <- fit$setup$matsetup
      #fit$ctstanmodelbase$pars
      indices <- ms$when %in% c(0, -1) & ms$param > 0 & ms$copyrow < 1
      pars <- data.frame(parnames = ms$parname[indices], parindices = ms$param[indices])
      pars <- pars[!duplicated(pars$parnames), ]
      pars <- pars[order(pars$parindices), ]
      parNames <- pars[pars$parindices > 0, 1]; parNames

      #targetCols <- grep("ov", parNames); targetCols
      #targetRaws <- fit$stanfit$rawposterior[, targetCols]
      #colnames(targetRaws) <- parNames[grep("ov", parNames)]
      #rawT0varTmp <- targetRaws
      #head(rawT0varTmp)

      # tform T0variances = random intercept covariances - just as a check that this is identical to ctmaInit
      n.studies2 <- length(ctmaFitObject$argumentList$primaryStudyList$studyNumbers); n.studies2
      T0COVCoeff <- T0COVCoeffMean <- T0COVCoeffSD <- list()
      n.mod.values.to.plot <- length(colnames(fit$data$tipredsdata)[1:(n.studies2-1)])+1; n.mod.values.to.plot
      #
      modPos <- 1:(n.studies2-1); modPos
      TIpred.values <- matrix(fit$data$tipredsdata[, modPos], ncol=length(modPos))
      #head(TIpred.values)
      effectCodingWeights <- unique(TIpred.values); effectCodingWeights
      #
      tmp1 <- which(!(is.na(fit$ctstanmodelbase$pars$param))); tmp1
      n.main.effects <- length(tmp1); n.main.effects # used later
      tmpPars <- fit$ctstanmodelbase$pars[tmp1,]; tmpPars
      T0varPos <- which(tmpPars$matrix == "T0VAR"); T0varPos
      rawT0varTmp <- fit$stanfit$rawposterior[ , T0varPos]
      targetNames <-fit$ctstanmodelbase$pars[fit$ctstanmodelbase$pars$matrix=="T0VAR", "param"]; targetNames
      targetNames <- targetNames[!is.na(targetNames)]; targetNames
      colnames(rawT0varTmp) <- targetNames
      #head(rawT0varTmp)
      #
      tmpNames <- paste0("RI Covriances for Study No ", unlist(lapply(ctmaFitObject$studyList, function(x) x$originalStudyNo)), "."); tmpNames
      #
      alternative <- 0
      if (alternative == 0) {
        TIpredEffTmp <- fit$stanfit$transformedparsfull$TIPREDEFFECT[,T0varPos, modPos]; head(TIpredEffTmp)
      }

      # alternative to above
      if (alternative == 1) {
        tmp2 <- grep("_effect", colnames(fit$ctstanmodelbase$pars)); tmp2 # all TIpred effects pos
        tmp2 <- tmp2[modPos]; tmp2                                        # mod TIpred effects pos
        tmp3 <- fit$ctstanmodelbase$pars[fit$ctstanmodelbase$pars$matrix == "T0VAR", tmp2]; tmp3 # mod TIpred effects logical
        n.tipred.params <- length(which(tmp3 == TRUE)); n.tipred.params
        n.all.params <- dim(fit$stanfit$rawposterior)[2]; n.all.params
        tmpData <- fit$stanfit$rawposterior[, (n.all.params-n.tipred.params+1):n.all.params] # all TIpred effects finishsapmples
        finishsamples <- nrow(tmpData); finishsamples
        TIpredEffTmp <- array(data=fit$stanfit$rawposterior[, (n.all.params-n.tipred.params+1):n.all.params],
                              dim=c(finishsamples, n.tipred.params/length(modPos), length(modPos)))
        #(TIpredEffTmp[1:5, 1, 2])
        #head(fit$stanfit$rawposterior[1:5,46])
      }

      #  n.studies2 should be > 2 for regular computation because otherwise there are more effectCodingWeights
      if (alternative != 1) {
        for (d in 1:nrow(effectCodingWeights)) {
          if (n.studies2 == 2) {
            tmp2 <- apply(TIpredEffTmp %*% t(effectCodingWeights[d, ]), 1, sum); tmp2
          } else {
            tmp2 <- apply(TIpredEffTmp %*% effectCodingWeights[d, ], 1, sum); tmp2
          }
          rawT0varTmp
          tmp1 <- t(apply(rawT0varTmp, 1 , function(x) x+tmp2))
          T0COVCoeff[[tmpNames[d]]] <- tmp1
        }
      }

      # alternative
      if (alternative == 1) {
        for (d in 1:nrow(effectCodingWeights)) {
          if (n.studies2 == 2) {
            tmp2 <- apply(TIpredEffTmp[1,,] * t(effectCodingWeights[d, ]), 1, sum); tmp2
            for (k in 2:finishsamples) {
              tmp2 <- rbind(tmp2, apply(TIpredEffTmp[k,,] * t(effectCodingWeights[d, ]), 1, sum)); tmp2
            }
          } else {
            tmp2 <- apply(TIpredEffTmp[1,,] %*% effectCodingWeights[d, ], 1, sum); tmp2
            for (k in 2:finishsamples) {
              tmp2 <- rbind(apply(TIpredEffTmp[k,,] * effectCodingWeights[d, ], 1, sum)); tmp2
            }
          }
          rawT0varTmp
          tmp1 <- t(apply(rawT0varTmp, 1 , function(x) x+tmp2))
          T0COVCoeff[[tmpNames[d]]] <- tmp1
        }
      }

      #str(T0COVCoeff)
      #savedOlds <- lapply(T0COVCoeff, function(x) apply(x, 2, mean))
      #round(lapply(T0COVCoeff, function(x) apply(x, 2, mean))[[3]], 4)
      #round(savedOlds[[3]], 4)
      #initFit_MM$studyFitList[[3]]$ctstanmodelbase$pars # 12-14
      #round(ctCollapse(initFit_MM$studyFitList[[3]]$stanfit$rawposterior,1, mean), 4)#[12:14]
      #s3 <- summary(initFit_MM$studyFitList[[3]])
      #s3$popsd

      tmp1a <- fit$ctstanmodelbase$pars[, "transform"]; tmp1a
      tmp1b <- fit$ctstanmodelbase$pars[, "param"]; tmp1b
      tmp1c <- grep("ov_", tmp1b); tmp1c
      transforms <- tmp1a[tmp1c]; transforms
      # compute tformed T0COVCoeff
      for (k in 1:(length(T0COVCoeff))) {
        counter <- 0
        for (l in 1:(n.latent^2)) {
          for (m in 1:l) {
            counter <- counter + 1
            paramTmp <- T0COVCoeff[[k]][, counter]; paramTmp
            for (p in 1:length(paramTmp)) {
              param <- paramTmp[p]
              T0COVCoeff[[k]][p , counter] <- eval(parse(text=transforms[counter])); T0COVCoeff[[k]][p, counter]
            }
          }
        }
      }
      #lapply(T0COVCoeff, function(x) apply(x, 2, mean))
      #s3$popsd

      CINTpopcov_mean[[i]] <- CINTpopcov_sd[[i]] <- CINTpopcov_T[[i]] <-
        CINTpopcov_2.5[[i]] <- CINTpopcov_97.5[[i]] <- list()
      CINTpopcor_mean[[i]] <- CINTpopcor_sd[[i]] <- CINTpopcor_T[[i]] <-
        CINTpopcor_2.5[[i]] <- CINTpopcor_97.5[[i]] <- list()
      MMpopcov_mean[[i]] <- MMpopcov_sd[[i]] <- MMpopcov_T[[i]] <-
        MMpopcov_2.5[[i]] <- MMpopcov_97.5[[i]] <- list()
      MMpopcor_mean[[i]] <- MMpopcor_sd[[i]] <- MMpopcor_T[[i]] <-
        MMpopcor_2.5[[i]] <- MMpopcor_97.5[[i]] <- list()

    } # end if ( (arguments$randomIntercepts == "MANIFEST") | (arguments$randomIntercepts == "CINT") )
    #apply(T0COVCoeff[[1]], 2, mean)

    ############### SPECIFIC COMPUTATIONS FOR REAL CoTiMA randomIntercept MANIFEST MODELS #################

    # This one is mainly taken from ctmaFit and adapted
    ( mmRI & (arguments$randomIntercepts == "MANIFEST") )
    if ( mmRI & (arguments$randomIntercepts == "MANIFEST") ) { # if random intercepts are modelled as manifest means instead cint
      for (j in 1:n.studies2) {
        #j <- 1
        # get the MM estimates
        popcov_tmp <- T0COVCoeff[[j]]
        popcov <- array(NA, dim=c(nrow(popcov_tmp), n.latent^2, n.latent^2))
        popcov_est <- popcov
        popcov_est[!(is.na(popcov_est))] <- NA
        popcor <- popcov_est
        popcor_est <- popcov_est
        for (k in 1:nrow(popcov_tmp)) { # non-transformed covariances & correlations
          upperTri <- popcov_tmp[k,]; upperTri
          popcov[k,,][upper.tri(popcov[k,,], diag=T)] <- upperTri
          popcov[k,,][lower.tri(popcov[k,,])] <- t(popcov[k,,][upper.tri(popcov[k,,])])
          popcor[k,,] <- stats::cov2cor(matrix(popcov[k,,], n.latent^2, n.latent^2))
        }
        # save the MM estimate
        MMpopcov_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, mean), digits)
        MMpopcov_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, sd), digits)
        MMpopcov_T[[i]][[j]] <- round(MMpopcov_mean[[i]][[j]]/MMpopcov_sd[[i]][[j]], digits)
        MMpopcov_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.025), digits)
        MMpopcov_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.975), digits)
        #
        MMpopcor_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, mean), digits)
        MMpopcor_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, sd), digits)
        MMpopcor_T[[i]][[j]] <- round(MMpopcor_mean[[i]][[j]]/MMpopcor_sd[[i]][[j]], digits)
        MMpopcor_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, stats::quantile, probs=.025), digits)
        MMpopcor_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, stats::quantile, probs=.975), digits)

        # now the transformations MMtoCINT
        popcov_est <- popcor_est <- array(NA, dim=c(nrow(popcov), n.latent^2, n.latent^2))
        UL <- UR <- diag(1, n.latent, n.latent); UL
        LL <- matrix(0, n.latent, n.latent); LL
        for (k in 1:(dim(popcov_est)[1])) {
          LR <- -e$pop_DRIFT[k,1:n.latent,1:n.latent]; LR
          trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
          popcov_est[k, , ] <- trans %*% popcov[k,,] %*% t(trans)
        }
        for (k in 1:(dim(popcov_est)[1])) {
          popcor_est[k,,] <- stats::cov2cor(matrix(popcov_est[k,,], n.latent^2, n.latent^2))
        }
        message <- "Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) were inferred from a model with individually varying manifest means instead of Cints."
        CINTpopcov_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, mean), digits)
        CINTpopcov_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, sd), digits)
        CINTpopcov_T[[i]][[j]] <- round(CINTpopcov_mean[[i]][[j]]/CINTpopcov_sd[[i]][[j]], digits)
        CINTpopcov_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, stats::quantile, probs=.025), digits)
        CINTpopcov_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, stats::quantile, probs=.975), digits)
        #
        CINTpopcor_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, mean), digits)
        CINTpopcor_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, sd), digits)
        CINTpopcor_T[[i]][[j]] <- round(CINTpopcor_mean[[i]][[j]]/CINTpopcor_sd[[i]][[j]], digits)
        CINTpopcor_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, stats::quantile, probs=.025), digits)
        CINTpopcor_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, stats::quantile, probs=.975), digits)
      }
    }
    CINTpopcor_97.5

    # This one is ALSO mainly taken from ctmaFit and adapted
    ( cintRI & (arguments$randomIntercepts == "CINT") )
    if ( cintRI & (arguments$randomIntercepts == "CINT") ) { # if random intercepts are modelled as manifest means instead cint
      for (j in 1:n.studies2) {
        #j <- 1
        # get CINT estimates
        popcov_tmp <- T0COVCoeff[[j]]
        popcov <- array(NA, dim=c(nrow(popcov_tmp), n.latent^2, n.latent^2))
        popcov_est <- popcov
        popcov_est[!(is.na(popcov_est))] <- NA
        popcor <- popcov_est
        popcor_est <- popcov_est
        for (k in 1:nrow(popcov_tmp)) { # non-transformed covariances & correlations
          upperTri <- popcov_tmp[k,]; upperTri
          popcov[k,,][upper.tri(popcov[k,,], diag=T)] <- upperTri
          #popcov[k,,][lower.tri(popcov[k,,])] <- t(popcov[k,,][upper.tri(popcov[k,,])])
          popcov[k,,][lower.tri(popcov[k,,])]  <- t(popcov[k,,])[lower.tri(popcov[k,,])]
          popcor[k,,] <- stats::cov2cor(matrix(popcov[k,,], n.latent^2, n.latent^2))
        }
        # save CINT estimates
        CINTpopcov_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, mean), digits)
        CINTpopcov_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, sd), digits)
        CINTpopcov_T[[i]][[j]] <- round(CINTpopcov_mean[[i]][[j]]/CINTpopcov_sd[[i]][[j]], digits)
        CINTpopcov_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.025), digits)
        CINTpopcov_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov, 1, stats::quantile, probs=.975), digits)
        #
        CINTpopcor_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, mean), digits)
        CINTpopcor_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, sd), digits)
        CINTpopcor_T[[i]][[j]] <- round(CINTpopcor_mean[[i]][[j]]/CINTpopcor_sd[[i]][[j]], digits)
        CINTpopcor_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, stats::quantile, probs=.025), digits)
        CINTpopcor_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor, 1, stats::quantile, probs=.975), digits)
        #
        #popcov_est <- popcor_est <- array(NA, dim=c(nrow(popcov), n.latent^2, n.latent^2))
        UL <- UR <- diag(1, n.latent, n.latent); UL
        LL <- matrix(0, n.latent, n.latent); LL
        for (k in 1:(dim(popcov_est)[1])) {  # transformed covariances & correlations
          LR <- -e$pop_DRIFT[k,1:n.latent,1:n.latent]; LR
          trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
          popcovTmp <- popcov[k,,]
          popcov_est[k,,] <- solve(trans) %*% popcovTmp %*% t(solve(trans))
          if (any(diag(popcov_est[k,,]) < 0)) {
            ErrorMsg <- c(paste0("Encountered a negative variances estimate for Study #", j," and finishssample #", k,". Have to stop!"))
            stop(ErrorMsg)
          }
          popcor_est[k,,] <- stats::cov2cor(matrix(popcov_est[k,,], n.latent^2, n.latent^2))
        }

        # save MM estimates
        message <- "MANIFEST MEANS (traits), T0means (initial means), and T0covs (initial (co-)vars) were inferred from a model with individually varying Cints."
        MMpopcov_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, mean), digits)
        MMpopcov_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, sd), digits)
        MMpopcov_T[[i]][[j]] <- round(MMpopcov_mean[[i]][[j]]/MMpopcov_sd[[i]][[j]], digits)
        MMpopcov_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, stats::quantile, probs=.025), digits)
        MMpopcov_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcov_est, 1, stats::quantile, probs=.975), digits)
        #
        MMpopcor_mean[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, mean), digits)
        MMpopcor_sd[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, sd), digits)
        MMpopcor_T[[i]][[j]] <- round(MMpopcor_mean[[i]][[j]]/MMpopcor_sd[[i]][[j]], digits)
        MMpopcor_2.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, stats::quantile, probs=.025), digits)
        MMpopcor_97.5[[i]][[j]] <- round(ctsem::ctCollapse(popcor_est, 1, stats::quantile, probs=.975), digits)
      }
    }
  }

  return(list(CINTpopcov_mean=CINTpopcov_mean, CINTpopcov_sd=CINTpopcov_sd, CINTpopcov_T=CINTpopcov_T,
              CINTpopcov_2.5=CINTpopcov_2.5, CINTpopcov_97.5=CINTpopcov_97.5,
              CINTpopcor_mean=CINTpopcor_mean, CINTpopcor_sd=CINTpopcor_sd, CINTpopcor_T=CINTpopcor_T,
              CINTpopcor_2.5=CINTpopcor_2.5, CINTpopcor_97.5=CINTpopcor_97.5,
              MMpopcov_mean=MMpopcov_mean, MMpopcov_sd=MMpopcov_sd, MMpopcov_T=MMpopcov_T,
              MMpopcov_2.5=MMpopcov_2.5, MMpopcov_97.5=MMpopcov_97.5,
              MMpopcor_mean=MMpopcor_mean, MMpopcor_sd=MMpopcor_sd, MMpopcor_T=MMpopcor_T,
              MMpopcor_2.5=MMpopcor_2.5, MMpopcor_97.5=MMpopcor_97.5,
              message1=message,
              message2=c("MM refers to latent and manifest mean modeling (traits), CINT refers to modeling intercepts.")))
}
