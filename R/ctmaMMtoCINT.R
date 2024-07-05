#' ctmaMMtoCINT
#'
#' @description Compute covariance of CINT-based random intercepts obtained with a MANIFESTMEANS model specification
#'
#' @param ctmaFitObject fit object created with ctmaInit or ctmaFitObject
#' @param digits digits used for rounding
#' @param undoTimeScaling if FALSE (default), results will correspond to the $randomIntercepts part of the summary of the equivalent model with indVarying="CINT"
#'
#' @importFrom ctsem ctCollapse
#' @importFrom stats quantile cov2cor
#'
#' @examples
#' \donttest{
#' RI_cov <- ctmaMMtoCINT(ctmaFitObject=CoTiMAFullFit_3)
#' print(RI_cov)
#' }
#'
#' @export ctmaMMtoCINT
#'
#' @return returns covariance of CINT-based random intercepts.
#'
ctmaMMtoCINT <- function(ctmaFitObject=NULL, undoTimeScaling=FALSE, digits=4) {
  # if ctStanFit instead of CoTiMA fit object is provided
  if (class(ctmaFitObject) == "ctStanFit") ctmaFitObject$studyFitList <- ctmaFitObject
  # if CoTiMA fit object contains one or more singleStudyFits
  if (class(ctmaFitObject$studyFitList[[1]]) == "ctStanFit") {
    n.studies <- length(ctmaFitObject$studyFitList)
  } else {
    n.studies <- 1
  }
  # if CoTiMA fit object is provided
  if (class(ctmaFitObject) == "CoTiMAFit") {
    #arguments <- ctmaFitObject$ctModel
    arguments <- ctmaFitObject$argumentList
    n.latent <- ctmaFitObject$ctModel$n.latent; n.latent
    if ( (arguments$randomIntercepts == "MANIFEST") |
         (arguments$randomIntercepts == "CINT") ) n.latent <- n.latent/2
    n.manifest <- ctmaFitObject$ctModel$n.manifest; n.manifest
    #digits <- arguments$digits; digits
  }
  if (class(ctmaFitObject) == "ctStanFit") {
    arguments <- ctmaFitObject$ctstanmodelbase
    arguments$scaleTime <- 1
    n.latent <- arguments$n.latent; n.latent
    n.manifest <- arguments$n.manifest; n.manifest
    #digits <- 4; digits
    #
    pars <- arguments$pars
    #if (all(pars[pars$matrix=="MANIFESTMEANS", "indvarying"] == TRUE)) arguments$indVarying <- TRUE
    if (all(pars[pars$matrix=="MANIFESTMEANS", "indvarying"] == TRUE)) arguments$indVarying <- "MANIFEST"
    if (all(pars[pars$matrix=="CINT", "indvarying"] == TRUE)) arguments$indVarying <- "CINT"
    if (is.null(arguments$indVarying)) {
      ErrorMsg <- "The fit object provided used neither CINT nor MANIFESTMEANS to model random intercepts."
      stop(ErrorMsg)
    }
  }

  popcov_mean <- popcov_sd <- popcov_2.5 <- popcov_97.5 <- popcov_T <- list()
  popcor_mean <- popcor_sd <- popcor_2.5 <- popcor_97.5 <- popcor_T <- list()
  if (!is.null(arguments$scaleTime)) scaleTime <- arguments$scaleTime else scaleTime <- 1
  if ( (arguments$indVarying == TRUE) | (arguments$indVarying == "MANIFEST") |
       (arguments$randomIntercepts == "MANIFEST") ) mmRI <- TRUE else mmRI <- FALSE
  if ( (arguments$indVarying == "CINT") | (arguments$randomIntercepts == "CINT") ) {
    ErrorMsg <- "The fit object provided used CINT rather the MANIFESTMEANS to model random intercepts."
    stop(ErrorMsg)
  }

  for (i in 1:n.studies) {
    #i <- 1; n.studies
    if (class(ctmaFitObject$studyFitList[[i]]) == "ctStanFit") {
      fit <- ctmaFitObject$studyFitList[[i]]
    } else {
      fit <- ctmaFitObject$studyFitList
    }
    e <- ctsem::ctExtract(fit)
    if (undoTimeScaling == TRUE) e$pop_DRIFT <- e$pop_DRIFT * scaleTime

    # get random intercept stats
    tmp1 <- ctsem::ctCollapse(e$pop_CINT, 1, mean); tmp1
    tmp2 <- ctsem::ctCollapse(e$pop_CINT, 1, sd); tmp2
    tmp3 <- ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, mean); tmp3
    tmp4 <- ctsem::ctCollapse(e$pop_MANIFESTMEANS, 1, sd); tmp4

    if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
      e$pop_T0MEANS_est <- e$pop_T0MEANS[ ,1:n.latent, ]        # 1.n.latent => eliminates last diminsion, which is not required
      e$pop_T0MEANS_est[!(is.na(e$pop_T0MEANS_est))] <- NA
      e$pop_T0MEANS <- e$pop_T0MEANS[ ,1:n.latent, ]
      e$pop_MANIFESTMEANS <- e$pop_MANIFESTMEANS[ ,1:n.latent,]
      for (j in 1:(dim(e$pop_T0MEANS_est)[1])) e$pop_T0MEANS_est[j,] <- e$pop_T0MEANS[j,] + e$pop_MANIFESTMEANS[j,]
      e$pop_T0MEANS <- e$pop_T0MEANS_est
      #e$pop_MANIFESTMEANS_backup <- e$pop_MANIFESTMEANS
      e$pop_MANIFESTMEANS[e$pop_MANIFESTMEANS != 0] <- 0
    }
    #
    if ( mmRI & (arguments$randomIntercepts != "MANIFEST") ) { # if random intercepts are modelled as manifest means instead cint
      #### IDEA: Transformation matrix describing popcov_MM into popcov_cint transformations and then popcov_cint <- trans %*% popcov_mm %*% t(trans) (see: #https://stats.stackexchange.com/questions/113700/covariance-of-a-random-vector-after-a-linear-transformation)
      #print("Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) are calculated based on a model with individually varying manifest means instead of Cints.")
      e$popcov_est <- e$popcov
      e$popcov_est[!(is.na(e$popcov_est))] <- NA
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
      popcov_mean[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits)
      popcov_sd[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, sd), digits)
      popcov_T[[i]] <- round(popcov_mean[[i]]/popcov_sd[[i]], digits)
      popcov_2.5[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.025), digits)
      popcov_97.5[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.975), digits)
      #
      popcor_mean[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, mean), digits)
      popcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, sd), digits)
      popcor_T[[i]] <- round(popcor_mean[[i]]/popcor_sd[[i]], digits)
      popcor_2.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.025), digits)
      popcor_97.5[[i]] <- round(ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.975), digits)
    }
    # This one is mainly taken from ctmaFit and adapted
    ( mmRI & (arguments$randomIntercepts == "MANIFEST") )
    if ( mmRI & (arguments$randomIntercepts == "MANIFEST") ) { # if random intercepts are modelled as manifest means instead cint
      # parNames <- ctsem:::getparnames(fitStanctModel); parNames
      # since getparnames is not exported, I took part fo the function and replicated it here
      ms <- fit$setup$matsetup
      indices <- ms$when %in% c(0, -1) & ms$param > 0 & ms$copyrow < 1
      pars <- data.frame(parnames = ms$parname[indices], parindices = ms$param[indices])
      pars <- pars[!duplicated(pars$parnames), ]
      pars <- pars[order(pars$parindices), ]
      parNames <- pars[pars$parindices > 0, 1]

      targetCols <- grep("ov", parNames); targetCols
      targetRaws <- fit$stanfit$rawposterior[, targetCols]
      colnames(targetRaws) <- parNames[grep("ov", parNames)]
      rawT0varTmp <- targetRaws
      # tform T0variances = random intercept covariances - just as a check that this is identical to ctmaInit
      {
        n.studies2 <- length(ctmaFitObject$argumentList$primaryStudyList$studyNumbers); n.studies2
        T0COVCoeff <- T0COVCoeffMean <- T0COVCoeffSD <- list()
        n.mod.values.to.plot <- length(colnames(fit$data$tipredsdata)[1:(n.studies2-1)])+1; n.mod.values.to.plot
        modPos <- 1:(n.studies2-1); modPos
        TIpred.values <- matrix(fit$data$tipredsdata[, modPos], ncol=length(modPos))
        effectCodingWeights <- unique(TIpred.values); effectCodingWeights
        #
        tmp1 <- which(!(is.na(fit$ctstanmodelbase$pars$param))); tmp1
        tmpPars <- fit$ctstanmodelbase$pars[tmp1,]; tmpPars
        T0varPos <- which(tmpPars$matrix == "T0VAR"); T0varPos
        rawT0varTmp <- fit$stanfit$rawposterior[ , T0varPos]
        #
        tmpNames <- paste0("RI Covriances for Study No ", unlist(lapply(ctmaFitObject$studyList, function(x) x$originalStudyNo)), "."); tmpNames
        #
        TIpredEffTmp <- fit$stanfit$transformedparsfull$TIPREDEFFECT[,T0varPos, modPos]; TIpredEffTmp
        #  n.studies2 should be > 2 for regular computation because otherwise there are more effectCodingWeights
        for (d in 1:nrow(effectCodingWeights)) {
          if (n.studies2 == 2) {
            tmp2 <- apply(TIpredEffTmp %*% t(effectCodingWeights[d, ]), 1, sum); tmp2
          } else {
            tmp2 <- apply(TIpredEffTmp %*% effectCodingWeights[d, ], 1, sum); tmp2
          }
          tmp1 <- t(apply(rawT0varTmp, 1 , function(x) x+tmp2))
          T0COVCoeff[[tmpNames[d]]] <- tmp1
        }
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
        # now the transformations MMtoCINT
        for (j in 1:n.studies2) {
          popcov <- T0COVCoeff[[j]]
          popcov_est <- popcor_est <- array(NA, dim=c(nrow(popcov), n.latent^2, n.latent^2))
          UL <- UR <- diag(1, n.latent, n.latent); UL
          LL <- matrix(0, n.latent, n.latent); LL
          for (k in 1:(dim(popcov_est)[1])) {
            LR <- -e$pop_DRIFT[k,1:n.latent,1:n.latent]; LR
            trans <- rbind(cbind(UL, UR), cbind(LL, LR)); trans
            popcovTmp <- matrix(NA, n.latent^2, n.latent^2); popcovTmp
            popcovTmp[upper.tri(popcovTmp, diag=T)] <- popcov[k,]; popcovTmp
            popcovTmp[lower.tri(popcovTmp)]  <- t(popcovTmp)[lower.tri(popcovTmp)]
            popcov_est[k, , ] <- trans %*% popcovTmp %*% t(trans)
          }
          for (k in 1:(dim(popcor)[1])) {
            popcor_est[k,,] <- stats::cov2cor(matrix(popcov_est[k,,], n.latent^2, n.latent^2))
          }
          message <- "Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) were inferred from a model with individually varying manifest means instead of Cints."
          popcov_mean[[j]] <- round(ctsem::ctCollapse(popcov_est, 1, mean), digits)
          popcov_sd[[j]] <- round(ctsem::ctCollapse(popcov_est, 1, sd), digits)
          popcov_T[[j]] <- round(popcov_mean[[j]]/popcov_sd[[j]], digits)
          popcov_2.5[[j]] <- round(ctsem::ctCollapse(popcov_est, 1, stats::quantile, probs=.025), digits)
          popcov_97.5[[j]] <- round(ctsem::ctCollapse(popcov_est, 1, stats::quantile, probs=.975), digits)
          #
          popcor_mean[[j]] <- round(ctsem::ctCollapse(popcor_est, 1, mean), digits)
          popcor_sd[[j]] <- round(ctsem::ctCollapse(popcor_est, 1, sd), digits)
          popcor_T[[j]] <- round(popcor_mean[[j]]/popcor_sd[[j]], digits)
          popcor_2.5[[j]] <- round(ctsem::ctCollapse(popcor_est, 1, stats::quantile, probs=.025), digits)
          popcor_97.5[[j]] <- round(ctsem::ctCollapse(popcor_est, 1, stats::quantile, probs=.975), digits)
        }
        # do the same for the average covariance ()
        counter <- 0
        T0varMean <- matrix(NA, nrow=(n.latent*(n.latent+1)+n.latent^2), ncol=length(paramTmp)); dim(T0varMean)
        for (l in 1:(n.latent^2)) {
          for (m in 1:l) {
            counter <- counter + 1
            paramTmp <- rawT0varTmp[,counter]; param
            for (p in 1:length(paramTmp)) {
              param <- paramTmp[p]; param
              T0varMean[counter, p] <- eval(parse(text=transforms[counter]))
            }
          }
        }
        # make matrices out of vector
        for (k in 1:(length(T0COVCoeff))) {
          tmpMatMean <- tmpMatSD <- matrix(NA, n.latent^2, n.latent^2)
          counter <- 0
          for (l in 1:(n.latent^2)) {
            for (m in 1:l) {
              counter <- counter + 1
              tmpMatMean[l,m] <- mean(T0COVCoeff[[k]][, counter])
              tmpMatSD[l,m] <- sd(T0COVCoeff[[k]][, counter])
            }
          }
          tmpMatMean[upper.tri(tmpMatMean, diag = T)] <- t(tmpMatMean)[upper.tri(tmpMatMean, diag=T)]
          tmpMatSD[upper.tri(tmpMatSD, diag = T)] <- t(tmpMatSD)[upper.tri(tmpMatSD, diag=T)]
          T0COVCoeffMean[[k]] <- tmpMatMean
          T0COVCoeffSD[[k]] <- tmpMatSD
        }
        # do the same for the average matrix
        counter <- 0
        tmpMatMean <- tmpMatSD <- matrix(NA, n.latent^2, n.latent^2)
        for (l in 1:(n.latent^2)) {
          for (m in 1:l) {
            counter <- counter + 1
            tmpMatMean[l,m] <- mean(T0varMean[counter, ])
            tmpMatSD[l,m] <- mean(T0varMean[counter, ])
          }
        }
        tmpMatMean[upper.tri(tmpMatMean, diag = T)] <- t(tmpMatMean)[upper.tri(tmpMatMean, diag=T)]
        tmpMatSD[upper.tri(tmpMatSD, diag = T)] <- t(tmpMatSD)[upper.tri(tmpMatSD, diag=T)]
        T0varMeanMean <- tmpMatMean
        T0varMeanSD <- tmpMatSD
        T0COVCoeffMean
        model_popsd <- "Random intercepts were estimated per primary study. It is recommended to compare them with ctmaInit results to ensure the present results are accurate. Currently, thes should correspond to the rawpopcovbase-slot in the ctsem fit object"
        model_popcov_m <- T0COVCoeffMean
        model_popcov_sd <- T0COVCoeffSD
      }
    }
  }
  return(list(popcov_mean=popcov_mean, popcov_sd=popcov_sd, popcov_T=popcov_T,
              popcov_2.5=popcov_2.5, popcov_97.5=popcov_97.5,
              popcor_mean=popcor_mean, popcor_sd=popcor_sd, popcor_T=popcor_T,
              popcor_2.5=popcor_2.5, popcor_97.5=popcor_97.5,
              message=message))
}
