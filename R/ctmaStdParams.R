#' ctmaStdParams
#'
#' @description Computes standardized drift effects from a CoTiMA or ctsem fit object. Can only handle CLPM or RI-CLPM fit objects.
#'
#' @param fit CoTiMA or ctsem fit object with or without random intercepts
#' @param times scalar (1 by defualt) or vector of scalars defining the discrete time lags for which standardized drift effects are computed.
#' @param digits rounding (4 by default)
#' @param standardize logical. TRUE (default) or FALSE (does not standardize and just computes discrete time effects)
#' @param oneTailed logical. FALSE (default) or TRUE. If TRUE, one-tailed CIs will be reported

#'
#' @importFrom  ctsem ctExtract ctCollapse
#' @importFrom  OpenMx expm
#' @importFrom  stats quantile sd
#'
#' @export ctmaStdParams
#'
#' @examples
#' \dontrun{
#' ctmaStdParams(CoTiMAFullFit_3_orig, times=c(.1, 1, 2), digits=6, standardize=TRUE)
#' }
#'
#'
#' @return ctmaStdParams returns a list of standardized discrete time drift matrices for different time intervals.
#'
ctmaStdParams <- function(fit=NULL, times=1, digits=4, standardize=TRUE, oneTailed=FALSE) {

  Msg <- "Function ctmaStdParams is under development. There is no guarantee effects will be correct."
  message(Msg)

  if (is.null(fit)) {
    ErrorMsg <- "\nA fitted ctsem or CoTiMA object has to be supplied to compute standardized effects. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if ((is(fit) != "ctStanFit") & (is(fit) != "CoTiMAFit"))  {
    ErrorMsg <- "\nThe fit object is neither of class ctsem nor of class CoTiMA. Please provide ctsem or CoTiMA fit object. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (length(times) == 1 & times[1] == 1 ) {
    Msg <- "Standardized effects will only be computed for time interval = 1. Provide more time intervals if needed."
    message(Msg)
  }

  if (standardize == FALSE) {
    Msg <- cat("Standardized effects will NOT be computed. Just discrete time effects across intervals = ",  times)
    message(Msg)
  }

  if (is(fit) == "CoTiMAFit") {
    tmp <- fit$studyFitList$ctstanmodelbase$pars
    if (is.null(tmp)) {
      ErrorMsg <- "\nFitted ctsem or CoTiMA object has been modified. Cannot find specification of params. \nGood luck for the next try!"
      stop(ErrorMsg)
    }
    fit <- fit$studyFitList
    fit$ctstanmodel$pars <- tmp
  }

  n.latent <- fit$ctstanmodel$n.latent; n.latent

  indvarying <- which(fit$ctstanmodelbase$pars$indvarying); indvarying
  tmp1 <- grep("T0MEANS", fit$ctstanmodelbase$pars$matrix[indvarying]); tmp1
  tmp1 <- c(tmp1, grep("CINT", fit$ctstanmodelbase$pars$matrix[indvarying])); tmp1
  tmp1 <- c(tmp1, grep("MANIFESTMEANS", fit$ctstanmodelbase$pars$matrix[indvarying])); tmp1
  if (length(tmp1) > n.latent) riT0 <- tmp1[1:n.latent] else riT0 <- c()
  if (length(tmp1) > n.latent) riTt <- tmp1[(n.latent+1):(2*n.latent) ] else riTt <- c()
  riT0; riTt
  # distingusih RI cint from RI manifestmeans
  if (length(grep("MANIFESTMEANS", fit$ctstanmodelbase$pars$matrix[indvarying])) > 0) power <- 2 else power <- 1
  power

  # report one- vs. twotailed
  if (oneTailed == FALSE) testText <- "two-tailed" else testText <- "one-tailed"

  # drift
  fitsum <- summary(fit)
  drift <- fitsum$parmatrices[, c("matrix", "Mean")]; drift
  drift <- matrix(drift[drift$matrix=="DRIFT", 2], n.latent, n.latent, byrow=T); drift

  drift_s <- fit$stanfit$transformedpars$pop_DRIFT; drift_s[1,,]

  # Q
  Q <- fitsum$parmatrices[, c("matrix", "Mean")]; Q
  Q <- matrix(Q[Q$matrix=="DIFFUSIONcov", 2], n.latent, n.latent, byrow=T); Q
  #t(chol(Q))

  Q_s <- fit$stanfit$transformedpars$pop_DIFFUSIONcov; Q_s[1,,]

  # drift hatch
  drift_hatch <- drift %x% diag(1, ncol=n.latent, nrow=n.latent) +
    diag(1, ncol=n.latent, nrow=n.latent) %x% drift; drift_hatch
  drift_hatch_solve <- solve(drift_hatch); drift_hatch_solve

  drift_hatch_s <- drift_hatch_solve_s <- array(data=NA, dim=c(dim(drift_s)[1], n.latent^2, n.latent^2))
  for (i in 1:dim(drift_s)[1]) {
    drift_hatch_s[i,,] <- drift_s[i,,] %x% diag(1, ncol=n.latent, nrow=n.latent) +
      diag(1, ncol=n.latent, nrow=n.latent) %x% drift_s[i,,]; drift_hatch_s[i,,]
    drift_hatch_solve_s[i,,] <- solve(drift_hatch_s[i,,]); drift_hatch_solve_s[i,,]
  }

  # useful for several estimates
  e <- ctsem::ctExtract(fit)

  # T0cov
  T0cov <- ctsem::ctCollapse(e$pop_T0cov, 1, mean)[1:n.latent,1:n.latent]; T0cov; t(chol(T0cov)/10)
  T0cov_s <- e$pop_T0cov[, 1:n.latent,1:n.latent]; T0cov_s[1,,]

  IVs <- 1:n.latent; IVs
  DVs <- (n.latent+1):(2*n.latent); DVs

  # RI covs
  if (!(is.null(e$popcov))) {
    #T0cov
    RIcov <- ctsem::ctCollapse(e$popcov, 1, mean); RIcov
    #RIcov[1:2,3:4] <- RIcov[1:2,3:4]/10
    #RIcov[3:4,] <- RIcov[3:4, ]/10
    #RIcov
    # partial (co-)variance of random intercepts
    # https://stats.stackexchange.com/questions/557855/partial-covariance-matrix-after-linear-transformations
    # https://link.springer.com/content/pdf/10.1007/BF02294097.pdf (Eq. 4)
    SXX <- RIcov[IVs, IVs]; SXX
    SXY <- RIcov[IVs, DVs]; SXY
    SYY <- RIcov[DVs, DVs]; SYY
    pRIcov <- SXX - SXY %*% solve(SYY) %*% t(SXY); pRIcov
    # try
    SYX <- RIcov[DVs, IVs]; SYX
    pRIcov <- SYY - SYX %*% solve(SXX) %*% t(SXY); pRIcov
    #if (power == 2) pRIcov <- pRIcov^2
    # https://psyarxiv.com/u2zdv/ p. 4
    #SXY <- RIcov[DVs, DVs]; SX
    #SXZ <- RIcov[DVs, IVs]; SXZ
    #SZZ <- RIcov[IVs, IVs]; SZZ
    #SZY <- RIcov[IVs, DVs]; SZY
    #pRIcov <- SXY - SXZ %*% solve(SZZ) %*% SZY; pRIcov_s[[i]]
  } else {
    pRIcov <- matrix(0, n.latent, n.latent); pRIcov
  }
  #RIcov
  #pRIcov

  # partial (co-)variance of random intercepts
  # https://stats.stackexchange.com/questions/557855/partial-covariance-matrix-after-linear-transformations
  pRIcov_s <- list() # array(data=NA, dim=c(dim(drift_s)[1], n.latent^2, n.latent^2))
  for (i in 1:dim(drift_s)[1]) {
    if (!(is.null(e$popcov))) {
      #RIcov_s <- e$popcov[,(n.latent+1):(2*n.latent), (n.latent+1):(2*n.latent)]#; RIcov_s[1,,]
      RIcov_s <- e$popcov#; RIcov_s[1,,]
      SXX <- RIcov_s[i, IVs, IVs]; SXX
      SXY <- RIcov_s[i, IVs, DVs]; SXY
      SYY <- RIcov_s[i,  DVs, DVs]; SYY
      #SXX - SXY %*% solve(SYY) %*% t(SXY)
      pRIcov_s[[i]] <- SXX - SXY %*% solve(SYY) %*% t(SXY); pRIcov_s[[i]]
      #if (power == 2) pRIcov_s[[i]] <- pRIcov_s[[i]]^2
      SYX <- RIcov_s[i, DVs, IVs]; SYX
      pRIcov_s[[i]] <- SYY - SYX %*% solve(SXX) %*% t(SXY); pRIcov_s[[i]]
      # https://psyarxiv.com/u2zdv/ p. 4
      #SXY <- RIcov_s[i, DVs, DVs]; SX
      #SXZ <- RIcov_s[i, DVs, IVs]; SXZ
      #SZZ <- RIcov_s[i, IVs, IVs]; SZZ
      #SZY <- RIcov_s[i, IVs, DVs]; SZY
      #pRIcov_s[[i]] <- SXY - SXZ %*% solve(SZZ) %*% SZY; pRIcov_s[[i]]
    } else {
    pRIcov_s[[i]] <- matrix(0, n.latent, n.latent)
    }
  }
  #RIcov_s[1,,]
  #pRIcov_s[[1]]

  # get drift names
  driftNames <- fit$ctstanmodelbase$latentNames; driftNames
  driftEffects <- fit$ctstanmodelbase$pars[fit$ctstanmodelbase$pars$matrix == "DRIFT", ]$param; driftEffects
  driftEffects <- c(t(matrix(driftEffects, n.latent, n.latent))); driftEffects

  # use popmeans
  standardized_effects <- list()
  standardized_effects_table <- list()
  for (time in times) {
    #time <- times[1]; time
    Ttvar <- OpenMx::expm(drift * time) %*% T0cov %*% t(OpenMx::expm(drift * time)); Ttvar
    Psi <- drift_hatch_solve %*% (OpenMx::expm(drift_hatch %x% time) - diag(1, n.latent^2, n.latent^2) )  %*% c(Q); Psi
    Psi <- matrix(Psi, n.latent, n.latent); Psi
    #T0cov
    #Ttvar + Psi
    #pRIcov
    #Ttvar + Psi + pRIcov
    #RIcov
    #fitsum$popsd$mean^2
    #Ttvar <- Ttvar + Psi + RIcov[riTt,riTt]; Ttvar
    if (power == 2) pRIcov <- 0
    Ttvar <- Ttvar + Psi + pRIcov; Ttvar
    standardized_effects[[paste0("time = ", time, " ", testText)]] <- matrix(NA, n.latent, n.latent)
    for (dv in 1:n.latent) {
      for (iv in 1:n.latent) {
        if (standardize == TRUE) {
          standardized_effects[[paste0("time = ", time, " ", testText)]][dv, iv] <- round(OpenMx::expm(drift * time)[dv,iv] * T0cov[iv,iv]^.5 / Ttvar[dv,dv]^.5, digits)
        } else {
          standardized_effects[[paste0("time = ", time, " ", testText)]][dv, iv] <- round(OpenMx::expm(drift * time)[dv,iv], digits)
        }
      }
    }
    dimnames(standardized_effects[[paste0("time = ", time, " ", testText)]]) <- list(c(driftNames), c(driftNames))
    standardized_effects_table[[paste0("time = ", time, " ", testText)]] <- matrix(c(standardized_effects[[paste0("time = ", time, " ", testText)]]), ncol=1)
    rownames(standardized_effects_table[[paste0("time = ", time, " ", testText)]]) <- driftEffects
  }
  standardized_effects_popmeans <- standardized_effects
  standardized_effects_table_popmeans <- standardized_effects_table
  #standardized_effects_popmeans
  #standardized_effects_table_popmeans

  # use all samples
  standardized_effects_s <- list()
  for (time in times) {
    for (i in 1:dim(drift_s)[1]) {
      if (time == times[1]) {
        standardized_effects_s[[i]] <- list()
      }
      Ttvar <- OpenMx::expm(drift_s[i,,] * time) %*% T0cov_s[i,,] %*% t(OpenMx::expm(drift_s[i,,] * time)); Ttvar
      Psi <- drift_hatch_solve_s[i,,] %*% (OpenMx::expm(drift_hatch_s[i,,] %x% time) - diag(1, n.latent^2, n.latent^2) )  %*% c(Q_s[i,,]); Psi
      Psi <- matrix(Psi, n.latent, n.latent); Psi
      if (power == 2) pRIcov_s[[i]] <- 0
      Ttvar <- Ttvar + Psi + pRIcov_s[[i]]; Ttvar
      standardized_effects_s[[i]][[paste0("time = ", time, " ", testText)]] <- matrix(NA, n.latent, n.latent)
      for (dv in 1:n.latent) {
        for (iv in 1:n.latent) {
          if (standardize == TRUE) {
            standardized_effects_s[[i]][[paste0("time = ", time, " ", testText)]][dv, iv] <- round(OpenMx::expm(drift_s[i,,] * time)[dv,iv] * T0cov_s[i,,][iv,iv]^.5 / Ttvar[dv,dv]^.5, digits)
          } else {
            standardized_effects_s[[i]][[paste0("time = ", time, " ", testText)]][dv, iv] <- round(OpenMx::expm(drift_s[i,,] * time)[dv,iv], digits)
          }
        }
      }
      dimnames(standardized_effects_s[[i]][[paste0("time = ", time, " ", testText)]]) <- list(c(driftNames), c(driftNames))
    } # end i
  }

  standardized_effects_table_s <- list()
  counter <- 0
  means <- sds <- ul95 <- ll95 <- ul99 <- ll99 <- medians <- matrix(NA, nrow=length(times), ncol=(n.latent^2))
  standardized_effects_m <- standardized_effects_sd <- list()
  standardized_effects_ul95 <- standardized_effects_ll95 <- list()
  standardized_effects_ul99 <- standardized_effects_ll99 <- list()
  standardized_effects_medians <- list()
  for (time in times) {
    #time <- times[1]
    counter <- counter + 1
    tmp <- lapply(standardized_effects_s, function(x) x[[paste0("time = ", time, " ", testText)]])
    for (i in 1:(n.latent^2)) {
      means[counter, i] <- mean(unlist(lapply(tmp, function(x) x[i])))
      sds[counter, i] <- sd(unlist(lapply(tmp, function(x) x[i])))
      tmp2 <- unlist(lapply(tmp, function(x) x[i]))
      if (oneTailed == FALSE) {
        ul95[counter, i] <- stats::quantile(tmp2, .975)
        ll95[counter, i] <- stats::quantile(tmp2, .025)
        ul99[counter, i] <- stats::quantile(tmp2, .995)
        ll99[counter, i] <- stats::quantile(tmp2, .005)
      } else {
        ul95[counter, i] <- stats::quantile(tmp2, .95)
        ll95[counter, i] <- stats::quantile(tmp2, .05)
        ul99[counter, i] <- stats::quantile(tmp2, .99)
        ll99[counter, i] <- stats::quantile(tmp2, .01)
      }
      medians[counter, i] <- stats::quantile(tmp2, .50)
    }
    standardized_effects_m[[paste0("time = ", time, " ", testText)]] <- round(means[counter,], digits)
    standardized_effects_sd[[paste0("time = ", time, " ", testText)]] <- round(sds[counter,], digits)
    standardized_effects_ul95[[paste0("time = ", time, " ", testText)]] <- round(ul95[counter,], digits)
    standardized_effects_ll95[[paste0("time = ", time, " ", testText)]] <- round(ll95[counter,], digits)
    standardized_effects_ul99[[paste0("time = ", time, " ", testText)]] <- round(ul99[counter,], digits)
    standardized_effects_ll99[[paste0("time = ", time, " ", testText)]] <- round(ll99[counter,], digits)
    standardized_effects_medians[[paste0("time = ", time, " ", testText)]] <- round(medians[counter,], digits)
    tmp <- cbind(standardized_effects_m[[paste0("time = ", time, " ", testText)]], standardized_effects_sd[[paste0("time = ", time, " ", testText)]],
          standardized_effects_ll99[[paste0("time = ", time, " ", testText)]], standardized_effects_ll95[[paste0("time = ", time, " ", testText)]],
          standardized_effects_medians[[paste0("time = ", time, " ", testText)]],
          standardized_effects_ul95[[paste0("time = ", time, " ", testText)]], standardized_effects_ul99[[paste0("time = ", time, " ", testText)]])
    rownames(tmp) <- driftEffects
    colnames(tmp) <- c("mean", "sd", "LL99", "LL95", "median", "UL95", "UL99"); tmp
    standardized_effects_table_s[[paste0("time = ", time, " ", testText)]] <- tmp
    }

  return(list(popmeans=standardized_effects_table_popmeans, sampled=standardized_effects_table_s))
}

