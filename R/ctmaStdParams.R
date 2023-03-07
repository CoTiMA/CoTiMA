#' ctmaStdParams
#'
#' @description Computes standardized drift effects from a CoTiMA or ctsem fit object. Can only handle CLPM or RI-CLPM fit objects.
#'
#' @param fit CoTiMA or ctsem fit object with or without random intercepts
#' @param times scalar (1 by defualt) or vector of scalars defining the discrete time lags for which standardized drift effects are computed.
#' @param digits rounding (4 by default)
#' @param standardize logical. TRUE (default) or FALSE (does not standardize and just computes discrete time effects)

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
ctmaStdParams <- function(fit=NULL, times=1, digits=4, standardize=TRUE) {

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
  if (length(tmp1) > n.latent) riT0 <- tmp1[1:n.latent] else riT0 <- c()
  if (length(tmp1) > n.latent) riTt <- tmp1[(n.latent+1):(2*n.latent) ] else riTt <- c()
  riT0; riTt

  # drift
  fitsum <- summary(fit)
  drift <- fitsum$parmatrices[, c("matrix", "Mean")]; drift
  drift <- matrix(drift[drift$matrix=="DRIFT", 2], n.latent, n.latent, byrow=T); drift

  drift_s <- fit$stanfit$transformedpars$pop_DRIFT; drift_s[1,,]

  # Q
  Q <- fitsum$parmatrices[, c("matrix", "Mean")]; Q
  Q <- matrix(Q[Q$matrix=="DIFFUSIONcov", 2], n.latent, n.latent, byrow=T); Q

  #Q_s <- fit$stanfit$transformedpars$pop_DIFFUSION; Q_s[1,,]
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

  # T0cov
  T0cov <- fitsum$parmatrices[, c("matrix", "Mean")]; T0cov
  T0cov <- matrix(T0cov[T0cov$matrix=="T0cov", 2], n.latent, n.latent, byrow=T); T0cov

  T0cov_s <- fit$stanfit$transformedpars$pop_T0cov[, riT0, riT0] #; T0cov_s[1,,]
  #if (dim(T0cov_s)[2] == 0) T0cov_s <- array(0, dim=c(dim(drift_s)[1], n.latent, n.latent))
  if (dim(T0cov_s)[2] == 0) T0cov_s <- fit$stanfit$transformedpars$pop_T0cov[, 1:n.latent, 1:n.latent] #; T0cov_s[1,,]

  # partial (co-)variance of random intercepts
  # https://stats.stackexchange.com/questions/557855/partial-covariance-matrix-after-linear-transformations
  Cov_s <- fit$stanfit$transformedpars$pop_T0cov[, c(riT0, riTt), c(riT0, riTt)]; Cov_s[1,,]
  #if (dim(Cov_s)[2] == 0) Cov_s <- array(0, dim=c(dim(drift_s)[1], n.latent, n.latent))
  pRIcov_s <- list() # array(data=NA, dim=c(dim(drift_s)[1], n.latent^2, n.latent^2))
  IVs <- 1:n.latent; IVs
  DVs <- (n.latent+1):(2*n.latent); DVs
  for (i in 1:dim(drift_s)[1]) {
    if ((dim(Cov_s)[2] != 0) & (length(riTt) != 0)) {
      SX <- Cov_s[i,IVs,IVs]; SX
      SXY <- Cov_s[i, IVs, DVs]; SXY
      SY <- Cov_s[i, DVs, DVs]; SY
      pRIcov_s[[i]] <- SX - SXY %*% solve(SY) %*% t(SXY); pRIcov_s[[i]]
    } else {
      pRIcov_s[[i]] <- matrix(0, n.latent, n.latent)
    }
  }
  #pRIcov_s

  # RI covs
  e <- ctsem::ctExtract(fit)
  dim(e$popcov)
  if (!(is.null(e$popcov))) { # e$popcov == Cov_s above
    RIcov <- ctsem::ctCollapse(e$popcov, 1, mean); RIcov
    # partial (co-)variance of random intercepts
    # https://stats.stackexchange.com/questions/557855/partial-covariance-matrix-after-linear-transformations
    SX <- RIcov[c(IVs), c(IVs)]; SX
    SXY <- RIcov[c(IVs), c(DVs)]; SXY
    SY <- RIcov[c(DVs), c(DVs)]; SY
    pRIcov <- SX - SXY %*% solve(SY) %*% t(SXY); pRIcov
    RIcov <- RIcov[riTt, riTt]; RIcov
    if (dim(RIcov)[1] == 0) {
      RIcov <- matrix(0, n.latent, n.latent); RIcov
      pRIcov <- matrix(0, n.latent, n.latent); pRIcov
    }
    RIcov_s <- fit$stanfit$transformedpars$pop_T0cov[, riTt, riTt]; RIcov_s[1,,]
  } else {
    RIcov <- matrix(0, n.latent, n.latent); RIcov
    pRIcov <- matrix(0, n.latent, n.latent); pRIcov
    RIcov_s <- array(data=0, dim=c(dim(drift_s)[1], n.latent, n.latent))
  }

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
    #Ttvar <- Ttvar + Psi + RIcov; Ttvar
    Ttvar <- Ttvar + Psi + pRIcov; Ttvar
    standardized_effects[[paste0("time = ", time)]] <- matrix(NA, n.latent, n.latent)
    for (dv in 1:n.latent) {
      for (iv in 1:n.latent) {
        if (standardize == TRUE) {
          standardized_effects[[paste0("time = ", time)]][dv, iv] <- round(OpenMx::expm(drift * time)[dv,iv] * T0cov[iv,iv]^.5 / Ttvar[dv,dv]^.5, digits)
        } else {
          standardized_effects[[paste0("time = ", time)]][dv, iv] <- round(OpenMx::expm(drift * time)[dv,iv], digits)
        }
      }
    }
    dimnames(standardized_effects[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    standardized_effects_table[[paste0("time = ", time)]] <- matrix(c(standardized_effects[[paste0("time = ", time)]]), ncol=1)
    rownames(standardized_effects_table[[paste0("time = ", time)]]) <- driftEffects
    #standardized_effects_table[[paste0("time = ", time)]]
  }
  standardized_effects_popmeans <- standardized_effects
  standardized_effects_table_popmeans <- standardized_effects_table
  #standardized_effects_popmeans
  #standardized_effects_table_popmeans

  # use all samples
  standardized_effects_s <- list()
  #standardized_effects_table_s <- list()
  for (time in times) {
    for (i in 1:dim(drift_s)[1]) {
      if (time == times[1]) {
        standardized_effects_s[[i]] <- list()
        #standardized_effects_table_s[[i]] <- list()
      }
      Ttvar <- OpenMx::expm(drift_s[i,,] * time) %*% T0cov_s[i,,] %*% t(OpenMx::expm(drift_s[i,,] * time)); Ttvar
      Psi <- drift_hatch_solve_s[i,,] %*% (OpenMx::expm(drift_hatch_s[i,,] %x% time) - diag(1, n.latent^2, n.latent^2) )  %*% c(Q_s[i,,]); Psi
      Psi <- matrix(Psi, n.latent, n.latent); Psi
      #Ttvar <- Ttvar + Psi + RIcov_s[i,,]; Ttvar
      Ttvar <- Ttvar + Psi + pRIcov_s[[i]]; Ttvar
      standardized_effects_s[[i]][[paste0("time = ", time)]] <- matrix(NA, n.latent, n.latent)
      #standardized_effects_table_s[[i]][[paste0("time = ", time)]] <- matrix(NA, nrow=n.latent^2, ncol=1)
      for (dv in 1:n.latent) {
        for (iv in 1:n.latent) {
          if (standardize == TRUE) {
            standardized_effects_s[[i]][[paste0("time = ", time)]][dv, iv] <- round(OpenMx::expm(drift_s[i,,] * time)[dv,iv] * T0cov_s[i,,][iv,iv]^.5 / Ttvar[dv,dv]^.5, digits)
          } else {
            standardized_effects_s[[i]][[paste0("time = ", time)]][dv, iv] <- round(OpenMx::expm(drift_s[i,,] * time)[dv,iv], digits)
          }
        }
      }
      dimnames(standardized_effects_s[[i]][[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
      #standardized_effects_table_s[[i]][[paste0("time = ", time)]] <- matrix(c(standardized_effects_s[[i]][[paste0("time = ", time)]]), ncol=1)
      #rownames(standardized_effects_table_s[[i]][[paste0("time = ", time)]]) <- driftEffects

    } # end i
  }

  standardized_effects_table_s <- list()
  #standardized_effects_s_backup <- standardized_effects_s
  counter <- 0
  means <- sds <- ul95 <- ll95 <- ul99 <- ll99 <- medians <- matrix(NA, nrow=length(times), ncol=(n.latent^2))
  #means_table <- sds_table <- ul95_table <- ll95_table <- ul99_table <- ll99_table <- medians_table <- matrix(NA, nrow=length(times), ncol=(n.latent^2))
  standardized_effects_m <- standardized_effects_sd <- list()
  standardized_effects_ul95 <- standardized_effects_ll95 <- list()
  standardized_effects_ul99 <- standardized_effects_ll99 <- list()
  standardized_effects_medians <- list()
  #standardized_effects_table_m <- standardized_effects_table_sd <- list()
  #standardized_effects_table_ul95 <- standardized_effects_table_ll95 <- list()
  #standardized_effects_table_ul99 <- standardized_effects_table_ll99 <- list()
  #standardized_effects_table_medians <- list()
  for (time in times) {
    #time <- times[1]
    counter <- counter + 1
    tmp <- lapply(standardized_effects_s, function(x) x[[paste0("time = ", time)]])
    #tmp2 <- lapply(standardized_effects_table_s, function(x) x[[paste0("time = ", time)]])
    for (i in 1:(n.latent^2)) {
      means[counter, i] <- mean(unlist(lapply(tmp, function(x) x[i])))
      sds[counter, i] <- sd(unlist(lapply(tmp, function(x) x[i])))
      tmp2 <- unlist(lapply(tmp, function(x) x[i]))
      ul95[counter, i] <- stats::quantile(tmp2, .95)
      ll95[counter, i] <- stats::quantile(tmp2, .05)
      ul99[counter, i] <- stats::quantile(tmp2, .99)
      ll99[counter, i] <- stats::quantile(tmp2, .01)
      medians[counter, i] <- stats::quantile(tmp2, .50)
    }
    #standardized_effects_m[[paste0("time = ", time)]] <- round(matrix(means[counter,], n.latent, n.latent), digits)
    #standardized_effects_sd[[paste0("time = ", time)]] <- round(matrix(sds[counter,], n.latent, n.latent), digits)
    #standardized_effects_ul95[[paste0("time = ", time)]] <- round(matrix(ul95[counter,], n.latent, n.latent), digits)
    #standardized_effects_ll95[[paste0("time = ", time)]] <- round(matrix(ll95[counter,], n.latent, n.latent), digits)
    #standardized_effects_ul99[[paste0("time = ", time)]] <- round(matrix(ul99[counter,], n.latent, n.latent), digits)
    #standardized_effects_ll99[[paste0("time = ", time)]] <- round(matrix(ll99[counter,], n.latent, n.latent), digits)
    #standardized_effects_medians[[paste0("time = ", time)]] <- round(matrix(medians[counter,], n.latent, n.latent), digits)
    #dimnames(standardized_effects_m[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #dimnames(standardized_effects_sd[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #dimnames(standardized_effects_ul95[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #dimnames(standardized_effects_ll95[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #dimnames(standardized_effects_ul99[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #dimnames(standardized_effects_ll99[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #dimnames(standardized_effects_medians[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
    #
    standardized_effects_m[[paste0("time = ", time)]] <- round(means[counter,], digits)
    standardized_effects_sd[[paste0("time = ", time)]] <- round(sds[counter,], digits)
    standardized_effects_ul95[[paste0("time = ", time)]] <- round(ul95[counter,], digits)
    standardized_effects_ll95[[paste0("time = ", time)]] <- round(ll95[counter,], digits)
    standardized_effects_ul99[[paste0("time = ", time)]] <- round(ul99[counter,], digits)
    standardized_effects_ll99[[paste0("time = ", time)]] <- round(ll99[counter,], digits)
    standardized_effects_medians[[paste0("time = ", time)]] <- round(medians[counter,], digits)
    tmp <- cbind(standardized_effects_m[[paste0("time = ", time)]], standardized_effects_sd[[paste0("time = ", time)]],
          standardized_effects_ll99[[paste0("time = ", time)]], standardized_effects_ll95[[paste0("time = ", time)]],
          standardized_effects_medians[[paste0("time = ", time)]],
          standardized_effects_ul95[[paste0("time = ", time)]], standardized_effects_ul99[[paste0("time = ", time)]])
    rownames(tmp) <- driftEffects
    colnames(tmp) <- c("mean", "sd", "LL99", "LL95", "median", "UL95", "UL99"); tmp
    standardized_effects_table_s[[paste0("time = ", time)]] <- tmp
    }
  #standardized_effects_sampled <- list(standardized_effects_mean=standardized_effects_m,
  #                                     standardized_effects_median=standardized_effects_medians,
  #                                     standardized_effects_sd=standardized_effects_sd,
  #                                     standardized_effects_ll95=standardized_effects_ll95,
  #                                     standardized_effects_ul95=standardized_effects_ul95,
  #                                     standardized_effects_ll99=standardized_effects_ll99,
  #                                     standardized_effects_ul99=standardized_effects_ul99
  #                                     )
  #standardized_effects_sampled
  return(list(popmeans=standardized_effects_table_popmeans, sampled=standardized_effects_table_s))
}

