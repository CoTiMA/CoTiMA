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

  tmp1 <- which(fit$ctstanmodel$pars$indvarying); tmp1
  riT0 <- grep("T0m", fit$ctstanmodel$pars[tmp1,"param"]); riT0
  riTt <- grep("cin", fit$ctstanmodel$pars[tmp1,"param"]); riTt
  riTt <- c(riTt, grep("mm_", fit$ctstanmodel$pars[tmp1,"param"])); riTt

  # drift
  fitsum <- summary(fit)
  drift <- fitsum$parmatrices[, c("matrix", "Mean")]; drift
  drift <- matrix(drift[drift$matrix=="DRIFT", 2], n.latent, n.latent, byrow=T); drift
  # Q
  Q <- fitsum$parmatrices[, c("matrix", "Mean")]; Q
  Q <- matrix(Q[Q$matrix=="DIFFUSIONcov", 2], n.latent, n.latent, byrow=T); Q
  # drift hatch
  drift_hatch <- drift %x% diag(1, ncol=n.latent, nrow=n.latent) +
    diag(1, ncol=n.latent, nrow=n.latent) %x% drift; drift_hatch
  drift_hatch_solve <- solve(drift_hatch); drift_hatch_solve
  # T0cov
  T0cov <- fitsum$parmatrices[, c("matrix", "Mean")]; T0cov
  T0cov <- matrix(T0cov[T0cov$matrix=="T0cov", 2], n.latent, n.latent, byrow=T); T0cov
  # RI covs
  e <- ctsem::ctExtract(fit)
  if (!(is.null(e$popcov))) {
    RIcov <- ctsem::ctCollapse(e$popcov, 1, mean); RIcov
    RIcov <- RIcov[riTt, riTt]; RIcov
    if (dim(RIcov)[1] == 0) RIcov <- matrix(0, n.latent, n.latent); RIcov
  } else {
    RIcov <- matrix(0, n.latent, n.latent); RIcov
  }

  # get drift names
  driftNames <- fit$ctstanmodelbase$latentNames; driftNames

  standardized_effects <- list()
  Ttvar <- list()
  for (time in times) {
    #time <- times[1]; time
    Ttvar[[paste0("time = ",time)]] <- OpenMx::expm(drift * time) %*% T0cov %*% t(OpenMx::expm(drift * time)); Ttvar[[paste0("time = ",time)]]
    OpenMx::expm(drift_hatch %x% time)
    #Psi <- drift_hatch_solve %*% (OpenMx::expm(drift_hatch %x% time) - diag(1, 2*n.latent, 2*n.latent) )  %*% c(Q); Psi
    Psi <- drift_hatch_solve %*% (OpenMx::expm(drift_hatch %x% time) - diag(1, n.latent^2, n.latent^2) )  %*% c(Q); Psi
    Psi <- matrix(Psi, n.latent, n.latent); Psi
    Ttvar[[paste0("time = ", time)]] <- Ttvar[[paste0("time = ", time)]] + Psi + RIcov; Ttvar[[paste0("time = ", time)]]
    standardized_effects[[paste0("time = ", time)]] <- matrix(NA, n.latent, n.latent)
    for (dv in 1:n.latent) {
      for (iv in 1:n.latent) {
        if (standardize == TRUE) {
          standardized_effects[[paste0("time = ", time)]][dv, iv] <- round(OpenMx::expm(drift * time)[dv,iv] * T0cov[iv,iv]^.5 / Ttvar[[paste0("time = ", time)]][dv,dv]^.5, digits)
        } else {
          standardized_effects[[paste0("time = ", time)]][dv, iv] <- round(OpenMx::expm(drift * time)[dv,iv], digits)
        }
      }
    }
    dimnames(standardized_effects[[paste0("time = ", time)]]) <- list(c(driftNames), c(driftNames))
  }

  return(standardized_effects)
}

