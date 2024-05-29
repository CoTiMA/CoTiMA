#' ctmaMMtoCINT
#'
#' @description Compute covariance of CINT-based random intercepts obtained with a MANIFESTMEANS model specification
#'
#' @param ctmaFitObject fit object created with ctmaInit or ctmaFitObject
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
ctmaMMtoCINT <- function(ctmaFitObject=NULL) {
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
    arguments <- ctmaFitObject$
    n.latent <- arguments$n.latent; n.latent
    n.manifest <- arguments$n.manifest; n.manifest
    digits <- arguments$digits; digits
  }
  if (class(ctmaFitObject) == "ctStanFit") {
    arguments <- ctmaFitObject$ctstanmodelbase
    arguments$scaleTime <- 1
    n.latent <- arguments$n.latent; n.latent
    n.manifest <- arguments$n.manifest; n.manifest
    digits <- 4; digits
    #
    pars <- arguments$pars
    if (all(pars[pars$matrix=="MANIFESTMEANS", "indvarying"] == TRUE)) arguments$indVarying <- TRUE
    if (all(pars[pars$matrix=="CINT", "indvarying"] == TRUE)) arguments$indVarying <- "CINT"
    if (is.null(arguments$indVarying)) {
      ErrorMsg <- "The fit object provided used neither CINT nor MANIFESTMEANS to model random intercepts."
      stop(ErrorMsg)
    }
  }

  popcov_mean <- popcov_sd <- popcov_2.5 <- popcov_97.5 <- popcov_T <- list()
  popcor_mean <- popcor_sd <- popcor_2.5 <- popcor_97.5 <- popcor_T <- list()
  if (!is.null(arguments$scaleTime)) scaleTime <- arguments$scaleTime else scaleTime <- 1
  if (arguments$indVarying == TRUE) mmRI <- TRUE else mmRI <- FALSE
  if (arguments$indVarying == "CINT") {
    ErrorMsg <- "The fit object provided used CINT rather the MANIFESTMEANS to model random intercepts."
    stop(ErrorMsg)
  }

  for (i in 1:n.studies) {
    if (class(ctmaFitObject$studyFitList[[i]]) == "ctStanFit") {
      fit <- ctmaFitObject$studyFitList[[i]]
    } else {
      fit <- ctmaFitObject$studyFitList
    }
    e <- ctsem::ctExtract(fit)
    e$pop_DRIFT <- e$pop_DRIFT * scaleTime

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
    if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
      #### IDEA: Transformation matrix describing popcov_MM into popciv_cint transformations and then popcov_cint <- trans %*% popcov_mm %*% t(trans) (see: #https://stats.stackexchange.com/questions/113700/covariance-of-a-random-vector-after-a-linear-transformation)
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
      popcov_mean[[i]] <- ctsem::ctCollapse(e$popcov, 1, mean)
      popcov_sd[[i]] <- ctsem::ctCollapse(e$popcov, 1, sd)
      popcov_T[[i]] <- popcov_mean[[i]]/popcov_sd[[i]]
      popcov_2.5[[i]] <- ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.025)
      popcov_97.5[[i]] <- ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.975)
      #
      popcor_mean[[i]] <- ctsem::ctCollapse(e$popcor, 1, mean)
      popcor_sd[[i]] <- ctsem::ctCollapse(e$popcor, 1, sd)
      popcor_T[[i]] <- popcor_mean[[i]]/popcor_sd[[i]]
      popcor_2.5[[i]] <- ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.025)
      popcor_97.5[[i]] <- ctsem::ctCollapse(e$popcor, 1, stats::quantile, probs=.975)
    }
  }
  return(list(popcov_mean=popcov_mean, popcov_sd=popcov_sd, popcov_T=popcov_T,
              popcov_2.5=popcov_2.5, popcov_97.5=popcov_97.5,
              popcor_mean=popcor_mean, popcor_sd=popcor_sd, popcor_T=popcor_T,
              popcor_2.5=popcor_2.5, popcor_97.5=popcor_97.5,
              message=message))
}
