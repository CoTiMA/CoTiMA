#' ctmaMMtoCINT
#'
#' @description Compute covariance of CINT-based random intercepts obtained with a MANIFESTMEANS model specification
#'
#' @param ctmaFit fit object created with ctmaInit or ctmaFit
#'
#' @importFrom ctsem stats
#'
#' @examples
#' \donttest{
#' RI_cov <- ctmaMMtoCINT(ctmaFit=CoTiMAFullFit_3)
#' print(RI_cov)
#' }
#'
#' @export ctmaMMtoCINT
#'
#' @return returns covariance of CINT-based random intercepts.
#'
MMRItoCINTRI <- function(ctmaFit=NULL) {
  ### TRANSFORM RI modeled as manifest means into cint-based estimates
  if (class(ctmaFit$studyFitList[[1]]) == "ctStanFit") {
    n.studies <- length(ctmaFit$studyFitList)
  } else {
    n.studies <- 1
  }
  #fit <- ctmaFit$studyFitList
  arguments <- ctmaFit$argumentList
  n.latent <- arguments$n.latent; n.latent
  n.manifest <- arguments$n.manifest; n.manifest
  digits <- arguments$digits; digits
  popcov_mean <- popcov_sd <- popcov_2.5 <- popcov_97.5 <- popcov_T <- list()
  if (!is.null(arguments$scaleTime)) scaleTime <- arguments$scaleTime else scaleTime <- 1
  if (arguments$indVarying == TRUE) mmRI <- TRUE else mmRI <- FALSE

  for (i in 1:n.studies) {
    #i <- 1
    if (class(ctmaFit$studyFitList[[i]]) == "ctStanFit") {
      fit <- ctmaFit$studyFitList[[i]]
    } else {
      fit <- ctmaFit$studyFitList
    }
    e <- ctsem::ctExtract(fit)
    #dim(e$pop_DRIFT)
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
      e$pop_MANIFESTMEANS_backup <- e$pop_MANIFESTMEANS
      e$pop_MANIFESTMEANS[e$pop_MANIFESTMEANS != 0] <- 0
    }
    #
    #if ( mmRI ) { # if random intercepts are modelled as manifest means instead cint
    #  initialMeans <- round(ctsem::ctCollapse(e$pop_T0MEANS, 1, mean), digits = digits); initialMeans
    #  initialMeansSD <- round(ctsem::ctCollapse(e$pop_T0MEANS, 1, stats::sd), digits = digits); initialMeansSD
    #  initialMeansLL <- ctsem::ctCollapse(e$pop_T0MEANS, 1, function(x) stats::quantile(x, .025)); initialMeansLL
    #  initialMeansUL <- ctsem::ctCollapse(e$pop_T0MEANS, 1, function(x) stats::quantile(x, .975)); initialMeansUL
    #}
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
      e$popcov <- e$popcov_est
      message <- "Cints (slope means), T0means (initial means), and T0covs (initial (co-)vars) were inferred from a model with individually varying manifest means instead of Cints."
      popcov_mean[[i]] <- ctsem::ctCollapse(e$popcov, 1, mean)
      popcov_sd[[i]] <- ctsem::ctCollapse(e$popcov, 1, sd)
      popcov_T[[i]] <- popcov_mean[[i]]/popcov_sd[[i]]
      popcov_2.5[[i]] <- ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.025)
      popcov_97.5[[i]] <- ctsem::ctCollapse(e$popcov, 1, stats::quantile, probs=.975)
    }
  }
  return(list(popcov_mean=popcov_mean, popcov_sd=popcov_sd, popcov_T=popcov_T,
              popcov_2.5=popcov_2.5, popcov_97.5=popcov_97.5,
              message=message))
}
