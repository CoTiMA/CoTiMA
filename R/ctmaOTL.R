#' ctmaOTL
#'
#' @description Numerically determines the optimal time lag (largest effect magnitude)
#'
#' @param ctmaFitFit fit object created with ctmaFit
#' @param timeRange time range across which to search for the optimal lag
#' @param driftMat drift matrix. Either a ctmaFit object or a drift matrix has to be supplied
#' @param undoTimeScaling undos time scaling in case the the scaleTime argument was used with ctmaFit
#' @param digits digits used for rounding
#'
#' @examples
#' \dontrun{
#' OTL <- ctmaOTL(ctmaFitFit=CoTiMA::CoTiMAFullFit_6_new)
#' print(OTL)
#' }
#'
#' @importFrom  OpenMx expm
#'
#' @export ctmaOTL
#'
#' @return A corrected correlation matrix (corEmpcov). Corrections leading to r > 1.0 are set to 1.0.
#'
ctmaOTL <- function(ctmaFitFit=NULL, timeRange=NULL, driftMat=NULL, undoTimeScaling=NULL, digits=4) {

  #################################################################################################

  if ((is.null(ctmaFitFit)) & (is.null(driftMat))) {
    ErrorMsg <- "\nEither a fit object created with ctmaFit or a drift matrix has to be supplied!"
    stop(ErrorMsg)
  }

  if (!(is.null(ctmaFitFit)) & !(is.null(driftMat))) {
    ErrorMsg <- "\nEither a fit object created with ctmaFit OR a drift matrix has to be supplied. You supplied both. Delete one of them!"
    stop(ErrorMsg)
  }

  if (!(is.null(driftMat)) & (is.null(timeRange))) {
    ErrorMsg <- "\nA drift matrix was to  supplied. I also need the argument timeRange=c(start, end, stepWidth)!"
    stop(ErrorMsg)
  }

  if (!(is.null(driftMat)) ) {
    if (!(is.matrix(driftMat))) {
      ErrorMsg <- "\nThe argument driftMat was to used, but the object supplied is not a matrix!"
      stop(ErrorMsg)
    }
  }

  if (!(is.null(ctmaFitFit)) ) {
    if (is(ctmaFitFit) != "CoTiMAFit") {
      ErrorMsg <- "\nThe argument ctmaFitFit was to used, but the object supplied was not created with ctmaFit!"
      stop(ErrorMsg)
    }
  }

  if (!(is.null(driftMat)) ) {
    n.latent <- ncol(driftMat); n.latent
  } else {
    n.latent <- ctmaFitFit$n.latent; n.latent
  }

  #driftMat <- matrix(c(-1.21669, 0.05623, 0.04629, -0.33088), 2, 2, byrow=T) ; driftMat
  #is.matrix(driftMat)

  if (!(is.null(ctmaFitFit)) & (is.null(undoTimeScaling))) {
    undoTimeScaling <- TRUE
    scaleTime <- 1
    Msg <- "\nThe argument undoTimeScaling ist set to TRUE. Set undoTimeScaling=FALSE to prevent this!"
    message(Msg)
  }

  if ( !(is.null(ctmaFitFit)) & (undoTimeScaling != TRUE)) {
    scaleTime <- ctmaFitFit$argumentList$scaleTime; scaleTime
  }

  if (!(is.null(ctmaFitFit))) {
    if (is.null(timeRange)) {
      allDeltas <- ctmaFitFit$statisticsList$allDeltas; allDeltas
      maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
      if (maxDelta > 1) maxDelta <- round(maxDelta)
      timeRange <- seq(0, 3*maxDelta, (3*maxDelta)/100); timeRange
      Msg <- paste0("\nThe argument timeRange was not supplied. I take the default, which is from 0 to ", round(3*maxDelta, digits), " in 100 Steps.")
      message(Msg)
    }
  }

  if (!(is.null(ctmaFitFit)) & (undoTimeScaling == TRUE) ) {
    driftMat <- matrix(ctmaFitFit$modelResults$DRIFToriginal_time_scale, n.latent, n.latent, byrow=T); driftMat
  }
  if (!(is.null(ctmaFitFit)) & (undoTimeScaling == FALSE) ) {
    driftMat <- matrix(ctmaFitFit$modelResults$DRIFT, n.latent, n.latent, byrow=T); driftMat
  }
  Msg <- paste0("\nThe following drift matrix is used for computing the optimal time lag!")
  message(Msg)
  print(driftMat)

  #################################################################################################

  optimalCrossLag <- matrix(NA, n.latent, n.latent)
  maxCrossEffect <- matrix(NA, n.latent, n.latent)
  targetParameters <- list()
  OTL <- function(timeRange) { OpenMx::expm(driftMat * timeRange)[targetRow, targetCol] }
  tmp1 <- 0
  if (0 %in% timeRange) tmp1 <- 1
  counter <- 0
  for (j in 1:n.latent) {
    for (h in 1:n.latent) {
      if (j != h) {
        targetRow <- j
        targetCol <- h
        if (driftMat[j, h] != 0) { # an effect that is zero has no optimal lag
          counter <- counter + 1
          targetParameters[[counter]] <- sapply(timeRange, OTL); targetParameters[[counter]]
          max(abs(targetParameters[[counter]]))[1]
          maxCrossEffect[j,h] <- max(abs(targetParameters[[counter]]))[1]; maxCrossEffect[j,h]
          tmp <- which(abs(targetParameters[[counter]])==maxCrossEffect[j,h])[1] - tmp1
          targetParameters
          optimalCrossLag[j,h] <- timeRange[tmp]
        } else {
          optimalCrossLag[j,h] <- NA
        }
      }
    }
  }
  maxCrossEffect <- round(maxCrossEffect, digits); maxCrossEffect
  optimalCrossLag <- round(optimalCrossLag, digits)

return(list(timeRangeUsed=timeRange, driftMatrixUsed=driftMat, timeScaleUsed=scaleTime,
            maxCrossEffect=maxCrossEffect, optimalCrossLag=optimalCrossLag))
}
