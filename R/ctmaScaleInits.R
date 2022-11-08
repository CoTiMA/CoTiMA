#' ctmaScaleInits
#'
#' @description This function rescales inits for drifts and sets all other inits to 0 (because it is too complicated to re-scale inits for diffusions).It uses the internal trasnformations of \code{\link{ctStanFit}} (i.e., tforms) to transform the raw estimates, then re-scale them, and finally use the inverse of tfrom to supplie raw estimates as inits.
#'
#' @param CoTiMAFit Fit object created with \code{\link{ctmaFit}}
#' @param ctsemFit Fit object created with \code{\link{ctStanFit}}
#' @param newTimeScale New Time scale  \code{\link{ctStanFit}}
#' @param autoRefit Whether to automatically refit the original model using the new inits
#
#' @importFrom stats uniroot
#' @importFrom ctsem ctStanFitUpdate ctExtract
#'

ctmaScaleInits <- function(CoTiMAFit=NULL, ctsemFit=NULL, newTimeScale=NULL, autoRefit=FALSE) {
  if (is.null(ctsemFit) & is.null(CoTiMAFit) ) {
    ErrorMsg <- "\nNo fit object has been supplied. Good luck for the next try!"
    stop(ErrorMsg)
  }

  if (is.null(ctsemFit)) ctsemFit <- CoTiMAFit$studyFitList

  if (is.null(newTimeScale) ) {
    ErrorMsg <- "\nNo newTimeScale has been supplied. Good luck for the next try!"
    stop(ErrorMsg)
  }

  inits <- ctsemFit$stanfit$rawest

  # this is a copy of ctsem's priorchecker, which is not exported by ctsem
  priorchecker <- function (sf, pars = c("rawpopmeans", "rawpopsdbase", "tipredeffectparams"),
            digits = 2)
  {
    e = ctExtract(sf)
    funcs <- c(base::mean, stats::sd)
    pars = unlist(lapply(pars, function(x) if (!is.null(dim(e[[x]])))
      x))
    out = round(do.call(cbind, lapply(funcs, function(fn) do.call(c,
                                                                  lapply(pars, function(obji) apply(e[[obji]], 2, fn, na.rm = TRUE))))),
                digits)
    rownames(out) = do.call(c, c(lapply(pars, function(obji) paste0(obji,
                                                                    "_", 1:ncol(e[[obji]])))))
    out = data.frame(out, do.call(c, c(lapply(pars, function(obji) 1:ncol(e[[obji]])))))
    out = data.frame(out, do.call(c, c(lapply(pars, function(obji) rep(obji,
                                                                       ncol(e[[obji]]))))), stringsAsFactors = FALSE)
    colnames(out) = c("mean", "sd", "param", "object")
    rownames(out) = getparnamesfromraw(priorcheck = out, sf = sf)
    return(out)
  }

  # the is ctsem's getparnamesfromraw
  getparnamesfromraw <- function (priorcheck, sf)
  {
    newnames = rownames(priorcheck)
    for (ni in 1:nrow(priorcheck)) {
      if (priorcheck$object[ni] %in% "rawpopmeans") {
        newnames[ni] = paste0("rawpop_", sf$setup$popsetup$parname[sf$setup$popsetup$param %in%
                                                                     priorcheck$param[ni]][1])
      }
      if (priorcheck$object[ni] %in% "tipredeffectparams") {
        newnames[ni] = paste0("rawtipredeffect_", paste0(which(sf$standata$TIPREDEFFECTsetup ==
                                                                 priorcheck$param[ni], arr.ind = TRUE), collapse = "_"))
      }
    }
    return(newnames)
  }

  indices <- priorchecker(ctsemFit)
  tmp1 <- grep("toV",rownames(indices)); tmp1
  allInitPos <- 1:length(inits)

  tmp2 <- which(ctsemFit$ctstanmodelbase$pars$matrix == "DRIFT")
  transforms <- ctsemFit$ctstanmodelbase$pars[tmp2, "transform"]; transforms
  # inverse of tforms
  transform_inverse <- list()
  inverse <- function (f, lower = -100, upper = 100) {
    function (y) stats::uniroot((function (x) f(x) - y), lower = lower, upper = upper, extendInt = "no", maxiter = 100000)[1]
  }


  for (i in 1:length(transforms)) {
    transform_inverse[[i]] <- inverse(function (param) eval(parse(text=transforms[i])), -100, 100)
  }

  initsTmp <- inits[c(tmp1)]; initsTmp

  rawToTrans <- c()
  for (i in 1:length(initsTmp)) {
    param <- initsTmp[i]
    rawToTrans[i] <- eval(parse(text=transforms[i]));
  }

  rawToTrans <- rawToTrans * 1/newTimeScale; rawToTrans
  initsTmp2 <- c()
  for (i in 1:length(rawToTrans)) {
    param <- rawToTrans[i]
    initsTmp2[i] <- unlist(transform_inverse[[i]](param)$root)
  }

  #initsTmp2
  # added 13 Sep 2022
  initsTmp <- inits;initsTmp
  initsTmp[tmp1] <- initsTmp2

  #n <- length(inits) - length(tmp1); n
  #initsTmp[!(allInitPos %in% tmp1)] <- runif(n, -.01, .01) # random
  # removed 13 Sep 2022
  #initsTmp[!(allInitPos %in% tmp1)] <- 0

  newFit <- list(comment="autoRefit was set to FALSE")
  if  (autoRefit == TRUE & (!(is.null(CoTiMAFit))) ) {
    #ctsemFit$data$time <- ctsemFit$data$time * newTimeScale
    ctsemFit$standata$time <-  ctsemFit$standata$time * newTimeScale
    newFit <- ctStanFitUpdate(oldfit=ctsemFit, refit=TRUE, inits=initsTmp)
  }
  return(list(initsTmp, newFit))
}



