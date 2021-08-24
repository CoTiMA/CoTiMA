#' ctmaModDrift
#'
#' @description Uses raw estimates of drift and TI effects to compute a set of different drift matrices for different moderator values
#'
#' @param CoTiMAFit A CoTiMA Fitobject.
#' @param ctStanFit A ctStanFit Fitobject.
#' @param transform A string vector containing the desired transformations. Usually NULL.
#' Could be used replace transforms contained in ctStanFit$pars.
#' @param mod.type A string vector containing either "cat" or "cont".
#' Could be used for older CoTiMA fit objects where this information is not included. Use with car.
#' @param mod.number The number of the moderator that should be used if multiple moderators are supplied
#' @param modValues A vector containing the desired values for the moderator. Usually -2:2.
#' Could be used replace possible categorical moderator values contained in CoTiMAFit.
#' @param scaleTime A scalar used to rescale the drift matrix. default = 1.
#'
#' @examples
#' \dontrun{
#' modDriftMatrices <- ctmaModDrift(CoTiMAFitObject)
#' }
#'
#' @importFrom  ctsem ctExtract
#'
#' @export ctmaModDrift
#'
#' @return A list of drift matrices for different moderator values.
#'
#'
ctmaModDrift <- function(CoTiMAFit=NULL, ctStanFit=NULL, transform=NULL, modValues=NULL, mod.type=NULL, mod.number=1, scaleTime=1) {
  if ((is.null(CoTiMAFit)) & is.null(ctStanFit) ) {
    ErrorMsg <- "\nNo CoTiMA Fit object found but required.  \nGood luck for the next try!"
    stop(ErrorMsg)
  }
  if ( (!(is.null(CoTiMAFit))) & (!(class(CoTiMAFit) == "CoTiMAFit")) ) {
    ErrorMsg <- "\nNo CoTiMA Fit provided is not of class \"CoTiMA\", which is required.  \nGood luck for the next try!"
    stop(ErrorMsg)
  }
  if ( (!(is.null(CoTiMAFit))) & (!(is.null(ctStanFit))) ) {
    ErrorMsg <- "\nBoth a CoTiMA Fit object and a ctStanFit object were provided. I am confused. Delete one of them. \nGood luck for the next try!"
    stop(ErrorMsg)
  }
  if (is.null(ctStanFit)) {
    if (is.null(CoTiMAFit$modelResults$MOD)) {
      ErrorMsg <- "\nNo The CoTiMA Fit object provided does not contain moderators.  \nGood luck for the next try!"
      stop(ErrorMsg)
    }
    if ( (is.null(CoTiMAFit$mod.type)) & (!(is.null(mod.type))) ) CoTiMAFit$mod.type <- mod.type
    if ( (is.null(CoTiMAFit$mod.type)) & is.null(mod.type) ) {
      ErrorMsg <- "\nCannot identify the type of moderator. Possible old CoTiMA Fit object? Please use mod.type argument. \nGood luck for the next try!"
      stop(ErrorMsg)
    }
    if (is.null(CoTiMAFit$studyFitList)) {
      ErrorMsg <- "\nNo The CoTiMA Fit does not contain a ctStanFit object. Possibly a reduced fit file. Cannot handle this one.  \nGood luck for the next try!"
      stop(ErrorMsg)
    }
    n.latent <- CoTiMAFit$n.latent; n.latent
    driftNames <- CoTiMAFit$parameterNames$DRIFT; driftNames
    n.mod <- n.TIpreds <- CoTiMAFit$n.moderators; n.mod; n.TIpreds
    n.TIpreds <- CoTiMAFit$n.studies-1; n.TIpreds
    mod.values <- list() # a list just for compatability with ctmaPlot code
    if ( (((all(length(modValues)==length(-2:2)) && all(modValues==-2:2)))) & (CoTiMAFit$mod.type == "cat") ) {
      mod.values[[1]] <- c(1, unique(as.numeric(substr(rownames(CoTiMAFit$summary$mod.effects), 1,2))))
    } else {
      mod.values[[1]] = modValues
    }
  }

  if (!(is.null(ctStanFit))) {
    n.latent <- ctStanFit$ctstanmodelbase$n.latent; n.latent
    driftNames <- c()
    for (i in 1:(n.latent)) {
      for (j in 1:(n.latent)) {
        driftNames <- c(driftNames, paste0("V",i,"toV", j))
      }
    }
    #driftNames
    n.mod <- ctStanFit$ctstanmodelbase$n.TIpred; n.mod
    mod.type <- "cont"
    n.TIpreds <- 0 # all TIpred are moderators in ctSatFit objects
    mod.values.tmp <- mod.values <- list()
    for (i in 1:n.mod) mod.values.tmp[[i]] <- unique(sort(ctStanFit$data$tipredsdata[, paste0("TI", i)]))
    mod.values[[1]] <- mod.values.tmp[[mod.number]]
    #mod.values
  } else {
    ctStanFit <- CoTiMAFit$studyFitList
    mod.type <- CoTiMAFit$mod.type
  }

  # the following is slightly adapted from ctmaPlot
  DRIFTCoeff <- list() # a list just for compatibility with ctmaPlot code
  counter <- 1

  for (i in mod.values[[1]]) {
    #i <- mod.values[[1]]; i
    tmp1 <- ctStanFit$stanfit$rawest[1:(n.latent^2)]; tmp1
    tmp1 <- matrix(tmp1, n.latent, byrow=TRUE); tmp1 # main effect
    # moderator effects (could be partial)
    e <- ctsem::ctExtract(ctStanFit)
    tmp2 <- apply(e$TIPREDEFFECT[,1:(n.latent^2),((n.TIpreds+1):(n.TIpreds+n.mod))], 2, mean); tmp2
    #tmp3 <- rownames(CoTiMAFit$modelResults$MOD); tmp3
    tmp3 <- paste0("Moderator_on_", driftNames)
    tmp4 <- c()
    for (l in 1:length(driftNames)) {
      tmp5 <- grep(unlist(driftNames[l]), tmp3); tmp4
      if (length(tmp5) == 0) tmp4 <- c(tmp4, NA) else tmp4 <- c(tmp4, tmp5)
    }
    tmp4[!(is.na(tmp4))] <- tmp2
    tmp4[(is.na(tmp4))] <- 0
    tmp2 <- matrix(tmp4, n.latent, byrow=TRUE); tmp2 # raw moderator effect to be added to raw main effect (followed by tform)

    DRIFTCoeff[[counter]] <- matrix(NA, n.latent, n.latent)
    if (mod.type == "cont") {
      DRIFTCoeff[[counter]] <- tmp1 + unlist(mod.values[[1]])[counter] * tmp2; DRIFTCoeff[[1]]
      #names(DRIFTCoeff[[1]]) <- paste0("Moderator Value = ", mod.values, " SD from mean if standardized (default setting)")
    }
    if (mod.type == "cat") {
      if (i == 1) {
        DRIFTCoeff[[counter]] <- tmp1 # copy main effects (= comparison group)
      } else {
        tmp2 <- CoTiMAFit$summary$mod.effects[,1]; tmp2
        tmp2 <- tmp2[((i - 2) * n.latent^2 + 1): ((i - 2) * n.latent^2 + 0 + n.latent^2)]; tmp2
        tmp2 <- matrix(tmp2, n.latent, n.latent, byrow=TRUE); tmp2
        DRIFTCoeff[[counter]] <- tmp1 + (i-1) * tmp2; DRIFTCoeff[[counter]]
      }
      #names(DRIFTCoeff[[DRIFTCoeff]][[counter]]) <- paste0("Moderator Value = ", mod.values[[g]][counter])
    }
    counter <- counter +1
  }

  ## new: apply tform to drift elements that should be tformed (extracted into tansforms)
  tmp1a <- ctStanFit$ctstanmodelbase$pars[, "transform"]; tmp1a
  tmp1b <- ctStanFit$ctstanmodelbase$pars[, "param"]; tmp1b
  transforms <- tmp1a[grep("toV", tmp1b)]; transforms
  if (length(transforms) == 0) transforms <- tmp1a[grep("drift", tmp1b)]; transforms # in case drift labels are 'auto'

  for (k in 1:(length(DRIFTCoeff))) {
    counter <- 0
    for (l in 1:(n.latent)) {
      for (m in 1:(n.latent)) {
        counter <- counter + 1
        param <- DRIFTCoeff[[k]][l,m]; param
        DRIFTCoeff[[k]][l,m] <- eval(parse(text=transforms[counter])) * scaleTime
      }
    }
  }

  #DRIFTCoeff
  invisible(DRIFTCoeff)
}

