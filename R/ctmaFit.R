#' ctmaFit
#'
#' @description Fits a ctsem model with invariant drift effects across primary studies, possible multiple moderators (but all of them of the
#' the same type, either "cont" or "cat"), and possible cluster (e.g., countries where primary studies were conducted).
#'
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param activeDirectory  defines another active directory than the one used in ctmaInitFit
#' @param allInvModel estimates a model with all parameters invariant (DRIFT, DIFFUSION, T0VAR) if set TRUE (defautl = FALSE)
#' @param binaries which manifest is a binary. Still experimental
#' @param catsToCompare when performing contrasts for categorical moderators, the categories (values, not positions) for which effects are set equal
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param cint default 'auto' (= 0). Are set free if random intercepts model with varying cints is requested (by indvarying='cint')
#' @param cluster  vector with cluster variables (e.g., countries). Has to be set up carfully. Will be included in \code{\link{ctmaPrep}} in later 'CoTiMA' versions.
#' @param coresToUse if negative, the value is subtracted from available cores, else value = cores to use
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function
#' @param ctmaInitFit object to which all single ctsem fits of primary studies has been assigned to (i.e., what has been returned by \code{\link{ctmaInit}})
#' @param customPar logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param digits Number of digits used for rounding (in outputs)
#' @param drift labels for drift effects. Have to be either of the type 'V1toV2' or '0' for effects to be excluded.
#' @param driftsToCompare when performing contrasts for categorical moderators, the (subset of) drift effects analyzed
#' @param equalDrift Constrains all listed effects to be equal (e.g., equalDrift = c("V1toV2", "V2toV1")). Note that this is not required for testing the assumption that two effects are equal in the population. Use the invariantDrift argument and then \code{\link{ctmaEqual}})
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param fit TRUE (default) fits the requested model. FALSE returns the \code{\link{ctsem}} code CoTiMA uses to set up the model, the ctsemmodelbase which can be modified to match users requirements, and the data set (in long format created). The model can then be fitted using \code{\link{ctStanFit}})
#' @param ind.mod.names vector of names for individual level (!) moderators used in output
#' @param ind.mod.number which in the vector of individual level (!) moderator values shall be used (e.g., 2 for a single moderator or 1:3 for 3 moderators simultaneously)
#' @param ind.mod.type 'cont' or 'cat' of the individual level (!) moderators (mixing them in a single model not yet possible)
#' @param indVarying allows continuous time intercepts to vary at the individual level (random intercepts model, accounts for unobserved heterogeneity)
#' @param indVaryingT0 (default = NULL). Automatically set to TRUE if not set to FALSE if indVarying ist set TRUE. indVaryingT0=TRUE forces T0MEANS (T0 scores) to vary interindividually, which undos the nesting of T0(co-)variances in primary studies. Was standard until Aug. 2022. Could provide better estimates if set to FALSE.
#' @param inits vector of start values
#' @param invariantDrift  drift labels for drift effects that are set invariant across primary studies (default = all drift effects).
#' @param iter number of iterations (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param manifestMeans default = 0. Are automatically set free is indvarying is set to TRUE. Can be assigned labels to estimate them freely.
#' @param manifestVars define the error variances (default = 0) of the manifests with a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param mod.names vector of names for moderators used in output
#' @param mod.number which in the vector of moderator values shall be used (e.g., 2 for a single moderator or 1:3 for 3 moderators simultaneously)
#' @param mod.type 'cont' or 'cat' (mixing them in a single model not yet possible)
#' @param moderatedDrift labels for drift effects that are moderated (default = all drift effects)
#' @param modsToCompare when performing contrasts for categorical moderators, the moderator numbers (position in mod.number) that is used
#' @param nopriors Deprecated, but still working. If TRUE, any priors are disabled – sometimes desirable for optimization
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param primaryStudyList  could be a list of primary studies compiled with \code{\link{ctmaPrep}} that defines the subset of studies in ctmaInitFit that should actually be used
#' @param priors if FALSE, any priors are disabled – sometimes desirable for optimization
#' @param sameInitialTimes Only important for raw data. If TRUE (default=FALSE), T0MEANS occurs for every subject at the same time, rather than just at the earliest observation.
#' @param scaleClus scale vector of cluster indicators - TRUE (default) yields avg. drift estimates, FALSE yields drift estimates of last cluster
#' @param scaleMod scale moderator variables - TRUE (default) recommended for continuous and categorical moderators, to separate withing and betwen efeccts
#' @param scaleTI scale TI predictors - not recommended until version 0.5.3.1. Does not change aggregated results anyways, just interpretation of effects for dummies representing primary studies.
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param T0var (default = 'auto')
#' @param transfMod more general option to change moderator values. A vector as long as number of moderators analyzed (e.g., c("mean(x)", "x - median(x)"))
#' @param useSampleFraction to speed up debugging. Provided as fraction (e.g., 1/10).
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging

#'
#' @importFrom  RPushbullet pbPost
#' @importFrom  parallel detectCores
#' @importFrom  ctsem ctWideToLong ctDeintervalise ctModel ctStanFit ctCollapse
#' @importFrom  OpenMx vech2full expm
#' @importFrom openxlsx addWorksheet writeData createWorkbook openXL saveWorkbook
#' @importFrom  stats cov2cor quantile sd
#'
#' @export ctmaFit
#'
#' @examples
#' \dontrun{
#' # Example 1. Fit a CoTiMA to all primary studies previously fitted one by one
#' # with the fits assigned to CoTiMAInitFit_6
#' CoTiMAFullFit_6 <- ctmaFit(ctmaInitFit=CoTiMAInitFit_6)
#' summary(CoTiMAFullFit_6)
#' }
#'
#' @examples
#' \dontrun{
#' # Example 2. Fit a CoTiMA with only 2 cross effects invariant (not the auto
#' # effects) to all primary studies previously fitted one by one with the fits
#' # assigned to CoTiMAInitFit_6
#' CoTiMAInitFit_6$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAFullInv23Fit_6 <- ctmaFit(ctmaInitFit=CoTiMAInitFit_6,
#'                         invariantDrift=c("V1toV2", "V2toV1"))
#' summary(CoTiMAFullInv23Fit_6)
#' }
#'
#' @examples
#' \dontrun{
#' # Example 3. Fit a moderated CoTiMA
#' CoTiMAInitFit_6$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAMod1onFullFit_6 <- ctmaFit(ctmaInitFit=CoTiMAInitFit_6,
#'                                  mod.number=1, mod.type="cont",
#'                                  mod.names=c("Control"))
#' summary(CoTiMAMod1onFullFit_6)
#' }
#'
#' @return ctmaFit returns a list containing somearguments supplied, the fitted model, different elements summarizing the main results,
#' model type, and the type of plot that could be performed with the returned object. The arguments in the returned object are activeDirectory,
#' coresToUse, moderator names (mod.names), and moderator type (mod.type). Further arguments, which are just copied from the init-fit object
#' supplied, are, n.latent, studyList, parameterNames, and statisticsList. The fitted model is found in studyFitList, which is a large list
#' with many elements (e.g., the ctsem model specified by CoTiMA, the rstan model created by ctsem, the fitted rstan model etc.). Further
#' results returned are n.studies = 1 (required for proper plotting), data (created pseudo raw data), and a list with modelResults (i.e.,
#' DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef, MOD=modTI_Coeff,  and
#' CLUS=clusTI_Coeff). Possible invariance constraints are included in invariantDrift. The number of moderators simultaneously analyzed are
#' included in ' n.moderators. The most important new results are returned as the list element "summary", which is printed if the summary
#' function is applied to the returned object. The summary list element comprises "estimates" (the aggregated effects), possible
#' randomEffects (not yet fully working),  the minus2ll value and its n.parameters, the opt.lag sensu Dormann & Griffin (2015) and the
#' max.effects that occur at the opt.lag, clus.effects and mod.effects, and possible warning messages (message). Plot type is
#' plot.type=c("drift") and model.type="stanct" ("omx" was deprecated).
#'
ctmaFit <- function(
    activateRPB=FALSE,
    activeDirectory=NULL,
    allInvModel=FALSE,
    binaries=NULL,
    catsToCompare=NULL,
    chains=NULL,
    cint=0,
    cluster=NULL,
    coresToUse=c(2),
    CoTiMAStanctArgs=NULL,
    ctmaInitFit=NULL,
    customPar=FALSE,
    digits=4,
    drift=NULL,
    driftsToCompare=NULL,
    equalDrift=NULL,
    finishsamples=NULL,
    fit=TRUE,
    ind.mod.names=NULL,
    ind.mod.number=NULL,
    ind.mod.type="cont",
    indVarying=FALSE,
    indVaryingT0=NULL,
    inits=NULL,
    invariantDrift=NULL,
    iter=NULL,
    lambda=NULL,
    manifestMeans=0,
    manifestVars=0,
    mod.names=NULL,
    mod.number=NULL,
    mod.type="cont",
    moderatedDrift=NULL,
    modsToCompare=NULL,
    nopriors=TRUE,
    optimize=TRUE,
    primaryStudyList=NULL,
    priors=FALSE,
    sameInitialTimes=FALSE,
    scaleClus=TRUE,
    scaleMod=NULL,
    scaleTI=TRUE,
    scaleTime=NULL,
    T0means=0,
    T0var='auto',
    transfMod=NULL,
    useSampleFraction=NULL,
    verbose=NULL
)


{  # begin function definition (until end of file)

  # adapt display of information during model fit
  if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMAStanctArgs$verbose}

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object has to be supplied to plot something. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (!(is.null(invariantDrift))) { # added 12.7.2023
    # check if invariantDrift == 'none', which is used to mimic ctmaInit
    if ( ((invariantDrift[1] == "none") | (invariantDrift[1] == "None") | (invariantDrift[1] == "NONE"))  & (scaleTI == TRUE) ) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      Msg <- "The argument invariantDrift=\'none\' was used. I assume you want to mimic ctmaInit with all drift effects varying across primary studis.
    Therefore, I set the argument scaleTI=FALSE."
      message(Msg)
      scaleTI <- FALSE
    }
  }

  # CHD added Aug 2023 because on github nopriors was replaced by priors argument
  tmp1 <- formals(ctsem::ctStanFit)
  if (is.na(tmp1$nopriors)) {
    nopriors <-NA
    CoTiMAStanctArgs$nopriors <- NA
  }


  { # set fitting params
    # Added 17. Aug 2022
    tmp1 <- names(CoTiMA::CoTiMAStanctArgs) %in% names(CoTiMAStanctArgs); tmp1
    tmp2 <- CoTiMA::CoTiMAStanctArgs
    if (!(is.null(CoTiMAStanctArgs))) tmp2[tmp1] <- CoTiMAStanctArgs
    CoTiMAStanctArgs <- tmp2

    if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
    if (!(is.null(scaleClus))) CoTiMAStanctArgs$scaleClus <- scaleClus
    if (!(is.null(scaleMod))) CoTiMAStanctArgs$scaleMod <- scaleMod
    if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
    if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
    if ( (!(is.null(nopriors))) & (!(is.null(nopriors))) ) CoTiMAStanctArgs$nopriors <- nopriors # changed Aug 2023
    if (!(is.null(priors))) CoTiMAStanctArgs$priors <- priors # added Aug 2023
    if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
    if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
    if (!(is.null(iter))) CoTiMAStanctArgs$iter <- iter
    if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
  }

  { # check of catsToCompare or catsToCompare is used only with mod.type = "cat"
    if ( ((!(is.null(modsToCompare))) & (mod.type != "cat")) |
         ((!(is.null(catsToCompare))) & (mod.type != "cat")) ) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      ErrorMsg <- "The arguments modsToCompare or catsToCompare are only allowed for mod.typ=\"cat\" "
      stop(ErrorMsg)
    }
  }
  if ( (!(is.null(catsToCompare))) & (is.null(modsToCompare)) ) modsToCompare <- 1


  { # check if scaleMod is not used in combination with transfMod
    if ( (!(is.null(scaleMod))) & (!(is.null(transfMod))) ) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      ErrorMsg <- "The arguments scaleMod and transfMod cannot be used in combination. Set one of them NULL (leave out)."
      stop(ErrorMsg)
    }
  }


  #######################################################################################################################
  ####### Copy/Change INIT File based on information delivered by different PREP files (e.g., moderator studies ) #######
  #######################################################################################################################

  # if primary study list is provided in addition to initfit-object, take primary study information from  primary study
  if (!(is.null(primaryStudyList))) {
    ctmaTempFit <- ctmaInitFit
    targetStudyNumbers <- unlist(primaryStudyList$studyNumbers); targetStudyNumbers; length(targetStudyNumbers)
    if (any(is.na(targetStudyNumbers))) targetStudyNumbers <- targetStudyNumbers[-which(is.na(targetStudyNumbers))]
    #targetStudyNumbers <- targetStudyNumbers[-which(is.na(targetStudyNumbers))]; targetStudyNumbers
    for (i in (length(ctmaTempFit$studyFitList)):1) {
      if (!(ctmaTempFit$studyList[[i]]$originalStudyNo %in% targetStudyNumbers)) {
        ctmaTempFit$studyList[[i]] <- NULL
        ctmaTempFit$studyFitList[[i]] <- NULL
        ctmaTempFit$emprawList[[i]] <- NULL
        ctmaTempFit$statisticsList$originalStudyNumbers[i] <- NA
        ctmaTempFit$statisticsList$allSampleSizes[i+1] <- NA
        ctmaTempFit$statisticsList$allTpoints[i] <- NA
        ctmaTempFit$modelResults[[1]][[i]] <- NULL
        ctmaTempFit$modelResults[[2]][[i]] <- NULL
        ctmaTempFit$modelResults[[3]][[i]] <- NULL
      }
    }
    ctmaTempFit$n.studies <- length(targetStudyNumbers); ctmaTempFit$n.studies
    ctmaTempFit$statisticsList$allDeltas <- unlist(lapply(ctmaTempFit$studyList, function(extract) extract$delta_t))
    ctmaTempFit$statisticsList$minDelta <- min(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$maxDelta <- max(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$meanDelta <- mean(ctmaTempFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaTempFit$statisticsList$overallSampleSize <- sum(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$meanSampleSize <- mean(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$maxSampleSize <- max(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$minSampleSize <- min(ctmaTempFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaTempFit$statisticsList$overallTpoints <- sum(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$meanTpoints <- mean(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$maxTpoints <- max(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$statisticsList$minTpoints <- min(ctmaTempFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaTempFit$summary$model <- "Moderator Model (for details see model summary)"
    tmpStudyNumber <- as.numeric(gsub("Study No ", "", rownames(ctmaTempFit$summary$estimates))); tmpStudyNumber
    targetRows <- which(tmpStudyNumber %in% targetStudyNumbers); targetRows; length(targetRows)
    ctmaTempFit$summary$estimates <- ctmaTempFit$summary$estimates[targetRows, ]
    ctmaTempFit$summary$confidenceIntervals <- ctmaTempFit$summary$confidenceIntervals[targetRows, ]
    ctmaTempFit$summary$n.parameters <- ctmaTempFit$studyFitList[[1]]$resultsSummary$npars * length(targetRows)
    ctmaTempFit$statisticsList$originalStudyNumbers <-
      ctmaTempFit$statisticsList$originalStudyNumbers[which(!(is.na(ctmaTempFit$statisticsList$originalStudyNumbers)))]
    ctmaTempFit$statisticsList$allSampleSizes <-
      ctmaTempFit$statisticsList$allSampleSizes[which(!(is.na(ctmaTempFit$statisticsList$allSampleSizes)))]
    ctmaTempFit$statisticsList$allTpoints <-
      ctmaTempFit$statisticsList$allTpoints[which(!(is.na(ctmaTempFit$statisticsList$allTpoints)))]

    ctmaInitFit <- ctmaTempFit
  }


  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################

  {
    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
    }

    if (coresToUse >= parallel::detectCores()) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- parallel::detectCores() - 1
      Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash."
      message(Msg)
    }
  }

  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {
    if (is.null(ctmaInitFit$n.latent)) {
      n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    } else {
      n.latent <- ctmaInitFit$n.latent
    }
    if (!(is.null(ctmaInitFit$n.manifest))) n.manifest <- ctmaInitFit$n.manifest else n.manifest <- n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory

    binaries.orig <- binaries
    if (is.null(binaries)) {
      binaries <- rep(0, max(n.manifest, n.latent))
    }
    if (length(binaries) != max(n.manifest, n.latent)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nThe number of binaries provided is incorrect! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ( (!(is.null(binaries.orig))) & (indVarying != 'cint') & (indVarying != 'Cint') & (indVarying != 'CINT') ) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nYou specified binary variables. You also need to specify \"indvarying=\'cint\'\". \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.null(binaries.orig))) message("Effects of binaries on cints not implemented yet.")

    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    usedTimeRange <- seq(0, 3*maxDelta, 1); usedTimeRange # new 8.7.2022
    lambda <- ctmaInitFit$statisticsList$lambda; lambda
  }


  # check match between cluster vector and n.studies
  if (!(is.null(cluster))) {
    if (length(cluster) != n.studies) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nThe vector of cluster numbers does not match the number of primary studies.\nGood luck for the next try!"
      stop(ErrorMsg)
    }
  }

  # check moderator information
  {
    #if (!(is.null(ind.mod.number))) {
    if (length(ind.mod.number) == 1) { # CHD 30.8.2023
      if (ind.mod.number == 0) {  # CHD 16.8.2023
        ind.mod.number <- NULL
        n.ind.moderators <- 0
      }
    }
    n.ind.moderators <- length(ind.mod.number); n.ind.moderators
    #if (ind.mod.number == 0 )  n.ind.moderators <- 0 # CHD 27.6. 2023
    #if ( is.null(ind.mod.number == 0 ) ) n.ind.moderators <- 0
    if (n.ind.moderators > 0) {
      mod.number <- NULL
      mod.type=ind.mod.type
      mod.names <- NULL

    }
    #n.ind.moderators

    if (n.ind.moderators == 0) { # proceed if only moderators at the study level are used
      n.moderators <- length(mod.number); n.moderators
      { # check if transfMod is as long as n.moderators
        if ( (!(is.null(transfMod))) ) {
          if ( length(transfMod) != n.moderators ) {
            if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
            ErrorMsg <- "More transformations for moderators (transfMod) provided than moderators."
            stop(ErrorMsg)
          }
        }
      }
      if (n.moderators > 0) {
        currentModerators <- matrix(as.numeric(unlist(lapply(ctmaInitFit$studyList, function(extract) extract$moderators[mod.number]))), ncol=n.moderators); currentModerators
        if (!(is.null(primaryStudyList))) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudyList$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE)
        #if (!is.null(primaryStudyList)) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudyList$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE); currentModerators
        if (is.na((currentModerators[length(currentModerators)])[[1]][1])) currentModerators <- currentModerators[-dim(currentModerators)[1],]; currentModerators
        if (is.null(dim(currentModerators)[1])) currentModerators <- matrix(currentModerators, ncol=1); currentModerators

        if (any(is.na(currentModerators)) == TRUE) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          ErrorMsg <- "\nAt least one of the primary studies does not have a valid value for the requested moderator. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
        if (var(currentModerators) == 0) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          ErrorMsg <- "\nModerator is constant across cases.\nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
    }

    if (n.ind.moderators > 0) { #
      n.moderators <- n.ind.moderators; n.moderators
      tmpMods <- lapply(ctmaInitFit$ind.mod.List, function(x) x)
      if(!(is.matrix(tmpMods[[1]]))) tmpMods <- lapply(tmpMods, function(x) matrix(x, ncol=1))
      tmp1 <- unlist(lapply(tmpMods, function(x) {
        if (is.null(dim(x))) return(1) else return(ncol(x))
        #nrow(x)[2]
        } )); tmp1
      if (length(ind.mod.number) > min(tmp1)) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
        ErrorMsg <- "\nIndividual level moderation model requested. At least one study has fewer individual level moderators in the raw data file than the number
          specified via ind.mod.number . Check your data, redo ctmaPrep and ctmaInit again. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (!(is.null(primaryStudyList))) {
        #if (length(primaryStudyList$deltas) != n.studies) {
        if (length(primaryStudyList$deltas) > n.studies) { # CHD changed Aug 2023
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
          ErrorMsg <- "\nYou provided a list to the argument primaryStudyList. The number of studies in this list is not identical to the number of
          studies fitted with ctmaInit as provided in the argument ctmaInitFit. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
      currentModerators <- tmpMods[[1]]
      for ( i in 2:length(tmpMods)) currentModerators <- rbind(currentModerators, tmpMods[[i]])
      #currentModerators <- currentModerators[, ind.mod.number]
      currentModerators <- as.matrix(currentModerators)
      colnames(currentModerators) <- paste0("mod", 1:dim(currentModerators)[2])
    }
  }

  #######################################################################################################################
  ################################################# data preparation ####################################################
  #######################################################################################################################

  {
    # manifests & latent names
    if (n.manifest > n.latent) {
      manifestNames <- paste0("y", 1:n.manifest); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    } else {
      manifestNames <- paste0("V", 1:n.latent); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    }

    # combine pseudo raw data
    {
      if (n.moderators > 0) {
        if (n.ind.moderators != 0) {
          #listOfStudyFits=ctmaInitFit
          tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit)
          casesToDelete <- which(is.na(currentModerators[, ind.mod.number])); casesToDelete
          if (length(casesToDelete) > 0)  {
            currentModerators <- as.matrix(currentModerators[-casesToDelete,])
            tmp$groups <- tmp$groups[-casesToDelete]; length(tmp$alldata)
            tmp$alldata <- tmp$alldata[-casesToDelete,]; dim(tmp$alldata)
          }
          tmp$moderatorGroups <- currentModerators; dim(tmp$moderatorGroups)
        } else {
          tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit, moderatorValues= currentModerators)
        }
      } else {
        tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit)
      }
      datawide_all <- tmp$alldata
      groups <- tmp$groups # vector of study IDs
    }

    # delete cases with missing moderator values (does not apply for ind level mods because this is already done before )
    if (n.moderators > 0) {
      casesToDelete <- tmp$casesToDelete; casesToDelete
      if (!(is.null(casesToDelete))) {
        datawide_all <- datawide_all[-casesToDelete, ]
        groups <- groups[-casesToDelete]
        n.studies <- length(unique(groups))
        currentModerators <- currentModerators[-casesToDelete,]
      }
    }

    # make data matrix with moderators
    if (n.moderators > 0) {
      moderatorGroups <- tmp$moderatorGroups
      if (!(is.matrix(moderatorGroups))) moderatorGroups <- matrix(moderatorGroups, ncol=1)
      colnames(moderatorGroups) <- paste0("mod", 1:(dim(currentModerators)[2])); colnames(moderatorGroups)

      if (n.ind.moderators != 0 ) {
        tmp1 <- as.matrix(moderatorGroups[, ind.mod.number], ncol=1)
        colnames(tmp1) <- colnames(moderatorGroups)[1]
        moderatorGroups <- tmp1
        #tmp$moderatorGroups <- tmp$moderatorGroups
        tmp$moderatorGroups <- moderatorGroups # CHD 31. AUG 2023
      }

      # change category values for making specific contrasts
      if (!(is.null(catsToCompare))) { # change category values to make the cats to compare the largest ones

        ### check if comparison of requested categories is possible
        tmp1 <- as.numeric(names(table(moderatorGroups))); tmp1
        counter <- 0
        for (c1 in catsToCompare) if (c1 %in% tmp1) counter <- counter + 1
        if (counter < length(catsToCompare)){
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          ErrorMsg <- "\nAt least one of the categories that should be compared is not available in the data. \nGood luck for the next try!"
          stop(ErrorMsg)
        }

        # check if each requested moderators is available in more than 1 study
        for (c1 in catsToCompare) if (!(c1 %in% as.numeric(names(table(currentModerators))))) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          ErrorMsg <- "\nOne of the categories that should be compared exists in a single primary study only. This model is not identified. Eliminate the primary study and redo ctmaPrep.. \nGood luck for the next try!"
          stop(ErrorMsg)
        } #

        for (c1 in 1:dim(moderatorGroups)[2]) {
          if (c1 %in% modsToCompare) {
            #maxModeratorGroups <- max(moderatorGroups[, c1]); maxModeratorGroups
            maxAbsModeratorGroups <- max(abs(moderatorGroups[, c1])) ; maxAbsModeratorGroups
            minAbsModeratorGroups <- min(abs(moderatorGroups[, c1])) ; minAbsModeratorGroups
            moderatorGroupsOffset <- maxAbsModeratorGroups + minAbsModeratorGroups; moderatorGroupsOffset
            for (c2 in catsToCompare) {
              #moderatorGroups[moderatorGroups[ ,c1] == c2] <- moderatorGroups[moderatorGroups[ ,c1] == c2] + maxModeratorGroups
              moderatorGroups[moderatorGroups[ ,c1] == c2] <- moderatorGroups[moderatorGroups[ ,c1] == c2] - moderatorGroupsOffset
            }
          }
        }
      }  # END if (!(is.null(catsToCompare)))

    }

    # label group vectors
    names(groups) <- c("Study_No_"); groups
    groupsNamed <- (paste0("Study_No_", groups)); groupsNamed

    # augment pseudo raw data by group ID and moderators
    if (n.moderators > 0) {
      dataTmp <- cbind(datawide_all, groups, moderatorGroups)
    } else {
      dataTmp <- cbind(datawide_all, groups)
    }


    # make TI out of group membership
    for (i in 1:(n.studies-1)) {
      tmp <- matrix(0, nrow=nrow(dataTmp)); tmp
      colnames(tmp) <- paste0("TI", i); tmp
      dataTmp <- cbind(dataTmp, tmp); dim(dataTmp)
      tmp <- which(dataTmp[,"groups"] == i); tmp
      dataTmp[tmp, ncol(dataTmp)] <- 1
      if (CoTiMAStanctArgs$scaleTI == TRUE) dataTmp[ , ncol(dataTmp)] <- scale(dataTmp[ , ncol(dataTmp)])
    }
    targetCols <- which(colnames(dataTmp) == "groups"); targetCols
    dataTmp <- dataTmp[ ,-targetCols]

    # make TI out of moderators
    modTIstartNum <- n.studies; modTIstartNum
    if (n.moderators > 0) {
      # make TI predictors as replacements for each moderator
      if (mod.type=="cont") {
        tmp1 <- paste0("mod", 1:n.moderators); tmp1
        if (length(tmp1) == 1) tmp <- matrix(dataTmp[ , tmp1], ncol=length(tmp1)) else tmp <- dataTmp[ , tmp1]
        if (CoTiMAStanctArgs$scaleMod == TRUE) tmp[ , 1:ncol(tmp)] <- scale(tmp[ , 1:ncol(tmp)])
        if (!(is.null(transfMod))) {
          tmp2 <- tmp #[ , 1:ncol(tmp)]
          for (t in 1:length(transfMod)) {
            x <- tmp2[, t]
            tmp2[, t] <- as.numeric(eval(parse(text=transfMod[t])))
          }
          tmp[ , 1:ncol(tmp)] <- tmp2
        }
        currentStartNumber <- modTIstartNum; currentStartNumber
        currentEndNumber <- currentStartNumber + n.moderators-1; currentEndNumber
        colnames(tmp) <- paste0("TI", currentStartNumber:currentEndNumber); tmp
        dataTmp <- cbind(dataTmp, tmp); dim(dataTmp)
        dataTmp <- dataTmp[ ,-grep("mod", colnames(dataTmp))]
      }
      if ((mod.type=="cat") | (ind.mod.type=="cat")) {
        tmp1 <- paste0("mod", 1:n.moderators); tmp1
        if (length(tmp1) == 1) tmp <- matrix(dataTmp[ , tmp1], ncol=length(tmp1)) else tmp <- dataTmp[ , tmp1]
        if (n.moderators > 1) {
          unique.mod <- list()
          for (i in 1:n.moderators) unique.mod[[i]] <- sort(c(unique(tmp[,i])))
        } else {
          unique.mod <- sort(c(unique(tmp)))
        }


        # determine number of required dummies
        if (n.moderators > 1) {
          catCounter <- 0
          for (i in 1:length(unique.mod)) catCounter <- catCounter + length(unique.mod[[i]]) -1
        } else {
          catCounter <- length(unique.mod) -1
        }

        # create dummies
        tmpTI <- matrix(0, dim(tmp)[1], catCounter)
        counter <- 0
        if (n.moderators > 1) {
          for (i in 1:n.moderators) {
            for (j in 2:length(unique.mod[[i]])) {
              counter <- counter + 1
              tmp2 <- which(tmp[, i] == unique.mod[[i]][j]); tmp2
              tmpTI[tmp2, counter] <- 1
            }
          }
        } else {
          for (i in 2:length(unique.mod)) {
            counter <- counter + 1; counter
            tmp2 <- which(tmp == unique.mod[i]); tmp2
            tmpTI[tmp2, counter] <- 1
          }
        }

        if (CoTiMAStanctArgs$scaleMod == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)], scale=FALSE)
        currentStartNumber <- modTIstartNum; currentStartNumber
        currentEndNumber <- currentStartNumber + ncol(tmpTI)-1; currentEndNumber
        colnames(tmpTI) <- paste0("TI", currentStartNumber:currentEndNumber)

        dataTmp <- cbind(dataTmp, tmpTI); dim(dataTmp)
        dataTmp <- dataTmp[ ,-grep("mod", colnames(dataTmp))]
        head(dataTmp)
      }
    } # END if (n.moderators > 0)

    # add clusters as dummy moderators
    if (!(is.null(cluster))) {
      # determine number of required dummies
      targetCluster <- which(table(cluster) > 1); targetCluster  # no cluster if only one study is included
      targetCluster <- names(targetCluster); targetCluster
      clusCounter <- length(targetCluster); clusCounter
      if (clusCounter == length(table(cluster))) clusCounter <- clusCounter -1 # if all cluster exits leave the last one as reference cluster
      # create dummies
      tmpTI <- matrix(0, dim(dataTmp)[1], clusCounter)
      for (i in 1:clusCounter) {
        targetGroups <- which(cluster == targetCluster[i]); targetGroups
        tmp2 <- which(groups %in% targetGroups); length(tmp2)
        tmpTI[tmp2, i] <- 1
      }
      # CHD changed 7.6.2023
      #if (CoTiMAStanctArgs$scaleClus == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)], scale=FALSE)
      if (scaleClus == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)], scale=FALSE)
      cluster.weights <- cluster.sizes <- matrix(NA, nrow=ncol(tmpTI), ncol=2)
      for (l in 1:ncol(tmpTI)) {
        cluster.weights[l,] <- round(as.numeric(names(table(tmpTI[ , l]))), digits)
        cluster.sizes[l,] <- table(tmpTI[ , l])
      }
      rownames(cluster.weights) <- paste0(1:ncol(tmpTI), "_on_")
      colnames(cluster.weights) <- c("non Members", "Cluster Member")
      rownames(cluster.sizes) <- rep("N", ncol(tmpTI))
      colnames(cluster.sizes) <- c("non Members", "Cluster Member")
      cluster.note <- capture.output(cat("The weights represent standardized cluster dummies. They are used to multiply a cluster's TI effect ",
                                         "and this product is then added to the average effect shown in $estimates, which overall yields",
                                         "the effects within a cluster as shown in $cluster.specific.effect.",
                                         sep="\n"))
      currentStartNumber <- n.studies; currentStartNumber
      currentEndNumber <- currentStartNumber + clusCounter -1; currentEndNumber
      colnames(tmpTI) <- paste0("TI", currentStartNumber:currentEndNumber); colnames(tmpTI)
      dataTmp <- cbind(dataTmp, tmpTI); dim(dataTmp)
    } else {
      clusCounter <- 0
      cluster.weights <- c()
      cluster.sizes <- c()
      cluster.note <- c()
    }

    # params for ctsem model specification
    {
      tmp1 <- 0
      n.all.moderators <- 0
      if (n.moderators > 0) {
        if (mod.type == "cont") tmp1 <- n.moderators
        if (mod.type == "cat") tmp1 <- catCounter
        n.all.moderators <- tmp1
      }
      if (!(is.null(cluster))) tmp1 <- tmp1 + clusCounter

      n.var <- max(n.manifest, n.latent); n.var
    }

    # possible subsample selection
    if (!(is.null(useSampleFraction))) {
      N <- dim(dataTmp)[1]; N
      stepwidth <- 100/useSampleFraction
      targetCases <- round(seq(1, N, stepwidth), 0); targetCases
      dataTmp <- dataTmp[targetCases, ]
      #groups <- groups[targetCases]
    }


    # make long data format
    {
      dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.var, n.TIpred = (n.studies-1+tmp1),
                                      manifestNames=manifestNames)

      dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
      dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
    }

    # eliminate rows where ALL latents are NA
    {
      if (n.manifest > n.latent) namePart <- "y" else namePart <- "V"
      dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0(namePart, 1:n.var)], 1, function(x) sum(is.na(x)) != n.var ), ]
      datalong_all <- dataTmp3
    }

    datalong_all <- as.data.frame(datalong_all)
  }

  # TI-identifiers for groups, moderators, and clusters
  {
    groupTIs <- paste0("TI", 1:(length(unique(groups))-1)); groupTIs
    tmp1 <- (length(unique(groups))); tmp1
    if (n.moderators > 0) modTIs <- paste0("TI", tmp1:(tmp1+n.all.moderators-1))
    tmp1 <- (tmp1+n.all.moderators); tmp1
    if (!(is.null(cluster))) clusTIs <- paste0("TI", tmp1:(tmp1+clusCounter-1))
  }


  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################


  # define Names (labels in compiled output) and Params (= labels for ctsem models)
  if (is.null(moderatedDrift) & ( (!(is.null(mod.number))) | (!(is.null(ind.mod.number)))) ) moderatedDrift <- "all" # will be changed by ctmaLabels

  namesAndParams <- ctmaLabels(
    n.latent=n.latent,
    n.manifest=n.manifest,
    lambda=lambda,
    drift=drift,
    invariantDrift=invariantDrift,
    moderatedDrift=moderatedDrift,
    equalDrift=equalDrift,
    T0means=T0means,
    manifestMeans=manifestMeans,
    manifestVars=manifestVars)

  driftNames <- namesAndParams$driftNames; driftNames
  driftFullNames <- namesAndParams$driftFullNames; driftFullNames
  driftParams <- namesAndParams$driftParams; driftParams
  diffNames <- namesAndParams$diffNames; diffNames
  diffParams <- namesAndParams$diffParams; diffParams
  diffFullNames <- namesAndParams$diffFullNames; diffFullNames
  invariantDriftNames <- namesAndParams$invariantDriftNames; invariantDriftNames
  invariantDriftParams <- namesAndParams$invariantDriftParams; invariantDriftParams
  # CHD added 28.6.2023
  if (!(is.null(invariantDrift))) { # added 12.7.2023
    if ( (invariantDrift[1] == "none") | (invariantDrift[1] == "None") | (invariantDrift[1] == "NONE")  ) {
      invariantDriftNames <- invariantDriftParams <- 'none'
    }
    if ( (invariantDrift[1] == "all") | (invariantDrift[1] == "All") | (invariantDrift[1] == "ALL")  ){
      invariantDriftNames <- invariantDriftParams <- driftFullNames
    }
  }
  moderatedDriftNames <- namesAndParams$moderatedDriftNames; moderatedDriftNames
  equalDriftNames <- namesAndParams$equalDriftNames; equalDriftNames
  equalDriftParams <- namesAndParams$equalDriftParams; equalDriftParams
  lambdaParams <- namesAndParams$lambdaParams; lambdaParams
  manifestMeansParams <- namesAndParams$manifestMeansParams; manifestMeansParams
  T0meansParams <- namesAndParams$T0meansParams; T0meansParams
  manifestVarsParams <- namesAndParams$manifestVarsParams; manifestVarsParams
  # CHD 9.6.2023
  #T0VARParams <- namesAndParams$T0VARParams; T0VARParams
  T0VARParams <- T0var; T0VARParams
  CINTParams <- cint; CINTParams

  if (is.null(invariantDriftNames)) invariantDriftNames <- driftNames

  if ( (!(is.null(invariantDrift))) &
       ( (indVarying != TRUE) | (indVarying != 'CINT') | (indVarying != 'cint' | (indVarying != 'Cint') ))  ) { # added 12.7.2023, added AUG 2023
    if ( (invariantDrift[1] == "none") | (invariantDrift[1] == "None") | (invariantDrift[1] == "NONE")  ) {

      # CHD 31. AUG 2023
      RIPos <- NA

      # CHD added AUG 2023: Start values for mimicing ctmaInit (RI not cannot be done because cov among random effects do not exist at study level)
      ctStanFitObject <- ctmaInitFit$studyFitList[[n.studies]]
      # get parameter positions in rawest vector
      tmp1 <- which(!(is.na(ctStanFitObject$ctstanmodelbase$pars$param))); tmp1
      tmpPars <- ctStanFitObject$ctstanmodelbase$pars[tmp1,]; tmpPars
      RIPos <- which(tmpPars$indvarying == TRUE); RIPos
      #
      if (length(RIPos) > 0){
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        #ErrorMsg <- "\nYou cannot mimic an cmtaInit using ctmaFit because your model has random intercepts. They cannot be estimated at the study-level with ctmaFit. \nGood luck for the next try!"
        #stop(ErrorMsg)
        # CHD 30. Aug. 2023
        Msg <- "\nYou cannot mimic an cmtaInit using ctmaFit because your model has random intercepts.
        \nCovariance among random effects is only estimated once - they do not exist at study level.
        \nI will fit the model anyway, but take care you correctly interpret the results."
        message(Msg)

      }


      # CHD 31. AUG 2023
      if (is.na(RIPos[1])) {

        print(paste0("#################################################################################"))
        print(paste0("###### Computing start values to improve convergence in mimicing ctmaInit. ######"))
        print(paste0("#################################################################################"))

        inits <- ctStanFitObject$stanfit$rawest; round(inits, 3)
        #
        T0meansPos <- which(tmpPars$matrix == "T0MEANS"); T0meansPos
        driftPos <- which(tmpPars$matrix == "DRIFT"); driftPos
        T0varPos <- which(tmpPars$matrix == "T0VAR"); T0varPos
        diffPos <- which(tmpPars$matrix == "DIFFUSION"); diffPos
        mmPos <- which(tmpPars$matrix == "MANIFESTMEANS"); mmPos
        cintPos <- which(tmpPars$matrix == "CINT"); cintPos

        rawDrift <-  modDrift <-  rawT0var <-  modT0var <-  rawDiff <- modDiff <- list()
        rawT0means <-  modT0means <-  rawMM <- modMM <- rawCint <- modCint <- list()
        # get rawest of reference study
        rawDrift[[n.studies]] <- ctStanFitObject$stanfit$rawest[driftPos]; rawDrift[[n.studies]]
        rawT0var[[n.studies]] <- ctStanFitObject$stanfit$rawest[T0varPos]; #rawT0var[[n.studies]]
        rawDiff[[n.studies]] <- ctStanFitObject$stanfit$rawest[diffPos]; rawDiff[[n.studies]]
        if (length(T0meansPos) > 0) rawT0means[[n.studies]] <- ctStanFitObject$stanfit$rawest[T0meansPos]; #rawT0means
        if (length(mmPos) > 0) rawMM[[n.studies]] <- ctStanFitObject$stanfit$rawest[mmPos]; #rawMM[[n.studies]]
        if (length(cintPos) > 0) rawCint[[n.studies]] <- ctStanFitObject$stanfit$rawest[cintPos]; #rawCint[[n.studies]]
        # rermaining studies
        for (k in 1:(n.studies-1)) {
          #k <- 1
          ctStanFitObject <- ctmaInitFit$studyFitList[[k]]
          rawDrift[[k]] <- ctStanFitObject$stanfit$rawest[driftPos]; rawDrift[[k]]
          rawT0var[[k]] <- ctStanFitObject$stanfit$rawest[T0varPos]; #rawT0var[[k]]
          rawDiff[[k]] <- ctStanFitObject$stanfit$rawest[diffPos]; rawDiff[[k]]
          if (length(T0meansPos) > 0) {rawT0means[[k]] <- ctStanFitObject$stanfit$rawest[T0meansPos]; rawT0means[[k]]}
          if (length(mmPos) > 0) {rawMM[[k]] <- ctStanFitObject$stanfit$rawest[mmPos]; rawMM[[k]]}
          if (length(cintPos) > 0) {rawCint[[k]] <- ctStanFitObject$stanfit$rawest[cintPos]; rawCint[[k]]}

          modDrift[[k]] <- rawDrift[[k]] - rawDrift[[n.studies]]; modDrift[[k]] # order of values is correct
          modT0var[[k]] <- rawT0var[[k]] - rawT0var[[n.studies]]; modT0var[[k]]
          modDiff[[k]] <- rawDiff[[k]] - rawDiff[[n.studies]]; modDiff[[k]]
          if (length(T0meansPos) > 0) modT0means[[k]] <- rawT0means[[k]] - rawT0means[[n.studies]]
          if (length(mmPos) > 0) modMM[[k]] <- rawMM[[k]] - rawMM[[n.studies]]
          if (length(cintPos) > 0) modCint[[k]] <- rawCint[[k]] - rawCint[[n.studies]]

        }

        for (l in 1:length(T0meansPos)) inits <- c(inits, unlist(lapply(modT0means, function(x) x[l])))
        for (l in 1:length(driftPos)) inits <- c(inits, unlist(lapply(modDrift, function(x) x[l])))
        for (l in 1:length(diffPos)) inits <- c(inits, unlist(lapply(modDiff, function(x) x[l])))
        for (l in 1:length(mmPos)) inits <- c(inits, unlist(lapply(modMM, function(x) x[l])))
        for (l in 1:length(cintPos)) inits <- c(inits, unlist(lapply(modCint, function(x) x[l])))
        for (l in 1:length(T0varPos)) inits <- c(inits, unlist(lapply(modT0var, function(x) x[l])))
      }
    }
  }



  if (allInvModel) {
    allInvModelFit <- ctmaAllInvFit(ctmaInitFit=ctmaInitFit,
                                    activeDirectory=activeDirectory,
                                    activateRPB=activateRPB,
                                    digits=digits,
                                    drift=drift,
                                    coresToUse=coresToUse,
                                    scaleTime=scaleTime,
                                    optimize=optimize,
                                    nopriors=nopriors,
                                    priors=priors,
                                    finishsamples=finishsamples,
                                    iter=iter,
                                    chains=chains,
                                    verbose=verbose,
                                    indVarying = indVarying,
                                    indVaryingT0 = indVaryingT0,
                                    customPar = customPar)
    fitStanctModel <- allInvModelFit$studyFitList[[1]]
    fitStanctModel_summary <- summary(fitStanctModel)
  } else {
    n.TIpred <- (n.studies-1+n.all.moderators+clusCounter); n.TIpred
    # scale Drift to cover changes in ctsem 3.4.1 (this would be for ctmaFit/ctmaModFit, but for Init individual study modification is done later)
    driftParamsTmp <- driftParams; driftParamsTmp
    diffParamsTmp  <- diffParams
    meanLag <- mean(allDeltas, na.rm=TRUE); meanLag
    if (customPar) {
      counter <- 0
      for (h in 1:(n.latent)) {
        for (j in 1:(n.latent)) {
          counter <- counter + 1
          if (h == j) {
            driftParamsTmp[counter] <- paste0(driftParamsTmp[counter], paste0("|-log1p_exp(-param *.1 -2)"))
            diffParamsTmp[counter] <- paste0(diffParamsTmp[counter], paste0("|log1p_exp(param *.1 +2)"))
          }
        }
      }
    }

    # Make model
    {
      # CHD 9.6.2023
      if ((indVarying == 'cint') | (indVarying == 'Cint')) indVarying <- 'CINT'

      # CHD 14. Jun 2023
      if ((indVarying == TRUE) & (is.null(indVaryingT0))) indVaryingT0 <- TRUE
      if ((indVarying == 'CINT') & (is.null(indVaryingT0))) indVaryingT0 <- TRUE
      if (is.null(indVaryingT0)) indVaryingT0 <- FALSE

      #if ((allInvModel == FALSE) & ((indVarying == TRUE) | (indVarying == 'cint') | (indVarying == 'CINT') ) ) {
      if (allInvModel == FALSE)  {

        # CHD 9.6.2023
        if ( (indVarying == 'CINT') & (indVaryingT0 == TRUE) ) {
          print(paste0("#################################################################################"))
          print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
          print(paste0("#################################################################################"))

          print(paste0("#################################################################################"))
          print(paste0("# T0means are set to \'auto\'. T0(co-)variances not modelled nested in primaries.##"))
          print(paste0("#################################################################################"))
          T0meansParams <- 'auto'

          print(paste0("#################################################################################"))
          print(paste0("####################### CT intercepts are set free.  ############################"))
          print(paste0("#################################################################################"))

          CINTParams <- c()
          for (c in 1:n.latent) {
            CINTParams <- c(CINTParams, paste0("cintV", c))
          }
        }
        #CINTParams

        if ( (indVarying == 'CINT') & (indVaryingT0 == FALSE) ) {
          print(paste0("#################################################################################"))
          print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
          print(paste0("#################################################################################"))

          print(paste0("#################################################################################"))
          print(paste0("### T0means are set to 0. T0(co-)variances are modelled nested in primaries. ####"))
          print(paste0("#################################################################################"))
          T0meansParams <- 0

          print(paste0("#################################################################################"))
          print(paste0("####################### CT intercepts are set free.  ########################"))
          print(paste0("#################################################################################"))

          CINTParams <- c()
          for (c in 1:n.latent) {
            CINTParams <- c(CINTParams, paste0("cintV", c))
          }
        }

        if ( (indVarying == TRUE) & (indVaryingT0 == TRUE) ) {
          print(paste0("#################################################################################"))
          print(paste0("###### Just a note: Individually varying manifest means model requested.  #######"))
          print(paste0("#################################################################################"))

          print(paste0("#################################################################################"))
          print(paste0("### T0means set to \'auto\'. T0(co-)variances not modelled nested in primaries. ###"))
          print(paste0("##### Consider setting \'indVaryingT0 = FALSE\' if estimation problems occur, #####"))
          print(paste0("###### however, be aware that this is not the regular RI model anymore then. ####"))
          print(paste0("#################################################################################"))
          T0meansParams <- 'auto'

          print(paste0("#################################################################################"))
          print(paste0("######### Manifest means (as replacement for intercepts) are set free.  #########"))
          print(paste0("#################################################################################"))

          manifestMeansParams <- 'auto'
        }


        if ( (indVarying == TRUE) & (indVaryingT0 == FALSE) ) {
          print(paste0("#################################################################################"))
          print(paste0("###### Just a note: Individually varying manifest means model requested.  #######"))
          print(paste0("#################################################################################"))

          print(paste0("#################################################################################"))
          print(paste0("### T0means are set to 0. T0(co-)variances are modelled nested in primaries. ####"))
          print(paste0("#################################################################################"))
          T0meansParams <- 0

          print(paste0("#################################################################################"))
          print(paste0("######### Manifest means (as replacement for intercepts) are set free.  #########"))
          print(paste0("#################################################################################"))

          manifestMeansParams <- 'auto'
        }

        if (!(is.null(binaries.orig))) {
          # check if really cints rather than manifest means are modelled
          # set TIpredeffects on cints to TRUE
        }

        stanctModel <- suppressMessages(
          ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, #Tpoints=maxTpoints,
                         manifestNames=manifestNames,
                         DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                         DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                         LAMBDA=lambdaParams,
                         # CHD 9.6.2023; AUG 2023
                         #CINT=CINTParams, #matrix(0, nrow=n.latent, ncol=1),
                         CINT=matrix(CINTParams, nrow=n.latent, ncol=1),
                         #T0MEANS = T0meansParams,
                         T0MEANS = matrix(T0meansParams, nrow=n.latent, ncol=1),
                         #MANIFESTMEANS = manifestMeansParams,
                         MANIFESTMEANS = matrix(manifestMeansParams, nrow=n.latent, ncol=1),
                         MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var),
                         T0VAR = T0VARParams,
                         type = 'stanct',
                         n.TIpred = n.TIpred,
                         TIpredNames = paste0("TI", 1:n.TIpred))
        )


        #CHD 13.6.2023
        if (indVaryingT0 == TRUE) {
          stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS','indvarying'] <- TRUE
        } else {
          stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS','indvarying'] <- FALSE
        }
        # CHD 13.6.2023
        if (indVarying == 'CINT') {
          stanctModel$pars[stanctModel$pars$matrix %in% 'CINT','indvarying'] <- TRUE
        } else {
          stanctModel$pars[stanctModel$pars$matrix %in% 'CINT','indvarying'] <- FALSE
        }
        if (indVarying == TRUE) {
          stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- TRUE
        } else {
          stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- FALSE
        }
      } else {
        stanctModel <- suppressMessages(
          ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, #Tpoints=maxTpoints,
                         manifestNames=manifestNames,
                         DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                         DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                         LAMBDA=lambdaParams,
                         # CHD 9.6.2023; Aug 2023
                         # CINT=CINTParams, #matrix(0, nrow=n.latent, ncol=1),
                         CINT=matrix(CINTParams, nrow=n.latent, ncol=1),
                         #T0MEANS = T0meansParams,
                         T0MEANS = matrix(T0meansParams, nrow=n.latent, ncol=1),
                         #MANIFESTMEANS = manifestMeansParams,
                         MANIFESTMEANS = matrix(manifestMeansParams, nrow=n.latent, ncol=1),
                         MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var),
                         T0VAR = T0VARParams,
                         type = 'stanct',
                         n.TIpred = n.TIpred,
                         TIpredNames = paste0("TI", 1:n.TIpred))
        )

      }
    }
    #stanctModel$pars

    # general setting for params
    stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'LAMBDA',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'CINT',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTMEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTVAR',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE


    if (!(is.null(cluster))) {
      stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames[(n.studies++n.all.moderators):(n.studies+n.all.moderators+clusCounter-1)],'_effect')] <- TRUE
    }

    if (n.moderators > 0) {
      tmp1 <- which(stanctModel$pars$matrix == "DRIFT"); tmp1
      tmp2 <- which((stanctModel$pars[tmp1, "param"] %in% moderatedDriftNames)); tmp2
      targetCols <- (n.studies):(n.studies-1+n.all.moderators); targetCols
      stanctModel$pars[ , paste0(stanctModel$TIpredNames[targetCols],'_effect')] <- FALSE
      stanctModel$pars[tmp1[tmp2] , paste0(stanctModel$TIpredNames[targetCols],'_effect')] <- TRUE

      # change model to allow for comparison of categories (-2ll is the only important result)
      if (!(is.null(catsToCompare))) {
        currentStartNumber <- modTIstartNum; currentStartNumber
        for (c1 in 1:modsToCompare) {
          #c1 <-1
          if (n.moderators > 1) {
            targetCols2 <- c()
            for (c3 in 1:length(unique.mod)) {
              targetCols2 <- c(targetCols2, currentStartNumber:(currentStartNumber+length(catsToCompare)-2))
              currentStartNumber <- currentStartNumber + length(unique.mod[[c3]]) - 1
            }
          } else {
            targetCols2 <- currentStartNumber:(currentStartNumber+length(catsToCompare)-2); targetCols2
          }
        }
        if (is.null(driftsToCompare)) driftsToCompare <- driftFullNames
        targetRows2 <- which(stanctModel$pars[, "param"] %in% driftsToCompare); targetRows2
        targetCols3 <- c()
        for (c4 in targetCols2) targetCols3 <- c(targetCols3, grep(c4, colnames(stanctModel$pars)))
        stanctModel$pars[targetRows2, targetCols3] <- FALSE


      }
    }

    # the target effects
    tmp1 <- which(stanctModel$pars$matrix == "DRIFT"); tmp1
    # CHD changes 28.6.2023; 13. =ct 2023
    tmp2 <- which(stanctModel$pars[tmp1, "param"] %in% invariantDriftParams); tmp2
    # CHD changes 12.7.2023 to include drift effects that were set to 0.0
    #varyingDrifts <- tmp1[!(tmp1 %in% tmp1[tmp2])]; varyingDrifts
    tmp3 <- which(is.na(stanctModel$pars[tmp1, "param"])); tmp3
    tmp4 <- sort(unique(c(tmp2, tmp3))); tmp4
    varyingDrifts <- tmp1[!(tmp1 %in% tmp1[tmp4])]; varyingDrifts



    if (length(varyingDrifts) > 0) stanctModel$pars[varyingDrifts, paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- TRUE
    #stanctModel$pars

    if (!(optimize)) {
      customPar <- FALSE
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      tmp1a <- paste0(" Bayesian sampling was selected, which does require appropriate scaling of time. ")
      tmp2 <- nchar(tmp1a); tmp2
      tmp3 <- (81 - tmp2)/2; tmp3
      tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
      tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
      tmp6a <- paste0(tmp4, tmp1a, tmp5); tmp6a

      tmp1b <- paste0(" See the end of the summary output ")
      tmp2 <- nchar(tmp1b); tmp2
      tmp3 <- (81 - tmp2)/2; tmp3
      tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
      tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
      tmp6b <- paste0(tmp4, tmp1b, tmp5); tmp6b

      Msg <- paste0("################################################################################# \n", tmp6a, "\n", tmp6b, "\n#################################################################################")
      message(Msg)
      #Msg <- "Bayesian sampling was selected, which does require appropriate scaling of time. See the end of the summary output\n"
      #message(Msg)
    }

    # CHD 31. AUG 2023 (not really necessary)
    tmp1 <- which(is.na(stanctModel$pars$param)); tmp1
    tmp2 <- grep("_effect", colnames(stanctModel$pars)); tmp2
    stanctModel$pars[tmp1, tmp2] <- FALSE

    stanctModel$manifesttype <- binaries

    if ( (indVarying == 'CINT') & (!(is.null(binaries.orig))) ) {
      tmp1 <- grep("_effect", colnames(stanctModel$pars)); tmp1
      tmp2 <- which(binaries.orig == 1); tmp2
      stanctModel$pars[(stanctModel$pars$matrix %in% 'CINT'), ][tmp2, tmp1] <- TRUE
    }

    fitStanctModel <- suppressMessages(ctsem::ctStanFit(
      fit=fit,
      datalong = datalong_all,
      ctstanmodel = stanctModel,
      sameInitialTimes=sameInitialTimes,
      savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
      stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
      iter=CoTiMAStanctArgs$iter,
      intoverstates=CoTiMAStanctArgs$intoverstates,
      binomial=CoTiMAStanctArgs$binomial,
      #fit=CoTiMAStanctArgs$fit,
      intoverpop=CoTiMAStanctArgs$intoverpop,
      stationary=CoTiMAStanctArgs$stationary,
      plot=CoTiMAStanctArgs$plot,
      #derrind=CoTiMAStanctArgs$derrind, # CHD deprecated, deleted Aug 2023
      optimize=CoTiMAStanctArgs$optimize,
      optimcontrol=CoTiMAStanctArgs$optimcontrol,
      nlcontrol=CoTiMAStanctArgs$nlcontrol,
      nopriors=CoTiMAStanctArgs$nopriors,
      priors=CoTiMAStanctArgs$priors, # added Aug 2023
      chains=CoTiMAStanctArgs$chains,
      forcerecompile=CoTiMAStanctArgs$forcerecompile,
      savescores=CoTiMAStanctArgs$savescores,
      gendata=CoTiMAStanctArgs$gendata,
      control=CoTiMAStanctArgs$control,
      verbose=CoTiMAStanctArgs$verbose,
      warmup=CoTiMAStanctArgs$warmup,
      cores=coresToUse,
      inits=inits))

    ### resample in parcels to avoid memory crash and speed up
    if (fit == TRUE) {
      if (!(is.null(CoTiMAStanctArgs$resample))) {
        fitStanctModel <- ctmaStanResample(ctmaFittedModel=fitStanctModel)
        #saveRDS(fitStanctModel, paste0(activeDirectory, "fitStanctModel.rds"))
        #fitStanctModel <- readRDS(paste0(activeDirectory, "fitStanctModel.rds"))
      }
      if (is.null(fitStanctModel$standata$priors)) fitStanctModel$standata$priors <- 0 # CHD added Sep 2023)
      fitStanctModel_summary <- summary(fitStanctModel, digits=2*digits, parmatrices=TRUE, residualcov=FALSE)
    } # end if (fit == TRUE)

    if (fit == FALSE) {
      print(paste0("#################################################################################"))
      print(paste0("#############  No model is fitted, only data and code are generated. ############"))
      print(paste0("#################################################################################"))

    }
  } # end if else (allInvModel)

  # Extract estimates & statistics
  if (fit == TRUE) { # CHD 16. Oct 2023 (end ~line 1534)
    # CHD 9.6.2023
    if ( (indVarying == TRUE) | (indVarying == 'cint') | (indVarying == 'CINT') ) {
      e <- ctsem::ctExtract(fitStanctModel)
      # CHD 12.6.2023
      model_popsd <- fitStanctModel_summary$popsd
      #if (indVaryingT0 == TRUE) {
      #if ( (indVaryingT0 == TRUE) & (T0meansParams[1] != 0) ) {
      if (dim(model_popsd)[1] != n.latent) {
        model_popsd <- fitStanctModel_summary$popsd
        model_popcov_m <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits)
        model_popcov_sd <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits)
        model_popcov_T <- round(ctsem::ctCollapse(e$popcov, 1, mean)/ctsem::ctCollapse(e$popcov, 1, stats::sd), digits)
        model_popcov_025 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025))
        model_popcov_50 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .50))
        model_popcov_975 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975))
        # convert to correlations and do the same (array to list then list to array)
        e$popcor <- lapply(seq(dim(e$popcov)[1]), function(x) e$popcov[x , ,])
        e$popcor <- lapply(e$popcor, stats::cov2cor)
        e$popcor <- array(unlist(e$popcor), dim=c(n.latent*2, n.latent*2, length(e$popcor)))
        model_popcor_m <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
        model_popcor_sd <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
        model_popcor_T <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
        model_popcor_025 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .025))
        model_popcor_50 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .50))
        model_popcor_975 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .975))
        #model_popcor <- stats::cov2cor(model_popcov_m)
      } else {
        model_popsd <- fitStanctModel_summary$popsd; model_popsd
        model_popcov_m <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits); model_popcov_m
        model_popcov_sd <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits)
        model_popcov_T <- round(ctsem::ctCollapse(e$popcov, 1, mean)/ctsem::ctCollapse(e$popcov, 1, stats::sd), digits)
        model_popcov_025 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .025))
        model_popcov_50 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .50))
        model_popcov_975 <- ctsem::ctCollapse(e$popcov, 1, function(x) stats::quantile(x, .975))
        # convert to correlations and do the same (array to list then list to array)
        e$popcor <- lapply(seq(dim(e$popcov)[1]), function(x) e$popcov[x , ,])
        e$popcor <- lapply(e$popcor, stats::cov2cor)
        e$popcor <- array(unlist(e$popcor), dim=c(n.latent   , n.latent  , length(e$popcor)))
        model_popcor_m <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
        model_popcor_sd <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
        model_popcor_T <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
        model_popcor_025 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .025))
        model_popcor_50 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .50))
        model_popcor_975 <- ctsem::ctCollapse(e$popcor, 3, function(x) stats::quantile(x, .975))
      }
    } else {
      model_popsd <- "no random effects estimated"
      model_popcov_m <- model_popcov_sd <- model_popcov_T <- model_popcov_025 <- model_popcov_50 <- model_popcov_975 <- "no random effects estimated"
      model_popcor_m <- model_popcor_sd <- model_popcor_T <- model_popcor_025 <- model_popcor_50 <- model_popcor_975 <- "no random effects estimated"
    }
    #model_popcov_m

    # account for changes in ctsem 3.4.1
    if ("matrix" %in% colnames(fitStanctModel_summary$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
    tmpMean <- grep("ean", colnames(fitStanctModel_summary$parmatrices)); tmpMean
    tmpSd <- tmpMean+1; tmpSd
    Tvalues <- fitStanctModel_summary$parmatrices[,tmpMean]/fitStanctModel_summary$parmatrices[,tmpSd]; Tvalues
    invariantDrift_Coeff <- cbind(fitStanctModel_summary$parmatrices, Tvalues); invariantDrift_Coeff
    invariantDrift_Coeff[, tmpMean:(dim(invariantDrift_Coeff)[2])] <- round(invariantDrift_Coeff[, tmpMean:(dim(invariantDrift_Coeff)[2])], digits); invariantDrift_Coeff

    # (create & ) re-label rownames
    {
      if (ctsem341) {
        tmp1 <- which(invariantDrift_Coeff[, "matrix"] == "DRIFT"); tmp1
        # replace row numbers by newly created names
        rownames(invariantDrift_Coeff) <- paste0(invariantDrift_Coeff[, c("matrix")], "_",
                                                 invariantDrift_Coeff[, c("row")], "_",
                                                 invariantDrift_Coeff[, c("col")])
      } else {
        tmp1 <- which(rownames(invariantDrift_Coeff) == "DRIFT")
      }
      rownames(invariantDrift_Coeff)[tmp1] <- driftFullNames; invariantDrift_Coeff
      # relabel if param is invariant
      tmp2 <- which(rownames(invariantDrift_Coeff) %in% invariantDriftNames); tmp2
      tmp3 <- paste0("DRIFT ", rownames(invariantDrift_Coeff)[tmp2] , " (invariant)"); tmp3
      rownames(invariantDrift_Coeff)[tmp2] <- tmp3; invariantDrift_Coeff
      # relabel if param is equal
      tmp2b <- which(rownames(invariantDrift_Coeff) %in% equalDriftNames); tmp2b
      if (length(tmp2b) > 0) {
        tmp3b <- paste0("DRIFT ", rownames(invariantDrift_Coeff)[tmp2b] , " (equal)"); tmp3b
        rownames(invariantDrift_Coeff)[tmp2b] <- tmp3b; invariantDrift_Coeff
      }

      tmp4 <- tmp1[which(!(tmp1 %in% tmp2))]; tmp4 # change to "DRIFT " for later extraction
      rownames(invariantDrift_Coeff)[tmp4] <- paste0("DRIFT ", driftFullNames[which(!(tmp1 %in% tmp2))]); invariantDrift_Coeff

      if (allInvModel == TRUE) {
        tmp5 <- (grep("DIFFUSIONcov", rownames(invariantDrift_Coeff))); tmp5
        tmp6 <- (grep("asymDIFFUSIONcov", rownames(invariantDrift_Coeff))); tmp6
        targetRows <- tmp5[which(!(tmp5 %in% tmp6))]; targetRows
        rownames(invariantDrift_Coeff)[targetRows] <- paste0(rownames(invariantDrift_Coeff)[targetRows], " (invariant)")

        targetRows <- (grep("T0cov", rownames(invariantDrift_Coeff))); targetRows
        rownames(invariantDrift_Coeff)[targetRows] <- paste0(rownames(invariantDrift_Coeff)[targetRows], " (invariant)")

        targetRows <- (grep("T0MEANS", rownames(invariantDrift_Coeff))); targetRows
        rownames(invariantDrift_Coeff)[targetRows] <- paste0(rownames(invariantDrift_Coeff)[targetRows], " (invariant)")

        targetRows <- (grep("MANIFESTMEANS", rownames(invariantDrift_Coeff))); targetRows
        rownames(invariantDrift_Coeff)[targetRows] <- paste0(rownames(invariantDrift_Coeff)[targetRows], " (invariant)")
      }

    }

    # extract fit
    invariantDrift_Minus2LogLikelihood  <- -2*fitStanctModel_summary$loglik; invariantDrift_Minus2LogLikelihood
    invariantDrift_estimatedParameters  <- fitStanctModel_summary$npars; invariantDrift_estimatedParameters
    invariantDrift_df <- "deprecated"

    # extract params
    {
      model_Drift_Coef <- invariantDrift_Coeff[(grep("DRIFT ", rownames(invariantDrift_Coeff))), tmpMean]; model_Drift_Coef
      tmp <- grep("DRIFT ", rownames(invariantDrift_Coeff)); tmp
      names(model_Drift_Coef) <- rownames(invariantDrift_Coeff)[tmp]; model_Drift_Coef


      model_Diffusion_Coef <- invariantDrift_Coeff[(rownames(invariantDrift_Coeff) == "DIFFUSIONcov"), tmpMean]; model_Diffusion_Coef
      if (length(model_Diffusion_Coef) < 1) {
        model_Diffusion_Coef <- invariantDrift_Coeff[invariantDrift_Coeff[, "matrix"] == "DIFFUSIONcov",  tmpMean]; model_Diffusion_Coef
      } else {
        model_Diffusion_Coef <- c(OpenMx::vech2full(model_Diffusion_Coef)); model_Diffusion_Coef
      }
      names(model_Diffusion_Coef) <- diffFullNames; model_Diffusion_Coef

      model_T0var_Coef <- invariantDrift_Coeff[(rownames(invariantDrift_Coeff) == "T0VAR"), 3]; model_T0var_Coef
      if (length(model_T0var_Coef) < 1) {
        model_T0var_Coef <- invariantDrift_Coeff[invariantDrift_Coeff[, "matrix"] == "T0cov",  tmpMean]; model_T0var_Coef
      } else {
        model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
      }
      names(model_T0var_Coef) <- driftFullNames; model_T0var_Coef
    }


    ## moderator effects
    modTI_Coeff <- NULL
    if (n.moderators > 0) {
      tmp1 <- c()
      for (i in modTIs) tmp1 <- c(tmp1, grep(i, rownames(fitStanctModel_summary$tipreds)))
      Tvalues <- fitStanctModel_summary$tipreds[tmp1, ][,6]; Tvalues
      modTI_Coeff <- round(cbind(fitStanctModel_summary$tipreds[tmp1, ], Tvalues), digits); modTI_Coeff

      # re-label
      if (!(is.null(mod.names))) {
        if (mod.type == "cont") {
          counter <- 0
          for (i in modTIs) {
            counter <- counter + 1
            targetNamePart <- paste0("tip_", modTIs[counter]); targetNamePart
            rownames(modTI_Coeff) <- sub(targetNamePart, paste0(mod.names[counter], "_on_"), rownames(modTI_Coeff))
          }
        }

        if (mod.type == "cat") {
          counter <- 0
          modNameCounter <- 1
          for (j in modTIs) {
            if (n.moderators == 1) unique.mod.tmp <- unique.mod else unique.mod.tmp <- unique.mod[[counter+1]]
            if (!(is.null(catsToCompare))) {
              origCats <- c(unique.mod.tmp[1:2] + moderatorGroupsOffset, unique.mod.tmp[-c(1:2)]); origCats
            } else {
              origCats <- unique.mod.tmp
            }
            for (i in 1:(length(unique.mod.tmp)-1)) {
              counter <- counter + 1; counter
              current.mod.names <- mod.names[modNameCounter]; current.mod.names
              targetNamePart <- paste0("tip_", modTIs[i]); targetNamePart
              tmp1 <- grep(targetNamePart, rownames(modTI_Coeff)); tmp1
              rownames(modTI_Coeff) <- sub(targetNamePart, paste0(origCats[counter+1], "  (category value) of ", mod.names[modNameCounter], "_on"), rownames(modTI_Coeff))
            }
            counter <- 0
            modNameCounter <- modNameCounter + 1
          }
        }
      }
      # eliminate z
      modTI_Coeff[, "z"] <- NULL; modTI_Coeff
    }

    ## cluster effects
    if (!(is.null(cluster))) {
      cluster.specific.effect <- matrix(NA, clusCounter, n.latent^2)
      tmp <- grep("DRIFT ", rownames(invariantDrift_Coeff)); tmp
      colnames(cluster.specific.effect) <- rownames(invariantDrift_Coeff)[tmp]; cluster.specific.effect
      rownames(cluster.specific.effect) <- paste0("Cluster No. ", seq(1,clusCounter, 1)); cluster.specific.effect
      tmp1 <- c()
      for (i in clusTIs) tmp1 <- c(tmp1, (grep(i, rownames(fitStanctModel_summary$tipreds))))
      Tvalues <- fitStanctModel_summary$tipreds[tmp1, ][,6]; Tvalues
      clusTI_Coeff <- round(cbind(fitStanctModel_summary$tipreds[tmp1, ], Tvalues), digits); clusTI_Coeff
      # re-label
      targetCluster <- which(table(cluster) > 1); targetCluster
      for (i in 1:clusCounter) {
        #i <- 2
        targetNamePart <- paste0("tip_", clusTIs[i]); targetNamePart
        rownames(clusTI_Coeff) <- sub(targetNamePart, paste0("Cluster_", targetCluster[i], "_on_"), rownames(clusTI_Coeff))
        for (j in 1:length(driftNames)) {
          if (driftNames[j] != 0) {
            tmp0 <- driftNames[j]; tmp0
            tmp0 <- gsub(" \\(invariant\\)", "", tmp0); tmp0
            tmp1 <- grep(tmp0, rownames((clusTI_Coeff))); tmp1
            tmp2 <- grep(tmp0, names((model_Drift_Coef))); tmp2
            # CHD changed 7.6.2023
            # cluster.specific.effect[i,j] <- round(model_Drift_Coef[tmp2] + clusTI_Coeff[tmp1, 1] * cluster.weights[i, 2], digits)
            cluster.specific.effect[i,j] <- round(model_Drift_Coef[tmp2] + clusTI_Coeff[tmp1[i], 1] * cluster.weights[i, 2], digits)
          }
        }
      }
    } else {
      clusTI_Coeff <- NULL
    }

    if (ctsem341) invariantDrift_Coeff[, "matrix"] <- NULL

    ### Numerically compute Optimal Time lag sensu Dormann & Griffin (2015)
    if (ctsem341) {
      driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=TRUE); driftMatrix # byrow set because order is different compared to mx model
    } else {
      driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=FALSE); driftMatrix # byrow set because order is different compared to mx model
    }

    if (length(invariantDriftNames) == length(driftNames)) {
      OTL <- function(timeRange) {
        OpenMx::expm(tmpDriftMatrix * timeRange)[targetRow, targetCol]}
      # use original time scale
      if (!(is.null(scaleTime))) {
        tmpDriftMatrix <- driftMatrix * scaleTime
      } else {
        tmpDriftMatrix <- driftMatrix
      }
      # loop through all cross effects
      tmp1 <- 0
      if (0 %in% usedTimeRange) tmp1 <- 1
      optimalCrossLag <- matrix(NA, n.latent, n.latent)
      maxCrossEffect <- matrix(NA, n.latent, n.latent)
      for (j in 1:n.latent) {
        for (h in 1:n.latent) {
          if (j != h) {
            targetRow <- j
            targetCol <- h
            if (tmpDriftMatrix[j, h] != 0) { # an effect that is zero has no optimal lag
              targetParameters <- sapply(usedTimeRange, OTL); targetParameters
              maxCrossEffect[j,h] <- max(abs(targetParameters))[1]; maxCrossEffect[j,h]
              optimalCrossLag[j,h] <- which(abs(targetParameters)==maxCrossEffect[j,h])[1]*1 - tmp1 # first targetParam is calculated for lag=0
            } else {
              optimalCrossLag[j,h] <- NA
            }
          }
        }
      }
      maxCrossEffect <- round(maxCrossEffect, digits); maxCrossEffect
    } else {
      optimalCrossLag <- "Drift Matrix is only partially invariant. (Generalizable) optimal intervall cannot be computed."
      maxCrossEffect <- "Drift Matrix is only partially invariant. (Generalizable) Largest effects cannot be computed."
    }

    # CHD 12.7.23
    #if (!(is.null(scaleTime))) {
    if ( (!(is.null(scaleTime))) & (length(invariantDriftNames) == length(driftNames)) ) {
      optimalCrossLag_scaledTime <- optimalCrossLag * CoTiMAStanctArgs$scaleTime
    } else {
      optimalCrossLag_scaledTime <- optimalCrossLag
    }

    #######################################################################################################################

    #end.time <- Sys.time()
    #time.taken <- end.time - start.time
    #st <- paste0("Computation started at: ", start.time); st
    #et <- paste0("Computation ended at: ", end.time); et
    #tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt

    meanDeltas <- mean(ctmaInitFit$statisticsList$allDeltas, na.rm=TRUE); meanDeltas
    largeDelta <- which(ctmaInitFit$statisticsList$allDeltas >= meanDeltas); largeDelta
    tmp1 <- table(ctmaInitFit$statisticsList$allDeltas[largeDelta]); tmp1
    tmp2 <- which(names(tmp1) == (max(as.numeric(names(tmp1))))); tmp2
    suggestedScaleTime <- as.numeric(names(tmp1[tmp2])); suggestedScaleTime
    message <- c()
    #if (meanDeltas > 3) {
    # CHD AUG 2023
    if ((meanDeltas * CoTiMAStanctArgs$scaleTime) > 3) { #xxx
      tmp2 <- paste0("Mean time interval was ", meanDeltas, "."); tmp2
      tmp3 <- paste0("scaleTime=1/", suggestedScaleTime); tmp3
      tmp4 <- paste0("It is recommended to fit the model again using the arguments ", tmp3, " and customPar=FALSE. "); tmp4
      message <- paste(tmp2, tmp4, "If the model fit (-2ll) is better (lower), continue using, e.g.,", tmp3, "in all subsequent models.", collapse="\n"); message
    }


    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

    if (ctsem341) { # new 8.7.2022
      tmp1 <- grep("CINT", fitStanctModel_summary$parmatrices[,"matrix"]); tmp1
      tmp2 <- grep("asym", fitStanctModel_summary$parmatrices[,"matrix"]); tmp2
      tmp3 <- grep("dt", fitStanctModel_summary$parmatrices[,"matrix"]); tmp3
      tmp4 <- tmp1[(tmp1 %in% c(tmp2, tmp3)) == FALSE]; tmp4
      model_Cint_Coef <- fitStanctModel_summary$parmatrices[tmp4, 4]; model_Cint_Coef
    } else {
      tmp1 <- grep("CINT", rownames(fitStanctModel_summary$parmatrices)); tmp1
      tmp2 <- grep("asym", rownames(fitStanctModel_summary$parmatrices)); tmp2
      tmp3 <- grep("dt", rownames(fitStanctModel_summary$parmatrices)); tmp3
      tmp4 <- tmp1[(tmp1 %in% c(tmp2, tmp3)) == FALSE]; tmp4
      model_Cint_Coef <- fitStanctModel_summary$parmatrices[tmp4, 3]; model_Cint_Coef
    }

    if (!(is.null(cluster))) {
      clus.effects=list(effects=clusTI_Coeff, weights=cluster.weights, sizes=cluster.sizes,
                        cluster.specific.effect=cluster.specific.effect, note=cluster.note)
    } else {
      clus.effects <- NULL
    }

    if (is.null(primaryStudyList)) primaryStudies <- ctmaInitFit$primaryStudyList else primaryStudies <- primaryStudyList

    if (is.null(scaleTime)) scaleTime2 <- 1 else scaleTime2 <- scaleTime

    if (!(is.null(scaleTime))) {
      model_Drift_Coef_original_time_scale <- model_Drift_Coef * scaleTime
      model_Diffusion_Coef_original_time_scale <- model_Diffusion_Coef * scaleTime
      if (!(is.null(modTI_Coeff))) {
        modTI_Coeff_original_time_scale <- modTI_Coeff * scaleTime
        modTI_Coeff_original_time_scale[, "Tvalues"] <- modTI_Coeff[, "Tvalues"]
      } else {
        modTI_Coeff_original_time_scale <- NULL
      }
      if (!(is.null(clusTI_Coeff))) {
        clusTI_Coeff_original_time_scale <- clusTI_Coeff * scaleTime
        clusTI_Coeff_original_time_scale[, "Tvalues"] <- clusTI_Coeff[, "Tvalues"]
      } else {
        clusTI_Coeff_original_time_scale <- NULL
      }
      tmp1<- invariantDrift_Coeff
      tmp2 <- grep("toV", rownames(tmp1))
      tmp3 <- grep("DIFFUSIONcov", rownames(tmp1))
      tmp4 <- grep("asym", rownames(tmp1))
      tmp3 <- tmp3[!(tmp3%in% tmp4)]
      tmp1 <- tmp1[c(tmp2, tmp3),]
      tmp1[, c("Mean", "sd", "2.5%", "50%", "97.5%")] <- tmp1[, c("Mean", "sd", "2.5%", "50%", "97.5%")] * scaleTime
      estimates_original_time_scale <- tmp1
      mod_effects_original_time_scale <- modTI_Coeff_original_time_scale
      if (!(is.null(clusTI_Coeff))) {
        clus_effects_original_time_scale <- round(clusTI_Coeff_original_time_scale, digits)
      } else {
        clus_effects_original_time_scale <- NULL
      }
    } else {
      model_Drift_Coef_original_time_scale <- round(model_Drift_Coef, digits); model_Drift_Coef_original_time_scale
      model_Diffusion_Coef_original_time_scale <- round(model_Diffusion_Coef, digits); model_Diffusion_Coef_original_time_scale
      if (!(is.null(modTI_Coeff))) { # new 9.7.20222
        modTI_Coeff_original_time_scale <- round(modTI_Coeff, digits)
        mod_effects_original_time_scale <- round(modTI_Coeff, digits)  # doubled (why?)
      } else {
        modTI_Coeff_original_time_scale <- NULL
        mod_effects_original_time_scale <- NULL
      }

      { # new 8.7.2022
        tmp1<- invariantDrift_Coeff
        tmp2 <- grep("toV", rownames(tmp1))
        tmp3 <- grep("DIFFUSIONcov", rownames(tmp1))
        tmp4 <- grep("asym", rownames(tmp1))
        tmp3 <- tmp3[!(tmp3%in% tmp4)]
        tmp1 <- tmp1[c(tmp2, tmp3),]
        tmp1[, c("Mean", "sd", "2.5%", "50%", "97.5%")] <- tmp1[, c("Mean", "sd", "2.5%", "50%", "97.5%")] * 1
        estimates_original_time_scale <- round(tmp1, digits)
      }
      clus_effects_original_time_scale <- NULL
      if (!(is.null(clusTI_Coeff))) {
        clusTI_Coeff_original_time_scale <- round(clusTI_Coeff, digits)
      } else {
        clusTI_Coeff_original_time_scale <- NULL
      }
    }

    results <- list(#activeDirectory=activeDirectory,
      plot.type="drift",  model.type="stanct",
      n.studies=1,
      n.latent=n.latent,
      n.moderators=length(mod.number),
      studyList=ctmaInitFit$studyList,
      studyFitList=fitStanctModel,
      data=datalong_all,
      statisticsList=ctmaInitFit$statisticsList,
      argumentList=list(ctmaInitFit=deparse(substitute(ctmaInitFit)),
                        primaryStudyList=primaryStudies,
                        cluster=cluster,
                        activeDirectory=activeDirectory,
                        activateRPB=activateRPB,
                        digits=digits,
                        drift=drift,
                        invariantDrift=invariantDrift,
                        moderatedDrift=moderatedDrift,
                        equalDrift=equalDrift,
                        ind.mod.names=ind.mod.names,
                        ind.mod.number=ind.mod.number,
                        ind.mod.type=ind.mod.type,
                        mod.number=mod.number,
                        mod.type=mod.type,
                        mod.names=mod.names,
                        indVarying=indVarying,
                        coresToUse=coresToUse,
                        scaleTI=scaleTI,
                        scaleMod=scaleMod,
                        transfMod=transfMod,
                        sameInitialTimes=sameInitialTimes,
                        scaleClus=scaleClus,
                        scaleTime=scaleTime,
                        optimize=optimize,
                        nopriors=nopriors,
                        finishsamples=finishsamples,
                        iter=iter,
                        chains=chains,
                        verbose=verbose,
                        allInvModel=allInvModel,
                        customPar=customPar,
                        inits=inits,
                        modsToCompare=modsToCompare,
                        catsToCompare=catsToCompare,
                        driftsToCompare=driftsToCompare,
                        useSampleFraction=useSampleFraction,
                        T0means=T0means,
                        manifestMeans=manifestMeans,
                        CoTiMAStanctArgs=CoTiMAStanctArgs),
      modelResults=list(DRIFToriginal_time_scale=model_Drift_Coef_original_time_scale,
                        DIFFUSIONoriginal_time_scale=model_Diffusion_Coef_original_time_scale,
                        T0VAR=model_T0var_Coef,
                        CINT=model_Cint_Coef,
                        MOD=modTI_Coeff_original_time_scale,
                        CLUS=clusTI_Coeff_original_time_scale,
                        DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef),
      parameterNames=ctmaInitFit$parameterNames,
      CoTiMAStanctArgs=CoTiMAStanctArgs,
      invariantDrift=invariantDrift,
      summary=list(model="Model name not specified", # paste(invariantDrift, "unequal but invariant across samples", collapse=" "),
                   scaledTime=scaleTime2,
                   estimates=invariantDrift_Coeff,
                   # changed 17. Aug. 2022
                   randomEffects=list(popsd=model_popsd,
                                      popcov_mean=model_popcov_m, model_popcov_sd=model_popcov_sd,
                                      model_popcov_T=model_popcov_T, model_popcov_025=model_popcov_025,
                                      model_popcov_50=model_popcov_50, model_popcov_975=model_popcov_975,
                                      popcor_mean=model_popcor_m, model_popcor_sd=model_popcor_sd,
                                      model_popcor_T=model_popcor_T, model_popcor_025=model_popcor_025,
                                      model_popcor_50=model_popcor_50, model_popcor_975=model_popcor_975),
                   minus2ll= invariantDrift_Minus2LogLikelihood,
                   n.parameters = invariantDrift_estimatedParameters,
                   #optimalLagInfo = "Optimal lag and effect was calculated for original time scale (i.e., ignoring possible scaleTime argument).",
                   opt.lag.orig.time = optimalCrossLag,
                   opt.lag.scaled.time = optimalCrossLag_scaledTime,
                   max.effects = maxCrossEffect,
                   clus.effects=clus.effects,
                   mod.effects=modTI_Coeff,
                   message=message,
                   estimates_original_time_scale =estimates_original_time_scale,
                   mod_effects_original_time_scale=mod_effects_original_time_scale,
                   clus_effects_original_time_scale=clus_effects_original_time_scale)
      # excel workbook is added later
    )

  } #end if (fit == TRUE) { # CHD 16. Oct 2023 (end ~line 1115)


  if (fit == FALSE) {
    results <- list(summary=c("No model was fitted, only data and code were generated. See $data & $ctModel section."),
                    data = datalong_all,
                    ctModel = fitStanctModel$ctstanmodelbase)
  }

  class(results) <- "CoTiMAFit"

  ### prepare Excel Workbook with several sheets ################################################################
  if (fit == TRUE) {
    {
      wb <- openxlsx::createWorkbook()
      sheet1 <- openxlsx::addWorksheet(wb, sheetName="model")
      sheet2 <- openxlsx::addWorksheet(wb, sheetName="modelResults")
      sheet3 <- openxlsx::addWorksheet(wb, sheetName="estimates")
      sheet7 <- openxlsx::addWorksheet(wb, sheetName="mod.effects")
      sheet4 <- openxlsx::addWorksheet(wb, sheetName="clus.effects")
      sheet5 <- openxlsx::addWorksheet(wb, sheetName="randomEffects")
      sheet6 <- openxlsx::addWorksheet(wb, sheetName="stats")
      openxlsx::writeData(wb, sheet1, results$summary$model)

      ### modelResults
      # DRIFT
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, matrix(driftNames, nrow=1), colNames = FALSE)
      startCol <- 2; startCol
      startRow <- 2; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          t(c(t(matrix(unlist(results$modelResults$DRIFT),
                                       nrow=n.latent, ncol=n.latent)))))
      startCol <- 1; startCol
      startRow <- 2; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix("Aggregated Drift Effects", ncol=1))
      # DIFFUSION
      offset <-  1 + 1
      startCol <- 2; startCol
      startRow <- 2 + offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow,
                          matrix(names(results$modelResults$DIFFUSION), nrow=1), colNames = FALSE)
      startCol <- 2; startCol
      startRow <- 2 + offset + 1# offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(results$modelResults$DIFFUSION,
                                 nrow=1, byrow=TRUE))
      results$modelResults$DIFFUSION
      startCol <- 1; startCol
      startRow <- 2 + offset + 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix("Diffusion Population Means (not aggregated)", nrow=1, byrow=TRUE))
      # T0Var
      offset <-  offset + 1 + 2
      startCol <- 2; startCol
      startRow <- 2 + offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow,
                          matrix(names(results$modelResults$T0VAR), nrow=1), colNames = FALSE)
      startCol <- 2; startCol
      startRow <- 2 + offset + 1# offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(results$modelResults$T0VAR), nrow=1, byrow=TRUE))
      startCol <- 1; startCol
      startRow <- 2 + offset + 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix("T0VAR Population Mean (not aggregated)", nrow=1, byrow=TRUE))
      ### estimates
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet3, startCol=startCol, startRow = startRow + 1,
                          rownames(results$summary$estimates), colNames = FALSE)
      openxlsx::writeData(wb, sheet3, startCol=startCol+1, startRow = startRow, results$summary$estimates,
                          colNames = TRUE)
      ### moderator Effects
      startCol <- 1; startCol
      startRow <- 1; startRow
      results$summary$mod.effects
      openxlsx::writeData(wb, sheet7, startCol=startCol, startRow = startRow + 1,
                          results$summary$mod.effects, colNames = TRUE, rowNames = TRUE)

      ### cluster Effects
      if (!(is.null(cluster))) {
        startCol <- 2; startCol
        startRow <- 1; startRow
        tmp1 <- dim(results$summary$clus.effects$effects); tmp1
        openxlsx::writeData(wb, sheet4, startCol=startCol, startRow = startRow + 1,
                            rownames(results$summary$clus.effects$effects), colNames = FALSE)
        openxlsx::writeData(wb, sheet4, startCol=startCol+1, startRow = startRow, results$summary$clus.effects$effects,
                            colNames = TRUE)
        openxlsx::writeData(wb, sheet4, startCol=startCol+1, startRow = startRow+tmp1[1] + 2,
                            results$summary$clus.effects$cluster.specific.effect,
                            colNames = TRUE, rowNames = TRUE)
      }
      ### random Effects
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet5, startCol=startCol, startRow = startRow, results$summary$randomEffects, colNames = FALSE)
      ### stats
      startCol <- 2; startCol
      startRow <- 1; startRow
      tmp <- cbind("-2ll = ", results$summary$minus2ll, "Number of Parameters = ", results$summary$n.parameters)
      openxlsx::writeData(wb, sheet6, startCol=startCol, startRow = startRow, t(tmp), colNames = FALSE)
    }

    results$excelSheets <- wb
  } #end if (fit == TRUE) {

  invisible(results)

} ### END function definition
