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
#' @param cint default 'auto' (= 0). Are set free if random intercepts model with varying cints is requested (by indVarying='cint')
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
#' @param indVaryingT0 deprecated. Automatically set to NULL.
#' @param inits vector of start values
#' @param invariantDrift  drift labels for drift effects that are set invariant across primary studies (default = all drift effects).
#' @param iter number of iterations (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param manifestMeans default = 0. Are automatically set free if indVarying is set to TRUE. Can be assigned labels to estimate them freely.
#' @param manifestVars define the error variances (default = 0) of the manifests with a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param mod.names vector of names for moderators used in output
#' @param mod.number which in the vector of moderator values shall be used (e.g., 2 for a single moderator or 1:3 for 3 moderators simultaneously)
#' @param mod.type 'cont' or 'cat' (mixing them in a single model not yet possible)
#' @param moderatedDrift labels for drift effects that are moderated (default = all drift effects)
#' @param modsToCompare when performing contrasts for categorical moderators, the moderator numbers (position in mod.number) that is used
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param primaryStudyList  could be a list of primary studies compiled with \code{\link{ctmaPrep}} that defines the subset of studies in ctmaInitFit that should actually be used
#' @param priors if FALSE, any priors are disabled – sometimes desirable for optimization
#' @param randomIntercepts (default = FALSE) Experimental. Overrides ctsem's default mode for modelling indVarying cints.
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
#' @param WEC (default = FALSE) Experimental. Uses weighted effect coding of TIpred representing the dummies of the primary studies. Returns drift matrices for all primary studies.
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
    #nopriors=TRUE,
    optimize=TRUE,
    primaryStudyList=NULL,
    priors=FALSE,
    randomIntercepts=FALSE,
    sameInitialTimes=FALSE,
    scaleClus=TRUE,
    scaleMod=TRUE,
    scaleTI=TRUE,
    scaleTime=NULL,
    T0means=0,
    T0var='auto',
    transfMod=NULL,
    useSampleFraction=NULL,
    verbose=0,
    WEC=FALSE
)
{  # begin function definition (until end of file)

  {
    ctmaInitFitName <- deparse(substitute(ctmaInitFit))

    if (is.null(scaleTime)) scaleTime <- 1

    # indVaryingT0 previously allowed the T0cov to vary across primaries, the cints to covary randomly for the entire sample, and
    # the cints NOT to covary with the latents at T0. The next three line prvent this, but can be deleted to make it work again.
    if (!(is.null(indVaryingT0))) {
      Msg <- "The argument \"indVaryingT0\" was deprecated and set to NULL. Try \"indVarying = TRUE\" or \"randomIntercepts = TRUE\" instead.\n"
      message(Msg)
    }

    { # adaptations to account for new arguments introduces
      if (is.null(T0var)) T0var <- 'auto'
      if (is.null(cint)) cint <- 0
      if (is.null(fit)) fit <- TRUE
      if (is.null(WEC)) WEC <- FALSE
      if (is.null(CoTiMAStanctArgs)) CoTiMAStanctArgs <- CoTiMA::CoTiMAStanctArgs
    }

    {
      if ( (indVarying == "CINT") | (indVarying == "Cint") | (indVarying == "cint")) indVarying <- "CINT"
      if ( (indVarying == "MANIFEST") | (indVarying == "Manifest") | (indVarying == "manifest")) indVarying <- TRUE
      #
      if (is.null(randomIntercepts)) randomIntercepts <- FALSE
      if ( (randomIntercepts == 'CINT') | (randomIntercepts == 'cint')  | (randomIntercepts == 'Cint')) randomIntercepts <- TRUE
      if ( (randomIntercepts == "Manifest") | (randomIntercepts == "manifest") | (randomIntercepts == "MANIFEST")) randomIntercepts <- "MANIFEST"
      if ( (randomIntercepts == "MANIFEST") | (randomIntercepts == TRUE) ) {
        indVarying <- FALSE
        indVaryingT0 <- NULL
      }
      randomInterceptsSettings <- randomIntercepts
    }

    # adapt display of information during model fit
    if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMAStanctArgs$verbose}

    # check if fit object is specified
    if (is.null(ctmaInitFit)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "\nA fitted CoTiMA object has to be supplied to plot something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    # CHD 18.12.2023
    if (WEC == TRUE) {
      Msg <- "The argument WEC=TRUE was used. I assume you want to mimic ctmaInit with all drift effects varying across primary studies (using weighted effect coding).
    Therefore, I set the argument invariantDrift=\'none\', so that all drift effects vary across primary studies."
      message(Msg)
      invariantDrift[1] <- 'none'
      scaleTI <- FALSE
    }

    if (!(is.null(invariantDrift))) { # added 12.7.2023
      if ( ((invariantDrift[1] == "none") | (invariantDrift[1] == "None") | (invariantDrift[1] == "NONE"))  & (scaleTI == TRUE) ) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
        Msg <- "The argument invariantDrift=\'none\' was used. I assume you want to mimic ctmaInit with all drift effects varying across primary studis.
    Therefore, I set the argument scaleTI=FALSE."
        message(Msg)
        scaleTI <- FALSE
      }
    }

    { # set fitting params

      #if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
      #if (!(is.null(scaleClus))) CoTiMAStanctArgs$scaleClus <- scaleClus
      #if (!(is.null(scaleMod))) CoTiMAStanctArgs$scaleMod <- scaleMod
      #if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
      if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
      #if  (!(is.null(nopriors))) CoTiMAStanctArgs$nopriors <- nopriors # changed Aug 2023
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
        Msg <- "The arguments scaleMod (default = TRUE) and transfMod cannot be used in combination. ScalMod was set to FALSE."
        message(Msg)
      }
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
        # CHD added 25.1.2024
        ctmaTempFit$ind.mod.List[[i]] <- NULL
        ctmaTempFit$modelResults$CINT[[i]] <- NULL
        ctmaTempFit$modelResults$DRIFToriginal_time_scale[[i]] <- NULL
        ctmaTempFit$modelResults$DIFFUSIONoriginal_time_scale[[i]] <- NULL
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
    # CHD 25.1.2024
    ctmaTempFit$summary$model <- "not specified"
    #ctmaTempFit$summary$model <- "Moderator Model (for details see model summary)"
    tmpStudyNumber <- as.numeric(gsub("Study No ", "", rownames(ctmaTempFit$summary$estimates))); tmpStudyNumber
    targetRows <- which(tmpStudyNumber %in% targetStudyNumbers); targetRows; #length(targetRows)
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

    # CHD 15.1.2024
    exclude <- which(!unlist(ctmaInitFit$primaryStudyList$studyNumbers) %in% targetStudyNumbers); exclude # study position
    if (length(exclude) > 0) {
      for (i in sort(exclude, decreasing = TRUE)) {
        ctmaInitFit$ind.mod.List[[i]] <- NULL
      }
    }

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

  {
    #start.time <- Sys.time(); start.time

    {
      if (is.null(ctmaInitFit$n.latent)) {
        n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
      } else {
        n.latent <- ctmaInitFit$n.latent
      }
      if (!(is.null(ctmaInitFit$n.manifest))) n.manifest <- ctmaInitFit$n.manifest else n.manifest <- n.latent
      if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory

      if (!(is.null(drift))) {
        if (!(is.matrix(drift))) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
          ErrorMsg <- paste0("The argument \"drift = c(", paste(drift, collapse=", "), ")\" was provided, but drift should be a ", n.latent, " x ", n.latent, " matrix.")
          stop(ErrorMsg)
        }
      }

      binaries.orig <- binaries
      if (is.null(binaries)) {
        binaries <- rep(0, max(n.manifest, n.latent))
      }
      if (length(binaries) != max(n.manifest, n.latent)) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        ErrorMsg <- "\nThe number of binaries provided is incorrect! \nGood luck for the next try!"
        stop(ErrorMsg)
      }

      if ( (!(is.null(binaries.orig))) & (indVarying != 'CINT') & (!(is.null(binaries.orig))) & (!all(binaries == 0)) ) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        ErrorMsg <- "\nYou specified binary variables. You also need to specify \"indVarying=\'CINT\'\". \nGood luck for the next try!"
        stop(ErrorMsg)
      }

      if ( (!(is.null(binaries.orig))) & !all(binaries == 0) ) message("Effects of binaries on cints not implemented yet.")

      n.studies <- unlist(ctmaInitFit$n.studies); n.studies
      allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
      maxTpoints <- max(allTpoints); maxTpoints
      allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
      maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
      usedTimeRange <- seq(0, 3*maxDelta, 1); usedTimeRange # new 8.7.2022
      lambda <- ctmaInitFit$statisticsList$lambda
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
      if (n.ind.moderators > 0) {
        mod.number <- NULL
        mod.type=ind.mod.type
        mod.names <- NULL
      }

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
        # CHD 15.1.2024
        validStudies <- unlist(lapply(tmpMods, function(x) !is.null(x))); validStudies
        if (any(validStudies == FALSE)) {
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
          ErrorMsg <- "\nIndividual level moderation model requested. At least one study has fewer individual level moderators (NULL?) in the
        raw data file than the number specified via ind.mod.number . Check your data, redo ctmaPrep and ctmaInit again. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
        #
        if (is.data.frame(tmpMods[[1]])) {
          tmpMods <- do.call(rbind, tmpMods) # rbinds all data frames
        } else {
          if (!(is.matrix(tmpMods[[1]]))) tmpMods <- lapply(tmpMods, function(x) matrix(x, ncol=1))
          if (is.list(tmpMods)) {
            tmpModstmp <- matrix(unlist(tmpMods[[1]]), ncol=ncol(tmpMods[[1]]))
            for (t in 2:length(tmpMods)) tmpModstmp <- rbind(tmpModstmp, matrix(unlist(tmpMods[[t]]), ncol=ncol(tmpMods[[t]])))
            tmpMods <- as.matrix(tmpModstmp, ncol=1)
          }
          # the following might be skipped (until skip end)
          tmp1 <- unlist(lapply(tmpMods, function(x) {
            if (is.null(dim(x))) return(1) else return(ncol(x))
          } )); tmp1
          if (length(ind.mod.number) > min(tmp1)) {
            if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
            ErrorMsg <- "\nIndividual level moderation model requested. At least one study has fewer individual level moderators in the raw data file than the number
          specified via ind.mod.number . Check your data, redo ctmaPrep and ctmaInit again. \nGood luck for the next try!"
            stop(ErrorMsg)
          } # skip end
        }
        #
        if (!(is.null(primaryStudyList))) {
          if (length(primaryStudyList$deltas) > n.studies) { # CHD changed Aug 2023
            if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
            ErrorMsg <- "\nYou provided a list to the argument primaryStudyList. The number of studies in this list is not identical to the number of
          studies fitted with ctmaInit as provided in the argument ctmaInitFit. \nGood luck for the next try!"
            stop(ErrorMsg)
          }
        }
        currentModerators <- tmpMods
        colnames(currentModerators) <- paste0("mod", 1:dim(currentModerators)[2])
      }
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
        }
        #
        for (c1 in 1:dim(moderatorGroups)[2]) {
          if (c1 %in% modsToCompare) {
            maxAbsModeratorGroups <- max(abs(moderatorGroups[, c1])) ; maxAbsModeratorGroups
            minAbsModeratorGroups <- min(abs(moderatorGroups[, c1])) ; minAbsModeratorGroups
            moderatorGroupsOffset <- maxAbsModeratorGroups + minAbsModeratorGroups; moderatorGroupsOffset
            for (c2 in catsToCompare) {
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
    if (is.null(dim(groups))) groups <- matrix(groups,  ncol=1)
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
      if (scaleTI == TRUE) dataTmp[ , ncol(dataTmp)] <- scale(dataTmp[ , ncol(dataTmp)])
    }
    targetCols <- which(colnames(dataTmp) == "groups"); targetCols
    dataTmp <- dataTmp[ ,-targetCols]

    # use weighted effect coding
    if (WEC == TRUE) {
      targetCols <- grep("TI", colnames(dataTmp)); targetCols
      targetRows <- which(apply(dataTmp[, targetCols], 1, sum) == 0); targetRows
      targetN <- length(targetRows); targetN
      for (i in 1:length(targetCols)) {
        tmp1 <- sum(dataTmp[, targetCols[i]])
        dataTmp[targetRows, targetCols[i]] <- -tmp1/targetN
      }
    }

    # make TI out of moderators
    modTIstartNum <- n.studies; modTIstartNum
    if (n.moderators > 0) {
      # make TI predictors as replacements for each moderator
      if (mod.type=="cont") {
        tmp1 <- paste0("mod", 1:n.moderators); tmp1
        if (length(tmp1) == 1) tmp <- matrix(dataTmp[ , tmp1], ncol=length(tmp1)) else tmp <- dataTmp[ , tmp1]
        if (scaleMod == TRUE) tmp[ , 1:ncol(tmp)] <- scale(tmp[ , 1:ncol(tmp)])
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

        if (scaleMod == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)], scale=FALSE)
        currentStartNumber <- modTIstartNum; currentStartNumber
        currentEndNumber <- currentStartNumber + ncol(tmpTI)-1; currentEndNumber
        colnames(tmpTI) <- paste0("TI", currentStartNumber:currentEndNumber)

        dataTmp <- cbind(dataTmp, tmpTI); dim(dataTmp)
        dataTmp <- dataTmp[ ,-grep("mod", colnames(dataTmp))]
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
      #
      n.var <- max(n.manifest, n.latent); n.var
    }

    # possible subsample selection
    if (!(is.null(useSampleFraction))) {
      N <- dim(dataTmp)[1]; N
      stepwidth <- 100/useSampleFraction
      targetCases <- round(seq(1, N, stepwidth), 0); targetCases
      dataTmp <- dataTmp[targetCases, ]
    }

    # make long data format
    {
      dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.var, n.TIpred = (n.studies-1+tmp1),
                                      manifestNames=manifestNames)
      dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
      dataTmp3[, "time"] <- dataTmp3[, "time"] * scaleTime
    }

    # eliminate rows where ALL latents are NA
    {
      if (n.manifest > n.latent) namePart <- "y" else namePart <- "V"
      dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0(namePart, 1:n.var)], 1, function(x) sum(is.na(x)) != n.var ), ]
      datalong_all <- dataTmp3
    }

    datalong_all <- as.data.frame(datalong_all)

    # TI-identifiers for groups, moderators, and clusters
    {
      groupTIs <- paste0("TI", 1:(length(unique(groups))-1)); groupTIs
      tmp1 <- (length(unique(groups))); tmp1
      if (n.moderators > 0) modTIs <- paste0("TI", tmp1:(tmp1+n.all.moderators-1))
      tmp1 <- (tmp1+n.all.moderators); tmp1
      if (!(is.null(cluster))) clusTIs <- paste0("TI", tmp1:(tmp1+clusCounter-1))
    }

  }

  #######################################################################################################################
  ############################################### Define Parameter Names ################################################
  #######################################################################################################################

  {
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
    manifestVarsParams <- namesAndParams$manifestVarsParams
    # CHD 9.6.2023
    T0VARParams <- T0var
    CINTParams <- cint

    if (is.null(invariantDriftNames)) invariantDriftNames <- driftNames
  }

  #######################################################################################################################
  ######################### All-invariant Model (used for calculation of statistical power) #############################
  #######################################################################################################################

  if (allInvModel == TRUE) {
    allInvModelFit <- ctmaAllInvFit(ctmaInitFit=ctmaInitFit,
                                    activeDirectory=activeDirectory,
                                    activateRPB=activateRPB,
                                    digits=digits,
                                    drift=drift,
                                    coresToUse=coresToUse,
                                    scaleTime=scaleTime,
                                    optimize=optimize,
                                    priors=priors,
                                    finishsamples=finishsamples,
                                    iter=iter,
                                    chains=chains,
                                    verbose=verbose,
                                    indVarying = indVarying,
                                    indVaryingT0 = indVaryingT0,
                                    customPar = customPar)
    stanctModel <- allInvModelFit$ctModel
    fitStanctModel <- allInvModelFit$studyFitList[[1]]
    fitStanctModel_summary <- summary(fitStanctModel)
  }

  #######################################################################################################################
  ################################################ CoTiMA Setup #########################################################
  #######################################################################################################################

  if (allInvModel == FALSE) {
    n.TIpred <- (n.studies-1+n.all.moderators+clusCounter); n.TIpred
    driftParamsTmp <- driftParams; driftParamsTmp
    diffParamsTmp  <- diffParams
    meanLag <- mean(allDeltas, na.rm=TRUE); meanLag
    if (customPar == TRUE) {
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
    # CHD 14. Jun 2023
    if ((indVarying == TRUE) & (is.null(indVaryingT0))) indVaryingT0 <- TRUE
    if ((indVarying == 'CINT') & (is.null(indVaryingT0))) indVaryingT0 <- TRUE
    if (is.null(indVaryingT0)) indVaryingT0 <- FALSE

    if ( (randomIntercepts != TRUE) & (randomIntercepts != "MANIFEST") ) {
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
    }
    if (!(is.null(binaries.orig))) {
      # check if really cints rather than manifest means are modelled
      # set TIpredeffects on cints to TRUE
    }

    #stanctModel <- suppressMessages(
    stanctModel <- (
      ctsem::ctModel(n.latent=n.latent, n.manifest=n.var,
                     manifestNames=manifestNames,
                     DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                     DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                     LAMBDA=lambdaParams,
                     CINT=matrix(CINTParams, nrow=n.latent, ncol=1),
                     T0MEANS = matrix(T0meansParams, nrow=n.latent, ncol=1),
                     MANIFESTMEANS = matrix(manifestMeansParams, nrow=n.latent, ncol=1),
                     MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var),
                     T0VAR = T0VARParams,
                     type = 'stanct',
                     n.TIpred = n.TIpred,
                     TIpredNames = paste0("TI", 1:n.TIpred))
    )

    if (indVaryingT0 == TRUE) {
      stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS','indvarying'] <- TRUE
    } else {
      stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS','indvarying'] <- FALSE
    }
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
    tmp2 <- which(stanctModel$pars[tmp1, "param"] %in% invariantDriftParams); tmp2
    tmp3 <- which(is.na(stanctModel$pars[tmp1, "param"])); tmp3
    tmp4 <- sort(unique(c(tmp2, tmp3))); tmp4
    varyingDrifts <- tmp1[!(tmp1 %in% tmp1[tmp4])]; varyingDrifts

    if (length(varyingDrifts) > 0) stanctModel$pars[varyingDrifts, paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- TRUE

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

    if ((randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") ) {
      print(paste0("#################################################################################"))
      print(paste0("#### Note: Correct random intercept model instead of rstan default requested ####"))
      print(paste0("#################################################################################"))

      nullMat <- matrix(0, n.latent, n.latent); nullMat
      DIFFUSIONtmp <- matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent); DIFFUSIONtmp
      DIFFUSIONtmp <- rbind(cbind(DIFFUSIONtmp, nullMat),
                            cbind(nullMat, nullMat)); DIFFUSIONtmp
      DRIFTtmp <- matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent); DRIFTtmp
      DRIFTtmp <- rbind(cbind(DRIFTtmp, diag(n.latent)),
                        cbind(nullMat, nullMat)); DRIFTtmp
      if (randomIntercepts == "MANIFEST") {
        DRIFTtmp <- matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent); DRIFTtmp
        DRIFTtmp <- rbind(cbind(DRIFTtmp, nullMat),
                          cbind(nullMat, nullMat)); DRIFTtmp
      }
      T0VARtmp = gsub("diff", "T0cov", matrix(diffParamsTmp, n.latent, n.latent)); T0VARtmp
      cint_cint_cov <- gsub("eta", "cint", T0VARtmp); cint_cint_cov
      cint_cint_cov <- gsub("T0c", "C", cint_cint_cov); cint_cint_cov
      T0eta_cint_cov <- matrix(NA, n.latent, n.latent); T0eta_cint_cov
      for (ii in 1:n.latent) {
        for (ll in 1:n.latent) {
          T0eta_cint_cov[ii, ll] <- paste0("Cov_T0eta", ll, "_cint", ii)
        }
      }
      T0VARtmp <- rbind(cbind(T0VARtmp, nullMat),
                        cbind(T0eta_cint_cov, cint_cint_cov)); T0VARtmp
      LAMBDAtmp <- cbind(lambdaParams, nullMat); LAMBDAtmp
      if (randomIntercepts == "MANIFEST") {
        LAMBDAtmp <- cbind(lambdaParams, lambdaParams); LAMBDAtmp
      }
      manifestNamesTmp <- manifestNames; manifestNamesTmp
      latentNamesTmp <- c(latentNames, paste0(latentNames, "_cint")); latentNamesTmp
      T0MEANStmp <- matrix(paste0("Mean", latentNamesTmp), n.latent*2, 1); T0MEANStmp
      TIpredNames <- paste0("TI", 1:n.TIpred); TIpredNames
      stanctModel <- (
        ctsem::ctModel(n.latent=n.latent*2,
                       n.manifest=n.var,
                       manifestNames=manifestNamesTmp,
                       latentNames = latentNamesTmp,
                       DIFFUSION=DIFFUSIONtmp,
                       DRIFT=DRIFTtmp,
                       LAMBDA=LAMBDAtmp,
                       CINT=matrix(0, nrow=n.latent*2, ncol=1),
                       T0MEANS = T0MEANStmp,
                       MANIFESTMEANS = matrix(manifestMeansParams, nrow=n.latent, ncol=1),
                       MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var),
                       T0VAR = T0VARtmp,
                       type = 'stanct',
                       n.TIpred = n.TIpred,
                       TIpredNames = paste0("TI", 1:n.TIpred))
      )
      stanctModel$pars$indvarying <- FALSE
      # general setting for params
      # not in this model stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
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
      tmp2 <- which(stanctModel$pars[tmp1, "param"] %in% invariantDriftParams); tmp2
      tmp3 <- which(is.na(stanctModel$pars[tmp1, "param"])); tmp3
      tmp4 <- sort(unique(c(tmp2, tmp3))); tmp4
      varyingDrifts <- tmp1[!(tmp1 %in% tmp1[tmp4])]; varyingDrifts
      if (length(varyingDrifts) > 0) stanctModel$pars[varyingDrifts, paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- TRUE

      # not really necessary
      tmp1 <- which(is.na(stanctModel$pars$param)); tmp1
      tmp2 <- grep("_effect", colnames(stanctModel$pars)); tmp2
      stanctModel$pars[tmp1, tmp2] <- FALSE
      stanctModel$manifesttype <- binaries
      if ( (indVarying == 'CINT') & (!(is.null(binaries.orig))) ) {
        tmp1 <- grep("_effect", colnames(stanctModel$pars)); tmp1
        tmp2 <- which(binaries.orig == 1); tmp2
        stanctModel$pars[(stanctModel$pars$matrix %in% 'CINT'), ][tmp2, tmp1] <- TRUE
      }
    } # end if ((randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )

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
  } # end if (allInvModel == FALSE)

  #stanctModel$pars

    #######################################################################################################################
    ################################################## CoTiMA Fit #########################################################
    #######################################################################################################################

  if (allInvModel == FALSE) {
    #fitStanctModel <- suppressMessages(ctsem::ctStanFit(
    fitStanctModel <- (ctsem::ctStanFit(
      fit=fit,
      datalong = datalong_all,
      ctstanmodel = stanctModel,
      sameInitialTimes=sameInitialTimes,
      savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
      stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
      iter=CoTiMAStanctArgs$iter,
      intoverstates=CoTiMAStanctArgs$intoverstates,
      binomial=CoTiMAStanctArgs$binomial,
      intoverpop=CoTiMAStanctArgs$intoverpop,
      stationary=CoTiMAStanctArgs$stationary,
      plot=CoTiMAStanctArgs$plot,
      optimize=CoTiMAStanctArgs$optimize,
      optimcontrol=CoTiMAStanctArgs$optimcontrol,
      nlcontrol=CoTiMAStanctArgs$nlcontrol,
      priors=CoTiMAStanctArgs$priors, # added Aug 2023
      chains=CoTiMAStanctArgs$chains,
      forcerecompile=CoTiMAStanctArgs$forcerecompile,
      savescores=CoTiMAStanctArgs$savescores,
      gendata=CoTiMAStanctArgs$gendata,
      control=CoTiMAStanctArgs$control,
      #verbose=CoTiMAStanctArgs$verbose,
      verbose=verbose,
      warmup=CoTiMAStanctArgs$warmup,
      cores=coresToUse,
      inits=inits))

    if (is.null(fitStanctModel$standata$priors)) fitStanctModel$standata$priors <- FALSE # CHD added Sep 2023
    fitStanctModel_summary <- summary(fitStanctModel, digits=2*digits, parmatrices=TRUE, residualcov=FALSE)

    if (fit == FALSE) {
      print(paste0("#################################################################################"))
      print(paste0("#############  No model is fitted, only data and code are generated. ############"))
      print(paste0("#################################################################################"))
    }

  } # end if (allInvModel == FALSE)

  #######################################################################################################################
  ####################################### Extract estimates & statistics ################################################
  #######################################################################################################################

  if (fit == TRUE) { # CHD 16. Oct 2023 (end ~line 1900)
    # CHD 22.1.2024
    if ( ( (indVarying == TRUE) | (indVarying == 'CINT') ) & ( (randomIntercepts != TRUE) | (randomIntercepts != "MANIFETS") ) ) {
      e <- ctsem::ctExtract(fitStanctModel)
      # CHD 12.6.2023
      model_popsd <- fitStanctModel_summary$popsd
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
      }
    }

    if ( (indVarying != TRUE) & (indVarying != 'CINT') & (randomIntercepts != TRUE) & (randomIntercepts != "MANIFEST") ) {
      model_popsd <- "no random Intercepts estimated"
      model_popcov_m <- model_popcov_sd <- model_popcov_T <- model_popcov_025 <- model_popcov_50 <- model_popcov_975 <- "no random intercepts estimated"
      model_popcor_m <- model_popcor_sd <- model_popcor_T <- model_popcor_025 <- model_popcor_50 <- model_popcor_975 <- "no random intercepts estimated"
    }

    if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST")  ) {
      #parNames <- ctsem:::getparnames(fitStanctModel); parNames
      # since getparnames is not exported, I took part fo the function and replicated it here # CHD 26.1.2024
      ms <- fitStanctModel$setup$matsetup
      indices <- ms$when %in% c(0, -1) & ms$param > 0 & ms$copyrow < 1
      pars <- data.frame(parnames = ms$parname[indices], parindices = ms$param[indices])
      pars <- pars[!duplicated(pars$parnames), ]
      pars <- pars[order(pars$parindices), ]
      parNames <- pars[pars$parindices > 0, 1]

      targetCols <- grep("ov", parNames); targetCols
      targetRaws <- fitStanctModel$stanfit$rawposterior[, targetCols]
      colnames(targetRaws) <- parNames[grep("ov", parNames)]
      rawT0varTmp <- targetRaws
      # tform T0variances = random intercept covariances - just as a check that this is identical to ctmaInit
      {
        T0COVCoeff <- T0COVCoeffMean <- T0COVCoeffSD <- list()
        n.mod.values.to.plot <- length(colnames(fitStanctModel$data$tipredsdata)[1:(n.studies-1)])+1; n.mod.values.to.plot
        modPos <- 1:(n.studies-1); modPos
        TIpred.values <- matrix(fitStanctModel$data$tipredsdata[, modPos], ncol=length(modPos)); head(TIpred.values)
        effectCodingWeights <- unique(TIpred.values); effectCodingWeights
        #
        tmp1 <- which(!(is.na(fitStanctModel$ctstanmodelbase$pars$param))); tmp1
        tmpPars <- fitStanctModel$ctstanmodelbase$pars[tmp1,]; tmpPars
        T0varPos <- which(tmpPars$matrix == "T0VAR"); T0varPos
        rawT0varTmp <- fitStanctModel$stanfit$rawposterior[ , T0varPos]; head(rawT0varTmp)
        #
        tmpNames <- paste0("RI Covriances for Study No ", unlist(lapply(ctmaInitFit$studyList, function(x) x$originalStudyNo)), "."); tmpNames
        #
        TIpredEffTmp <- fitStanctModel$stanfit$transformedparsfull$TIPREDEFFECT[,T0varPos, modPos]; TIpredEffTmp
        #  n.studies should be > 2 for regular computation because otherwise there are more effectCodingWeights
        for (d in 1:nrow(effectCodingWeights)) {
          if (n.studies == 2) {
            tmp2 <- apply(TIpredEffTmp %*% t(effectCodingWeights[d, ]), 1, sum); tmp2
          } else {
            tmp2 <- apply(TIpredEffTmp %*% effectCodingWeights[d, ], 1, sum); tmp2
          }
          tmp1 <- t(apply(rawT0varTmp, 1 , function(x) x+tmp2))
          T0COVCoeff[[tmpNames[d]]] <- tmp1
        }
        T0COVCoeff
        tmp1a <- fitStanctModel$ctstanmodelbase$pars[, "transform"]; tmp1a
        tmp1b <- fitStanctModel$ctstanmodelbase$pars[, "param"]; tmp1b
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
        # do the same vor the average matrix
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
      }
      #T0varMeanMean

      model_popsd <- "Random intercepts were estimated per primary study. It is recommended to compare them with ctmaInit results to ensure the present results are accurate. Currently, thes should correspond to the rawpopcovbase-slot in the ctsem fit object"
      model_popcov_m <- T0COVCoeffMean
      model_popcov_sd <- T0COVCoeffSD
      #model_popcov_m <- model_popcov_sd <- model_popcov_T <- model_popcov_025 <- model_popcov_50 <- model_popcov_975 <- "random intercepts were estimated but output is not yet available"
      #model_popcor_m <- model_popcor_sd <- model_popcor_T <- model_popcor_025 <- model_popcor_50 <- model_popcor_975 <- "random intercepts were estimated but output is not yet available"
    } # end if (randomIntercepts == TRUE)...

    # account for changes in ctsem 3.4.1
    if ("matrix" %in% colnames(fitStanctModel_summary$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
    tmpMean <- grep("ean", colnames(fitStanctModel_summary$parmatrices)); tmpMean
    tmpSd <- tmpMean+1; tmpSd
    Tvalues <- fitStanctModel_summary$parmatrices[,tmpMean]/fitStanctModel_summary$parmatrices[,tmpSd]; Tvalues
    invariantDrift_Coeff <- cbind(fitStanctModel_summary$parmatrices, Tvalues); invariantDrift_Coeff
    invariantDrift_Coeff[, tmpMean:(dim(invariantDrift_Coeff)[2])] <- round(invariantDrift_Coeff[, tmpMean:(dim(invariantDrift_Coeff)[2])], digits); invariantDrift_Coeff
    #invariantDrift_Coeff

    # (create & ) re-label rownames
    {
      if (ctsem341 == TRUE) {
        tmp1 <- which(invariantDrift_Coeff[, "matrix"] == "DRIFT"); tmp1
        if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") ) {
          tmp1 <- tmp1[which( (invariantDrift_Coeff[tmp1, "row"] <= n.latent) & (invariantDrift_Coeff[tmp1, "col"] <= n.latent) )]; tmp1
        }
        # replace row numbers by newly created names
        rownames(invariantDrift_Coeff) <- paste0(invariantDrift_Coeff[, c("matrix")], "_",
                                                 invariantDrift_Coeff[, c("row")], "_",
                                                 invariantDrift_Coeff[, c("col")])
      } else {
        tmp1 <- which(rownames(invariantDrift_Coeff) == "DRIFT"); tmp1
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
      #invariantDrift_Coeff
      #fitStanctModel_summary

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
      if ( (randomIntercepts == TRUE) |  (randomIntercepts == "MANIFEST") ) {
        tmp <- which(model_Diffusion_Coef != 0)
        model_Diffusion_Coef <- model_Diffusion_Coef[tmp]
        names(model_Diffusion_Coef) <- diffFullNames[tmp]; model_Diffusion_Coef
      } else {
        names(model_Diffusion_Coef) <- diffFullNames; model_Diffusion_Coef
      }

      model_T0var_Coef <- invariantDrift_Coeff[(rownames(invariantDrift_Coeff) == "T0VAR"), 3]; model_T0var_Coef
      if (length(model_T0var_Coef) < 1) {
        model_T0var_Coef <- invariantDrift_Coeff[invariantDrift_Coeff[, "matrix"] == "T0cov",  tmpMean]; model_T0var_Coef
      } else {
        model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
      }
      if ( (randomIntercepts == TRUE) |  (randomIntercepts == "MANIFEST") ) model_T0var_Coef <- c(T0varMean[1:n.latent, 1:n.latent])
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

    if (ctsem341 == TRUE) invariantDrift_Coeff[, "matrix"] <- NULL

    ### Numerically compute Optimal Time lag sensu Dormann & Griffin (2015)
    if (ctsem341 == TRUE) {
      driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=TRUE); driftMatrix # byrow set because order is different compared to mx model
    } else {
      driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=FALSE); driftMatrix # byrow set because order is different compared to mx model
    }

    if (length(invariantDriftNames) == length(driftNames)) {
      OTL <- function(timeRange) {
        OpenMx::expm(tmpDriftMatrix * timeRange)[targetRow, targetCol]}
      # use original time scale
      #if (!(is.null(scaleTime))) {
        tmpDriftMatrix <- driftMatrix * scaleTime
      #} else {
      #  tmpDriftMatrix <- driftMatrix
      #}
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
    #if ( (!(is.null(scaleTime))) & (length(invariantDriftNames) == length(driftNames)) ) {
    if (length(invariantDriftNames) == length(driftNames))  {
      optimalCrossLag_scaledTime <- optimalCrossLag * scaleTime
    }
    else {
      optimalCrossLag_scaledTime <- 'Not available for this model.'
    }

    #######################################################################################################################


    meanDeltas <- mean(ctmaInitFit$statisticsList$allDeltas, na.rm=TRUE); meanDeltas
    largeDelta <- which(ctmaInitFit$statisticsList$allDeltas >= meanDeltas); largeDelta
    tmp1 <- table(ctmaInitFit$statisticsList$allDeltas[largeDelta]); tmp1
    tmp2 <- which(names(tmp1) == (max(as.numeric(names(tmp1))))); tmp2
    suggestedScaleTime <- as.numeric(names(tmp1[tmp2])); suggestedScaleTime
    suggestedScaleTime <- round(suggestedScaleTime, digits); suggestedScaleTime
    message <- c()
    #if (meanDeltas > 3) {
    # CHD AUG 2023
    #if ((meanDeltas * CoTiMAStanctArgs$scaleTime) > 3) { #
    if ((meanDeltas * scaleTime) > 3) { #
      tmp2 <- paste0("Mean time interval was ", meanDeltas, "."); tmp2
      tmp3 <- paste0("scaleTime=1/", suggestedScaleTime); tmp3
      tmp4 <- paste0("It is recommended to fit the model again using the arguments ", tmp3, " and customPar=FALSE. "); tmp4
      message <- paste(tmp2, tmp4, "If the model fit (-2ll) is better (lower), continue using, e.g.,", tmp3, "in all subsequent models.", collapse="\n"); message
    }

    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

    if (ctsem341 == TRUE) { # new 8.7.2022
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

    # new 18.12.2023
    DRIFTCoeff <- list()
    DRIFTCoeffMean <- DRIFTCoeffSD <- list()
    if (WEC == TRUE) {
      n.mod.values.to.plot <- length(colnames(fitStanctModel$data$tipredsdata)[1:(n.studies-1)])+1; n.mod.values.to.plot
      modPos <- 1:(n.studies-1); modPos
      TIpred.values <- matrix(fitStanctModel$data$tipredsdata[, modPos], ncol=length(modPos)); head(TIpred.values)
      effectCodingWeights <- unique(TIpred.values); effectCodingWeights
      #
      tmp1 <- which(!(is.na(fitStanctModel$ctstanmodelbase$pars$param))); tmp1
      tmpPars <- fitStanctModel$ctstanmodelbase$pars[tmp1,]; tmpPars
      driftPos <- which(tmpPars$matrix == "DRIFT"); driftPos
      rawDriftTmp <- fitStanctModel$stanfit$rawposterior[ , driftPos]
      #str(rawDriftTmp)
      #
      tmpNames <- paste0("Drift for Study No ", unlist(lapply(ctmaInitFit$studyList, function(x) x$originalStudyNo)), "."); tmpNames
      #
      TIpredEffTmp <- fitStanctModel$stanfit$transformedparsfull$TIPREDEFFECT[,driftPos, modPos]; TIpredEffTmp
      for (d in 1:nrow(effectCodingWeights)) {
        tmp2 <- apply(TIpredEffTmp %*% effectCodingWeights[d, ], 1, sum); tmp2
        tmp1 <- t(apply(rawDriftTmp , 1 , function(x) x + tmp2))
        DRIFTCoeff[[tmpNames[d]]] <- tmp1
      }

      tmp1a <- fitStanctModel$ctstanmodelbase$pars[, "transform"]; tmp1a
      tmp1b <- fitStanctModel$ctstanmodelbase$pars[, "param"]; tmp1b
      tmp1c <- tmp1b %in% driftNames; tmp1c
      transforms <- tmp1a[tmp1c]; transforms
      #
      # compute tformed drift effects
      for (k in 1:(length(DRIFTCoeff))) {
        counter <- 0
        for (l in 1:(n.latent)) {
          for (m in 1:(n.latent)) {
            counter <- counter + 1
            paramTmp <- DRIFTCoeff[[k]][, counter]; paramTmp
            for (p in 1:length(paramTmp)) {
              param <- paramTmp[p]
              DRIFTCoeff[[k]][p , counter] <- eval(parse(text=transforms[counter]))
            }
            DRIFTCoeff[[k]][, counter] <- DRIFTCoeff[[k]][ ,counter] * scaleTime2
          }
        }
      }
      for (p in 1:length(DRIFTCoeff)) {
        DRIFTCoeffMean[[p]] <- matrix(apply(DRIFTCoeff[[p]], 2, mean), n.latent, n.latent, byrow=T)
        DRIFTCoeffSD[[p]] <- matrix(apply(DRIFTCoeff[[p]], 2, sd), n.latent, n.latent, byrow=T)
      }
    } # end if (WEC == TRUE)

    if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") ) {
      randomIntercepts <- list(popsd=model_popsd,
                               popcov_mean=model_popcov_m,
                               popcov_sd=model_popcov_sd)
      # move T0Means of phantom latents to cints
      tmp1 <- grep("T0MEAN", rownames(invariantDrift_Coeff)); tmp1
      tmp1 <- tmp1[(n.latent+1):(n.latent^2)]; tmp1
      tmp1 <- invariantDrift_Coeff[tmp1, ]; tmp1
      tmp2 <- grep("CINT", rownames(invariantDrift_Coeff)); tmp2
      invariantDrift_Coeff[tmp2[1:n.latent],] <- tmp1
      # delete irrelevant rows
      toDelete <- c()
      for (i in (n.latent+1):(n.latent^2)) {
        toDelete <- c(toDelete, grep(i, rownames(invariantDrift_Coeff)))
      }
      invariantDrift_Coeff <- invariantDrift_Coeff[-toDelete, ]
      toDelete <- c()
      for (i in (n.latent+1):(n.latent^2)) {
        toDelete <- c(toDelete, grep(i, rownames(estimates_original_time_scale)))
      }
      estimates_original_time_scale <- estimates_original_time_scale[-toDelete, ]
    } else {
      if ( (indVarying == 'CINT') | (indVarying == TRUE)  | (indVarying != FALSE)){
        randomIntercepts <- list(popsd=model_popsd,
                                 popcov_mean=model_popcov_m, model_popcov_sd=model_popcov_sd,
                                 model_popcov_T=model_popcov_T, model_popcov_025=model_popcov_025,
                                 model_popcov_50=model_popcov_50, model_popcov_975=model_popcov_975,
                                 popcor_mean=model_popcor_m, model_popcor_sd=model_popcor_sd,
                                 model_popcor_T=model_popcor_T, model_popcor_025=model_popcor_025,
                                 model_popcor_50=model_popcor_50, model_popcor_975=model_popcor_975)
      } else {
        randomIntercepts <- "Not estimated."
      }
    }

    results <- list(
      plot.type="drift",  model.type="stanct",
      n.studies=1,
      n.latent=n.latent,
      n.moderators=length(mod.number),
      studyList=ctmaInitFit$studyList,
      studyFitList=fitStanctModel,
      data=datalong_all,
      statisticsList=ctmaInitFit$statisticsList,
      argumentList=list(ctmaInitFit=ctmaInitFitName,
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
                        #nopriors=nopriors,
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
                        manifestVars=manifestVarsParams,
                        WEC=WEC,
                        randomIntercepts=randomInterceptsSettings,
                        CoTiMAStanctArgs=CoTiMAStanctArgs,
                        T0var=T0var,
                        priors=priors,
                        binaries=binaries,
                        T0var=T0var,
                        cint=cint,
                        indVaryingT0=indVaryingT0,
                        fit=fit
      ),
      modelResults=list(DRIFToriginal_time_scale=model_Drift_Coef_original_time_scale,
                        DIFFUSIONoriginal_time_scale=model_Diffusion_Coef_original_time_scale,
                        T0VAR=model_T0var_Coef,
                        CINT=model_Cint_Coef,
                        MOD=modTI_Coeff_original_time_scale,
                        CLUS=clusTI_Coeff_original_time_scale,
                        DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef,
                        DRIFTCoeffViaWECMean=DRIFTCoeffMean,
                        DRIFTCoeffViaWECSD=DRIFTCoeffSD),
      ctModel = fitStanctModel$ctstanmodelbase,
      parameterNames=ctmaInitFit$parameterNames,
      #CoTiMAStanctArgs=CoTiMAStanctArgs,
      invariantDrift=invariantDrift,
      summary=list(model="Model name not specified", # paste(invariantDrift, "unequal but invariant across samples", collapse=" "),
                   scaledTime=scaleTime2,
                   estimates=invariantDrift_Coeff,
                   # changed 17. Aug. 2022
                   randomIntercepts=randomIntercepts,
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

  } #end if (fit == TRUE)


  if (fit == FALSE) {
    results <- list(summary=c("No model was fitted, only data and code were generated. See $data & $ctModel section."),
                    data = datalong_all,
                    ctModel = fitStanctModel$ctstanmodelbase)
  }

  class(results) <- "CoTiMAFit"

  #######################################################################################################################
  ####################################### Excel Workbook with several sheets ############################################
  #######################################################################################################################


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
