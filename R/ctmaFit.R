#' ctmaFit
#'
#' @description Fits a ctsem model with invariant drift effects across primary studies, possible multiple moderators (but all of them of the
#' the same type, either "cont§ or "cat"), and possible cluster (e.g., countries where primary studies were conducted).
#'
#' @param ctmaInitFit object to which all single ctsem fits of primary studies has been assigned to (i.e., what has been returned by ctmaInit)
#' @param primaryStudyList  could be a list of primary studies compiled with ctmaPrep that defines the subset of studies in ctmaInitFit that should actually be used
#' @param cluster  vector with cluster variables (e.g., countries). Has to be set up carfully. Will be included in ctmaPrep later.
#' @param activeDirectory  defines another active directory than the one used in ctmaInitFit
#' @param activateRPB  set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param digits Number of digits used for rounding (in outputs)
#' @param invariantDrift  drift Labels for drift effects that are set invariant across primary studies (default = all drift effects).
#' @param drift Labels for drift effects. Have to be either of the type V1toV2 or 0 for effects to be excluded, which is usually not recommended)
#' @param moderatedDrift drift Labels for drift effects that are moderated (default = all drift effects) (PRELIMINARY TEST)
#' @param mod.number which in the vector of moderator values to use (e.g., 2 for a single moderator or 1:3 for 3 moderators simultaneously)
#' @param mod.type "cont" or "cat" (mixing them in a single model not yet possible)
#' @param mod.names vector of names for moderators used in output
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#' @param indVarying Allows ct intercepts to vary at the individual level (random effects model, accounts for unobserved heterogeneity)
#' @param scaleTI scale TI predictors - not recommended if TI are dummies representing primary studies as probably in most instances
#' @param scaleClus scale vector of cluster indicators - TRUE (default) yields avg. drift estimates, FALSE yields drift estimates of last cluster
#' @param scaleMod scale moderator variables - TRUE (default) highly recommended
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param nopriors if TRUE, any priors are disabled – sometimes desirable for optimization
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param iter number of interation (defaul = 1000). Sometimes larger values could be reqiured fom Baysian estimation
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#' @param allInvModel estimates a model with all parameters invariant (DRIFT, DIFFUSION, T0VAR)
#' @param customPar logical. Leverages the first pass using priors and ensure that the drift diagonal cannott easily go too negative (could help with ctsem > 3.4)
#' @param equalDrift Not enabled
#' @param inits vector of start values
#'
#' @importFrom  RPushbullet pbPost
#' @importFrom  crayon red
#' @importFrom  parallel detectCores
#' @importFrom  ctsem ctWideToLong ctDeintervalise ctModel ctStanFit
#' @importFrom  OpenMx vech2full expm
#' @importFrom openxlsx addWorksheet writeData createWorkbook openXL saveWorkbook
#'
#' @export ctmaFit
#'
#' @examples
#' \donttest{
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
  ctmaInitFit=NULL,
  primaryStudyList=NULL,
  cluster=NULL,
  activeDirectory=NULL,
  activateRPB=FALSE,
  digits=4,
  drift=NULL,
  invariantDrift=NULL,
  moderatedDrift=NULL,
  equalDrift=NULL,
  mod.number=NULL,
  mod.type="cont",
  mod.names=NULL,
  #n.manifest=0,
  indVarying=FALSE,
  coresToUse=c(1),
  scaleTI=NULL,
  scaleMod=NULL,
  scaleClus=NULL,
  scaleTime=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  finishsamples=NULL,
  iter=NULL,
  chains=NULL,
  verbose=NULL,
  allInvModel=FALSE,
  customPar=TRUE,
  inits=NULL
)


{  # begin function definition (until end of file)

  # adapt display of information during model fit
  if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0} else {verbose <- CoTiMAStanctArgs$verbose}

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
    stop("Good luck for the next try!")
  }

  { # set fitting params
    if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
    if (!(is.null(scaleClus))) CoTiMAStanctArgs$scaleClus <- scaleClus
    if (!(is.null(scaleMod))) CoTiMAStanctArgs$scaleMod <- scaleMod
    if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
    if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
    if (!(is.null(nopriors))) CoTiMAStanctArgs$nopriors <- nopriors
    if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
    if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
    if (!(is.null(iter))) CoTiMAStanctArgs$iter <- iter
    if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
  }


  #######################################################################################################################
  ####### Copy/Change INIT File based on information delivered by different PREP files (e.g., moderator studies ) #######
  #######################################################################################################################

  # if primary study list is provided in addition to initfit-object, take primary study information from  primary study
  if (!(is.null(primaryStudyList))) {
    ctmaTempFit <- ctmaInitFit
    targetStudyNumbers <- unlist(primaryStudyList$studyNumbers); targetStudyNumbers; length(targetStudyNumbers)
    targetStudyNumbers <- targetStudyNumbers[-which(is.na(targetStudyNumbers))]; targetStudyNumbers
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
      cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
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
    if (!(is.null(ctmaInitFit$n.manifest))) n.manifest <- ctmaInitFit$n.manifest
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    usedTimeRange <- seq(0, 1.5*maxDelta, 1)
    lambda <- ctmaInitFit$statisticsList$lambda; lambda
  }


  # check match between cluster vector and n.studies
  if (!(is.null(cluster))) {
    if (length(cluster) != n.studies) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("The vector of cluster numbers does not match the number of primary studies.\n"))
      stop("Good luck for the next try!")
    }
  }

  # check moderator information
  {
    n.moderators <- length(mod.number); n.moderators
    if (n.moderators > 0) {
      currentModerators <- matrix(as.numeric(unlist(lapply(ctmaInitFit$studyList, function(extract) extract$moderators[mod.number]))), ncol=n.moderators); currentModerators
      if (!is.null(primaryStudyList)) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudyList$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE); currentModerators
      if (is.na((currentModerators[length(currentModerators)])[[1]][1])) currentModerators <- currentModerators[-dim(currentModerators)[1],]; currentModerators
      if (is.null(dim(currentModerators)[1])) currentModerators <- matrix(currentModerators, ncol=1); currentModerators

      if (any(is.na(currentModerators)) == TRUE) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        cat(crayon::red$bold("At least one of the primary studies does not have a valid value for the requested moderator.  \n"))
        stop("Good luck for the next try!")
      }
      if (var(currentModerators) == 0) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        cat(crayon::red$bold("Moderator is constant across cases.", sep="\n"))
        stop("Good luck for the next try!")
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
        tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit, moderatorValues= currentModerators)
      } else {
        tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit)
      }
      datawide_all <- tmp$alldata
      groups <- tmp$groups # vector of study IDs
    }

    # delete cases with missing moderator values
    if (n.moderators > 0) {
      casesToDelete <- tmp$casesToDelete; casesToDelete
      if (!(is.null(casesToDelete))) {
        datawide_all <- datawide_all[-casesToDelete, ]
        groups <- groups[-casesToDelete]
        n.studies <- length(unique(groups))
      }
    }

    # make data matrix with moderators
    if (n.moderators > 0) {
      moderatorGroups <- tmp$moderatorGroups
      if (!(is.matrix(moderatorGroups))) moderatorGroups <- matrix(moderatorGroups, ncol=1)
      colnames(moderatorGroups) <- paste0("mod", 1:(dim(currentModerators)[2])); colnames(moderatorGroups)
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
    #modTIstartNum <- n.studies+1; modTIstartNum
    modTIstartNum <- n.studies; modTIstartNum
    if (n.moderators > 0) {
      # make TI predictors as replacements for each moderator
      if (mod.type=="cont") {
        tmp1 <- paste0("mod", 1:n.moderators); tmp1; dim(tmp1)
        if (length(tmp1) == 1) tmp <- matrix(dataTmp[ , tmp1], ncol=length(tmp1)) else tmp <- dataTmp[ , tmp1]
        if (CoTiMAStanctArgs$scaleMod == TRUE) tmp[ , 1:ncol(tmp)] <- scale(tmp[ , 1:ncol(tmp)])
        currentStartNumber <- modTIstartNum; currentStartNumber
        currentEndNumber <- currentStartNumber + n.moderators-1; currentEndNumber
        colnames(tmp) <- paste0("TI", currentStartNumber:currentEndNumber); tmp
        dataTmp <- cbind(dataTmp, tmp); dim(dataTmp)
        dataTmp <- dataTmp[ ,-grep("mod", colnames(dataTmp))]
      }
      if (mod.type=="cat") {
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
            counter <- counter + 1
            tmp2 <- which(tmp == unique.mod[i]); tmp2
            tmpTI[tmp2, counter] <- 1
          }
        }

        if (CoTiMAStanctArgs$scaleMod == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)])
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
      # create dummies
      tmpTI <- matrix(0, dim(dataTmp)[1], clusCounter)
      for (i in 1:clusCounter) {
        targetGroups <- which(cluster == targetCluster[i]); targetGroups
        tmp2 <- which(groups %in% targetGroups); length(tmp2)
        tmpTI[tmp2, i] <- 1
      }
      if (CoTiMAStanctArgs$scaleClus == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)])
      cluster.weights <- cluster.sizes <- matrix(NA, nrow=ncol(tmpTI), ncol=2)
      for (l in 1:ncol(tmpTI)) {
        cluster.weights[l,] <- round(as.numeric(names(table(tmpTI[ , l]))), digits)
        cluster.sizes[l,] <- table(tmpTI[ , l])
      }
      rownames(cluster.weights) <- paste0(1:ncol(tmpTI), "_on__")
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
    #head(datalong_all)
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
  if (is.null(moderatedDrift) & (!(is.null(mod.number)))) moderatedDrift <- "all" # will be changed by ctmaLabels

  namesAndParams <- ctmaLabels(
    n.latent=n.latent,
    n.manifest=n.manifest,
    lambda=lambda,
    drift=drift,
    invariantDrift=invariantDrift,
    moderatedDrift=moderatedDrift,
    equalDrift=equalDrift
  )
  driftNames <- namesAndParams$driftNames; driftNames
  driftFullNames <- namesAndParams$driftFullNames; driftFullNames
  driftParams <- namesAndParams$driftParams; driftParams
  diffNames <- namesAndParams$diffNames; diffNames
  diffParams <- namesAndParams$diffParams; diffParams
  diffFullNames <- namesAndParams$diffFullNames; diffFullNames
  invariantDriftNames <- namesAndParams$invariantDriftNames; invariantDriftNames
  invariantDriftParams <- namesAndParams$invariantDriftParams; invariantDriftParams
  moderatedDriftNames <- namesAndParams$moderatedDriftNames; moderatedDriftNames
  equalDriftNames <- namesAndParams$equalDriftNames; equalDriftNames
  equalDriftParams <- namesAndParams$equalDriftParams; equalDriftParams
  lambdaParams <- namesAndParams$lambdaParams; lambdaParams
  T0VARParams <- namesAndParams$T0VARParams; T0VARParams
  manifestmeansParams <- namesAndParams$manifestMeansParams; manifestmeansParams
  manifestVarParams <- namesAndParams$manifestVarParams; manifestVarParams

  if (is.null(invariantDriftNames)) invariantDriftNames <- driftNames

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
                                    finishsamples=finishsamples,
                                    iter=iter,
                                    chains=chains,
                                    verbose=verbose,
                                    indVarying = indVarying,
                                    customPar = customPar)
    fitStanctModel <- allInvModelFit$studyFitList[[1]]
    invariantDriftStanctFit <- summary(fitStanctModel)
  } else {

    n.TIpred <- (n.studies-1+n.all.moderators+clusCounter); n.TIpred
    # scale Drift to cover changes in ctsem 3.4.1 (this would be for ctmaFit/ctmaModFit, but for Init individual study modification is done later)
    driftParamsTmp <- driftParams; driftParamsTmp
    diffParamsTmp  <- diffParams
    meanLag <- mean(allDeltas, na.rm=TRUE); meanLag
    if ((meanLag > 6) & (customPar)) {
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

    # Make model with max time points
    {
      if ((allInvModel == FALSE) & (indVarying == TRUE)) {
        print(paste0("#################################################################################"))
        print(paste0("######## Just a note: Individually varying intercepts model requested.  #########"))
        print(paste0("#################################################################################"))

        #manifestNames <- paste0("y", 1:n.manifest); manifestNames
        #latentNames <- paste0("V", 1:2); latentNames
        #if (n.manifest == n.latent) manifestNames <- latentNames
        manifestmeansParams <- manifestNames; manifestmeansParams
        #if (n.manifest == 0) n.manifestTmp <- n.latent else n.manifestTmp <- n.manifest
        manifestVarParams <- c()
        for (u in 1:n.var) {
          for (v in 1:n.var) {
            if ( v < u) manifestVarParams <- c(manifestVarParams, "0") else manifestVarParams <- c(manifestVarParams, paste0("var_", v, "_", u))
          }
        }
        manifestVarParams

        stanctModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=maxTpoints, manifestNames=manifestNames,
                                      DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                                      DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                                      LAMBDA=lambdaParams,
                                      CINT=matrix(0, nrow=n.latent, ncol=1),
                                      T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                      MANIFESTMEANS = matrix(manifestmeansParams, nrow = n.var, ncol = 1),
                                      MANIFESTVAR=matrix(manifestVarParams, nrow=n.var, ncol=n.var),
                                      type = 'stanct',
                                      n.TIpred = n.TIpred,
                                      TIpredNames = paste0("TI", 1:n.TIpred),
                                      #TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+clusCounter+n.all.moderators))
                                      TIPREDEFFECT = matrix(0, n.latent, n.TIpred))

      } else {
        stanctModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=maxTpoints, manifestNames=manifestNames,
                                      DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                                      DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                                      LAMBDA=lambdaParams,
                                      CINT=matrix(0, nrow=n.latent, ncol=1),
                                      T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                      MANIFESTMEANS = matrix(manifestmeansParams, nrow = n.var, ncol = 1),
                                      MANIFESTVAR=matrix(manifestVarParams, nrow=n.var, ncol=n.var),
                                      type = 'stanct',
                                      n.TIpred = n.TIpred,
                                      #TIpredNames = paste0("TI", 1:(n.studies-1+clusCounter+n.all.moderators)),
                                      TIpredNames = paste0("TI", 1:n.TIpred),
                                      #TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+clusCounter+n.all.moderators
                                      TIPREDEFFECT = matrix(0, n.latent, n.TIpred))
        stanctModel$pars[, "indvarying"] <- FALSE
      }
    }
    #stanctModel$pars

    # general setting for params
    stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'LAMBDA',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'CINT',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTMEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTVAR',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    #n.studies+n.all.moderators+clusCounter-1
    #stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames[(n.studies+n.all.moderators):(n.studies+n.all.moderators+clusCounter-1)],'_effect')]
    if (!(is.null(cluster))) {
      stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames[(n.studies++n.all.moderators):(n.studies+n.all.moderators+clusCounter-1)],'_effect')] <- TRUE
    }

    if (n.moderators > 0) {
      tmp1 <- which(stanctModel$pars$matrix == "DRIFT"); tmp1
      #tmp2 <- which(!(stanctModel$pars[tmp1, "param"] %in% moderatedDrift)); tmp2
      tmp2 <- which((stanctModel$pars[tmp1, "param"] %in% moderatedDriftNames)); tmp2
      #targetCols <- (n.studies-1+clusCounter+1):(n.studies-1+clusCounter+n.moderators); targetCols
      #targetCols <- (n.studies-1+clusCounter+1):(n.studies-1+clusCounter+n.all.moderators); targetCols
      targetCols <- (n.studies):(n.studies-1+n.all.moderators); targetCols
      #stanctModel$pars[ , paste0(stanctModel$TIpredNames[targetCols],'_effect')]
      stanctModel$pars[ , paste0(stanctModel$TIpredNames[targetCols],'_effect')] <- FALSE
      stanctModel$pars[tmp1[tmp2] , paste0(stanctModel$TIpredNames[targetCols],'_effect')] <- TRUE
    }
    stanctModel$pars


    # the target effects
    tmp1 <- which(stanctModel$pars$matrix == "DRIFT"); tmp1
    tmp2 <- which(stanctModel$pars[tmp1, "param"] %in% invariantDriftNames); tmp2
    stanctModel$pars[tmp1[tmp2], paste0(stanctModel$TIpredNames[1:(n.studies-1)],'_effect')] <- FALSE
    #stanctModel$pars

    if (!(optimize)) {
      customPar <- FALSE
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      cat(crayon::red("Bayesian sampling was selected, which does require appropriate scaling of time. See the end of the summary output","\n"))
    }
    #stanctModel$pars

    fitStanctModel <- ctsem::ctStanFit(
      datalong = datalong_all,
      ctstanmodel = stanctModel,
      savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
      stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
      iter=CoTiMAStanctArgs$iter,
      intoverstates=CoTiMAStanctArgs$intoverstates,
      binomial=CoTiMAStanctArgs$binomial,
      fit=CoTiMAStanctArgs$fit,
      intoverpop=CoTiMAStanctArgs$intoverpop,
      stationary=CoTiMAStanctArgs$stationary,
      plot=CoTiMAStanctArgs$plot,
      derrind=CoTiMAStanctArgs$derrind,
      optimize=CoTiMAStanctArgs$optimize,
      optimcontrol=CoTiMAStanctArgs$optimcontrol,
      nlcontrol=CoTiMAStanctArgs$nlcontrol,
      nopriors=CoTiMAStanctArgs$nopriors,
      chains=CoTiMAStanctArgs$chains,
      forcerecompile=CoTiMAStanctArgs$forcerecompile,
      savescores=CoTiMAStanctArgs$savescores,
      gendata=CoTiMAStanctArgs$gendata,
      #control=CoTiMAStanctArgs$control,
      verbose=CoTiMAStanctArgs$verbose,
      warmup=CoTiMAStanctArgs$warmup,
      cores=coresToUse,
      inits=inits)

    ### resample in parcels to avoid memory crash and speed up
    if (!(is.null(CoTiMAStanctArgs$resample))) {
      fitStanctModel <- ctmaStanResample(ctmaFittedModel=fitStanctModel) #, CoTiMAStanctArgs=CoTiMAStanctArgs)
    }
        invariantDriftStanctFit <- summary(fitStanctModel, digits=2*digits, parmatrices=TRUE, residualcov=FALSE)
  } # end if else (allInvModel)

  # Extract estimates & statistics
  # account for changes in ctsem 3.4.1
  if ("matrix" %in% colnames(invariantDriftStanctFit$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
  tmpMean <- grep("ean", colnames(invariantDriftStanctFit$parmatrices)); tmpMean
  tmpSd <- tmpMean+1; tmpSd
  Tvalues <- invariantDriftStanctFit$parmatrices[,tmpMean]/invariantDriftStanctFit$parmatrices[,tmpSd]; Tvalues
  invariantDrift_Coeff <- cbind(invariantDriftStanctFit$parmatrices, Tvalues); invariantDrift_Coeff
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
      #driftFullNames <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); driftFullNames
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
  #invariantDrift_Coeff

  # extract fit
  invariantDrift_Minus2LogLikelihood  <- -2*invariantDriftStanctFit$loglik; invariantDrift_Minus2LogLikelihood
  invariantDrift_estimatedParameters  <- invariantDriftStanctFit$npars; invariantDrift_estimatedParameters
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
  #modTIendNum <- 0
  modTI_Coeff <- NULL
  if (n.moderators > 0) {
    #tmp1 <- tmp2 <- c()

    #if (mod.type == "cont") tmp2 <- length(mod.number)-1
    #if (mod.type == "cat") tmp2 <- length(unlist(unique.mod))-n.moderators
    #modTIendNum <- modTIstartNum+tmp2-1; modTIendNum

    #for (i in modTIstartNum:modTIendNum) tmp1 <- c(tmp1, (grep(i, rownames(invariantDriftStanctFit$tipreds))))
    #for (i in modTIstartNum:modTIendNum) tmp1 <- c(tmp1, (grep(paste0("TI", i), rownames(invariantDriftStanctFit$tipreds))))
    tmp1 <- c()
    invariantDriftStanctFit$tipreds
    for (i in modTIs) tmp1 <- c(tmp1, grep(i, rownames(invariantDriftStanctFit$tipreds)))
    Tvalues <- invariantDriftStanctFit$tipreds[tmp1, ][,6]; Tvalues
    modTI_Coeff <- round(cbind(invariantDriftStanctFit$tipreds[tmp1, ], Tvalues), digits); modTI_Coeff

    # re-label
    if (!(is.null(mod.names))) {
      if (mod.type == "cont") {
        #for (i in 1:n.moderators) {
        counter <- 0
        for (i in modTIs) {
          #i <- 1
          counter <- counter + 1
          #targetNamePart <- paste0("tip_TI", n.studies+i-1); targetNamePart
          targetNamePart <- paste0("tip_", modTIs[counter]); targetNamePart
          rownames(modTI_Coeff) <- sub(targetNamePart, paste0(mod.names[counter], "_on_"), rownames(modTI_Coeff))
        }
      }
      if (mod.type == "cat") {
        counter <- 0
        modNameCounter <- 1
        #for (j in 1:n.moderators) {
        for (j in modTIs) {
          if (n.moderators == 1) unique.mod.tmp <- unique.mod else unique.mod.tmp <- unique.mod[[counter+1]]
          for (i in 1:(length(unique.mod.tmp)-1)) {
            counter <- counter + 1
            current.mod.names <- mod.names[modNameCounter]; current.mod.names
            #targetNamePart <- paste0("tip_TI", n.studies+i-1); targetNamePart
            targetNamePart <- paste0("tip_", modTIs[i]); targetNamePart
            tmp1 <- grep(targetNamePart, rownames(modTI_Coeff)); tmp1
            rownames(modTI_Coeff) <- sub(targetNamePart, paste0(counter+1, ". smallest value (category) of ", mod.names[modNameCounter], "_on"), rownames(modTI_Coeff))
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
    for (i in clusTIs) tmp1 <- c(tmp1, (grep(i, rownames(invariantDriftStanctFit$tipreds))))
    Tvalues <- invariantDriftStanctFit$tipreds[tmp1, ][,6]; Tvalues
    clusTI_Coeff <- round(cbind(invariantDriftStanctFit$tipreds[tmp1, ], Tvalues), digits); clusTI_Coeff
    # re-label
    targetCluster <- which(table(cluster) > 1); targetCluster
    for (i in 1:clusCounter) {
      targetNamePart <- paste0("tip_", clusTIs[i]); targetNamePart
      rownames(clusTI_Coeff) <- sub(targetNamePart, paste0("Cluster_", targetCluster[i], "_on_"), rownames(clusTI_Coeff))
      for (j in 1:length(driftNames)) {
        if (driftNames[j] != 0) {
          tmp0 <- driftNames[j]; tmp0
          tmp0 <- gsub(" \\(invariant\\)", "", tmp0); tmp0
          tmp1 <- grep(tmp0, rownames((clusTI_Coeff))); tmp1
          tmp2 <- grep(tmp0, names((model_Drift_Coef))); tmp2
          cluster.specific.effect[i,j] <- round(model_Drift_Coef[tmp2] + clusTI_Coeff[tmp1, 1] * cluster.weights[i, 2], digits)
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
      OpenMx::expm(driftMatrix * timeRange)[targetRow, targetCol]}
    # loop through all cross effects
    optimalCrossLag <- matrix(NA, n.latent, n.latent)
    maxCrossEffect <- matrix(NA, n.latent, n.latent)
    for (j in 1:n.latent) {
      for (h in 1:n.latent) {
        if (j != h) {
          targetRow <- j
          targetCol <- h
          if (driftMatrix[j, h] != 0) { # an effect that is zero has no optimal lag
            targetParameters <- sapply(usedTimeRange, OTL)
            maxCrossEffect[j,h] <- max(abs(targetParameters))
            optimalCrossLag[j,h] <- which(abs(targetParameters)==maxCrossEffect[j,h])*1+0
          } else {
            optimalCrossLag[j,h] <- NA
          }
        }
      }
    }
    maxCrossEffect <- round(maxCrossEffect, digits)
  } else {
    optimalCrossLag <- "Drift Matrix is only partially invariant. (Generalizable) optimal intervall cannot be computed."
    maxCrossEffect <- "Drift Matrix is only partially invariant. (Generalizable) Largest effects cannot be computed."
  }
  #} ## END  fit stanct model


  #######################################################################################################################

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  st <- paste0("Computation started at: ", start.time); st
  et <- paste0("Computation ended at: ", end.time); et
  tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt


  meanDeltas <- mean(ctmaInitFit$statisticsList$allDeltas, na.rm=TRUE); meanDeltas
  largeDelta <- which(ctmaInitFit$statisticsList$allDeltas >= meanDeltas); largeDelta
  tmp1 <- table(ctmaInitFit$statisticsList$allDeltas[largeDelta]); tmp1
  tmp2 <- which(tmp1 == (max(tmp1))); tmp2
  suggestedScaleTime <- as.numeric(names(tmp1[tmp2])); suggestedScaleTime
  message <- c()
  if (meanDeltas > 3) {
    tmp2 <- paste0("Mean time interval was ", meanDeltas, "."); tmp2
    tmp3 <- paste0("scaleTime=1/", suggestedScaleTime); tmp3
    tmp4 <- paste0("It is recommended to fit the model again using the arguments ", tmp3, " and customPar=FALSE. "); tmp4
    message <- paste(tmp2, tmp4, "If the model fit (-2ll) is better (lower), continue using, e.g.,", tmp3, "in all subsequent models.", collapse="\n"); message
  }


  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}


  tmp1 <- grep("CINT", rownames(invariantDriftStanctFit$parmatrices)); tmp1
  tmp2 <- grep("asym", rownames(invariantDriftStanctFit$parmatrices)); tmp2
  tmp3 <- grep("dt", rownames(invariantDriftStanctFit$parmatrices)); tmp3
  tmp4 <- tmp1[(tmp1 %in% c(tmp2, tmp3)) == FALSE]; tmp4
  model_Cint_Coef <- invariantDriftStanctFit$parmatrices[tmp4, 3]; model_Cint_Coef

  if (!(is.null(cluster))) {
    clus.effects=list(effects=clusTI_Coeff, weights=cluster.weights, sizes=cluster.sizes,
                      cluster.specific.effect=cluster.specific.effect, note=cluster.note)
  } else {
    clus.effects <- NULL
  }

  results <- list(activeDirectory=activeDirectory,
                  plot.type="drift",  model.type="stanct",
                  coresToUse=coresToUse,
                  n.studies=1,
                  n.latent=n.latent,
                  n.moderators=length(mod.number),
                  mod.names=mod.names,
                  mod.type=mod.type,
                  studyList=ctmaInitFit$studyList,
                  studyFitList=fitStanctModel,
                  data=datalong_all,
                  statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef,
                                    CINT=model_Cint_Coef, MOD=modTI_Coeff, CLUS=clusTI_Coeff),
                  parameterNames=ctmaInitFit$parameterNames,
                  CoTiMAStanctArgs=CoTiMAStanctArgs,
                  invariantDrift=invariantDrift,
                  summary=list(model="Model name not specified", # paste(invariantDrift, "unequal but invariant across samples", collapse=" "),
                               estimates=invariantDrift_Coeff,
                               randomEffects=invariantDriftStanctFit$popsd,
                               minus2ll= invariantDrift_Minus2LogLikelihood,
                               n.parameters = invariantDrift_estimatedParameters,
                               #df= invariantDrift_df,
                               opt.lag = optimalCrossLag,
                               max.effects = maxCrossEffect,
                               clus.effects=clus.effects,
                               mod.effects=modTI_Coeff,
                               message=message)
                  # excel workboo is added later
                  )

  class(results) <- "CoTiMAFit"

  ### prepare Excel Workbook with several sheets ################################################################
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
                        matrix(unlist(results$modelResults$DRIFT),
                               nrow=1, ncol=n.latent^2, byrow=TRUE))
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
    results$modelResults
    openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                        matrix(names(results$modelResults$DIFFUSION),
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

  invisible(results)

} ### END function definition
