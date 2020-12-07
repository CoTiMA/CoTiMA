############################################ CoTiMA Moderator STANCT ##################################################
#######################################################################################################################
#' ctmaModFit
#'
#' @description Performs moderator analysis of one or more drift coefficients moderated by one or more time independent moderators
#'
#' @param ctmaInitFit object to which all single ctsem fits of primary studies has been assigned to (i.e., what has been returned by ctmaInit)
#' @param primaryStudyList  could be a list of primary studies compiled with ctmaPrep that defines the subset of studies in ctmaInitFit that should actually be used
#' @param cluster  vector with cluster variables (e.g., countries). Has to be set up carfully. Will be included in ctmaPrep later.
#' @param activeDirectory  defines another active directory than the one used in ctmaInitFit
#' @param activateRPB  set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param mod.number which in the vector of moderator values to use (e.g., 2 for a single moderator or 1:3 for 3 moderators simultaneously)
#' @param mod.type "cont" or "cat" (mixing them in a single model not yet possible)
#' @param mod.names vector of names for moderators used in output
#' @param digits Number of digits used for rounding (in outputs)
#' @param moderatedDrift drift Labels for drift effects that are moderated (default = all drift effects)
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#' @param indVarying Allows ct intercepts to vary at the individual level (random effects model, accounts for unobserved heteregeneity)
#' @param scaleTI scale TI predictors - not recommended if TI are dummies representing primary studies as probably in most instances
#' @param scaleMod scale moderators - useful in particular if continuous moderators are used
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param optimize if set to FALSE, Stan’s Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param nopriors if TRUE, any priors are disabled – sometimes desirable for optimization
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red
#' @importFrom parallel detectCores
#' @importFrom ctsem ctWideToLong ctDeintervalise ctModel ctStanFit
#' @importFrom OpenMx vech2full
#'
#' @export ctmaModFit
#'
#' @examples
#' # Fit a moderated CoTiMA to all primary studies previously fitted
#' # one by one with the fits assigned to CoTiMAInitFit_Ex1
#' \dontrun{
#' CoTiMAModlFit_Ex1 <- ctmaModFit(ctmaInitFit=CoTiMAInitFit_Ex1,
#' mod.number=2, mod.type="cat", mod.names=c("H not stated"))
#'
#' saveRDS(CoTiMAModlFit_Ex1, file=paste0(activeDirectory, "CoTiMAModlFit_Ex1.rds"))
#' summary(CoTiMAModlFit_Ex1)
#' }
#'
ctmaModFit <- function(
  ctmaInitFit=NULL,
  primaryStudyList=NULL,
  cluster=NULL,
  activeDirectory=NULL,
  mod.number=1,
  mod.type="cont",
  mod.names=NULL,
  activateRPB=FALSE,
  digits=4,
  moderatedDrift=c(),
  coresToUse=c(1),
  indVarying=FALSE,
  scaleTI=NULL,
  scaleMod=NULL,
  scaleTime=NULL,
  optimize=TRUE,
  nopriors=TRUE,
  finishsamples=NULL,
  chains=NULL,
  verbose=NULL
)
{  # begin function definition (until end of file)

  { # fitting params
    if (!(is.null(scaleTI))) CoTiMAStanctArgs$scaleTI <- scaleTI
    if (!(is.null(scaleMod))) CoTiMAStanctArgs$scaleMod <- scaleMod
    if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime
    if (!(is.null(optimize))) CoTiMAStanctArgs$optimize <- optimize
    if (!(is.null(nopriors))) CoTiMAStanctArgs$nopriors <- nopriors
    if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
    if (!(is.null(chains))) CoTiMAStanctArgs$chains <- chains
    if (!(is.null(verbose))) CoTiMAStanctArgs$verbose <- verbose
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
  ######### Copy/Change INIT File based on information delivered by the PREP file of moderator analysis #################
  #######################################################################################################################

  if (!(is.null(primaryStudyList))) {
    ctmaModFit <- ctmaInitFit
    targetStudyNumbers <- unlist(primaryStudyList$studyNumbers); targetStudyNumbers; length(targetStudyNumbers)
    for (i in (length(ctmaModFit$studyFitList)):1) {
      if (!(ctmaModFit$studyList[[i]]$originalStudyNo %in% targetStudyNumbers)) {
        ctmaModFit$studyList[[i]] <- NULL
        ctmaModFit$studyFitList[[i]] <- NULL
        ctmaModFit$emprawList[[i]] <- NULL
        ctmaModFit$statisticsList$originalStudyNumbers[i] <- NA
        ctmaModFit$statisticsList$allSampleSizes[i+1] <- NA
        ctmaModFit$statisticsList$allTpoints[i] <- NA
        (ctmaModFit$modelResults[[1]])
        ctmaModFit$modelResults[[1]][[i]] <- NULL
        ctmaModFit$modelResults[[2]][[i]] <- NULL
        ctmaModFit$modelResults[[3]][[i]] <- NULL
      }
    }
    ctmaModFit$n.studies <- length(targetStudyNumbers); ctmaModFit$n.studies
    ctmaModFit$statisticsList$allDeltas <- unlist(lapply(ctmaModFit$studyList, function(extract) extract$delta_t))
    ctmaModFit$statisticsList$minDelta <- min(ctmaModFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaModFit$statisticsList$maxDelta <- max(ctmaModFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaModFit$statisticsList$meanDelta <- mean(ctmaModFit$statisticsList$allDeltas, na.rm=TRUE)
    ctmaModFit$statisticsList$overallSampleSize <- sum(ctmaModFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaModFit$statisticsList$meanSampleSize <- mean(ctmaModFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaModFit$statisticsList$maxSampleSize <- max(ctmaModFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaModFit$statisticsList$minSampleSize <- min(ctmaModFit$statisticsList$allSampleSizes, na.rm = TRUE)
    ctmaModFit$statisticsList$overallTpoints <- sum(ctmaModFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaModFit$statisticsList$meanTpoints <- mean(ctmaModFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaModFit$statisticsList$maxTpoints <- max(ctmaModFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaModFit$statisticsList$minTpoints <- min(ctmaModFit$statisticsList$allTpoints, na.rm = TRUE)
    ctmaModFit$summary$model <- "Moderator Model (for details see model summary)"
    tmpStudyNumber <- as.numeric(gsub("Study No ", "", rownames(ctmaModFit$summary$estimates))); tmpStudyNumber
    targetRows <- which(tmpStudyNumber %in% targetStudyNumbers); targetRows; length(targetRows)
    ctmaModFit$summary$estimates <- ctmaModFit$summary$estimates[targetRows, ]
    ctmaModFit$summary$confidenceIntervals <- ctmaModFit$summary$confidenceIntervals[targetRows, ]
    ctmaModFit$summary$n.parameters <- ctmaModFit$studyFitList[[1]]$resultsSummary$npars * length(targetRows)
    ctmaModFit$statisticsList$originalStudyNumbers <-
      ctmaModFit$statisticsList$originalStudyNumbers[which(!(is.na(ctmaModFit$statisticsList$originalStudyNumbers)))]
    ctmaModFit$statisticsList$allSampleSizes <-
      ctmaModFit$statisticsList$allSampleSizes[which(!(is.na(ctmaModFit$statisticsList$allSampleSizes)))]
    ctmaModFit$statisticsList$allTpoints <-
      ctmaModFit$statisticsList$allTpoints[which(!(is.na(ctmaModFit$statisticsList$allTpoints)))]

    ctmaInitFit <- ctmaModFit
  }




  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################
  {
    start.time <- Sys.time(); start.time
    n.latent <- ctmaInitFit$n.latent; n.latent
    maxTpoints <- ctmaInitFit$statisticsList$maxTpoints; maxTpoints
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints; length(allTpoints)
    allTpoints <- allTpoints[which(!(is.na(allTpoints)))]; allTpoints
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies

    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    allDeltas <- allDeltas[which(!(is.na(allDeltas)))]; allDeltas
    maxDelta <- max(allDeltas); maxDelta
    usedTimeRange <- seq(0, 1.5*maxDelta, 1); usedTimeRange
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    driftNames <- c(t(matrix(ctmaInitFit$parameterNames$DRIFT, n.latent))); driftNames

    if (length(moderatedDrift) < 1) moderatedDrift <- driftNames

    usedTimeRange <- seq(0, 1.5*maxDelta, 1)

    n.moderators <- length(mod.number); n.moderators
    currentModerators <- matrix(as.numeric(unlist(lapply(ctmaInitFit$studyList, function(extract) extract$moderators[mod.number]))), ncol=n.moderators); currentModerators
    #if (!is.null(primaryStudies)) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudies$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE); currentModerators
    if (!is.null(primaryStudyList)) currentModerators <- matrix(as.numeric(unlist(lapply(primaryStudyList$moderators, function(extract) extract[mod.number]))), ncol=n.moderators, byrow=TRUE); currentModerators
    ##currentModerators <- unlist(lapply(primaryStudyList, function(extract) extract$moderators[mod.number])); currentModerators
    if (is.na((currentModerators[length(currentModerators)])[[1]][1])) currentModerators <- currentModerators[-dim(currentModerators)[1],]; currentModerators
    if (is.null(dim(currentModerators)[1])) currentModerators <- matrix(currentModerators, ncol=1)

    if (var(currentModerators) == 0) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("Moderator is constant across cases.", sep="\n"))
      stop("Good luck for the next try!")

    }

    }

    #######################################################################################################################
    ################################################# data preparation ####################################################
    #######################################################################################################################
    {
      # combine pseudo raw data for model
      tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit, moderatorValues = currentModerators)
      datawide_all <- tmp$alldata
      groups <- tmp$groups
      names(groups) <- c("Study_No_"); groups
      groupsNamed <- (paste0("Study_No_", groups)); groupsNamed
      moderatorGroups <- tmp$moderatorGroups
      colnames(moderatorGroups) <- paste0("mod", 1:(dim(currentModerators)[2]))

      # augment pseudo raw data for stanct model
      modTIstartNum <- n.studies; modTIstartNum
      dataTmp <- cbind(datawide_all, groups, moderatorGroups)
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

      # make TI predictors as replacements for each moderator
      if (mod.type=="cont") {
        tmp1 <- paste0("mod", 1:n.moderators); tmp1; dim(tmp1)
        if (length(tmp1) == 1) tmp <- matrix(dataTmp[ , tmp1], ncol=length(tmp1)) else tmp <- dataTmp[ , tmp1]
        if (CoTiMAStanctArgs$scaleMod == TRUE) tmp[ , 1:ncol(tmp)] <- scale(tmp[ , 1:ncol(tmp)])
        tmp[1000:1360]
        scale(tmp[ , 1:ncol(tmp)])
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

      # add clusters as dummy moderators
      if (!(is.null(cluster))) {
        # determine number of required dummies
        targetCluster <- which(table(cluster) > 1); targetCluster  # no cluster if only one study is included
        targetCluster <- names(targetCluster); targetCluster
        clusCounter <- length(targetCluster); clusCounter
        #tmp1 <- which(cluster %in% targetCluster); tmp1

        # create dummies
        tmpTI <- matrix(0, dim(dataTmp)[1], clusCounter)
        for (i in 1:clusCounter) {
          targetGroups <- which(cluster == targetCluster[i]); targetGroups
          tmp2 <- which(groups %in% targetGroups); length(tmp2)
          tmpTI[tmp2, i] <- 1
        }
        if (CoTiMAStanctArgs$scaleClus == TRUE) tmpTI[ , 1:ncol(tmpTI)] <- scale(tmpTI[ , 1:ncol(tmpTI)])
        currentStartNumber <- as.numeric(gsub("TI", "", colnames(dataTmp)[max(grep("TI", colnames(dataTmp)))]))+1; currentStartNumber
        currentEndNumber <- currentStartNumber + clusCounter -1; currentEndNumber
        colnames(tmpTI) <- paste0("TI", currentStartNumber:currentEndNumber); colnames(tmpTI)
        dataTmp <- cbind(dataTmp, tmpTI); dim(dataTmp)
        #head(dataTmp)
      }

      if (mod.type == "cont") tmp1 <- n.moderators
      if (mod.type == "cat") tmp1 <- catCounter
      if (!(is.null(cluster))) tmp1 <- tmp1 + clusCounter
      dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.latent, n.TIpred = (n.studies-1+tmp1),
                               manifestNames=manifestNames)

      dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
      dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
      # eliminate rows where ALL latents are NA
      dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
      datalong <- dataTmp3

      if (CoTiMAStanctArgs$scaleLongData == TRUE) {
        tmp <- grep ("TI", colnames(datalong)); tmp
        sampleCols <- tmp[1:(n.studies-1)]; sampleCols
        modCols <- tmp[!(tmp %in% sampleCols)]; modCols
        if (CoTiMAStanctArgs$scaleTI == TRUE) datalong[ , sampleCols] <- scale(datalong[ , sampleCols])
        if (CoTiMAStanctArgs$scaleMod == TRUE) datalongl[ , modCols] <- scale(datalong[ , modCols])
      }

    }

    #######################################################################################################################
    ########################################### Model with moderator effects  #############################################
    #######################################################################################################################
    {
      print(paste0("#################################################################################"))
      print(paste0("######## Fit Moderator Model with the following Drift Effects Moderated: ########"))
      print(paste0("### ",moderatedDrift, " ###"))
      print(paste0("#################################################################################"))

      # Make model with most time points
      n.all.moderators <- length(colnames(datalong)[grep("TI", colnames(datalong))])-n.studies+1; n.all.moderators
      stanctModModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                                DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent, byrow=TRUE),
                                LAMBDA=diag(n.latent),
                                CINT=matrix(0, nrow=n.latent, ncol=1),
                                T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                                MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent),
                                MANIFESTTRAITVAR = 'auto',
                                type = 'stanct',
                                n.TIpred = (n.studies-1+n.all.moderators),
                                TIpredNames = paste0("TI", 1:(n.studies-1+n.all.moderators)),
                                TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+n.all.moderators)))
      stanctModModel$pars[stanctModModel$pars$matrix %in% 'DRIFT',paste0(stanctModModel$TIpredNames,'_effect')] <- TRUE
      stanctModModel$pars[stanctModModel$pars$matrix %in% 'T0MEANS',paste0(stanctModModel$TIpredNames,'_effect')] <- FALSE
      stanctModModel$pars[stanctModModel$pars$matrix %in% 'LAMBDA',paste0(stanctModModel$TIpredNames,'_effect')] <- FALSE
      stanctModModel$pars[stanctModModel$pars$matrix %in% 'CINT',paste0(stanctModModel$TIpredNames,'_effect')] <- FALSE
      stanctModModel$pars[stanctModModel$pars$matrix %in% 'MANIFESTMEANS',paste0(stanctModModel$TIpredNames,'_effect')] <- FALSE
      stanctModModel$pars[stanctModModel$pars$matrix %in% 'MANIFESTVAR',paste0(stanctModModel$TIpredNames,'_effect')] <- FALSE

      # target effects (to be aggregated):
      targetNames <- driftNames; targetNames
      tmp1 <- which(stanctModModel$pars$matrix == "DRIFT"); tmp1
      tmp2 <- which(stanctModModel$pars[tmp1, "param"] %in% targetNames); tmp2
      stanctModModel$pars[tmp1[tmp2], paste0(stanctModModel$TIpredNames[1:(n.studies-1)],'_effect')] <- FALSE
      tmp2 <- which(!(stanctModModel$pars[tmp1, "param"] %in% moderatedDrift)); tmp2
      stanctModModel$pars[tmp1[tmp2], paste0(stanctModModel$TIpredNames,'_effect')] <- FALSE


      # DIFFUSION effects (not to be moderated):
      tmp1 <- which(stanctModModel$pars$matrix == "DIFFUSION"); tmp1
      stanctModModel$pars[tmp1, paste0(stanctModModel$TIpredNames[(n.studies):(n.studies+n.all.moderators-1)],'_effect')] <- FALSE

      # T0var effects (not to be moderated):
      tmp1 <- which(stanctModModel$pars$matrix == "T0VAR"); tmp1
      stanctModModel$pars[tmp1, paste0(stanctModModel$TIpredNames[(n.studies):(n.studies+n.all.moderators-1)],'_effect')] <- FALSE

      #stanctModModel$pars[, 1:4]

      fitStanctModModel <- ctsem::ctStanFit(
        datalong = datalong,
        ctstanmodel = stanctModModel,
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
        control=CoTiMAStanctArgs$control,
        verbose=CoTiMAStanctArgs$verbose,
        cores=coresToUse)


      if (!(is.null(CoTiMAStanctArgs$resample))) {
        fitStanctModModel <- ctmaStanResample(ctmaFittedModel=fitStanctModModel) #, CoTiMAStanctArgs=CoTiMAStanctArgs)
      }

      stanctModFit <- summary(fitStanctModModel, digits=2*digits, parmatrices=TRUE, residualcov=FALSE)

    }

    # Extract estimates & statistics
    tmp1 <- which(rownames(stanctModFit$parmatrices) %in% c("T0VAR")); tmp1
    tmp2 <- which(rownames(stanctModFit$parmatrices) %in% c("DRIFT")); tmp2
    tmp2 <- c(matrix(tmp2, n.latent, byrow=TRUE)); tmp2 # correct order
    tmp3 <- which(rownames(stanctModFit$parmatrices) %in% c("DIFFUSIONcov")); tmp3
    targetRows <- c(tmp1, tmp2, tmp3); targetRows
    Tvalues <- stanctModFit$parmatrices[targetRows,3]/stanctModFit$parmatrices[targetRows, 4]; Tvalues
    modDrift_Coeff <- round(cbind(stanctModFit$parmatrices[targetRows,], Tvalues), digits); modDrift_Coeff
    # relabel rownames
    for (i in 1:dim(modDrift_Coeff)[1]) {
      rownames(modDrift_Coeff)[i] <- paste0(rownames(modDrift_Coeff)[i], "_V",
                                            modDrift_Coeff[i, 2], "toV", modDrift_Coeff[i, 1])
    }
    #modDrift_Coeff

    ## moderator effects
    tmp1 <- tmp2 <- c()
    if (mod.type == "cont") tmp2 <- length(mod.number)-1
    if (mod.type == "cat") tmp2 <- length(unlist(unique.mod))-n.moderators

    for (i in modTIstartNum:(modTIstartNum+tmp2)) tmp1 <- c(tmp1, (grep(i, rownames(stanctModFit$tipreds))))
    Tvalues <- stanctModFit$tipreds[tmp1, ][,6]; Tvalues
    modTI_Coeff <- round(cbind(stanctModFit$tipreds[tmp1, ], Tvalues), digits); modTI_Coeff
    # re-label
    if (!(is.null(mod.names))) {
      if (mod.type == "cont") {
        for (i in 1:n.moderators) {
          targetNamePart <- paste0("tip_TI", n.studies+i-1); targetNamePart
          rownames(modTI_Coeff) <- sub(targetNamePart, paste0(mod.names[i], "_on_"), rownames(modTI_Coeff))
        }
      }
      if (mod.type == "cat") {
        counter <- 0
        modNameCounter <- 1
        for (j in 1:n.moderators) {
          if (n.moderators == 1) unique.mod.tmp <- unique.mod else unique.mod.tmp <- unique.mod[[j]]
          for (i in 1:(length(unique.mod.tmp)-1)) {
            counter <- counter + 1
            current.mod.names <- mod.names[modNameCounter]; current.mod.names
            targetNamePart <- paste0("tip_TI", n.studies+i-1); targetNamePart
            tmp1 <- grep(targetNamePart, rownames(modTI_Coeff)); tmp1
            rownames(modTI_Coeff) <- sub(targetNamePart, paste0(counter+1, ". smallest value (categorie) of ", mod.names[modNameCounter], "_on"), rownames(modTI_Coeff))
            #modTI_Coeff
          }
          counter <- 0
          modNameCounter <- modNameCounter + 1
        }
      }
    }

    ## cluster effects
    if (!(is.null(cluster))) {
      tmp1 <- c()
      for (i in (modTIstartNum+length(mod.number)):((modTIstartNum+length(mod.number))+clusCounter-1)) tmp1 <- c(tmp1, (grep(i, rownames(stanctModFit$tipreds))))
      Tvalues <- stanctModFit$tipreds[tmp1, ][,6]; Tvalues
      clusTI_Coeff <- round(cbind(stanctModFit$tipreds[tmp1, ], Tvalues), digits); clusTI_Coeff
      # re-label
      for (i in 1:clusCounter) {
        targetNamePart <- paste0("tip_TI", n.studies+i+n.moderators-1); targetNamePart
        newNamePart <- paste0(targetCluster[i], "_on_"); newNamePart
        rownames(clusTI_Coeff) <- sub(targetNamePart, paste0(targetCluster[i], "_on_"), rownames(clusTI_Coeff))
      }
    } else {
      clusTI_Coeff <- NULL
    }

    modDrift_Minus2LogLikelihood  <- -2*stanctModFit$loglik; modDrift_Minus2LogLikelihood
    modDrift_estimatedParameters  <- stanctModFit$npars; modDrift_estimatedParameters
    modDrift_df <- NULL
    #modDrift_Coeff
    #model_Drift_Coef <- modDrift_Coeff["DRIFT" %in% rownames(modDrift_Coeff)), 3]; model_Drift_Coef
    model_Drift_Coef <- modDrift_Coeff[grep("DRIFT", rownames(modDrift_Coeff)), 3]; model_Drift_Coef
    # re-sort
    #model_Drift_Coef <- c(matrix(model_Drift_Coef, n.latent, n.latent, byrow=FALSE)); model_Drift_Coef
    names(model_Drift_Coef) <- driftNames; model_Drift_Coef

    #model_Diffusion_Coef <- modDrift_Coeff[(rownames(modDrift_Coeff) == "DIFFUSIONcov"), 3]; model_Diffusion_Coef
    model_Diffusion_Coef <- modDrift_Coeff[grep("DIFFUSIONcov", rownames(modDrift_Coeff)), 3]; model_Diffusion_Coef
    model_Diffusion_Coef <- c(OpenMx::vech2full(model_Diffusion_Coef)); model_Diffusion_Coef
    names(model_Diffusion_Coef) <- driftNames; model_Diffusion_Coef

    model_T0var_Coef <- modDrift_Coeff[grep("T0VAR", rownames(modDrift_Coeff)), 3]; model_T0var_Coef
    model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
    names(model_T0var_Coef) <- driftNames; model_Diffusion_Coef

    allResults <- list(estimates=modDrift_Coeff, Minus2LL=modDrift_Minus2LogLikelihood,
                       n.parameters=modDrift_estimatedParameters, df=modDrift_df,
                       mod.effects=modTI_Coeff, clus.effects=clusTI_Coeff)

    model_Cint_Coef <- NULL

    #######################################################################################################################
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    st <- paste0("Computation started at: ", start.time); st
    et <- paste0("Computation ended at: ", end.time); et
    tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt

    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

    tmp1 <- grep("CINT", rownames(stanctModFit$parmatrices)); tmp1
    tmp2 <- grep("asym", rownames(stanctModFit$parmatrices)); tmp2
    tmp3 <- grep("dt", rownames(stanctModFit$parmatrices)); tmp3
    tmp4 <- tmp1[(tmp1 %in% c(tmp2, tmp3)) == FALSE]; tmp4
    model_Cint_Coef <- stanctModFit$parmatrices[tmp4, 3]; model_Cint_Coef


    results <- list(activeDirectory=activeDirectory,
                    time=list(start.time=start.time, end.time=end.time, time.taken=time.taken),
                    plot.type="drift", #model.type="stanct",
                    coresToUse=coresToUse, n.studies=n.studies,
                    n.latent=n.latent, n.moderators=length(mod.number), mod.names=mod.names, mod.type=mod.type,
                    studyList=ctmaInitFit$studyList, studyFitList=fitStanctModModel,
                    data=datalong, statisticsList=ctmaInitFit$statisticsList,
                    modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef,
                                      CINT=model_Cint_Coef, MOD=modTI_Coeff),
                    parameterNames=ctmaInitFit$parameterNames,
                    summary=list(model="All Drift Effects Moderated",
                                 estimates=modDrift_Coeff,
                                 randomEffects=stanctModFit$popsd,
                                 minus2ll= modDrift_Minus2LogLikelihood,
                                 n.parameters = modDrift_estimatedParameters,
                                 df=modDrift_df,
                                 mod.effects=modTI_Coeff,
                                 clus.effects=clusTI_Coeff))

    class(results) <- "CoTiMAFit"

    invisible(results)

  } ### END function definition
