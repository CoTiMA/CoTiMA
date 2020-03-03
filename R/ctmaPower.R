#######################################################################################################################
############################################ CoTiMA Statistical Power #################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   useMBESS=TRUE
#   activeDirectory="/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/ARTIKEL/COTIMA BURNOUT/NEU 2019/PB SUBMISSION (BUL-2019-1709)/REVISION/_TRYOUT/"
#   resultsFilePrefix="ctmaPower"
#   saveFilePrefix="ctmaPower"
#   CoTiMAInitFit=CoTiMAInitFit1
#   datawide_all=datawide_all
#   activateRPB=FALSE
#   silentOverwrite=FALSE
#   allInvStartValues=NULL
#   compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2]
#   loadAllInvFit=c()
#   saveAllInvFit=c("ctmaAllInvMx")
#   loadAllInvWOSingFit=c()
#   saveAllInvWOSingFit=c("ctmaAllInvWOSingMx")
#   statisticalPower=c(.90, .80)                   # beta, e.g., 95. Computes required sample sizes for discrete time drift paramters across range of time lags. Sets testAllInvariantModel=TRUE
#   ##### ENTER DEBUG INFO HERE #######
#   silentOverwrite=FALSE
#   ctmaInitFit=ctmaInitFit1
#   timeRange=NULL
# }
# debug <- 0

#' ctmaPower
#'
#' @param ctmaInitFit ?
#' @param activeDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param statisticalPower ?
#' @param allInvStartValues ?
#' @param timeRange ?
#' @param compSVMethod ?
#' @param useMBESS ?
#' @param type ?
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param fullDriftStartValues ?
#' @param activateRPB ?
#' @param silentOverwrite ?
#' @param digits ?
#' @param saveStatPower ?
#' @param loadStatPower ?
#' @param loadAllInvFit ?
#' @param saveAllInvFit ?
#' @param loadAllInvWOSingFit ?
#' @param saveAllInvWOSingFit ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaPower <- function(
  # Primary Study Fits
  ctmaInitFit=NULL,                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,
  #sourceDirectory= NULL,
  resultsFilePrefix="CoTiMBbias",
  saveFilePrefix="ctmaBias",

  statisticalPower=c(),
  allInvStartValues=NULL,
  timeRange=NULL,                     # e.g., seq(0, 120, 1)
  compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  useMBESS=FALSE,

  # Fitting Parameters
  type="mx",
  retryattempts=30,
  refits=1,
  NPSOL=FALSE,
  coresToUse=c(-1),
  fullDriftStartValues=NULL,


  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  silentOverwrite=FALSE,

  digits=4,

  saveStatPower="statPowerMx",
  loadStatPower=NULL,
  loadAllInvFit=c(),
  saveAllInvFit=c("ctmaAllInvMx"),
  loadAllInvWOSingFit=c(),
  saveAllInvWOSingFit=c("ctmaAllInvWOSingMx")

)

{  # begin function definition (until end of file)

  { ### CHECKS
    if (is.null(activeDirectory)) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("No working directory has been specified!", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat("Press 'q' to quit and change or 'c' to continue with the active directory of the ctmaInitFit file. Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
      activeDirectory <- ctmaInitFit$activeDirectory
    }

    # check if fit object is specified
    if (is.null(ctmaInitFit)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("A fitted CoTiMA object (\"ctmaInitFit\") has to be supplied to analyse something. \n"))
      stop("Good luck for the next try!")
    }

    if (length(statisticalPower) < 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold(" ", sep="\n"))
      cat(crayon::red$bold("Levels (betas) for statistical power (\"statisticalPower\") have not been suppiled \n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat("Press 'q' to quit and change or 'c' to continue with the statisticalPower=c(.90, .80). Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
      statisticalPower <- c(.90, .80)
    }

    if (length(statisticalPower) > 0) {
      for (i in 1:length(statisticalPower)) {
        if (statisticalPower[i] < 0 | statisticalPower[i] > 1){
          if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          cat(crayon::red$bold("Parameter for statistical poweranalysis outside the allowed interval!\n"))
          cat(crayon::red("Values have to be between 0 and 1\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    if (resultsFilePrefix=="ctmaPower") {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat("The default results file prefix (ctmaPower) has been chosen.", "\n")
      cat("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards ", "\n")
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }
  } # END Checks


  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  {
    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  {
        if (.Platform$OS.type == "unix") {
          coresToUse <- parallel::detectCores() + coresToUse
          #coresToUse <- detectCores() + coresToUse
          cat("     #################################################################################\n")
          cat("     # You decided using multiple CPU cores for parallel fitting CoTiMA models. This #\n")
          cat("     # may work well in many cases, but we highly recommend not to run this script   #\n")
          cat("     # in RStudio. Instead, go to the console (e.g., console.app on MAC OS), then    #\n")
          cat("     # change to the direcrtory where your R-script is located (by typing, e.g.,     #\n")
          cat("     # \'cd \"/Users/trump/only/temp\"\'), and then, e.g., \'Rscript \"CoTiMA1.R\"\'.#\n")
          cat("     #################################################################################")
          Sys.sleep(1)
        } else {
          coresToUse <- 1
        }
      }
    } else {
      coresToUse <- 1
    }
    coresToUse
    #CoTiMAStanctArgs$cores <- coresToUse

    if ((-1) * coresToUse >= parallel::detectCores()) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- parallel::detectCores() - 1
      cat(crayon::red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
    }

    numOfThreads <- round((coresToUse/refits - .4), 0); numOfThreads
    if(numOfThreads < 1) numOfThreads <- 1
    tmp <- paste0("No of Threads set to ", numOfThreads); tmp
    cat("","\n")
    cat(tmp,"\n")
  }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {
    tmp <- ctmaCombPRaaw(listOfStudyFits=ctmaInitFit$studyFitList)
    datawide_all <- tmp$alldata
    groups <- tmp$groups
    names(groups) <- c("Study_No_"); groups
    groupsNamed <- (paste0("Study_No_", groups)); groupsNamed

    n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    all_Coeff <- (lapply(ctmaInitFit$studyFitList, function(extract) extract$mxobj$output$estimate)); all_Coeff
    all_SE <- (lapply(ctmaInitFit$studyFitList, function(extract) extract$mxobj$output$standardErrors)); all_SE
    allSampleSizes <- unlist(lapply(ctmaInitFit$studyList, function(extract) extract$sampleSize)); allSampleSizes
    allSampleSizes <- allSampleSizes[-length(allSampleSizes)]; allSampleSizes
    tmp1 <- which(ctmaInitFit$statisticsList$allTpoints == max(ctmaInitFit$statisticsList$allTpoints)); tmp1
    ctsemModel <- ctmaInitFit$studyFitList[[tmp1]]$ctmodelobj
    if (is.null(timeRange)) usedTimeRange <- seq(0, 1.5*maxDelta, 1) else usedTimeRange <- timeRange
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas); maxDelta
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints
    driftNames <- ctmaInitFit$parameterNames$DRIFT; driftNames

    stepWidth <- max(usedTimeRange)/(length(usedTimeRange)-1)

  }



  #######################################################################################################################
  ##################################### Analyses of Statistical Power ###################################################
  #######################################################################################################################

  print(paste0("#################################################################################"))
  print(paste0("########################## Statistical Power Analysis ###########################"))
  print(paste0("#################################################################################"))

  print(paste0("#################################################################################"))
  print(paste0("######## Fitting all fixed CoTiMA - ALL parameters equal across groups ##########"))
  print(paste0("#################################################################################"))

  # make models
  allFixedModel <- ctsemModel
  allFixedModel$DRIFT <- matrix("groupfixed", n.latent, n.latent)
  allFixedModel$DIFFUSION[grep("diffusion", allFixedModel$DIFFUSION)] <- "groupfixed"
  allFixedModel$T0VAR[grep("T0var", allFixedModel$T0VAR)] <- "groupfixed"
  if (any(ctsemModel$CINT != 0) ) {
    allFixedModel$CINT <- matrix("groupfixed", n.latent, 1)
  } else {
    allFixedModel$CINT <- matrix(0, n.latent, 1)
  }
  ctsemModel <- ctsemModel
  allFixedModel <- allFixedModel
  modModel <- NULL
  if (!(is.null(allInvStartValues))) {
    startValues <- allInvStartValues
  } else {
    startValues <- ctmaCompSV(singleStudyFits = ctmaInitFit$studyFitList,
                          ctModelObj = ctsemModel, fixedModel = allFixedModel,
                          modModel = modModel, moderatorValues = NULL,
                          compSVMethod=compSVMethod)
  }

  # FIT (NOT LOAD)
  OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
  results <- parallel::mclapply(seq(1, refits, by=1),
                      function(x) ctsem::ctFit(dat=datawide_all,
                                        retryattempts=retryattempts,
                                        ctmodelobj = ctsemModel,
                                        omxStartValues=startValues),
                      mc.cores=coresToUse)
  OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
  # Select model with best fit
  allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
  homAllFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
  homAllFit <- homAllFit$mxobj

  # LOAD
  if (length(loadAllInvFit) > 0) {
    x1 <- paste0(activeDirectory, loadAllInvFit[1], " homAllFit.rds"); x1
    homAllFit <- readRDS(file=x1)
  }

  # SAVE
  if (length(saveAllInvFit) > 0)  {
    x1 <- paste0(saveAllInvFit[1], "Fit.rds"); x1
    x2 <- paste0(activeDirectory); x2
    ctmaSaveFile(activateRPB, "", homAllFit, x1, x2, silentOverwrite=silentOverwrite)
  }

  ### Extract estimates & statistics
  allStandardErrors <- sqrt(2*(diag(solve(homAllFit$output$calculatedHessian)))); allStandardErrors
  ## Drift
  driftRows <- seq(1, n.latent^2, 1); driftRows
  homAll_Drift_Coef <- homAllFit$output$estimate[driftRows]; homAll_Drift_Coef
  homAll_Drift_SE <- homAllFit$output$standardErrors[driftRows]; homAll_Drift_SE
  homAll_Drift_Tvalue <- (homAll_Drift_Coef/homAll_Drift_SE); homAll_Drift_Tvalue
  ## Diffusion (there are 4 'Normal' coefficients, which are 3 (lower triagonal) transformed base coefficients)
  diffusionRows <- (max(driftRows)+1):(max(driftRows)+(n.latent*(n.latent+1)/2)); diffusionRows
  homAll_Diffusion_Coef <- c(homAllFit$output$algebras$ctsem.DIFFUSION); homAll_Diffusion_Coef
  homAll_Diffusion_SE <- c(OpenMx::vech2full(allStandardErrors[diffusionRows])); homAll_Diffusion_SE
  homAll_Diffusion_Tvalue <- (homAll_Diffusion_Coef/homAll_Diffusion_SE); homAll_Diffusion_Tvalue
  ## T0Var (there are 4 'Normal' coefficients, which are 3 (lower triagonal) transformed base coefficients)
  T0VarRows <- (max(diffusionRows)+1):(max(diffusionRows)+(n.latent*(n.latent+1)/2)); T0VarRows
  # Select the non-transformed ctparameters
  homAll_T0Var_Coef <- c(homAllFit$output$algebras$ctsem.T0VAR); homAll_T0Var_Coef
  homAll_T0Var_SE <- c(OpenMx::vech2full(allStandardErrors[T0VarRows])); homAll_T0Var_SE
  homAll_T0Var_Tvalue <- (homAll_T0Var_Coef/homAll_T0Var_SE); homAll_T0Var_Tvalue
  ## Extract Model Fit
  homAll_Minus2LogLikelihood  <- homAllFit$output$Minus2LogLikelihood; homAll_Minus2LogLikelihood
  homAll_estimatedParameters  <- length(homAllFit$output$estimate); homAll_estimatedParameters
  homAll_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
    homAll_estimatedParameters; homAll_df

  # Combine summary information
  homAll_effects <- matrix(t(cbind((homAll_Drift_Coef), (homAll_Drift_SE),
                                   (homAll_Drift_Tvalue))), 1, 3*length(driftRows), byrow=T); homAll_effects
  homAll_effects <- rbind(homAll_effects,
                          matrix(t(cbind((homAll_Diffusion_Coef), (homAll_Diffusion_SE),
                                         (homAll_Diffusion_Tvalue) )), 1, 3*length(driftRows), byrow=T)); homAll_effects
  homAll_effects <- rbind(homAll_effects,
                          matrix(t(cbind((homAll_T0Var_Coef), (homAll_T0Var_SE),
                                         (homAll_T0Var_Tvalue) )), 1, 3*length(driftRows), byrow=T)); homAll_effects
  # Label summary table
  rownames(homAll_effects) <- c("Fixed Effects Drift", "Fixed Effects Diffusion", "Fixed Effects T0Var")
  newColNames <- c()
  for (j in 1:n.latent) {
    for (h in 1:n.latent) {
      newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)", "Tvalue")
    }
  }
  colnames(homAll_effects) <- newColNames; homAll_effects

  DRIFT_Coef <- matrix(homAll_Drift_Coef, n.latent, n.latent); DRIFT_Coef
  DIFFUSION_Coef <- matrix(homAll_T0Var_Coef, n.latent, n.latent); DIFFUSION_Coef
  T0VAR_Coef <- matrix(homAll_T0Var_Coef, n.latent, n.latent); T0VAR_Coef

  print(paste0("#################################################################################"))
  print(paste0("############ Fitting model with single cross effects fixed to 0.0 ###############"))
  print(paste0("#################################################################################"))

  # Make series of models. In each one there is one drift effects set to zero
  # ctmaCompSV is not used because effects are extracted from previous hom model
  fullWOCross <- list()
  counter <- 0
  for (j in 1:n.latent) {
    for (h in 1:n.latent) {
      counter <- counter +1
      fullWOCross[[counter]] <- ctsemModel
      fullWOCross[[counter]]$DRIFT[h,j] <- 0
    }
  }

  # fit models
  fullWOSingleFit <- list()
  fullWOSingle_Diffusion_Coef  <- list()
  fullWOSingle_Diffusion_SE <- list()
  fullWOSingle_Diffusion_Tvalue <- list()
  fullWOSingle_Diffusion_Coef_Normal <- list()
  fullWOSingle_T0Var_Coef  <- list()
  fullWOSingle_T0Var_SE <- list()
  fullWOSingle_T0Var_Tvalue <- list()
  fullWOSingle_T0Var_Coef_Normal <- list()
  fullWOSingle_Drift_Coef <- list()
  fullWOSingle_Drift_SE <- list()
  fullWOSingle_Drift_Tvalue <- list()

  DRIFT.without.j <- list()             # coefficient in matrix formant with missing coeff set to 0 (used in subsequent power computations)
  DIFFUSION.without.j <- list()         # coefficient in matrix formant (used in subsequent power computations)
  T0VAR.without.j <- list()             # coefficient in matrix formant (used in subsequent power computations)

  driftRows <- 1:(n.latent^2-1); driftRows
  diffusionRows <- (max(driftRows)+1):(max(driftRows)+(n.latent*(n.latent+1)/2)); diffusionRows
  T0VarRows <- (max(diffusionRows)+1):(max(diffusionRows)+(n.latent*(n.latent+1)/2)); T0VarRows

  counter <- 0
  allMinus2LogLikelihood <- fullWOSingleFit <- list()

  # Sequence of cross effects in drift matrices (only model without cross effects are computed to save time)
  coeffSeq <- seq(1, n.latent^2, 1)[!(seq(1, n.latent^2, 1) %in% seq(1, n.latent^2, (n.latent+1)))]; coeffSeq

  for (j in 1:(n.latent^2)) {
    counter <- counter + 1
    if (j %in% coeffSeq) {
      newStartValues <- OpenMx::omxGetParameters(homAllFit)[-counter]; newStartValues
      fullWOCross[[counter]]

      # LOAD
      if (length(loadAllInvWOSingFit) > 0)  {
        x1 <- paste0(loadAllInvWOSingFit[1], " WO", j, ".rds"); x1
        fullWOSingleFit[[j]] <- readRDS(file=x1)
      }

      # FIT (NOT LOAD)
      if (length(loadAllInvWOSingFit) < 1)  {
        OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
        results <- parallel::mclapply(seq(1, refits, by=1),
                            function(refits) ctsem::ctFit(dat=datawide_all, retryattempts = retryattempts,
                                                   omxStartValues=newStartValues,
                                                   ctmodelobj = fullWOCross[[counter]]),
                            mc.cores=coresToUse)
        OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        fullWOSingleFit[[j]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
      }

      # SAVE
      if (length(saveAllInvWOSingFit) > 0)  {
        x1 <- paste0(saveAllInvWOSingFit[1], " WO", j, ".rds"); x1
        x2 <- paste0(activeDirectory); x2
        ctmaSaveFile(activateRPB, "", fullWOSingleFit[[j]], x1, x2, silentOverwrite=silentOverwrite)
      }

      # Collect all results (not really necessary)
      fullWOSingle_Drift_Coef[[j]] <- fullWOSingleFit[[j]]$mxobj$output$estimate[driftRows]; fullWOSingle_Drift_Coef[[j]]
      fullWOSingle_Drift_SE[[j]] <- fullWOSingleFit[[j]]$mxobj$output$standardErrors[driftRows]; fullWOSingle_Drift_SE[[j]]
      fullWOSingle_Drift_Tvalue[[j]] <- (fullWOSingle_Drift_Coef[[j]]/fullWOSingle_Drift_SE[[j]]); fullWOSingle_Drift_Tvalue[[j]]
      fullWOSingle_Diffusion_Coef[[j]] <- fullWOSingleFit[[j]]$mxobj$output$estimate[diffusionRows]; fullWOSingle_Diffusion_Coef[[j]]
      fullWOSingle_Diffusion_SE[[j]] <- fullWOSingleFit[[j]]$mxobj$output$standardErrors[diffusionRows]; fullWOSingle_Diffusion_SE[[j]]
      fullWOSingle_Diffusion_Tvalue[[j]] <- (fullWOSingle_Diffusion_Coef[[j]]/fullWOSingle_Diffusion_SE[[j]]); fullWOSingle_Diffusion_Tvalue[[j]]
      fullWOSingle_Diffusion_Coef_Normal[[j]] <- c(fullWOSingleFit[[j]]$mxobj$algebras$DIFFUSION$result); fullWOSingle_Diffusion_Coef_Normal[[j]]
      fullWOSingle_T0Var_Coef[[j]] <- fullWOSingleFit[[j]]$mxobj$output$estimate[T0VarRows]; fullWOSingle_T0Var_Coef[[j]]
      fullWOSingle_T0Var_SE[[j]] <- fullWOSingleFit[[j]]$mxobj$output$standardErrors[T0VarRows]; fullWOSingle_T0Var_SE[[j]]
      fullWOSingle_T0Var_Tvalue[[j]] <- (fullWOSingle_T0Var_Coef[[j]]/fullWOSingle_T0Var_SE[[j]]); fullWOSingle_T0Var_Tvalue[[j]]
      fullWOSingle_T0Var_Coef_Normal[[j]] <- c(fullWOSingleFit[[j]]$mxobj$algebras$T0VAR$result); fullWOSingle_T0Var_Coef_Normal[[j]]
      # Add one zero in each returned vector of drift coefficents to make a full drift matrix later
      DRIFT.without.j[[j]]<-NA
      DIFFUSION.without.j[[j]]<-NA
      T0VAR.without.j[[j]]<-NA
      counter2 <- 0
      for (h in 1:(n.latent^2)) {
        counter2 <- counter2 + 1; counter2
        DRIFT.without.j[[j]][h] <- fullWOSingle_Drift_Coef[[j]][counter2]
        DIFFUSION.without.j[[j]][h] <- fullWOSingle_Diffusion_Coef_Normal[[j]][h]
        T0VAR.without.j[[j]][h] <- fullWOSingle_T0Var_Coef_Normal[[j]][h]
        if (j == h) {
          DRIFT.without.j[[j]][h] <- 0
          counter2 <- counter2 - 1; counter2
        }
      }

      ### Make drift.without.j, diffusion.without.j & T0VAR.without.j matrices,
      ## DRIFT
      DRIFT.without.j[[j]] <- matrix(DRIFT.without.j[[j]], n.latent, n.latent); DRIFT.without.j[[j]]

      ## DIFFUSION
      DIFFUSION.without.j[[j]] <- matrix(DIFFUSION.without.j[[j]], n.latent, n.latent); DIFFUSION.without.j[[j]]

      ## T0VAR
      T0VAR.without.j[[j]] <- matrix(T0VAR.without.j[[j]], n.latent, n.latent); T0VAR.without.j[[j]]
      T0VAR.without.j
    } # if (j %in% coeffSeq)
  } # for (j in 1:(n.latent^2))

  DRIFT <- matrix(homAll_Drift_Coef, n.latent, n.latent); DRIFT
  DIFFUSION <- matrix(homAll_Diffusion_Coef, n.latent, n.latent); DIFFUSION
  T0VAR <- matrix(homAll_T0Var_Coef, n.latent, n.latent); T0VAR

  print(paste0("#################################################################################"))
  print(paste0("########### Compute required sample size to achieve requested power #############"))
  print(paste0("#################################################################################"))


  # Fast function to calculate required sample sizes later (as optional replacement for MBESS::ss.power.reg.coef)
  nestedProbFunT <- function (fvalue, alpha=.05, power=.80, p=2, x) (1-
                                                                       stats::pt(
                                                                         stats::qt((1 - alpha/2), df = (x)-p-1,
                                                                            lower.tail = TRUE, log.p = FALSE),
                                                                         df = (x)-p-1, ncp = sqrt(x) * abs(fvalue),
                                                                         lower.tail = TRUE, log.p = FALSE)) - power

    # Create table: sampleSizes x deltas (of primary studies) for post hoc power calculations
    tableNxDeltas <- matrix(NA, nrow=n.studies, ncol=maxTpoints); tableNxDeltas
    tableNxDeltas[ ,1]  <- unlist(allSampleSizes); tableNxDeltas
    for (j in 1:n.studies) {
      for (h in 1:length(allDeltas[[j]])) {
        tableNxDeltas[j , (1+h)] <- allDeltas[[j]][h]
      }
    }
    tableNxDeltas[is.na(tableNxDeltas)] <- -99; tableNxDeltas
    tableNxPowerAlpha05 <- tableNxDeltas; tableNxPowerAlpha05
    tableNxPowerAlpha05[ , 2:maxTpoints] <- NA; tableNxPowerAlpha05
    tableNxPowerAlpha01 <- tableNxPowerAlpha05; tableNxPowerAlpha01

    listPowerAlpha05 <- list()
    listPowerAlpha01 <- list()

    # Define functions to compute discrete time effects
    discreteDriftFunction <- function(driftMatrix, timeScale, number) {
      discreteDriftValue <- expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }
    discreteDiffusionFunction <- function(diffusionMatrix, driftMatrix, timeScale, number) {
      driftHatch <- driftMatrix %x% diag(dim(diffusionMatrix)[1]) + diag(dim(diffusionMatrix)[1]) %x% driftMatrix
      discreteDiffusionValue <- solve(driftHatch) %*% (expm(timeScale %x% driftHatch) - diag(dim(driftHatch)[1])) %*% c(diffusionMatrix)
      discreteDiffusionValue[number] }

    # Loop through a range of lags to determine sample sizes (same parameters as for plotting the effects furter below)
    plotPairs <- array(dim=c(n.latent^2, length(statisticalPower), length(usedTimeRange), 2))  # all drift effects, all powers, time range, timePoint+SampleSize

    # Plot required sample size for cross effects.
    for (h in 1:length(statisticalPower)) {
      #h <- 1
      counter <- 0
      for (j1 in 1:(n.latent)) {
        for (j2 in 1:(n.latent)) {
          #j1 <- 1
          #j2 <- 2
          counter <- counter + 1; counter
          if (j1 != j2) {
            j <- counter; j
            for (k in 1:(length(usedTimeRange)-1)) {
              #k <- 1
              delta_t <- usedTimeRange[k+1]; delta_t
              plotPairs[j, h, k, 1] <- usedTimeRange[k+1] # time point

              # R2 in terms of Kelley & Maxwell 2008
              A <- expm(DRIFT %x% delta_t) %*% T0VAR %*% t(expm(DRIFT %x% delta_t)); A         # implied variance at later Tpoint
              S <- matrix(discreteDiffusionFunction(DIFFUSION, DRIFT, delta_t, 1:4),
                          n.latent, n.latent); S   # residual variance at later Tpoint
              R2 <- A[j2,j2]/( (A + S)[j2,j2] ); R2                                            # explained variance at later Tpoint

              # R2 without j (cross effect) in terms of Kelley & Maxwell 2008
              DRIFT.without.j[[counter]]
              A.j <- expm(DRIFT.without.j[[counter]] %x% delta_t) %*%
                T0VAR.without.j[[j]] %*%
                t(expm(DRIFT.without.j[[j]] %x% delta_t)); A.j           # implied variance without j (counter) as predictor at later Tpoint
              S.j <- matrix(
                discreteDiffusionFunction(DIFFUSION.without.j[[j]],      # implied residual variance without j (counter) as predictor at later Tpoint
                                          DRIFT.without.j[[j]],
                                          delta_t, 1:4), n.latent,n.latent); S.j
              R2.j <- A.j[j2,j2]/( (A.j + S.j)[j2,j2] ); R2.j           # explained variance without j (counter) as predictor at later Tpoint

              # Rounding to prevent endless computations
              powerRounding <- 6
              R2 <- round(R2, powerRounding); R2
              R2.j <- round(R2.j, powerRounding); R2.j

              # Skip in case that more variance is explained after predictor is removed
              if (R2 < R2.j) {
                plotPairs[j, h, k, 2] <- 100000
              } else {
                # The following is Kelley's function, which is replaced below by our own, which is > 50 times faster
                # Kelley, K. (2019). The MBESS R Package. R package version 4.6.0. Retrieved from:
                # https://cran.r-project.org/web/packages/MBESS/MBESS.pdf
                if (useMBESS == TRUE) {
                  plotPairs[j, h, k, 2] <- MBESS::ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                                             p = n.latent, desired.power = statisticalPower[h],
                                                             alpha.level = 0.05)[[1]] #
                } else {
                  # The following uses our own function
                  signalToNoiseRatios <- sqrt((R2-R2.j)/(1-R2)); signalToNoiseRatios
                  helper <- round(rootSolve::uniroot.all(nestedProbFunT, c(n.latent+2,999999999),
                                              fvalue=signalToNoiseRatios, alpha=.05,
                                              power=statisticalPower[h], p=n.latent) + .49999, 0)
                  if (length(helper) < 1) helper <- NA
                  plotPairs[j, h, k, 2] <- helper
                }

                # Post hoc power computations
                if ( (delta_t %in% tableNxDeltas[ ,-1]) & (h == 1) ){  # do only once (not for all a priori powers)
                  M <- tableNxDeltas[ ,-1] == delta_t; M # temp matrix used below
                  empiricalN <- matrix(tableNxDeltas[apply(M, 1, any), ], ncol=maxTpoints)[,1]; empiricalN
                  #empiricalN <- matrix(tableNxDeltas[M, ], ncol=maxTpoints)[,1]; empiricalN
                  empiricalN <- na.omit(empiricalN); empiricalN
                  for (l in empiricalN) {
                    p05 <- MBESS::ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                             p = n.latent, Specified.N = l, alpha.level = 0.05)[2]
                    p01 <- MBESS::ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                             p = n.latent, Specified.N = l, alpha.level = 0.01)[2]
                    for (m in 1:n.studies) { # poke power into matrices
                      if (tableNxDeltas[m,1] == l) { # if study has current empirical N ...
                        for (n in 2:(maxTpoints)) {  # ... loop through al lags
                          if (tableNxDeltas[m,n] == delta_t) { # if lag corresponds to current lag ...
                            tableNxPowerAlpha05[m,n] <- as.numeric(p05)  # ... poke
                            tableNxPowerAlpha01[m,n] <- as.numeric(p01)  # ... poke
                          }
                        } # (poke power into matrices) end maxTpoints loops
                      } # (poke power into matrices) end if tableNxDeltas[m,1] == l (delta)
                    } # (poke power into matrices) end n.studies loop
                  } # (compute power if necessary) end empiricalN loop
                } # (compute power if necessary) end delta_t %in% tableNxDeltas[-1,]
              } # (compute power if possible) ende else
            } # end usedTimeRange loop
            if (h == 1) {
              listPowerAlpha05[[counter]] <- tableNxPowerAlpha05
              listPowerAlpha01[[counter]] <- tableNxPowerAlpha01
            }
          } # end j1!=j2 condition
        } # end j2 loop
      } # end j1 loop
      print(paste0("#################################################################################"))
      print(paste0("###################### power calculation for ", statisticalPower[h], " completed ######################"))
      print(paste0("#################################################################################"))
    } # end h loop (length(statisticalPower))

    # Table of required sample sizes for range of different time lags (a priori power)
    requiredSampleSizes <- list()
    currentDriftNames <- c()
    counter1 <- counter2 <- 0
    for (j1 in 1:(n.latent)) {
      for (j2 in 1:(n.latent)) {
        counter1 <- counter1 + 1
        if (j1 != j2 ) {
          counter2 <- counter2 + 1
          requiredSampleSizes[[counter2]] <- plotPairs[counter1, , , 2]
          currentDriftNames <- c(currentDriftNames, driftNames[counter1])
          rowNames  <- plotPairs[counter1, 1, , 1]
        }
      }
    }

    # re-structure into a single table and replace 10000 by NA
    numberOfEffects <- n.latent^2 - n.latent; numberOfEffects
    tmp <- requiredSampleSizes[[1]]; tmp
    for (h in 2:(numberOfEffects)) tmp <- rbind(tmp, requiredSampleSizes[[h]])
    requiredSampleSizes <- t(tmp); requiredSampleSizes
    requiredSampleSizes[requiredSampleSizes==100000] <- NA
    # Label columns
    columnNames <- c()
    for (h in 1:numberOfEffects) {
      for (j in 1:length(statisticalPower)) {
        columnNames <- c(columnNames, paste0(currentDriftNames[h], " Power=", statisticalPower[j]))
      }
    }
    colnames(requiredSampleSizes) <- columnNames
    rownames(requiredSampleSizes) <- round(rowNames, digits)

    # Determine optimal time lag in terms of min samle size required
    rowNames <- c(rownames(requiredSampleSizes), "Min N", "Opt. Lag"); rowNames
    minN <- apply(requiredSampleSizes, 2, min, na.rm=TRUE); minN
    optimalCrossLagN <- c()
    for (h in 1:(dim(requiredSampleSizes)[2])) optimalCrossLagN[h] <- mean(which(requiredSampleSizes[ ,h ] == minN[h]))
    optimalCrossLagN <- optimalCrossLagN*stepWidth; optimalCrossLagN
    requiredSampleSizes <- rbind(requiredSampleSizes, minN, optimalCrossLagN)
    rownames(requiredSampleSizes) <- rowNames

    # Formatting of post hoc results
    postHocPowerList <- list()
    tableNxDeltas[tableNxDeltas == -99] <- NA; tableNxDeltas
    columnNames <- c("N", rep(c("Time Lag", "Power (\u03b1=.05)", "Power (\u03b1=.01)"), (maxTpoints-1))); columnNames

    # Remove empty elements form list and combine
    listPowerAlpha05 <- listPowerAlpha05[!sapply(listPowerAlpha05, is.null)]
    listPowerAlpha01 <- listPowerAlpha01[!sapply(listPowerAlpha01, is.null)]
    for (j in 1:length(currentDriftNames)) {
      tmp05 <- round(listPowerAlpha05[[j]], digits)
      tmp01 <- round(listPowerAlpha01[[j]], digits)
      postHocPower <- tableNxDeltas[, 1]; postHocPower
      for (k in 1:(maxTpoints-1)) { postHocPower <- cbind(postHocPower,
                                                          tableNxDeltas[,k+1],
                                                          tmp05[,k+1],
                                                          tmp01[,k+1])
      }
      colnames(postHocPower) <- columnNames; postHocPower

      # Compute mean & median power (for first lag only) with and without temporarily replacing NA with 0
      tmp <- postHocPower
      tmp[is.na(tmp)] <- 0; tmp
      targetCols <- grep("Power", columnNames); targetCols
      meanPower0 <- apply(tmp[, targetCols[1:2]], 2, mean, na.rm=TRUE); meanPower0
      meanPowerNA <- apply(postHocPower[, targetCols[1:2]], 2, mean, na.rm=TRUE); meanPowerNA
      medianPower0 <- apply(tmp[, targetCols[1:2]], 2, stats::median, na.rm=TRUE); medianPower0
      medianPowerNA <- apply(postHocPower[, targetCols[1:2]], 2, stats::median, na.rm=TRUE); medianPowerNA
      postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
      postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
      postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
      postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower # four times is correct
      postHocPower[dim(postHocPower)[1]-3, targetCols[1:2]] <- round(meanPower0, digits)
      postHocPower[dim(postHocPower)[1]-2, targetCols[1:2]] <- round(meanPowerNA, digits)
      postHocPower[dim(postHocPower)[1]-1, targetCols[1:2]] <- round(medianPower0, digits)
      postHocPower[dim(postHocPower)[1], targetCols[1:2]] <- round(medianPowerNA, digits)
      newNames <- c(paste0("Study_No_", 1:n.studies),
                    "Mean (NA = 0 Power)", "Mean (NA = NA)",
                    "Median (NA = 0 Power)", "Median (NA = NA)")
      rownames(postHocPower) <- newNames
      postHocPowerList[[j]] <- postHocPower
      names(postHocPowerList)[[j]] <- currentDriftNames[j]
    }

  results <- list(activeDirectory=activeDirectory, #sourceDirectory=sourceDirectory,
                  plot.type=c("power"), model.type="mx",
                  coresToUse=NULL, n.studies=1,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=list(homAllFit), #fullWOSingleFit)
                  emprawList=NULL,
                  statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=DRIFT, DIFFUSION=DIFFUSION, T0VAR=T0VAR, CINT=NULL,
                                    DRIFTwoJ=DRIFT.without.j, DIFFUSIONwoJ=DIFFUSION.without.j, T0VARwoJ=T0VAR.without.j),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="Analysis of Statistical Power and Required Sample Sizes",
                               estimates=list("Estimates of Model with all Effects Invariant"=round(homAll_effects, digits),
                                              "Requested Statistical Power"=statisticalPower,
                                              "Power (post hoc) for Drift Effects"=postHocPowerList,
                                              "Required Sample Sizes"=round(requiredSampleSizes, digits))))
  class(results) <- "CoTiMAFit"

  saveRDS(results, paste0(activeDirectory, "ctmaPower_allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
