#######################################################################################################################
############################################## CoTiMA FullDrift MX ####################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   datawide_all=datawide_all
#   groups=groups
#   groupsNamed=groupsNamed
#   activeDirectory=activeDirectory
#   sourceDirectory=sourceDirectory
#   resultsFilePrefix="ctmaFullMx"
#   saveFilePrefix="ctmaFullMx"
#   ctmaInitFit=ctmaInitFit
#   activateRPB=activateRPB
#   checkSingleStudyResults=TRUE
#   digits=4
#   retryattempts=retryattempts
#   refits=refits
#   NPSOL=NPSOL
#   coresToUse=coresToUse
#   fullDriftStartValues=fullDriftStartValues
#   compSVMethod=compSVMethod
#   useCTMultiGroupAlt=useCTMultiGroupAlt
#   confidenceIntervals=confidenceIntervals
#   saveFullDriftModelFit=saveFilePrefix
# }
# debug <- 0

#' ctmaFullMx
#'
#' @param datawide_all ?
#' @param groups ?
#' @param groupsNamed ?
#' @param activeDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param ctmaInitFit ?
#' @param activateRPB ?
#' @param silentOverwrite ?
#' @param checkSingleStudyResults ?
#' @param digits ?
#' @param n.latent ?
#' @param n.studies ?
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param numOfThreads ?
#' @param fullDriftStartValues ?
#' @param compSVMethod ?
#' @param useCTMultiGroupAlt ?
#' @param confidenceIntervals ?
#' @param saveFullDriftModelFit ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaFullMx <- function(
  #type="fullDrift",
  datawide_all=NULL,
  groups=NULL,
  groupsNamed=NULL,
  # Directory names and file names
  activeDirectory = NULL,
  #sourceDirectory = NULL,
  resultsFilePrefix = "ctmaFullMx",
  saveFilePrefix = "ctmaFullMx",
  ctmaInitFit = NULL,
  activateRPB = NULL,
  silentOverwrite = FALSE,
  checkSingleStudyResults = TRUE,
  digits = 4,
  n.latent = NULL,
  n.studies = NULL,
  retryattempts = 1,
  refits = 1,
  NPSOL = FALSE,
  coresToUse = 1,
  numOfThreads = 1,
  fullDriftStartValues = c(),
  compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  useCTMultiGroupAlt = TRUE,
  confidenceIntervals = FALSE,
  saveFullDriftModelFit = saveFilePrefix  )
{  # begin function definition (until end of file)

  # get NPSOL optimizer if requested
  if( NPSOL==TRUE) { # if not already installed but requested get NPSOL
    if ("NPSOL" %in% OpenMx::mxAvailableOptimizers()) {
      OpenMx::mxOption(model=NULL, key="Default optimizer", value="NPSOL", reset = FALSE)
    } else {
      cat(crayon::red$bold("You requested the NPSOL optimizer, which has to be installed manually!!", sep="\n"))
      cat(crayon::red$bold("\n", sep="\n"))
      cat(crayon::red$bold("\n", sep="\n"))
      cat(crayon::red$bold("1. Manually download from https://vipbg.vcu.edu/vipbg/OpenMx2/software/travis/OpenMx_latest.tgz", sep="\n"))
      cat(crayon::red$bold("\n", sep="\n"))
      cat(crayon::red$bold("2. Install.packages(\'/Users/cdormann/Downloads/OpenMx_latest.tar\'", sep="\n"))
      cat(crayon::red$bold("\n", sep="\n"))
      cat(crayon::red$bold("                 lib=\'/Library/Frameworks/R.framework/Versions/3.6/Resources/library\',repos = NULL, type=\"binary\")", "\n"))
      cat(crayon::red$bold("\n", sep="\n"))
      cat(crayon::red$bold("or: install.packages(\"https://vipbg.vcu.edu/vipbg/OpenMx2/software/bin/macosx/travis/OpenMx_latest.tgz\" ", sep="\n"))
      cat(crayon::red$bold("\n", sep="\n"))
      cat(crayon::red$bold("\n", sep="\n"))
    }        # ... version of OpenMx
  } else {
    OpenMx::mxOption(model=NULL, key="Default optimizer", value="CSOLNP", reset = FALSE)
  }

  #######################################################################################################################
  ############################################### attach further packages ###############################################
  #######################################################################################################################
  # {
  #   print(paste0("#################################################################################"))
  #   print(paste0("############################ Attach Further Packages ############################"))
  #   print(paste0("#################################################################################"))
  #
  #   # If OpenMx and ctsem are already attached, detaching them is required to enable OpenMx to run in parallel mode.
  #   #if ("OpenMx" %in% (.packages())) suppressWarnings(detach("package:OpenMx", force=TRUE)) #, unload=TRUE)
  #   #if ("ctsem" %in% (.packages())) suppressWarnings(detach("package:ctsem", force=TRUE)) #, unload=TRUE))
  #   #Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
  #   #library(ctsem)
  #   #OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores()) #now
  #
  #   # Attach all required packages that were not already attach with "library(ctsem)" before
  #   #if (!("MASS" %in% (.packages()))) library(MASS)
  #   #if (!("MBESS" %in% (.packages()))) library(MBESS)
  #   #if (!("rootSolve" %in% (.packages()))) library(rootSolve)
  #   #if (!("doParallel" %in% (.packages()))) library(doParallel)  # this is for parallel processing loops (not for internal parallel processing of OpenMx)
  #   ##if (!("crayon" %in% (.packages()))) library(crayon)
  #   #if (!("psych" %in% (.packages()))) library(psych)
  #
  #   #if (activateRPB==TRUE) {
  #   #  if("RPushbullet" %in% rownames(installed.packages()) == FALSE) {install.packages("RPushbullet")}
  #   #  if (!("RPushbullet" %in% (.packages()))) library(RPushbullet)
  #   #}
  # }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  #####################
  #######################################################################################################################
  {
    start.time <- Sys.time(); start.time
    n.latent <- ctmaInitFit$n.latent
    maxTpoints <- ctmaInitFit$statisticsList$maxTpoints
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas); maxDelta
    usedTimeRange <- seq(0, 1.5*maxDelta, 1)
  }



  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################

  # Make model with most time points
  hetModel <- ctmaInitFit$studyFitList[[1]]$ctmodelobj
  hetModel$Tpoints <- maxTpoints

  # get start values
  homDRIFTallCTmodelobj <- hetModel
  homDRIFTallFixedModel <- hetModel
  homDRIFTallFixedModel$DRIFT <- matrix("groupfixed", n.latent, n.latent)

  homDRIFTallModModel <- NULL
  if (!(is.null(fullDriftStartValues))) {
    homDRIFTallStartValues <- fullDriftStartValues
  } else {
    homDRIFTallStartValues <- ctmaCompSV(singleStudyFits = ctmaInitFit$studyFitList,
                                     ctModelObj = homDRIFTallCTmodelobj, fixedModel = homDRIFTallFixedModel,
                                     modModel = homDRIFTallModModel,
                                     compSVMethod=compSVMethod)
  }

  # FIT
  if (useCTMultiGroupAlt == FALSE) {
    OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
    homDRIFTallStartValues <- ctmaChangeSVLab(startValues = homDRIFTallStartValues, ctmodelobj = homDRIFTallCTmodelobj,
                                          fixedmodel = homDRIFTallFixedModel, n.Studies <- n.studies)
    results <- parallel::mclapply(seq(1, refits, by=1),
                        function(refits) ctsem::ctMultigroupFit(dat=datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                         omxStartValues=homDRIFTallStartValues,
                                                         ctmodelobj = homDRIFTallCTmodelobj,
                                                         fixedmodel = homDRIFTallFixedModel),
                        mc.cores=coresToUse)
    OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
    # Select model with best fit
    allMinus2LogLikelihood <- lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
    homDRIFTallFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
    homDRIFTallFit <- homDRIFTallFit$mxobj
  } else {
    method <- 2
    if (useCTMultiGroupAlt == "1") method <- 1
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("################ Fit CoTiMA refit times & try hard to find optimal fit #################.", "", sep=""))
    cat(paste("", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
    results <- parallel::mclapply(seq(1, refits, by=1),
                        function(refits) ctmaCTMultiGroupAlt(dat=datawide_all, groupings = groups,
                                                            retryattempts = retryattempts,
                                                            startValues = homDRIFTallStartValues,
                                                            ctmodelobj = homDRIFTallCTmodelobj,
                                                            fixedmodel = homDRIFTallFixedModel,
                                                            method = method),
                        mc.cores=coresToUse)
    OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
    # Select model with best fit
    allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
    homDRIFTallFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
    #allEstimatedParameters <- length(homDRIFTallFit$output$estimate); allEstimatedParameters

  } # END else

  if ((silentOverwrite=FALSE)  & (
    any(homDRIFTallFit$output$standardErrors > 1) |
    any(is.na(homDRIFTallFit$output$standardErrors)) |
    any(abs(homDRIFTallFit$output$estimate) > 3) |
    any(is.na(homDRIFTallFit$output$estimate)) ) ) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    print(cbind((homDRIFTallFit$output$estimate), (homDRIFTallFit$output$standardErrors)))
    cat(crayon::red$bold(" ", "", sep="\n"))
    cat(crayon::red$bold("Some estimates or their standard errors seem to be out of range!", sep="\n"))
    cat(crayon::red$bold("Optimal time lags will not be computed if you continue.", sep="\n"))
    cat(crayon::red$bold(" ", "", sep="\n"))
    cat(crayon::blue("Press 'q' to quit and specify or'c'to continue. Press ENTER afterwards ", "\n"))
    char <- readline(" ")
    while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
      cat((crayon::blue("Please press 'q' to quit and 'c' to continue without changes. Press ENTER afterwards.", "\n")))
      char <- readline(" ")
    }
    if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
  }

  # SAVE
  if (length(saveFullDriftModelFit) > 0) {
    x1 <- paste0(saveFullDriftModelFit[1], " ", (Sys.time()), ".rds"); x1
    x2 <- c()
    ctmaSaveFile(activateRPB, activeDirectory, homDRIFTallFit, x1, x2, silentOverwrite=silentOverwrite)
  }

  # Extract estimates & statistics
  driftRows <- seq(1, n.latent^2, 1); driftRows
  homDRIFTall_Drift_Coef <- homDRIFTallFit$output$estimate[driftRows]; homDRIFTall_Drift_Coef
  homDRIFTall_Drift_SE <- homDRIFTallFit$output$standardErrors[driftRows]; homDRIFTall_Drift_SE
  if (is.null(homDRIFTall_Drift_SE)) homDRIFTall_Drift_SE <- 1 # in case of fast fitting no SE are computed
  Tvalue <- (homDRIFTall_Drift_Coef/homDRIFTall_Drift_SE); Tvalue
  if (all(homDRIFTall_Drift_SE==1)) {
    Tvalue <- rep(NA, length(Tvalue)) # in case of fast fitting no SE are computed
    homDRIFTall_Drift_SE <- rep(NA, length(Tvalue))
  }
  homDRIFTall_Minus2LogLikelihood  <- homDRIFTallFit$output$Minus2LogLikelihood; homDRIFTall_Minus2LogLikelihood
  homDRIFTall_estimatedParameters  <- length(homDRIFTallFit$output$estimate); homDRIFTall_estimatedParameters
  homDRIFTall_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
    homDRIFTall_estimatedParameters; homDRIFTall_df


  # Compute confidence intervals
  if (confidenceIntervals == TRUE ) {
    print(paste0("#################################################################################"))
    print(paste0("############### Computing Confidence Intervals for CoTiMA Model  ################"))
    print(paste0("#################################################################################"))

    homDRIFTallFitCI <- homDRIFTallCI <- list()
    tmpModelMxobjFit <- list()
    for (k in 1:(n.latent^2)) tmpModelMxobjFit[[k]] <- homDRIFTallFit  # copy mxobj part of fitted models multiple times (for each drift coefficient)
    ci <- list()
    if (useCTMultiGroupAlt == TRUE) driftNames <- c(tmpModelMxobjFit[[1]]$DRIFT_1$labels) else driftNames <- c(tmpModelMxobjFit[[1]]$Study_No_1$DRIFT$labels)
    for (k in 1:(n.latent^2)) ci[[k]] <- OpenMx::mxCI(driftNames[k]) # make mxCI object for every drift coefficients
    tmpModelMxobj <- list()
    for (k in 1:(n.latent^2)) {
      tmpModelMxobj[[k]] <- OpenMx::mxModel(tmpModelMxobjFit[[k]], ci[[k]]) # make OpenMx Models for all drift coefficients
      #tmpModelMxobj[[k]] <- OpenMx::mxOption(tmpModelMxobj[[k]], "Calculate Hessian", "No")
      #tmpModelMxobj[[k]] <- OpenMx::mxOption(tmpModelMxobj[[k]], "Standard Errors"  , "No")
    }
    OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
    for (i in  1:(n.latent^2)) results[[i]] <- OpenMx::mxRun(tmpModelMxobj[[i]], intervals=TRUE)
    # parallelized version does not work for some reason (rsession have no computing time after the first 5 sec)
    #OpenMx::mxOption(NULL, 'Number of Threads', 1)
    #results <- parallel::mclapply(seq(1, (n.latent^2), by=1),
    #                    function(x) OpenMx::mxRun(tmpModelMxobj[[x]], intervals=TRUE),
    #                    mc.cores=coresToUse)
    #
    #OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())

    for (k in 1:(n.latent^2)) {
      homDRIFTallFitCI[[k]] <- results[[k]]
      homDRIFTallCI[[k]] <- homDRIFTallFitCI[[k]]$output$confidenceIntervals
    }

    # SAVE
    if (length(saveFullDriftModelFit) > 0) {
      for (k in 1:(n.latent^2)) {
        x1 <- paste0(saveFullDriftModelFit[1], " CI_", k, " ", (Sys.time()), ".rds"); x1
        x2 <- c()
        ctmaSaveFile(activateRPB, activeDirectory, homDRIFTallFitCI[k], x1, x2, silentOverwrite=silentOverwrite)

      }
    }
  } ## END confidence intervals

  # Combine summary information
  homDRIFTallDRIFT_effects <- matrix(t(cbind((homDRIFTall_Drift_Coef), (homDRIFTall_Drift_SE), Tvalue)), 1, 3*length(driftRows), byrow=T)

  # Label summary table
  rownames(homDRIFTallDRIFT_effects) <- c("Fixed Effects")
  newColNames <- c()
  for (j in 1:n.latent) {
    for (h in 1:n.latent) {
      newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)", "Tvalue")
    }
  }
  colnames(homDRIFTallDRIFT_effects) <- newColNames; homDRIFTallDRIFT_effects

  if (confidenceIntervals == TRUE) {
    limits <- matrix(NA, 2, (3*(n.latent^2)))
    tmp <- unlist(homDRIFTallCI); tmp
    tmpSeq1 <- seq(1, (3*(n.latent^2)), 3); tmpSeq1 # drift cols & lower limit cols
    tmpSeq2 <- seq(3, (3*(n.latent^2)), 3); tmpSeq2 # upper limit cols
    limits[1, tmpSeq1] <- tmp[tmpSeq2]
    limits[2, tmpSeq1] <- tmp[tmpSeq1]
    homDRIFTallDRIFT_effects <- rbind(homDRIFTallDRIFT_effects, limits); homDRIFTallDRIFT_effects
    rownames(homDRIFTallDRIFT_effects) <- rbind("Fixed Effects", "upper bound",  "lower bound"); homDRIFTallDRIFT_effects
  }


  ### Numerically compute Optimal Time lag sensu Dormann & Griffin (2015)
  driftMatrix <- matrix(homDRIFTall_Drift_Coef, n.latent, n.latent); driftMatrix
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
        targetParameters <- sapply(usedTimeRange, OTL)
        maxCrossEffect[j,h] <- max(abs(targetParameters))
        optimalCrossLag[j,h] <- which(abs(targetParameters)==maxCrossEffect[j,h])*1+0
      }
    }
  }


  #######################################################################################################################
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  st <- paste0("Computation started at: ", start.time); st
  et <- paste0("Computation ended at: ", end.time); et
  tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt


  ###########################################################################################################
  ############################################ SAVE RESULTS #################################################
  ###########################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("################################# Save Results ##################################"))
    print(paste0("#################################################################################"))

    resultsFileName <- paste0(activeDirectory, resultsFilePrefix, " ", Sys.time(), ".txt"); resultsFileName

    sink(file = resultsFileName, append = TRUE, type = c("output"), split = TRUE)

    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("##################### Continuous Time Meta Analysis (CoTiMA)  #########################", "", sep="\n"))
    cat(paste("####################  (cf. Dormann, Guthier, & Cortina, 2019) ##########################", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat(st, sep="\n")
    cat(et, sep="\n")
    cat(tt, sep="\n")
    cat(" ", "", sep="\n")
    cat(paste("########################### contact: cdormann@uni-mainz.de #############################", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("--------------------------------- CoTiMA Parameters ------------------------------------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Fitted model is of type: MX", sep="")
    cat(" ", "", sep="\n")
    cat("Number of latent variables (n.latent): ", n.latent,  sep="")
    cat(" ", "", sep="\n")
    cat("Number of iterations to be used by ctsem (retryattempts): ", retryattempts,  sep="")
    cat(" ", "", sep="\n")
    cat("Number of re-fits used by this CoTiMA Function (refits): ", refits,  sep="")
    cat(" ", "", sep="\n")
    cat("Request NPSOL optimizer (NPSOL): ", NPSOL,  sep="")
    cat(" ", "", sep="\n")
    cat("Allows parallel processing on unix-like computers (e.g., Mac;) using # cores: ", coresToUse,  sep="")
    cat(" ", "", sep="\n")
    cat("Rounding used in output (digits): ", digits,  sep="")
    cat(" ", "", sep="\n")
    cat("Displays estimates from single study ctsem models and waits for user inoput to continue (checkSingleStudyResults): ", checkSingleStudyResults,  sep="")
    cat(" ", "", sep="\n")
    cat("The active directory (activeDirectory) is: ", activeDirectory,  sep="\n")
    cat(" ", "", sep="\n")
    cat("The result file name that is created (resultsFileName): ", resultsFileName,  sep="")
    cat(" ", "", sep="\n")
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("---------- Multi-Sample Homogeneity Model (all Drift Coefficients Invariant) -----------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients and their Standard Errors (SE):")
    cat(" ", "", sep="\n")
    print(round(homDRIFTallDRIFT_effects, digits))
    cat(" ", "", sep="\n")
    cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
    cat("-2 Log Likelihood: ", round(homDRIFTall_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Overall Number of Estimated Parameters: ", round(homDRIFTall_estimatedParameters, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(homDRIFTall_df, digits), sep="")
    cat(" ", "", sep="\n")
    #cat("Comparison with Unconstrained Modell (All Samples Analyzed Separately (~Heterogeneity Model)):")
    #cat(" ", "", sep="\n")
    #cat("Delta(df): ", round(allStudies_homDRIFTall_df, digits), sep="")
    #cat(" ", "", sep="\n")
    #cat("Delta(-2LL) (= Delta(Chi-square)): ", round(allStudies_homDRIFTall_Minus2LogLikelihood, digits), sep="")
    #cat(" ", "", sep="\n")
    #cat("p-value (with double number of digits): ", round(allStudies_homDRIFTall_prob, 2*digits), sep="")
    #cat(" ", "", sep="\n")
    #cat("A sign. value (p < .05) indicates that entire process (i.e., the full Drift Matrix) varies among primary studies.")
    #cat(" ", "", sep="\n")
    cat("The following matrix contains the time lags for which the discrete time cross-lagged effects become largest.", sep="\n")
    cat("(i.e., the optimal time lags sensu Dormann & Griffin, 2015).", sep="\n")
    cat("Note that for autoregressive effects short time lags are always optimal because they steadily decline.")
    cat(" ", "", sep="\n")
    print(optimalCrossLag)
    cat(" ", "", sep="\n")
    cat("The following matrix contains the discrete time cross-lagged effects across the optimal time lags")
    cat(" ", "", sep="\n")
    print(round(maxCrossEffect, digits))
    cat(" ", "", sep="\n")
  }

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n", st, "\n", et, "\nAnalysis successfully completed. \nThank you for using CoTiMA.\nHave a nice day!"))}

  sink()

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}



  if (!(exists("homDRIFTallFitCI"))) homDRIFTallFitCI <- c()

  tmp <- grep("toV", colnames(homDRIFTallDRIFT_effects)); tmp
  model_Drift_Coef <- homDRIFTallDRIFT_effects[1, tmp]; model_Drift_Coef

  tmp <- grep("DIFFUSIONbase", names(homDRIFTallFit$output$matrices)); tmp
  model_Diffusion_Coef <- homDRIFTallFit$output$matrices[tmp]; model_Diffusion_Coef
  model_Diffusion_Coef <- lapply(model_Diffusion_Coef, ctmaBase2Full); model_Diffusion_Coef
  names(model_Diffusion_Coef) <- gsub("base", "", names(model_Diffusion_Coef)); model_Diffusion_Coef

  tmp <- grep("T0VARbase", names(homDRIFTallFit$output$matrices)); tmp
  model_T0var_Coef <- homDRIFTallFit$output$matrices[tmp]; model_T0var_Coef
  model_T0var_Coef <- lapply(model_T0var_Coef, ctmaBase2Full); model_T0var_Coef
  names(model_T0var_Coef) <- gsub("base", "", names(model_T0var_Coef)); model_T0var_Coef

  tmp <- grep("CINT", names(homDRIFTallFit$output$matrices)); tmp
  model_Cint_Coef <- homDRIFTallFit$output$matrices[tmp]; model_Cint_Coef

  results <- list(activeDirectory=activeDirectory, #sourceDirectory=sourceDirectory,
                  plot.type="drift", model.type="mx",
                  coresToUse=coresToUse, n.studies=1,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=list(homDRIFTallFit), # , homDRIFTallFitCI),
                  data=datawide_all, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="all drift fixed (hom. model)",
                               estimates=round(homDRIFTallDRIFT_effects, digits),
                               minus2ll= round(homDRIFTall_Minus2LogLikelihood, digits),
                               n.parameters = round(homDRIFTall_estimatedParameters, digits),
                               df= c(round(homDRIFTall_df, digits)),
                               opt.lag = optimalCrossLag,
                               max.effects = round(maxCrossEffect, digits)))
  class(results) <- "CoTiMAFit"

  saveRDS(results, paste0(activeDirectory, "ctmaFullMx_allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
