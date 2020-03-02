#######################################################################################################################
############################################# CoTiMA SingleDrift MX ###################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   single=single #c("V1toV2", "V2toV2")
#   datawide_all=datawide_all
#   groups=groups
#   groupsNamed=groupsNamed
#   activeDirectory=activeDirectory
#   sourceDirectory=sourceDirectory
#   resultsFilePrefix="ctmaStingleMx"
#   saveFilePrefix="ctmaStingleMx"
#   ctmaInitFit=ctmaInitFit1
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
#   saveSingleDriftModelFit=saveSingleDriftModelFit
#   confidenceIntervals=TRUE
# }
# debug <- 0

#' ctmaSingleMx
#'
#' @param single ?
#' @param datawide_all ?
#' @param groups ?
#' @param groupsNamed ?
#' @param activeDirectory ?
#' @param sourceDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param ctmaInitFit ?
#' @param activateRPB ?
#' @param silentOverwrite ?
#' @param checkSingleStudyResults ?
#' @param digits ?
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param numOfThreads ?
#' @param fullDriftStartValues ?
#' @param compSVMethod ?
#' @param useCTMultiGroupAlt ?
#' @param confidenceIntervals ?
#' @param saveSingleDriftModelFit ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaSingleMx <- function(
  single=c(),
  datawide_all=NULL,
  groups=NULL,
  groupsNamed=NULL,
  # Directory names and file names
  activeDirectory = NULL,
  sourceDirectory = NULL,
  resultsFilePrefix = "ctmaStingleMx",
  saveFilePrefix = "ctmaStingleMx",
  ctmaInitFit = NULL,
  activateRPB = NULL,
  silentOverwrite = FALSE,
  checkSingleStudyResults = TRUE,
  digits = 4,
  #n.latent = NULL,
  #n.studies = NULL,
  retryattempts = 1,
  refits = 1,
  NPSOL = FALSE,
  coresToUse = 1,
  numOfThreads = 1,
  fullDriftStartValues = c(),
  compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  useCTMultiGroupAlt = TRUE,
  confidenceIntervals = FALSE,
  saveSingleDriftModelFit = saveSingleDriftModelFit  )

{  # begin function definition (until end of file)

  # get NPSOL optimizer if requested
  if( NPSOL==TRUE) { # if not already installed but requested get NPSOL
    if ("NPSOL" %in% mxAvailableOptimizers()) {
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
  # #######################################################################################################################
  # {
  #   print(paste0("#################################################################################"))
  #   print(paste0("############################ Attach Further Packages ############################"))
  #   print(paste0("#################################################################################"))

    # If OpenMx and ctsem are already attached, detaching them is required to enable OpenMx to run in parallel mode.
    #if ("OpenMx" %in% (.packages())) suppressWarnings(detach("package:OpenMx", force=TRUE)) #, unload=TRUE)
    #if ("ctsem" %in% (.packages())) suppressWarnings(detach("package:ctsem", force=TRUE)) #, unload=TRUE))
    #Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
    #library(ctsem)
    #OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores()) #now

    # Attach all required packages that were not already attach with "library(ctsem)" before
    #if (!("MASS" %in% (.packages()))) library(MASS)
    #if (!("MBESS" %in% (.packages()))) library(MBESS)
    #if (!("rootSolve" %in% (.packages()))) library(rootSolve)
    #if (!("doParallel" %in% (.packages()))) library(doParallel)  # this is for parallel processing loops (not for internal parallel processing of OpenMx)
    #if (!("crayon" %in% (.packages()))) library(crayon)
    #if (!("psych" %in% (.packages()))) library(psych)

    #if (activateRPB==TRUE) {
    #  if("RPushbullet" %in% rownames(installed.packages()) == FALSE) {install.packages("RPushbullet")}
    #  if (!("RPushbullet" %in% (.packages()))) library(RPushbullet)
#    }
#  }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
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
  ######################### Series of homogeneity models with single drift effects invariant ############################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("##### Fit Series of Homogeneity Models with Single Drift Effects Invariant ######"))
    print(paste0("#################################################################################"))

    homSingleModel <- fixedSingleModel <- homDRIFTSingleFitCI <- fixedSingleModeratorModel <- list()
    homDRIFTSingleFit <- list()
    homDRIFTSingle_Drift_Coef <- homDRIFTSingle_Drift_SE <- list()
    homDRIFTSingle_Minus2LogLikelihood <- homDRIFTSingle_estimatedParameters <- list()
    homDRIFTSingle_df  <- homDRIFTSingleCI <- list()
    homDRIFTSingleDRIFT_effects <- list()
    # Vector of drift coefficients to be tested for invariance
    targetNames <- single; targetNames

    # get start values
    homDRIFTSingleCTmodelobj <- list()
    homDRIFTSingleFixedModel <- list()
    homDRIFTSingleModModel <- NULL
    homDRIFTSingleStartValues <- list()
    #for (k in 1:(n.latent^2)) {
    for (k in 1:length(targetNames)) {
      #k <- 1
      homDRIFTSingleCTmodelobj[[k]] <- ctmaInitFit$studyFitList[[1]]$ctmodelobj
      homDRIFTSingleCTmodelobj[[k]]$Tpoints <- maxTpoints
      homDRIFTSingleFixedModel[[k]] <- ctmaInitFit$studyFitList[[1]]$ctmodelobj
      #homDRIFTSingleFixedModel[[k]]$DRIFT[k] <- "groupfixed"
      homDRIFTSingleFixedModel[[k]]$DRIFT <- sub(targetNames[k], "groupfixed", homDRIFTSingleFixedModel[[k]]$DRIFT)
      homDRIFTSingleStartValues[[k]] <- ctmaCompSV(singleStudyFits = ctmaInitFit$studyFitList,
                                                           ctModelObj = homDRIFTSingleCTmodelobj[[k]],
                                                           fixedModel = homDRIFTSingleFixedModel[[k]],
                                                           modModel = homDRIFTSingleModModel[[k]], moderatorValues = NULL,
                                                           compSVMethod=compSVMethod)
    }

    # Loop through all CROSS coefficients to be tested for invariance
    for (i in 1:length(targetNames)) {
      #i <- 1
      coeffNumber <- i
      print(paste0("#################################################################################"))
      print(paste0("## Fitting model ", i," of ", length(targetNames), " with single drift effects invariant (homDRIFTSingleFit) #"))
      print(paste0("#################################################################################"))

      ## LOAD
      #if (length(loadDRIFTSingleModelFit) > 0 ) {
      #  x1 <- paste0(activeDirectory, loadDRIFTSingleModelFit, " homDRIFTSingleFits/", loadDRIFTSingleModelFit, " homDRIFTSingleFit", i, ".rds"); x1
      #  homDRIFTSingleFit[[i]] <- readRDS(file=x1)
      #  FitList$tmpName <- homDRIFTSingleFit[[i]]
      #  tmp <- names(FitList); tmp
      #  tmp[length(tmp)] <- paste0("homDRIFTSingleFit", i); tmp
      #  names(FitList) <- tmp; names(FitList)
      #}

      # FIT (NOT LOAD)
      #if (length(loadDRIFTSingleModelFit) < 1 ) {

      # FIT (NOT LOAD)
      if (useCTMultiGroupAlt == FALSE) {
        OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
        # change names of start values
        homDRIFTSingleStartValues[[i]] <- changeSVLab(startValues = homDRIFTSingleStartValues[[i]],
                                                                 ctmodelobj = homDRIFTSingleCTmodelobj[[i]],
                                                                 fixedmodel = homDRIFTSingleFixedModel[[i]], noOfStudies <- n.studies)
        results <- parallel::mclapply(seq(1, refits, by=1),
                            function(refits) ctMultigroupFit(dat=datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                             omxStartValues = homDRIFTSingleStartValues[[i]],
                                                             ctmodelobj = homDRIFTSingleCTmodelobj[[i]],
                                                             fixedmodel = homDRIFTSingleFixedModel[[i]]),
                            mc.cores=coresToUse)
        OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        homDRIFTSingleFit[[i]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
        homDRIFTSingleFit[[i]] <- homDRIFTSingleFit[[i]]$mxobj
      } else {
        method <- 2
        if (useCTMultiGroupAlt == "1") method <- 1
        OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
        results <- parallel::mclapply(seq(1, refits, by=1),
                            function(refits) ctmaCTMultiGroupAlt(dat=datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                startValues=homDRIFTSingleStartValues[[i]],
                                                                ctmodelobj = homDRIFTSingleCTmodelobj[[i]],
                                                                fixedmodel=homDRIFTSingleFixedModel[[i]],
                                                                #fastIterativeCoTiMA = fastIterativeCoTiMA,
                                                                method = method,
                                                                tryHard = TRUE),
                            mc.cores=coresToUse)
        OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
        homDRIFTSingleFit[[i]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
      }

      # STORE
      #FitList$tmpName <- homDRIFTSingleFit[[i]]
      #tmp <- names(FitList); tmp
      #tmp[length(tmp)] <- paste0("homDRIFTSingleFit", i); tmp
      #names(FitList) <- tmp; names(FitList)

      #} # END FIT (NOT LOAD)

      # SAVE
      #saveSingleDriftModelFit
      #paste0(saveSingleDriftModelFit[1], " ", targetNames[i], ".rds")
      if (length(saveSingleDriftModelFit) > 0)  {
        # x1 <- paste0(saveSingleDriftModelFit, " SingleDriftFit ", targetNames[i], ".rds"); x1
        x1 <- paste0(saveSingleDriftModelFit[1], "MX ", targetNames[i], " ", (Sys.time()), ".rds"); x1
        #x2 <-paste0(saveSingleDriftModelFit[1], "/"); x2
        x2 <- ""; x2
        ctmaSaveFile(activateRPB, activeDirectory, homDRIFTSingleFit[[i]], x1, x2, silentOverwrite=silentOverwrite)
      }

      # Extract estimates & statistics
      driftRows <- coeffNumber
      homDRIFTSingle_Drift_Coef[[i]] <- homDRIFTSingleFit[[i]]$output$estimate; homDRIFTSingle_Drift_Coef[[i]]
      homDRIFTSingle_Drift_SE[[i]] <- homDRIFTSingleFit[[i]]$output$standardErrors; homDRIFTSingle_Drift_SE[[i]]
      homDRIFTSingle_Minus2LogLikelihood[[i]] <- homDRIFTSingleFit[[i]]$output$Minus2LogLikelihood; homDRIFTSingle_Minus2LogLikelihood[[i]]
      homDRIFTSingle_estimatedParameters[[i]] <- length(homDRIFTSingleFit[[i]]$output$estimate); homDRIFTSingle_estimatedParameters[[i]]
      homDRIFTSingle_df[[i]] <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
        homDRIFTSingle_estimatedParameters[[i]]; homDRIFTSingle_df[[i]]

      # Combine summary information#
      homDRIFTSingleDRIFT_effects[[i]] <- cbind(homDRIFTSingle_Drift_Coef[[i]][i], homDRIFTSingle_Drift_SE[[i]][i])
      Tvalue <- (homDRIFTSingleDRIFT_effects[[i]][1]/homDRIFTSingleDRIFT_effects[[i]][2]); Tvalue

      # Label summary table
      newColNames <- c()
      for (j in 1:n.latent) {
        for (h in 1:n.latent) {
          newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)")
        }
      }
      rownames(homDRIFTSingleDRIFT_effects[[i]]) <- c("Fixed Effects")
      colnames(homDRIFTSingleDRIFT_effects[[i]]) <- newColNames[((driftRows-1)*2+1):((driftRows-1)*2+2)]
      homDRIFTSingleDRIFT_effects[[i]] <- cbind(homDRIFTSingleDRIFT_effects[[i]], Tvalue); homDRIFTSingleDRIFT_effects[[i]]

      # STORE
      #FitList$tmpName <- homDRIFTSingleDRIFT_effects[[i]]
      #tmp <- names(FitList); tmp
      #tmp[length(tmp)] <- paste0("homDRIFTSingleDRIFT_effects", i); tmp
      #names(FitList) <- tmp; names(FitList)

      # Compute confidence intervals
      if (confidenceIntervals == TRUE) {
        print(paste0("#################################################################################"))
        print(paste0(" Computing Confidence Intervalls for Models with Single Drift Effects Invariant  "))
        print(paste0("#################################################################################"))


        # LOAD
        #if (length(loadDRIFTSingleModelFit) > 0) {
        #  x1 <- paste0(activeDirectory, loadDRIFTSingleModelFit, " homDRIFTSingleFits/", loadDRIFTSingleModelFit, " homDRIFTSingleFitCI", i, ".rds"); x1
        #  homDRIFTSingleFitCI[[i]] <- readRDS(file=x1)
        #  homDRIFTSingleCI[[i]] <- homDRIFTSingleFitCI[[i]]$output$confidenceIntervals
        #}

        # FIT (NOT LOAD)
        #if (length(loadDRIFTSingleModelFit) < 1) {
        tmpModelMxobjFit <- homDRIFTSingleFit[[i]]
        #ci <- OpenMx::mxCI(driftNames[i]); ci
        ci <- OpenMx::mxCI(targetNames[i]); ci
        tmpModelMxobj <- OpenMx::mxModel(tmpModelMxobjFit, ci)
        results <- OpenMx::mxRun(tmpModelMxobj, intervals=TRUE)
        homDRIFTSingleFitCI[[i]] <- results
        homDRIFTSingleCI[[i]] <- homDRIFTSingleFitCI[[i]]$output$confidenceIntervals
        #}

        # STORE
        #FitList$tmpName <- homDRIFTSingleFitCI[[i]]
        #tmp <- names(FitList); tmp
        #tmp[length(tmp)] <- paste0("homDRIFTSingleFitCI", i); tmp
        #names(FitList) <- tmp; names(FitList)

        # SAVE
        if (length(saveSingleDriftModelFit) > 0)  {
          x1 <- paste0(saveSingleDriftModelFit[1], "MX CI_", targetNames[i], " ", (Sys.time()), ".rds"); x1
          #x2 <-paste0(saveSingleDriftModelFit[1], "/")
          x2 <- ""
          ctmaSaveFile(activateRPB, activeDirectory, homDRIFTSingleFitCI[[i]], x1, x2, silentOverwrite=silentOverwrite)
        } ## END save fits

      } # END confidence intervals
    } # END for (i in 1:length(targetNames)=

    #} # END testDRIFTSingleModel == TRUE
  } ##

  # Combine Results
  singleDRIFT_effects <- matrix(NA, length(targetNames), 6)
  for (k in 1:length(targetNames)) singleDRIFT_effects[k, ] <- c(homDRIFTSingleDRIFT_effects[[k]],
                                                                 homDRIFTSingle_Minus2LogLikelihood[[k]],
                                                                 homDRIFTSingle_estimatedParameters[[k]],
                                                                 homDRIFTSingle_df[[k]])
  #singleDRIFT_effects
  rownames(singleDRIFT_effects) <- targetNames
  colnames(singleDRIFT_effects) <- c("Estimate", "Std. Err.", "T-value", "-2LL", "n.params", "df")
  if (confidenceIntervals == TRUE) {
    singleDRIFT_effects2 <- matrix(NA, length(targetNames), 2)
    newColnames <- c(colnames(singleDRIFT_effects), "lbound", "ubound")
    for (k in 1:length(targetNames)) singleDRIFT_effects2[k,] <- homDRIFTSingleCI[[k]][c(1,3)]
    singleDRIFT_effects <- cbind(singleDRIFT_effects, singleDRIFT_effects2)
    colnames(singleDRIFT_effects) <- newColnames
  }
  if (confidenceIntervals == TRUE) singleDRIFT_effects <- singleDRIFT_effects[, c(1:3, 7:8, 4:6)] else singleDRIFT_effects <- singleDRIFT_effects
  singleDRIFT_effects


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
    cat(paste("- Multi-Sample Homogeneity Model (only a single Drift Coefficient Invariant per Model) -", "", sep=""))
    cat(sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients, their Standard Errors (SE), T-Values, 2-LL values and df:")
    cat(" ", "", sep="\n")
    print(round(singleDRIFT_effects, digits))
    cat(" ", "", sep="\n")
    cat("Note: If confidence intervals are reported, they do not necessarily have to be symmetric!")
  }

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n", st, "\n", et, "\nAnalysis successfully completed. \nThank you for using CoTiMA.\nHave a nice day!"))}

  sink()

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}



  if (!(exists("homDRIFTallFitCI"))) homDRIFTallFitCI <- c()

  model_Drift_Coef <- model_Diffusion_Coef <- model_T0var_Coef <- model_Cint_Coef <- list()
  for (i in 1:length(targetNames)) {
    tmp <- grep("toV", names(homDRIFTSingleFit[[i]]$output$estimate)); tmp
    model_Drift_Coef[[i]] <- homDRIFTSingleFit[[i]]$output$estimate[tmp]; model_Drift_Coef

    tmp <- grep("diffusion", names(homDRIFTSingleFit[[i]]$output$estimate)); tmp
    model_Diffusion_Coef[[i]] <- homDRIFTSingleFit[[i]]$output$estimate[tmp]; model_Diffusion_Coef

    tmp <- grep("T0var", names(homDRIFTSingleFit[[i]]$output$estimate)); tmp
    model_T0var_Coef[[i]] <- homDRIFTSingleFit[[i]]$output$estimate[tmp]; model_T0var_Coef

    tmp <- grep("CINT", names(homDRIFTSingleFit[[i]]$output$estimate)); tmp
    model_Cint_Coef[[i]] <- homDRIFTSingleFit[[i]]$output$estimate[tmp]; model_Cint_Coef
  }

  if (confidenceIntervals == TRUE) {
    estimates <- singleDRIFT_effects[, 1:5]; estimates
    minus2ll <- singleDRIFT_effects[, 6]; minus2ll
    n.parameters <- singleDRIFT_effects[, 7]; n.parameters
    df <- singleDRIFT_effects[, 8]; df
  } else {
    estimates <- singleDRIFT_effects[, 1:3]; estimates
    minus2ll <- singleDRIFT_effects[, 4]; minus2ll
    n.parameters <- singleDRIFT_effects[, 5]; n.parameters
    df <- singleDRIFT_effects[, 6]; df
  }

  results <- list(activeDirectory=activeDirectory, sourceDirectory=sourceDirectory,
                  plot.type="drift", model.type="mx",
                  coresToUse=coresToUse, n.studies=length(single),
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=homDRIFTSingleFit, # , homDRIFTallFitCI),
                  data=datawide_all, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model=c(paste0(targetNames, " fixed")),
                               estimates=round(estimates, digits), #[]
                               minus2ll= round(minus2ll, digits),
                               n.parameters = round(n.parameters, digits),
                               df= c(round(df, digits))))
  class(results) <- "CoTiMAFit"

  saveRDS(results, paste0(activeDirectory, "ctmaStingleMx_allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
