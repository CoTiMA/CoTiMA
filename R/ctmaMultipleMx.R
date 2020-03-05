#######################################################################################################################
############################################ CoTiMA MultipleDrift MX ##################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   multipleDrift=multipleDrift #c("V1toV2", "V2toV2")
#   datawide_all=datawide_all
#   groups=groups
#   groupsNamed=groupsNamed
#   activeDirectory=activeDirectory
#   sourceDirectory=sourceDirectory
#   resultsFilePrefix="ctmaMultipleMx"
#   saveFilePrefix="ctmaMultipleMx"
#   ctmaInitFit=ctmaInitFit1
#   activateRPB=activateRPB
#   checkMultipleStudyResults=TRUE
#   digits=4
#   retryattempts=retryattempts
#   refits=refits
#   NPSOL=NPSOL
#   coresToUse=coresToUse
#   fullDriftStartValues=fullDriftStartValues
#   compSVMethod=compSVMethod
#   useCTMultiGroupAlt=useCTMultiGroupAlt
#   confidenceIntervals=confidenceIntervals
#   saveMultipleDriftModelFit=saveMultipleDriftModelFit
#   confidenceIntervals=TRUE
# }
# debug <- 0

#' ctmaMultipleMx
#'
#' @param multipleDrift ?
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
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param numOfThreads ?
#' @param fullDriftStartValues ?
#' @param compSVMethod ?
#' @param useCTMultiGroupAlt ?
#' @param confidenceIntervals ?
#' @param saveMultipleDriftModelFit ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaMultipleMx <- function(
  multipleDrift=c(),
  datawide_all=NULL,
  groups=NULL,
  groupsNamed=NULL,
  # Directory names and file names
  activeDirectory = NULL,
  #sourceDirectory = NULL,
  resultsFilePrefix = "ctmaMultipleMx",
  saveFilePrefix = "ctmaMultipleMx",
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
  saveMultipleDriftModelFit = saveMultipleDriftModelFit  )

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
  #{
  #  print(paste0("#################################################################################"))
  #  print(paste0("############################ Attach Further Packages ############################"))
  #  print(paste0("#################################################################################"))
  #
  #  # If OpenMx and ctsem are already attached, detaching them is required to enable OpenMx to run in parallel mode.
  #  #if ("OpenMx" %in% (.packages())) suppressWarnings(detach("package:OpenMx", force=TRUE)) #, unload=TRUE)
  #  #if ("ctsem" %in% (.packages())) suppressWarnings(detach("package:ctsem", force=TRUE)) #, unload=TRUE))
  #  #Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
  #  #library(ctsem)
  #  #OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores()) #now
  #
  #  # Attach all required packages that were not already attach with "library(ctsem)" before
  #  #if (!("MASS" %in% (.packages()))) library(MASS)
  #  #if (!("MBESS" %in% (.packages()))) library(MBESS)
  #  #if (!("rootSolve" %in% (.packages()))) library(rootSolve)
  #  #if (!("doParallel" %in% (.packages()))) library(doParallel)  # this is for parallel processing loops (not for internal parallel processing of OpenMx)
  #  if (!("crayon" %in% (.packages()))) library(crayon)
  #  #if (!("psych" %in% (.packages()))) library(psych)
  #
    # if (activateRPB==TRUE) {
    #   if("RPushbullet" %in% rownames(installed.packages()) == FALSE) {install.packages("RPushbullet")}
    #   if (!("RPushbullet" %in% (.packages()))) library(RPushbullet)
    # }
  #}


  #######################################################################################################################
  ############## Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  ######################
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
  ################### Series of models with multiple drift effects invariant across primary studies #####################
  #######################################################################################################################

    print(paste0("#################################################################################"))
    print(paste0("######### Fit Homogeneity Model with Multiple Drift Effects Invariant ###########"))
    print(paste0("#################################################################################"))

    # Vector of drift coefficients to be tested for invariance
    targetNames <- multipleDrift; targetNames

    # get start values
    homDRIFTMultipleModModel <- NULL
      homDRIFTMultipleCTmodelobj <- ctmaInitFit$studyFitList[[1]]$ctmodelobj
      homDRIFTMultipleCTmodelobj$Tpoints <- maxTpoints
      homDRIFTMultipleFixedModel <- ctmaInitFit$studyFitList[[1]]$ctmodelobj
      for (j in 1:length(targetNames)) homDRIFTMultipleFixedModel$DRIFT <- sub(targetNames[j], "groupfixed", homDRIFTMultipleFixedModel$DRIFT)
      homDRIFTMultipleStartValues <- ctmaCompSV(singleStudyFits = ctmaInitFit$studyFitList,
                                                           ctModelObj = homDRIFTMultipleCTmodelobj,
                                                           fixedModel = homDRIFTMultipleFixedModel,
                                                           modModel = homDRIFTMultipleModModel, moderatorValues = NULL,
                                                           compSVMethod=compSVMethod)
      homDRIFTMultipleStartValues

      # FIT (NOT LOAD)
      if (useCTMultiGroupAlt == FALSE) {
        OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
        # change names of start values
        homDRIFTMultipleStartValues <- ctmaChangeSVLab(startValues = homDRIFTMultipleStartValues,
                                                                 ctmodelobj = homDRIFTMultipleCTmodelobj,
                                                                 fixedmodel = homDRIFTMultipleFixedModel, noOfStudies <- n.studies)
        results <- parallel::mclapply(seq(1, refits, by=1),
                            function(refits) ctsem::ctMultigroupFit(dat=datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                             omxStartValues = homDRIFTMultipleStartValues,
                                                             ctmodelobj = homDRIFTMultipleCTmodelobj,
                                                             fixedmodel = homDRIFTMultipleFixedModel),
                            mc.cores=coresToUse)
        OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        homDRIFTMultipleFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
        homDRIFTMultipleFit <- homDRIFTMultipleFit$mxobj
        #summary(homDRIFTMultipleFit)
      } else {
        method <- 2
        if (useCTMultiGroupAlt == "1") method <- 1
        OpenMx::mxOption(NULL, 'Number of Threads', numOfThreads)
        results <- parallel::mclapply(seq(1, refits, by=1),
                            function(refits) ctmaCTMultiGroupAlt(dat=datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                startValues=homDRIFTMultipleStartValues,
                                                                ctmodelobj = homDRIFTMultipleCTmodelobj,
                                                                fixedmodel=homDRIFTMultipleFixedModel,
                                                                #fastIterativeCoTiMA = fastIterativeCoTiMA,
                                                                method = method,
                                                                tryHard = TRUE),
                            mc.cores=coresToUse)
        OpenMx::mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
        homDRIFTMultipleFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
      }

      # SAVE
      if (length(saveMultipleDriftModelFit) > 0)  {
        tmp <- paste(targetNames, collapse=" "); tmp
        x1 <- paste0(saveMultipleDriftModelFit[1], "MX ", tmp, " ", (Sys.time()), ".rds"); x1
        x2 <- ""; x2
        ctmaSaveFile(activateRPB, activeDirectory, homDRIFTMultipleFit, x1, x2, silentOverwrite=silentOverwrite)
      }

      # Extract estimates & statistics
      homDRIFTMultiple_Drift_Coef <- homDRIFTMultipleFit$output$estimate; homDRIFTMultiple_Drift_Coef
      homDRIFTMultiple_Drift_SE <- homDRIFTMultipleFit$output$standardErrors; homDRIFTMultiple_Drift_SE
      homDRIFTMultiple_Minus2LogLikelihood <- homDRIFTMultipleFit$output$Minus2LogLikelihood; homDRIFTMultiple_Minus2LogLikelihood
      homDRIFTMultiple_estimatedParameters <- length(homDRIFTMultipleFit$output$estimate); homDRIFTMultiple_estimatedParameters
      homDRIFTMultiple_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
        homDRIFTMultiple_estimatedParameters; homDRIFTMultiple_df

      # Combine summary information
      homDRIFTMultipleDRIFT_effects <- cbind(homDRIFTMultiple_Drift_Coef, homDRIFTMultiple_Drift_SE); homDRIFTMultipleDRIFT_effects
      Tvalue <- (homDRIFTMultipleDRIFT_effects[,1]/homDRIFTMultipleDRIFT_effects[,2]); Tvalue
      homDRIFTMultipleDRIFT_effects <- cbind(homDRIFTMultipleDRIFT_effects, Tvalue)
      homDRIFTMultipleDRIFT_effects <- homDRIFTMultipleDRIFT_effects[which(rownames(homDRIFTMultipleDRIFT_effects) %in% targetNames), ]; homDRIFTMultipleDRIFT_effects

      # Label summary table
      paste0(homDRIFTMultipleDRIFT_effects, " (fixed effect)")
      rownames(homDRIFTMultipleDRIFT_effects) <- paste0(rownames(homDRIFTMultipleDRIFT_effects), " (fixed effect)"); homDRIFTMultipleDRIFT_effects
      colnames(homDRIFTMultipleDRIFT_effects) <- c("estimate", "(SE)", "Tvalue"); homDRIFTMultipleDRIFT_effects


      # Compute confidence intervals
      if (confidenceIntervals == TRUE) {
        print(paste0("#################################################################################"))
        print(paste0(" Computing Confidence Intervalls for Models with Multiple Drift Effects Invariant  "))
        print(paste0("#################################################################################"))


        # FIT (NOT LOAD)
        tmpModelMxobjFit <- homDRIFTMultipleFit
        ci <- OpenMx::mxCI(targetNames); ci
        tmpModelMxobj <- OpenMx::mxModel(tmpModelMxobjFit, ci)
        results <- OpenMx::mxRun(tmpModelMxobj, intervals=TRUE)
        homDRIFTMultipleFitCI <- results
        homDRIFTMultipleCI <- homDRIFTMultipleFitCI$output$confidenceIntervals


        # SAVE
        if (length(saveMultipleDriftModelFit) > 0)  {
          tmp <- paste0("CI_", targetNames, collapse=" "); tmp
          x1 <- paste0(saveMultipleDriftModelFit[1], "MX ", tmp, " ", (Sys.time()), ".rds"); x1
          x2 <- ""
          ctmaSaveFile(activateRPB, activeDirectory, homDRIFTMultipleFitCI, x1, x2, silentOverwrite=silentOverwrite)
        }

        multipleDRIFT_effects <- cbind(homDRIFTMultipleDRIFT_effects, homDRIFTMultipleCI[,c(1,3)]); multipleDRIFT_effects

      } # END confidence intervals

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
    #cat("Displays estimates from single study ctsem models and waits for user inoput to continue (checkMultipleStudyResults): ", checkMultipleStudyResults,  sep="")
    #cat(" ", "", sep="\n")
    cat("The active directory (activeDirectory) is: ", activeDirectory,  sep="\n")
    cat(" ", "", sep="\n")
    cat("The result file name that is created (resultsFileName): ", resultsFileName,  sep="")
    cat(" ", "", sep="\n")
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("#### - Multiple Drift Homogeneity Model: " , length(targetNames), " Coefficients Simultaneously Invariant - #####", "", sep=""))
    cat(sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients, their Standard Errors (SE), T-Values, 2-LL values and df:")
    cat(" ", "", sep="\n")
    print(round(homDRIFTMultiple_Drift_Coef, digits))
    cat(" ", "", sep="\n")
    cat("Note: If confidence intervals are reported, they do not necessarily have to be symmetric!")
  }

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n", st, "\n", et, "\nAnalysis successfully completed. \nThank you for using CoTiMA.\nHave a nice day!"))}

  sink()

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}



  if (!(exists("homDRIFTallFitCI"))) homDRIFTallFitCI <- c()

  model_Drift_Coef <- model_Diffusion_Coef <- model_T0var_Coef <- model_Cint_Coef <- list()

  if (useCTMultiGroupAlt == TRUE) {
    tmp1 <- grep("toV", names(homDRIFTMultiple_Drift_Coef)); tmp1
    all_model_Drift_Coef <-homDRIFTMultiple_Drift_Coef[tmp1]; all_model_Drift_Coef
    tmp1 <- grep("_G", names(all_model_Drift_Coef)); tmp1
    group_model_Drift_Coef <- all_model_Drift_Coef[tmp1]; group_model_Drift_Coef
    fixed_model_Drift_Coef <- all_model_Drift_Coef[!(all_model_Drift_Coef %in% group_model_Drift_Coef)]; fixed_model_Drift_Coef
    for (i in 1:n.studies) {
     tmp1 <- grep(paste0("_G", i), names(group_model_Drift_Coef)); tmp1
     tmp2 <- c(group_model_Drift_Coef[tmp1], fixed_model_Drift_Coef)
     model_Drift_Coef[[i]] <- matrix(tmp2[sort(names(tmp2))], n.latent, n.latent); model_Drift_Coef[[i]]

     tmp1 <- grep(paste0("_G", i), names(homDRIFTMultiple_Drift_Coef)); tmp1
     tmp2 <- homDRIFTMultiple_Drift_Coef[tmp1]; tmp2
     tmp3 <- grep("diffusion", names(tmp2)); tmp3
     model_Diffusion_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_Diffusion_Coef[[i]]
     #model_Diffusion_Coef[[i]] <- base2Full(model_Diffusion_Coef[[i]]); model_Diffusion_Coef[[i]]

     tmp3 <- grep("T0var", names(tmp2)); tmp3
     model_T0var_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_T0var_Coef[[i]]
     #model_T0var_Coef[[i]] <- base2Full(model_T0var_Coef[[i]]); model_T0var_Coef[[i]]
    }
  } else {
    tmp1 <- grep("toV", names(homDRIFTMultiple_Drift_Coef)); tmp1
    all_model_Drift_Coef <-homDRIFTMultiple_Drift_Coef[tmp1]; all_model_Drift_Coef
    tmp1 <- grep("Study_No", names(all_model_Drift_Coef)); tmp1
    group_model_Drift_Coef <- all_model_Drift_Coef[tmp1]; group_model_Drift_Coef
    fixed_model_Drift_Coef <- all_model_Drift_Coef[!(all_model_Drift_Coef %in% group_model_Drift_Coef)]; fixed_model_Drift_Coef
    for (i in 1:n.studies) {
      tmp1 <- grep(paste0("Study_No_", i), names(group_model_Drift_Coef)); tmp1
      tmp2 <- c(group_model_Drift_Coef[tmp1], fixed_model_Drift_Coef); tmp2
      names(tmp2) <- gsub(paste0("Study_No_", i, "_"), "", names(tmp2)); tmp2
      model_Drift_Coef[[i]] <- matrix(tmp2[sort(names(tmp2))], n.latent, n.latent); model_Drift_Coef[[i]]

      tmp1 <- grep(paste0("Study_No_", i), names(homDRIFTMultiple_Drift_Coef)); tmp1
      tmp2 <- homDRIFTMultiple_Drift_Coef[tmp1]; tmp2
      tmp3 <- grep("diffusion", names(tmp2)); tmp3
      model_Diffusion_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_Diffusion_Coef[[i]]
      #model_Diffusion_Coef[[i]] <- base2Full(model_Diffusion_Coef[[i]]); model_Diffusion_Coef[[i]]

      tmp3 <- grep("T0var", names(tmp2)); tmp3
      model_T0var_Coef[[i]] <- OpenMx::vech2full(tmp2[tmp3]); model_T0var_Coef[[i]]
      #model_T0var_Coef[[i]] <- base2Full(model_T0var_Coef[[i]]); model_T0var_Coef[[i]]
    }
  }



  results <- list(activeDirectory=activeDirectory, #sourceDirectory=sourceDirectory,
                  plot.type=NULL, model.type="mx",
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=homDRIFTMultipleFit, # , homDRIFTMultipleFitCI),
                  data=NULL, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model=paste(targetNames, "fixed", collapse=" & "),
                               estimates=round(homDRIFTMultiple_Drift_Coef, digits), #[]
                               minus2ll= round(homDRIFTMultiple_Minus2LogLikelihood, digits),
                               n.parameters = round(homDRIFTMultiple_estimatedParameters, digits),
                               df= c(round(homDRIFTMultiple_df, digits))))
  class(results) <- "CoTiMAFit"

  tmp <- paste(targetNames, collapse="_"); tmp
  saveRDS(results, paste0(activeDirectory, "ctmaMultipleMx ", tmp, " allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
