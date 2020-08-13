#######################################################################################################################
########################################### CoTiMA FullDrift stanct ###################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   datalong_all=datalong_all
#   #groups=groups
#   #groupsNamed=groupsNamed
#   activeDirectory=activeDirectory
#   sourceDirectory=sourceDirectory
#   resultsFilePrefix="ctmaFullStanct"
#   saveFilePrefix="CoTiMADFullDriftStanct"
#   # Primary Study Fits
#   ctmaInitFit=ctmaInitFit1
#   activateRPB=activateRPB
#   #checkSingleStudyResults=TRUE
#   digits=4
#   #type="mx",                           # this option switsches between mx/ctsem & stanct
#   retryattempts=retryattempts
#   refits=refits
#   NPSOL=NPSOL
#   coresToUse=coresToUse
#   fullDriftStartValues=fullDriftStartValues
#   #compSVMethod=compSVMethod
#   #useCTMultiGroupAlt=useCTMultiGroupAlt
#   #confidenceIntervals=confidenceIntervals
#   saveFullDriftModelFit=saveFullDriftModelFit
#   CoTiMAStanctArgs=list(test=TRUE, scaleTI=TRUE, scaleTime=1/1,
#                         savesubjectmatrices=FALSE, verbose=1,
#                         datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
#                         iter=1000, intoverstates=TRUE,
#                         binomial=FALSE, fit=TRUE,
#                         intoverpop=FALSE, stationary=FALSE,
#                         plot=FALSE, derrind="all",
#                         optimize=TRUE, optimcontrol=list(is=F, stochastic=FALSE),
#                         nlcontrol=list(),
#                         nopriors=TRUE,
#                         chains=2,
#                         #cores=numOfThreads,
#                         inits=NULL, forcerecompile=FALSE,
#                         savescores=FALSE, gendata=FALSE,
#                         control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
#                         verbose=0)
# }
# debug <- 0

#' ctmaFullStanct
#'
#' @param datalong_all ?
#' @param activeDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param ctmaInitFit ?
#' @param activateRPB ?
#' @param silentOverwrite ?
#' @param digits ?
#' @param n.latent ?
#' @param n.studies ?
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param numOfThreads ?
#' @param fullDriftStartValues ?
#' @param saveFullDriftModelFit ?
#' @param CoTiMAStanctArgs ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaFullStanct <- function(
  #type="fullDrift",
  datalong_all=NULL,
  #groups=NULL,
  #groupsNamed=NULL,
  # Directory names and file names
  activeDirectory = NULL,
  #sourceDirectory = NULL,
  resultsFilePrefix = "ctmaFullStanct",
  saveFilePrefix = "ctmaFullStanct",
  ctmaInitFit = NULL,
  activateRPB = NULL,
  silentOverwrite = FALSE,
  #checkSingleStudyResults = TRUE,
  digits = 4,
  #type = "mx",
  n.latent = NULL,
  n.studies = NULL,
  retryattempts = 1,
  refits = 1,
  NPSOL = FALSE,
  coresToUse = 1,
  numOfThreads = 1,
  fullDriftStartValues = c(),
  #compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  #useCTMultiGroupAlt = TRUE,
  #confidenceIntervals = FALSE,
  saveFullDriftModelFit = c(),
  CoTiMAStanctArgs=list(test=TRUE, scaleTI=TRUE, scaleTime=1/1,
                        savesubjectmatrices=FALSE, verbose=1,
                        datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
                        iter=1000, intoverstates=TRUE,
                        binomial=FALSE, fit=TRUE,
                        intoverpop=FALSE, stationary=FALSE,
                        plot=FALSE, derrind="all",
                        optimize=TRUE, optimcontrol=list(is=F, stochastic=FALSE),
                        nlcontrol=list(),
                        nopriors=TRUE,
                        chains=2,
                        #cores=numOfThreads,
                        inits=NULL, forcerecompile=FALSE,
                        savescores=FALSE, gendata=FALSE,
                        control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
                        verbose=0)
)
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
    #}
  #}


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
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctmodelobj$manifestNames; manifestNames
    #driftNames <- c(ctmaInitFit$studyFitList[[1]]$ctmodelobj$DRIFT); driftNames
    driftNames <- ctmaInitFit$parameterNames$DRIFT; driftNames
  }

  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################

  { # fit stanct model


    # Make model with most time points
    stanctModel <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                           DRIFT=matrix(driftNames, nrow=n.latent, ncol=n.latent),
                           LAMBDA=diag(n.latent),
                           CINT=matrix(0, nrow=n.latent, ncol=1),
                           T0MEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                           MANIFESTMEANS = matrix(c(0), nrow = n.latent, ncol = 1),
                           MANIFESTVAR=matrix(0, nrow=n.latent, ncol=n.latent),
                           MANIFESTTRAITVAR = 'auto',
                           type = 'stanct',
                           n.TIpred = (n.studies-1),
                           TIpredNames = paste0("TI", 1:(n.studies-1)),
                           TIPREDEFFECT = matrix(0, n.latent, (n.studies-1)))

    # only relevant for moderator analysis
    #ctsemStanctModModel <- ctsemStanctModel
    #ctsemStanctModModel$n.TIpred <- (n.studies-1+length(moderatorNumber)); ctsemStanctModModelTemplate$n.TIpred
    #ctsemStanctModModel$TIpredNames <- paste0("TI", 1:(n.studies-1+length(moderatorNumber))); ctsemStanctModModelTemplate$TIpredNames
    #ctsemStanctModModel$TIPREDEFFECT = matrix(0, n.latent, (n.studies-1+length(moderatorNumber))); ctsemStanctModModelTemplate$TIPREDEFFECT

    stanctModel$pars[stanctModel$pars$matrix %in% 'DRIFT',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'T0MEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'LAMBDA',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'CINT',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTMEANS',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    stanctModel$pars[stanctModel$pars$matrix %in% 'MANIFESTVAR',paste0(stanctModel$TIpredNames,'_effect')] <- FALSE
    # set diag elements of T0var to FALSE
    #tmp1 <- which(stanctModel$pars$matrix == 'T0VAR'); tmp1
    #tmp2 <- which(stanctModel$pars$row == stanctModel$pars$col); tmp2
    #tmp3 <- tmp1[which(tmp1 %in% tmp2)]; tmp3
    #stanctModel$pars[tmp3, paste0(stanctModel$TIpredNames,'_effect')] <- FALSE

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
      control=CoTiMAStanctArgs$control,
      verbose=CoTiMAStanctArgs$verbose,
      cores=coresToUse)

    fullDriftStanctFit <- summary(fitStanctModel, digits=digits)

    # SAVE
    if (length(saveFullDriftModelFit) > 0) {
      x1 <- paste0(saveFullDriftModelFit[1], " ", (Sys.time()), ".rds"); x1
      x2 <- c()
      ctmaSaveFile(activateRPB, activeDirectory, fitStanctModel, x1, x2, silentOverwrite=silentOverwrite)
    }

    # Extract estimates & statistics
    Tvalues <- fullDriftStanctFit$popmeans[,1]/fullDriftStanctFit$popmeans[,2]; Tvalues
    fullDrift_Coeff <- round(cbind(fullDriftStanctFit$popmeans, Tvalues), digits); fullDrift_Coeff
    fullDrift_Minus2LogLikelihood  <- -2*fullDriftStanctFit$logprob; fullDrift_Minus2LogLikelihood
    fullDrift_estimatedParameters  <- fullDriftStanctFit$npars; fullDrift_estimatedParameters
    fullDrift_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
      fullDrift_estimatedParameters; fullDrift_df
    model_Drift_Coef <- fullDrift_Coeff[grep("toV", rownames(fullDrift_Coeff)), 1]; model_Drift_Coef
    names(model_Drift_Coef) <- rownames(fullDrift_Coeff)[grep("toV", rownames(fullDrift_Coeff))]; model_Drift_Coef
    model_Diffusion_Coef <- fullDrift_Coeff[grep("diffusion", rownames(fullDrift_Coeff)), 1]; model_Diffusion_Coef
    names(model_Diffusion_Coef) <- rownames(fullDrift_Coeff)[grep("diffusion", rownames(fullDrift_Coeff))]; model_Diffusion_Coef
    model_T0var_Coef <- fullDrift_Coeff[grep("T0var", rownames(fullDrift_Coeff)), 1]; model_T0var_Coef
    names(model_T0var_Coef) <- rownames(fullDrift_Coeff)[grep("T0var", rownames(fullDrift_Coeff))]; model_Diffusion_Coef

    ## LOAD
    #if ((testDRIFTallModel == TRUE) & (length(loadFullDriftModelFit) > 0) ) {
    #  x1 <- paste0(loadFullDriftModelFit[1], " fullDriftFit.rds"); x1
    #  fullDriftFit <- readRDS(file=x1)
    #  FitList$FullDriftModelFit <- fullDriftFit
    #}


    ### Numerically compute Optimal Time lag sensu Dormann & Griffin (2015)
    #driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent); driftMatrix
    driftMatrix <- matrix(model_Drift_Coef, n.latent, n.latent, byrow=T); driftMatrix # byrow set because order is different compared to mx model
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
  } ## END  fit stanct model

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
    cat(" ", "", sep="\n")
    cat("Fitted model is of type: stanct !!!!", sep="")
    cat(" ", "", sep="\n")
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
    #cat("Displays estimates from single study ctsem models and waits for user inoput to continue (checkSingleStudyResults): ", checkSingleStudyResults,  sep="")
    #cat(" ", "", sep="\n")
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
    print(round(fullDrift_Coeff, digits))
    cat(" ", "", sep="\n")
    cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
    cat("-2 Log Likelihood: ", round(fullDrift_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Overall Number of Estimated Parameters: ", round(fullDrift_estimatedParameters, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(fullDrift_df, digits), sep="")
    cat(" ", "", sep="\n")
    #cat("Comparison with Unconstrained Modell (All Samples Analyzed Separately (~Heterogeneity Model)):")
    #cat(" ", "", sep="\n")
    #cat("Delta(df): ", round(allStudies_fullDrift_df, digits), sep="")
    #cat(" ", "", sep="\n")
    #cat("Delta(-2LL) (= Delta(Chi-square)): ", round(allStudies_fullDrift_Minus2LogLikelihood, digits), sep="")
    #cat(" ", "", sep="\n")
    #cat("p-value (with double number of digits): ", round(allStudies_fullDrift_prob, 2*digits), sep="")
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


  tmp1 <- grep("CINT", rownames(fullDriftStanctFit$parmatrices)); tmp1
  tmp2 <- grep("asym", rownames(fullDriftStanctFit$parmatrices)); tmp2
  tmp3 <- grep("dt", rownames(fullDriftStanctFit$parmatrices)); tmp3
  tmp4 <- tmp1[(tmp1 %in% c(tmp2, tmp3)) == FALSE]; tmp4
  model_Cint_Coef <- fullDriftStanctFit$parmatrices[tmp4, 3]; model_Cint_Coef

  results <- list(activeDirectory=activeDirectory, #sourceDirectory=sourceDirectory,
                  plot.type="drift",  model.type="stanct",
                  coresToUse=coresToUse, n.studies=1,
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=list(fitStanctModel), # , fullDriftFitCI),
                  data=datalong_all, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model="all drift fixed (hom. model)",
                               estimates=round(fullDrift_Coeff, digits),
                               minus2ll= round(fullDrift_Minus2LogLikelihood, digits),
                               n.parameters = round(fullDrift_estimatedParameters, digits),
                               df= c(round(fullDrift_df, digits)),
                               opt.lag = optimalCrossLag,
                               max.effects = round(maxCrossEffect, digits)))
  class(results) <- "CoTiMAFit"

  saveRDS(results, paste0(activeDirectory, "ctmaFullStanct_allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
