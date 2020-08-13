#######################################################################################################################
########################################### CoTiMA SingleDrift STANCT #################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   single=c("V1toV2", "V2toV2")
#   datalong_all=datalong_all
#   #groups=groups
#   #groupsNamed=groupsNamed
#   activeDirectory=activeDirectory
#   sourceDirectory=sourceDirectory
#   resultsFilePrefix="ctmaSingleStanct"
#   saveFilePrefix="ctmaSingleStanct"
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
#   saveSingleDriftModelFit=saveSingleDriftModelFit
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

#' ctmaSingleStanct
#'
#' @param single ?
#' @param datalong_all ?
#' @param activeDirectory ?
#' @param sourceDirectory ?
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
#' @param saveSingleDriftModelFit ?
#' @param CoTiMAStanctArgs ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaSingleStanct <- function(
  single=c(),
  datalong_all=NULL,
  #groups=NULL,
  #groupsNamed=NULL,
  # Directory names and file names
  activeDirectory = NULL,
  sourceDirectory = NULL,
  resultsFilePrefix = "ctmaSingleStanct",
  saveFilePrefix = "ctmaSingleStanct",
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
  saveSingleDriftModelFit = c(),
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
    #print(paste0("#################################################################################"))
    #print(paste0("############################ Attach Further Packages ############################"))
    #print(paste0("#################################################################################"))

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

   # if (activateRPB==TRUE) {
  #    if("RPushbullet" %in% rownames(installed.packages()) == FALSE) {install.packages("RPushbullet")}
  #    if (!("RPushbullet" %in% (.packages()))) library(RPushbullet)
  #  }
  #}


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
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctmodelobj$manifestNames; manifestNames
    driftNames <- c(ctmaInitFit$studyFitList[[1]]$ctmodelobj$DRIFT)
  }

  #######################################################################################################################
  ######################### Series of homogeneity models with single drift effects invariant ############################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("##### Fit Series of Homogeneity Models with Single Drift Effects Invariant ######"))
    print(paste0("#################################################################################"))

    stanctModel <- fitStanctModel <- singleStanctFit <- list()
    targetNames <- single

    for (i in 1:length(single)) {
      #i <- 1
      # Make model with most time points
      stanctModel[[i]] <- ctsem::ctModel(n.latent=n.latent, n.manifest=n.latent, Tpoints=maxTpoints, manifestNames=manifestNames,    # 2 waves in the template only
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

      stanctModel[[i]]$pars[stanctModel[[i]]$pars$matrix %in% 'DRIFT',paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- TRUE
      stanctModel[[i]]$pars[stanctModel[[i]]$pars$matrix %in% 'T0MEANS',paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- FALSE
      stanctModel[[i]]$pars[stanctModel[[i]]$pars$matrix %in% 'LAMBDA',paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- FALSE
      stanctModel[[i]]$pars[stanctModel[[i]]$pars$matrix %in% 'CINT',paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- FALSE
      stanctModel[[i]]$pars[stanctModel[[i]]$pars$matrix %in% 'MANIFESTMEANS',paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- FALSE
      stanctModel[[i]]$pars[stanctModel[[i]]$pars$matrix %in% 'MANIFESTVAR',paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- FALSE

      # the target effect:
      targetNames[i]
      tmp1 <- which(stanctModel[[i]]$pars$matrix == "DRIFT"); tmp1
      tmp2 <- which(stanctModel[[i]]$pars[tmp1, "param"] == targetNames[i]); tmp2

      #stanctModel[[i]]$pars[tmp, paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- TRUE
      #stanctModel[[i]]$pars[tmp1[tmp2],]
      stanctModel[[i]]$pars[tmp1[tmp2], paste0(stanctModel[[i]]$TIpredNames,'_effect')] <- FALSE
      stanctModel[[i]]$pars

      # set diag elements of T0var to FALSE
      #tmp1 <- which(stanctModel$pars$matrix == 'T0VAR'); tmp1
      #tmp2 <- which(stanctModel$pars$row == stanctModel$pars$col); tmp2
      #tmp3 <- tmp1[which(tmp1 %in% tmp2)]; tmp3
      #stanctModel$pars[tmp3, paste0(stanctModel$TIpredNames,'_effect')] <- FALSE

      fitStanctModel[[i]] <- ctsem::ctStanFit(
        datalong = datalong_all,
        ctstanmodel = stanctModel[[i]],
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

      singleStanctFit[[i]] <- summary(fitStanctModel[[i]], digits=digits)

      # SAVE
      if (length(saveSingleDriftModelFit) > 0) {
        #x1 <- paste0(saveSingleDriftModelFit[1], "STANCT", " ", (Sys.time()), ".rds"); x1
        x1 <- paste0(saveSingleDriftModelFit[1], "STANCT", " ", targetNames[i], " ", (Sys.time()), ".rds"); x1
        x2 <- c()
        ctmaSaveFile(activateRPB, activeDirectory, fitStanctModel[[i]], x1, x2, silentOverwrite=silentOverwrite)
      }
    }

    # Extract estimates & statistics
    Tvalues <- single_Coeff <- single_Minus2LogLikelihood <- single_estimatedParameters <- single_df <- list()
    model_Drift_Coef <- model_Diffusion_Coef <- model_T0var_Coef <- model_Cint_Coef <- allResults <- list()

    for (i in 1:length(targetNames)){
      Tvalues[[i]] <- singleStanctFit[[i]]$popmeans[,1]/singleStanctFit[[i]]$popmeans[,2]; Tvalues[[i]]
      single_Coeff[[i]] <- round(cbind(singleStanctFit[[i]]$popmeans, Tvalues[[i]]), digits); single_Coeff[[i]]
      # re-label
      tmp <- colnames(single_Coeff[[i]]); tmp
      tmp[length(tmp)] <- "Tvalues"
      colnames(single_Coeff[[i]]) <- tmp; single_Coeff[[i]]
      tmp1 <- which(rownames(single_Coeff[[i]]) == targetNames[i]); tmp1
      tmp2 <- paste0(rownames(single_Coeff[[i]])[tmp1] , " (fixed)"); tmp2
      rownames(single_Coeff[[i]])[tmp1] <- tmp2; single_Coeff[[i]]
      #
      single_Minus2LogLikelihood[[i]]  <- -2*singleStanctFit[[i]]$logprob; single_Minus2LogLikelihood[[i]]
      single_estimatedParameters[[i]]  <- singleStanctFit[[i]]$npars; single_estimatedParameters[[i]]
      #single_estimatedParameters[[i]]  <- length(fitStanctModel[[i]]$stanfit$optimfit$par); single_estimatedParameters[[i]]
      single_df[[i]] <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
        single_estimatedParameters[[i]]; single_df[[i]]
      model_Drift_Coef[[i]] <- single_Coeff[[i]][grep("toV", rownames(single_Coeff[[i]])), 1]; model_Drift_Coef[[i]]
      #names(model_Drift_Coef[[i]]) <- driftNames; model_Drift_Coef[[i]]
      names(model_Drift_Coef[[i]]) <- rownames(single_Coeff[[i]])[grep("toV", rownames(single_Coeff[[i]]))]; model_Drift_Coef
      model_Diffusion_Coef[[i]] <- single_Coeff[[i]][grep("diffusion", rownames(single_Coeff[[i]])), 1]; model_Diffusion_Coef[[i]]
      names(model_Diffusion_Coef[[i]]) <- rownames(single_Coeff[[i]])[grep("diffusion", rownames(single_Coeff[[i]]))]; model_Diffusion_Coef[[i]]
      model_T0var_Coef[[i]] <- single_Coeff[[i]][grep("T0var", rownames(single_Coeff[[i]])), 1]; model_T0var_Coef[[i]]
      names(model_T0var_Coef[[i]]) <- rownames(single_Coeff[[i]])[grep("T0var", rownames(single_Coeff[[i]]))]; model_T0var_Coef[[i]]

      allResults[[i]] <- list(estimates=single_Coeff[[i]], Minus2LL=single_Minus2LogLikelihood[[i]],
                              n.parameters=single_estimatedParameters[[i]], df=single_df[[i]])

      model_Cint_Coef[[i]] <- NULL
    }

    ## LOAD
    #if ((testDRIFTallModel == TRUE) & (length(loadFullDriftModelFit) > 0) ) {
    #  x1 <- paste0(loadFullDriftModelFit[1], " homDRIFTallFit.rds"); x1
    #  homDRIFTallFit <- readRDS(file=x1)
    #  FitList$FullDriftModelFit <- homDRIFTallFit
    #}
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
    cat("Fitted model is of type: STANCT", sep="")
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
    cat(paste("- Multi-Sample Homogeneity Model (only a single Drift Coefficient Invariant per Model) -", "", sep=""))
    cat(sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients, their Standard Errors (SE), T-Values, 2-LL values and df:")
    cat(" ", "", sep="\n")
    print(allResults)
    cat(" ", "", sep="\n")
    cat("Note: Confidence intervals do not necessarily have to be symmetric!")
  }

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n", st, "\n", et, "\nAnalysis successfully completed. \nThank you for using CoTiMA.\nHave a nice day!"))}

  sink()

  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","CoTiMA has finished!"))}

  results <- list(activeDirectory=activeDirectory, sourceDirectory=sourceDirectory,
                  plot.type="drift", model.type="stanct",
                  coresToUse=coresToUse, n.studies=length(single),
                  n.latent=n.latent,
                  studyList=ctmaInitFit$studyList, studyFitList=fitStanctModel, # , homDRIFTallFitCI),
                  data=datalong_all, statisticsList=ctmaInitFit$statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef),
                  parameterNames=ctmaInitFit$parameterNames,
                  summary=list(model=c(paste0(targetNames, " fixed")),
                    estimates=lapply(allResults, function(extract) extract$estimates),
                               minus2ll= lapply(allResults, function(extract) extract$Minus2LL),
                               n.parameters = lapply(allResults, function(extract) extract$n.parameters),
                               df= lapply(allResults, function(extract) extract$df)))
  class(results) <- "CoTiMAFit"

  saveRDS(results, paste0(activeDirectory, "ctmaSingleStanct_allResults", " ", Sys.time(), ".rds"))

  invisible(results)

} ### END function definition
