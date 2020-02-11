debug <- 0
if (debug == 1) {
  ######################### ENTER MODEL-SPECIFIC SETTINGS AT THE END #########################
  # Directory names and file names
  workingDirectory= NULL
  sourceDirectory= NULL
  resultsFileName="CoTiMA.txt"
  filePrefix="CoTiMA"
  # Primary Study Information
  primaryStudies=NULL
  # Workflow
  activateRPB=FALSE
  checkSingleStudyResults=TRUE
  # General Model Setup
  nlatents=NULL
  manifesttraitvar=FALSE
  # Load/save Raw Data
  saveRawData=list()
  # Fitting Parameters
  retryattempts=30
  refits=5
  extraRefits=c()
  factorExtraRefits=1
  NPSOL=FALSE
  coresToUse=c(1)
  useCTMultigroupFitAlt=FALSE
  # Figure Parameters
  plotCrossEffects=TRUE
  plotAutoEffects=TRUE
  timeUnit="timeUnit"
  timeRange=c()
  yLimitsForEffects=c()
  digits=4
  # Specific Model Tests
  testHeterogeneityModel=FALSE
  testDRIFTallModel=TRUE
  testDRIFTallWOdtModel=FALSE
  testDRIFTSingleModel=FALSE
  testDRIFTCrossModel=FALSE
  testDRIFTAutoModel=FALSE
  testAllInvariantModel=FALSE
  testModeratorModel=FALSE
  moderatorNumber=c()
  userSpecifiedModel=list()
  # Further Analyses
  confidenceIntervals=FALSE
  statisticalPower=c()
  statisticalPowerRefits=1
  yLimitsForPower=c()
  fixedAndRandomEffects=TRUE
  publicationBias=FALSE
  # Save computed model fits or load previous fits (then provide (old) file prefix)
  saveSingleStudyModelFit=c()
  loadSingleStudyModelFit=c()
  saveHeterogeneityModelFit=c()
  loadHeterogeneityModelFit=c()
  saveDRIFTAllModelFit=c()
  loadDRIFTAllModelFit=c()
  saveDRIFTSingleModelFit=c()
  loadDRIFTSingleModelFit=c()
  saveDRIFTCrossModelFit=c()
  loadDRIFTCrossModelFit=c()
  saveDRIFTAutoModelFit=c()
  loadDRIFTAutoModelFit=c()
  saveDRIFTallWOdtModelFit=c()
  loadDRIFTallWOdtModelFit=c()
  saveAllInvariantModelFit=c()
  loadAllInvariantModelFit=c()
  saveHomAllWOSingleDriftModelFit=c()
  loadHomAllWOSingleDriftModelFit=c()
  saveModeratorModelFit=c()
  loadModeratorModelFit=c()
  loadRawData=list()                     #DEPRECATED
  saveSingleStudyFit=c()                 #DEPRECATED
  loadSingleStudyFit=c()                 #DEPRECATED
  saveHomogeneityStudyFit=FALSE          #DEPRECATED
  loadHomogeneityStudyFit=FALSE          #DEPRECATED
  saveHeterogeneityStudyFit=FALSE        #DEPRECATED
  loadHeterogeneityStudyFit=FALSE        #DEPRECATED
  loadFilePrefix=NULL                    #DEPRECATED
  unixLike=FALSE                         #DEPRECATED
  figureFilePrefix=c()

  ######################### ENTER MODEL-SPECIFIC SETTINGS HERE #################################
  #testUserSpecifiedModel=list(singleStudyModelFits = 1:5,
  #                        moderatorValues = moderatorValues1,
  #                        listOfFixedModelMatrices = listOfFixedModelMatrices1,
  #                        listOfModeratorModelMatrices = listOfModeratorModelMatrices1)

  ##### no WD in Package Build!?
  # workingDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/"
  # sourceDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/"


   resultsFileName="CoTiMA 4.txt"
  filePrefix="CoTiMA 4"
  # Primary Study information
  primaryStudies=firstStudyList
  nlatents=2
  # Workflow
  checkSingleStudyResults=FALSE
  # Figure Parameters
  timeUnit="Months"
  # Specific Model Tests
  testModeratorModel=FALSE
  moderatorNumber = 1
  # Fitting Parameters
  refits=1
  confidenceIntervals=TRUE
  testDRIFTallModel=TRUE
  #loadSingleStudyModelFit=c("CoTiMA 2", 1:5)
  #saveUserSpecifiedModel=c("CoTiMA 2")
  loadUserSpecifiedModel=c()
}
debug <- 0

### STOP ###


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' Title
#'
#' @param workingDirectory workingDirectory
#' @param sourceDirectory sourceDirectory
#' @param resultsFileName the prefix for the result file that is created
#' @param filePrefix the prefix for the figure files and the model fit files that are created
#' @param primaryStudies list of lists: list(deltas, sampleSizes, empcovs, moderators, startValues, studyNumbers)
#' @param activateRPB set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param checkSingleStudyResults displays estimates from single study ctsem models and waits for user inoput to continue
#' @param nlatents n latent variables
#' @param manifesttraitvar test model with unit effects/unobserved heterogeneity (requires more then 2 time points could still fail then)
#' @param saveRawData TBN
#' @param retryattempts number of iterations to be used by ctsem
#' @param refits number of re-fits used by this CoTiMA Function
#' @param extraRefits (vector of) study number(s) that should be re-fitted even more (difficult models)
#' @param factorExtraRefits factor by which "refits" are multiplied for difficult models
#' @param NPSOL request NPSOL optimizer (sometimes better estimates, e.g., for confidence intervals)
#' @param coresToUse could be set to -n (usually -1) which is subtracted from the overall number available (on non Windows machines)
#' @param useCTMultigroupFitAlt TBN
#' @param plotCrossEffects plot cross effects across a range of discrete intervals
#' @param plotAutoEffects plot auto effects across a range of discrete intervals
#' @param timeUnit timelag unit to lable x-axis of plots
#' @param timeRange used for Plotting and Poweranalysis. Populated by 0 to 1.5*maxDelta (Steptwidth = 1) if not specified as list(min,max,stepWidth)
#' @param yLimitsForEffects used the y-axis of Drift-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
#' @param digits rounding used in output
#' @param testHeterogeneityModel multigoup model with all parameters freely estimated across groups
#' @param testDRIFTallModel multigoup model (CoTiMA!) with all DRIFT parameters invariantly estimated across groups
#' @param testDRIFTallWOdtModel TBN
#' @param testDRIFTSingleModel multigoup models with every single DRIFT parameters invariantly estimated across groups
#' @param testDRIFTCrossModel multigoup model to test equality of cross effects
#' @param testDRIFTAutoModel multigoup model to test equality of auto effects
#' @param testAllInvariantModel multigoup model (CoTiMA!) with ALL DRIFT paramters invariant. Required for Statistical Power Computations
#' @param testModeratorModel performs moderator analysis. Fails if no (vector of) moderator values are provided (e.g. moderator1 <- c(13))
#' @param moderatorNumber selects from the vector of moderator values (the 1st is used if not otehrwiseotherwise specificed)
#' @param testUserSpecifiedModel analysis of a user specified model
#' @param confidenceIntervals estimate confidence intervals (for all requested models)
#' @param statisticalPower beta, e.g., 95. Computes required sample sizes for discrete time drift paramters across range of time lags. Sets testAllInvariantModel=TRUE
#' @param statisticalPowerRefits TBN
#' @param yLimitsForPower used the y-axis of RequiredSampleSize-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
#' @param fixedAndRandomEffects performs fixed and random effects analysis
#' @param publicationBias perform Egger's test, funnel plots and PET-PEESE estimates
#' @param saveSingleStudyModelFit save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#' @param loadSingleStudyModelFit load the fit of single study ctsem models
#' @param saveHeterogeneityModelFit save the fit of the heterogeneity model
#' @param loadHeterogeneityModelFit load the fit of the heterogeneity model
#' @param saveDRIFTAllModelFit save the fit of the full CoTiMA model
#' @param loadDRIFTAllModelFit load the fit of the full CoTiMA model
#' @param saveDRIFTSingleModelFit save the fit of the CoTiMA models with single invariant Drift Parameters
#' @param loadDRIFTSingleModelFit load the fit of the CoTiMA models with single invariant Drift Parameters
#' @param saveDRIFTCrossModelFit TBN
#' @param loadDRIFTCrossModelFit TBN
#' @param saveDRIFTAutoModelFit TBN
#' @param loadDRIFTAutoModelFit TBN
#' @param saveDRIFTallWOdtModelFit TBN
#' @param loadDRIFTallWOdtModelFit TBN
#' @param saveAllInvariantModelFit TBN
#' @param loadAllInvariantModelFit TBN
#' @param saveHomAllWOSingleDriftModelFit save models required for calculation of statistical power (only relevant is statisticalPower is requested)
#' @param loadHomAllWOSingleDriftModelFit load models required for calculation of statistical power (only relevant is statisticalPower is requested)
#' @param saveModeratorModelFit save moderator models (the moderator number will be added)
#' @param loadModeratorModelFit load moderator models
#' @param saveUserSpecifiedModel TBN
#' @param loadUserSpecifiedModel TBN
#' @param loadRawData DEPRECATED: information which raw data to load is now included in list of primary studies
#' @param saveSingleStudyFit DEPRECATED: save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#' @param loadSingleStudyFit DEPRECATED: load the fit of single study ctsem models
#' @param saveHomogeneityStudyFit DEPRECATED: save the fit of the CoTiMA model
#' @param loadHomogeneityStudyFit load the fit of the CoTiMA model
#' @param saveHeterogeneityStudyFit save the fit of the heterogeneity model
#' @param loadHeterogeneityStudyFit load the fit of the heterogeneity model
#' @param loadFilePrefix DEPRECATED: the prefix of (older) files that contain model fits to be be loaded (useful if figureFilePrefix has changed since then)
#' @param unixLike DEPRECATED: instead the machine type is detected and parallelisation of loops with mclapply on unix machine activated
#' @param figureFilePrefix TBN
#'
#' @return
#' @export
#'

CoTiMA <- function(
  # Directory names and file names
  workingDirectory= NULL,
  sourceDirectory= NULL,
  resultsFileName="CoTiMA.txt",           # the prefix for the result file that is created
  filePrefix="CoTiMA",                    # the prefix for the figure files and the model fit files that are created

  # Primary Study Information
  primaryStudies=NULL,                    #NEW: list of lists: list(deltas, sampleSizes, empcovs, moderators, startValues, studyNumbers)

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,                      #set to TRUE to receive push messages with CoTiMA notifications on your phone
  checkSingleStudyResults=TRUE,          # displays estimates from single study ctsem models and waits for user inoput to continue

  # General Model Setup
  nlatents=NULL,
  manifesttraitvar=FALSE,                 # test model with unit effects/unobserved heterogeneity (requires more then 2 time points could still fail then)

  # Load/save Raw Data
  saveRawData=list(),

  # Fitting Parameters
  retryattempts=30,                       # number of iterations to be used by ctsem
  refits=5,                               # number of re-fits used by this CoTiMA Function
  extraRefits=c(),                        # (vector of) study number(s) that should be re-fitted even more (difficult models)
  factorExtraRefits=1,                    # factor by which "refits" are multiplied for difficult models
  NPSOL=FALSE,                            # request NPSOL optimizer (sometimes better estimates, e.g., for confidence intervals)
  coresToUse=c(1),                        # could be set to -n (usually -1) which is subtracted from the overall number available (on non Windows machines)
  useCTMultigroupFitAlt=FALSE,

  # Figure Parameters
  plotCrossEffects=TRUE,                  # plot cross effects across a range of discrete intervals
  plotAutoEffects=TRUE,                   # plot auto effects across a range of discrete intervals
  timeUnit="timeUnit",                    # timelag unit to lable x-axis of plots
  timeRange=c(),                          # used for Plotting and Poweranalysis. Populated by 0 to 1.5*maxDelta (Steptwidth = 1) if not specified as list(min,max,stepWidth)
  yLimitsForEffects=c(),                  # used the y-axis of Drift-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
  digits=4,                               # rounding used in output

  # Specific Model Tests
  testHeterogeneityModel=FALSE,           # multigoup model with all parameters freely estimated across groups
  testDRIFTallModel=TRUE,                 # multigoup model (CoTiMA!) with all DRIFT parameters invariantly estimated across groups
  testDRIFTallWOdtModel=FALSE,
  testDRIFTSingleModel=FALSE,             # multigoup models with every single DRIFT parameters invariantly estimated across groups
  testDRIFTCrossModel=FALSE,              # multigoup model to test equality of cross effects
  testDRIFTAutoModel=FALSE,               # multigoup model to test equality of auto effects
  testAllInvariantModel=FALSE,            # multigoup model (CoTiMA!) with ALL DRIFT paramters invariant. Required for Statistical Power Computations
  testModeratorModel=FALSE,               # performs moderator analysis. Fails if no (vector of) moderator values are provided (e.g. moderator1 <- c(13)  )
  moderatorNumber=c(),                    # selects from the vector of moderator values (the 1st is used if not otehrwiseotherwise specificed)
  testUserSpecifiedModel=list(singleStudyModelFits = NULL,
                              moderatorValues = NULL,
                              listOfFixedModelMatrices = list(T0VAR=NULL, DIFFUSION=NULL, DRIFT=NULL),
                              listOfModeratorModelMatrices = list(T0VAR=NULL, DIFFUSION=NULL, DRIFT=NULL)),              # analysis of a user specified model

  # Further Analyses
  confidenceIntervals=FALSE,              # estimate confidence intervals (for all requested models)
  statisticalPower=c(),                   # beta, e.g., 95. Computes required sample sizes for discrete time drift paramters across range of time lags. Sets testAllInvariantModel=TRUE
  statisticalPowerRefits=1,
  yLimitsForPower=c(),                    # used the y-axis of RequiredSampleSize-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
  fixedAndRandomEffects=TRUE,             # performs fixed and random effects analysis
  publicationBias=FALSE,                  # perform Egger's test, funnel plots and PET-PEESE estimates

  ## Save computed model fits or load previous fits (then provide (old) file prefix)
  # ctsem models for all primasry studies (always required for getting good starting values (and model fit for heterogeneity model without fitting it))
  saveSingleStudyModelFit=c(),            # save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
  loadSingleStudyModelFit=c(),            # load the fit of single study ctsem models
  #startValues=list(),
  # all heterogneity model (only included for completeness; the same is achieved with single model fitting, which is always done anyway)
  saveHeterogeneityModelFit=c(),          # save the fit of the heterogeneity model
  loadHeterogeneityModelFit=c(),          # load the fit of the heterogeneity model
  # CoTiMA (all homogeneity model; model with all drift effects simultaneously estimated as invariant across primary studies)
  saveDRIFTAllModelFit=c(),               # save the fit of the full CoTiMA model.
  loadDRIFTAllModelFit=c(),               # load the fit of the full CoTiMA model.
  # Single drift homogeneity model; model with each drift effects (one after another) estimated as invariant across primary studies
  saveDRIFTSingleModelFit=c(),            # save the fit of the CoTiMA models with single invariant Drift Parameters
  loadDRIFTSingleModelFit=c(),            # load the fit of the CoTiMA models with single invariant Drift Parameters
  # Model testing for invariance AND equality of cross effects respective auto effects
  saveDRIFTCrossModelFit=c(),
  loadDRIFTCrossModelFit=c(),
  saveDRIFTAutoModelFit=c(),
  loadDRIFTAutoModelFit=c(),
  # Model in which time lag values (dT) are nullified (replaced by overall mean dT). Usefull to descriptively compare -2LL with DRIFTAllModelFit
  saveDRIFTallWOdtModelFit=c(),
  loadDRIFTallWOdtModelFit=c(),
  # Model that adds further invariance constraint to the CoTiMA model (diffusion co-/variances; T0 co-/variances). Also required for Power Calculations (then automatically requested)
  saveAllInvariantModelFit=c(),
  loadAllInvariantModelFit=c(),
  # Models required for Power Calculations (Note: Only computed for drift effects - usually very very poor results for auto effects (where power declines with time and Power Analysis not really needed))
  saveHomAllWOSingleDriftModelFit=c(),    # save models required for calculation of statistical power (only relevant is statisticalPower is requested)
  loadHomAllWOSingleDriftModelFit=c(),
  # Moderator Models
  saveModeratorModelFit=c(),              # save moderator models (the moderator number will be added)
  loadModeratorModelFit=c(),              # load moderator models
  # User-Specified Models
  saveUserSpecifiedModel=c(),
  loadUserSpecifiedModel=c(),

  loadRawData=list(),                     #DEPRECATED: information which raw data to load is now included in list of primary studies
  saveSingleStudyFit=c(),                 #DEPRECATED: save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
  loadSingleStudyFit=c(),                 #DEPRECATED: load the fit of single study ctsem models
  saveHomogeneityStudyFit=FALSE,          #DEPRECATED: save the fit of the CoTiMA model.
  loadHomogeneityStudyFit=FALSE,          #DEPRECATED: load the fit of the CoTiMA model
  saveHeterogeneityStudyFit=FALSE,        #DEPRECATED: save the fit of the heterogeneity model
  loadHeterogeneityStudyFit=FALSE,        #DEPRECATED: load the fit of the heterogeneity model
  loadFilePrefix=NULL,                    #DEPRECATED: the prefix of (older) files that contain model fits to be be loaded (useful if figureFilePrefix has changed since then)
  unixLike=FALSE,                         #DEPRECATED: instead the machine type is detected and parallelisation of loops with mclapply on unix machine activated
  figureFilePrefix=c()
)
{  # begin function definition (until end of file)

  #######################################################################################################################
  ############################################ load required subfunctions ###############################################
  #######################################################################################################################

  {
    #print(paste0("#################################################################################"))
    #print(paste0("########################### Load Required Functions #############################"))
    #print(paste0("#################################################################################"))

    #source(paste0(sourceDirectory, "CoTiMASaveFile.R"))
    #source(paste0(sourceDirectory, "CoTiMAfittingSSMF.R"))
    #source(paste0(sourceDirectory, "CoTiMActMultigroupFitAlt.R"))
    #source(paste0(sourceDirectory, "CoTiMAcomputeStartValues.R"))
    #source(paste0(sourceDirectory, "CoTiMAchangeStartValueLabels.R"))
    #source(paste0(sourceDirectory, "CoTiMAcompileListOfPrimaryStudies.R"))
    #source(paste0(sourceDirectory, "CoTiMApseudoRawData.R"))
    #source(paste0(sourceDirectory, "CoTiMAimportStartValues.R"))
    #source(paste0(sourceDirectory, "CoTiMAcombinePseudoRawData.R"))
    #source(paste0(sourceDirectory, "CoTiMAfitUserSpecifiedModel.R"))
  } ### END load required subfunctions ###



  #######################################################################################################################
  ####################################### install and attach further packages ###########################################
  #######################################################################################################################
  #{
    #print(paste0("#################################################################################"))
    #print(paste0("##################### Install and Attach Further Packages #######################"))
    #print(paste0("#################################################################################"))

    #if("MASS" %in% rownames(installed.packages()) == FALSE) {install.packages("MASS")} # for computing pseudo raw data
    #if("MBESS" %in% rownames(installed.packages()) == FALSE) {install.packages("MBESS")} # for computing statistical power
    #if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel")} # also for PP
    #if("rootSolve" %in% rownames(installed.packages()) == FALSE) {install.packages("rootSolve")} # for fast power calculation
    #if("crayon" %in% rownames(installed.packages()) == FALSE) {install.packages("crayon")}
    #if("psych" %in% rownames(installed.packages()) == FALSE) {install.packages("psych")}

    ## get NPSOL optimizer if requested
    #if( ("OpenMx" %in% rownames(installed.packages()) == TRUE) & (NPSOL==TRUE) ) { # if not already installed but requested get NPSOL
    #  library(OpenMx)
    #  if ("NPSOL" %in% mxAvailableOptimizers() == FALSE) { # if not already installed but requested get NPSOL
    #    source('https://openmx.ssri.psu.edu/software/getOpenMx.R')
    #  }        # ... version of OpenMx
    #}
    #if( ("OpenMx" %in% rownames(installed.packages()) == FALSE) & (NPSOL==TRUE) ) { # if not already installed but requested get NPSOL
    #  source('https://openmx.ssri.psu.edu/software/getOpenMx.R')
    #}        # ... version of OpenMx
    #if( ("OpenMx" %in% rownames(installed.packages()) == FALSE) & (NPSOL==FALSE) ) {install.packages("OpenMx")}   # if not already installed get standard version of OpenMx
    ## get ctsem (after OpenMx)
    #if("ctsem" %in% rownames(installed.packages()) == FALSE) {install.packages("ctsem")}

    ## If OpenMx and ctsem are already attached, detaching them is required to enable OpenMx to run in parallel mode.
    #if ("OpenMx" %in% (.packages())) suppressWarnings(detach("package:OpenMx", force=TRUE)) #, unload=TRUE)
    #if ("ctsem" %in% (.packages())) suppressWarnings(detach("package:ctsem", force=TRUE)) #, unload=TRUE))
    #Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
    #library(ctsem)
    #mxOption(key='Number of Threads', value=parallel::detectCores()) #now

    ## Attach all required packages that were not already attach with "library(ctsem)" before
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
  #} ### END install and attach further packages ###


  #######################################################################################################################
  ########################################### Check Model Specification #################################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("########################## Check Model Specification ############################"))
    print(paste0("#################################################################################"))


    if (is.null(primaryStudies)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold(" List with lists of primary study information not specified!", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat(red$bold("Should I try to get primary study information (e.g., empcov1, delta_t22) from the global environment?", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat(red$bold("(This is the way how earlier version of CoTiMA were run but no longer recommended.)", sep="\n"))
      cat(red$bold(" ", " ", sep="\n"))
      cat(blue("Press 'q' to quit and specify or 'c' to give it a try. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and specify primaryStudies or 'c' to try reading from global environment Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') {
        stop("Good luck for the next try!")
      } else {
        if (!(is.numeric(tryCatch(get(paste("sampleSize", 1, sep = "")), error = function(e) e)))) {
          cat(red$bold("Getting primary study information from the global environment failed", sep="\n"))
          cat(red$bold(" ", " ", sep="\n"))
          cat(blue("To test, I searched for sampleSize1 and could not find it!", "\n"))
          cat(red$bold(" ", " ", sep="\n"))
          stop("Good luck for the next try!")
        } else {
          cat(blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n"))
          char <- as.numeric(readline(""))
          while ((is.na(char))) {
            cat((blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n")))
            char <- as.numeric(readline(""))
          }
          primaryStudies <- compileListOfPrimaryStudies(selectedStudies=1:char)
        } # END else (catch failed)
      } # END else
    } #


    if (is.null(nlatents)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("Number of variables (nlatents) not specified!", sep="\n"))
      stop("Good luck for the next try!")
    }


    if (is.null(nlatents)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("Number of variables (nlatents) not specified!", sep="\n"))
      stop("Good luck for the next try!")
    }

    if (is.null(workingDirectory)) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("No working directory has been specified!", sep="\n"))
      stop("Good luck for the next try!")
    }

    if (length(statisticalPower) > 0) {
      for (i in 1:length(statisticalPower)) {
        if (statisticalPower[i] < 0 | statisticalPower[i] > 1){
          if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
          cat(red$bold("Parameter for statistical poweranalysis outside the allowed interval!\n"))
          cat(red("Values have to be between 0 and 1\n"))
          stop("Good luck for the next try!")
        }
      }
    }

    if (timeUnit=="timeUnit") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default time lag (timeUnit) has been chosen.", "\n"))
      cat(blue("Press 'q' to quit and specify or'c'to continue. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and specify time lag or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    if (resultsFileName=="CoTiMA.txt") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default output file name (CoTiMA.txt) has been chosen.", "\n"))
      cat(blue("Press 'q' to quit and change or'c'to continue. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change filename or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    if (filePrefix=="CoTiMA") {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The default figure file prefix (CoTiMA) has been chosen.", "\n"))
      cat(blue("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards "))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change filename or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    if (is.null(loadFilePrefix)) {ncharLoadFilePrefix <- 0} else {ncharLoadFilePrefix <- nchar(loadFilePrefix)} # check maximum path + filename length
    ncharResultsFileName <- nchar(resultsFileName)
    ncharFilePrefix <- nchar(filePrefix)
    ncharWorkingDirectory <- nchar(workingDirectory)
    ncharModelFits <- max( (nchar(saveSingleStudyModelFit[1])+2),
                           nchar(loadSingleStudyModelFit),
                           nchar(saveDRIFTAllModelFit),
                           nchar(loadDRIFTAllModelFit),
                           nchar(saveDRIFTSingleModelFit),
                           nchar(loadDRIFTSingleModelFit),
                           nchar(saveDRIFTCrossModelFit),
                           nchar(loadDRIFTCrossModelFit),
                           nchar(saveDRIFTAutoModelFit),
                           nchar(loadDRIFTAutoModelFit),
                           nchar(saveHeterogeneityModelFit),
                           nchar(loadHeterogeneityModelFit),
                           1) # 1 added to prevent empty vector
    if (nchar("_CoTiMA results for moderated auto-regressive (all effects in the model were moderated) effects of V2.png") +
        max(ncharLoadFilePrefix, ncharResultsFileName, ncharFilePrefix, ncharModelFits) +
        ncharWorkingDirectory > 259) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("The overall length of the path name (working directory) + the file names to be created exceeds 260.", "\n"))
      cat(blue("Shorten the file names or prefixes, or better work in a different working directory.", sep="\n"))
      stop("Good luck for the next try!")
    }

    if ( (testModeratorModel != FALSE) & is.null(moderatorNumber) ) {
      moderatorNumber <- 1
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("No moderator number has been specified. The first moderator will be used.", "\n"))
      cat(blue("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((blue("Please press 'q' to quit and change or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    if (length(saveSingleStudyFit) > 1 | length(loadSingleStudyFit) > 1 | saveHeterogeneityStudyFit == TRUE | loadHeterogeneityStudyFit == TRUE
        | saveHomogeneityStudyFit == TRUE | loadHomogeneityStudyFit == TRUE | !(is.null(loadFilePrefix))) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(red$bold("DEPRECATED parameters detected!","\n"))
      cat(red("You're using parameters of an older Version of CoTiMA (prior v1.6)","\n"))
      stop(red("Please adjust your file to the newest version or use an older function"))
    }

    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  {
        if (.Platform$OS.type == "unix") {
          coresToUse <- detectCores() + coresToUse
          cat("     #################################################################################\n")
          cat("     # You decided using multiple CPU cores for parallel fitting CoTiMA models. This #\n")
          cat("     # may work well in many cases, but we highly recommend not to run this script   #\n")
          cat("     # in RStudio. Instead, go to the console (e.g., console.app on MAC OS), then    #\n")
          cat("     # change to the direcrtory where your R-script is located (by typing, e.g.,     #\n")
          cat("     # \'cd \"/Users/trump/only/temp\"\'), and then, e.g., \'Rscript \"CoTiMA1.R\"\'.#\n")
          cat("     #################################################################################")
          Sys.sleep(5)
        } else {
          coresToUse <- 1
        }
      }
    } else {
      coresToUse <- 1
    }

    if ((-1) * coresToUse >= detectCores()) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Attention!"))}
      coresToUse <- detectCores() - 1
      cat(red("No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.","\n"))
    }
  } ### END Check Model Specification ###



  #######################################################################################################################
  ##################### Read user provided data and create list with all study information ##############################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("###### Read user provided data and create list with all study information #######"))
    print(paste0("#################################################################################"))

    #setwd(workingDirectory)
    start.time <- Sys.time(); start.time

    studyList <- list()
    maxLengthModeratorVector <- 0
    loadRawDataStudyNumbers <- unlist(lapply(primaryStudies$rawData,
                                             function(extract) extract$studyNumbers)); loadRawDataStudyNumbers

    # resorting primary study information
    if (!(exists("moderatorNumber"))) moderatorNumber <- 1; moderatorNumber
    for (i in 1:length(unlist(primaryStudies$studyNumbers))) {
      studyList[[i]] <- list(studyNumber=i, empcov=primaryStudies$empcovs[[i]], delta_t=primaryStudies$deltas[[i]],
                             sampleSize=primaryStudies$sampleSizes[[i]], originalStudyNo=primaryStudies$studyNumber[[i]],
                             timePoints=(length(primaryStudies$deltas[[i]]) + 1), moderators=primaryStudies$moderators[[i]],
                             maxModerators=length(primaryStudies$moderators[[i]]), startValues=primaryStudies$startValues[[i]],
                             rawData=primaryStudies$rawData[[i]], pairwiseN=primaryStudies$pairwiseNs[[i]])
      if (length(primaryStudies$moderators[[i]]) > maxLengthModeratorVector) maxLengthModeratorVector <- length(primaryStudies$moderators[[i]])
    }

    # check matrix symmetry if matrix is provided
    if (!(primaryStudies$studyNumbers[i] %in% loadRawDataStudyNumbers)) {
      if (isSymmetric(primaryStudies$empcovs[[i]]) == FALSE) {
        cat(red$bold("The correlation matrix of study no.", counter, "is not symmetric. Check and re-start!", "\n"))
        stop("Good luck for the next try!")
      }
    }

    # determine number of studies (if list is created by CoTiMAcompileListOfPrimaryStudies.R it is 1 element too long)
    tmp <- length(unlist(primaryStudies$studyNumbers)); tmp
    if ( is.na(primaryStudies$deltas[tmp]) & is.na(primaryStudies$sampleSizes[tmp]) & is.na(primaryStudies$deltas[tmp]) &
         is.null(dim(primaryStudies$pairwiseNs[tmp])) & is.null(dim(primaryStudies$emcovs[tmp])) &
         is.na(primaryStudies$moderators[tmp]) & is.null(primaryStudies$rawData$fileName[tmp]) ) {
      n.studies <- tmp-1
      primaryStudies$studyNumbers[[tmp]] <- NULL
    } else {
      n.studies <- tmp
    }
    n.studies

    if (testModeratorModel != FALSE) {
      if (maxLengthModeratorVector < moderatorNumber )  {
        if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        moderatorNumber <- 1
        cat(red$bold("The 'moderatorNumber' was automatically set to 1 because it was larger than the longest vector of moderator values provided.", "", sep="\n"))
        cat(blue("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards.", "\n"))
        char <- readline(" ")
        while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
          cat((blue("Please press 'q' to quit and change or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
      }
    }

    ### create pseudo raw data for all studies or load raw data if available & specified
    empraw <- lags <- moderators <- emprawMod <- allSampleSizes <- lostN <- overallNDiff <- relativeNDiff <- list()
    allTpoints <-lapply(studyList, function(extract) extract$timePoints); allTpoints
    #allSampleSizes <-lapply(studyList, function(extract) extract$sampleSize); allSampleSizes
    #allSampleSizes[length(allSampleSizes)] <- NULL
    #allSampleSizes <- unlist(allSampleSizes); allSampleSizes
    #if (length(allSampleSizes) > n.studies) allSampleSizes[tmp] <- NULL

    for (i in 1:n.studies) {
      if (!(i %in% loadRawDataStudyNumbers)) {
        currentVarnames <- c()
        currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        currentEmpcov <- (lapply(studyList, function(extract) extract$empcov))[[i]]; currentEmpcov
        currentLags <- (lapply(studyList, function(extract) extract$delta_t))[[i]]; currentLags
        currentPairwiseN <- (lapply(studyList, function(extract) extract$pairwiseN))[[i]]; currentPairwiseN
        if (testModeratorModel == TRUE) {
          currentModerators <- (lapply(studyList, function(extract) extract$moderators))[[i]]; currentModerators
        }

        for (j in 1:(currentTpoints)) {
          for (h in 1:nlatents) {
            currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
          }
        }

        tmp <- suppressWarnings(pseudoRawData(empCovMat=currentEmpcov, empN=currentSampleSize,
                                                      empNMat=currentPairwiseN))
        empraw[[i]] <- tmp$data
        lostN[[i]] <- tmp$lostN
        overallNDiff[[i]] <- tmp$overallLostN
        relativeNDiff[[i]] <- tmp$relativeLostN
        lags[[i]] <- matrix(currentLags, nrow=dim(empraw[[i]])[1], ncol=currentTpoints-1, byrow=TRUE)
        empraw[[i]] <- cbind(empraw[[i]], lags[[i]]);
        colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))
        empraw[[i]] <- as.data.frame(empraw[[i]])
        ## START correction of current lags if entire time point is missing for a case
        # wide to long
        emprawLong <- ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=nlatents, manifestNames=paste0("V", 1:nlatents))
        emprawLong <- suppressMessages(ctDeintervalise(datalong = emprawLong, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        emprawLong <- emprawLong[, ][ apply(emprawLong[, paste0("V", 1:nlatents)], 1, function(x) sum(is.na(x)) != nlatents ), ]
        # eliminate rows where time is NA
        emprawLong <- emprawLong[which(!(is.na(emprawLong[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctLongToWide(emprawLong, id='id', time='time', manifestNames=paste0("V", 1:nlatents)))
        # inrervalise
        emprawWide <- suppressMessages(ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=nlatents, manifestNames=paste0("V", 1:nlatents)))
        # restore
        empraw[[i]] <- as.data.frame(emprawWide)
        # END correction

        if (testModeratorModel == TRUE) {
          moderators[[i]] <- matrix(currentModerators, nrow=dim(empraw[[i]])[1], ncol=length(currentModerators), byrow=TRUE)
          emprawMod[[i]] <- cbind(empraw[[i]], moderators[[i]]);
          colnames(emprawMod[[i]]) <- c(c(currentVarnames,
                                          paste0("dT", seq(1:(currentTpoints-1))),
                                          paste0("mod", seq(1:length(currentModerators)))))
          emprawMod[[i]] <- as.data.frame(emprawMod[[i]])
        }
      }

      # load raw data on request
      if ( i %in% loadRawDataStudyNumbers ) {
        x1 <- studyList[[i]]$rawData$fileName; x1
        tmpData <- read.table(file=studyList[[i]]$rawData$fileName,
                              header=studyList[[i]]$rawData$header,
                              dec=studyList[[i]]$rawData$dec,
                              sep=studyList[[i]]$rawData$sep)

        # replace mising values
        tmpData <- as.matrix(tmpData) # important: line below will not work without having data as a matrix
        tmpData[tmpData %in% studyList[[i]]$rawData$missingValues] <- NA
        empraw[[i]] <- as.data.frame(tmpData)
        ## START correction of current lags if entire time point is missing for a case
        # change variable names
        tmp1 <- dim(empraw[[i]])[2]; tmp1
        currentTpoints <- (tmp1 + 1)/(nlatents+1); currentTpoints
        colnames(empraw[[i]])[1:(currentTpoints * nlatents)] <- paste0(paste0("V", 1:nlatents), "_T", rep(0:(currentTpoints-1), each=nlatents))
        # wide to long
        emprawLong <- ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=nlatents, manifestNames=paste0("V", 1:nlatents))
        emprawLong <- suppressMessages(ctDeintervalise(datalong = emprawLong, id='id', dT='dT'))
        # eliminate rows where ALL latents are NA
        emprawLong <- emprawLong[, ][ apply(emprawLong[, paste0("V", 1:nlatents)], 1, function(x) sum(is.na(x)) != nlatents ), ]
        # eliminate rows where time is NA
        emprawLong <- emprawLong[which(!(is.na(emprawLong[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctLongToWide(emprawLong, id='id', time='time', manifestNames=paste0("V", 1:nlatents)))
        # inrervalise
        emprawWide <- suppressMessages(ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=nlatents, manifestNames=paste0("V", 1:nlatents)))
        # restore
        empraw[[i]] <- as.data.frame(emprawWide)
        # END correction

        # Change the NAs provided for deltas if raw data are loaded
        for (h in 1:(currentTpoints-1)) studyList[[i]]$delta_t[h] <- mean(empraw[[i]][, paste0("dT", 1)], na.rm=TRUE)

      }

      # change sample size if entire cases were deleted
      studyList[[i]]$sampleSize <- (dim(empraw[[i]]))[1]
      allSampleSizes[[i]] <- dim(empraw[[i]])[1]

      currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
      currentTpoints <- allTpoints[[i]]; currentTpoints
      currentVarnames <- c()
      for (j in 1:(currentTpoints)) {
        for (h in 1:nlatents) {
          currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
        }
      }
      colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))

      # standardize (variables - not time lags) if option is chosen
      if (studyList[[i]]$rawData$standardize == TRUE) empraw[[i]][, currentVarnames] <- scale(empraw[[i]][, currentVarnames])

      # replace missing values for time lags dTi by .000001 (has to be so because dTi are definition variables)
      tmpData <- empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))]
      tmpData[is.na(tmpData)] <- .001
      empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))] <- tmpData

      if (testModeratorModel == TRUE) {
        #currentModerators <- studyList[[i]]$maxModerators; currentModerators
        currentModerators <- studyList[[i]]$moderators; currentModerators
        moderators[[i]] <- matrix(currentModerators, nrow=dim(empraw[[i]])[1], ncol=length(currentModerators), byrow=TRUE)
        emprawMod[[i]] <- as.data.frame(cbind(as.matrix(empraw[[i]]), as.matrix(moderators[[i]])))
        #colnames(emprawMod[[i]]) <- c(colnames(empraw[[i]]), paste0("mod", seq(1:currentModerators)))
        colnames(emprawMod[[i]]) <- c(colnames(empraw[[i]]), paste0("mod", seq(1:length(currentModerators))))
      }

      # Save raw data  on request
      if ( i %in% saveRawData$studyNumbers ) {
        x1 <- paste0(saveRawData$fileName, i, ".dat"); x1
        write.table(empraw[[i]], file=x1, row.names=saveRawData$row.names, col.names=saveRawData$col.names,
                    sep=saveRawData$sep, dec=saveRawData$dec)
      }
    }
  } ### END Read user provided data and create list with all study information ###


  #######################################################################################################################
  ################################################### Some Statistics ###################################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("################# Compute Summary Statistics of Primary Studies #################"))
    print(paste0("#################################################################################"))

    ### some stats
    # Sample size
    allSampleSizes <-unlist(lapply(studyList, function(extract) extract$sampleSize)); allSampleSizes
    allSampleSizes <- allSampleSizes[-length(allSampleSizes)]; allSampleSizes
    overallSampleSize <- sum(allSampleSizes, na.rm=TRUE); overallSampleSize
    meanSampleSize <- mean(allSampleSizes, na.rm=TRUE); meanSampleSize
    maxSampleSize <- max(allSampleSizes, na.rm=TRUE); maxSampleSize
    minSampleSize <- min(allSampleSizes, na.rm=TRUE); minSampleSize
    # lags
    # Will be computed again below if raw data are loaded #
    allDeltas <-unlist(lapply(studyList, function(extract) extract$delta_t)); allDeltas
    allDeltas <- allDeltas[-length(allDeltas)]; allDeltas
    meanDelta <- mean(allDeltas, na.rm=TRUE); meanDelta
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    minDelta <- min(allDeltas, na.rm=TRUE); minDelta
    # Time points
    allTpoints <-unlist(lapply(studyList, function(extract) extract$timePoints)); allTpoints
    overallTpoints <- sum(allTpoints, na.rm=TRUE); overallTpoints
    meanTpoints  <- mean(allTpoints, na.rm=TRUE); meanTpoints
    maxTpoints  <- max(allTpoints, na.rm=TRUE); maxTpoints
    minTpoints  <- min(allTpoints, na.rm=TRUE); minTpoints
    # Moderators
    if (testModeratorModel == TRUE) {
      allModerators <-unlist(lapply(studyList, function(extract) extract$moderators[moderatorNumber])); allModerators
      allModerators <- allModerators[-length(allModerators)]; allModerators
      numberOfModerators <- length(table(allModerators)); numberOfModerators
      categoriesOfModerators <- names(table(allModerators)); categoriesOfModerators
      maxModerators <-unlist(lapply(studyList, function(extract) extract$maxModerators))
      maxModerators <- maxModerators[-length(maxModerators)]; maxModerators
      maxModerators  <- max(maxModerators); maxModerators
    }
  } ### END Some Statistics ###



  #######################################################################################################################
  ############################# Create ctsem model template to fit all primary studies ##################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("############ Create ctsem Model Template to fit all Primary Studies #############"))
    print(paste0("#################################################################################"))

    manifestNames <- c()
    for (i in 1:nlatents) manifestNames <- c(manifestNames, paste0("V",i))
    driftNames <- c()
    for (i in 1:(nlatents)) {
      for (j in 1:(nlatents)) {
        driftNames <- c(driftNames, paste0("V",i,"toV", j))
      }
    }

    cintNames <- c()
    for (i in 1:(nlatents)) {
      cintNames <- c(cintNames, paste0("cint",i))
    }

    # general ctsem model template
    ctsemModelTemplate <- ctModel(n.latent=nlatents, n.manifest=nlatents, Tpoints=2, manifestNames=manifestNames,    # 2 waves in the template only
                                  DRIFT=matrix(driftNames, nrow=nlatents, ncol=nlatents),
                                  LAMBDA=diag(nlatents),
                                  #CINT=matrix(cintNames, nrow=nlatents, ncol=1),
                                  CINT=matrix(0, nrow=nlatents, ncol=1),
                                  T0MEANS = matrix(c(0), nrow = nlatents, ncol = 1),
                                  MANIFESTMEANS = matrix(c(0), nrow = nlatents, ncol = 1),
                                  MANIFESTVAR=matrix(0, nrow=nlatents, ncol=nlatents))
    if (manifesttraitvar == TRUE) {
      ctsemModelTemplate <- ctModel(n.latent=nlatents, n.manifest=nlatents, Tpoints=2, manifestNames=manifestNames,    # 2 waves in the template only
                                    DRIFT=matrix(driftNames, nrow=nlatents, ncol=nlatents),
                                    LAMBDA=diag(nlatents),
                                    #CINT=matrix(cintNames, nrow=nlatents, ncol=1),
                                    CINT=matrix(0, nrow=nlatents, ncol=1),
                                    T0MEANS = matrix(c(0), nrow = nlatents, ncol = 1),
                                    MANIFESTMEANS = matrix(c(0), nrow = nlatents, ncol = 1),
                                    MANIFESTVAR=matrix(0, nrow=nlatents, ncol=nlatents),
                                    MANIFESTTRAITVAR = 'auto')
    }

    # ctsem models for each primary study with the correct number of time points
    ctsemModel <- list()
    counter <- 1
    for (i in unique(unlist(allTpoints))) {
      helper <- ctsemModelTemplate
      helper$Tpoints <- i
      ctsemModel[[counter]] <- helper
      counter <- counter +1
    }
  } ### END Create ctsem model template to fit all primary studies ###



  #######################################################################################################################
  ##################################### Check Specification of Primary Studies ##########################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("#################### Check Specification of Primary Studies #####################"))
    print(paste0("#################################################################################"))

    if (length(saveSingleStudyModelFit) == 1){
      {
        if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        cat(red$bold("You have indicated that you want to save SingleStudyModelFits, but have not selected any study/studies to save the fit for!","\n"))
        cat(red("Would you like to save the SingleStudyModelFits for ALL studies??","\n"))
        cat(blue("Press 'y' to save ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(blue("select the study/studies you whish to save the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((blue("Please press 'y' to save ALL, 's' to specify the study/studies to save, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in 1:n.studies) {
            counter <- i
            if (length(saveSingleStudyModelFit > 0)){
              saveSingleStudyModelFit <- c(saveSingleStudyModelFit, counter)
              i <- i+1
            }
          }
        } else if (char == 's' | char == 'S') {
          cat(blue("Which SingleStudyModelFits would you like to save?", "\n"))
          cat(blue("Please enter the no. of study/studies separated by commas!", "\n"))
          chars <- readline(" ")
          chars <- gsub(" ", "", chars, fixed = TRUE)
          chars <- unlist(strsplit(chars, ","))
          chars
          saveSingleStudyModelFit <- c(saveSingleStudyModelFit, chars)
        } else if (char == 'q' | char == 'Q'){
          stop("Good luck for the next try!")
        }
      }

      if (length(loadSingleStudyModelFit) == 1){
        if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        cat(red$bold("You have indicated that you want to load SingleStudyModelFits, but have not selected any study/studies to load the fit for!","\n"))
        cat(red("Would you like to load the SingleStudyModelFits for ALL studies??","\n"))
        cat(blue("Press 'y' to load ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(blue("select the study/studies you whish to load the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((blue("Please press 'y' to load ALL, 's' to specify the study/studies to load, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in 1:n.studies) {
            counter <- i
            if (length(loadSingleStudyModelFit > 0)){
              loadSingleStudyModelFit <- c(loadSingleStudyModelFit, counter)
              i <- i+1
            }
          }
        } else if (char == 's' | char == 'S') {
          cat(blue("Which SingleStudyModelFits would you like to load?", "\n"))
          cat(blue("Please enter the no. of study/studies separated by commas!", "\n"))
          chars <- readline(" ")
          chars <- gsub(" ", "", chars, fixed = TRUE)
          chars <- unlist(strsplit(chars, ","))
          loadSingleStudyModelFit <- c(loadSingleStudyModelFit, chars)
        } else if (char == 'q' | char == 'Q'){
          stop("Good luck for the next try!")
        }
      }
    }
  } ### END check specification of primary studies ###



  #######################################################################################################################
  ##################################### Fit ctsem Model to each Primary Study ###########################################
  #######################################################################################################################
  {
    # loop through all primary studies
    FitList=list()

    studyFit <- studyFitCI <- studyFit_Minus2LogLikelihood <- studyFit_estimatedParameters <- list()
    study_Drift_Coef <- study_Drift_SE <- list()
    study_Diffusion_Coef <- study_Diffusion_SE <- list()
    study_T0var_Coef <- study_T0var_SE <- list()
    study_Cint_Coef <- study_Cint_SE <- list()
    origRefits <- refits; origRefits

    for (i in 1:n.studies) {
      if ( (length(loadSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) {
        print(paste0("#################################################################################"))
        print(paste0("######################## LOADING SingleStudyFit ", i, " of ", n.studies, " ##########################"))
        print(paste0("#################################################################################"))
        x1 <- paste0(workingDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        studyFit[[i]] <- readRDS(file=x1)

        # STORE
        FitList$tmpName <-studyFit[[i]]
        tmp <- names(FitList); tmp
        tmp[length(tmp)] <- paste0("SingleStudyFit", i); tmp
        names(FitList) <- tmp; names(FitList)
      }

      if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1])) {
        print(paste0("#################################################################################"))
        print(paste0("################### Fitting SingleStudyModel ", i, " of ", n.studies, " #######################"))
        print(paste0("#################################################################################"))
        if (studyList[[i]]$originalStudyNo %in% extraRefits) refits <- factorExtraRefits * refits
        if (!(studyList[[i]]$originalStudyNo %in% extraRefits)) refits <- origRefits
        # select correct template
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        modelToSelect <- which(unique(allTpoints) == currentTpoints); modelToSelect
        currentModel <- ctsemModel[[modelToSelect]]; currentModel

        # Set Start Values
        currentStartValues <- (lapply(studyList, function(extract) extract$startValues))[[i]]; currentStartValues
        if (is.na(currentStartValues[1])) currentStartValues <- NULL
        studyList

        # FIT (NOT LOAD)
        mxOption(NULL, 'Number of Threads', 1)
        results <- CoTiMAfittingSSMF(coresToUse=coresToUse, empraw=empraw[[i]], currentModel=currentModel,
                                     refits=refits, retryattempts=retryattempts,
                                     singleModelStartValues=currentStartValues)
        mxOption(key='Number of Threads', value=parallel::detectCores())

        # Extract estimates & statistics
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        studyFit[[i]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]

        # STORE
        FitList$tmpName <-studyFit[[i]]
        tmp <- names(FitList); tmp
        tmp[length(tmp)] <- paste0("SingleStudyFit", studyList[[i]]$originalStudyNo); tmp
        names(FitList) <- tmp; names(FitList)

        # SAVE
        if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
          x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
          x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/")
          CoTiMASaveFile(activateRPB, workingDirectory, studyFit[[i]], x1, x2)
        }
      }

      studyFit_Minus2LogLikelihood[[i]] <- (summary(studyFit[[i]]))$omxsummary$Minus2LogLikelihood
      studyFit_estimatedParameters[[i]] <- (summary(studyFit[[i]]))$omxsummary$estimatedParameters
      study_Drift_Coef[[i]] <- (summary(studyFit[[i]]))$ctparameters[1:(nlatents^2), 1]
      study_Drift_SE[[i]] <- (summary(studyFit[[i]]))$ctparameters[1:(nlatents^2), 3]
      study_Diffusion_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("diffusion", names(studyFit[[1]]$mxobj$output$estimate))]
      study_Diffusion_SE[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("diffusion", rownames(studyFit[[1]]$mxobj$output$standardErrors))]
      study_T0var_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("T0var", names(studyFit[[1]]$mxobj$output$estimate))]
      study_T0var_SE[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("T0var", rownames(studyFit[[1]]$mxobj$output$standardErrors))]
      study_Cint_Coef[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("cint", names(studyFit[[1]]$mxobj$output$estimate))]
      study_Cint_SE[[i]] <- studyFit[[i]]$mxobj$output$estimate[grep("cint", rownames(studyFit[[1]]$mxobj$output$standardErrors))]
    }

    # Compute confidence intervals
    if (confidenceIntervals == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("######################## Computing Confidence Intervals #########################"))
      print(paste0("#################################################################################"))

      # LOAD
      if ( (length(loadSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) {
        for (i in 1:n.studies) {
          print(paste0("#################################################################################"))
          print(paste0("############ LOADING SingleStudyFit with confidence intervals ", i, " of ", n.studies, " ############"))
          print(paste0("#################################################################################"))
          x1 <- paste0(workingDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFitCI", studyList[[i]]$originalStudyNo, ".rds"); x1
          studyFitCI[[i]] <- readRDS(file=x1)

          # STORE
          FitList$tmpName <-studyFitCI[[i]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("SingleStudyFitCI", i); tmp
          names(FitList) <- tmp; names(FitList)
        }
      }

      # FIT (NOT LOAD)
      if ( (length(loadSingleStudyModelFit) < 1) ) {
        mxOption(NULL, 'Number of Threads', 1)
        results <- mclapply(seq(1, n.studies, by=1),
                            function(studyNo) ctCI(ctfitobj=studyFit[[studyNo]], confidenceintervals = driftNames),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())

        # STORE
        for (k in 1:n.studies) {
          # fits
          #FitList$tmpName <- results[[k]]
          #tmp <- names(FitList); tmp
          #tmp[length(tmp)] <- paste0("SingleStudyFitCI", k); tmp
          #names(FitList) <- tmp; names(FitList)
          # just confidence intervals
          FitList$tmpName <- results[[k]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("SingleStudyFitCI", k); tmp
          names(FitList) <- tmp; names(FitList)
        }

        # SAVE
        studyFitCI <- results; studyFitCI
        for (h in 1:n.studies) {
          if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[h]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
            x1 <- paste0(saveSingleStudyModelFit[1], " studyFitCI", studyList[[h]]$originalStudyNo, ".rds"); x1
            x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/")
            CoTiMASaveFile(activateRPB, workingDirectory, results[[h]], x1, x2)
          }
        }
      }
    } # END Computing Confidence Intervals

    # Combine summary information and fit statistics
    allStudies_Minus2LogLikelihood <- sum(unlist(studyFit_Minus2LogLikelihood)); allStudies_Minus2LogLikelihood
    allStudies_estimatedParameters <- sum(unlist(studyFit_estimatedParameters)); allStudies_estimatedParameters
    allStudies_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
      allStudies_estimatedParameters; allStudies_df
    if (confidenceIntervals == TRUE) allDriftCI <-lapply(studyFitCI, function(extract) extract$mxobj$output$confidenceIntervals)
    allStudiesDRIFT_effects <- matrix(t(cbind(unlist(study_Drift_Coef), unlist(study_Drift_SE)) ), n.studies, 2*nlatents^2, byrow=T)

    # Label summary table
    rownames(allStudiesDRIFT_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    newColNames <- c()
    for (j in 1:nlatents) {
      for (h in 1:nlatents) {
        newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)")
      }
    }
    colnames(allStudiesDRIFT_effects) <- newColNames; round(allStudiesDRIFT_effects, digits)
    FitList$allStudiesDRIFT_effects <- allStudiesDRIFT_effects

    # check single study results
    if (checkSingleStudyResults == TRUE) {
      print(round(allStudiesDRIFT_effects, digits))
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(blue(" Press 'q' to quit or any other key to continue. Press ENTER afterwards."))
      char <- readline(" ")
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    # Extract parameter estimates and paramter names to be used as start values later
    est <- unlist(lapply(studyFit, function(extract) extract$mxobj$output$estimate)); est
    startValues <- c(t(est)); startValues
    noOfParams <- length(studyFit[[1]]$mxobj$output$estimate); noOfParams
    headNames <- c(matrix(paste0(rep(c(paste0("Study_No_", primaryStudies$studyNumbers, "_")), noOfParams)),
                          nrow=noOfParams, ncol=n.studies, byrow=T)); headNames
    names(startValues) <- paste0(headNames, names(est)); startValues

    DRIFTCoeff <- matrix(unlist(study_Drift_Coef), n.studies, nlatents^2, byrow=TRUE); DRIFTCoeff
    DRIFTSE <- matrix(unlist(study_Drift_SE), n.studies, nlatents^2, byrow=TRUE); DRIFTSE

    if (n.studies < 2) {
      if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(blue("Only a single primary study was handed over to CoTiMA. No further (meta-) analyses can be conducted."))
      cat(" ", sep="\n")
      cat(blue("I guess this stop is intended! You could ignore further warning messages such as \"sqrt(n-2): NaNs produced\""))
      cat(" ", sep="\n")
      stop("Hopefully the solution is reasonable. Set refits to a higher value if this is not the case.")
      cat(" ", sep="\n")

    }

  } ### END fit ctsem model to each primary study


  #######################################################################################################################
  ################## Specification of Parameters for Plotting, Statistical Power, Optimal Lags  #########################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("## Specification of Parameters for Plotting, Statistical Power, & Optimal Lags ##"))
    print(paste0("#################################################################################"))

    ## timeRange, stepWidth, & noOfSteps
    if (length(timeRange) < 1) {
      stepWidth <- 1
      usedTimeRange <- seq(0, 1.5*maxDelta, stepWidth)
      # add empirical lags not yet included
      usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
      noOfSteps <- length(usedTimeRange)
    } else {
      stepWidth <- timeRange[3]
      usedTimeRange <- seq(timeRange[1], timeRange[2], stepWidth)
      # add empirical lags not yet included
      usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
      noOfSteps <- length(usedTimeRange)
    }

    ## yLimitsForEffects
    # discrete effects across time range
    allSingleStudyDiscreteDrift <- array(dim=c(n.studies, noOfSteps, nlatents^2))  # primary studies
    for (h in 1:n.studies) {
      for (i in usedTimeRange[1]:noOfSteps){
        timeValue <- i * stepWidth; timeValue
        allSingleStudyDiscreteDrift[h, i, 1:(nlatents^2)] <- c(expm(matrix(study_Drift_Coef[[h]], nlatents, nlatents) %x% timeValue))
      }
    }
    # min, max, & yLimitsForEffects
    coeffSeq <- seq(1, nlatents^2, 1)[!(seq(1, nlatents^2, 1) %in% seq(1, nlatents^2, (nlatents+1)))]; coeffSeq
    if (length(yLimitsForEffects) < 1) {
      yLimitsForEffects[1] <- round(min(allSingleStudyDiscreteDrift[, , coeffSeq]) - .10, 1); yLimitsForEffects[1]
      yLimitsForEffects[2] <- round(max(allSingleStudyDiscreteDrift[, , coeffSeq]) + .10, 1); yLimitsForEffects[2]
    }

    ## yLimitsForPower
    if (length(yLimitsForPower) < 1) {
      yLimitsForPower[1] <- 0
      yLimitsForPower[2] <- 3*10^(nchar(maxSampleSize)-1)
    }
  } ### END Specification of Parameters for Plotting, Statistical Power, Optimal Lags ###



  #######################################################################################################################
  ####################################### Data Management for CoTiMA ####################################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("########################## Data Management for CoTiMA ###########################"))
    print(paste0("#################################################################################"))

    empraw2 <- emprawMod2 <- list()
    currentVarnames <- c()
    for (j in 1:(maxTpoints)) {
      for (h in 1:nlatents) {
        currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
      }
    }

    for (j in 1:(maxTpoints-1)) currentVarnames <- c(currentVarnames, paste0("dT", j))
    if (testModeratorModel == TRUE) {
      currentVarnamesMod <- currentVarnames
      for (j in 1:(maxModerators)) currentVarnamesMod <- c(currentVarnamesMod, paste0("mod", j))
    }

    for (i in 1:n.studies) {
      addTpoints <- maxTpoints - allTpoints[[i]]; addTpoints
      addVariables <- matrix(NA, nrow=allSampleSizes[[i]], ncol=(2*addTpoints)); dim(addVariables)
      addLags <- matrix(rep(.00001, addTpoints), nrow=allSampleSizes[[i]], ncol=addTpoints); dim(addLags)
      empraw2[[i]]  <- cbind(empraw[[i]][ , 1:(nlatents*allTpoints[[i]])], addVariables,
                             empraw[[i]][ , (nlatents*allTpoints[[i]]+1):(nlatents*allTpoints[[i]]+(allTpoints[[i]]-1))], addLags)
      colnames(empraw2[[i]]) <- currentVarnames
      if (testModeratorModel == TRUE) {
        addModerators <- matrix(moderators[[i]], nrow=allSampleSizes[[i]], ncol=dim(moderators[[i]])[2]); addModerators
        emprawMod2[[i]] <- cbind(empraw2[[i]], addModerators)
        colnames(emprawMod2[[i]]) <- currentVarnamesMod
      }
    }

    ### tryout
    # make data & grouping variabe
    if (!(exists("allModerators"))) allModerators <- 1
    tmp <- combinePseudoRawData(listOfStudyFits=studyFit, moderatorValues=allModerators)
    datawide_all <- tmp$alldata
    groups <- tmp$groups
    names(groups) <- c("Study_No_"); groups
    groupsNamed <- (paste0("Study_No_", groups)); groupsNamed
    moderatorGroups <- tmp$moderatorGroups

    ### try skipping and test tryout
    trySkip <- 1
    if (trySkip != 1)  {
      # combine datasets
      datawide_all <- matrix(NA, 0, ncol=dim(empraw2[[1]])[2]); datawide_all
      if (testModeratorModel == TRUE) {
        datawideMod_all <- matrix(NA, 0, ncol=dim(emprawMod2[[1]])[2]); datawideMod_all
      }
      for (i in 1:n.studies) {
        datawide_all <- rbind(datawide_all, empraw2[[i]])
        if (testModeratorModel == TRUE) datawideMod_all <- rbind(datawideMod_all, emprawMod2[[i]])
      }

      colnames(datawide_all) <- colnames(empraw[[min(which(max(unique(unlist(allTpoints)))==unlist(allTpoints)))]])
      if (testModeratorModel == TRUE) {
        colnames(datawideMod_all) <- colnames(emprawMod[[min(which(max(unique(unlist(allTpoints)))==unlist(allTpoints)))]])
      }

      # create grouping variable
      #studyNo <-unlist(rep(primaryStudies$studyNumbers, allSampleSizes))
      #if (testModeratorModel == TRUE) datawideMod_all <- as.data.frame(cbind(datawideMod_all, studyNo))
      groups <- c(rep(primaryStudies$studyNumbers, allSampleSizes))
      groups <- unlist(groups)
      names(groups) <- c("Study_No_"); groups
      groupsNamed <- (paste0("Study_No_", groups)); groupsNamed

      if (testModeratorModel == TRUE) groupsMod <- unlist(rep(allModerators, allSampleSizes))
    }

  } ### END Data Management for CoTiMA ###



  #######################################################################################################################
  ##################################### Create Model Templates for Multigroup Fit #######################################
  #######################################################################################################################
  {
    print(paste0("#################################################################################"))
    print(paste0("################## Create Model Templates for Multigroup Fit ####################"))
    print(paste0("#################################################################################"))

    # General preparation for heterogeneity and homogeneity models
    hetModel <- ctsemModel[[which(unique(unlist(allTpoints)) == max(unique(unlist(allTpoints))))]]

    allFixedModel <- hetModel
    allFixedModel$DRIFT <- matrix("groupfixed", nlatents, nlatents)
    allFixedModel$DIFFUSION[grep("diffusion", allFixedModel$DIFFUSION)] <- "groupfixed"
    allFixedModel$T0VAR[grep("T0var", allFixedModel$T0VAR)] <- "groupfixed"
    if (any(hetModel$CINT != 0) ) {
      allFixedModel$CINT <- matrix("groupfixed", nlatents, 1)
    } else {
      allFixedModel$CINT <- matrix(0, nlatents, 1)
    }

  } ### END Create Model Templates for Multigroup Fit ###


  #######################################################################################################################
  ############################################### Heterogeneity Model ###################################################
  #######################################################################################################################
  {
    if ((testHeterogeneityModel == TRUE) & (length(loadHeterogeneityModelFit) < 1)) {
      print(paste0("#################################################################################"))
      print(paste0("######### Fitting Heterogeneity Model (Not Recommended - May Take Ages) #########"))
      print(paste0("#################################################################################"))

      # get start values
      hetModelCTmodelobj <- hetModel
      hetModelFixedModel <- NULL
      hetModelModModel <- NULL
      hetModelStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                ctModelObj = hetModelCTmodelobj, fixedModel = hetModelFixedModel,
                                                modModel = hetModelModModel)
      # FIT (NOT LOAD)
      if (useCTMultigroupFitAlt == FALSE) {
        mxOption(NULL, 'Number of Threads', 1)
        hetModelStartValues <- changeStartValueLabels(startValues = hetModelStartValues, ctmodelobj = hetModelCTmodelobj,
                                                      fixedmodel = hetModelCTmodelobj, noOfStudies <- n.studies)
        results <- mclapply(seq(1, refits, by=1),
                            function(refits) ctMultigroupFit(dat=datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                             omxStartValues=hetModelStartValues,
                                                             ctmodelobj = hetModelCTmodelobj),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        hetModelFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
        hetModelFit <- hetModelFit$mxobj
      } else {
        mxOption(NULL, 'Number of Threads', 1)
        results <- mclapply(seq(1, refits, by=1),
                            function(refits) ctMultigroupFitAlt(dat=datawide_all, groupings = groups,
                                                                retryattempts = retryattempts,
                                                                #startValues = studyFit,
                                                                startValues = hetModelStartValues,
                                                                ctmodelobj = hetModelCTmodelobj),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
        hetModelFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
      }
      FitList$HeterogeneityModelFit <- hetModelFit
    }

    # LOAD
    if ((testHeterogeneityModel == TRUE) & (length(loadHeterogeneityModelFit) > 0) ) {
      x1 <- paste0(loadHeterogeneityModelFit[1], " hetFit.rds"); x1
      hetModelFit <- readRDS(file=x1)
      FitList$HeterogeneityModelFit <- hetModelFit
    }

    # SAVE
    if ((testHeterogeneityModel == TRUE) & (length(saveHeterogeneityModelFit) > 0) ) {
      x1 <- paste0(saveHeterogeneityModelFit[1], " hetFit.rds"); x1
      x2 <- paste0(workingDirectory)
      CoTiMASaveFile(activateRPB, workingDirectory, hetModelFit, x1, x2)
    }

    if (testHeterogeneityModel == TRUE) {
      hetModel_Drift_Coef <- hetModelFit$output$estimate[grep("toV", names(hetModelFit$output$estimate))]
      hetModel_Drift_SE <- hetModelFit$output$estimate[grep("toV", rownames(hetModelFit$output$standardErrors))]
      hetModel_Minus2LogLikelihood <- hetModelFit$output$Minus2LogLikelihood; hetModel_Minus2LogLikelihood
      hetModel_estimatedParameters <- length(hetModelFit$output$estimate); hetModel_estimatedParameters
      hetModel_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
        hetModel_estimatedParameters; hetModel_df

      # Compute confidence intervals (may take ages)
      if (confidenceIntervals == TRUE) {
        hetDRIFTallFitCI <- hetDRIFTallCI <- list()
        print(paste0("#################################################################################"))
        print(paste0(" Computing Confidence Intervals for Heterogeneity Model (Not Recommended Either) "))
        print(paste0("#################################################################################"))

        # LOAD
        if (length(loadHeterogeneityModelFit) > 0) {
          hetModelFitCI <- list()
          x1 <- paste0(workingDirectory, loadHeterogeneityModelFit[1], " hetFitCI.rds"); x1
          hetModelFitCI <- readRDS(file=x1)
          hetModelCI <- hetModelFitCI$output$confidenceIntervals
        }

        # FIT (NOT LOAD)
        if (length(loadHeterogeneityModelFit) < 1) {
          hetModelCItmp <- hetModelFit
          ci <- mxCI(names(hetModelCItmp$output$estimate)[grep("to", names(hetModelCItmp$output$estimate))]); ci
          hetModelCItmp <- mxModel(hetModelCItmp, ci)
          hetModelFitCI <- mxRun(hetModelCItmp, intervals=TRUE)
          hetModelCI <- hetModelFitCI$output$confidenceIntervals
        }

        # STORE
        FitList$hetModelFitCI <- hetModelFitCI
        FitList$hetModelCI <- hetModelCI

        # SAVE
        if (length(saveHeterogeneityModelFit) > 0) {
          x1 <- paste0(saveHeterogeneityModelFit, loadHeterogeneityModelFit[1], " hetFitCI.rds"); x1
          x2 <- paste0(workingDirectory); x2
          CoTiMASaveFile(activateRPB, workingDirectory, hetModelFitCI, x1, x2)
        }

      } ## END confidence intervals

      # Combine summary information and fit statistics
      hetModelDRIFT_effects <- matrix(t(cbind((hetModel_Drift_Coef), (hetModel_Drift_SE)) ), n.studies, 2*nlatents^2, byrow=T)
      # Label summary table
      rownames(hetModelDRIFT_effects) <- paste0("Study No ", 1:n.studies)
      newColNames <- c()
      for (j in 1:nlatents) {
        for (h in 1:nlatents) {
          newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)")
        }
      }
      colnames(hetModelDRIFT_effects) <- newColNames; hetModelDRIFT_effects
      FitList$hetModelDRIFT_effects <- hetModelDRIFT_effects
    }
  } ### END Heterogeneity Model ###



  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################
  {
    if ((testDRIFTallModel == TRUE) & (length(loadDRIFTAllModelFit) < 1)) {
      print(paste0("#################################################################################"))
      print(paste0("## Fitting CoTiMA Model with all Drift Coefficients Invariant (DRIFTallModel) ###"))
      print(paste0("#################################################################################"))

      # get start values
      homDRIFTallCTmodelobj <- hetModel
      homDRIFTallFixedModel <- hetModel
      homDRIFTallFixedModel$DRIFT <- matrix("groupfixed", nlatents, nlatents)
      homDRIFTallModModel <- NULL
      homDRIFTallStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                   ctModelObj = homDRIFTallCTmodelobj, fixedModel = homDRIFTallFixedModel,
                                                   modModel = homDRIFTallModModel)

      # FIT (NOT LOAD)
      if (useCTMultigroupFitAlt == FALSE) {
        mxOption(NULL, 'Number of Threads', 1)
        homDRIFTallStartValues <- changeStartValueLabels(startValues = homDRIFTallStartValues, ctmodelobj = homDRIFTallCTmodelobj,
                                                         fixedmodel = homDRIFTallFixedModel, noOfStudies <- n.studies)
        results <- mclapply(seq(1, refits, by=1),
                            function(refits) ctMultigroupFit(dat=datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                             omxStartValues=homDRIFTallStartValues,
                                                             ctmodelobj = homDRIFTallCTmodelobj,
                                                             fixedmodel = homDRIFTallFixedModel),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
        homDRIFTallFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
        homDRIFTallFit <- homDRIFTallFit$mxobj
        FitList$DRIFTAllModelFit <- homDRIFTallFit
      } else {
        mxOption(NULL, 'Number of Threads', 1)
        results <- mclapply(seq(1, refits, by=1),
                            function(refits) ctMultigroupFitAlt(dat=datawide_all, groupings = groups,
                                                                retryattempts = retryattempts,
                                                                startValues = homDRIFTallStartValues,
                                                                ctmodelobj = homDRIFTallCTmodelobj,
                                                                fixedmodel = homDRIFTallFixedModel),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
        homDRIFTallFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
        FitList$DRIFTAllModelFit <- homDRIFTallFit
      }
    }

    # LOAD
    if ((testDRIFTallModel == TRUE) & (length(loadDRIFTAllModelFit) > 0) ) {
      x1 <- paste0(loadDRIFTAllModelFit[1], " homDRIFTallFit.rds"); x1
      homDRIFTallFit <- readRDS(file=x1)
      FitList$DRIFTAllModelFit <- homDRIFTallFit
    }

    if ((testDRIFTallModel == TRUE) | (length(loadDRIFTAllModelFit) > 0)) {
      if ( any(homDRIFTallFit$output$standardErrors > 1) |
           any(is.na(homDRIFTallFit$output$standardErrors)) |
           any(abs(homDRIFTallFit$output$estimate) > 3) |
           any(is.na(homDRIFTallFit$output$estimate)) ) {
        if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        print(cbind((homDRIFTallFit$output$estimate), (homDRIFTallFit$output$standardErrors)))
        cat(red$bold(" ", "", sep="\n"))
        cat(red$bold("Some estimates or their standard errors seem to be out of range!", sep="\n"))
        cat(red$bold("Optimal time lags will not be computed if you continue.", sep="\n"))
        cat(red$bold(" ", "", sep="\n"))
        cat(blue("Press 'q' to quit and specify or'c'to continue. Press ENTER afterwards ", "\n"))
        char <- readline(" ")
        while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
          cat((blue("Please press 'q' to quit and 'c' to continue without changes. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
      }
    }


    # SAVE
    if ((testDRIFTallModel == TRUE) & (length(saveDRIFTAllModelFit) > 0) ) {
      x1 <- paste0(saveDRIFTAllModelFit[1], " homDRIFTallFit", ".rds"); x1
      x2 <- paste0(workingDirectory)
      CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTallFit, x1, x2)
    }

    if (testDRIFTallModel == TRUE) {
      # Extract estimates & statistics
      driftRows <- seq(1, nlatents^2, 1); driftRows
      homDRIFTall_Drift_Coef <- homDRIFTallFit$output$estimate[driftRows]; homDRIFTall_Drift_Coef
      homDRIFTall_Drift_SE <- homDRIFTallFit$output$standardErrors[driftRows]; homDRIFTall_Drift_SE
      Tvalue <- (homDRIFTall_Drift_Coef/homDRIFTall_Drift_SE); Tvalue
      homDRIFTall_Minus2LogLikelihood  <- homDRIFTallFit$output$Minus2LogLikelihood; homDRIFTall_Minus2LogLikelihood
      homDRIFTall_estimatedParameters  <- length(homDRIFTallFit$output$estimate); homDRIFTall_estimatedParameters
      homDRIFTall_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
        homDRIFTall_estimatedParameters; homDRIFTall_df

      # Compute confidence intervals
      if (confidenceIntervals == TRUE) {
        print(paste0("#################################################################################"))
        print(paste0("######## Computing Confidence Intervals for CoTiMA Model (DRIFTallModel) ########"))
        print(paste0("#################################################################################"))

        homDRIFTallFitCI <- homDRIFTallCI <- list()

        # LOAD
        if (length(loadDRIFTAllModelFit) > 0) {
          homDRIFTallFitCI <- list()
          for (i in 1:(nlatents^2)) {
            x1 <- paste0(workingDirectory, loadDRIFTAllModelFit[1], " homDRIFTallFitCI", i, ".rds"); x1
            homDRIFTallFitCI[[i]] <- readRDS(file=x1)
            homDRIFTallCI[[i]] <- homDRIFTallFitCI[[i]]$output$confidenceIntervals
          }
        }

        # FIT (NOT LOAD)
        if (length(loadDRIFTAllModelFit) < 1) {
          tmpModelMxobjFit <- list()
          for (k in 1:(nlatents^2)) tmpModelMxobjFit[[k]] <- homDRIFTallFit  # copy mxobj part of fitted models multiple times (for each drift coefficient)
          ci <- list()
          for (k in 1:(nlatents^2)) ci[[k]] <- mxCI(driftNames[k]) # make mxCI object for every drift coefficients
          tmpModelMxobj <- list()
          for (k in 1:(nlatents^2)) tmpModelMxobj[[k]] <- mxModel(tmpModelMxobjFit[[k]], ci[[k]]) # make OpenMx Models for all drift coefficients
          mxOption(NULL, 'Number of Threads', 1)
          results <- mclapply(seq(1, (nlatents^2), by=1),
                              function(x) mxRun(tmpModelMxobj[[x]], intervals=TRUE),
                              mc.cores=coresToUse)
          mxOption(key='Number of Threads', value=parallel::detectCores())

          for (k in 1:(nlatents^2)) {
            homDRIFTallFitCI[[k]] <- results[[k]]
            homDRIFTallCI[[k]] <- homDRIFTallFitCI[[k]]$output$confidenceIntervals
          }
        }

        # STORE
        for (k in 1:(nlatents^2)) {
          # fits
          FitList$tmpName <- homDRIFTallFitCI[[k]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("homDRIFTallFitCI", k); tmp
          names(FitList) <- tmp; names(FitList)
          # just confidence intervals
          FitList$tmpName <- homDRIFTallCI[[k]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("homDRIFTallCI", k); tmp
          names(FitList) <- tmp; names(FitList)
        }

        # SAVE
        if (length(saveDRIFTAllModelFit) > 0) {
          for (k in 1:(nlatents^2)) {
            x1 <- paste0(workingDirectory, saveDRIFTAllModelFit[1], " homDRIFTallFitCI", k, ".rds"); x1
            x2 <- paste0(workingDirectory); x2
            CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTallFitCI[[k]], x1, x2)
          }
        }

      } ## END confidence intervals

      # Combine summary information
      homDRIFTallDRIFT_effects <- matrix(t(cbind((homDRIFTall_Drift_Coef), (homDRIFTall_Drift_SE), Tvalue)), 1, 3*length(driftRows), byrow=T)

      # Label summary table
      rownames(homDRIFTallDRIFT_effects) <- c("Fixed Effects")
      newColNames <- c()
      for (j in 1:nlatents) {
        for (h in 1:nlatents) {
          newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)", "Tvalue")
        }
      }
      colnames(homDRIFTallDRIFT_effects) <- newColNames; homDRIFTallDRIFT_effects
      if (confidenceIntervals == TRUE) {
        limits <- matrix(NA, 2, (3*(nlatents^2)))
        tmp <- unlist(homDRIFTallCI); tmp
        tmpSeq1 <- seq(1, (3*(nlatents^2)), 3); tmpSeq1 # drift cols & lower limit cols
        tmpSeq2 <- seq(3, (3*(nlatents^2)), 3); tmpSeq2 # upper limit cols
        limits[1, tmpSeq1] <- tmp[tmpSeq2]
        limits[2, tmpSeq1] <- tmp[tmpSeq1]
        homDRIFTallDRIFT_effects <- rbind(homDRIFTallDRIFT_effects, limits); homDRIFTallDRIFT_effects
        rownames(homDRIFTallDRIFT_effects) <- rbind("Fixed Effects", "upper bound",  "lower bound")
      }
      FitList$homDRIFTallDRIFT_effects <- homDRIFTallDRIFT_effects

      ### Numerically compute Optimal Time lag sensu Dormann & Griffin (2015)
      driftMatrix <- matrix(homDRIFTall_Drift_Coef, nlatents, nlatents); driftMatrix
      OTL <- function(timeRange) {
        expm(driftMatrix * timeRange)[targetRow, targetCol]}
      # loop through all cross effects
      optimalCrossLag <- matrix(NA, nlatents, nlatents)
      maxCrossEffect <- matrix(NA, nlatents, nlatents)
      for (j in 1:nlatents) {
        for (h in 1:nlatents) {
          if (j != h) {
            targetRow <- j
            targetCol <- h
            targetParameters <- sapply(usedTimeRange, OTL)
            maxCrossEffect[j,h] <- max(abs(targetParameters))
            optimalCrossLag[j,h] <- which(abs(targetParameters)==maxCrossEffect[j,h])*stepWidth+0
          }
        }
      }
      #} # END (compute OTL only if CoTiMA results are in a reasonable range )
    }
  } ## END PART 1 CoTiMA (ctsem multigroup): ESTIMATING DRIFTAllModelFit  ##



  ############################################# CoTiMA wit all time lags made equal #####################################
  {
    if (testDRIFTallWOdtModel == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("### Fitting Model without dT (DRIFTallWOdtModel; all dT set to mean(delta_t)) ###"))
      print(paste0("#################################################################################"))

      # Create additional data sets with time lags equalized
      datawide_all_eq_dt <- datawide_all
      targetCols <- grep("dT", colnames(datawide_all_eq_dt)); targetCols
      homDRIFTallWOdtFit_timeLags <- meanDelta
      helper <- datawide_all_eq_dt[,targetCols]
      helper[helper > .001] <- homDRIFTallWOdtFit_timeLags
      datawide_all_eq_dt[,targetCols] <- helper

      # get start values
      homDRIFTallWOdtCTmodelobj <- hetModel
      homDRIFTallWOdtFixedModel <- hetModel
      homDRIFTallWOdtFixedModel$DRIFT <- matrix("groupfixed", nlatents, nlatents)
      homDRIFTallWOdtModModel <- NULL
      homDRIFTallWOdtStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                       ctModelObj = homDRIFTallWOdtCTmodelobj, fixedModel = homDRIFTallWOdtFixedModel,
                                                       modModel = homDRIFTallWOdtModModel, moderatorValues = NULL)
      # LOAD
      if ((testDRIFTallWOdtModel == TRUE) & (length(loadDRIFTallWOdtModelFit) > 0) ) {
        x1 <- paste0(loadDRIFTallWOdtModelFit[1], " homDRIFTallWOdtFit.rds"); x1
        homDRIFTallWOdtFit <- readRDS(file=x1)
      }

      # FIT (NOT LOAD)
      if ((testDRIFTallWOdtModel == TRUE) & (length(loadDRIFTallWOdtModelFit) < 1) ) {
        if (useCTMultigroupFitAlt == FALSE) {
          mxOption(NULL, 'Number of Threads', 1)
          homDRIFTallWOdtStartValues <- changeStartValueLabels(startValues = homDRIFTallWOdtStartValues, ctmodelobj = homDRIFTallWOdtCTmodelobj,
                                                               fixedmodel = homDRIFTallWOdtFixedModel, noOfStudies <- n.studies)
          results <- mclapply(seq(1, refits, by=1),
                              function(refits) ctMultigroupFit(dat=datawide_all_eq_dt, groupings = groupsNamed, retryattempts = retryattempts,
                                                               omxStartValues=homDRIFTallWOdtStartValues,
                                                               ctmodelobj = homDRIFTallWOdtCTmodelobj,
                                                               fixedmodel = homDRIFTallWOdtFixedModel),
                              mc.cores=coresToUse)
          mxOption(key='Number of Threads', value=parallel::detectCores())
          # Select model with best fit
          allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
          homDRIFTallWOdtFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          homDRIFTallWOdtFit <- homDRIFTallWOdtFit$mxobj
        } else {
          mxOption(NULL, 'Number of Threads', 1)
          results <- mclapply(seq(1, refits, by=1),
                              function(refits) ctMultigroupFitAlt(dat=datawide_all_eq_dt, groupings = groups, retryattempts = retryattempts,
                                                                  startValues = homDRIFTallWOdtStartValues,
                                                                  ctmodelobj = homDRIFTallWOdtCTmodelobj,
                                                                  fixedmodel = homDRIFTallWOdtFixedModel),
                              mc.cores=coresToUse)
          mxOption(key='Number of Threads', value=parallel::detectCores())
          # Select model with best fit
          allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
          homDRIFTallWOdtFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
        }
      }

      # STORE
      FitList$DRIFTAllWOdtModelFit <- homDRIFTallWOdtFit

      # SAVE
      if ((testDRIFTallWOdtModel == TRUE) & (length(saveDRIFTallWOdtModelFit) > 0) ) {
        x1 <- paste0(saveDRIFTallWOdtModelFit[1], " homDRIFTallWOdtFit", ".rds"); x1
        x2 <- paste0(workingDirectory); x2
        CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTallWOdtFit, x1, x2)
      }
    }
  } ## END PART 2 CoTiMA (ctsem multigroup): ESTIMATING DRIFTAllModelFit without dt ##

  ################################## FULL CoTiMA (ALL parameters equal across groups) ###################################

  {
    # 5.08.2019 CHD: Test Model with all parameters invariant accross primary studies (on request or if statisticalPower is requested )
    if ((testAllInvariantModel == TRUE) | length(statisticalPower) > 0 ) {

      print(paste0("#################################################################################"))
      print(paste0("## Fitting FULL CoTiMA (AllInvariantModel; ALL parameters equal across groups) ##"))
      print(paste0("#################################################################################"))

      # get start values
      homAllCTmodelobj <- hetModel
      homAllFixedModel <- allFixedModel
      homAllModModel <- NULL
      homAllStartValues <- computeStartValues(singleStudyFits = studyFit,
                                              ctModelObj = homAllCTmodelobj, fixedModel = homAllFixedModel,
                                              modModel = homAllModModel, moderatorValues = NULL)
      # FIT (NOT LOAD)
      mxOption(NULL, 'Number of Threads', 1)
      results <- mclapply(seq(1, refits, by=1),
                          function(x) ctFit(dat=datawide_all, #homModel,
                                            retryattempts=retryattempts,
                                            ctmodelobj = homAllCTmodelobj, #fixedmodel = homAllFixedModel,
                                            omxStartValues=homAllStartValues),
                          mc.cores=coresToUse)
      mxOption(key='Number of Threads', value=parallel::detectCores())
      # Select model with best fit
      allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
      homAllFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
      homAllFit <- homAllFit$mxobj

      # LOAD
      if (length(loadAllInvariantModelFit) > 0) {
        x1 <- paste0(workingDirectory, loadAllInvariantModelFit[1], " homAllFit.rds"); x1
        homAllFit <- readRDS(file=x1)
      }

      # STORE
      FitList$AllInvariantModelFit <- homAllFit

      # SAVE
      if ( ((testAllInvariantModel == TRUE) | length(statisticalPower) > 0 ) & (length(saveAllInvariantModelFit) > 0) ) {
        x1 <- paste0(saveAllInvariantModelFit[1], " homAllFit", ".rds"); x1
        x2 <- paste0(workingDirectory); x2
        CoTiMASaveFile(activateRPB, workingDirectory, homAllFit, x1, x2)
      }

      ### Extract estimates & statistics
      allStandardErrors <- sqrt(2*(diag(solve(homAllFit$output$calculatedHessian)))); allStandardErrors
      ## Drift
      driftRows <- seq(1, nlatents^2, 1); driftRows
      homAll_Drift_Coef <- homAllFit$output$estimate[driftRows]; homAll_Drift_Coef
      homAll_Drift_SE <- homAllFit$output$standardErrors[driftRows]; homAll_Drift_SE
      homAll_Drift_Tvalue <- (homAll_Drift_Coef/homAll_Drift_SE); homAll_Drift_Tvalue
      ## Diffusion (there are 4 'Normal' coefficients, which are 3 (lower triagonal) transformed base coefficients)
      diffusionRows <- (max(driftRows)+1):(max(driftRows)+(nlatents*(nlatents+1)/2)); diffusionRows
      homAll_Diffusion_Coef <- c(homAllFit$output$algebras$ctsem.DIFFUSION); homAll_Diffusion_Coef
      homAll_Diffusion_SE <- allStandardErrors[diffusionRows]; homAll_Diffusion_SE
      homAll_Diffusion_SE_mat <- matrix(NA, nlatents, nlatents)
      counter <- 0
      for (k in 1:nlatents) {
        for (j in k:nlatents){
          counter <- counter + 1
          homAll_Diffusion_SE_mat[k,j] <- homAll_Diffusion_SE[counter]
          if (k != j) homAll_Diffusion_SE_mat[j,k] <- homAll_Diffusion_SE[counter]
        }
      }
      homAll_Diffusion_SE <- c(homAll_Diffusion_SE_mat)
      homAll_Diffusion_Tvalue <- (homAll_Diffusion_Coef/homAll_Diffusion_SE); homAll_Diffusion_Tvalue
      ## T0Var (there are 4 'Normal' coefficients, which are 3 (lower triagonal) transformed base coefficients)
      T0VarRows <- (max(diffusionRows)+1):(max(diffusionRows)+(nlatents*(nlatents+1)/2)); T0VarRows
      # Select the non-transformed ctparameters
      homAll_T0Var_Coef <- c(homAllFit$output$algebras$ctsem.T0VAR); homAll_T0Var_Coef
      homAll_T0Var_SE <- allStandardErrors[T0VarRows]; homAll_T0Var_SE
      homAll_T0Var_SE_mat <- matrix(NA, nlatents, nlatents)
      counter <- 0
      for (k in 1:nlatents) {
        for (j in k:nlatents){
          counter <- counter + 1
          homAll_T0Var_SE_mat[k,j] <- homAll_T0Var_SE[counter]
          if (k != j) homAll_T0Var_SE_mat[j,k] <- homAll_T0Var_SE[counter]
        }
      }
      homAll_T0Var_SE <- c(homAll_T0Var_SE_mat)
      homAll_T0Var_Tvalue <- (homAll_T0Var_Coef/homAll_T0Var_SE); homAll_T0Var_Tvalue
      ## Extract Model Fit
      homAll_Minus2LogLikelihood  <- homAllFit$output$Minus2LogLikelihood; homAll_Minus2LogLikelihood
      homAll_estimatedParameters  <- length(homAllFit$output$estimate); homAll_estimatedParameters
      homAll_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
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
      for (j in 1:nlatents) {
        for (h in 1:nlatents) {
          newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)", "Tvalue")
        }
      }

      colnames(homAll_effects) <- newColNames; homAll_effects
      FitList$homAll_effects <- homAll_effects
    }
  } ## END PART 3 CoTiMA (ctsem multigroup): ESTIMATING FULL CoTiMA (ALL parameters equal across groups) ##



  #######################################################################################################################
  ############################################ STATISTICAL POWER ANALYSES ###############################################
  #######################################################################################################################
  {
    current.start.time <- Sys.time(); current.start.time

    #### Estimate models without single drift effects (required for sample size computations) ###########################
    if (length(statisticalPower) > 0 ) {

      print(paste0("#################################################################################"))
      print(paste0("########################## Statistical Power Analysis ###########################"))
      print(paste0("#################################################################################"))

      # Make series of models. In each one there is one drift effects set to zero
      # computeStartValues is not used because effects are extracted from previous hom model
      homWOCrossModel <- list()
      counter <- 0
      for (j in 1:nlatents) {
        for (h in 1:nlatents) {
          counter <- counter +1
          homWOCrossModel[[counter]] <- hetModel
          homWOCrossModel[[counter]]$DRIFT[h,j] <- 0
        }
      }

      # fit models
      homAllWOSingleDriftFit <- list()
      homAllWOSingleDrift_Diffusion_Coef  <- list()
      homAllWOSingleDrift_Diffusion_SE <- list()
      homAllWOSingleDrift_Diffusion_Tvalue <- list()
      homAllWOSingleDrift_Diffusion_Coef_Normal <- list()
      homAllWOSingleDrift_T0Var_Coef  <- list()
      homAllWOSingleDrift_T0Var_SE <- list()
      homAllWOSingleDrift_T0Var_Tvalue <- list()
      homAllWOSingleDrift_T0Var_Coef_Normal <- list()
      homAllWOSingleDrift_Drift_Coef <- list()
      homAllWOSingleDrift_Drift_SE <- list()
      homAllWOSingleDrift_Drift_Tvalue <- list()

      DRIFT.without.j <- list()             # coefficient in matrix formant with missing coeff set to 0 (used in subsequent power computations)
      DIFFUSION.without.j <- list()         # coefficient in matrix formant (used in subsequent power computations)
      T0VAR.without.j <- list()             # coefficient in matrix formant (used in subsequent power computations)

      driftRows <- 1:(nlatents^2-1); driftRows
      diffusionRows <- (max(driftRows)+1):(max(driftRows)+(nlatents*(nlatents+1)/2)); diffusionRows
      T0VarRows <- (max(diffusionRows)+1):(max(diffusionRows)+(nlatents*(nlatents+1)/2)); T0VarRows

      counter <- 0
      allMinus2LogLikelihood <- homAllWOSingleDriftFit <- list()

      # Sequence of cross effects in drift matrices (only model without cross effects are computed to save time)
      coeffSeq <- seq(1, nlatents^2, 1)[!(seq(1, nlatents^2, 1) %in% seq(1, nlatents^2, (nlatents+1)))]; coeffSeq

      for (j in 1:(nlatents^2)) {
        counter <- counter + 1
        if (j %in% coeffSeq) {
          newStartValues <- omxGetParameters(homAllFit)[-counter]; newStartValues
          homWOCrossModel[[counter]]

          # LOAD
          if (length(loadHomAllWOSingleDriftModelFit) > 0)  {
            x1 <- paste0(loadHomAllWOSingleDriftModelFit[1], " homAllWOSingleDriftFit", j, ".rds"); x1
            homAllWOSingleDriftFit[[j]] <- readRDS(file=x1)
          }

          # FIT (NOT LOAD)
          if (length(loadHomAllWOSingleDriftModelFit) < 1)  {
            mxOption(NULL, 'Number of Threads', 1)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctFit(dat=datawide_all, retryattempts = retryattempts,
                                                       omxStartValues=newStartValues,
                                                       ctmodelobj = homWOCrossModel[[counter]]),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homAllWOSingleDriftFit[[j]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          }

          # STORE
          FitList$tmpName <-homAllWOSingleDriftFit[[j]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("homAllWOSingleDriftFit", j); tmp
          names(FitList) <- tmp; names(FitList)

          # SAVE
          if (length(saveHomAllWOSingleDriftModelFit) > 0)  {
            x1 <- paste0(saveHomAllWOSingleDriftModelFit[1], " homAllWOSingleDriftFit", j, ".rds"); x1
            x2 <- paste0(workingDirectory); x2
            CoTiMASaveFile(activateRPB, workingDirectory, homAllWOSingleDriftFit[[j]], x1, x2)
          }

          # Collect all results (not really necessary)
          homAllWOSingleDrift_Drift_Coef[[j]] <- homAllWOSingleDriftFit[[j]]$mxobj$output$estimate[driftRows]; homAllWOSingleDrift_Drift_Coef[[j]]
          homAllWOSingleDrift_Drift_SE[[j]] <- homAllWOSingleDriftFit[[j]]$mxobj$output$standardErrors[driftRows]; homAllWOSingleDrift_Drift_SE[[j]]
          homAllWOSingleDrift_Drift_Tvalue[[j]] <- (homAllWOSingleDrift_Drift_Coef[[j]]/homAllWOSingleDrift_Drift_SE[[j]]); homAllWOSingleDrift_Drift_Tvalue[[j]]
          homAllWOSingleDrift_Diffusion_Coef[[j]] <- homAllWOSingleDriftFit[[j]]$mxobj$output$estimate[diffusionRows]; homAllWOSingleDrift_Diffusion_Coef[[j]]
          homAllWOSingleDrift_Diffusion_SE[[j]] <- homAllWOSingleDriftFit[[j]]$mxobj$output$standardErrors[diffusionRows]; homAllWOSingleDrift_Diffusion_SE[[j]]
          homAllWOSingleDrift_Diffusion_Tvalue[[j]] <- (homAllWOSingleDrift_Diffusion_Coef[[j]]/homAllWOSingleDrift_Diffusion_SE[[j]]); homAllWOSingleDrift_Diffusion_Tvalue[[j]]
          homAllWOSingleDrift_Diffusion_Coef_Normal[[j]] <- c(homAllWOSingleDriftFit[[j]]$mxobj$algebras$DIFFUSION$result); homAllWOSingleDrift_Diffusion_Coef_Normal[[j]]
          homAllWOSingleDrift_T0Var_Coef[[j]] <- homAllWOSingleDriftFit[[j]]$mxobj$output$estimate[T0VarRows]; homAllWOSingleDrift_T0Var_Coef[[j]]
          homAllWOSingleDrift_T0Var_SE[[j]] <- homAllWOSingleDriftFit[[j]]$mxobj$output$standardErrors[T0VarRows]; homAllWOSingleDrift_T0Var_SE[[j]]
          homAllWOSingleDrift_T0Var_Tvalue[[j]] <- (homAllWOSingleDrift_T0Var_Coef[[j]]/homAllWOSingleDrift_T0Var_SE[[j]]); homAllWOSingleDrift_T0Var_Tvalue[[j]]
          homAllWOSingleDrift_T0Var_Coef_Normal[[j]] <- c(homAllWOSingleDriftFit[[j]]$mxobj$algebras$T0VAR$result); homAllWOSingleDrift_T0Var_Coef_Normal[[j]]
          # Add one zero in each returned vector of drift coefficents to make a full drift matrix later
          DRIFT.without.j[[j]]<-NA
          DIFFUSION.without.j[[j]]<-NA
          T0VAR.without.j[[j]]<-NA
          counter2 <- 0
          for (h in 1:(nlatents^2)) {
            counter2 <- counter2 + 1; counter2
            DRIFT.without.j[[j]][h] <- homAllWOSingleDrift_Drift_Coef[[j]][counter2]
            DIFFUSION.without.j[[j]][h] <- homAllWOSingleDrift_Diffusion_Coef_Normal[[j]][h]
            T0VAR.without.j[[j]][h] <- homAllWOSingleDrift_T0Var_Coef_Normal[[j]][h]
            if (j == h) {
              DRIFT.without.j[[j]][h] <- 0
              counter2 <- counter2 - 1; counter2
            }
          }

          ### Make drift.without.j, diffusion.without.j & T0VAR.without.j matrices,
          ## DRIFT
          DRIFT.without.j[[j]] <- matrix(DRIFT.without.j[[j]], nlatents, nlatents); DRIFT.without.j[[j]]

          # STORE
          FitList$tmpName <- DRIFT.without.j[[j]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("DRIFT.without.j", j); tmp
          names(FitList) <- tmp; names(FitList)

          ## DIFFUSION
          DIFFUSION.without.j[[j]] <- matrix(DIFFUSION.without.j[[j]], nlatents, nlatents); DIFFUSION.without.j[[j]]

          # STORE
          FitList$tmpName <- DIFFUSION.without.j[[j]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("DIFFUSION.without.j", j); tmp
          names(FitList) <- tmp; names(FitList)

          ## T0VAR
          T0VAR.without.j[[j]] <- matrix(T0VAR.without.j[[j]], nlatents, nlatents); T0VAR.without.j[[j]]

          # STORE
          FitList$tmpName <- T0VAR.without.j[[j]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("T0VAR.without.j", j); tmp
          names(FitList) <- tmp; names(FitList)
        }
      }

      DRIFT <- matrix(homAll_Drift_Coef, nlatents, nlatents); DRIFT
      DIFFUSION <- matrix(homAll_Diffusion_Coef, nlatents, nlatents); DIFFUSION
      T0VAR <- matrix(homAll_T0Var_Coef, nlatents, nlatents); T0VAR
    }
  } ## END PART 1 of STATISTICAL POWER ANALYSES

  #### Compute Required Sample Sizes ###################################################################################
  {
    # Fast function to calculate required sample sizes later (as optional replacement for ss.power.reg.coef)
    nestedProbFunT <- function (fvalue, alpha=.05, power=.80, p=2, x) (1-
                                                                         pt(
                                                                           qt((1 - alpha/2), df = (x)-p-1,
                                                                              lower.tail = TRUE, log.p = FALSE),
                                                                           df = (x)-p-1, ncp = sqrt(x) * abs(fvalue),
                                                                           lower.tail = TRUE, log.p = FALSE)) - power

    if (length(statisticalPower) > 0 ) {
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
      plotPairs <- array(dim=c(nlatents^2, length(statisticalPower), length(usedTimeRange), 2))  # all drift effects, all powers, time range, timePoint+SampleSize

      # Plot required sample size for cross effects.
      for (h in 1:length(statisticalPower)) {
        counter <- 0
        for (j1 in 1:(nlatents)) {
          for (j2 in 1:(nlatents)) {
            counter <- counter + 1; counter
            if (j1 != j2) {
              j <- counter
              for (k in 1:(length(usedTimeRange)-1)) {
                delta_t <- usedTimeRange[k+1]; delta_t
                plotPairs[j, h, k, 1] <- usedTimeRange[k+1] # time point

                # R2 in terms of Kelley & Maxwell 2008
                A <- expm(DRIFT %x% delta_t) %*% T0VAR %*% t(expm(DRIFT %x% delta_t)); A         # implied variance at later Tpoint
                S <- matrix(discreteDiffusionFunction(DIFFUSION, DRIFT, delta_t, 1:4), 2,2); S   # residual variance at later Tpoint
                R2 <- A[j2,j2]/( (A + S)[j2,j2] ); R2                                            # explained variance at later Tpoint

                # R2 without j (cross effect) in terms of Kelley & Maxwell 2008
                DRIFT.without.j[[counter]]
                A.j <- expm(DRIFT.without.j[[counter]] %x% delta_t) %*%
                  T0VAR.without.j[[j]] %*%
                  t(expm(DRIFT.without.j[[j]] %x% delta_t)); A.j           # implied variance without j (counter) as predictor at later Tpoint
                S.j <- matrix(
                  discreteDiffusionFunction(DIFFUSION.without.j[[j]],      # implied residual variance without j (counter) as predictor at later Tpoint
                                            DRIFT.without.j[[j]],
                                            delta_t, 1:4), 2,2); S.j
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
                  #plotPairs[j, h, k, 2] <- ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                  #                                           p = nlatents, desired.power = statisticalPower[h],
                  #                                           alpha.level = 0.05)[[1]] #

                  # The following uses our own function
                  signalToNoiseRatios <- sqrt((R2-R2.j)/(1-R2))
                  helper <- round(uniroot.all(nestedProbFunT, c(nlatents+2,999999999),
                                              fvalue=signalToNoiseRatios, alpha=.05,
                                              power=statisticalPower[h], p=nlatents) + .49999, 0)
                  if (length(helper) < 1) helper <- NA
                  plotPairs[j, h, k, 2] <- helper

                  # Post hoc power computations
                  if ( (delta_t %in% tableNxDeltas[ ,-1]) & (h == 1) ){  # do only once (not for all a priori powers)
                    M <- tableNxDeltas[ ,-1] == delta_t; M # temp matrix used below
                    empiricalN <- matrix(tableNxDeltas[apply(M, 1, any), ], ncol=maxTpoints)[,1]; empiricalN
                    empiricalN <- na.omit(empiricalN); empiricalN
                    for (l in empiricalN) {
                      p05 <- ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                               p = nlatents, Specified.N = l, alpha.level = 0.05)[2]
                      p01 <- ss.power.reg.coef(Rho2.Y_X = R2, Rho2.Y_X.without.j = R2.j,
                                               p = nlatents, Specified.N = l, alpha.level = 0.01)[2]
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
      for (j1 in 1:(nlatents)) {
        for (j2 in 1:(nlatents)) {
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
      numberOfEffects <- nlatents^2 - nlatents; numberOfEffects
      tmp <- requiredSampleSizes[[1]]
      for (h in 2:(numberOfEffects)) tmp <- rbind(tmp, requiredSampleSizes[[h]])
      requiredSampleSizes <- t(tmp)
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
      FitList$requiredSampleSizes <- requiredSampleSizes

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
        tmp05 <- listPowerAlpha05[[j]]
        tmp01 <- listPowerAlpha01[[j]]
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
        medianPower0 <- apply(tmp[, targetCols[1:2]], 2, median, na.rm=TRUE); medianPower0
        medianPowerNA <- apply(postHocPower[, targetCols[1:2]], 2, median, na.rm=TRUE); medianPowerNA
        postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
        postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
        postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower
        postHocPower <- rbind(postHocPower, c(c(NA), rep(NA, 3*(maxTpoints-1)))); postHocPower # four times is correct
        postHocPower[dim(postHocPower)[1]-3, targetCols[1:2]] <- meanPower0
        postHocPower[dim(postHocPower)[1]-2, targetCols[1:2]] <- meanPowerNA
        postHocPower[dim(postHocPower)[1]-1, targetCols[1:2]] <- medianPower0
        postHocPower[dim(postHocPower)[1], targetCols[1:2]] <- medianPowerNA
        newNames <- c(paste0("Study_No_", 1:n.studies),
                      "Mean (NA = 0 Power)", "Mean (NA = NA)",
                      "Median (NA = 0 Power)", "Median (NA = NA)")
        rownames(postHocPower) <- newNames
        postHocPowerList[[j]] <- postHocPower

        # STORE
        FitList$tmpName <- postHocPowerList[[j]]
        tmp <- names(FitList); tmp
        tmp[length(tmp)] <- paste0("postHocPowerList", j); tmp
        names(FitList) <- tmp; names(FitList)
      }
      postHocPowerListNames <- currentDriftNames; postHocPowerListNames
    }
  } ### END STATISTICAL POWER ANALYSES ###



  #######################################################################################################################
  ######################### Series of homogeneity models with single drift effects invariant ############################
  #######################################################################################################################
  {
    if (testDRIFTSingleModel == TRUE) {

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
      targetNames <- c()
      for (j in 1:nlatents) {
        for (h in 1:nlatents) {
          targetNames <- c(targetNames, paste0("V",j,"toV", h))
        }
      }

      # get start values
      homDRIFTSingleCTmodelobj <- list()
      homDRIFTSingleFixedModel <- list()
      homDRIFTSingleModModel <- NULL
      homDRIFTSingleStartValues <- list()
      for (k in 1:(nlatents^2)) {
        homDRIFTSingleCTmodelobj[[k]] <- hetModel
        homDRIFTSingleFixedModel[[k]] <- hetModel
        homDRIFTSingleFixedModel[[k]]$DRIFT[k] <- "groupfixed"
        homDRIFTSingleStartValues[[k]] <- computeStartValues(singleStudyFits = studyFit,
                                                             ctModelObj = homDRIFTSingleCTmodelobj[[k]],
                                                             fixedModel = homDRIFTSingleFixedModel[[k]],
                                                             modModel = homDRIFTSingleModModel[[k]], moderatorValues = NULL)
      }

      # Loop through all CROSS coefficients to be tested for invariance
      for (i in 1:(nlatents^2)) {
        coeffNumber <- i
        j <- nlatents^2
        print(paste0("#################################################################################"))
        print(paste0("## Fitting models ", i," of ", j, " with single drift effects invariant (homDRIFTSingleFit) #"))
        print(paste0("#################################################################################"))

        # LOAD
        if (length(loadDRIFTSingleModelFit) > 0 ) {
          x1 <- paste0(workingDirectory, loadDRIFTSingleModelFit, " homDRIFTSingleFits/", loadDRIFTSingleModelFit, " homDRIFTSingleFit", i, ".rds"); x1
          homDRIFTSingleFit[[i]] <- readRDS(file=x1)
          FitList$tmpName <- homDRIFTSingleFit[[i]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("homDRIFTSingleFit", i); tmp
          names(FitList) <- tmp; names(FitList)
        }

        # FIT (NOT LOAD)
        if (length(loadDRIFTSingleModelFit) < 1 ) {

          # FIT (NOT LOAD)
          if (useCTMultigroupFitAlt == FALSE) {
            mxOption(NULL, 'Number of Threads', 1)
            # change names of start values
            homDRIFTSingleStartValues[[i]] <- changeStartValueLabels(startValues = homDRIFTSingleStartValues[[i]],
                                                                     ctmodelobj = homDRIFTSingleCTmodelobj[[i]],
                                                                     fixedmodel = homDRIFTSingleFixedModel[[i]], noOfStudies <- n.studies)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFit(dat=datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                                 omxStartValues = homDRIFTSingleStartValues[[i]],
                                                                 ctmodelobj = homDRIFTSingleCTmodelobj[[i]],
                                                                 fixedmodel = homDRIFTSingleFixedModel[[i]]),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTSingleFit[[i]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
            homDRIFTSingleFit[[i]] <- homDRIFTSingleFit[[i]]$mxobj
          } else {
            mxOption(NULL, 'Number of Threads', 1)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFitAlt(dat=datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                    startValues=homDRIFTSingleStartValues[[i]],
                                                                    ctmodelobj = homDRIFTSingleCTmodelobj[[i]],
                                                                    fixedmodel=homDRIFTSingleFixedModel[[i]]),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTSingleFit[[i]] <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          }

          # STORE
          FitList$tmpName <- homDRIFTSingleFit[[i]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("homDRIFTSingleFit", i); tmp
          names(FitList) <- tmp; names(FitList)

        } # END FIT (NOT LOAD)

        # SAVE
        if (length(saveDRIFTSingleModelFit) > 0)  {
          x1 <- paste0(saveDRIFTSingleModelFit, " homDRIFTSingleFit", i, ".rds"); x1
          x2 <-paste0(saveDRIFTSingleModelFit, " homDRIFTSingleFits/")
          CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTSingleFit[[i]], x1, x2)
        }

        # Extract estimates & statistics
        driftRows <- coeffNumber
        homDRIFTSingle_Drift_Coef[[i]] <- homDRIFTSingleFit[[i]]$output$estimate; homDRIFTSingle_Drift_Coef[[i]]
        homDRIFTSingle_Drift_SE[[i]] <- homDRIFTSingleFit[[i]]$output$standardErrors; homDRIFTSingle_Drift_SE[[i]]
        homDRIFTSingle_Minus2LogLikelihood[[i]] <- homDRIFTSingleFit[[i]]$output$Minus2LogLikelihood; homDRIFTSingle_Minus2LogLikelihood[[i]]
        homDRIFTSingle_estimatedParameters[[i]] <- length(homDRIFTSingleFit[[i]]$output$estimate); homDRIFTSingle_estimatedParameters[[i]]
        homDRIFTSingle_df[[i]] <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
          homDRIFTSingle_estimatedParameters[[i]]; homDRIFTSingle_df[[i]]

        # Combine summary information#
        homDRIFTSingleDRIFT_effects[[i]] <- cbind(homDRIFTSingle_Drift_Coef[[i]][i], homDRIFTSingle_Drift_SE[[i]][i])
        Tvalue <- (homDRIFTSingleDRIFT_effects[[i]][1]/homDRIFTSingleDRIFT_effects[[i]][2]); Tvalue

        # Label summary table
        newColNames <- c()
        for (j in 1:nlatents) {
          for (h in 1:nlatents) {
            newColNames <- c(newColNames, paste0("V",j,"toV", h), "(SE)")
          }
        }
        rownames(homDRIFTSingleDRIFT_effects[[i]]) <- c("Fixed Effects")
        colnames(homDRIFTSingleDRIFT_effects[[i]]) <- newColNames[((driftRows-1)*2+1):((driftRows-1)*2+2)]
        homDRIFTSingleDRIFT_effects[[i]] <- cbind(homDRIFTSingleDRIFT_effects[[i]], Tvalue); homDRIFTSingleDRIFT_effects[[i]]

        # STORE
        FitList$tmpName <- homDRIFTSingleDRIFT_effects[[i]]
        tmp <- names(FitList); tmp
        tmp[length(tmp)] <- paste0("homDRIFTSingleDRIFT_effects", i); tmp
        names(FitList) <- tmp; names(FitList)

        # Compute confidence intervals
        print(paste0("#################################################################################"))
        print(paste0(" Computing Confidence Intervalls for Models with Single Drift Effects Invariant  "))
        print(paste0("#################################################################################"))

        if (confidenceIntervals == TRUE) {
          # LOAD
          if (length(loadDRIFTSingleModelFit) > 0) {
            x1 <- paste0(workingDirectory, loadDRIFTSingleModelFit, " homDRIFTSingleFits/", loadDRIFTSingleModelFit, " homDRIFTSingleFitCI", i, ".rds"); x1
            homDRIFTSingleFitCI[[i]] <- readRDS(file=x1)
            homDRIFTSingleCI[[i]] <- homDRIFTSingleFitCI[[i]]$output$confidenceIntervals
          }

          # FIT (NOT LOAD)
          if (length(loadDRIFTSingleModelFit) < 1) {
            tmpModelMxobjFit <- homDRIFTSingleFit[[i]]
            ci <- mxCI(driftNames[i]); ci
            tmpModelMxobj <- mxModel(tmpModelMxobjFit, ci)
            results <- mxRun(tmpModelMxobj, intervals=TRUE)
            homDRIFTSingleFitCI[[i]] <- results
            homDRIFTSingleCI[[i]] <- homDRIFTSingleFitCI[[i]]$output$confidenceIntervals
          }

          # STORE
          FitList$tmpName <- homDRIFTSingleFitCI[[i]]
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("homDRIFTSingleFitCI", i); tmp
          names(FitList) <- tmp; names(FitList)

          # SAVE
          if (length(saveDRIFTSingleModelFit) > 0)  {
            x1 <- paste0(saveDRIFTSingleModelFit, " homDRIFTSingleFitCI", i, ".rds"); x1
            x2 <-paste0(saveDRIFTSingleModelFit, " homDRIFTSingleFits/")
            CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTSingleFitCI[[i]], x1, x2)
          } ## END save fits

        } # END confidence intervals
      } # END for (i in 1:(nlatents^2))
    } # END testDRIFTSingleModel == TRUE
  } ## END Series of homogeneity models with single drift effects invariant #



  #######################################################################################################################
  ######################################### Test Homogeneity of Cross Effects ###########################################
  #######################################################################################################################
  {
    if (testDRIFTCrossModel == TRUE & nlatents > 1) {

      print(paste0("#################################################################################"))
      print(paste0("####################### Test Homogeneity of Cross Effects #######################"))
      print(paste0("#################################################################################"))

      ## MODEL 2 ##################################################################################
      {
        # LOAD
        if ((testDRIFTCrossModel == TRUE) & (length(loadDRIFTCrossModelFit > 0))) {
          x1 <- paste0(loadDRIFTCrossModelFit[1], " homDRIFTCrossFit", ".rds"); x1
          homDRIFTCross2Fit <- readRDS(file=x1)
        }

        # FIT (NOT LOAD)
        if ((testDRIFTCrossModel == TRUE) & (length(loadDRIFTCrossModelFit) < 1)) {

          homDRIFTCross2CTmodelobj <- hetModel
          homDRIFTCross2FixedModel <- hetModel
          for (i in 1:nlatents) {
            for (j in 1:nlatents) {
              if (i != j) homDRIFTCross2FixedModel$DRIFT[i,j] <- "groupfixed"
            }
          }
          homDRIFTCross2ModModel <- NULL
          homDRIFTCross2StartValues <- computeStartValues(singleStudyFits = studyFit,
                                                          ctModelObj = homDRIFTCross2CTmodelobj,
                                                          fixedModel = homDRIFTCross2FixedModel,
                                                          modModel = homDRIFTCross2ModModel, moderatorValues = NULL)
          # fitting
          print(paste0("#################################################################################"))
          print(paste0("### Fitting Model with Invariant but Unequal Cross Effects (homDRIFTCross2Fit) ##"))
          print(paste0("#################################################################################"))
          if (useCTMultigroupFitAlt == FALSE) {
            mxOption(NULL, 'Number of Threads', 1)
            # change names of start values
            homDRIFTCross2StartValues <- changeStartValueLabels(startValues = homDRIFTCross2StartValues,
                                                                ctmodelobj = homDRIFTCross2CTmodelobj,
                                                                fixedmodel = homDRIFTCross2FixedModel, noOfStudies <- n.studies)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFit(dat= datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                                 #omxStartValues=newStartValues,
                                                                 omxStartValues=homDRIFTCross2StartValues,
                                                                 ctmodelobj = homDRIFTCross2CTmodelobj, fixedmodel = homDRIFTCross2FixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTCross2Fit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
            homDRIFTCross2Fit <- homDRIFTCross2Fit$mxobj
          } else {
            mxOption(NULL, 'Number of Threads', 1)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFitAlt(dat= datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                    startValues = homDRIFTCross2StartValues,
                                                                    ctmodelobj = homDRIFTCross2CTmodelobj,
                                                                    fixedmodel = homDRIFTCross2FixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTCross2Fit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          }
        } # END FIT (NOT LOAD)

        # STORE
        FitList$homDRIFTCrossFit <- homDRIFTCross2Fit

        # SAVE
        if ((testDRIFTCrossModel == TRUE) & (length(saveDRIFTCrossModelFit) > 0)) {
          x1 <- paste0(saveDRIFTCrossModelFit, " homDRIFTCrossFit", ".rds"); x1
          x2 <- paste0(workingDirectory)
          CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTCross2Fit, x1, x2)
        }

        # Extract estimates & statistics
        homDRIFTCross2_Minus2LogLikelihood <- homDRIFTCross2Fit$output$Minus2LogLikelihood
        homDRIFTCross2_estimatedParameters <- length(homDRIFTCross2Fit$output$estimate)
        homDRIFTCross2_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
          homDRIFTCross2_estimatedParameters; homDRIFTCross2_df
      } ## END MODEL 2

      ## MODEL 1 ##################################################################################
      {
        # LOAD
        if ((testDRIFTCrossModel == TRUE) & (length(loadDRIFTCrossModelFit > 0))) {
          x1 <- paste0(loadDRIFTCrossModelFit[1], " homEquDRIFTCrossFit", ".rds"); x1
          homDRIFTCrossFit <- readRDS(file=x1)
        }

        # FIT (NOT LOAD)
        if ((testDRIFTCrossModel == TRUE) & (length(loadDRIFTCrossModelFit) < 1)) {

          homDRIFTCrossCTmodelobj <- hetModel
          for (i in 1:nlatents) {
            for (j in 1:nlatents) {
              if (i != j) homDRIFTCrossCTmodelobj$DRIFT[i,j] <- "V1toV2"
            }
          }
          homDRIFTCrossCTmodelobj$DRIFT[1,2] <- homDRIFTCrossCTmodelobj$DRIFT[2,1]
          homDRIFTCrossFixedModel <- hetModel
          for (i in 1:nlatents) {
            for (j in 1:nlatents) {
              if (i != j) homDRIFTCrossFixedModel$DRIFT[i,j] <- "groupfixed"
            }
          }
          homDRIFTCrossModModel <- NULL
          homDRIFTCrossStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                         ctModelObj = homDRIFTCrossCTmodelobj,
                                                         fixedModel = homDRIFTCrossFixedModel,
                                                         modModel = homDRIFTCrossModModel, moderatorValues = NULL)
          homDRIFTCrossStartValuesTmp <- homDRIFTCrossStartValues # used again further below
          homDRIFTCrossStartValues <- homDRIFTCrossStartValues[unique(names(homDRIFTCrossStartValuesTmp))]

          # fitting
          print(paste0("#################################################################################"))
          print(paste0("### Fitting Model with Invariant & Equal Cross Effects (homEquDRIFTCrossFit) ####"))
          print(paste0("#################################################################################"))
          if (useCTMultigroupFitAlt == FALSE) {
            mxOption(NULL, 'Number of Threads', 1)
            # change names of start values
            homDRIFTCrossStartValues <- changeStartValueLabels(startValues = homDRIFTCrossStartValuesTmp,
                                                               ctmodelobj = homDRIFTCrossCTmodelobj,
                                                               fixedmodel = homDRIFTCrossFixedModel, noOfStudies <- n.studies)
            homDRIFTCrossStartValues <- homDRIFTCrossStartValues[unique(names(homDRIFTCrossStartValues))]
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFit(dat= datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                                 omxStartValues=homDRIFTCrossStartValues,
                                                                 ctmodelobj = homDRIFTCrossCTmodelobj,
                                                                 fixedmodel = homDRIFTCrossFixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select Model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTCrossFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
            homDRIFTCrossFit <- homDRIFTCrossFit$mxobj
          } else {
            mxOption(NULL, 'Number of Threads', 1)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFitAlt(dat= datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                    startValues = homDRIFTCrossStartValues,
                                                                    ctmodelobj = homDRIFTCrossCTmodelobj,
                                                                    fixedmodel = homDRIFTCrossFixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select Model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTCrossFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          }
        } # END FIT (NOT LOAD)

        # STORE
        FitList$homEquDRIFTCrossFit <- homDRIFTCrossFit

        # SAVE
        if ((testDRIFTCrossModel == TRUE) & (length(saveDRIFTCrossModelFit) > 0)) {
          x1 <- paste0(saveDRIFTCrossModelFit, " homEquDRIFTCrossFit", ".rds"); x1
          x2 <- paste0(workingDirectory)
          CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTCrossFit, x1, x2)
        }

        # Extract estimates & statistics
        homDRIFTCross_Minus2LogLikelihood <- homDRIFTCrossFit$output$Minus2LogLikelihood
        homDRIFTCross_estimatedParameters <- length(homDRIFTCrossFit$output$estimate)
        driftRows <- which(names(homDRIFTCrossFit$output$estimate) == "V1toV2")
        #driftRows <- which(table(homDRIFTCrossCTmodelobj$DRIFT) > 1); driftRows
        homDRIFTCross_Drift_Coef <- homDRIFTCrossFit$output$estimate[driftRows]; homDRIFTCross_Drift_Coef
        names(homDRIFTCross_Drift_Coef) <- "AllCrossEffects"
        homDRIFTCross_Drift_SE <- homDRIFTCrossFit$output$standardErrors[driftRows]; homDRIFTCross_Drift_SE
        homDRIFTCross_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
          homDRIFTCross_estimatedParameters; homDRIFTCross_df

        # Compute confidence intervals
        if (confidenceIntervals == TRUE) {

          print(paste0("#################################################################################"))
          print(paste0("###### Computing Confidence Intervalls for Invariant & Equal Cross Effects ######"))
          print(paste0("#################################################################################"))

          # LOAD
          if (length(loadDRIFTCrossModelFit) > 0) {
            x1 <- paste0(loadDRIFTCrossModelFit[1], " homEquDRIFTCrossFitCI", ".rds"); x1
            homDRIFTCrossFitCI <- readRDS(file=x1)
            homDRIFTCrossCI <- homDRIFTCrossFitCI$output$confidenceIntervals
          }

          # FIT (NOT LOAD)
          if (length(loadDRIFTCrossModelFit) < 1) {
            tmpModelMxobjFit <- homDRIFTCrossFit  # copy mxobj part of fitted models
            ci <- mxCI("V1toV2") # make mxCI object for cross effect
            tmpModelMxobj <- mxModel(tmpModelMxobjFit, ci) # make OpenMx Model where constrained parameter matches CI-request
            results <- mxRun(tmpModelMxobj, intervals=TRUE) #mc.cores=coresToUse)
            homDRIFTCrossFitCI <- results
            homDRIFTCrossCI <- homDRIFTCrossFitCI$output$confidenceIntervals
            rownames(homDRIFTCrossCI) <- "AllCrossEffects"; homDRIFTCrossCI
          }

          # STORE
          FitList$homDRIFTCrossFitCI <- homDRIFTCrossFitCI
          FitList$homDRIFTCrossCI <- homDRIFTCrossCI

          # SAVE
          if (length(saveDRIFTCrossModelFit) > 0)  {
            x1 <- paste0(saveDRIFTCrossModelFit, " homEquDRIFTCrossFitCI", ".rds"); x1
            x2 <- paste0(workingDirectory)
            CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTCrossFitCI, x1, x2)
          }
        } # END confidence intervals

        # Combine summary information
        homDRIFTCrossDRIFT_effects <- matrix(t(cbind((homDRIFTCross_Drift_Coef), (homDRIFTCross_Drift_SE)) ), 1, 2*length(min(driftRows)), byrow=T)
        # Label summary table
        rownames(homDRIFTCrossDRIFT_effects) <- c("Fixed Effects")
        colnames(homDRIFTCrossDRIFT_effects) <- c("AllCrossEffects", "(SE)")
        Tvalue <- (homDRIFTCrossDRIFT_effects[1]/homDRIFTCrossDRIFT_effects[2]); Tvalue
        homDRIFTCrossDRIFT_effects <- cbind(homDRIFTCrossDRIFT_effects, Tvalue); homDRIFTCrossDRIFT_effects
        # add confidence intervals
        if (confidenceIntervals == TRUE) {
          homDRIFTCrossDRIFT_effects <- cbind(homDRIFTCrossDRIFT_effects, homDRIFTCrossCI[1], homDRIFTCrossCI[3])
          colnames(homDRIFTCrossDRIFT_effects) <- c(colnames(homDRIFTCrossDRIFT_effects)[1:3], "lower bound", "upper bound")
        }

        # STORE
        FitList$homDRIFTCrossDRIFT_effects <- homDRIFTCrossDRIFT_effects

      } ## END MODEL 1
    }
  } ### END Test Homogeneity of Cross Effects  ###

  #######################################################################################################################
  ######################################### Test Homogeneity of Auto Effects ###########################################
  #######################################################################################################################

  {
    if (testDRIFTAutoModel == TRUE & nlatents > 1) {

      print(paste0("#################################################################################"))
      print(paste0("####################### Test Homogeneity of Auto Effects ########################"))
      print(paste0("#################################################################################"))

      ## MODEL 2 ##################################################################################
      {
        # LOAD
        if (length(loadDRIFTAutoModelFit > 0)) {
          x1 <- paste0(loadDRIFTAutoModelFit[1], " homDRIFTAutoFit", ".rds"); x1
          homDRIFTAuto2Fit <- readRDS(file=x1)
        }

        # FIT (NOT LOAD)
        if (length(loadDRIFTAutoModelFit) < 1) {
          homDRIFTAuto2CTmodelobj <- hetModel
          homDRIFTAuto2FixedModel <- hetModel
          for (i in 1:nlatents) {
            for (j in 1:nlatents) {
              if (i == j) homDRIFTAuto2FixedModel$DRIFT[i,j] <- "groupfixed"
            }
          }
          homDRIFTAuto2ModModel <- NULL
          homDRIFTAuto2StartValues <- computeStartValues(singleStudyFits = studyFit,
                                                         ctModelObj = homDRIFTAuto2CTmodelobj,
                                                         fixedModel = homDRIFTAuto2FixedModel,
                                                         modModel = homDRIFTAuto2ModModel, moderatorValues = NULL)
          # fitting
          print(paste0("#################################################################################"))
          print(paste0("### Fitting Model with Invariant but Unequal Auto Effects (homDRIFTAuto2Fit) ####"))
          print(paste0("#################################################################################"))
          if (useCTMultigroupFitAlt == FALSE) {
            mxOption(NULL, 'Number of Threads', 1)
            homDRIFTAuto2StartValues <- changeStartValueLabels(startValues = homDRIFTAuto2StartValues,
                                                               ctmodelobj = homDRIFTAuto2CTmodelobj,
                                                               fixedmodel = homDRIFTAuto2FixedModel, noOfStudies <- n.studies)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFit(dat= datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                                 omxStartValues = homDRIFTAuto2StartValues,
                                                                 ctmodelobj = homDRIFTAuto2CTmodelobj,
                                                                 fixedmodel = homDRIFTAuto2FixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select Model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTAuto2Fit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
            homDRIFTAuto2Fit <- homDRIFTAuto2Fit$mxobj
          } else {
            mxOption(NULL, 'Number of Threads', 1)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFitAlt(dat= datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                    startValues = homDRIFTAuto2StartValues,
                                                                    ctmodelobj = homDRIFTAuto2CTmodelobj,
                                                                    fixedmodel = homDRIFTAuto2FixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select Model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTAuto2Fit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          }
        } # END FIT (NOT LOAD)

        # STORE
        FitList$homDRIFTAutoFit <- homDRIFTAuto2Fit

        # SAVE
        if ((testDRIFTAutoModel == TRUE) & (length(saveDRIFTAutoModelFit) > 0)) {
          x1 <- paste0(saveDRIFTAutoModelFit, " homDRIFTAutoFit", ".rds"); x1
          x2 <- paste0(workingDirectory)
          CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTAuto2Fit, x1, x2)
        }

        # Extract estimates & statistics
        homDRIFTAuto2_Minus2LogLikelihood <- homDRIFTAuto2Fit$output$Minus2LogLikelihood
        homDRIFTAuto2_estimatedParameters <- length(homDRIFTAuto2Fit$output$estimate)
        homDRIFTAuto2_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
          homDRIFTAuto2_estimatedParameters; homDRIFTAuto2_df

      } ## END MODEL 2

      ## MODEL 1 ##################################################################################
      {
        # LOAD
        if (length(loadDRIFTAutoModelFit > 0)) {
          x1 <- paste0(loadDRIFTAutoModelFit[1], " homEquDRIFTAutoFit", ".rds"); x1
          homDRIFTAutoFit <- readRDS(file=x1)
        }

        if (length(loadDRIFTAutoModelFit) < 1) {
          homDRIFTAutoCTmodelobj <- hetModel

          for (i in 1:nlatents) {
            for (j in 1:nlatents) {
              if (i == j) homDRIFTAutoCTmodelobj$DRIFT[i,j] <- "V1toV1"
            }
          }
          homDRIFTAutoFixedModel <- hetModel
          for (i in 1:nlatents) {
            for (j in 1:nlatents) {
              if (i == j) homDRIFTAutoFixedModel$DRIFT[i,j] <- "groupfixed"
            }
          }
          homDRIFTAutoModModel <- NULL
          homDRIFTAutoStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                        ctModelObj = homDRIFTAutoCTmodelobj,
                                                        fixedModel = homDRIFTAutoFixedModel,
                                                        modModel = homDRIFTAutoModModel, moderatorValues = NULL)
          homDRIFTAutoStartValuesTmp <- homDRIFTAutoStartValues
          homDRIFTAutoStartValues <- homDRIFTAutoStartValues[unique(names(homDRIFTAutoStartValuesTmp))]

          # fitting
          print(paste0("#################################################################################"))
          print(paste0("#### Fitting Model with Invariant & Equal Auto Effects (homEquDRIFTAutoFit) #####"))
          print(paste0("#################################################################################"))
          if (useCTMultigroupFitAlt == FALSE) {
            mxOption(NULL, 'Number of Threads', 1)
            homDRIFTAutoStartValues <- changeStartValueLabels(startValues = homDRIFTAutoStartValuesTmp,
                                                              ctmodelobj = homDRIFTAutoCTmodelobj,
                                                              fixedmodel = homDRIFTAutoFixedModel, noOfStudies <- n.studies)
            homDRIFTAutoStartValues <- homDRIFTAutoStartValues[unique(names(homDRIFTAutoStartValues))]
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFit(dat= datawide_all, groupings = groupsNamed, retryattempts = retryattempts,
                                                                 omxStartValues = homDRIFTAutoStartValues,
                                                                 ctmodelobj = homDRIFTAutoCTmodelobj,
                                                                 fixedmodel = homDRIFTAutoFixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$mxobj$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTAutoFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
            homDRIFTAutoFit <- homDRIFTAutoFit$mxobj
          } else {
            mxOption(NULL, 'Number of Threads', 1)
            results <- mclapply(seq(1, refits, by=1),
                                function(refits) ctMultigroupFitAlt(dat= datawide_all, groupings = groups, retryattempts = retryattempts,
                                                                    startValues = homDRIFTAutoStartValues,
                                                                    ctmodelobj = homDRIFTAutoCTmodelobj,
                                                                    fixedmodel = homDRIFTAutoFixedModel),
                                mc.cores=coresToUse)
            mxOption(key='Number of Threads', value=parallel::detectCores())
            # Select model with best fit
            allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
            homDRIFTAutoFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          }

          # STORE
          FitList$homEquDRIFTAutoFit <- homDRIFTAutoFit

          # SAVE
          if (length(saveDRIFTAutoModelFit) > 0) {
            x1 <- paste0(saveDRIFTAutoModelFit, " homEquDRIFTAutoFit", ".rds"); x1
            x2 <- paste0(workingDirectory)
            CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTAutoFit, x1, x2)
          }
        }

        # Extract estimates & statistics
        homDRIFTAuto_Minus2LogLikelihood <- homDRIFTAutoFit$output$Minus2LogLikelihood
        homDRIFTAuto_estimatedParameters <- length(homDRIFTAutoFit$output$estimate)
        driftRows <- which(names(homDRIFTAutoFit$output$estimate) == "V1toV1")
        homDRIFTAuto_Drift_Coef <- homDRIFTAutoFit$output$estimate[driftRows]; homDRIFTAuto_Drift_Coef
        names(homDRIFTAuto_Drift_Coef) <- "AllAutoEffects"
        homDRIFTAuto_Drift_SE <- homDRIFTAutoFit$output$standardErrors[driftRows]; homDRIFTAuto_Drift_SE
        homDRIFTAuto_df <- ((nlatents * unlist(allTpoints)) %*% ((nlatents * unlist(allTpoints)) +1 )) / 2 -
          homDRIFTAuto_estimatedParameters; homDRIFTAuto_df

        # Compute confidence intervals
        if (confidenceIntervals == TRUE) {
          print(paste0("#################################################################################"))
          print(paste0("###### Computing Confidence Intervalls for Invariant & Equal Auto Effects #######"))
          print(paste0("#################################################################################"))

          # LOAD
          if (length(loadDRIFTAutoModelFit) > 0) {
            x1 <- paste0(loadDRIFTAutoModelFit[1], " homEquDRIFTAutoFitCI", ".rds"); x1
            homDRIFTAutoFitCI <- readRDS(file=x1)
            homDRIFTAutoCI <- homDRIFTAutoFitCI$output$confidenceIntervals
          }

          # FIT (NOT LOAD)
          if (length(loadDRIFTAutoModelFit) < 1) {
            tmpModelMxobjFit <- homDRIFTAutoFit  # copy mxobj part of fitted models
            ci <- mxCI("V1toV1") # make mxCI object for cross effect
            tmpModelMxobj <- mxModel(tmpModelMxobjFit, ci) # make OpenMx Model where constrained parameter matches CI-request
            results <- mxRun(tmpModelMxobj, intervals=TRUE) #mc.cores=coresToUse)
            homDRIFTAutoFitCI <-  results
            homDRIFTAutoCI <- homDRIFTAutoFitCI$output$confidenceIntervals
            rownames(homDRIFTAutoCI) <- "AllAutoEffects"
          }

          # STORE
          FitList$homDRIFTAutoFitCI <- homDRIFTAutoFitCI
          FitList$homDRIFTAutoCI <- homDRIFTAutoCI

          # SAVE
          if (length(saveDRIFTAutoModelFit) > 0)  {
            x1 <- paste0(saveDRIFTAutoModelFit, " homEquDRIFTAutoFitCI", ".rds"); x1
            x2 <- paste0(workingDirectory)
            CoTiMASaveFile(activateRPB, workingDirectory, homDRIFTAutoFitCI, x1, x2)
          }
        } # END confidence intervals

        # Combine summary information
        homDRIFTAutoDRIFT_effects <- matrix(t(cbind((homDRIFTAuto_Drift_Coef), (homDRIFTAuto_Drift_SE)) ), 1, 2*length(min(driftRows)), byrow=T)
        # Label summary table
        rownames(homDRIFTAutoDRIFT_effects) <- c("Fixed Effects")
        colnames(homDRIFTAutoDRIFT_effects) <- c("AutoEffects", "(SE)")
        Tvalue <- (homDRIFTAutoDRIFT_effects[1]/homDRIFTAutoDRIFT_effects[2]); Tvalue
        homDRIFTAutoDRIFT_effects <- cbind(homDRIFTAutoDRIFT_effects, Tvalue); homDRIFTAutoDRIFT_effects
        # Add confidence intervals
        if (confidenceIntervals == TRUE) {
          homDRIFTAutoDRIFT_effects <- cbind(homDRIFTAutoDRIFT_effects, homDRIFTAutoCI[1], homDRIFTAutoCI[3])
          colnames(homDRIFTAutoDRIFT_effects) <- c(colnames(homDRIFTAutoDRIFT_effects)[1:3], "lower bound", "upper bound")
        }

        # STORE
        FitList$homDRIFTAutoDRIFT_effects <- homDRIFTAutoDRIFT_effects

      } ## END MODEL 1
    }
  } ### END Test Homogeneity of Auto Effects  ###



  #######################################################################################################################
  ############################################ Modertor Analyses ########################################################
  #######################################################################################################################
  {
    if (testModeratorModel == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("############################ Fitting Moderator Model ############################"))
      print(paste0("#################################################################################"))

      # LOAD
      if (length(loadModeratorModelFit) > 0) {
        x1 <- paste0(loadModeratorModelFit[1], " modModel1mxFit_modNo", moderatorNumber ,".rds"); x1
        x2 <- paste0(workingDirectory)
        modModel1mxFit <- readRDS(file=paste0(x2, x1))
        oldDrift <- studyFit[[1]]$mxobj$matrices$DRIFT; oldDrift # should always exists
      }

      if (length(loadModeratorModelFit) < 1) {
        # FIT
        modModel1mxCTmodelobj <- hetModel
        modModel1mxFixedModel <- hetModel
        modModel1mxModModel <-  hetModel
        modModel1mxModModel$DRIFT <- matrix("moderated", nlatents, nlatents)
        moderatorValues <- unlist(lapply(allModerators, function(extract) extract[moderatorNumber])); moderatorValues
        modModel1mxStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                     ctModelObj = modModel1mxCTmodelobj,
                                                     fixedModel = modModel1mxFixedModel,
                                                     modModel = modModel1mxModModel,
                                                     moderatorValues = moderatorValues)
        modModel1mxStartValues
        mxOption(NULL, 'Number of Threads', 1)
        results <- mclapply(seq(1, refits, by=1),
                            function(refits) ctMultigroupFitAlt(dat=datawide_all,
                                                                groupings = groups,
                                                                startValues = modModel1mxStartValues,
                                                                ctmodelobj = modModel1mxCTmodelobj,
                                                                fixedmodel = modModel1mxFixedModel,
                                                                moderators = moderatorGroups,
                                                                modmodel = modModel1mxModModel,
                                                                extraTries=40),
                            mc.cores=coresToUse)
        mxOption(key='Number of Threads', value=parallel::detectCores())
        # Select model with best fit
        allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
        modModel1mxFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
      } # END if (length(loadModeratorModelFit) < 1 #

      # STORE
      FitList$modModel1mxFit <- modModel1mxFit

      # SAVE
      if (length(saveModeratorModelFit) > 0) {
        x1 <- paste0(saveModeratorModelFit[1], " modModel1mxFit_modNo", moderatorNumber ,".rds"); x1
        x2 <- paste0(workingDirectory)
        CoTiMASaveFile(activateRPB, workingDirectory, modModel1mxFit, x1, x2)
      }

      # extract Results
      modModel1mx_Minus2LogLikelihood <- modModel1mxFit$output$Minus2LogLikelihood; modModel1mx_Minus2LogLikelihood
      modModel1mx_estimatedParameters <- length(modModel1mxFit$output$estimate); modModel1mx_estimatedParameters
      targetDriftNames <- c()
      for (i in unique(moderatorValues)) {
        newLabels <- paste0(driftNames, "_M", i)         # change labels of drift elements (adding "_" and value of moderator)
        targetDriftNames <- c(targetDriftNames, newLabels); targetDriftNames
      }
      modModel1mx_Drift_Coef <- modModel1mxFit$output$estimate[c(targetDriftNames)]; modModel1mx_Drift_Coef
      modModel1mx_Drift_SE <- modModel1mxFit$output$standardErrors[c(targetDriftNames),]; modModel1mx_Drift_SE
      modModel1mx_Tvalue <- modModel1mx_Drift_Coef/modModel1mx_Drift_SE; modModel1mx_Tvalue

      modModel1DRIFT_effects <- matrix(t(cbind(modModel1mx_Drift_Coef,
                                               modModel1mx_Drift_SE,
                                               modModel1mx_Tvalue)), nrow=3)
      colnames(modModel1DRIFT_effects) <- targetDriftNames; modModel1DRIFT_effects
      rownames(modModel1DRIFT_effects) <- c( "all Grps.", "(SE)", "Tvalue"); modModel1DRIFT_effects

      newModModel1DRIFT_effects <- modModel1DRIFT_effects[, 1:(nlatents^2)]; newModModel1DRIFT_effects
      for (j in 2:numberOfModerators) {
        newModModel1DRIFT_effects <- rbind(newModModel1DRIFT_effects,
                                           modModel1DRIFT_effects[,(1+(j-1)*nlatents^2):((j)*nlatents^2) ])
      }
      targetNames <- c()
      for (j in unique(unlist(allModerators))) targetNames <- c(targetNames, paste0("Mod. Grp. ", j), "(SE)", "Tvalue")
      rownames(newModModel1DRIFT_effects) <- targetNames;
      colnames(newModModel1DRIFT_effects) <- c(studyFit[[1]]$mxobj$DRIFT$labels); newModModel1DRIFT_effects

      #STORE
      FitList$modModel1DRIFT_effects <- newModModel1DRIFT_effects

      ############### moderator analyses of single UNCONSTRAINED drift coefficients ############################

      print(paste0("#################################################################################"))
      print(paste0("################ Moderator Analyses of Single Drift Coefficients ################"))
      print(paste0("#################################################################################"))

      # LOAD
      if (length(loadModeratorModelFit) > 0) {
        for (i in 1:(nlatents^2)) {
          x1 <- paste0(loadModeratorModelFit[1], " modModel2mxFit", i, "_modNo", moderatorNumber ,".rds"); x1
          x2 <- paste0(workingDirectory)
          modModel2mxFit[[i]] <- readRDS(file=paste0(x2, x1))
        }
      }

      # CREATE MODELS & FIT (NOT LOAD)
      if (length(loadModeratorModelFit) < 1) {
        modModel2mxFit <- list()
        for (j in 1:(nlatents^2)) {
          targetCol <- seq(1, (nlatents^2))[-j]; targetCol
          modModel2mxCTmodelobj <- hetModel
          modModel2mxFixedModel <- hetModel
          modModel2mxFixedModel$DRIFT[targetCol] <- "groupfixed"
          modModel2mxModModel <-  hetModel
          modModel2mxModModel$DRIFT[j] <- "moderated"
          moderatorValues <- unlist(lapply(allModerators, function(extract) extract[moderatorNumber])); moderatorValues
          modModel2mxStartValues <- computeStartValues(singleStudyFits = studyFit,
                                                       ctModelObj = modModel2mxCTmodelobj,
                                                       fixedModel = modModel2mxFixedModel,
                                                       modModel = modModel2mxModModel, moderatorValues = moderatorValues)

          mxOption(NULL, 'Number of Threads', 1)
          results <- mclapply(seq(1, refits, by=1),
                              function(refits) ctMultigroupFitAlt(dat=datawide_all,
                                                                  groupings = groups, retryattempts = retryattempts,
                                                                  startValues = modModel2mxStartValues,
                                                                  ctmodelobj = modModel2mxCTmodelobj,
                                                                  fixedmodel = modModel2mxFixedModel,
                                                                  moderators = moderatorGroups,
                                                                  modmodel = modModel2mxModModel,
                                                                  extraTries=40),
                              mc.cores=coresToUse)
          mxOption(key='Number of Threads', value=parallel::detectCores())
          # Select model with best fit
          allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
          currentModModel2mxFit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
          FitList$tmpName <- currentModModel2mxFit
          tmp <- names(FitList); tmp
          tmp[length(tmp)] <- paste0("ModModel2mxFit", i); tmp
          names(FitList) <- tmp; names(FitList)
          modModel2mxFit[[j]] <- currentModModel2mxFit
        }
      }

      # SAVE
      if (length(saveModeratorModelFit) > 0) {
        for (i in 1:(nlatents^2)) {
          x1 <- paste0(saveModeratorModelFit[1], " modModel2mxFit", i, "_modNo", moderatorNumber ,".rds"); x1
          x2 <- paste0(workingDirectory)
          CoTiMASaveFile(activateRPB, workingDirectory, modModel2mxFit[[i]], x1, x2)
        }
      }

      # Required for parameter extraction
      oldLabels <- c(studyFit[[1]]$mxobj$DRIFT$labels); oldLabels
      allLabels <- list()
      for (j in unique(unlist(allModerators))) allLabels[[j]] <- c(paste0(oldLabels, "_M", j))
      x1 <- unlist(allLabels)
      allLabels <- c(oldLabels, x1); allLabels

      # Extract model fits
      modModel2mx_Minus2LogLikelihood <-  modModel2mx_estimatedParameters <- modModel2mx_Drift_Coef <- list()
      modModel2mx_Drift_SE <- modModel2mx_Tvalue <- list()
      for (i in 1:(nlatents^2)) {
        modModel2mx_Minus2LogLikelihood[[i]] <- modModel2mxFit[[i]]$output$Minus2LogLikelihood
        modModel2mx_estimatedParameters[[i]] <- length(modModel2mxFit[[i]]$output$estimate)
        modModel2mx_Drift_Coef[[i]] <- modModel2mxFit[[i]]$output$estimate[c(names(modModel2mxFit[[i]]$output$estimate) %in% allLabels)]
        modModel2mx_Drift_SE[[i]] <- modModel2mxFit[[i]]$output$standardErrors[c(rownames(modModel2mxFit[[i]]$output$standardErrors) %in% allLabels)]
        names(modModel2mx_Drift_SE[[i]]) <- names(modModel2mx_Drift_Coef[[i]])
        modModel2mx_Tvalue[[i]] <- modModel2mx_Drift_Coef[[i]]/modModel2mx_Drift_SE[[i]]
        names(modModel2mx_Tvalue[[i]]) <- names(modModel2mx_Drift_Coef[[i]])
      }

      # Summarize results
      modModel2mx_all_Drift_Coef <- matrix(NA, (nlatents^2), length(allLabels))
      modModel2mx_all_Drift_SE <- matrix(NA, (nlatents^2), length(allLabels))
      modModel2mx_all_Tvalue <- matrix(NA, (nlatents^2), length(allLabels))
      colnames(modModel2mx_all_Drift_Coef) <- allLabels
      colnames(modModel2mx_all_Drift_SE) <- allLabels
      colnames(modModel2mx_all_Tvalue) <- allLabels
      for (i in 1:(nlatents^2)) {
        modModel2mx_all_Drift_Coef[i , names(modModel2mx_Drift_Coef[[i]])] <- modModel2mx_Drift_Coef[[i]]
        modModel2mx_all_Drift_SE[i , names(modModel2mx_Drift_SE[[i]])] <- modModel2mx_Drift_SE[[i]]
        modModel2mx_all_Tvalue[i , names(modModel2mx_Tvalue[[i]])] <- modModel2mx_Tvalue[[i]]
      }
      selectCols <- c(rbind(
        matrix(1:ncol(modModel2mx_all_Drift_Coef),ncol=ncol(modModel2mx_all_Drift_Coef)),
        matrix(1:ncol(modModel2mx_all_Drift_SE),ncol=ncol(modModel2mx_all_Drift_SE)) +
          ncol(modModel2mx_all_Drift_Coef),
        matrix(1:ncol(modModel2mx_all_Tvalue),ncol=ncol(modModel2mx_all_Tvalue)) +
          ncol(modModel2mx_all_Drift_Coef) + ncol(modModel2mx_all_Drift_SE)
      )); selectCols
      modModel2mxAllDRIFT_effects <- cbind(modModel2mx_all_Drift_Coef,
                                           modModel2mx_all_Drift_SE,
                                           modModel2mx_all_Tvalue)[ , selectCols]
      modModel2mxAllDRIFT_effects
      newColNames <- colnames(modModel2mxAllDRIFT_effects); newColNames
      newColNames[seq(2, length(newColNames), 3)] <- "(SE)"; newColNames
      newColNames[seq(3, length(newColNames), 3)] <- "Tvalue"; newColNames
      colnames(modModel2mxAllDRIFT_effects) <- newColNames
      modModel2mxAllDRIFT_effects <- round(t(modModel2mxAllDRIFT_effects), digits)
      colnames(modModel2mxAllDRIFT_effects) <- driftNames

      newResults1 <- (modModel2mxAllDRIFT_effects[-(1:(3*(nlatents^2))),]); newResults1
      for (i in unique((unlist(allModerators)))){
        targetRowNames <- grep(paste0("_",i), rownames(newResults1)); targetRowNames
        rownames(newResults1)[targetRowNames] <- c(rep(paste0("Moderator=", i), nlatents^2))
      }
      modModel2mxAllDRIFT_effects <- newResults1
      FitList$modModel2mxAllDRIFT_effects <- modModel2mxAllDRIFT_effects
    }
  } ### END Modertor Analyses ###



  #######################################################################################################################
  ###################################### Analysis Of User-Specified Model ###############################################
  #######################################################################################################################
  {

      print(paste0("#################################################################################"))
      print(paste0("######################### Fitting User-Specified Model  #########################"))
      print(paste0("#################################################################################"))

      #backup <- testUserSpecifiedModel$singleStudyModelFits; backup
      #backup2 <- studyFit

      if(!(is.null(testUserSpecifiedModel$listOfFixedModelMatrices$DRIFT)
           & is.null(testUserSpecifiedModel$listOfFixedModelMatrices$DIFFUSION)
           & is.null(testUserSpecifiedModel$listOfFixedModelMatrices$T0VAR)
           & is.null(testUserSpecifiedModel$listOfModeratorModelMatrices$DRIFT)
           & is.null(testUserSpecifiedModel$listOfModeratorModelMatrices$DIFFUSION)
           & is.null(testUserSpecifiedModel$listOfModeratorModelMatrices$T0VAR)
           & is.null(testUserSpecifiedModel$moderatorValues))) { # anything specified?

        # FIT (NOT LOAD)
        if (length(loadUserSpecifiedModel) < 1) {

          if (is.null(testUserSpecifiedModel$singleStudyModelFits)) {
            testUserSpecifiedModel$singleStudyModelFits <- studyFit
          } else  {
            tmp <- list()
            for (i in testUserSpecifiedModel$singleStudyModelFit) {
              tmp[[i]] <- studyFit[i]
            }
            testUserSpecifiedModel$singleStudyModelFits <- tmp
          }

          testUserSpecifiedModel

          newList <- list()
          for (i in 1:length(testUserSpecifiedModel$singleStudyModelFits)) newList[[i]] <- testUserSpecifiedModel$singleStudyModelFits[[i]][[1]]

          UserSpecifiedModelFit <- fitUserSpecifiedModel(singleStudyModelFits = newList,
                                                         moderatorValues = testUserSpecifiedModel$moderatorValues ,
                                                         listOfFixedModelMatrices = testUserSpecifiedModel$listOfFixedModelMatrices,
                                                         listOfModeratorModelMatrices = testUserSpecifiedModel$listOfModeratorModelMatrices,
                                                         coresToUse=coresToUse,
                                                         refits=refits,
                                                         confidenceIntervals=confidenceIntervals,
                                                         digits=digits)
          FitList$UserSpecifiedModelFit <- UserSpecifiedModelFit
        } # END FIT (NOT LOAD)

        # SAVE
        if (length(saveUserSpecifiedModel) > 0)  {
          x1 <- paste0(saveUserSpecifiedModel, "userSpecifiedModelFit.rds"); x1
          x2 <- paste0(workingDirectory)
          CoTiMASaveFile(activateRPB, workingDirectory, UserSpecifiedModelFit, x1, x2)
        }

      } # END anything specified?

      # LOAD
      if (length(loadUserSpecifiedModel) > 0)  {
        x1 <- paste0(loadUserSpecifiedModel, "userSpecifiedModelFit.rds"); x1
        UserSpecifiedModelFit <- readRDS(file=x1)
        FitList$UserSpecifiedModelFit <- UserSpecifiedModelFit
      }
    }


  #######################################################################################################################
  ################## Model Comparisons (-2LL Difference Test ~ Chi-square Difference Tests) #############################
  #######################################################################################################################

  {
    print(paste0("#################################################################################"))
    print(paste0("#################### Model Comparisons between Nested Models ####################"))
    print(paste0("#################################################################################"))

    if (testDRIFTallModel == TRUE) {
      allStudies_homDRIFTall_Minus2LogLikelihood <- homDRIFTall_Minus2LogLikelihood -
        allStudies_Minus2LogLikelihood; allStudies_homDRIFTall_Minus2LogLikelihood
      allStudies_homDRIFTall_df <- homDRIFTall_df -
        allStudies_df; allStudies_homDRIFTall_df
      allStudies_homDRIFTall_prob <- 1-pchisq(abs(allStudies_homDRIFTall_Minus2LogLikelihood),
                                              allStudies_homDRIFTall_df); allStudies_homDRIFTall_prob
    }

    if (testDRIFTSingleModel == TRUE) {
      allStudies_homDRIFTSingle_Minus2LogLikelihood <- allStudies_homDRIFTSingle_df <- allStudies_homDRIFTSingle_prob <- list()

      for (i in 1:(nlatents^2)) {
        allStudies_homDRIFTSingle_Minus2LogLikelihood[[i]] <- homDRIFTSingle_Minus2LogLikelihood[[i]] -
          allStudies_Minus2LogLikelihood; allStudies_homDRIFTSingle_Minus2LogLikelihood[[i]]
        allStudies_homDRIFTSingle_df[[i]] <- homDRIFTSingle_df[[i]] -
          allStudies_df; allStudies_homDRIFTSingle_df[[i]]
        allStudies_homDRIFTSingle_prob[[i]] <- 1-pchisq(abs(allStudies_homDRIFTSingle_Minus2LogLikelihood[[i]]),
                                                        allStudies_homDRIFTSingle_df[[i]]); allStudies_homDRIFTSingle_prob[[i]]
      }
    }

    if ( (testDRIFTallWOdtModel == TRUE) &  (testDRIFTallModel == TRUE) ) {
      homDRIFTall_DRIFTallWOdt_Minus2LogLikelihood <- homDRIFTall_Minus2LogLikelihood -
        homDRIFTallWOdtFit$mxobj$output$Minus2LogLikelihood; homDRIFTall_DRIFTallWOdt_Minus2LogLikelihood
    }

    if (testDRIFTCrossModel == TRUE & nlatents > 1) {
      homDRIFTCross_homDRIFTCross2_Minus2LogLikelihood <- homDRIFTCross2_Minus2LogLikelihood -
        homDRIFTCross_Minus2LogLikelihood; homDRIFTCross_homDRIFTCross2_Minus2LogLikelihood
      homDRIFTCross_homDRIFTCross2_df <- homDRIFTCross_df -
        homDRIFTCross2_df; homDRIFTCross_homDRIFTCross2_df
      homDRIFTCross_homDRIFTCross2_prob <- 1-pchisq(abs(homDRIFTCross_homDRIFTCross2_Minus2LogLikelihood),
                                                    homDRIFTCross_homDRIFTCross2_df); homDRIFTCross_homDRIFTCross2_prob
    }

    if (testDRIFTAutoModel == TRUE  & nlatents > 1) {
      homDRIFTAuto_homDRIFTAuto2_Minus2LogLikelihood <- homDRIFTAuto2_Minus2LogLikelihood -
        homDRIFTAuto_Minus2LogLikelihood; homDRIFTAuto_homDRIFTAuto2_Minus2LogLikelihood
      homDRIFTAuto_homDRIFTAuto2_df <- homDRIFTAuto2_df -
        homDRIFTAuto_df; homDRIFTAuto_homDRIFTAuto2_df
      homDRIFTAuto_homDRIFTAuto2_prob <- 1-pchisq(abs(homDRIFTAuto_homDRIFTAuto2_Minus2LogLikelihood),
                                                  abs(homDRIFTAuto_homDRIFTAuto2_df)); homDRIFTAuto_homDRIFTAuto2_prob
    }

    if ( (testModeratorModel == TRUE) & (testDRIFTallModel == TRUE) ) {
      homDRIFTall_modModel1_Minus2LogLikelihood <- homDRIFTall_Minus2LogLikelihood -
        modModel1mx_Minus2LogLikelihood; homDRIFTall_modModel1_Minus2LogLikelihood
      homDRIFTall_modModel1_df <- modModel1mx_estimatedParameters -
        homDRIFTall_estimatedParameters; homDRIFTall_modModel1_df
      homDRIFTall_modModel1_prob <- 1-pchisq(abs(homDRIFTall_modModel1_Minus2LogLikelihood),
                                             homDRIFTall_modModel1_df); homDRIFTall_modModel1_prob

      homDRIFTall_modModel2_Minus2LogLikelihood <- homDRIFTall_modModel2_df <- homDRIFTall_modModel2_prob <- list()

      for (i in 1:(nlatents^2)) {
        homDRIFTall_modModel2_Minus2LogLikelihood[[i]] <- homDRIFTall_Minus2LogLikelihood -
          modModel2mx_Minus2LogLikelihood[[i]]; homDRIFTall_modModel2_Minus2LogLikelihood[[i]]
        homDRIFTall_modModel2_df[[i]] <- modModel2mx_estimatedParameters[[i]] -
          homDRIFTall_estimatedParameters ; homDRIFTall_modModel2_df[[i]]
        homDRIFTall_modModel2_prob[[i]] <- 1-pchisq(abs(homDRIFTall_modModel2_Minus2LogLikelihood[[i]]),
                                                    homDRIFTall_modModel2_df[[i]]); homDRIFTall_modModel2_prob[[i]]
      }
    }

  } ### END Model Comparisons (-2LL Difference Test ~ Chi-square Difference Tests) ###

  #######################################################################################################################
  ##################################### Analyses of Publication Bias ####################################################
  #######################################################################################################################

  if (publicationBias == TRUE) {

    print(paste0("#################################################################################"))
    print(paste0("######################### Analysis of Publication Bias ##########################"))
    print(paste0("#################################################################################"))

    targetNames <- c()
    for (j in 1:nlatents) {
      for (h in 1:nlatents) {
        targetNames <- c(targetNames, paste0("V",j,"toV", h))
      }
    }
    colnames(DRIFTCoeff) <- targetNames; DRIFTCoeff
    targetNames <- c()
    for (j in 1:nlatents) {
      for (h in 1:nlatents) {
        targetNames <- c(targetNames, paste0("SE(V",j,"toV", h, ")"))
      }
    }
    colnames(DRIFTSE) <- targetNames; DRIFTSE

    DRIFTCoeffSND <- DRIFTCoeff / DRIFTSE; DRIFTCoeffSND
    DRIFTPrecision <- c(rep(1, nlatents^2))/(DRIFTSE); DRIFTPrecision
    colnames(DRIFTPrecision) <- colnames(DRIFTCoeffSND); DRIFTPrecision

    message1 <- "The pos. & sign. intercept indicates that SMALLER studies produced more positive (or less negative) effects"
    message2 <- "The neg. & sign. intercept indicates that LARGER studies produced more positive (or less negative) effects"
    tmp <- c()

    eggerDrift <- list()
    for (j in 1:(nlatents^2)) {
      eggerDrift[[j]] <- lm(DRIFTCoeffSND[,j]~DRIFTPrecision[,j]) # This is identical to a weighted regression of drift on se ...
      if (summary(eggerDrift[[j]])$coefficients[1,1] > 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message1
      }
      if (summary(eggerDrift[[j]])$coefficients[1,1] < 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message2
      }
    }

    FREAResults <- list()

    FREAResults[[1]] <- "############# Eggers Test for DRIFT Parameter Estimates  ###############################"
    FREACounter <- 1
    for (j in 1:(nlatents^2)) {
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- paste0("-------------------------------- Eggers Test for ",
                                           colnames(DRIFTCoeff)[j], "--------------------------------")
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- eggerDrift[[j]]$message
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- summary(eggerDrift[[j]])
    }

    # Unweighted means
    DriftMeans <- colMeans(DRIFTCoeff); DriftMeans
  }



  #######################################################################################################################
  ################################### Fixed & Random Effects Analyses ###################################################
  #######################################################################################################################

  if (fixedAndRandomEffects == TRUE) {

    print(paste0("#################################################################################"))
    print(paste0("######### Fixed and Random Effect Analysis of Single Drift Coefficients #########"))
    print(paste0("#################################################################################"))

    # FIXED EFFECTS ANALYSIS ###############################################################################
    DriftMeans <- colMeans(DRIFTCoeff); DriftMeans
    DRIFTPrecision <- c(rep(1, nlatents^2))/(DRIFTSE); DRIFTPrecision
    # Sum of within weights  and weight * effect size
    T_DriftWeights <- colSums(DRIFTPrecision^2); T_DriftWeights
    T_DriftMeans <- colSums(DRIFTCoeff * DRIFTPrecision^2); T_DriftMeans
    names(T_DriftMeans) <- names(T_DriftWeights); T_DriftMeans
    # Fixed effects results
    FixedEffect_Drift <- T_DriftMeans/T_DriftWeights; FixedEffect_Drift
    FixedEffect_DriftVariance <- 1/T_DriftWeights; FixedEffect_DriftVariance
    FixedEffect_DriftSE <- FixedEffect_DriftVariance^.5; FixedEffect_DriftSE
    FixedEffect_DriftUpperLimit <- FixedEffect_Drift + 1.96*FixedEffect_DriftSE; FixedEffect_DriftUpperLimit
    FixedEffect_DriftLowerLimit <- FixedEffect_Drift - 1.96*FixedEffect_DriftSE; FixedEffect_DriftLowerLimit
    FixedEffect_DriftZ <- FixedEffect_Drift/FixedEffect_DriftSE; FixedEffect_DriftZ
    FixedEffect_DriftProb <- round(1-pnorm(abs(FixedEffect_DriftZ),
                                           mean=c(rep(0, (nlatents^2))), sd=c(rep(1, (nlatents^2))), log=F), digits=digits); FixedEffect_DriftProb
    Q_Drift <- colSums(DRIFTPrecision^2 * DRIFTCoeff^2)- (colSums(DRIFTPrecision^2 * DRIFTCoeff))^2 / colSums(DRIFTPrecision^2); Q_Drift
    H2_Drift <- Q_Drift/(n.studies-1); H2_Drift
    I2_Drift <- (H2_Drift-1)/H2_Drift*100; I2_Drift
    # Tau square
    T2_DriftWeights <- colSums(DRIFTPrecision^2); T2_DriftWeights
    cDrift <- T_DriftWeights-T2_DriftWeights/T_DriftWeights; cDrift
    tau2Drift <- (Q_Drift-(n.studies-1))/cDrift; tau2Drift
    SElnHDrift <- c()
    SElnHDrift[] <- 0
    for (j in 1:(nlatents^2)) {
      if (Q_Drift[j] > n.studies) SElnHDrift[j] <- 1/2*(log(Q_Drift[j])-log(n.studies-1))/((2*Q_Drift[j])^.5-(2*(n.studies-1)-1)^.5)
      if (Q_Drift[j] <= n.studies) SElnHDrift[j] <-  (1/(2*(n.studies-2)) * (1 - 1/(3*(n.studies-2)^.5)) )^.5
    }

    H2DriftUpperLimit <- exp(log(H2_Drift) + 1.96*SElnHDrift); H2DriftUpperLimit
    H2DriftLowerLimit <- exp(log(H2_Drift) - 1.96*SElnHDrift); H2DriftLowerLimit
    L <- exp(0.5*log(Q_Drift/(n.studies-1))-1.96*SElnHDrift)
    U <- exp(0.5*log(Q_Drift/(n.studies-1))+1.96*SElnHDrift)
    I2DriftUpperLimit <- (U^2-1)/U^2 * 100; I2DriftUpperLimit
    I2DriftLowerLimit <- (L^2-1)/L^2 * 100; I2DriftLowerLimit

    MeanOfDriftValues <- DriftMeans
    fixedEffectDriftResults <- rbind(MeanOfDriftValues, FixedEffect_Drift, FixedEffect_DriftVariance, FixedEffect_DriftSE,
                                     FixedEffect_DriftUpperLimit, FixedEffect_DriftLowerLimit,
                                     FixedEffect_DriftZ, FixedEffect_DriftProb, tau2Drift, Q_Drift, H2_Drift,
                                     H2DriftUpperLimit, H2DriftLowerLimit, I2_Drift,
                                     I2DriftUpperLimit, I2DriftLowerLimit)

    # RANDOM EFFECTS ANALYSIS ###############################################################################
    # Total variance weighting
    Ttot_DriftWeights <- 0
    Ttot_DriftMeans <- 0
    tau2DriftExtended <- do.call(rbind, replicate(n.studies, tau2Drift, simplify=FALSE))
    Ttot_DriftWeights <-colSums(1/ (DRIFTSE^2 + tau2DriftExtended)); Ttot_DriftWeights
    Ttot_DriftMeans <- colSums(DRIFTCoeff * 1/ (DRIFTSE^2 + tau2DriftExtended)); Ttot_DriftMeans
    # Random effects results
    RandomEffecttot_Drift <- Ttot_DriftMeans/Ttot_DriftWeights; RandomEffecttot_Drift
    RandomEffecttot_DriftVariance <- 1/Ttot_DriftWeights; RandomEffecttot_DriftVariance
    RandomEffecttot_DriftSE <- RandomEffecttot_DriftVariance^.5; RandomEffecttot_DriftSE
    RandomEffecttot_DriftUpperLimit <- RandomEffecttot_Drift + 1.96*RandomEffecttot_DriftSE; RandomEffecttot_DriftUpperLimit
    RandomEffecttot_DriftLowerLimit <- RandomEffecttot_Drift - 1.96*RandomEffecttot_DriftSE; RandomEffecttot_DriftLowerLimit
    RandomEffecttot_DriftZ <- RandomEffecttot_Drift/RandomEffecttot_DriftSE; RandomEffecttot_DriftZ
    RandomEffecttot_DriftProb <- round(1-pnorm(abs(RandomEffecttot_DriftZ),
                                               mean=c(rep(0, (nlatents^2))), sd=c(rep(1, (nlatents^2))), log=F), digits=digits); RandomEffecttot_DriftProb
    RandomEffectDriftResults <- rbind(RandomEffecttot_Drift, RandomEffecttot_DriftVariance, RandomEffecttot_DriftSE,
                                      RandomEffecttot_DriftUpperLimit, RandomEffecttot_DriftLowerLimit,
                                      RandomEffecttot_DriftZ, RandomEffecttot_DriftProb)
    RandomEffecttot_DriftUpperLimitPI <- RandomEffecttot_Drift + 1.96*(tau2Drift^.5); RandomEffecttot_DriftUpperLimitPI
    RandomEffecttot_DriftLowerLimitPI <- RandomEffecttot_Drift - 1.96*(tau2Drift^.5); RandomEffecttot_DriftLowerLimitPI
    RandomEffectDriftResults <- rbind(RandomEffecttot_Drift, RandomEffecttot_DriftVariance, RandomEffecttot_DriftSE,
                                      RandomEffecttot_DriftUpperLimit, RandomEffecttot_DriftLowerLimit,
                                      RandomEffecttot_DriftZ, RandomEffecttot_DriftProb,
                                      RandomEffecttot_DriftUpperLimitPI, RandomEffecttot_DriftLowerLimitPI)

    ### PET, PEESE & WLS approaches to correct for bias
    PETDrift_fit <- list()
    PEESEDrift_fit <- list()
    WLSDrift_fit <- list()
    PET_PEESEDrift_fit <- list()
    Egger2Drift_fit <- list()

    sampleSizes <- unlist(allSampleSizes)
    for (ii in 1:(nlatents^2)) {
      driftCoeff <- allStudiesDRIFT_effects[ , (ii*2-1)]; driftCoeff
      driftSE <- allStudiesDRIFT_effects[ , (ii*2)]; driftSE

      # PET
      IV <- driftSE; IV
      DV <- driftCoeff; DV
      currentWeigths <- (1/(driftSE^2)); currentWeigths
      PETDrift_fit[[ii]] <- lm(DV ~ IV, weights=currentWeigths); PETDrift_fit[[ii]]

      # Egger's Test (alternative but algebraically identical model)
      Egger2Drift_fit[[ii]] <- t(c(summary(PETDrift_fit[[ii]])$coefficients[2,1:4]))

      # PEESE
      IV <- driftSE^2; IV
      DV <- driftCoeff
      currentWeigths <- (1/(driftSE^2)); currentWeigths
      PEESEDrift_fit[[ii]] <- lm(DV ~ IV, weights=currentWeigths); PEESEDrift_fit[[ii]]

      # PET-PEESE
      if ( (summary(PETDrift_fit[[ii]]))$coefficients[1,4] > .10) {
        PET_PEESEDrift_fit[[ii]] <- PETDrift_fit[[ii]]
        PET_PEESEDrift_fit[[ii]]
      }
      if ( (summary(PETDrift_fit[[ii]]))$coefficients[1,4] <= .10) {
        PET_PEESEDrift_fit[[ii]] <- PEESEDrift_fit[[ii]]
        PET_PEESEDrift_fit[[ii]]
      }

      # WLS
      DV <- driftCoeff/driftSE
      IV <- 1/driftSE
      WLSDrift_fit[[ii]] <- lm(DV ~ IV + 0); WLSDrift_fit[[ii]]; FixedEffect_Drift[[ii]] # should be identical to FixedEffect_Drift
      WLSDriftSE_fit <- summary(WLSDrift_fit[[ii]])$coefficients[2]; WLSDriftSE_fit; FixedEffect_DriftSE[[ii]] # should outperform FixedEffect_DriftSE
    }

    # Combine results
    PET_Drift <-unlist(lapply(PETDrift_fit, function(extract) extract$coefficients))[seq(1, 2*nlatents^2, 2)]; PET_Drift
    PET_SE <- c()
    for (k in 1:(nlatents^2)) PET_SE <- c(PET_SE, summary(PETDrift_fit[[k]])$coefficients[1,2])

    PEESE_Drift <-unlist(lapply(PEESEDrift_fit, function(extract) extract$coefficients))[seq(1, 2*nlatents^2, 2)]; PEESE_Drift
    PEESE_SE <- c()
    for (k in 1:(nlatents^2)) PEESE_SE <- c(PEESE_SE, summary(PEESEDrift_fit[[k]])$coefficients[1,2])

    PET_PEESE_Drift <-unlist(lapply(PET_PEESEDrift_fit, function(extract) extract$coefficients))[seq(1, 2*nlatents^2, 2)]; PET_PEESE_Drift
    PET_PEESE_SE <- c()
    for (k in 1:(nlatents^2)) PET_PEESE_SE <- c(PET_PEESE_SE, summary(PET_PEESEDrift_fit[[k]])$coefficients[1,2])

    WLS_Drift <- unlist(lapply(WLSDrift_fit, function(extract) extract$coefficients)); WLS_Drift
    WLS_SE <- c()
    for (k in 1:(nlatents^2)) WLS_SE <- c(WLS_SE, summary(WLSDrift_fit[[k]])$coefficients[1,2])

    Egger2Drift_results <- matrix(unlist(Egger2Drift_fit), ncol=nlatents^2, nrow=4); Egger2Drift_results

    PET_PEESE_DRIFTresults <- rbind(PET_Drift, PET_SE,
                                    PEESE_Drift, PEESE_SE,
                                    PET_PEESE_Drift, PET_PEESE_SE,
                                    WLS_Drift, WLS_SE,
                                    Egger2Drift_results)
    colnames(PET_PEESE_DRIFTresults) <- colnames(DRIFTCoeff)
    rownames(PET_PEESE_DRIFTresults) <- c(rownames(PET_PEESE_DRIFTresults)[1:8], "Egger's b0", "SE(b0)", "T", "p")

  } ### END Fixed & Random Effects Analyses ###



  #######################################################################################################################
  ################################################## Plotting ###########################################################
  #######################################################################################################################

  {
    # Function to compute discrete parameters using drift parameters and time-scaling factors
    discreteDrift <-function(driftMatrix, timeScale, number) {
      discreteDriftValue <- expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }

    if (plotCrossEffects == TRUE || plotAutoEffects == TRUE) {

      print(paste0("#################################################################################"))
      print(paste0("################################### Plotting ####################################"))
      print(paste0("#################################################################################"))

      ##################################### SELECT DRIFT MATRICES ########################################

      # Drift matrix used for plotting effects of primary studies
      DriftForPlot <- array(dim=c(n.studies, nlatents, nlatents))
      for (h in 1:n.studies) {
        DriftForPlot[h, 1:nlatents,1:nlatents] <- matrix(DRIFTCoeff[h, ], nlatents, nlatents)
      }

      # Drift matrix used for plotting effects of moderator models
      if (testModeratorModel == TRUE) {

        DriftForPlot.M.D <- array(dim=c(length(unique(unlist(allModerators))), (nlatents^2) ))
        for (j in 1:length(unique(unlist(allModerators)))) {
          DriftForPlot.M.D[j, ] <- modModel1mx_Drift_Coef[((j-1)*(nlatents^2)+1) : ((j-1)*(nlatents^2)+(nlatents^2))]
        }

        DriftForPlot.M <- array(dim=c((nlatents^2), length(unique(unlist(allModerators))), (nlatents^2)))
        latentCounter <- 0
        for (h in 1:(nlatents^2)) {
          latentCounter <- latentCounter + 1
          moderatorCounter <- 0
          for (j in unique(unlist(allModerators))) {
            moderatorCounter <- moderatorCounter + 1
            if (j == unique(unlist(allModerators))[1]) {
              DriftForPlot.M[h, j, ] <- modModel2mx_Drift_Coef[[latentCounter]][1:(nlatents^2)]
            }
            if (j != unique(unlist(allModerators))[1]) {
              DriftForPlot.M[h, j, ] <- DriftForPlot.M[h, 1, ]
              DriftForPlot.M[h, j, latentCounter] <- modModel2mx_Drift_Coef[[latentCounter]][nlatents^2 + moderatorCounter-1]
            }
          }
        }

      } # END Drift Matrix used for plotting effects of moderator models

      # Drift matrix used for plotting effects of CoTiMA model
      if (testDRIFTallModel == TRUE) DriftForPlot.D <- matrix(homDRIFTall_Drift_Coef, nlatents, nlatents)

      # Create objects to store the dots for plotting
      plotPairs <- array(dim=c(n.studies, noOfSteps, 2+nlatents^2))
      dotPlotPairs <- array(dim=c(n.studies, noOfSteps, 2+nlatents^2))
      if (testDRIFTallModel == TRUE) plotPairs.D <- array(dim=c(noOfSteps, 2+nlatents^2))
      if (testModeratorModel == TRUE) {
        plotPairs.M <- array(dim=c(nlatents^2, length(unique(unlist(allModerators))), noOfSteps, 2+nlatents^2))
        plotPairs.M.D <- array(dim=c(length(unique(unlist(allModerators))), noOfSteps, 2+nlatents^2))
      }

      ##################################### COMPUTE DOTS FOR PLOTTING ########################################
      ## Loop through all noOfSteps (full time range)
      for (currentTimeScale in 0:noOfSteps){

        # All primary studies
        for (h in 1:n.studies) {
          timeValue <- usedTimeRange[currentTimeScale+1]; timeValue
          plotPairs[h,currentTimeScale,1] <- currentTimeScale; plotPairs[h,currentTimeScale,1]
          plotPairs[h,currentTimeScale,2] <- timeValue; plotPairs[h,currentTimeScale,2]
          for (j in 1:(nlatents^2)) {
            plotPairs[h,currentTimeScale,(2+j)] <- discreteDrift((DriftForPlot[h,1:nlatents,1:nlatents]),timeValue, j)
            if (timeValue %in% studyList[[h]]$delta_t) {
              dotPlotPairs[h, currentTimeScale, 1] <- currentTimeScale
              dotPlotPairs[h, currentTimeScale, 2] <- timeValue
              dotPlotPairs[h, currentTimeScale, (2+j)] <- discreteDrift((DriftForPlot[h,1:nlatents,1:nlatents]),timeValue, j)
            }
          }
        } # END All primary studies

        # CoTiMA (DRIFTallModel)
        if (testDRIFTallModel == TRUE) {
          plotPairs.D[currentTimeScale, 1] <- currentTimeScale
          plotPairs.D[currentTimeScale, 2] <- timeValue
          for (j in 1:(nlatents^2)) {
            plotPairs.D[currentTimeScale,(2+j)] <- discreteDrift((DriftForPlot.D[1:nlatents,1:nlatents]),timeValue, j)
          }
        } # END (DRIFTallModel)

        # Moderator models
        if (testModeratorModel == TRUE) {
          for (h in 1:length(unique(unlist(allModerators)))) {
            timeValue <- currentTimeScale * stepWidth; timeValue
            plotPairs.M[, h,currentTimeScale,1] <- currentTimeScale
            plotPairs.M[, h,currentTimeScale,2] <- timeValue
            for (j in 1:(nlatents^2)) {
              plotPairs.M[j, h,currentTimeScale,(2+j)] <- discreteDrift(matrix(DriftForPlot.M[j, h, ], nlatents, nlatents), timeValue, j)
            }
            plotPairs.M.D[h, currentTimeScale,1] <- currentTimeScale; plotPairs.M.D[h,currentTimeScale,1]
            plotPairs.M.D[h, currentTimeScale,2] <- timeValue; plotPairs.M.D[h,currentTimeScale,2]
            plotPairs.M.D[h, currentTimeScale, 3:(3+nlatents^2-1)] <- discreteDrift(matrix(DriftForPlot.M.D[h, ], nlatents, nlatents), timeValue, 1:(nlatents^2))
          }
        } # END Moderator models
      } ## END Loop through all noOfSteps (full time range)

      if (!(exists("plotPairs.D"))) plotPairs.D <- array(dim=c(2, 2+nlatents^2))

      ##################################### PLOTTING PARAMETERS ##########################################
      yMin <- yLimitsForEffects[1]; yMin
      yMax <- yLimitsForEffects[2]; yMax
      xMax <- max(usedTimeRange); xMax
      xMin <- usedTimeRange[1]; xMin
      targetRows <- max(usedTimeRange)/stepWidth; targetRows

      ############################################ PLOTTING ##############################################

      ## PLOT (auto effects) ##############################
      if (plotAutoEffects == TRUE) {
        figureFileNameAuto <- c()

        #### Auto-regressive Effects ###
        counter <- 0
        for (j in seq(1, nlatents^2, (nlatents+1))) {
          counter <- counter + 1
          for (h in 1:n.studies) {
            plot(plotPairs[h, ,2], plotPairs[h, ,2+j], type="l", col="grey", lwd=1.5,
                 xlim = c(xMin, xMax), ylim = c(yMin, 1),
                 xaxt='n', yaxt='n', ann=FALSE)
            par(new=T)
            plot(dotPlotPairs[h, ,2], plotPairs[h, ,2+j], type="b", col="black", lwd=.5, pch=16,
                 xlim = c(xMin, xMax), ylim = c(yMin, 1),
                 xaxt='n', yaxt='n', ann=FALSE)
            par(new=T)
          }
          if (testDRIFTallModel == TRUE) {
            plot(plotPairs.D[,2], plotPairs.D[,2+j], type="l", lty=2, col="black", lwd=2.5,
                 xlim = c(xMin, xMax), ylim = c(yMin, 1),
                 xaxt='n', yaxt='n', ann=FALSE)
            par(new=T) }
          plot(c(0,0), type="l", col="white", lwd=1.5,
               xlim = c(xMin, xMax), ylim = c(yMin, 1),
               xaxt='n',
               ann=FALSE)

          # Correct differences in axis length
          atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
          labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
          if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
          if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
          axis(side=1, at = atSeq*stepWidth, labels=labelsSeq)
          # Add labels and title
          title(main = paste0("Auto-regressive Effects of V", counter), sub = NULL,
                xlab=paste0("Time Lag in ", timeUnit), ylab = "Auto-regressive Beta")

          par(new=F)

          figureFileNameAuto[counter] <- paste0(filePrefix,"_CoTiMA results for auto-regressive effects of V", counter, ".png")
          dev.copy(png, figureFileNameAuto[counter], width = 8, height = 8, units = 'in', res = 300)
          dev.off()
        }
      } ## END PLOT (auto effects)

      ## PLOT (moderated auto effects)
      if (plotAutoEffects == TRUE & testModeratorModel == TRUE) {

        figureFileNameAutoModeratedSingle <- figureFileNameAutoModeratedFull <- c()

        #### MODERATED Auto-regressive Effects of X
        ## (only a single effect moderated in each model)
        counter <- 0
        for (j in seq(1, nlatents^2, (nlatents+1))) {
          counter <- counter + 1
          for (h in 1:length(unique(unlist(allModerators)))) {
            modColor <- h
            modLWD <- h # h/2
            plot(plotPairs.M[j, h, ,2], plotPairs.M[j, h, ,2+j], type="l", col=modColor, lwd=modLWD,
                 xlab=paste0("Time Lag in ", timeUnit), ylab ="",
                 xlim = c(xMin,xMax), ylim=c(yMin,yMax),
                 xaxt='n',
                 ann=TRUE)

            legend(0, yMin+.1, legend=paste0("Moderator Value = ", unique(unlist(allModerators))),
                   col=1:length(unique(unlist(allModerators))),
                   lwd=1:length(unique(unlist(allModerators))), cex=0.8)

            # CHD 26.08.2019: correct differences in axis length
            atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
            labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
            if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
            if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
            axis(side=1, at = atSeq*stepWidth, labels=labelsSeq)

            # CHD 26.08.2019 add labels and title
            title(main = paste0("Moderated Auto-regressive Effects (only this effects was moderated) of V", counter), sub = NULL,
                  xlab=paste0("Time Lag in ", timeUnit), ylab = "Auto-regressive Beta")

            par(new=T)

            figureFileNameAutoModeratedSingle[counter] <- paste0(filePrefix,"_CoTiMA results for moderated auto-regressive (only this effect moderated) effects of V", counter, ".png")
            dev.copy(png, figureFileNameAutoModeratedSingle[counter], width = 8, height = 8, units = 'in', res = 300)
            dev.off()
          }
          par(new=F)
        }

        ## (all effects moderated in the model)
        counter <- 0
        for (j in seq(1, nlatents^2, (nlatents+1))) {
          counter <- counter + 1
          for (h in 1:length(unique(unlist(allModerators)))) {
            dim(plotPairs.M.D)
            modColor <- h
            modLWD <- h
            plot(plotPairs.M.D[h, ,2], plotPairs.M.D[h, ,2+j], type="l", col=modColor, lwd=modLWD,
                 xlab=paste0("Time Lag in ", timeUnit), ylab ="",
                 xlim = c(xMin,xMax), ylim=c(yMin,yMax),
                 xaxt='n',
                 ann=TRUE)

            legend(0, yMin+.1, legend=paste0("Moderator Value = ", unique(unlist(allModerators))),
                   col=1:length(unique(unlist(allModerators))),
                   lwd=1:length(unique(unlist(allModerators))), cex=0.8)

            # Correct differences in axis length
            atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
            labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
            if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
            if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
            axis(side=1, at = atSeq*stepWidth, labels=labelsSeq)

            # Add labels and title
            title(main = paste0("Moderated Auto-regressive Effects (all effects were moderated) of V", counter), sub = NULL,
                  xlab=paste0("Time Lag in ", timeUnit), ylab = "Auto-regressive Beta")

            par(new=T)

            figureFileNameAutoModeratedFull[counter] <- paste0(filePrefix,"_CoTiMA results for moderated auto-regressive (all effects in the model were moderated) effects of V", counter, ".png")
            dev.copy(png, figureFileNameAutoModeratedFull[counter], width = 8, height = 8, units = 'in', res = 300)
            dev.off()
          }
          par(new=F)

        }
      } # END PLOT (moderated auto effects)


      ## PLOT (cross effects) ##############################
      if (plotCrossEffects == TRUE & nlatents > 1) {

        figureFileNameCross <- c()

        coeffSeq <- seq(1, nlatents^2, 1)[!(seq(1, nlatents^2, 1) %in% seq(1, nlatents^2, (nlatents+1)))]; coeffSeq
        counter <- 0
        for (j in coeffSeq) {
          counter <- counter + 1
          for (h in 1:n.studies) {
            plot(plotPairs[h, ,2], plotPairs[h, ,2+j], type="l", col="grey", lwd=1.5,
                 xlim = c(xMin, xMax), ylim = c(yMin, yMax),
                 xaxt='n', yaxt='n', ann=FALSE)
            par(new=T)
            plot(dotPlotPairs[h, ,2], plotPairs[h, ,2+j], type="b", col="black", lwd=5, pch=16,
                 xlim = c(xMin, xMax), ylim = c(yMin, yMax),
                 xaxt='n', yaxt='n', ann=FALSE)
            par(new=T)
            abline(v=meanDelta, lty=2)
            par(new=T)
          }
          if (testDRIFTallModel == TRUE) {
            plot(plotPairs.D[,2], plotPairs.D[,2+j], type="l", lty=2, col="black", lwd=2.5,
                 xlim = c(xMin, xMax), ylim = c(yMin, yMax),
                 xaxt='n', yaxt='n', ann=FALSE)
            par(new=T) }
          plot(c(0,0), type="l", col="white", lwd=1.5,
               xlim = c(xMin, xMax), ylim = c(yMin, yMax),
               xaxt='n',
               ann=FALSE)
          par(new=F)

          # Correct differences in axis length
          atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
          labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
          if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
          if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
          axis(side=1, at = atSeq*stepWidth, labels=labelsSeq)
          # Add labels and title
          title(main = paste0("Cross-lagged Effects of ", driftNames[j]), sub = NULL,
                xlab=paste0("Time Lag in ", timeUnit), ylab = "Cross-lagged Beta")

          par(new=F)

          if (plotCrossEffects == TRUE) {
            figureFileNameCross[counter] <- paste0(filePrefix,"_CoTiMA results for cross-lagged effects of ", driftNames[j], ".png")
            dev.copy(png, figureFileNameCross[counter], width = 8, height = 8, units = 'in', res = 300)
            dev.off()
          }
        }
      }

      if (publicationBias == TRUE) {
        # funnel plots #########################################################

        stretch <- 1.2
        figureName <- c()
        driftPlots <- list()
        for (i in 1:(nlatents^2)) {
          # Determine range of X- and Y-axes
          pairs <- cbind(DRIFTCoeff[,i], DRIFTSE[,i]); pairs
          minX <- min(DRIFTCoeff[,i]); minX
          maxX <- max(DRIFTCoeff[,i]); maxX
          maxAbsX <- max(abs(c(minX, maxX))); maxAbsX
          avgX <-FixedEffect_Drift[i]; avgX
          yMax2 <- stretch * max(DRIFTSE[,i]); yMax2
          # Determine pseudo confidence intervals (http://www.metafor-project.org/doku.php/plots:funnel_plot_variations)
          lowLeft <- FixedEffect_Drift[i] - 1.96 * yMax2; lowLeft
          lowRight <- FixedEffect_Drift[i] + 1.96 * yMax2; lowRight
          xMin2 <- 0 - max(abs(lowLeft), abs(lowRight), abs(DRIFTCoeff[,i])); xMin2
          xMax2 <- 0 + max(abs(lowLeft), abs(lowRight), abs(DRIFTCoeff[,i])); xMax2

          plot.new()
          plot(pairs, xlab=colnames(DRIFTCoeff)[i], ylab="Standard Error", xlim = c(xMin2, xMax2), ylim = c(yMax2, 0))
          polygon(c(lowLeft, avgX, lowRight), c(stretch * yMax2, 0, stretch * yMax2), lty=3)
          abline(v=avgX, col="black", lwd=1.5, lty=1)

          figureContent <- colnames(DRIFTCoeff)[i]
          figureName[i] <- paste0(filePrefix, "_",  "Funnel Plot for ", figureContent, ".png")
          dev.copy(png, figureName[i], width = 8, height = 8, units = 'in', res = 300)
          dev.off()
        }
      }

      # Plot moderator effects #########################################################

      if (plotCrossEffects == TRUE & nlatents > 1 & testModeratorModel == TRUE) {

        figureFileNameCrossModeratedSingle <- figureFileNameCrossModeratedFull <-  c()
        counter <- c()
        coeffSeq <- seq(1, nlatents^2, 1)[!(seq(1, nlatents^2, 1) %in% seq(1, nlatents^2, (nlatents+1)))]; coeffSeq

        #### Moderated Cross-lagged Effects
        ## (only a single effect moderated in each model)
        counter <- 0
        for (j in coeffSeq) {
          counter <- counter + 1
          for (h in 1:length(unique(unlist(allModerators)))) {
            modColor <- h
            modLWD <- h # h/2
            plot(plotPairs.M[j, h, ,2], plotPairs.M[j, h, ,2+j], type="l", col=modColor, lwd=modLWD,
                 xlab="", ylab ="",
                 xlim = c(xMin,xMax), ylim=c(yMin,yMax),
                 xaxt='n',
                 ann=TRUE)

            legend(0, yMax, legend=paste0("Moderator Value = ", unique(unlist(allModerators))),
                   col=1:length(unique(unlist(allModerators))),
                   lwd=1:length(unique(unlist(allModerators))), cex=0.8)

            # Correct differences in axis length
            atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
            labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
            if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
            if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
            axis(side=1, at = atSeq*stepWidth, labels=labelsSeq)

            # Add labels and title
            title(main = paste0("Cross-lagged Effects (only this effects was moderated) of ", driftNames[j]), sub = NULL,
                  xlab=paste0("Time Lag in ", timeUnit), ylab = "Cross-lagged Beta")

            par(new=T)

            figureFileNameCrossModeratedSingle[counter] <- paste0(filePrefix,"_CoTiMA results for moderated cross-lagged (only this effect moderated) effects of V", counter, ".png")
            dev.copy(png, figureFileNameCrossModeratedSingle[counter], width = 8, height = 8, units = 'in', res = 300)
            dev.off()
          }
          par(new=F)
        }

        ## (all effects moderated in the model)
        counter <- 0
        for (j in coeffSeq) {
          counter <- counter + 1
          for (h in 1:length(unique(unlist(allModerators)))) {
            modColor <- h
            modLWD <- h
            plot(plotPairs.M[j, h, ,2], plotPairs.M[j, h, ,2+j], type="l", col=modColor, lwd=modLWD,
                 xlab="", ylab ="",
                 xlim = c(xMin,xMax), ylim=c(yMin,yMax),
                 xaxt='n',
                 ann=TRUE)

            legend(0, yMax, legend=paste0("Moderator Value = ", unique(unlist(allModerators))),
                   col=1:length(unique(unlist(allModerators))),
                   lwd=1:length(unique(unlist(allModerators))), cex=0.8)

            # Correct differences in axis length
            atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
            labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
            if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
            if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
            axis(side=1, at = atSeq*stepWidth, labels=labelsSeq)

            # Add labels and title
            title(main = paste0("Cross-lagged Effects (all effects were moderated) of ", driftNames[j]), sub = NULL,
                  xlab=paste0("Time Lag in ", timeUnit), ylab = "Cross-lagged Beta")

            par(new=T)

            figureFileNameCrossModeratedFull[counter] <- paste0(filePrefix,"_CoTiMA results for moderated cross-lagged (all effects in the model were moderated) effects of V", counter, ".png")
            dev.copy(png, figureFileNameCrossModeratedFull[counter], width = 8, height = 8, units = 'in', res = 300)
            dev.off()
          }
          par(new=F)
        }
      }
    } ## END if (plotCrossEffects == TRUE || plotAutoEffects == TRUE) ##

    #### Required Sample Sizes ##################################################################################

    if (length(statisticalPower) > 0 ) {

      figureFileNameStatisticalPower <- c()

      yMin <-  yLimitsForPower[1]; yMin
      yMax <-  yLimitsForPower[2]; yMax
      yStep <- round(abs(yMin-yMax)/100, 0); yStep

      targetRows <- max(usedTimeRange)/stepWidth; targetRows

      xMax <- targetRows; xMax
      xMin <- usedTimeRange[1]; xMin

      currentColour <- rep(c("black", "grey"), length(statisticalPower)); currentColour
      currentLWD <- c(3, 2, 1); currentLWD

      plot.new()
      for (j in 1:length(currentDriftNames)) {
        offset <- (j-1)*length(statisticalPower); offset
        par(new=F)
        for (h in 1:length(statisticalPower)) {
          plot(requiredSampleSizes[1:targetRows, (h+offset)],
               type="l", col=currentColour[h], lwd=currentLWD[h],
               main=paste0("Required Sample Size For the Effect of ", currentDriftNames[j]),
               xlab=paste0("Time Lag in ", timeUnit),
               ylab="Required Sample Size",
               xlim = c(xMin, xMax),
               ylim = c(yMin, yMax),
               xaxt='n' #, yaxt='n', ann=FALSE
          )
          par(new=T)
        }

        # Correct differences in axis length
        atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
        labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
        if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]
        if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]
        axis(side=1, at = atSeq, labels=labelsSeq)
        legend('bottomright', legend=statisticalPower, lty=1, col=currentColour, lwd=currentLWD, bty='n', cex=.75)

        figureFileNameStatisticalPower[j] <- paste0(filePrefix,"_RequiredSampleSizesFor", currentDriftNames[j], ".png")
        dev.copy(png, figureFileNameStatisticalPower[j], width = 8, height = 8, units = 'in', res = 300)
        dev.off()

      }

    } # end if (length(statisticalPower) > 0 )
  } ### END Plotting ###


  end.time <- Sys.time()
  time.taken <- end.time - start.time
  st <- paste0("Computation started at: ", start.time); st
  et <- paste0("Computation ended at: ", end.time); et
  tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt



  ###########################################################################################################
  ###################################### SAVE RESULTS AND FIGURES ###########################################
  ###########################################################################################################

  print(paste0("#################################################################################"))
  print(paste0("########################### Save Results and Figures ############################"))
  print(paste0("#################################################################################"))

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
  cat("Number of latent variables (nlatents): ", nlatents,  sep="")
  cat(" ", "", sep="\n")
  cat("Number of iterations to be used by ctsem (retryattempts): ", retryattempts,  sep="")
  cat(" ", "", sep="\n")
  cat("Number of re-fits used by this CoTiMA Function (refits): ", refits,  sep="")
  cat(" ", "", sep="\n")
  if(is.null(extraRefits)) {
    x <- "empty"
  } else {
    x <- extraRefits
  }
  cat("(Vector of) study number(s) that should be re-fitted even more (difficult models) (extraRefits): ", x,  sep="")
  cat(" ", "", sep="\n")
  cat("Factor by which refits are multiplied for difficult models (factorExtraRefits): ", factorExtraRefits,  sep="")
  cat(" ", "", sep="\n")
  cat("Multigoup model with all parameters freely estimated across groups (testHeterogeneityModel): ", testHeterogeneityModel,  sep="")
  cat(" ", "", sep="\n")
  cat("Multigoup model (CoTiMA!) with all DRIFT parameters invariantly estimated across groups (testDRIFTallModel): ", testDRIFTallModel,  sep="")
  cat(" ", "", sep="\n")
  cat("Multigoup models with every single DRIFT parameters invariantly estimated across groups (testDRIFTSingleModel): ", testDRIFTSingleModel,  sep="")
  cat(" ", "", sep="\n")
  cat("Multigoup model to test equality of cross effects (testDRIFTCrossModel): ", testDRIFTCrossModel,  sep="")
  cat(" ", "", sep="\n")
  cat("Multigoup model to test equality of auto effects (testDRIFTAutoModel): ", testDRIFTAutoModel,  sep="")
  cat(" ", "", sep="\n")
  cat("Estimate confidence intervals (for all requested models) (confidenceIntervals): ", confidenceIntervals,  sep="")
  cat(" ", "", sep="\n")
  cat("Performs fixed and random effects analysis, funnel plots, Egger's test (fixedAndRandomEffects): ", fixedAndRandomEffects,  sep="")
  cat(" ", "", sep="\n")
  cat("Plot cross effects across a range of discrete intervals (plotCrossEffects): ", plotCrossEffects,  sep="")
  cat(" ", "", sep="\n")
  cat("Plot auto effects across a range of discrete intervals (plotAutoEffects): ", plotAutoEffects,  sep="")
  cat(" ", "", sep="\n")
  cat("Request NPSOL optimizer (NPSOL): ", NPSOL,  sep="")
  cat(" ", "", sep="\n")
  cat("Allows parallel processing on unix-like computers (e.g., Mac;) using # cores: ", coresToUse,  sep="")
  cat(" ", "", sep="\n")
  cat("Use alternative multigroup fitting procedure (useCTMultigroupFitAlt): ", useCTMultigroupFitAlt,  sep="")
  cat(" ", "", sep="\n")
  cat("Performs moderator analysis (testModeratorModel): ", testModeratorModel,  sep="")
  cat(" ", "", sep="\n")
  cat("Selects from the vector of moderator values (the 1st is used if not otherwise specificed) (moderatorNumber): ", moderatorNumber,  sep="")
  cat(" ", "", sep="\n")
  cat("Rounding used in output (digits): ", digits,  sep="")
  cat(" ", "", sep="\n")
  cat("Displays estimates from single study ctsem models and waits for user inoput to continue (checkSingleStudyResults): ", checkSingleStudyResults,  sep="")
  cat(" ", "", sep="\n")
  cat("Test model with unit effects/unobserved heterogeneity (manifesttraitvar): ", manifesttraitvar,  sep="")
  cat(" ", "", sep="\n")
  if (is.null(saveSingleStudyModelFit)) x <- "no" else {
    x <- paste(rep(saveSingleStudyModelFit[1], length(saveSingleStudyModelFit)-1),
               " studyFit", saveSingleStudyModelFit[-1], ".rds", sep="") }
  cat("Save the fit of single study ctsem models (saveSingleStudyModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadSingleStudyModelFit)) x <- "no" else {
    x <- paste(rep(loadSingleStudyModelFit[1], length(loadSingleStudyModelFit)-1),
               " studyFit", loadSingleStudyModelFit[-1], ".rds", sep="") }
  cat("load the fit of single study ctsem models (loadSingleStudyModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(saveDRIFTAllModelFit)) x <- "no" else {
    x <- paste(saveDRIFTAllModelFit, " homDRIFTallFit.rds", sep="") }
  cat("Save the fit of the full CoTiMA model (saveDRIFTAllModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadDRIFTAllModelFit)) x <- "no" else {
    x <- paste(loadDRIFTAllModelFit, " homDRIFTallFit.rds", sep="") }
  cat("Load the fit of the full CoTiMA model (loadDRIFTAllModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(saveDRIFTSingleModelFit)) x <- "no" else {
    x <- paste(rep(saveDRIFTSingleModelFit, nlatents^2),
               " homDRIFTSingleFit", seq(1, nlatents^2, 1), ".rds", sep="") }
  cat("Save the fit of the CoTiMA models with single invariant Drift Parameters (saveDRIFTSingleModelFit): \n", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadDRIFTSingleModelFit)) x <- "no" else {
    x <- paste(rep(loadDRIFTSingleModelFit, nlatents^2),
               " homDRIFTSingleFit", seq(1, nlatents^2, 1), ".rds", sep="") }
  cat("Load the fit of the CoTiMA models with single invariant Drift Parameters (loadDRIFTSingleModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(saveDRIFTCrossModelFit)) x <- "no" else {
    x <- paste(saveDRIFTCrossModelFit, " homDRIFTCrossModelFit.rds", "\n",
               saveDRIFTCrossModelFit, " homEquDRIFTCrossModelFit.rds", sep="") }
  cat("Save the fit of the CoTiMA models testing the equality of the cross effects (saveDRIFTCrossModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadDRIFTCrossModelFit)) x <- "no" else {
    x <- paste(loadDRIFTCrossModelFit, " homDRIFTCrossModelFit.rds", "\n",
               loadDRIFTCrossModelFit, " homEquDRIFTCrossModelFit.rds", sep="") }
  cat("Load the fit of the CoTiMA models testing the equality of the cross effects (loadDRIFTCrossModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(saveDRIFTAutoModelFit)) x <- "no" else {
    x <- paste(saveDRIFTAutoModelFit, " homDRIFTAutoModelFit.rds", "\n",
               saveDRIFTAutoModelFit, " homEquDRIFTAutoModelFit.rds", sep="") }
  cat("Save the fit of the CoTiMA models testing the equality of the auto effects (saveDRIFTAutoModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadDRIFTAutoModelFit)) x <- "no" else {
    x <- paste(loadDRIFTAutoModelFit, " homDRIFTAutoModelFit.rds", "\n",
               loadDRIFTAutoModelFit, " homEquDRIFTAutoModelFit.rds", sep="") }
  cat("Load the fit of the CoTiMA models testing the equality of the auto effects (loadDRIFTAutoModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(saveHeterogeneityModelFit)) x <- "no" else {
    x <- paste(saveHeterogeneityModelFit, " hetModelFit.rds", sep="") }
  cat("Save the fit of the heterogeneity model (saveHeterogeneityModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  if(is.null(loadHeterogeneityModelFit)) x <- "no" else {
    x <- paste(loadHeterogeneityModelFit, " hetModelFit.rds", sep="") }
  cat("Load the fit of the heterogeneity model (loadHeterogeneityModelFit): ", x, sep="\n")
  cat(" ", "", sep="\n")
  cat("The working directory (workingDirectory) is: ", workingDirectory,  sep="\n")
  cat(" ", "", sep="\n")
  cat("The prefix for the result file that is created (resultsFileName): ", resultsFileName,  sep="")
  cat(" ", "", sep="\n")
  cat("The prefix for the figure files and the model fit files that are created (filePrefix): ", filePrefix,  sep="")
  cat(" ", "", sep="\n")
  cat(" ", "", sep="\n")
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("------------------------------ Primary Study Statistics---------------------------------", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat(paste("-------------------------------- Sample Sizes ------------------------------------------", "", sep="\n"))
  cat("Sample Sizes (k):       ", unlist(allSampleSizes), "", sep="\n")
  cat(paste("Reported sample sizes should be treated with some caution. Primary studies with         ", "", sep="\n"))
  cat(paste("correlations matrices with pairwise deleted missing values or raw data with missing     ", "", sep="\n"))
  cat(paste("do not have a single number representing the sample size. However, the internal         ", "", sep="\n"))
  cat(paste("algorithm used to transform the correlation matrix with pairwise deleted missing values ", "", sep="\n"))
  cat(paste("into pseudo raw data estimates how large the sample size had to be at least in order to ", "", sep="\n"))
  cat(paste("reproduce the pattern of pairwise N provided for the respective primary studies.", "", sep="\n"))
  cat(" ", "", sep="\n")
  cat("Overall Sample Size (N):", overallSampleSize)
  cat(" ", "", sep="\n")
  for (i in 1:n.studies) {
    if (dim(primaryStudies$pairwiseNs[[i]])[1] > 0) {
      cat(paste("Primary study No.", i, "had a correlation matrix with pairwise deleted cases. N was:    ", "", sep=" "))
      print(round(primaryStudies$pairwiseNs[[i]], digits))
      cat(" ", "", sep="\n")
      cat(paste("Some cases were lost during pseudo raw data generation. Lost N was                      ", "", sep="\n"))
      print(round(lostN[[i]], digits))
      cat(paste("Overall lost N was: ", overallNDiff[[i]], "", sep="\n"))
      cat(paste("Relative lost N was: ", round(relativeNDiff[[i]], digits), "", sep="\n"))
      cat(" ", "", sep="\n")
    }
  }
  cat(paste("---------------------------------- Time Lags -------------------------------------------", "", sep="\n"))
  cat("All Time Lags (deltas): ", sort(unlist(allDeltas)))
  cat(" ", "", sep="\n")
  cat("Mean Time Lag (delta):  ", meanDelta)
  cat(" ", "", sep="\n")
  cat(paste("--------------------------------- Time Points ------------------------------------------", "", sep="\n"))
  cat("All Time Points (Waves):", sort(unlist(allTpoints)))
  cat(" ", "", sep="\n")
  cat("Sum of all Time Points: ", overallTpoints)
  cat(" ", "", sep="\n")
  cat("Mean No. of Time Points:", round(meanTpoints, digits))
  cat(" ", "", sep="\n")
  if (testModeratorModel == TRUE) {
    cat(paste("--------------------------------- Moderators -------------------------------------------", "", sep="\n"))
    cat("Moderator values of all studies:", round(unlist(allModerators), digits))
    cat(" ", "", sep="\n")
    cat("Number of different moderator values:", round(numberOfModerators, digits))
    cat(" ", "", sep="\n")
    cat("Moderator values found:", categoriesOfModerators)
    cat(" ", "", sep="\n")
    x <- table(unlist(allModerators))
    x <- cbind((paste0(rep("value = ", length(names(x))), names(x))), (table(unlist(allModerators)))[2])
    colnames(x) <- c("Moderator Value", "Frequency")
    cat("Frequencies of moderator values: ", sep="\n")
    print(x)
    cat(" ", "", sep="\n")
  }
  cat(" ", "", sep="\n")
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("--------------- All Studies Analyzed Separately (~ Heterogeneity Model) ----------------", "", sep="\n"))
  cat(paste("########################################################################################", "", sep="\n"))
  cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
  cat("All Drift Coefficients and their Standard Errors (SE):")
  cat(" ", "", sep="\n")
  print(round(allStudiesDRIFT_effects, digits))
  cat(" ", "", sep="\n")
  if (confidenceIntervals == TRUE) {
    cat("All Drift Coefficients and their 95% Confidence Intervals (Not Necessarily Symmetric):")
    cat(" ", "", sep="\n")
    for (i in 1:n.studies) {print(paste0("Study No. ", i)); print(round(allDriftCI[[i]][, 1:3], 3))}
  }
  cat(" ", "", sep="\n")
  cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
  cat("-2 Log Likelihood: ", round(allStudies_Minus2LogLikelihood, digits), sep="")
  cat(" ", "", sep="\n")
  cat("Overall Number of Estimated Parameters: ", round(allStudies_estimatedParameters, digits), sep="")
  cat(" ", "", sep="\n")
  cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(allStudies_df, digits), sep="")
  cat(" ", "", sep="\n")
  cat(" ", "", sep="\n")

  if (testHeterogeneityModel == TRUE) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("--------- Actual Heterogeneity Model (Results Should be Identical to Those Above) ------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat("All Drift Coefficients and their Standard Errors (SE):")
    cat(" ", "", sep="\n")
    print(round(hetModelDRIFT_effects, digits))
    cat(" ", "", sep="\n")
    if (confidenceIntervals == TRUE) {
      cat("All Drift Coefficients and their 95% Confidence Intervals (Not Necessarily Symmetric):")
      cat(" ", "", sep="\n")
      print(round(hetModelCI[, 1:3], digits))
    }
    cat(" ", "", sep="\n")
    cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
    cat("-2 Log Likelihood: ", round(hetModel_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Overall Number of Estimated Parameters: ", round(hetModel_estimatedParameters, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(hetModel_df, digits), sep="")
    cat(" ", "", sep="\n")
  }

  if (testDRIFTallModel == TRUE) {
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
    cat("Comparison with Unconstrained Modell (All Samples Analyzed Separately (~Heterogeneity Model)):")
    cat(" ", "", sep="\n")
    cat("Delta(df): ", round(allStudies_homDRIFTall_df, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Delta(-2LL) (= Delta(Chi-square)): ", round(allStudies_homDRIFTall_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("p-value (with double number of digits): ", round(allStudies_homDRIFTall_prob, 2*digits), sep="")
    cat(" ", "", sep="\n")
    cat("A sign. value (p < .05) indicates that entire process (i.e., the full Drift Matrix) varies among primary studies.")
    cat(" ", "", sep="\n")
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

  if (testModeratorModel == TRUE) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("------ Full Drift Moderator Model (all Drift Coefficients Vary Among Moderators) -------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    print(round(newModModel1DRIFT_effects, digits))
    cat(" ", "", sep="\n")
    cat("-2 Log Likelihood: ", round(modModel1mx_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Overall Number of Estimated Parameters: ", round(modModel1mx_estimatedParameters, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Comparison with Multi-Sample Homogeneity Model (all Drift Coefficients Invariant):")
    cat(" ", "", sep="\n")
    cat("Delta(df): ", round(homDRIFTall_modModel1_df, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Delta(-2LL) (= Delta(Chi-square)): ", round(homDRIFTall_modModel1_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("p-value (with double number of digits): ", round(homDRIFTall_modModel1_prob, 2*digits), sep="")
    cat(" ", "", sep="\n")
    cat("A sign. value (p < .05) indicates that entire process (i.e., the full Drift Matrix) varies ")
    cat("conditional on moderator values.")
    cat(" ", "", sep="\n")

    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-- Single Drift Moderator Models (only one Coefficient Varies within a Model Column) --", "", sep=""))
    cat(sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients and their Standard Errors (SE):")
    cat(" ", "", sep="\n")
    print(modModel2mxAllDRIFT_effects)
    cat(" ", "", sep="\n")
    cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
    currentMatrix <- matrix(NA, 1, (nlatents^2))
    colnames(currentMatrix) <- colnames(modModel2mxAllDRIFT_effects)
    currentMatrix[1,] <- c(round(unlist(modModel2mx_Minus2LogLikelihood), digits))
    rownames(currentMatrix) <- c("-2 Log Likelihoods:                     ")
    print(currentMatrix)
    cat(" ", "", sep="\n")
    currentMatrix[1,] <- c(round(unlist(modModel2mx_estimatedParameters), digits))
    rownames(currentMatrix) <- c("Overall Number of Estimated Parameters: ")
    print(currentMatrix)
    cat(" ", "", sep="\n")
    currentMatrix[1,] <- c(round(unlist(homDRIFTall_modModel2_df), digits))
    rownames(currentMatrix) <- c("Delta(df):                              ")
    print(currentMatrix)
    cat(" ", "", sep="\n")
    currentMatrix[1,] <- c(round(unlist(homDRIFTall_modModel2_Minus2LogLikelihood), digits))
    rownames(currentMatrix) <- c("Delta(-2LL) (= Delta(Chi-square)):      ")
    print(currentMatrix)
    cat(" ", "", sep="\n")
    currentMatrix[1,] <- c(round(unlist(homDRIFTall_modModel2_prob), 2*digits))
    rownames(currentMatrix) <- c("p-value (with double number of digits): ")
    print(currentMatrix)
    cat(" ", "", sep="\n")
    cat("A sign. value (p < .05) indicates that the Drift Coefficients vary among moderator groups):")
    cat(" ", "", sep="\n")
  }

  targetNames <- c()
  for (j in 1:nlatents) {
    for (h in 1:nlatents) {
      targetNames <- c(targetNames, paste0("V",j,"toV", h))
    }
  }

  if (testDRIFTSingleModel == TRUE) {
    for (j in 1:length(targetNames)) {
      cat(" ", "", sep="\n")
      cat(paste("########################################################################################", "", sep="\n"))
      cat(paste("---------- Multi-Sample Homogeneity Model (only ", targetNames[j], " Coefficient Invariant) ----------", "", sep=""))
      cat(sep="\n")
      cat(paste("########################################################################################", "", sep="\n"))
      cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
      cat(" ", "", sep="\n")
      cat("Synthesized/Aggregated Drift Coefficients and their Standard Errors (SE):")
      cat(" ", "", sep="\n")
      print(round(homDRIFTSingleDRIFT_effects[[j]], digits))
      cat(" ", "", sep="\n")
      if (confidenceIntervals == TRUE) {
        cat("Synthesized/Aggregated Drift Coefficients and their 95% Confidence Intervals (Not Necessarily Symmetric):")
        cat(" ", "", sep="\n")
        homDRIFTSingleCI
        toPrint <- homDRIFTSingleCI[[j]][1, 1:3]
        toPrint <- matrix(toPrint, nrow=1, dimnames=list("Fixed Effects", names(toPrint)))
        print(round(toPrint, digits))
      }
      cat(" ", "", sep="\n")
      cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
      cat("-2 Log Likelihood: ", round(homDRIFTSingle_Minus2LogLikelihood[[j]], digits), sep="")
      cat(" ", "", sep="\n")
      cat("Overall Number of Estimated Parameters: ", round(homDRIFTSingle_estimatedParameters[[j]], digits), sep="")
      cat(" ", "", sep="\n")
      cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(homDRIFTSingle_df[[j]], digits), sep="")
      cat(" ", "", sep="\n")
      cat("Comparison with Unconstrained Modell (All Samples Analyzed Separately (~Heterogeneity Model)):")
      cat(" ", "", sep="\n")
      cat("Delta(df): ", round(allStudies_homDRIFTSingle_df[[j]], digits), sep="")
      cat(" ", "", sep="\n")
      cat("Delta(-2LL) (= Delta(Chi-square)): ", round(allStudies_homDRIFTSingle_Minus2LogLikelihood[[j]], digits), sep="")
      cat(" ", "", sep="\n")
      cat("p-value (with double number of digits): ", round(allStudies_homDRIFTSingle_prob[[j]], 2*digits), sep="")
      cat(" ", "", sep="\n")
      cat("A sign. value (p < .05) indicates that the", targetNames[j] , "effects vary among primary studies.")
      cat(" ", "", sep="\n")
    }
  }

  if (testDRIFTCrossModel == TRUE & nlatents > 1) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("------ Multi-Sample Homogeneity Model with all Cross Effects Identical & Invariant -----", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients and their Standard Errors (SE):")
    cat(" ", "", sep="\n")
    print(round(homDRIFTCrossDRIFT_effects, digits))
    cat(" ", "", sep="\n")
    if (confidenceIntervals == TRUE) {
      cat("Synthesized/Aggregated Drift Coefficients and their 95% Confidence Intervals (Not Necessarily Symmetric):")
      cat(" ", "", sep="\n")
      toPrint <- homDRIFTCrossCI[which(rownames(homDRIFTCrossCI) == "AllCrossEffects"), 1:3]
      toPrint <- matrix(toPrint, nrow=1, dimnames=list("Fixed Effects", names(toPrint)))
      print(round(toPrint, digits))
    }
    cat(" ", "", sep="\n")
    cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
    cat("-2 Log Likelihood: ", round(homDRIFTCross_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Overall Number of Estimated Parameters: ", round(homDRIFTCross_estimatedParameters, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(homDRIFTCross_df, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Comparison with Unconstrained Modell (See Note below)")
    cat(" ", "", sep="\n")
    cat("Delta(df): ", round(homDRIFTCross_homDRIFTCross2_df, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Delta(-2LL) (= Delta(Chi-square)): ", round(abs(homDRIFTCross_homDRIFTCross2_Minus2LogLikelihood), digits), sep="")
    cat(" ", "", sep="\n")
    cat("p-value (with double number of digits): ", round(homDRIFTCross_homDRIFTCross2_prob, 2*digits), sep="")
    cat(" ", "", sep="\n")
    cat("Two models were fitted:")
    cat("One model (not reported) with every cross effect invariant across primary studies (but the",
        nlatents^2-nlatents, "cross effects different).")
    cat("A second model (reported above) with all cross effects being equal.")
    cat(" ", "", sep="\n")
    cat("A sign. value (p < .05) indicates that the synthesized/aggregated cross effects are probably NOT equal.")
    cat(" ", "", sep="\n")
  }

  if (testDRIFTAutoModel == TRUE & nlatents > 1) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("------ Multi-Sample Homogeneity Model with all Auto Effects Identical & Invariant ------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Coefficients ------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Synthesized/Aggregated Drift Coefficients and their Standard Errors (SE):")
    cat(" ", "", sep="\n")
    print(round(homDRIFTAutoDRIFT_effects, digits))
    cat(" ", "", sep="\n")
    if (confidenceIntervals == TRUE) {
      cat("Synthesized/Aggregated Drift Coefficients and their 95% Confidence Intervals (Not Necessarily Symmetric):")
      cat(" ", "", sep="\n")
      toPrint <- homDRIFTAutoCI [1, 1:3]
      toPrint <- matrix(toPrint, nrow=1, dimnames=list("Fixed Effects", names(toPrint)))
      print(round(toPrint, digits))
    }
    cat(" ", "", sep="\n")
    cat(paste("------------------------------- Fit Statistics -----------------------------------------", "", sep="\n"))
    cat("-2 Log Likelihood: ", round(homDRIFTAuto_Minus2LogLikelihood, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Overall Number of Estimated Parameters: ", round(homDRIFTAuto_estimatedParameters, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Degrees of Freedom (Standard SEM type. NOT OpenMx/CTSEM type): ", round(homDRIFTAuto_df, digits), sep="")
    cat(" ", "", sep="\n")
    cat("Comparison with Unconstrained Modell (See Note below)")
    cat(" ", "", sep="\n")
    cat("Delta(df): ", round(abs(homDRIFTAuto_homDRIFTAuto2_df), digits), sep="")
    cat(" ", "", sep="\n")
    cat("Delta(-2LL) (= Delta(Chi-square)): ", round(abs(homDRIFTAuto_homDRIFTAuto2_Minus2LogLikelihood), digits), sep="")
    cat(" ", "", sep="\n")
    cat("p-value (with double number of digits): ", round(homDRIFTAuto_homDRIFTAuto2_prob, 2*digits), sep="")
    cat(" ", "", sep="\n")
    cat("Two models were fitted:")
    cat("One model (not reported) with AutoX and AutoY invariant across primary studies (but the", nlatents, "auto effects different).")
    cat("A second model (reported above) with the", nlatents, "auto effects being equal.")
    cat(" ", "", sep="\n")
    cat("A sign. value (p < .05) indicates that the synthesized/aggregated Auto effects are probably NOT equal.")
    cat(" ", "", sep="\n")
  }

  if (exists("UserSpecifiedModelFit")) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------- Summary of user-specified Model. ---------------------------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    print(UserSpecifiedModelFit$results)
    cat(" ", "", sep="\n")
  }

  if (fixedAndRandomEffects == TRUE) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("---- Fixed Effects Analysis of Drift Coefficients (Taken From Heterogeneity Model) -----", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------- Analysis of Fixed And Random Effects ---------------------------", "", sep="\n"))
    cat("-------- i.e., All Drift Effects in All Studies Analyzed Separately (not a CoTiMA) --------")
    cat(" ", "", sep="\n")
    round(allStudiesDRIFT_effects, digits=digits)
    cat(" ", "", sep="\n")
    cat("Fixed Effects Analysis (Borenstein et al., 2007)")
    cat(" ", "", sep="\n")
    fixedResults <-round(fixedEffectDriftResults, digits=digits)
    rownames(fixedResults) <- c("Mean", "Fixed Effect", "Variance of Fixed Effect", "Standard Error of Fixed Effect",
                                "Upper Limit of Fixed Effect", "Lower Limit of Fixed Effect", "Z-Value", "p(Z)",
                                "tau square", "Q", "H square (H2)", "Upper Limit of H2", "Lower Limit of H2",
                                "I square (I2)", "Upper Limit of I2", "Lower Limit of I2")
    colnames(fixedResults) <- c(targetNames)
    print(fixedResults)
    cat(" ", "", sep="\n")
    cat("Random Effects Analysis (Borenstein et al., 2007)")
    randomResults <- round(RandomEffectDriftResults, digits=digits)
    cat(" ", "", sep="\n")
    rownames(randomResults) <- c("Random Effect", "Variance of Random Effect", "Standard Error of Random Effect",
                                 "Upper Limit of Random Effect", "Lower Limit of Random Effect", "Z-Value", "p(Z)",
                                 "Upper Limit of Prediction Interval", "Lower Limit of Prediction Interval")
    colnames(randomResults) <- c(targetNames)
    print(randomResults)
    cat(" ", "", sep="\n")
    cat("Note: If NaN (not available) values appear above, this could indicate negative between study variance (tau-squared).", "", sep="\n")
    cat("\"If the number of studies is very small, then it may be impossible to estimate the ", "", sep="")
    cat("between-studies variance (tau-squared) with any precision. In this case, the fixed effect ", "", sep="")
    cat("model may be the only viable option.\" (Borenstein, Hedges, & Rothstein, 2007, p 114).", "")
  }

  if (publicationBias == TRUE) {
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-- Analysis of Publication Bias (Egger's tests and PET-PEESE conditional estimators) ---", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("-------------------------------- Eggers Test -------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    print(FREAResults)
    cat("There are two established versions of the Egger test. The one above regresses the SNDs of Drift  (instead of Drift) on precision.")
    cat(" ", "", sep="\n")
    cat(paste("--------------------------------- PET-PEESE --------------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    print(round(PET_PEESE_DRIFTresults, digits))
    cat(" ", "", sep="\n")
    cat(paste("Note: ", "", sep=""))
    cat("There are two established versions of the Egger test.")
    cat(" ", "", sep="\n")
    cat("The one above regresses Drift (instead of SNDs of Drift) on precision (or precision squared) .")
    cat(" ", "", sep="\n")
    cat(" ", "", sep="\n")
    cat(paste("-------------------------------- Funnel Plots  ----------------------------------------", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Funnel Plots can be found in ", workingDirectory)
    cat(" ", "", sep="\n")
    cat("The files are labelled :")
    cat(" ", "", sep="\n")
    for (j in 1:length(figureName)) cat(figureName[j], "", sep="\n")
    cat(" ", "", sep="\n")
  }


  if (length(statisticalPower) > 0) {
    tmp <- getOption("max.print"); tmp
    options(max.print = 999999999)
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("------------- Values for Required Sample Sizes for Different Time Lags & Power ---------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    print(round(requiredSampleSizes, digits))
    cat(" ", "", sep="\n")
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("------------- Post Hoc Power of the Primary Studies and their Time Lags Power ----------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat("Calculations are based on the MBESS R Package by Ken Kelley using the ss.power.reg.coef function.", "", sep="\n")
    cat(" ", "", sep="\n")
    cat("R-square values used for sample size calculations were rounded to ", powerRounding, " digits.", "", sep="")
    cat("This substantially reduces (still long) computation times, but may set the required sample size ", "", sep=" ")
    cat("to the default 'largest' value (i.e., 100.000) if they would be large anyway (e.g., > 20.000).", "", sep="\n")
    cat(" ", "", sep="\n")
    for (k in 1:length(postHocPowerList)) {
      cat(" ", "", sep="\n")
      cat("Post Hoc Power for the Effect of ", postHocPowerListNames[k], sep="")
      cat(" ", "", sep="\n")
      print(round(postHocPowerList[[k]], digits))
      cat("The fourth but last row shows the mean power with NAs equated with 0.0 (assuming no power at all).", sep="\n")
      cat("The third row shows the mean power with NAs excluded (assuming estimation problems).", sep="\n")
      cat("The second but last row shows the median power with NAs equated with 0.0 (assuming no power at all).", sep="\n")
      cat("The last row shows the median power with NAs excluded (assuming estimation problems).", sep="\n")
      cat(" ", "", sep="\n")
    }
    cat(" ", "", sep="\n")
    options(max.print=tmp)
  }


  if (plotCrossEffects == TRUE  & nlatents > 1) {

    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("------------------------ Discrete Time Plots of Cross-lagged Effects -------------------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Plots showing the discrete time cross-lagged effects can be found in ", workingDirectory)
    cat(" ", "", sep="\n")
    cat("The files are labelled :")
    cat(" ", "", sep="\n")
    for (j in 1:length(figureFileNameCross)) cat(figureFileNameCross[j], "", sep="\n")
    cat(" ", "", sep="\n")
    cat("Grey and solid lines represent the expected discrete time cross-lagged effects of the primary studies.")
    cat(" ", "", sep="\n")
    cat("The black and dashed line represents the CoTiMA-based expected discrete time cross-lagged effects.")
    cat(" ", "", sep="\n")

    if (testModeratorModel == TRUE) {
      cat(" ", "", sep="\n")
      cat(" ", "", sep="\n")
      cat(" ", "", sep="\n")
      cat("Files with plots showing the moderated discrete time cross-lagged effects are labelled :")
      cat(" ", "", sep="\n")
      for (j in 1:length(figureFileNameCrossModeratedSingle)) cat(figureFileNameCrossModeratedSingle[j], "", sep="\n")
      cat("and ", "", sep="\n")
      for (j in 1:length(figureFileNameCrossModeratedFull)) cat(figureFileNameCrossModeratedFull[j], "", sep="\n")
      cat(" ", "", sep="\n")
    }
  }

  if (plotAutoEffects == TRUE) {

    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("----------------------- Discrete Time Plots of Autoregressive Effects ------------------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Plots showing the discrete time autoregressive effects can be found in ", workingDirectory)
    cat(" ", "", sep="\n")
    cat("The files are labelled :")
    cat(" ", "", sep="\n")
    for (j in 1:length(figureFileNameAuto)) cat(figureFileNameAuto[j], "", sep="\n")
    cat(" ", "", sep="\n")
    cat("Grey and solid lines represent the expected discrete time autoregressive effects of the primary studies.")
    cat(" ", "", sep="\n")
    cat("The black and dashed line represents the CoTiMA-based expected discrete time autoregressive effects.")
    cat(" ", "", sep="\n")

    if (testModeratorModel == TRUE) {
      cat(" ", "", sep="\n")
      cat(" ", "", sep="\n")
      cat(" ", "", sep="\n")
      cat("Files with plots showing the moderated discrete time autoregressive effects are labelled :")
      cat(" ", "", sep="\n")
      for (j in 1:length(figureFileNameAutoModeratedSingle)) cat(figureFileNameAutoModeratedSingle[j], "", sep="\n")
      cat("and ", "", sep="\n")
      for (j in 1:length(figureFileNameAutoModeratedFull)) cat(figureFileNameAutoModeratedFull[j], "", sep="\n")
      cat(" ", "", sep="\n")
    }

  }

  if (length(statisticalPower) > 0) {
    cat(paste("########################################################################################", "", sep="\n"))
    cat(paste("----------------------- Discrete Time Plots of Required Sample Sizes -------------------", "", sep="\n"))
    cat(paste("########################################################################################", "", sep="\n"))
    cat(" ", "", sep="\n")
    cat("Plots showing the required sample sizes to achieve a given statistical power across discrete  be found in ", workingDirectory)
    cat(" ", "", sep="\n")
    cat("The files are labelled :")
    cat(" ", "", sep="\n")
    for (j in 1:length(figureFileNameStatisticalPower)) cat(figureFileNameStatisticalPower[j], "", sep="\n")
    cat(" ", "", sep="\n")
    cat("The plots are based on a model in which all parameters (T0 (co-) variances, drift & diffusion effects ....")
    cat(" ", "", sep="\n")
    cat("... were invariant across primary studies.")
    cat(" ", "", sep="\n")
  }
  if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n", st, "\n", et, "\nAnalysis successfully completed. \nThank you for using CoTiMA.\nHave a nice day!"))}

  sink()
  #cat(green("The warning meassage \'package \'OpenMx\' is required by \'ctsem\', which may no longer work correctly\' could be ignored", "\n"))
  #cat(green("The warning meassage \'In (sqrt(n - 2) : NaNs produced\' could be ignored", "\n"))
  invisible(FitList)
}   # end function definition
