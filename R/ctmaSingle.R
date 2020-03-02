#######################################################################################################################
################################################ CoTiMA FullDrift #####################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   activeDirectory=NULL
#   single=c("V1toV2", "V2toV1")
#   sourceDirectory= NULL
#   resultsFilePrefix="ctmaSingle"
#   saveFilePrefix="ctmaSingle"
#   ctmaInitFit=CoTiMAInitFit1
#   activateRPB=FALSE
#   checkSingleStudyResults=TRUE
#   type="mx"
#   retryattempts=30
#   refits=1
#   NPSOL=FALSE
#   coresToUse=c(-1)
#   fullDriftStartValues=NULL
#   compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2]
#   useCTMultiGroupAlt=FALSE
#   confidenceIntervals=FALSE
#   saveSingleDriftModelFit=saveFilePrefix
#   ##### ENTER DEBUG INFO HERE #######
#   silentOverwrite=FALSE
#   useCTMultiGroupAlt=TRUE
#   confidenceIntervals=FALSE
#   digits=4
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
#                         cores=1,
#                         #numOfThreads=numOfThreads,
#                         inits=NULL, forcerecompile=FALSE,
#                         savescores=FALSE, gendata=FALSE,
#                         control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
#                         verbose=0)
# }
# debug <- 0

#' ctmaSingle
#'
#' @param activeDirectory ?
#' @param sourceDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param ctmaInitFit ?
#' @param activateRPB ?
#' @param checkSingleStudyResults ?
#' @param silentOverwrite ?
#' @param digits ?
#' @param single ?
#' @param type ?
#' @param retryattempts ?
#' @param refits ?
#' @param NPSOL ?
#' @param coresToUse ?
#' @param fullDriftStartValues ?
#' @param compSVMethod ?
#' @param useCTMultiGroupAlt ?
#' @param confidenceIntervals ?
#' @param saveSingleDriftModelFit ?
#' @param CoTiMAStanctArgs ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaSingle <- function(
  # Directory names and file names
  activeDirectory=NULL,
  sourceDirectory= NULL,
  resultsFilePrefix="ctmaSingle",
  saveFilePrefix="ctmaSingle",

  # Primary Study Fits
  ctmaInitFit=NULL,                    #list of lists: could be more than one fit object

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,
  checkSingleStudyResults=TRUE,
  silentOverwrite=FALSE,

  digits=4,

  # General Model Setup
  single=c(),

  # Fitting Parameters
  type="mx",
  retryattempts=30,
  refits=1,
  NPSOL=FALSE,
  coresToUse=c(-1),
  fullDriftStartValues=NULL,
  compSVMethod=c("mean", "fixed", "random", "rnw", "all")[2],
  useCTMultiGroupAlt=FALSE,
  confidenceIntervals=FALSE,
  saveSingleDriftModelFit=saveFilePrefix,
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
                        cores=1,
                        inits=NULL, forcerecompile=FALSE,
                        savescores=FALSE, gendata=FALSE,
                        control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
                        verbose=0)

)


{  # begin function definition (until end of file)

  # check if fit object is specified
  if (is.null(ctmaInitFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
    stop("Good luck for the next try!")
  }

  if (length(single) < 1) {
    n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    cat(crayon::red$bold("You have not specified the single drift parameters that you want to set free across groups in a series of models!","\n"))
    cat(crayon::red("Would you like to test if ALL single drift coefficients vary across primary studies?","\n"))
    cat(crayon::blue("Press 'y' to test ALL single drift coefficients, 's' to specify the drift coefficients to test,  or 'q' to quit.","\n"))
    char <- readline(" ")
    while (!(char == 'q') & !(char == 'Q') & !(char == 'y') & !(char == 'Y') & !(char == 's') & !(char == 'S')) {
      cat((crayon::blue("Please press 'y' to save ALL, 's' to specify the drift coefficients to test, or 'q' to quit. Press ENTER afterwards.", "\n")))
      char <- readline(" ")
    }
    if (char == 'y' | char == 'Y') {
      single <- c()
      for (i in 1:(n.latent)) {
        for (j in 1:(n.latent)) {
          single <- c(single, paste0("V",i,"toV", j))
        }
      }
    } else if (char == 's' | char == 'S') {
      cat(crayon::blue("Which drift coefficients would you like to test?", "\n"))
      cat(crayon::blue("Please enter them as V1toV1, V1toV2 etc. separated by commas!", "\n"))
      chars <- readline(" ")
      chars <- gsub(" ", "", chars, fixed = TRUE)
      chars <- unlist(strsplit(chars, ","))
      chars
      single <- c(single, chars)
    } else if (char == 'q' | char == 'Q'){
      stop("Good luck for the next try!")
    }
  }

  if (length(saveSingleDriftModelFit) == 1){
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    if (silentOverwrite==FALSE) {
      if ((activateRPB==TRUE) &  (silentOverwrite==FALSE)) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("You have indicated that you want to save singleModelFit, but have not selected which model to save the fit for!","\n"))
      cat(crayon::red("Would you like to save the singleModelFit for specified single drift effects??","\n"))
      cat(crayon::blue("Press 'y' to save ALL singleModelFits or 's' to continue and","\n"))
      cat(crayon::blue("select the singleModelFit you whish to save the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
      char <- readline(" ")
      while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'y' to save ALL, 's' to specify the singleModelFits to save, or 'q' to quit. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'y' | char == 'Y') {
        #for (i in 1:n.studies) {
        for (i in 1:length(single)) {
          counter <- i
          if (length(saveSingleDriftModelFit > 0)){
            #saveSingleDriftModelFit <- c(saveSingleDriftModelFit, counter)
            saveSingleDriftModelFit <- c(saveSingleDriftModelFit, paste0(" ", single[i], " "))
            i <- i+1
          }
        }
      } else if (char == 's' | char == 'S') {
        cat(crayon::blue("Which singleModelFit would you like to save?", "\n"))
        cat(crayon::blue("Please enter the drift effects separated by commas (e.g., V1toV2, V2toV2)!", "\n"))
        chars <- readline(" ")
        chars <- gsub(" ", "", chars, fixed = TRUE)
        chars <- unlist(strsplit(chars, ","))
        chars
        saveSingleDriftModelFit <- c(saveSingleDriftModelFit, chars)
      } else if (char == 'q' | char == 'Q'){
        stop("Good luck for the next try!")
      }
    } else {
      saveSingleDriftModelFit <- c(saveSingleDriftModelFit, single)
    }
    if (silentOverwrite==TRUE) {
      saveSingleDriftModelFit <- c(saveSingleDriftModelFit, single)
    }
  }

  if (resultsFilePrefix=="ctmaSingle") {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat("The default results file prefix (ctmaSingle) has been chosen.", "\n")
    cat("Press 'q' to quit and change or 'c' to continue. Press ENTER afterwards ", "\n")
    char <- readline(" ")
    while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
      cat((crayon::blue("Please press 'q' to quit and change prefix or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
      char <- readline(" ")
    }
    if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
  }



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
    CoTiMAStanctArgs$cores <- coresToUse

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
    #ctmaInitFit <- list(ctmaInitFit1, ctmaInitFit1)
    #ctmaInitFit <- ctmaInitFit1
    #ctmaInitFit$modelResults$DRIFT[[1]]
    n.latent <- length(ctmaInitFit$modelResults$DRIFT[[1]])^.5; n.latent
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    if (is.null(sourceDirectory)) sourceDirectory <- ctmaInitFit$sourceDirectory; sourceDirectory
    n.studies <- unlist(ctmaInitFit$n.studies); n.studies
    #drift_Coef <- ctmaInitFit$studyResults$DRIFT; drift_Coef
    allTpoints <- ctmaInitFit$statisticsList$allTpoints; allTpoints
    maxTpoints <- max(allTpoints); maxTpoints
    allDeltas <- ctmaInitFit$statisticsList$allDeltas; allDeltas
    maxDelta <- max(allDeltas); maxDelta
    manifestNames <- c(ctmaInitFit$studyFitList[[1]]$ctmodelobj$manifestNames); manifestNames

    # derived values
    #allStudies_Minus2LogLikelihood <- sum(unlist(lapply(ctmaInitFit$studyFitList,
    #                                                    function(extract) extract$mxobj$output$Minus2LogLikelihood)))
    #allStudies_Minus2LogLikelihood
    #allStudies_estimatedParameters <- length(unlist(lapply(ctmaInitFit$studyFitList,
    #                                                       function(extract) extract$mxobj$output$estimate)))
    #allStudies_estimatedParameters
    #allStudies_df <- ((n.latent * unlist(allTpoints)) %*% ((n.latent * unlist(allTpoints)) +1 )) / 2 -
    #  allStudies_estimatedParameters; allStudies_df

    usedTimeRange <- seq(0, 1.5*maxDelta, 1)

  }

  #######################################################################################################################
  ############################################ load required subfunctions ###############################################
  #######################################################################################################################

  #{
  #  print(paste0("#################################################################################"))
  #  print(paste0("########################## DEFINE Required Functions ############################"))
  #  print(paste0("#################################################################################"))
  #
  #  source(paste0(sourceDirectory, "CoTiMASaveFile.R"))
  #  source(paste0(sourceDirectory, "CoTiMAcombPRaw"))
  #  source(paste0(sourceDirectory, "CoTiMAcompSV.R"))
  #  source(paste0(sourceDirectory, "CoTiMAchangeSVLab.R"))
  #  source(paste0(sourceDirectory, "CoTiMActMultigroupFitAlt.R"))
  #  source(paste0(sourceDirectory, "CoTiMAbase2Full.R"))
  #  #source(paste0(sourceDirectory, "CoTiMAmx 1.0.0.0.R"))
  #} ### END DEFINE required subfunctions ###


  #######################################################################################################################
  ################################################# data preparation ####################################################
  #######################################################################################################################
  {
    # combine pseudo raw data for mx model
    tmp <- ctmaCombPRaw(listOfStudyFits=ctmaInitFit$studyFitList)
    datawide_all <- tmp$alldata
    groups <- tmp$groups
    names(groups) <- c("Study_No_"); groups
    groupsNamed <- (paste0("Study_No_", groups)); groupsNamed

    # auggment pseudo raw data for stanct model
    if (!(type == "mx" || type == "Mx" || type == "MX"|| type == "mX")) {
      dataTmp <- cbind(datawide_all, groups)
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
      #if (testModeratorModel == TRUE) {
      #  moderatorGroups <- moderators[[1]]; moderatorGroups
      #  for (l in 2:n.studies) moderatorGroups <- rbind(moderatorGroups, moderators[[l]])
      #  moderatorGroups <- moderatorGroups[, moderatorNumber]
      #  tmp <- as.data.frame(moderatorGroups)
      #  if (CoTiMActStanFit$scaleTI == TRUE) tmp <- scale(tmp)
      #  paste0("TI", ((n.studies):(n.studies-1+length(moderatorNumber))))
      #  colnames(tmp) <- paste0("TI", ((n.studies):(n.studies-1+length(moderatorNumber)))); dim(tmp)
      #  if (CoTiMActStanFit$scaleTI == TRUE) tmp <- scale(tmp)
      #  dataTmp <- cbind(dataTmp, tmp); dim(dataTmp)
      #}
      #if (testModeratorModel == TRUE) n.TIpred <- n.studies-1 + length(moderatorNumber) else n.TIpred <- (n.studies-1)
      dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=maxTpoints, n.manifest=n.latent, n.TIpred = (n.studies-1),
                               manifestNames=manifestNames)
      dataTmp3 <- ctsem::ctDeintervalise(dataTmp2)
      dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime
      # eliminate rows where ALL latents are NA
      dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent)], 1, function(x) sum(is.na(x)) != n.latent ), ]
      datalong_all <- dataTmp3
    }

  }

  #######################################################################################################################
  ############################################# CoTiMA (ctsem multigroup) ###############################################
  #######################################################################################################################

  if (type == "mx" || type == "Mx" || type == "MX"|| type == "mX") {
    results <- ctmaSingleMx(
      #type="single",
      single=single,
      datawide_all=datawide_all,
      groups=groups,
      groupsNamed=groupsNamed,
      activeDirectory=activeDirectory,
      sourceDirectory=sourceDirectory,
      resultsFilePrefix="ctmaSingleMx",
      saveFilePrefix="ctmaSingleMx",
      ctmaInitFit=ctmaInitFit,
      activateRPB=activateRPB,
      checkSingleStudyResults=checkSingleStudyResults,
      digits=4,
      retryattempts=retryattempts,
      refits=refits,
      NPSOL=NPSOL,
      coresToUse=coresToUse,
      fullDriftStartValues=fullDriftStartValues,
      compSVMethod=compSVMethod,
      useCTMultiGroupAlt=useCTMultiGroupAlt,
      confidenceIntervals=confidenceIntervals,
      saveSingleDriftModelFit=saveSingleDriftModelFit
    )
  } else {
    results <- ctmaSingleStanct(
      #type="single",
      single=single,
      ctmaInitFit=ctmaInitFit,
      datalong_all=datalong_all,
      activeDirectory=activeDirectory,
      sourceDirectory=sourceDirectory,
      resultsFilePrefix="ctmaSingleStanct",
      saveFilePrefix="ctmaSingleStanct",
      activateRPB=activateRPB,
      silentOverwrite=silentOverwrite,
      digits=digits,
      #n.latent=n.latent,
      #n.studies=n.studies,
      retryattempts=retryattempts,
      refits=refits,
      NPSOL=NPSOL,
      coresToUse=coresToUse,
      numOfThreads=numOfThreads,
      fullDriftStartValues=fullDriftStartValues,
      saveSingleDriftModelFit=saveSingleDriftModelFit,
      CoTiMAStanctArgs=CoTiMAStanctArgs
    )
  }

  #class(results) <- "CoTiMAFit"
  invisible(results)

} ### END function definition

