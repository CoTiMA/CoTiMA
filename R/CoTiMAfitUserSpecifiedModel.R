
debug <- 0
if (debug == 1) {
  singleStudyModelFits = testUserSpecifiedModel$singleStudyModelFits
  moderatorValues = testUserSpecifiedModel$moderatorValues
  listOfFixedModelMatrices = testUserSpecifiedModel$listOfFixedModelMatrices
  listOfModeratorModelMatrices = testUserSpecifiedModel$listOfModeratorModelMatrices
  }
debug <- 0

fitUserSpecifiedModel <- function(singleStudyModelFits = NULL,
                                  moderatorValues = NULL,
                                  userModel = NULL,
                                  listOfFixedModelMatrices = list(T0VAR=NULL, DIFFUSION=NULL, DRIFT=NULL),
                                  listOfModeratorModelMatrices = list(T0VAR=NULL, DIFFUSION=NULL, DRIFT=NULL),
                                  coresToUse=1,
                                  refits=1,
                                  confidenceIntervals=FALSE,
                                  digits=4) 
  
{
  
  if (is.null(singleStudyModelFits)) {
    cat(red$bold("A list of ftted primary study models has to be provided!", sep="\n")) 
    stop("Good luck for the next try!")
  }

  if (is.null(moderatorValues)) {
    cat(red$bold("No moderator values were provided. All were set to 1 (no moderator)!", sep="\n")) 
    stop("Good luck for the next try!")
  }
  
  n.studies <- length(singleStudyModelFits); n.studies
  
  studyFit <- singleStudyModelFits

  # Extract Models
  if (is.null(userModel)) {
    userModel1 <- studyFit[[1]]$ctmodelobj
  } else {
    userModel1 <- userModel
  }
  fixedUserModel1 <- userModel1
  moderatedUserModel1 <- userModel1
  
  # Change user model
  # tonot necessary because effects can be elimited by fixing them to 0 in the fixed Model.
  
  # Change fixed model
  if (!(is.null(listOfFixedModelMatrices$T0VAR))) {
    fixedUserModel1$T0VAR <- listOfFixedModelMatrices$T0VAR; fixedUserModel1$T0VAR
    tmp1 <- suppressWarnings(as.numeric(fixedUserModel1$T0VAR)); tmp1
    tmp2 <- which(!(is.na(tmp1)) & (tmp1 != "0") ); tmp2
    if (length(tmp2) > 0) fixedUserModel1$T0VAR[tmp2] <- tmp1[tmp2]
    dimnames(fixedUserModel1$T0VAR) <- dimnames(studyFit[[1]]$ctmodelobj$T0VAR)
    fixedUserModel1$T0VAR
  }
  if (!(is.null(listOfFixedModelMatrices$DIFFUSION)))  {
    fixedUserModel1$DIFFUSION <- listOfFixedModelMatrices$DIFFUSION
    tmp1 <- suppressWarnings(as.numeric(fixedUserModel1$DIFFUSION)); tmp1
    tmp2 <- which(!(is.na(tmp1)) & (tmp1 != "0") ); tmp2
    if (length(tmp2) > 0) fixedUserModel1$DIFFUSION[tmp2] <- tmp1[tmp2]
    dimnames(fixedUserModel1$DIFFUSION) <- dimnames(studyFit[[1]]$ctmodelobj$fixedUserModel1$DIFFUSION)
    fixedUserModel1$DIFFUSION
  }
  if (!(is.null(listOfFixedModelMatrices$DRIFT))) {
    fixedUserModel1$DRIFT <- listOfFixedModelMatrices$DRIFT
    tmp1 <- suppressWarnings(as.numeric(fixedUserModel1$DRIFT)); tmp1
    tmp2 <- which(!(is.na(tmp1))  ); tmp2
    if (length(tmp2) > 0) fixedUserModel1$DRIFT[tmp2] <- tmp1[tmp2]
    dimnames(fixedUserModel1$DRIFT) <- dimnames(studyFit[[1]]$ctmodelobj$fixedUserModel1$DRIFT)
    fixedUserModel1$DRIFT
  }

  # Change moderated model (fixing values for moderator groups not yet implemented)
  if (!(is.null(listOfModeratorModelMatrices$T0VAR))) moderatedUserModel1$T0VAR <- listOfModeratorModelMatrices$T0VAR
  if (!(is.null(listOfModeratorModelMatrices$DIFFUSION))) moderatedUserModel1$DIFFUSION <- listOfModeratorModelMatrices$DIFFUSION
  if (!(is.null(listOfModeratorModelMatrices$DRIFT))) moderatedUserModel1$DRIFT <- listOfModeratorModelMatrices$DRIFT
  if (!(is.null(listOfModeratorModelMatrices$CINT))) moderatedUserModel1$CINT <- listOfModeratorModelMatrices$CINT
  
  if (!(is.null(moderatedUserModel1))) {
    tmp1 <- which(!(is.na(suppressWarnings(as.numeric(fixedUserModel1$DRIFT))))); tmp1
    tmp2 <- which(moderatedUserModel1$DRIFT == "moderated"); tmp2
    if ( length(tmp1) > 0 & length(tmp2) > 0 ) tmp3 <- match(tmp1, tmp2) else tmp3 <- NA
    if (!(is.na(tmp3))) {
      tmp3 <- tmp2[tmp3]; tmp3
      moderatedUserModel1$DRIFT[tmp3] <- userModel1$DRIFT[tmp3]
    }
    userModel1
    moderatedUserModel1$DRIFT
    
    tmp1 <- which(!(is.na(suppressWarnings(as.numeric(fixedUserModel1$DIFFUSION))))); tmp1
    tmp2 <- which(moderatedUserModel1$DIFFUSION == "moderated"); tmp2
    if ( length(tmp1) > 0 & length(tmp2) > 0 ) tmp3 <- match(tmp1, tmp2) else tmp3 <- NA
    if (!(is.na(tmp3))) {
      tmp3 <- tmp3[-tmp3[is.na(tmp3)]]
      moderatedUserModel1$DIFFUSION[tmp3] <- userModel1$DIFFUSION[tmp3] # "something"
    }
    moderatedUserModel1$DIFFUSION
    
    tmp1 <- which(!(is.na(suppressWarnings(as.numeric(fixedUserModel1$T0VAR))))); tmp1
    tmp2 <- which(moderatedUserModel1$T0VAR == "moderated"); tmp2
    if ( length(tmp1) > 0 & length(tmp2) > 0 ) tmp3 <- match(tmp1, tmp2) else tmp3 <- NA
    if (!(is.na(tmp3))) {
      tmp3 <- tmp3[-tmp3[is.na(tmp3)]]
      moderatedUserModel1$T0VAR[tmp3] <- userModel1$T0VAR[tmp3] # <- "something"
    }
    moderatedUserModel1$T0VAR

    tmp1 <- which(!(is.na(suppressWarnings(as.numeric(fixedUserModel1$CINT))))); tmp1
    tmp2 <- which(moderatedUserModel1$CINT == "moderated"); tmp2
    if ( length(tmp1) > 0 & length(tmp2) > 0 ) tmp3 <- match(tmp1, tmp2) else tmp3 <- NA
    if (!(is.na(tmp3))) {
      tmp3 <- tmp3[-tmp3[is.na(tmp3)]]
      moderatedUserModel1$CINT[tmp3] <- userModel1$CINT[tmp3] # <- "something"
    }
    moderatedUserModel1$CINT
  }
  
  # make data & grouping variabe
  tmp <- combinePseudoRawData(listOfStudyFits=studyFit, moderatorValues=moderatorValues)
  alldata <- tmp$alldata
  groups <- tmp$groups
  moderatorGroups <- tmp$moderatorGroups
  
  # make model adaptations
  tpoints <- length(grep("dT", colnames(alldata)))+1; tpoints 
  userModel1$Tpoints <- tpoints
  fixedUserModel1$Tpoints <- tpoints
  moderatedUserModel1$Tpoints <- tpoints
  
  # Provide start values
  userModel1StartValues <- computeStartValues(singleStudyFits = studyFit,
                                              ctModelObj = userModel1, 
                                              fixedModel = fixedUserModel1, 
                                              modModel = moderatedUserModel1, 
                                              moderatorValues = moderatorValues)
  userModel1StartValues
  
  # Fit
  mxOption(NULL, 'Number of Threads', 1)
  results <- mclapply(seq(1, refits, by=1),
                      function(refits) ctMultigroupFitAlt(dat=alldata, 
                                                          groupings = groups,
                                                          startValues = userModel1StartValues,
                                                          ctmodelobj = userModel1, 
                                                          fixedmodel = fixedUserModel1, 
                                                          modmodel = moderatedUserModel1, 
                                                          moderators = moderatorGroups,
                                                          extraTries=40),
                      mc.cores=coresToUse)
  mxOption(key='Number of Threads', value=parallel::detectCores())
  # Select model with best fit
  allMinus2LogLikelihood <-lapply(results, function(extract) extract$output$Minus2LogLikelihood); allMinus2LogLikelihood
  userModel1Fit <- results[[min(which(allMinus2LogLikelihood==min(unlist(allMinus2LogLikelihood))))]]
  
  # Confidence Intervals
  userModel1FitCI <- list()
  if (confidenceIntervals == TRUE) { 
    tmpModelMxobjFit <- tmpModelMxobj <- ci <- userModel1FitCI <- userModel1CI <- list()
    driftNames <- names(userModel1Fit$output$estimate)[grep("to", names(userModel1Fit$output$estimate))]; driftNames
    for (k in 1:(length(driftNames))) tmpModelMxobjFit[[k]] <- userModel1Fit  # copy mxobj part of fitted models multiple times (for each drift coefficient)
    for (k in 1:(length(driftNames))) ci[[k]] <- mxCI(driftNames[k]) # make mxCI object for every drift coefficients
    for (k in 1:(length(driftNames))) tmpModelMxobj[[k]] <- mxModel(tmpModelMxobjFit[[k]], ci[[k]]) # make OpenMx Models for all drift coefficients
    mxOption(NULL, 'Number of Threads', 1)
    results <- mclapply(seq(1, (length(driftNames)), by=1),
                        function(x) mxRun(tmpModelMxobj[[x]], intervals=TRUE),
                        mc.cores=coresToUse)
    mxOption(key='Number of Threads', value=parallel::detectCores())
    for (k in 1:(length(driftNames))) {
      userModel1FitCI[[k]] <- results[[k]]
      userModel1CI[[k]] <- userModel1FitCI[[k]]$output$confidenceIntervals
    }
  }
  
  # Collect results
  userModel1Est <- userModel1Fit$output$estimate; userModel1Est
  userModel1SE <- userModel1Fit$output$standardErrors; userModel1SE
  userModel1Tvalues <- 1/(userModel1SE/userModel1Est); userModel1Tvalues
  combinedResults <- cbind(userModel1Est, userModel1SE, userModel1Tvalues); combinedResults
  # collect fixed values
  newRows <- matrix(NA, 0, 3); newRows
  # DRIFT
  tmp1 <- which(userModel1Fit$matrices$DRIFT_1$free == FALSE); tmp1
  tmp2 <- userModel1Fit$matrices$DRIFT_1$values[tmp1]; tmp2
  tmp3 <- userModel1Fit$matrices$DRIFT_1$labels[tmp1]; tmp3
  if (!(is.na(tmp3))) {
    newRows <- cbind(tmp2, "-", "fixed value"); newRows
    rownames(newRows) <- tmp3; newRows
  }
  # DIFFUSION
  tmp1 <- which(userModel1Fit$matrices$DIFFUSIONbase_1$free == FALSE); tmp1
  tmp2 <- userModel1Fit$matrices$DIFFUSIONbase_1$values[tmp1]; tmp2
  tmp3 <- userModel1Fit$matrices$DIFFUSIONbase_1$labels[tmp1]; tmp3
  if (!(is.na(tmp3))) {
    newRows2 <- cbind(tmp2, "-", "fixed value"); newRows2
    rownames(newRows2) <- tmp3
    newRows <- rbind(newRows, newRows2)
  }
  # T0VAR
  tmp1 <- which(userModel1Fit$matrices$T0VAR_1$free == FALSE); tmp1
  tmp2 <- userModel1Fit$matrices$T0VAR_1$values[tmp1]; tmp2
  tmp3 <- userModel1Fit$matrices$T0VAR_1$labels[tmp1]; tmp3
  if (!(is.na(tmp3))) {
    newRows2 <- cbind(tmp2, "-", "fixed value"); newRows2
    rownames(newRows2) <- tmp3
    newRows <- rbind(newRows, newRows2)
  }
  # will be combined after CI were added
  
  colnames(combinedResults) <- c("Est.", "SE", "T-value"); combinedResults
  if (confidenceIntervals == TRUE) { 
    combinedResults <- cbind(combinedResults, matrix(NA,dim(combinedResults)[1], 2)); combinedResults
    for (k in 1:(length(driftNames))) combinedResults[k, 4:5] <- userModel1CI[[k]][c(1,3)]
    colnames(combinedResults)[4:5] <- c("lbound", "ubound")
  }
  combinedResults <- round(combinedResults, digits)
  combinedResults
  
  # combine with user-fixed values
  if (dim(newRows)[1] > 0) {
    colnames(newRows) <- c("Est.", "SE", "Note")
    combinedResults <- list(estimatedParameters=combinedResults, fixedParameters=newRows)
  }  
  combinedResults
  
  currentFit <- userModel1Fit$output$Minus2LogLikelihood; currentFit
  names(currentFit) <- "Minus2LogLikelihood"; currentFit
  
  estimatedParameters <- length(userModel1Est); estimatedParameters
  names(estimatedParameters) <- "Number of estimated parameters"; estimatedParameters
  
  combinedResults$Minus2LogLikelihood <- currentFit
  combinedResults$numberOfEstimatedParameters <- estimatedParameters
  
  allResults <- list(userModelFit=userModel1Fit, userModelFitCI=userModel1FitCI, results=combinedResults)
  return(allResults)
}
