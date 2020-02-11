debug <- 0
if (debug == 1) {
  singleStudyFits = studyFit
  ctModelObj = homDRIFTallCTmodelobj
  fixedModel = homDRIFTallFixedModel 
  modModel = homDRIFTallModModel  
  moderatorValues = NULL
  }
debug <- 0

computeStartValues <- function(singleStudyFits = NULL, 
                               ctModelObj = NULL, 
                               fixedModel = NULL, 
                               modModel = NULL, 
                               moderatorValues = NULL) {     # moderatorValues short vector with unique moderator values
  if (is.null(ctModelObj)){
    cat(red$bold("ctmodelobj is not specified, which always has to be done!", sep="\n")) 
    stop("Good luck for the next try!")
  }  
  if (is.null(fixedModel)){
    cat(red$bold("fixedModel is not specified. Has been set to ctmodelobj (nothing fixed).", sep="\n")) 
    fixedModel <- ctModelObj
  }  
  if ((!(is.null(modModel))) & (is.null(moderatorValues)) ) {
    cat(red$bold("modModel is specified but moderatorValues are not provided. Has to be done!", sep="\n")) 
    stop("Good luck for the next try!")
  }  
  
  # define required objects
  noOfStudies <- length(singleStudyFits); noOfStudies
  nlatents <- dim(ctModelObj$DRIFT)[1]; nlatents
  
  # extract all start values & label them
  est <- lapply(singleStudyFits, function(extract) extract$mxobj$output$estimate); est
  SE <- lapply(singleStudyFits, function(extract) c(t(extract$mxobj$output$standardError))); SE
  for (i in 1:noOfStudies) names(SE[[i]]) <- names(est[[1]])
  
  # fixed and moderator pattern
  tmp1 <- ctModelObj$DIFFUSION[ctModelObj$DIFFUSION != "0"]; tmp1
  tmp2 <- ctModelObj$T0VAR[ctModelObj$T0VAR != "0"]; tmp2
  #labelPattern <- c(ctModelObj$DRIFT, ctModelObj$DIFFUSION, ctModelObj$T0VAR, ctModelObj$CINT); labelPattern 
  #labelPattern <- labelPattern[-which(labelPattern == "0")]; labelPattern
  labelPattern <- c(ctModelObj$DRIFT, tmp1, tmp2, ctModelObj$CINT); labelPattern 
  
  tmp1 <- fixedModel$DIFFUSION[lower.tri(fixedModel$DIFFUSION, diag=TRUE)]; tmp1
  tmp2 <- fixedModel$T0VAR[lower.tri(fixedModel$T0VAR, diag=TRUE)]; tmp2
  #fixedPattern <- c(fixedModel$DRIFT, fixedModel$DIFFUSION, fixedModel$T0VAR, fixedModel$CINT); fixedPattern 
  #fixedPattern <- fixedPattern[-which(fixedPattern == "0")]; fixedPattern
  fixedPattern <- c(fixedModel$DRIFT, tmp1, tmp2, fixedModel$CINT); fixedPattern
  
  if (!(is.null(modModel))) {
    tmp1 <- modModel$DIFFUSION[lower.tri(modModel$DIFFUSION, diag=TRUE)]; tmp1
    tmp2 <- modModel$T0VAR[lower.tri(modModel$T0VAR, diag=TRUE)]; tmp2
    #modPattern <- c(modModel$DRIFT, modModel$DIFFUSION, modModel$T0VAR, modModel$CINT); modPattern 
    #modPattern <- modPattern[-which(modPattern == "0")]; modPattern
    modPattern <- c(modModel$DRIFT, tmp1, tmp2, modModel$CINT); modPattern 
  }
  
  if (!(is.null(modModel))) {
    if (!( (length(labelPattern) == length(fixedPattern)) &
           (length(labelPattern) == length(modPattern)) &
           (length(fixedPattern) == length(modPattern)) )){
      cat(red$bold("ctmodelobj, fixedModel, and modModel do not macth", sep="\n")) 
      stop("Good luck for the next try!")
    }  
  } else {
    if (!( (length(labelPattern) == length(fixedPattern)) )){
      cat(red$bold("ctmodelobj and fixedModel do not macth", sep="\n")) 
      stop("Good luck for the next try!")
    }
  }
  
  ### All Start Values (and SE) in one Vector
  allStartValues <- est; allStartValues
  allStartValuesSE <- SE; allStartValuesSE
  
  # eliminate parameters not specified 
  labelPattern2 <- labelPattern[labelPattern != "0"]; labelPattern2 # new to eliminate "0" fo CINT
  
  allStartValues <- lapply(allStartValues, function(x) x[match(labelPattern2, names(x))]); allStartValues
  allStartValuesSE <- lapply(allStartValuesSE, function(x) x[match(labelPattern2, names(x))]); allStartValuesSE
  
  # assign names
  for (i in 1:noOfStudies) {
    names(allStartValues[[i]]) <- paste0(names(allStartValues[[i]]), "_G", i)
    names(allStartValuesSE[[i]]) <- paste0(names(allStartValuesSE[[i]]), "_G", i)
  }
  # if two or more (sets of) lables are identical replace by mean value and change lables
  for (i in which(duplicated(labelPattern2))) {
    toCombine <- which(labelPattern2 == labelPattern2[i]); toCombine
    newStartValues <- lapply(allStartValues, function(combine) mean(combine[toCombine])); newStartValues
    newStartValuesSE <- lapply(allStartValuesSE, function(combine) mean(combine[toCombine])); newStartValuesSE
    for (j in 1:length(allStartValues)) {
      allStartValues[[j]][toCombine] <- newStartValues[[j]]
      allStartValuesSE[[j]][toCombine] <- newStartValuesSE[[j]]
      names(allStartValues[[j]])[toCombine] <- names(allStartValues[[j]])[toCombine[1]]
      names(allStartValuesSE[[j]])[toCombine] <- names(allStartValues[[j]])[toCombine[1]]
    }
  }
  allStartValues <- unlist(allStartValues); allStartValues
  allStartValuesSE <- unlist(allStartValuesSE); allStartValuesSE
  
  # aggregate across all Groups (= fixed effects computation)
  aggregatedFixedStartValues <- list()
  precision <- 1/(matrix(unlist(allStartValuesSE), nrow=noOfStudies, byrow=TRUE)); precision
  weigths <- colSums(precision^2); weigths
  tmp <- matrix(unlist(allStartValues), nrow=noOfStudies, byrow=TRUE) 
  aggregatedFixedStartValues <- colSums( tmp * precision^2 ) / weigths; aggregatedFixedStartValues
  names(aggregatedFixedStartValues) <- names(est[[1]][names(est[[1]]) %in% labelPattern])
  aggregatedFixedStartValues <- aggregatedFixedStartValues[labelPattern[!(fixedPattern %in% labelPattern)]]
  if (length(aggregatedFixedStartValues) < 1) aggregatedFixedStartValues <- c()
  aggregatedFixedStartValues
  # replace values by fixed-by-number values
  aggregatedFixedStartValues2 <- aggregatedFixedStartValues
  tmp1 <- suppressWarnings(as.numeric(fixedModel$DRIFT)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$DRIFT[tmp2]; tmp3
  aggregatedFixedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$DRIFT[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$DRIFT
  aggregatedFixedStartValues2
  
  tmp1 <- suppressWarnings(as.numeric(fixedModel$DIFFUSION)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$DIFFUSION[tmp2]; tmp3
  aggregatedFixedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$DIFFUSION[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$DIFFUSION
  aggregatedFixedStartValues2
  
  tmp1 <- suppressWarnings(as.numeric(fixedModel$T0VAR)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$T0VAR[tmp2]; tmp3
  aggregatedFixedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$T0VAR[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$T0VAR
  aggregatedFixedStartValues2
  
  tmp1 <- suppressWarnings(as.numeric(fixedModel$CINT)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$CINT[tmp2]; tmp3
  aggregatedFixedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$CINT[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$CINT
  aggregatedFixedStartValues2
  
  aggregatedFixedStartValues <- aggregatedFixedStartValues2[names(aggregatedFixedStartValues2) %in% names(aggregatedFixedStartValues)]
  aggregatedFixedStartValues
  
  # aggregate across Moderators
  aggregatedModeratedStartValues <- list()
  if (!(is.null(modModel))) {
    for (i in 1:length(unique(unlist(moderatorValues)))) { 
      targetModValue <- unique(unlist(moderatorValues))[i]; targetModValue
      tmp1 <- matrix(unlist(allStartValues), nrow=noOfStudies, byrow=TRUE)[which(unlist(moderatorValues) == targetModValue), ]
      tmp2 <- (precision^2)[which(unlist(moderatorValues) == targetModValue), ]
      tmp3 <- dim(tmp1); tmp3
      if (is.null(dim(tmp1))) {
        aggregatedModeratedStartValues[[i]] <- sum( tmp1 * tmp2 ) / weigths; aggregatedModeratedStartValues[[i]]
      } else {
        aggregatedModeratedStartValues[[i]] <- colSums( tmp1 * tmp2 ) / weigths; aggregatedModeratedStartValues[[i]]
      }
      names(aggregatedModeratedStartValues[[i]]) <- paste0(names(est[[1]][names(est[[1]]) %in% labelPattern]), "_M", i)
      testPattern <- paste0(labelPattern[!(modPattern %in% labelPattern)], "_M", i); testPattern
      aggregatedModeratedStartValues[[i]] <- aggregatedModeratedStartValues[[i]][testPattern]
    }
    aggregatedModeratedStartValues <- unlist(aggregatedModeratedStartValues)
    if (length(aggregatedModeratedStartValues) < 1) aggregatedModeratedStartValues <- c()
  }
  # eliminate moderated effects inlcuded in fixed effects
  labelsToDelete <- c()
  if (length(aggregatedFixedStartValues) > 0) {
    for (k in 1:length(aggregatedFixedStartValues)) {
      tmp <- grep(names(aggregatedFixedStartValues)[k] , names(aggregatedModeratedStartValues)); tmp
      if (length(tmp) > 0) labelsToDelete <- c(labelsToDelete, tmp)
    }
  }
  if (!(is.null(labelsToDelete))) aggregatedModeratedStartValues <- aggregatedModeratedStartValues[-labelsToDelete]
  aggregatedModeratedStartValues
  
  # remove elements from allStartValues that are fixed across groups or moderators or do not exist 
  labelsToDelete <- c()
  
  if (!(is.null(modModel))) {
    labelsToDelete <- c(labelsToDelete, labelPattern[modPattern == "moderated"]); labelsToDelete
  }
  
  labelsToDelete <- c(labelsToDelete, labelPattern[fixedPattern == "groupfixed"]); labelsToDelete
  
  labelsToDelete <- c(labelsToDelete, labelPattern[suppressWarnings((!(is.na(as.numeric(fixedPattern)))))]); labelsToDelete # remove fixed-by-value entries
  
  labelsToDelete <- unique(labelsToDelete); labelsToDelete
  labelsToDelete <- labelsToDelete[labelsToDelete != "0"]; labelsToDelete
  
  if (length(labelsToDelete) > 0 ) {
    for (i in 1:length(labelsToDelete)) allStartValues <- allStartValues[-(grep(labelsToDelete[i], names(allStartValues)))]
  }
  # combine all start values
  allStartValues <- unlist(c(aggregatedFixedStartValues, aggregatedModeratedStartValues, allStartValues))
  
  return(allStartValues)
} # END Function definition

