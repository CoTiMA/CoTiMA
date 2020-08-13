# debug <- 0
# if (debug == 1) {
#   singleStudyFits = studyFit
#   ctModelObj = homDRIFTallCTmodelobj
#   fixedModel = homDRIFTallFixedModel
#   modModel = homDRIFTallModModel
#   compSVMethod=compSVMethod
#   moderatorValues = NULL             # moderatorValues short vector with unique moderator values
# }
# debug <- 0

#' ctmaCompSV
#'
#' @param singleStudyFits ?
#' @param ctModelObj ?
#' @param fixedModel ?
#' @param modModel ?
#' @param moderatorValues ?
#' @param compSVMethod ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaCompSV <- function(singleStudyFits = NULL,
                               ctModelObj = NULL,
                               fixedModel = NULL,
                               modModel = NULL,
                               moderatorValues = NULL,             # moderatorValues short vector with unique moderator values
                               compSVMethod="all") {     # "mean", "fixed", "random", "all" = mean of all three
  if (is.null(ctModelObj)){
    cat(crayon::red$bold("ctmodelobj is not specified, which always has to be done!", sep="\n"))
    stop("Good luck for the next try!")
  }
  if (is.null(fixedModel)){
    cat(crayon::red$bold("fixedModel is not specified. Has been set to ctmodelobj (nothing fixed).", sep="\n"))
    fixedModel <- ctModelObj
  }
  if ((!(is.null(modModel))) & (is.null(moderatorValues)) ) {
    cat(crayon::red$bold("modModel is specified but moderatorValues are not provided. Has to be done!", sep="\n"))
    stop("Good luck for the next try!")
  }

  if (!(compSVMethod %in% c("mean", "fixed", "random", "rnw", "all", "lottery") )) {
    cat(crayon::red$bold("Option compSVMethod has to be \"mean\", \"fixed\", \"random\", \"rnw\", or \"all\". Use lower case!", sep="\n"))
    stop("Good luck for the next try!")
  }

  # define required objects
  noOfStudies <- length(singleStudyFits); noOfStudies
  nlatents <- dim(ctModelObj$DRIFT)[1]; nlatents

  # extract all start values & label them
  #singleStudyFits
  est <- lapply(singleStudyFits, function(extract) extract$mxobj$output$estimate); est
  SE <- lapply(singleStudyFits, function(extract) c(t(extract$mxobj$output$standardError))); SE
  for (i in 1:noOfStudies) names(SE[[i]]) <- names(est[[1]])
  sampleSizes <- unlist(lapply(singleStudyFits, function(extract) dim(extract$mxobj$data$observed)[1])); sampleSizes

  # fixed and moderator pattern
  tmp1 <- ctModelObj$DIFFUSION[ctModelObj$DIFFUSION != "0"]; tmp1
  tmp2 <- ctModelObj$T0VAR[ctModelObj$T0VAR != "0"]; tmp2
  labelPattern <- c(ctModelObj$DRIFT, tmp1, tmp2, ctModelObj$CINT); labelPattern

  tmp1 <- fixedModel$DIFFUSION[lower.tri(fixedModel$DIFFUSION, diag=TRUE)]; tmp1
  tmp2 <- fixedModel$T0VAR[lower.tri(fixedModel$T0VAR, diag=TRUE)]; tmp2
  fixedPattern <- c(fixedModel$DRIFT, tmp1, tmp2, fixedModel$CINT); fixedPattern

  if (!(is.null(modModel))) {
    tmp1 <- modModel$DIFFUSION[lower.tri(modModel$DIFFUSION, diag=TRUE)]; tmp1
    tmp2 <- modModel$T0VAR[lower.tri(modModel$T0VAR, diag=TRUE)]; tmp2
    modPattern <- c(modModel$DRIFT, tmp1, tmp2, modModel$CINT); modPattern
  }

  if (!(is.null(modModel))) {
    if (!( (length(labelPattern) == length(fixedPattern)) &
           (length(labelPattern) == length(modPattern)) &
           (length(fixedPattern) == length(modPattern)) )){
      cat(crayon::red$bold("ctmodelobj, fixedModel, and modModel do not macth", sep="\n"))
      stop("Good luck for the next try!")
    }
  } else {
    if (!( (length(labelPattern) == length(fixedPattern)) )){
      cat(crayon::red$bold("ctmodelobj and fixedModel do not macth", sep="\n"))
      stop("Good luck for the next try!")
    }
  }

  ### All Start Values (and SE) in one Vector
  allStartValues <- est; allStartValues
  allStartValuesSE <- SE; allStartValuesSE

  # replace NAs for SE by average SE (happens if single models are under-identified)
  tmp <- apply(apply(simplify2array(allStartValuesSE), 1:2, mean), 1, mean, na.rm=TRUE); tmp
  if (any(is.na(unlist(allStartValuesSE)))) {
    for (i in 1:length(allStartValuesSE)) {
      x <- which(is.na(allStartValuesSE[[i]]))
      allStartValuesSE[[i]][x] <- tmp[x]
    }
  }
  allStartValuesSE

  # eliminate parameters not specified
  labelPattern2 <- labelPattern[labelPattern != "0"]; labelPattern2 # new to eliminate "0" fo CINT

  allStartValues <- lapply(allStartValues, function(x) x[match(labelPattern2, names(x))]); allStartValues
  allStartValuesSE <- lapply(allStartValuesSE, function(x) x[match(labelPattern2, names(x))]); allStartValuesSE

  # assign names
  for (i in 1:noOfStudies) {
    names(allStartValues[[i]]) <- paste0(names(allStartValues[[i]]), "_G", i)
    names(allStartValuesSE[[i]]) <- paste0(names(allStartValuesSE[[i]]), "_G", i)
  }
  # if two or more (sets of) labels are identical replace by mean value and change labels
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

  ### aggregate across all Groups
  aggregatedFixedStartValues <- list()
  #allStartValuesSE
  precision <- 1/(matrix(unlist(allStartValuesSE), nrow=noOfStudies, byrow=TRUE)); precision
  weights <- colSums(precision^2); weights
  coeff <- matrix(unlist(allStartValues), nrow=noOfStudies, byrow=TRUE); coeff

  #allStartValuesSE <- se <- c(.03, .03, .05, .01, .05, .02)^.5; se
  #coeff <- c(.1, .3, .35, .65, .45, .15); coeff
  #precision <- 1/(matrix(se, nrow=6, byrow=TRUE)); precision
  #weights <- (precision^2); weights
  ##weights <- colSums(precision^2); weights
  #CoeffxWeight <- coeff * weights; CoeffxWeight
  #sum(CoeffxWeight)
  #sum(weights)
  #sum(CoeffxWeight)/sum(weights)
  #(CoeffxWeight)/colSums(precision^2)
  #noOfStudies<-3

  # fixed effects computation
  aggregatedFixedStartValues <- colSums( coeff * precision^2 ) / weights; aggregatedFixedStartValues
  names(aggregatedFixedStartValues) <- names(est[[1]][names(est[[1]]) %in% labelPattern])
  aggregatedFixedStartValues <- aggregatedFixedStartValues[labelPattern[!(fixedPattern %in% labelPattern)]]
  if (length(aggregatedFixedStartValues) < 1) aggregatedFixedStartValues <- c()
  aggregatedFixedStartValues

  # NEW random effects
  Q <- colSums(precision^2 * coeff^2) - (colSums(precision^2 * coeff))^2 / colSums(precision^2); Q
  T_weights <- colSums(precision^2); T_weights
  T2_weights <- colSums(precision^2^2); T2_weights
  C <- T_weights-T2_weights/T_weights; C
  tau2 <- (Q -(noOfStudies-1))/C; tau2
  tau2Extended <- do.call(rbind, replicate(noOfStudies, tau2, simplify=FALSE)); tau2Extended
  Ttot_Weights <-colSums(1/ (allStartValuesSE^2 + tau2Extended)); Ttot_Weights
  Ttot_Means <- colSums(coeff * 1/ (allStartValuesSE^2 + tau2Extended)); Ttot_Means
  aggregatedRandomStartValues <- Ttot_Means/Ttot_Weights; aggregatedRandomStartValues
  names(aggregatedRandomStartValues) <- names(est[[1]][names(est[[1]]) %in% labelPattern])
  aggregatedRandomStartValues <- aggregatedRandomStartValues[labelPattern[!(fixedPattern %in% labelPattern)]]
  if (length(aggregatedRandomStartValues) < 1) aggregatedRandomStartValues <- c()
  aggregatedRandomStartValues
  # END NEW

  # NEW means
  aggregatedMeanStartValues <- colMeans(coeff); aggregatedMeanStartValues
  names(aggregatedMeanStartValues) <- names(est[[1]][names(est[[1]]) %in% labelPattern])
  aggregatedMeanStartValues <- aggregatedMeanStartValues[labelPattern[!(fixedPattern %in% labelPattern)]]
  if (length(aggregatedMeanStartValues) < 1) aggregatedMeanStartValues <- c()
  aggregatedMeanStartValues
  # END NEW

  # NEW RootNWeighted RNW
  aggregatedRNWStartValues <- colMeans(coeff*(sampleSizes^.5)) / mean(sampleSizes^.5); aggregatedRNWStartValues
  names(aggregatedRNWStartValues) <- names(est[[1]][names(est[[1]]) %in% labelPattern])
  aggregatedRNWStartValues <- aggregatedRNWStartValues[labelPattern[!(fixedPattern %in% labelPattern)]]
  if (length(aggregatedRNWStartValues) < 1) aggregatedRNWStartValues <- c()
  aggregatedRNWStartValues
  # END NEW


  # Chose starting values
  compSVMethod2 <- compSVMethod
  if (compSVMethod == "lottery") compSVMethod2 <- c("mean", "fixed", "random", "all")[stats::runif(1,1,4)]
  if (compSVMethod2 == "mean") {
    aggregatedStartValues <- aggregatedMeanStartValues
    tmp <- paste0("The method for computing start values is: ", compSVMethod2, " of single study effects" )
  }
  if (compSVMethod2 == "fixed") {
    aggregatedStartValues <- aggregatedFixedStartValues
    tmp <- paste0("The method for computing start values is: ", compSVMethod2, " effects estimates of single study effects" )
  }
  if (compSVMethod2 == "random") {
    aggregatedStartValues <- aggregatedRandomStartValues
    tmp <- paste0("The method for computing start values is: ", compSVMethod2, " effects estimates of single study effects" )
  }
  if (compSVMethod2 == "rnw") {
    aggregatedStartValues <- aggregatedRNWStartValues
    tmp <- paste0("The method for computing start values is: ", "Root Of N-weighted mean", " effects estimates of single study effects" )
  }
  if (compSVMethod2 == "all") {
    aggregatedStartValues <- colMeans(rbind(aggregatedMeanStartValues,
                                            aggregatedFixedStartValues,
                                            aggregatedRandomStartValues,
                                            aggregatedRNWStartValues))
    tmp <- paste0("The method for computing start values is: mean of all availble aggregation methods" )
  }
  print(tmp)

  # replace values by fixed-by-number values
  aggregatedStartValues2 <- aggregatedStartValues
  tmp1 <- suppressWarnings(as.numeric(fixedModel$DRIFT)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$DRIFT[tmp2]; tmp3
  aggregatedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$DRIFT[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$DRIFT
  aggregatedStartValues2

  tmp1 <- suppressWarnings(as.numeric(fixedModel$DIFFUSION)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$DIFFUSION[tmp2]; tmp3
  aggregatedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$DIFFUSION[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$DIFFUSION
  aggregatedStartValues2

  tmp1 <- suppressWarnings(as.numeric(fixedModel$T0VAR)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$T0VAR[tmp2]; tmp3
  aggregatedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$T0VAR[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$T0VAR
  aggregatedStartValues2

  tmp1 <- suppressWarnings(as.numeric(fixedModel$CINT)); tmp1
  tmp2 <- which(!(is.na(tmp1))); tmp2
  tmp3 <- ctModelObj$CINT[tmp2]; tmp3
  aggregatedStartValues2[tmp3] <- tmp1[(!(is.na(tmp1)))]
  ctModelObj$CINT[tmp2] <- tmp1[!(is.na(tmp1))]; ctModelObj$CINT
  aggregatedStartValues2

  aggregatedStartValues <- aggregatedStartValues2[names(aggregatedStartValues2) %in% names(aggregatedStartValues)]
  aggregatedStartValues

  # aggregate across Moderators (fixed effect only - needs to be extended)
  aggregatedModeratedStartValues <- list()
  if (!(is.null(modModel))) {
    for (i in 1:length(unique(unlist(moderatorValues)))) {
      targetModValue <- unique(unlist(moderatorValues))[i]; targetModValue
      tmp1 <- matrix(unlist(allStartValues), nrow=noOfStudies, byrow=TRUE)[which(unlist(moderatorValues) == targetModValue), ]
      tmp2 <- (precision^2)[which(unlist(moderatorValues) == targetModValue), ]
      tmp3 <- dim(tmp1); tmp3
      if (is.null(dim(tmp1))) {
        aggregatedModeratedStartValues[[i]] <- sum( tmp1 * tmp2 ) / weights; aggregatedModeratedStartValues[[i]]
      } else {
        aggregatedModeratedStartValues[[i]] <- colSums( tmp1 * tmp2 ) / weights; aggregatedModeratedStartValues[[i]]
      }
      names(aggregatedModeratedStartValues[[i]]) <- paste0(names(est[[1]][names(est[[1]]) %in% labelPattern]), "_M", i)
      testPattern <- paste0(labelPattern[!(modPattern %in% labelPattern)], "_M", i); testPattern
      aggregatedModeratedStartValues[[i]] <- aggregatedModeratedStartValues[[i]][testPattern]
    }
    aggregatedModeratedStartValues <- unlist(aggregatedModeratedStartValues)
    if (length(aggregatedModeratedStartValues) < 1) aggregatedModeratedStartValues <- c()
  }
  # eliminate moderated effects included in fixed effects
  labelsToDelete <- c()
  if (length(aggregatedStartValues) > 0) {
    for (k in 1:length(aggregatedStartValues)) {
      tmp <- grep(names(aggregatedStartValues)[k] , names(aggregatedModeratedStartValues)); tmp
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
  allStartValues <- unlist(c(aggregatedStartValues, aggregatedModeratedStartValues, allStartValues))

  return(allStartValues)
} # END Function definition

