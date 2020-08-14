#' ctmaprep
#'
#' @param selectedStudies "Vector of studynumbers to include in CoTiMA"
#' @param excludedElements "Vector of Studyelements to exclude from CoTiMA"
#' @param addElements "Vector of extra Studyelements to include in CoTiMA"
#'
#' @return "List of primary studies and parameters for the following CoTiMA"
#' @export ctmaPrep
#'
#'
ctmaPrep <- function(selectedStudies=NULL,
                     excludedElements=NULL, # vector that might include and of the following (just for cosmetic reasons:
                                            # "deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues",
                                            # "studyNumbers", "rawData", "empMeans", "empVars", "source",
                                            # "ageM", "malePercent", "occupation", "country", "alphas", "targetVariables"
                     addElements=NULL       # vector of names of additiona variables, (e.g. "quality") included in "PREP xxx.R"
                     ) {

    if (is.null(selectedStudies)) {
    cat(crayon::red$bold("Number of primary studies to combine in the list was not specified!", sep="\n"))
    stop("Good luck for the next try!")
  }

  deltas <- sampleSizes <- empcovs <- moderators <- startValues <- studyNumbers <- pairwiseNs <- rawData <- empMeans <- empVars <- source <- list()
  ageM <- malePercent <- occupation <- country <- alphas <- targetVariables <- list()

  if (!(is.null(addElements))) {
    addElementsList <- list()
    for (i in 1:length(addElements)) addElementsList[[addElements[i]]] <- list()
  }

  insideRawData <- list(NULL, NULL, -99, TRUE, FALSE, ".", " ")
  names(insideRawData) <- list("fileName", "studyNumbers", "missingValues", "standardize", "header", "dec", "sep")
  for (i in 1:(length(selectedStudies)+1)) {
    deltas[[i]] <- NA
    sampleSizes[[i]] <- NA
    empcovs[[i]] <- matrix(NA, 0, 0)
    moderators[[i]] <- NA
    startValues[[i]] <- NA
    studyNumbers[[i]] <- selectedStudies[i]
    pairwiseNs[[i]] <- matrix(NA, 0, 0)
    rawData[[i]] <- insideRawData
    empMeans[[i]] <- NA
    empVars[[i]] <- NA
    source[[i]] <- NA
    ageM[[i]] <- NA
    malePercent[[i]] <- NA
    occupation[[i]] <- NA
    country[[i]] <- NA
    alphas[[i]] <- NA
    targetVariables[[i]] <- NA
    if (!(is.null(addElements))) for (j in 1:length(addElements)) addElementsList[[j]][[i]] <- NA
  }

  for (i in 1:length(selectedStudies)) { # 'length' ensures consecutive numbering
    if (exists(paste0("delta_t", selectedStudies[i]))) deltas[[i]] <- get(paste0("delta_t", selectedStudies[i]))
    if (exists(paste0("empcov", selectedStudies[i]))) empcovs[[i]] <- get(paste0("empcov", selectedStudies[i]))
    if (exists(paste0("pairwiseN", selectedStudies[i]))) pairwiseNs[[i]] <- get(paste0("pairwiseN", selectedStudies[i]))
    if (exists(paste0("moderator", selectedStudies[i]))) moderators[[i]] <- get(paste0("moderator", selectedStudies[i]))
    if (exists(paste0("startValues", selectedStudies[i]))) startValues[[i]] <- get(paste0("startValues", selectedStudies[i]))
    if (exists(paste0("studyNumber", selectedStudies[i]))) studyNumbers[[i]] <- get(paste0("studyNumber", selectedStudies[i]))
    if (exists(paste0("sampleSize", selectedStudies[i]))) sampleSizes[[i]] <- get(paste0("sampleSize", selectedStudies[i]))
    if (exists(paste0("empMeans", selectedStudies[i]))) empMeans[[i]] <- get(paste0("empMeans", selectedStudies[i]))
    if (exists(paste0("empVars", selectedStudies[i]))) empVars[[i]] <- get(paste0("empVars", selectedStudies[i]))
    if (exists(paste0("source", selectedStudies[i]))) source[[i]] <- get(paste0("source", selectedStudies[i]))
    if (exists(paste0("ageM", selectedStudies[i]))) ageM[[i]] <- get(paste0("ageM", selectedStudies[i]))
    if (exists(paste0("malePercent", selectedStudies[i]))) malePercent[[i]] <- get(paste0("malePercent", selectedStudies[i]))
    if (exists(paste0("occupation", selectedStudies[i]))) occupation[[i]] <- get(paste0("occupation", selectedStudies[i]))
    if (exists(paste0("country", selectedStudies[i]))) country[[i]] <- get(paste0("country", selectedStudies[i]))
    if (exists(paste0("targetAlphas", selectedStudies[i]))) alphas[[i]] <- get(paste0("targetAlphas", selectedStudies[i]))
    if (exists(paste0("targetVariables", selectedStudies[i]))) targetVariables[[i]] <- get(paste0("targetVariables", selectedStudies[i]))
    if (exists(paste0("rawData", selectedStudies[i]))) {
      rawData[[i]] <- get(paste0("rawData", selectedStudies[i]))
      rawData[[i]]$studyNumbers <- selectedStudies[i]
    }

    if ( (is.na(sampleSizes[[i]]) & (is.null(dim(pairwiseNs[[i]]))) & (is.null(rawData[[i]])) ) ) {
      cat(crayon::red$bold("Neither sample size nor matrix of pairwise N nor rawData was provided for primary study ", i, sep=""))
      cat(" ", sep="\n")
      stop("Good luck for the next try!")
    }

    if (!(is.null(addElements))) {
      for (j in 1:length(addElements)) {
        if (exists(paste0(addElements[j], selectedStudies[i]))) addElementsList[[j]][[i]] <- get(paste0(addElements[[j]], selectedStudies[i]))
      }
    }

  }


  primaryStudies <- list(deltas, sampleSizes, pairwiseNs, empcovs, moderators, startValues,
                         studyNumbers, rawData, empMeans, empVars, source,
                         ageM, malePercent, occupation, country, alphas, targetVariables)

  if (!(is.null(addElements))) {
    for (i in 1:length(addElements)) primaryStudies[[length(primaryStudies)+1]] <- addElementsList[[i]]
  }

  tmpNames <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues",
                "studyNumbers", "rawData", "empMeans", "empVars", "source",
                "ageM", "malePercent", "occupation", "country", "alphas", "targetVariables")
  if (!(is.null(addElements))) {
    for (i in 1:length(addElements)) tmpNames <- c(tmpNames, addElements[i])
  }

  names(primaryStudies) <- tmpNames

  # exclude elements by setting them 0
  if (!(is.null(excludedElements))) {
    targetNames <- c()
    targetNames <- which(excludedElements == names(primaryStudies)); targetNames
    for (i in 1:(length(targetNames))) {
      for (j in 1:length(selectedStudies)) primaryStudies[[targetNames[i]]][j] <- 0
    }
  }

  primaryStudies$n.studies <- length(selectedStudies)

  return(primaryStudies)
}
