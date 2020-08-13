#' ctmaPrep
#'
#' @param selectedStudies ?
#' @param excludedElements ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaPrep <- function(selectedStudies=NULL, excludedElements=NULL) {
  #selectedStudies=1:5
  #excludedElements="moderators"
  #excludedElements=NULL

  if (is.null(selectedStudies)) {
    cat(crayon::red$bold("Number of primary studies to combine in the list was not specified!", sep="\n"))
    stop("Good luck for the next try!")
  }

  deltas <- sampleSizes <- empcovs <- moderators <- startValues <- studyNumbers <- pairwiseNs <- rawData <- empMeans <- empVars <- list()
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
    if (exists(paste0("rawData", selectedStudies[i]))) {
      rawData[[i]] <- get(paste0("rawData", selectedStudies[i]))
      rawData[[i]]$studyNumbers <- selectedStudies[i]
    }

    if ( (is.na(sampleSizes[[i]]) & (is.null(dim(pairwiseNs[[i]]))) & (is.null(rawData[[i]])) ) ) {
      cat(crayon::red$bold("Neither sample size nor matrix of pairwise N nor rawData was provided for primary study ", i, sep=""))
      cat(" ", sep="\n")
      stop("Good luck for the next try!")
    }
  }

  primaryStudies <- list(deltas, sampleSizes, pairwiseNs, empcovs, moderators, startValues,
                         studyNumbers, rawData, empMeans, empVars)
  names(primaryStudies) <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues",
                             "studyNumbers", "rawData", "empMeans", "empVars")

  # exclude elements by setting them 0
  if (!(is.null(excludedElements))) {
    targetNames <- c()
    targetNames <- which(excludedElements == names(primaryStudies)); targetNames
    for (i in 1:(length(targetNames))) {
      for (j in 1:length(selectedStudies)) primaryStudies[[targetNames[i]]][j] <- 0
    }
  }

  return(primaryStudies)
}
