#' ctmaPrep
#'
#' @description Combines information of primary studies into a list object and returns this list. This list is then used as input to fit CoTiMA models (ctmaFit, ctmaModFit), for analyses of publication bias (ctmaBias) or statistical power (ctmaPower), or for plotting the results (ctmaPlot).
#'              Primary study information is expected to be assigned to 'numbered' objects. Some of these objects are pre-defined (e.g., 'empcov', 'ageM'). Most of the pre-defined objects could be empty, or they could be dropped by entering their names in the excludedElements-object (e.g., excludedElements = c('ageM')), but dropping them is not really necessary. Additional elements could also be added, which could be useful to put together all information about primary studies at the convenience of the researcher.
#'
#' @param selectedStudies Vector of primary study numbers (numeric values with no leading 0; e.g., '2' but not '02')
#' @param excludedElements Vector of predefined objects used to code primary study information. Some predefined objects are strongly defined; they have to be used in a special way because they are actually used in subsequent analyses. Some other objects could be used at the researcher's convenience (information is just collected). Strongly predefined objects are 'delta_t' (vector of time intervals; the only mandatory requirement; should be of the type c(NA, NA) in cases when raw data are provided), 'sampleSize' (single number), 'pairwiseN' (matrix of pairwise N; could be used if correlation matrix is based on pairwise N), 'empcov' (correlation matrix), 'moderator' (vector of numbers; could be continuous or categorical), 'startValues' (vector of start values), 'rawData' (information about file name and structure of raw data), 'empMeans' (means for variables; usually 0), and 'empVars' (varainces for variables; usually 1). Weakly predefined objects are 'studyNumber' (intended as a special number used for the outputs of subsequently fitted CoTiMA models), 'source' (intended as vector of authors' names and publication year), 'ageM' (intended as
#'                          value indicating the mean age of participants in a primary study), 'malePercent' (intended as value indicating the percentage of male participants in a primary study), 'occupation' (intended as vector of character strings representing the occupations of participants in a primary study), 'country' (intended as single character string representing the country in which a primary study was conducted), 'alphas' (intended as vector of Cronbach's alphas of the variables of a primary study; not yet functional), and 'targetVariables' (intended as vector of character strings representing information about the variables used).'
#' @param addElements User-added objects that are handled as the weakly predefined objects. The major purpose is to collect information a researcher regards as important.
#'
#' @importFrom crayon red
#'
#' @return List of primary studies and parameters for the following CoTiMA
#' @export ctmaPrep
#'
#' @note The following example shows information a researcher has about two studies, which have the numbers '2' and '4'.
#' All information about these studies are stored in objects ending with '2' and '4', respectively.
#' In most instances, one relevant piece of information is the empirical correlation (or covariance) matrix reported in this study,
#' which is stored in the objects 'empcov2' and 'empcov4'. Note that full and symmetric matrices are required for ctmaPrep.
#' Usually, sample sizes ('sampleSize2', 'sampleSize4') and time lags ('delta_t2', 'delta_t4'), are required, too
#'
#' @examples # First Study
#' source2 <- c("Dollard", "& Bakker", "2010")
#' delta_t2 <- 12
#' sampleSize2 <- 209
#' empcov2 <- matrix(c(
#'  1, 0.55, 0.69, 0.37,
#'  0.55, 1, 0.43, 0.55,
#'  0.69, 0.43, 1, 0.58,
#'  0.37, 0.55, 0.58, 1), nrow = 4, ncol = 4)
#' moderator2 <- c(2, 2)
#' addedByResearcher2 <- "something you want to add"
#'
#  Second Study
#' delta_t4 <- c(12, 6)
#' sampleSize4 <- 261
#' empcov4 <- matrix(c(c(1.00, 0.44, 0.74, 0.36, 0.71, 0.32,
#'                       0.44, 1.00, 0.35, 0.66, 0.38, 0.65,
#'                       0.74, 0.35, 1.00, 0.43, 0.83, 0.35,
#'                       0.36, 0.66, 0.43, 1.00, 0.41, 0.71,
#'                       0.71, 0.38, 0.83, 0.41, 1.00, 0.44,
#'                       0.32, 0.65, 0.35, 0.71, 0.44, 1.00),
#'                       nrow=6, ncol=6))
#' moderator4 <- c(3, 1) #
#' addedByResearcher4 <- "another comment"
#' studyList_Ex1 <- ctmaPrep(selectedStudies = c(2, 4),
#'                           excludedElements = "ageM",
#'                           addElements ="addedByResearcher")
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
