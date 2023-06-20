#' ctmaPrep
#'
#' @description Combines information of primary studies into a list object and returns this list. This list is then used as input to
#' fit 'ctsem' models. Primary study information is expected to be assigned to 'numbered' objects. Some of these objects are pre-defined
#' (e.g., 'empcov', 'ageM'). Most of the pre-defined objects could be empty, or they could be dropped by entering their names in the
#' excludedElements-object (e.g., excludedElements = c('ageM')), but dropping them is not really necessary. Additional elements could
#' also be added, which could be useful to put together all information about primary studies at the convenience of the researcher.
#'
#' @param selectedStudies Vector of primary study numbers (numeric values with no leading 0; e.g., '2' but not '02')
#' @param excludedElements Vector of predefined objects used to code primary study information. Some predefined objects are strongly
#' defined; they have to be used in a special way because they are actually used in subsequent analyses. Some other objects could be
#' used at the researcher's convenience (information is just collected). Strongly predefined objects are 'delta_t' (vector of time
#' intervals; the only mandatory requirement; should be of the type c(NA, NA) in cases when raw data are provided), 'sampleSize'
#' (single number), 'pairwiseN' (matrix of pairwise N; could be used if correlation matrix is based on pairwise N), 'empcov' (correlation
#' matrix), 'moderator' (vector of numbers; could be continuous or categorical), 'startValues' (vector of start values), 'rawData'
#' (information about file name and structure of raw data), 'empMeans' (means for variables; usually 0), and 'empVars' (varainces for
#' variables; usually 1). Weakly predefined objects are 'studyNumber' (intended as a special number used for the outputs of subsequently
#' fitted CoTiMA models), 'source' (intended as vector of authors' names and publication year), 'ageM' (intended as value indicating the
#' mean age of participants in a primary study), 'malePercent' (intended as value indicating the percentage of male participants in a
#' primary study), 'occupation' (intended as vector of character strings representing the occupations of participants in a primary study),
#' 'country' (intended as single character string representing the country in which a primary study was conducted), 'alphas' (intended as
#' vector of Cronbach's alphas of the variables of a primary study; not yet functional), and 'targetVariables' (intended as vector of
#' character strings representing information about the variables used).'
#' @param addElements User-added objects that are handled as the weakly predefined objects. The major purpose is to collect information
#' a researcher regards as important.
#' @param digits Rounding used for summary function
#' @param moderatorLabels character vector of names
#' @param moderatorValues list of character vectors
#' @param summary if TRUE (default) creates summary table and xlsx sheets. Could be set to FALSE in case of errors.
#' @param activeDirectory Mandatory. If subsequent fitting is done using different folders or on different computers, it can be
#' changed so that raw data files can be loaded.
#'
#' @importFrom crayon red
#' @importFrom openxlsx addWorksheet writeData createWorkbook openXL saveWorkbook
#'
#' @return List of primary studies and parameters for the following CoTiMA (plus StudyInformation which could be saved to Excel)
#'
#' @export ctmaPrep
#'
#' @note The following example shows information a researcher has about three studies, which have the numbers '2', '4' and '17'.
#' All information about these studies are stored in objects ending with '2', '4', and '17', respectively. In most instances, one
#' relevant piece of information is the empirical correlation (or covariance) matrix reported in this study, which is stored in the
#' objects 'empcov2', 'empcov4', and 'empcov17'. Note that full and symmetric matrices are required for ctmaPrep. Usually, sample
#' sizes ('sampleSize2', 'sampleSize4', & 'sampleSize17') and time lags ('delta_t2', 'delta_t4', & 'delta_t17'), are required, too.
#'
#' @examples
#' # First Study
#' empcov2 <- matrix(c(1.00, 0.45, 0.57, 0.18,
#'                     0.45, 1.00, 0.31, 0.66,
#'                     0.57, 0.31, 1.00, 0.40,
#'                     0.18, 0.66, 0.40, 1.00), nrow=4, ncol=4)
#' delta_t2 <- 12
#' sampleSize2 <- 148
#' moderator2 <- c(1, 0.72)
#' source2 <- c("Houkes, I,", "Janssen, P, P, M,", "de Jonge, J",
#'               "& Bakker, A, B", "Study1", "2003")
#' addedByResearcher2 <- "something you want to add"
#'
#' # Second Study
#' empcov3 <- matrix(c(1.00, 0.43, 0.71, 0.37,
#'                     0.43, 1.00, 0.34, 0.69,
#'                     0.71, 0.34, 1.00, 0.50,
#'                     0.37, 0.69, 0.50, 1.00), nrow=4, ncol=4)
#' delta_t3 <- 12
#' sampleSize3 <- 88
#' moderator3 <- c(1, 0.72)
#' source3 <- c("Houkes, I,", "Janssen, P, P, M,", "de Jonge, J",
#'               "& Bakker, A, B", "Study2", "2003")
#' addedByResearcher3 <- ""
#'
#' # Third Study
#' empcov313 <- matrix(c(1.00, 0.38, 0.54, 0.34, 0.60, 0.28,
#'                       0.38, 1.00, 0.34, 0.68, 0.28, 0.68,
#'                       0.54, 0.34, 1.00, 0.47, 0.66, 0.39,
#'                       0.34, 0.68, 0.47, 1.00, 0.38, 0.72,
#'                       0.60, 0.28, 0.66, 0.38, 1.00, 0.38,
#'                       0.28, 0.68, 0.39, 0.72, 0.38, 1.00), nrow=6, ncol=6)
#' delta_t313 <- c(1.5, 1.5)
#' sampleSize313 <- 335
#' moderator313 <- c(0.8,	2.47)
#' source313 <- c("Demerouti", "Bakker", "& Bulters", "2004")
#' addedByResearcher313 <- "check correlation matrix"
#'
#' # Add Labels and Values for Moderators (just for optional excel tables)
#' moderatorLabels <- c("Control", "Social Support")
#' moderatorValues <- list("continuous", c("1 = very low", "2 = low",
#'                        "3 = medium", "4 = high", "5 = very high"))
#'
#' CoTiMAstudyList_3 <- ctmaPrep(selectedStudies = c(2, 3, 313),
#'                               activeDirectory="/user/",
#'                               excludedElements = "ageM",
#'                               addElements = "addedByResearcher",
#'                               moderatorLabels=moderatorLabels,
#'                               moderatorValues=moderatorValues)
#'
ctmaPrep <- function(selectedStudies=NULL,
                     excludedElements=NULL,
                     addElements=NULL,
                     digits=4,
                     moderatorLabels=NULL,
                     moderatorValues=NULL,
                     summary=TRUE,
                     activeDirectory=NULL
) {

  ctma <- globalenv()

  if (is.null(selectedStudies)) {
    ErrorMsg <- "Number of primary studies to combine in the list was not specified! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  deltas <- sampleSizes <- empcovs <- moderators <- startValues <- studyNumbers <- pairwiseNs <- rawData <- empMeans <- empVars <- source <- list()
  ageM <- malePercent <- occupation <- country <- alphas <- targetVariables <- list()
  recodeVariables <- combineVariables <- combineVariablesNames <- missingVariables  <- list()

  if (!(is.null(addElements))) {
    addElementsList <- list()
    for (i in 1:length(addElements)) addElementsList[[addElements[i]]] <- list()
  }

  if (is.null(activeDirectory)) {
    ErrorMsg <- "\nNo active directory has been specified! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  #insideRawData <- list(NULL, NULL, -99, TRUE, FALSE, ".", " ")
  #names(insideRawData) <- list("fileName", "studyNumbers", "missingValues", "standardize", "header", "dec", "sep")
  insideRawData <- list(fileName=NULL, studyNumbers=NULL, missingValues=-99,
                        standardize=TRUE, header=FALSE, dec =".", sep=" ",
                        n.ind.mod=0)

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
    recodeVariables[[i]] <- NA
    combineVariables[[i]] <- list()
    combineVariablesNames[[i]] <- NA
    missingVariables[[i]] <- NA

    if (!(is.null(addElements))) for (j in 1:length(addElements)) addElementsList[[j]][[i]] <- NA
  }

  for (i in 1:length(selectedStudies)) { # 'length' ensures consecutive numbering
    if (exists(paste0("delta_t", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) deltas[[i]] <- get(paste0("delta_t", selectedStudies[i]))
    if (exists(paste0("empcov", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) empcovs[[i]] <- get(paste0("empcov", selectedStudies[i]))
    if (exists(paste0("pairwiseN", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) pairwiseNs[[i]] <- get(paste0("pairwiseN", selectedStudies[i]))
    if (exists(paste0("moderator", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) moderators[[i]] <- get(paste0("moderator", selectedStudies[i]))
    if (exists(paste0("startValues", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) startValues[[i]] <- get(paste0("startValues", selectedStudies[i]))
    if (exists(paste0("studyNumber", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) studyNumbers[[i]] <- get(paste0("studyNumber", selectedStudies[i]))
    #if (exists(paste0("sampleSize", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) sampleSizes[[i]] <- get(paste0("sampleSize", selectedStudies[i]))

    if (exists(paste0("sampleSize", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      tmp1 <- get(paste0("sampleSize", selectedStudies[i]), envir =parent.frame(), inherits=FALSE); tmp1
      if (!(is.null(tmp1))) sampleSizes[[i]] <- get(paste0("sampleSize", selectedStudies[i]),
                                                    envir =parent.frame(), inherits=FALSE) else sampleSizes[[i]] <- NA
    } else {
      sampleSizes[[i]] <- NA
    }


    if (exists(paste0("empMeans", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) empMeans[[i]] <- get(paste0("empMeans", selectedStudies[i]))
    if (exists(paste0("empVars", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) empVars[[i]] <- get(paste0("empVars", selectedStudies[i]))
    if (exists(paste0("source", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) source[[i]] <- get(paste0("source", selectedStudies[i]))
    if (exists(paste0("ageM", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) ageM[[i]] <- get(paste0("ageM", selectedStudies[i]))
    if (exists(paste0("malePercent", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) malePercent[[i]] <- get(paste0("malePercent", selectedStudies[i]))
    if (exists(paste0("occupation", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) occupation[[i]] <- get(paste0("occupation", selectedStudies[i]))
    if (exists(paste0("country", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) country[[i]] <- get(paste0("country", selectedStudies[i]))
    if (exists(paste0("alphas", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) alphas[[i]] <- get(paste0("alphas", selectedStudies[i]))

    if (exists(paste0("targetVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      tmp1 <- get(paste0("targetVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE); tmp1
      if (!(is.null(tmp1))) targetVariables[[i]] <- get(paste0("targetVariables", selectedStudies[i]),
                                                        envir =parent.frame(), inherits=FALSE) else targetVariables[[i]] <- NA
    } else {
      targetVariables[[i]] <- NA
    }

    if (exists(paste0("recodeVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      tmp1 <- get(paste0("recodeVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE); tmp1
      if (!(is.null(tmp1))) recodeVariables[[i]] <- get(paste0("recodeVariables", selectedStudies[i]),
                                                        envir =parent.frame(), inherits=FALSE) else recodeVariables[[i]] <- NA
    } else {
      recodeVariables[[i]] <- NA
    }

    if (exists(paste0("combineVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      tmp1 <- get(paste0("combineVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE); tmp1
      if (length(tmp1) > 0) {
        tmp2 <- c()
        for (l in 1:length(tmp1)) {
          tmp3 <- c()
          for (m in 1:length(tmp1[[l]])) {
            tmp3 <- paste(tmp3, tmp1[[l]][m], sep=" + "); tmp3
          }
          tmp2[l] <- substring(tmp3, 4)
        }
        combineVariables[[i]] <- tmp2
      } else {
        combineVariables[[i]] <- NA
      }
    } else {
      combineVariables[[i]] <- NA
    }

    if (exists(paste0("combineVariablesNames", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      tmp1 <- get(paste0("combineVariablesNames", selectedStudies[i])); tmp1
      if (!(is.null(tmp1))) combineVariablesNames[[i]] <- get(paste0("combineVariablesNames", selectedStudies[i]),
                                                              envir =parent.frame(), inherits=FALSE) else combineVariablesNames[[i]] <- NA
    } else {
      combineVariablesNames[[i]] <- NA
    }

    if (exists(paste0("missingVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      tmp1 <- get(paste0("missingVariables", selectedStudies[i]), envir =parent.frame(), inherits=FALSE); tmp1
      if (!(is.null(tmp1))) missingVariables[[i]] <- get(paste0("missingVariables", selectedStudies[i]),
                                                         envir =parent.frame(), inherits=FALSE) else missingVariables[[i]] <- NA
    } else {
      missingVariables[[i]] <- NA
    }

    if (exists(paste0("rawData", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)) {
      rawData[[i]] <- get(paste0("rawData", selectedStudies[i]), envir =parent.frame(), inherits=FALSE)
      rawData[[i]]$studyNumbers <- selectedStudies[i]
      # CHD 19.6.2023 to allow separating ind level moderators from rest of data file in ctmaInit
      if (is.null(rawData[[i]]$n.ind.mod)) rawData[[i]]$n.ind.mod <- 0
    }
    if ( (is.na(sampleSizes[[i]]) & (is.null(dim(pairwiseNs[[i]]))) & (is.null(rawData[[i]])) ) ) {
      cat(crayon::red$bold("Neither sample size nor matrix of pairwise N nor rawData was provided for primary study ", i, sep=""))
      ErrorMsg <- "Good luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.null(addElements))) {
      for (j in 1:length(addElements)) {
        if (exists(paste0(addElements[j], selectedStudies[i]),
                   envir =parent.frame(), inherits=FALSE)) addElementsList[[j]][[i]] <- get(paste0(addElements[[j]], selectedStudies[i]),
                                                                                            envir =parent.frame(), inherits=FALSE)
      }
    }
  }

  primaryStudies <- list(deltas, sampleSizes, pairwiseNs, empcovs, moderators, startValues,
                         studyNumbers, rawData, empMeans, empVars, source,
                         ageM, malePercent, occupation, country, alphas, targetVariables,
                         recodeVariables, combineVariables, combineVariablesNames, missingVariables)

  if (!(is.null(addElements))) {
    for (i in 1:length(addElements)) primaryStudies[[length(primaryStudies)+1]] <- addElementsList[[i]]
  }

  tmpNames <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues",
                "studyNumbers", "rawData", "empMeans", "empVars", "source",
                "ageM", "malePercent", "occupation", "country", "alphas", "targetVariables",
                "recodeVariables", "combineVariables", "combineVariablesNames", "missingVariables")
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

  if (summary == TRUE) {
    # create summary
    # values required for printing matrix values in a single row
    primaryStudies2 <- primaryStudies
    n.studies  <- primaryStudies$n.studies
    maxWaves <- max(unlist(lapply(primaryStudies$deltas, length)))+1; maxWaves
    maxEmpcov <- max(unlist(lapply(primaryStudies$empcovs, length)))^.5; maxEmpcov
    #maxPairwiseNs <- max(unlist(lapply(primaryStudies$pairwiseNs, length)))^.5; maxPairwiseNs
    n.variables <- maxEmpcov/maxWaves; n.variables

    studyListCategories <- vector("list", length=length(names(primaryStudies2))); studyListCategories
    names(studyListCategories) <- names(primaryStudies2); studyListCategories
    studyListCategories$n.studies <- NULL # do not summarize n.studies (is constant)
    primaryStudies2$n.studies <- NULL #
    studyListCategories$startValues <- NULL # do not summarize start values
    primaryStudies2$startValues <- NULL
    studyListCategories$rawData <- NULL # do not summarize raw data
    primaryStudies2$rawData <- NULL #
    studyListCategories$empMeans <- NULL # do not summarize means
    primaryStudies2$empMeans <- NULL
    studyListCategories$empVars <- NULL # do not summarize variances
    primaryStudies2$empVars <- NULL

    summaryTable <- matrix(NA, nrow=n.studies, ncol=0); summaryTable
    for (i in 1:length(studyListCategories)) {
      if (any(!(is.na(primaryStudies2[[i]])))) {
        # check max length of list elements across studies
        maxLength <- max(unlist(lapply(primaryStudies2[[i]], length))); maxLength
        object <- "vector"
        if (names(studyListCategories)[i] %in% c("empcovs", "pairwiseNs")) object <- "matrix"
        if (names(studyListCategories)[i] %in% c("combineVariables")) object <- "list"
        if (object == "matrix") {
          maxLength <- maxLength ^.5; maxLength # correction if input is matrix
          maxDim <- maxLength # CHD 2.2.23
          maxLength <- maxLength * (maxLength - 1) / 2; maxLength
        }
        if (names(studyListCategories)[i] %in% c("alphas")) maxLength <- maxWaves * n.variables

        if (maxLength > 0) {
          tmpTable <- matrix(NA, nrow=n.studies, ncol=maxLength); tmpTable
          for (j in (1:n.studies)) {
            #j <- 1
            if (length(primaryStudies2[[i]][[j]]) > 0) {
              if (object == "matrix") {
                # currentLength <- length(primaryStudies2[[i]][[j]])^.5; currentLength
                # currentLength <- currentLength * (currentLength-1) / 2; currentLength
                # for (k in 1:currentLength) tmpTable[j, k] <- round(primaryStudies2[[i]][[j]][lower.tri(primaryStudies2[[i]][[j]])][k], digits)
                # new CHD 2.2.23
                tmp1 <- primaryStudies2[[i]][[j]]
                currentDim <- length(tmp1)^.5; currentDim
                tmp2 <- maxDim - currentDim; tmp2
                if (tmp2 > 0) {
                  tmp1 <- cbind(tmp1, matrix(NA, nrow=currentDim, ncol=tmp2))
                  tmp1 <- rbind(tmp1, matrix(NA, nrow=tmp2, ncol=maxDim))
                }
                currentLength <- length(tmp1)^.5; currentLength
                currentLength <- currentLength * (currentLength-1) / 2; currentLength
                tmp1[lower.tri(tmp1)]
                for (k in 1:currentLength) if (!(is.na(tmp1[lower.tri(tmp1)][k]))) tmpTable[j, k] <- round(tmp1[lower.tri(tmp1)][k], digits)
              }
              if (object == "vector") {
                for (k in 1:maxLength) tmpTable[j, k] <- primaryStudies2[[i]][[j]][k]
              }
              if (object == "list") {
                tmp1 <- primaryStudies2[[i]][[j]]; tmp1
                if (length(tmp1) > 0) {
                  tmp2 <- c()
                  for (l in 1:length(tmp1)) {
                    tmp3 <- c()
                    for (m in 1:length(tmp1[[l]])) {
                      tmp3 <- paste0(tmp3, tmp1[[l]][m]); tmp3
                    }
                    tmp2[l] <- tmp3
                  }
                  for (m in 1:maxLength) tmpTable[j, m] <- tmp2[m]
                }
              }
            }
          } # end for (j in (1:n.studies))
          #tmpTable

          if (names(studyListCategories)[i] %in% c("ageM", "ageSD", "malePercent")) tmpTable <- round(tmpTable, digits)
          tmpTableNames <- tmpTableNamesBackup <- gsub("$", "", names(studyListCategories[i])); tmpTableNames
          if (tmpTableNamesBackup == "deltas") tmpTableNames <- paste0("Delta", " Lag ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "moderators") tmpTableNames <- paste0("Moderator", " # ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "sampleSizes") tmpTableNames <- "N"; tmpTableNames
          if (tmpTableNamesBackup == "pairwiseNs") tmpTableNames <- "pairwise N"; tmpTableNames
          if (tmpTableNamesBackup == "studyNumbers") tmpTableNames <- "Orig. Study No."; tmpTableNames
          if (tmpTableNamesBackup == "source") tmpTableNames <- paste0("Source Info ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "occupation") tmpTableNames <- paste0("Occupation ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "targetVariables") tmpTableNames <- paste0("Variable ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "country") tmpTableNames <- paste0("Country ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "recodeVariables") tmpTableNames <- paste0("recodeVariables Nr.", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "combineVariables") tmpTableNames <- paste0("combineVariables ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "combineVariablesNames") tmpTableNames <- paste0("combineVariablesNames ", seq(1, maxLength, 1)); tmpTableNames
          if (tmpTableNamesBackup == "missingVariables") tmpTableNames <- paste0("missingVariables", seq(1, maxLength, 1)); tmpTableNames

          if (tmpTableNamesBackup == "alphas") {
            tmpTableNames <- c()
            for (k in 1:maxWaves) {
              for (l in 1:n.variables) {
                tmpTableNames <- c(tmpTableNames, paste0("alpha Y", l, "_T", k-1)); tmpTableNames
              }
            }
          }

          if (tmpTableNamesBackup == "empcovs") {
            tmpTableNames <- c()
            for (k in 1:maxWaves) {
              for (l in 1:n.variables) {
                for (m in 1:maxWaves) {
                  for (n in 1:n.variables) {
                    tmpTableNames <- c(tmpTableNames, paste0("r(Y", l, "_T", k-1, ") (Y", n, "_T", m-1, ")"))
                  }
                }
              }
            }
            tmpTableNamesMat <- matrix(tmpTableNames, n.variables*maxWaves, n.variables*maxWaves); tmpTableNamesMat
            tmpTableNames <- tmpTableNamesMat[lower.tri(tmpTableNamesMat)]; tmpTableNames
            tmpTableNames <- tmpTableNames[1:maxLength] # test
          }

          if (tmpTableNamesBackup == "pairwiseNs") {
            tmpTableNames <- c()
            for (k in 1:maxWaves) {
              for (l in 1:n.variables) {
                for (m in 1:maxWaves) {
                  for (n in 1:n.variables) {
                    tmpTableNames <- c(tmpTableNames, paste0("N(Y", l, "_T", k-1, ") (Y", n, "_T", m-1, ")"))
                  }
                }
              }
            }
            tmpTableNamesMat <- matrix(tmpTableNames, n.variables*maxWaves, n.variables*maxWaves); tmpTableNamesMat
            tmpTableNames <- tmpTableNamesMat[lower.tri(tmpTableNamesMat)]; tmpTableNames
            tmpTableNames <- tmpTableNames[1:maxLength] # test
          }

          # if less columnames than columns are present
          if (length(tmpTableNames) < dim(tmpTable)[2]) tmpTableNames <- rep(tmpTableNames[1], dim(tmpTable)[2])
          colnames(tmpTable) <- tmpTableNames

          if (tmpTableNamesBackup == "source") summaryTable <- cbind(tmpTable,summaryTable) else summaryTable <- cbind(summaryTable, tmpTable)
        } # end if (maxLength > 0)
      }
    }

    primaryStudies$summary <- as.data.frame(summaryTable)

    moderatorLabelsBackup <- moderatorLabels; moderatorLabelsBackup
    moderatorValuesBackup <- moderatorValues; moderatorValuesBackup
    if (is.null(moderatorLabels)) moderatorLabels <- NA
    if (is.null(moderatorValues)) moderatorValues <- NA
    primaryStudies$moderatorLabels <- moderatorLabels; primaryStudies$moderatorLabels
    primaryStudies$moderatorValues <- moderatorValues; primaryStudies$moderatorValues


    ### prepare Excel Workbook with several sheets ################################################################
    wb <- openxlsx::createWorkbook()
    sheet1 <- openxlsx::addWorksheet(wb, sheetName="All Primary Study Information")
    sheet2 <- openxlsx::addWorksheet(wb, sheetName="Deltas")
    sheet3 <- openxlsx::addWorksheet(wb, sheetName="Sample Sizes")
    sheet4 <- openxlsx::addWorksheet(wb, sheetName="Correlations")
    sheet5 <- openxlsx::addWorksheet(wb, sheetName="Moderators")
    sheet6 <- openxlsx::addWorksheet(wb, sheetName="Countries")
    sheet7 <- openxlsx::addWorksheet(wb, sheetName="Occupations")
    sheet8 <- openxlsx::addWorksheet(wb, sheetName="Demographics")
    sheet9 <- openxlsx::addWorksheet(wb, sheetName="Variable Information")

    openxlsx::writeData(wb, 1, primaryStudies$summary)

    tmp1 <- grep("Source", colnames(primaryStudies$summary)); tmp1
    tmp2 <- grep("Orig.", colnames(primaryStudies$summary)); tmp2
    tmp3 <- grep("Delta", colnames(primaryStudies$summary)); tmp3
    openxlsx::writeData(wb, sheet2, primaryStudies$summary[c(tmp1, tmp2, tmp3)])

    tmp3 <- which(colnames(primaryStudies$summary) == "N")
    tmp4 <- grep("N\\(", colnames(primaryStudies$summary)); tmp4
    openxlsx::writeData(wb, sheet3, primaryStudies$summary[c(tmp1, tmp2, tmp3, tmp4)])


    tmp3 <- grep("r\\(", colnames(primaryStudies$summary)); tmp3
    tmp4 <- grep("N\\(", colnames(primaryStudies$summary)); tmp4
    openxlsx::writeData(wb, sheet4, primaryStudies$summary[c(tmp1, tmp2, tmp3, tmp4)])

    tmp3 <- grep("Moderator", colnames(primaryStudies$summary)); tmp3
    tmp8 <- primaryStudies$summary[c(tmp1, tmp2, tmp3)]; tmp8
    openxlsx::writeData(wb, sheet5, tmp8)
    if (!(is.null(moderatorLabelsBackup))) {
      tmp5 <- grep("Moderator", colnames(tmp8)); tmp5
      tmp6 <- rep("", (min(tmp5)-1)); tmp6
      tmp7 <- c(tmp6, moderatorLabelsBackup); tmp7
      tmp8 <- rbind(tmp8, tmp7); tmp8
      openxlsx::writeData(wb, sheet5, tmp8)
    }
    if (!(is.null(moderatorValuesBackup))) {
      maxCategories <- max(unlist(lapply(moderatorValuesBackup, function(extract) length(extract)))); maxCategories
      tmp5 <- grep("Moderator", colnames(tmp8)); tmp5
      tmp6 <- rep("", (min(tmp5)-1)); tmp6   # empty leading columns
      tmp7 <- c()
      for (j in 1: maxCategories){
        for (i in 1:length(moderatorValuesBackup)){
          tmp7 <- cbind(tmp7, paste0(moderatorValuesBackup[[i]][j]))
        }
      }
      tmp7 <- matrix(tmp7, ncol=length(tmp5), byrow=TRUE); tmp7  # matrix with labels
      tmp6b <- tmp6; tmp6b
      #if ((dim(tmp7)[1]-1) > 1) {
      if ((dim(tmp7)[1]-1) > 0) {
        for (i in 1:(dim(tmp7)[1]-1)) tmp6b <- rbind(tmp6b, tmp6); tmp6b
        tmp7 <- cbind(tmp6b, tmp7); tmp7
      } else {
        tmp7 <- matrix(c(tmp6b,tmp7), nrow=1)
      }
      tmp7[which(tmp7 =="NA")] <- ""
      colnames(tmp7) <- colnames(tmp8)
      rownames(tmp7) <- NULL
      tmp8 <- rbind(tmp8, tmp7); tmp8
      openxlsx::writeData(wb, sheet5, tmp8)
    }
    tmp3 <- grep("Country", colnames(primaryStudies$summary)); tmp3
    openxlsx::writeData(wb, sheet6, primaryStudies$summary[c(tmp1, tmp2, tmp3)])
    tmp3 <- grep("Occupation", colnames(primaryStudies$summary)); tmp3
    openxlsx::writeData(wb, sheet7, primaryStudies$summary[c(tmp1, tmp2, tmp3)])
    tmp3 <- grep("age", colnames(primaryStudies$summary)); tmp3
    tmp3 <- c(tmp3, grep("male", colnames(primaryStudies$summary))); tmp3
    openxlsx::writeData(wb, sheet8, primaryStudies$summary[c(tmp1, tmp2, tmp3)])
    tmp3 <- grep("Variable ", colnames(primaryStudies$summary)); tmp3  # space behind Variable is important
    tmp3 <- c(tmp3, grep("alpha", colnames(primaryStudies$summary))); tmp3
    tmp3 <- c(tmp3, grep("recodeVariables", colnames(primaryStudies$summary))); tmp3
    tmp3 <- c(tmp3, grep("combineVariables ", colnames(primaryStudies$summary))); tmp3
    tmp3 <- c(tmp3, grep("combineVariablesNames", colnames(primaryStudies$summary))); tmp3
    openxlsx::writeData(wb, sheet9, primaryStudies$summary[c(tmp1, tmp2, tmp3)])

    primaryStudies$excelSheets <- wb
  } # end if (summary == TRUE)

  primaryStudies$plot.type="none"

  primaryStudies$activeDirectory <- activeDirectory

  class(primaryStudies) <-  "CoTiMAFit"

  rm(ctma)
  return(primaryStudies)
}
