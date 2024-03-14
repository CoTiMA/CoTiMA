#' ctmaPub
#'
#' @description Compute publication and citation scores for studies based on the (team of) authors' publication scores .
#'
#' @param getPubObj publication information compiled with \code{\link{ctmaGetPub}}
#' @param primaryStudyList vector with numbers of studies (e.g., c(1,3); requires source1 and source3 to be available)
#' @param yearsToExclude years to exclude from publications
#' @param targetYear year (default = last year) after which publications are ignored
#' @param recency years before targetYear that are considered for recency analysis
#' @param indFUN function (default = sum) how publications of each author within a collective (team) are summarized
#' @param colFUN function (default = mean) how publications all authors of collective (team) are summarized
#' @param addAsMod currently disabled. Add to existing moderator objects (or create them) in primaryStudyList, which is part of the returned object
#'
#' @importFrom crayon red
#'
#' @examples
#' \donttest{
#' pubResults_6 <- ctmaPub(getPubObj=pubList_8,
#'                         primaryStudyList=CoTiMAstudyList_6)
#' summary(pubResults_6)
#' }
#'
#' @export ctmaPub
#'
#' @return returns NEPP (= the \*number\* of studies published by the authors of the primary studies supplied UNTIL the year when the
#' primary study was published), NEPPRecency (like NEPP, but limited to the number of years before the publication as specified with the
#' recency argument), "Meaning of NEPP" and "Meaning of NEPPRecency" which explain what \*number\* exactly means (e.g., could be the mean
#' of the sum of each author's publication, or the sum of the maximum publications per year of the authors), and "primaryStudyList(full)",
#' which just returns the primaryStudyList supplied).
#'
ctmaPub <- function(getPubObj=NULL, primaryStudyList=NULL,
                    yearsToExclude=0, recency=5, targetYear=NULL,
                    indFUN="sum", colFUN="mean", addAsMod=FALSE) {

  if ( (!(indFUN %in% c("mean", "sum", "max", "min", "var"))) |
       (!(colFUN %in% c("mean", "sum", "max", "min", "var"))) ) {
    ErrorMsg <- "Unknown function to collect publication & citation information specified.  \nFunctions must be mean, sum, max, min, or var.  \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  # backup to use for labeling return objects later
  indFUNbackup <- indFUN
  colFUNbackup <- colFUN

  { # select specified function
    skip <- 0
    if (indFUN == "mean") {indFUN <- mean; skip <- 1}
    #if (skip == 0) if (indFUN == "sum") {indFUN <- function(x) sum(x); skip <- 1}
    if (skip == 0) if (indFUN == "sum") {indFUN <- sum; skip <- 1}
    if (skip == 0) if (indFUN == "max") {indFUN <- max; skip <- 1}
    if (skip == 0) if (indFUN == "min") {indFUN <- min; skip <- 1}
    if (skip == 0) if (indFUN == "var") {indFUN <- var; skip <- 1}

    skip <- 0
    if (colFUN == "mean") {colFUN <- mean; skip <- 1}
    #if (skip == 0) if (colFUN == "sum") {colFUN <- function(x) sum(x); skip <- 1}
    if (skip == 0) if (colFUN == "sum") {colFUN <- sum; skip <- 1}
    if (skip == 0) if (colFUN == "max") {colFUN <- max; skip <- 1}
    if (skip == 0) if (colFUN == "min") {colFUN <- min; skip <- 1}
    if (skip == 0) if (colFUN == "var") {colFUN <- var; skip <- 1}
  }

  source <- year <- NEPP <- list()
  messages1 <- messages2 <- messages3 <- c("all okay")
  counter1 <- counter2 <- counter3 <- 0
  NEPP <- NEPPrecency <- c() # NEed to Publish Papers (in recent years befor focal publication)

  for (i in 1:(length(primaryStudyList$source)-1)) {
    # get authors (sources) and years
    # source
    tmp1 <- primaryStudyList$source[[i]]; tmp1
    if (is.null(tmp1)) {
      ErrorMsg <- paste0("No source specified for Study no. ", primaryStudyList$studyNumbers[[i]], ".")
      stop(ErrorMsg)
    } else {
      source[[i]] <- tmp1; source[[i]]
    }
    # year
    year[[i]] <- as.numeric(source[[i]][length(source[[i]])]); year[[i]]
    if (is.na(year[[i]])) {
      ErrorMsg <- paste0("Year was not correctly specified for Study ", primaryStudyList$studyNumbers[[i]], ".")
      stop(ErrorMsg)
    }
    # target year
    if (is.null(targetYear)) targetYear2 <- as.numeric(year[[i]]) else targetYear2 <- targetYear; targetYear2
    n.authors <- n.authors2 <- which(source[[i]] == year[[i]]) - 1; n.authors
    pubCounter <- pubCounterRecency <- c()

    for (j in 1:n.authors) {
      # identify author record in results of ctmaGetPub
      tmpName <- source[[i]][j]; tmpName
      tmpName <- gsub("& ", "", tmpName); tmpName
      tmp <- grep(tolower(tmpName), tolower(getPubObj$Authors)); tmp
      if (length(tmp) > 1) {
        counter1 <- counter1 + 1
        messages1[counter1] <- paste0("Found multiple records for author ", tmpName, ". Locations were ", toString(tmp), ". Selected 1st occurence.")
        tmp <- tmp[1]
      }
      if (length(tmp) == 0) {
        counter3 <- counter3 + 1
        messages3[counter3] <- paste0("Found no matching name for author ", tmpName, ". in getPubObj. Check names. Currently this author is excluded.")
        tmp1a <- NA
      } else {
        # get cumulated pub frequenceis
        tmp1a <- getPubObj$`Publication Frequencies`[[grep(tolower(tmpName), tolower(getPubObj$Authors))[1]]]; tmp1a # grep 1st occurence in case of multiple entries
      }
      if (!(is.na(tmp1a[1]))) { # if there are publications for this author found
        # correction of cum pub frequencies if years should be excluded
        if (length(yearsToExclude) > 0)  {
          if (any(as.numeric(names(tmp1a)) %in% yearsToExclude))  { # if there are any publications to be excluded
            tmp2 <- max(tmp1a[as.numeric(names(tmp1a)) %in% yearsToExclude]); tmp2
            tmp1a[as.numeric(names(tmp1a)) %in% yearsToExclude] <- 0; tmp1a
            tmp1a[!(as.numeric(names(tmp1a)) %in% yearsToExclude)] <- tmp1a[!(as.numeric(names(tmp1a)) %in% yearsToExclude)] - tmp2
          }
        }
        tmp1a
        # cum pub freq across al years
        if (any(as.numeric(names(tmp1a)) < targetYear2)) {
          # single author's cum pub freq before current publication (pos in vector)
          tmp1b <- max(which(as.numeric(names(tmp1a)) < targetYear2)); tmp1b
          pubCounter[j] <- indFUN(tmp1a[1:tmp1b], na.rm=TRUE); pubCounter[j]
          # single author's cum pub freq in recent years before current publication (pos in vector)
          tmp1c <- which(as.numeric(names(tmp1a)) < (targetYear2-recency)); tmp1c
          if (!(is.na(tmp1c[1]))) tmp1c <- max(tmp1c) else tmp1c <- 0; tmp1c  #
          if (tmp1c > 0) pubCounterRecency[j] <- indFUN(tmp1a[1:tmp1b], na.rm=TRUE) - indFUN(tmp1a[1:tmp1c], na.rm=TRUE) else pubCounterRecency[j] <- NA; pubCounterRecency[j]
        } else {
          pubCounter[j] <- NA
          pubCounterRecency[j] <- NA
        }
      } else {
        counter2 <- counter2 + 1
        n.authors2 <- n.authors2 - 1
        messages2[counter2] <- paste0("No publications found for ", tmpName, ", who was therefor excluded from computations.")
        pubCounter[j] <- NA
        pubCounterRecency[j] <- NA
      }
    }
    if (n.authors2 > 0) {
      NEPP[i] <- colFUN(pubCounter, na.rm=TRUE); NEPP[i]
      NEPPrecency[i] <- colFUN(pubCounterRecency, na.rm=TRUE); NEPPrecency[i]
    } else {
      NEPP[i] <- NEPPrecency[i] <- NA
    }
  }

  ## create labels to explain output
  {
    if (is.null(targetYear)) tmp1 <- "the years of the groups' publication"  else tmp1 <- paste0("before ", targetYear); tmp1
    tmp2 <- ""
    if (indFUNbackup == "mean") tmp2 <- "mean in his/her yearly"; tmp2
    if (indFUNbackup == "sum") tmp2 <- "sum of all of his/her"; tmp2
    if (indFUNbackup == "var") tmp2 <- "variance in his/her yearly"; tmp2
    if (indFUNbackup == "max") tmp2 <- "maximum of his/her yearly"; tmp2
    if (indFUNbackup == "min") tmp2 <- "minimum of his/her yearly"; tmp2
    tmp3 <- ""
    if (colFUNbackup == "mean") tmp3 <- "mean"; tmp3
    if (colFUNbackup == "sum") tmp3 <- "sum"; tmp3
    if (colFUNbackup == "var") tmp3 <- "variance"; tmp3
    if (colFUNbackup == "max") tmp3 <- "maximum"; tmp3
    if (colFUNbackup == "min") tmp3 <- "minimum"; tmp3
    NEPPLabel <- paste0("NEPP represents the groups' ", tmp3, " of each study authors\' ", tmp2, " publications before ", tmp1); NEPPLabel
    if (is.null(targetYear)) tmp1 <- "the year of their study's publication"  else tmp1 <- paste0("before ", targetYear); tmp1
    NEPPRecencyLabel <- paste0("NEPPRecency is like NEPP, but is limited to the ", recency, " years before ", tmp1); NEPPRecencyLabel
  }

  ## add to existing moderator object (or create them) in primary study list
  if (addAsMod == TRUE) {
    counter <- 0
    NAcounter <- c()
    for (i in 1:length(NEPP)) {
      counter <- counter + 1
      if (!(is.na(NEPP[i]))) {
        NEPP_log <- log(NEPP[i]+1); NEPP_log
        NEPPrecency_log <- log(NEPPrecency[i]+1); NEPPrecency_log
        toAssign <- c(NEPP[i], NEPPrecency[i], NEPP_log, NEPPrecency_log); toAssign
        names(toAssign) <- c("NEPP", "NEPPrecency", "NEPP_log", "NEPPrecency_log")
        if (!(is.na(unlist(primaryStudyList$moderators[counter])[1]))) { # if other moderators are already present
          tmp1 <- names(primaryStudyList$moderators[[counter]])
          primaryStudyList$moderators[[counter]] <- c(unlist(primaryStudyList$moderators[[counter]]), toAssign)
          names(primaryStudyList$moderators[[counter]]) <- c(tmp1, "NEPP", "NEPPrecency", "NEPP_log", "NEPPrecency_log")
        } else {
          primaryStudyList$moderators[[counter]] <- toAssign
          names(primaryStudyList$moderators[[counter]]) <- c("NEPP", "NEPPrecency", "NEPP_log", "NEPPrecency_log")
        }
      } else {
        NAcounter <- c(NAcounter, counter)
        if (is.na(primaryStudyList$moderators[[counter]])) {
          primaryStudyList$moderators[[counter]] <- c(unlist(primaryStudyList$moderators[[counter]]), rep(NA, 3)) # 3 times because one NA is present
        } else {
          primaryStudyList$moderators[[counter]] <- c(unlist(primaryStudyList$moderators[[counter]]), rep(NA, 4)) # 4 times
        }
        names(primaryStudyList$moderators[[counter]]) <- c("NEPP", "NEPPrecency", "NEPP_log", "NEPPrecency_log")
      }
    }
  }

  useThis <- 0 # currently disabled because ctmaPrep only reads from the global environment
  if (useThis == 1) {
    # make new list without studies for which NEPP is NA
    objectNames <- c("delta_t", "sampleSize", "pairwiseN", "empcov", "moderator",
                     "studyNumber", "rawData", "empMeans", "empVars", "source",
                     "ageM",
                     #"ageSD",
                     "malePercent", "occupation", "country",  "alphas",
                     "targetVariables", "recodeVariables", "combineVariables", "combineVariablesNames",
                     "missingVariables"); objectNames
    namesInStudyList <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators",
                          #"startValues",
                          "studyNumbers", "rawData", "empMeans", "empVars", "source",
                          "ageM", "malePercent", "occupation", "country", "alphas",
                          "targetVariables", "recodeVariables", "combineVariables", "combineVariablesNames",
                          "missingVariables")
    # assign values to pre-defined object names
    counter <- 0
    for (i in namesInStudyList) {
      counter <- counter + 1; counter
      #i <- namesInStudyList[counter]; i
      # object names
      tmp1 <- rep(objectNames[counter], length(primaryStudyList[[i]])); tmp1
      tmp2 <- paste0(tmp1, primaryStudyList$studyNumbers); tmp2
      tmp2 <- tmp2[-length(tmp2)]; tmp2
      # values
      tmp3 <- primaryStudyList[[i]]; tmp3
      if ( length(tmp3) > length(tmp2) )  tmp3 <- tmp3[-length(tmp3)]; tmp3  # some lists do not have NA at the end
      for (j in 1:length(tmp2)) {
        if (!(is.na(NEPP[j]))) { assign(tmp2[j],tmp3[[j]], inherits=TRUE) } # inherits paste variables to the global environment (required for ctmaPrep)
        #if (!(is.na(NEPP[j]))) { assign(tmp2[j],tmp3[[j]], inherits=FALSE) } # inherits paste variables to the global environment (required for ctmaPrep)
      }
    }

    # make new (reduced) study list
    tmp1 <- which(!(is.na(NEPP))); tmp1
    studiesToSelect <- unlist(primaryStudyList$studyNumbers)[tmp1]; studiesToSelect

    if (exists("primaryStudyList$moderatorLabels")) {
      moderatorLabels <- c(primaryStudyList$moderatorLabels, "NEPP", "NEPPrecency", "NEPP_log", "NEPPrecency_log")
    } else {
      moderatorLabels <- c("NEPP", "NEPPrecency", "NEPP_log", "NEPPrecency_log")
    }
    if (exists("primaryStudyList$moderatorValues")) {
      moderatorValues <- c(primaryStudyList$moderatorValues, rep("countinuous", 4))
    } else {
      moderatorValues <- rep("countinuous", 4)
    }
    newStudyList <- ctmaPrep(selectedStudies=studiesToSelect,
                             moderatorLabels=moderatorLabels,
                             moderatorValues=moderatorValues)
    results <- list('NEPP'=NEPP,
                    'NEPPrecency'=NEPPrecency,
                    "Meaning of NEPP"=NEPPLabel,
                    "Meaning of NEPPrecency"=NEPPRecencyLabel,
                    "primaryStudyList(full)"=primaryStudyList,
                    "primaryStudyList(without NEPP=NA)"=newStudyList,
                    summary=list('NEPP'=NEPP,
                                 'NEPPrecency'=NEPPrecency,
                                 "Meaning of NEPP"=NEPPLabel,
                                 "Meaning of NEPPrecency"=NEPPRecencyLabel))
  } else {
    results <- list('NEPP'=NEPP,
                    'NEPPrecency'=NEPPrecency,
                    "Meaning of NEPP"=NEPPLabel,
                    "Meaning of NEPPrecency"=NEPPRecencyLabel,
                    "primaryStudyList(full)"=primaryStudyList,
                    summary=list('NEPP'=NEPP,
                                 'NEPPrecency'=NEPPrecency,
                                 "Meaning of NEPP"=NEPPLabel,
                                 "Meaning of NEPPrecency"=NEPPRecencyLabel))
  }

  class(results) <- "CoTiMAFit"

  invisible(results)
}

