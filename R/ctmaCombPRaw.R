#' ctmaCombPRaw
#'
#' @description Combine Pseudo Raw Data (extract them from 'CoTiMAFit object'$studyFitList)
#'
#' @param listOfStudyFits "Listobject of Studyfits"
#' @param moderatorValues "Moderators
#'
#' @return returns a pseudo raw data set that combines pseudo raw data and moderators of primary studies
#'
ctmaCombPRaw <- function(listOfStudyFits=NULL, moderatorValues=NULL) {

  allSampleSizes <- maxLatents <- allTpoints <- c()
  n.latent <- listOfStudyFits$n.latent; n.latent
  n.manifest <- listOfStudyFits$n.manifest; n.manifest
  if (is.null(n.manifest)) n.manifest <- 0

  allSampleSizes <- listOfStudyFits$statisticsList$allSampleSizes; allSampleSizes
  allTpoints <- listOfStudyFits$statisticsList$allTpoints; allTpoints; length(allTpoints)

  currentVarnames <- c()
  for (j in 1:(max(allTpoints))) {
    if (n.manifest == 0) {
      for (h in 1:n.latent) {
        currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
      }
    } else {
      for (h in 1:n.manifest) {
        currentVarnames <- c(currentVarnames, paste0("y",h,"_T", (j-1)))
      }
    }
  }
  targetColNames <- c(c(currentVarnames, paste0("dT", seq(1:(max(allTpoints)-1))))); targetColNames

  n.var <- max(n.latent, n.manifest); n.var
  alldata <- matrix(NA, 0, (n.var*max(allTpoints)+max(allTpoints)-1)); dim(alldata)
  colnames(alldata) <- targetColNames; alldata
  groups <- c()
  #moderatorValues
  if (!(is.null(moderatorValues))) {
    moderatorGroups <- matrix(NA, nrow=0, ncol=dim(moderatorValues)[[2]])
  }

  for (i in 1:length(listOfStudyFits$studyFitList)) {
    tmp <- listOfStudyFits$emprawList[[i]]
    # CHD Aug 2023
    # check of data file is loaded that has been created with old ctmaInit function (that did not save empraw data)
    tmp2 <- grep("V", colnames(tmp))
    #tmp3 <- mean(apply(tmp[tmp2], 2, mean, na.rm=TRUE))
    tmp3 <- all(tmp[tmp2] == 0)
    if (tmp3 == TRUE) {
      ErrorMsg <- "\nIt seems you used an outdated version of ctmaInit for initial fitting. I do not have access to pseudo raw data. Please update CoTiMA and run ctmaInit again."
      stop(ErrorMsg)
    }
    #
    tmp2 <- colnames(tmp); tmp2
    if (n.manifest > 0) colnames(tmp) <- gsub("V", "y", tmp2)
    tmp <- tmp[, colnames(tmp) %in% targetColNames]
    missingColNames <- colnames(alldata)[!(colnames(alldata) %in% colnames(tmp))]; missingColNames
    tmp2 <- matrix(NA, dim(tmp)[1], length(missingColNames));
    colnames(tmp2) <- missingColNames
    tmp2 <- cbind(tmp, tmp2)
    tmp2 <- as.data.frame(tmp2)
    # replace missing time values
    tmp3 <- grep("dT", colnames(tmp2)); tmp3
    tmp2[tmp3][is.na(tmp2[tmp3])] <- .00001
    alldata <- rbind(alldata, tmp2[, targetColNames])
    groups <- c(groups, rep(i, dim(tmp)[1]))
    #moderatorValues[i]
    if (!(is.null(moderatorValues))) {
      if (!(is.na(moderatorValues[i]))) {
      moderatorGroups <- rbind(moderatorGroups, (matrix(rep(unlist(moderatorValues[i, ]), dim(tmp2)[1]), nrow=dim(tmp2)[1], byrow=T)))
      } else {
        moderatorGroups <- rbind(moderatorGroups, (matrix(rep(NA, dim(tmp2)[1]), nrow=dim(tmp2)[1], byrow=T)))
      }
    }
  }
  casesToDelete <- c()
  if(exists("moderatorGroups")) {
    if (any(is.na(moderatorGroups))) {
      casesToDelete <- which(is.na(moderatorGroups)); casesToDelete
      alldata <- alldata[-casesToDelete]
      moderatorGroups <- moderatorGroups[-casesToDelete]
    }
  }

  if (!(is.null(moderatorValues))) {
    toReturn <- list(alldata=alldata, groups=groups, moderatorGroups=moderatorGroups, casesToDelete=casesToDelete)
  } else {
    toReturn <- list(alldata=alldata, groups=groups)
  }

  return(toReturn)
}
