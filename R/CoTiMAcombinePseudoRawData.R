
#listOfStudyFits <- list(studyFit1, studyFit2, studyFit3, studyFit4, studyFit5)
#moderatorValues <- moderatorValues

#' Title
#'
#' @param listOfStudyFits ?
#' @param moderatorValues ?
#'
#' @return
#' @export
#'
#' @examples
#'
combinePseudoRawData <- function(listOfStudyFits=NULL, moderatorValues=NULL) {

  allSampleSizes <- maxLatents <- allTpoints <- c()
  nlatents <- listOfStudyFits[[1]]$ctmodelobj$n.latent; nlatents
  for (i in 1:length(listOfStudyFits)) {
    allSampleSizes <- c(allSampleSizes, dim(listOfStudyFits[[i]]$mxobj$data$observed)[1])
    allTpoints <- c(allTpoints, listOfStudyFits[[i]]$ctmodelobj$Tpoints)
  }
  currentVarnames <- c()
  for (j in 1:(max(allTpoints))) {
    for (h in 1:nlatents) {
      currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
    }
  }
  currentVarnames
  targetColNames <- c(c(currentVarnames, paste0("dT", seq(1:(max(allTpoints)-1))))); targetColNames

  alldata <- matrix(NA, 0, (nlatents*max(allTpoints)+max(allTpoints)-1)); dim(alldata)
  colnames(alldata) <- targetColNames; alldata

  groups <- moderatorGroups <- c()

  for (i in 1:length(listOfStudyFits)) {
    tmp <- listOfStudyFits[[i]]$mxobj$data$observed
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

    moderatorGroups <- c(moderatorGroups, rep(moderatorValues[i], dim(tmp2)[1]))
  }

  toReturn <- list(alldata=alldata, groups=groups, moderatorGroups=moderatorGroups)

  return(toReturn)
}
