# Combine Pseudo Raw Data (extrect them from ctmaInitFit$studyFitList)
#' Title
#'
#' @param listOfStudyFits "Listobject of Studyfits"
#' @param moderatorValues "Moderators
#'
#' @return
#' @export
#'
ctmaCombPRaw <- function(listOfStudyFits=NULL, moderatorValues=NULL) {
  allSampleSizes <- maxLatents <- allTpoints <- c()
  n.latent <- listOfStudyFits$n.latent; n.latent

  allSampleSizes <- listOfStudyFits$statisticsList$allSampleSizes; allSampleSizes
  allTpoints <- listOfStudyFits$statisticsList$allTpoints; allTpoints; length(allTpoints)

  currentVarnames <- c()
  for (j in 1:(max(allTpoints))) {
    for (h in 1:n.latent) {
      currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
    }
  }
  targetColNames <- c(c(currentVarnames, paste0("dT", seq(1:(max(allTpoints)-1))))); targetColNames

  alldata <- matrix(NA, 0, (n.latent*max(allTpoints)+max(allTpoints)-1)); dim(alldata)
  colnames(alldata) <- targetColNames; alldata

  groups <- c()
  #if (!(is.null(moderatorValues))) moderatorGroups <- matrix(NA, nrow=0, ncol=length(moderatorValues[[1]]))
  if (!(is.null(moderatorValues))) {
    moderatorGroups <- matrix(NA, nrow=0, ncol=dim(moderatorValues)[[2]])
  }

  for (i in 1:length(listOfStudyFits$studyFitList)) {
    tmp <- listOfStudyFits$emprawList[[i]]
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
    dim(alldata)
    groups <- c(groups, rep(i, dim(tmp)[1]))
    if (!(is.null(moderatorValues))) {
      moderatorGroups <- rbind(moderatorGroups, (matrix(rep(unlist(moderatorValues[i, ]), dim(tmp2)[1]), nrow=dim(tmp2)[1], byrow=T)))
    }
  }
  if (!(is.null(moderatorValues))) {
    toReturn <- list(alldata=alldata, groups=groups, moderatorGroups=moderatorGroups)
  } else {
    toReturn <- list(alldata=alldata, groups=groups)
  }

  return(toReturn)
}
