#' ctmaPub
#'
#' @description Extract publication and citation scores for (multiple) authors in a list of studies
#'
#' @param activateRPB if TRUE, messages (warning, stops) could be send to smart phone (default = FALSE)
#' @param authorList list of authors and googe scholar addresses
#' @param yearsToExclude the years to be excluded (default = current year)
#' @param flush if TRUE, the cache will be cleared and the data reloaded from Google.
#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red
#' @importFrom stringi stri_split_fixed
#' @importFrom scholar get_profile get_publications
#'
#' @export ctmaGetPub
#'
#' @examples
#'
#' results <- ctmaGetPub(list( c("Wilmar B.; Schaufeli",
#'      "https://scholar.google.de/citations?hl=en&user=w1tHcj4AAAAJ"),
#'                             c("Maureen; Dollard",
#'      "user=J6oH3rgAAAAJ") ) )
#'
ctmaPub <- function(getPubObj=NULL, listOfStudies=NULL,
                    yearsToExclude=0, recency=5, targetYear=NULL,
                    indFUN=sum, colFUN=mean) {
  source <- year <- NEPP <- list()
  messages1 <- messages2 <- c("all okay")
  counter1 <- counter2 <- 0
  NEPP <- NEPPrecency <- c() # NEed to Publish Papers (in recent years befor focal publication)

  for (i in 1:length(listOfStudies)) {
    #i <- 1
    # get authors and years
    if (exists(paste0("source", listOfStudies[[i]]))) {
      source[[i]] <- get(paste0("source", listOfStudies[[i]])); source[[i]]
      year[[i]] <- get(paste0("year", listOfStudies[[i]])); year[[i]]
      if (is.null(targetYear)) targetYear <- as.numeric(year[[i]])
    } else {
      stop(paste0("No source specified for Study no. ", listOfStudies[[i]], "."))
    }
    n.authors <- n.authors2 <- which(source[[i]] == year[[i]]) - 1; n.authors
    pubCounter <- pubCounterRecency <- c()
    for (j in 1:n.authors) {
      #j <- 1
      # identify author record in results of ctmaGetPub
      tmpName <- source[[i]][j]; tmpName
      tmp <- grep(tolower(tmpName), tolower(getPubObj$Authors)); tmp
      if (length(tmp) > 1) {
        counter1 <- counter1 + 1
        messages1[counter1] <- paste0("Found multiple records for author ", tmpName, ". Locations were ", toString(tmp), ". Selected 1st occurence.")
        tmp <- tmp[1]
      }
      # get cumulated pub frequenceis
      tmp1a <- getPubObj$`Cum. Publication Frequencies`[[grep(tolower(tmpName), tolower(getPubObj$Authors))[1]]]; tmp1a # grep 1st occurence in case of multiple entries
      #tmp1a <- getPubObj$`Publication Frequencies`[[grep(tolower(tmpName), tolower(getPubObj$Authors))[1]]]; tmp1a2 # grep 1st occurence in case of multiple entries
      # correction of cum pub frequencies if years should be excluded
      if (length(excludeYears) > 0)  {
        if (any(as.numeric(names(tmp1a)) %in% excludeYears))  { # if there are any publications to be excluded
          tmp2 <- max(tmp1a[as.numeric(names(tmp1a)) %in% excludeYears]); tmp2
          tmp1a[as.numeric(names(tmp1a)) %in% excludeYears] <- 0; tmp1a
          tmp1a[!(as.numeric(names(tmp1a)) %in% excludeYears)] <- tmp1a[!(as.numeric(names(tmp1a)) %in% excludeYears)] - tmp2
        }
      }
      tmp1a
      # cum pub freq across al years
      if (!(is.null(tmp1a))) {
        if (any(as.numeric(names(tmp1a)) < targetYear)) {
        #if (any(as.numeric(names(tmp1a)) < as.numeric(year[[i]]))) {
          # single author's cum pub freq before current publication
          ##tmp1b <- max(tmp1a[which(as.numeric(names(tmp1a)) < as.numeric(year[[i]]))]); tmp1b
          #tmp1b <- max(which(as.numeric(names(tmp1a)) < as.numeric(year[[i]]))); tmp1b
          tmp1b <- max(which(as.numeric(names(tmp1a)) < targetYear)); tmp1b
          # single author's cum pub freq in recent years before current publication
          ##tmp1c <- max(tmp1a[which(as.numeric(names(tmp1a)) < as.numeric(year[[i]])-recency)]); tmp1c
          #tmp1c <- (which(as.numeric(names(tmp1a)) < as.numeric(year[[i]])-recency)); tmp1c
          tmp1c <- which(as.numeric(names(tmp1a)) < (targetYear-recency)); tmp1c
          if (!(is.na(tmp1c[1]))) tmp1c <- max(tmp1c) else tmp1c <- 0  #
          #pubCounter <- pubCounter + tmp1a[tmp1b]; pubCounter
          #pubCounter <- pubCounter + indFUN(tmp1a[1:tmp1b]); pubCounter
          pubCounter[j] <- indFUN(tmp1a[1:tmp1b]); pubCounter[j]
          #if (tmp1c > 0) pubCounterRecency <- pubCounterRecency  + tmp1a[tmp1b] - tmp1a[tmp1c] else pubCounterRecency <- pubCounterRecency + 0
          #if (tmp1c > 0) pubCounterRecency <- pubCounterRecency  + indFUN(tmp1a[1:tmp1b]) - indFUN(tmp1a[1:tmp1c]) else pubCounterRecency <- pubCounterRecency + 0
          if (tmp1c > 0) pubCounterRecency[j] <- indFUN(tmp1a[1:tmp1b]) - indFUN(tmp1a[1:tmp1c]) else pubCounterRecency[j] <- NA
          pubCounterRecency[j]
        } else {
          #pubCounter <- 0
          pubCounter[j] <- NA
        }
      } else {
        counter2 <- counter2 + 1
        n.authors2 <- n.authors2 - 1
        messages2[counter2] <- paste0("No publications found for ", tmpName, ", who was therefor excluded from computations.")
      }
      indFUN(tmp1a[1:tmp1b]) # cum pub freq until year before current publication
      indFUN(tmp1a[1:tmp1b]) - indFUN(tmp1a[1:tmp1c])  # cum pub freq in recency years before current publication
      pubCounter
      pubCounterRecency
    }
    if (n.authors2 > 0) {
      #NEPP[i] <- pubCounter/n.authors2
      NEPP[i] <- colFUN(pubCounter); NEPP[i]
      #NEPPrecency[i] <- pubCounterRecency/n.authors2
      NEPPrecency[i] <- colFUN(pubCounterRecency); NEPPrecency[i]
    } else {
      NEPP[i] <- NEPPrecency[i] <- NA
    }
  }
}




  xxx <- list('inconsisten Names'=nameInconsitencies,
                      'Range of Publications Activities'=years,
                      'Publication Frequencies'=pubFreqs,
                      'Cum. Publication Frequencies'=cumPubFreqs,
                      'Citation Frequencies'=citeFreqs,
                      'Cum. Citation Frequencies'=cumCiteFreqs,
                      'All Publication Info'=pubs)

  invisible(pubAnalysis)
}
