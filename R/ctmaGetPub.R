#' ctmaGetPub
#'
#' @description Retrieves publication and citation information from google scholar
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
ctmaGetPub <- function(authorList=NULL,
                       activateRPB=FALSE,
                       flush=FALSE, # if TRUE, the cache will be cleared and data is reloaded from Google
                       yearsToExclude=as.integer(format(Sys.Date(), "%Y")) # exclude current year (default)
)
{
  researchers <- years <- pubs <- pubFreqs <- cumPubFreqs <- citeFreqs <- cumCiteFreqs <- list()
  nameInconsitencies <- list()
  researcherID <- c()

  if ( is.null(yearsToExclude) | is.na(yearsToExclude) ) yearsToExclude <- 0

  if (length(authorList) == 0) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(crayon::red$bold("At least one author name is required (e.g., list(\"Christian; Dormann\", NA))  \n"))
    stop("Good luck for the next try!")
  }

  if (all(is.na(lapply(authorList, function(x) x[2])))) {
    cat(crayon::red$bold("No google scholar IDs (i.e., https) were provided. Not much to do. Processing continues anyway.   \n"))
  }

  cat(crayon::red$bold("                       BE CAREFUL                                   \n"))
  cat(crayon::red$bold("Google Scholar might lock you out for a day or so if you retrieve \n"))
  cat(crayon::red$bold("a particular author's information too frequently in short time.\n"))
  cat(crayon::red$bold("                       BE CAREFUL                                   \n"))

  counter <- 0
  for (i in 1:length(authorList)) {

    cat(paste0("Retrieving publication information of author no. ", i, " in the list provided."), "\n")

    researchers[[i]] <- (unlist(stringi::stri_split_fixed(authorList[[i]][1], "; "))); researchers[[i]]   # user provided name
    researcherID[i] <- unlist(stringi::stri_split_fixed(authorList[[i]][2], "user="))[2]; researcherID[i] # gs ID
    if (!(is.na(researcherID[i]))) {
      # use package scholar
      tmp1 <- unlist(stringi::stri_split_fixed(scholar::get_profile(researcherID[i])$name, " ")); tmp1 # get names (last = surname) from gs
      tmp2 <- grep(tolower(tmp1[length(tmp1)]), tolower(researchers[[i]][2])); tmp2  # check if gs name = user provided name
      if (length(tmp2) < 1) {
        counter <- counter +1; counter
        nameInconsitencies[[counter]] <- paste0(
          cat("Researcher's surname provided does not match surname extracted from google scholar.", "\n"),
          cat("You provided the name", researchers[[i]][2], "with google scholar ID =", researcherID[i],
              " but I found the following name:", tmp1[length(tmp1)], ".",  "\n")
        )
      }
      # load author record
      #pubs[[i]] <- get_publications(researcherID[i], sortby = "year", flush=TRUE)
      pubs[[i]] <- scholar::get_publications(researcherID[i], sortby = "citation", flush=flush)
      # publications
      pubFreqs[[i]] <- pubs[[i]]$year; pubFreqs[[i]]   # all years extracted from all publications (multiple years possible)
      pubFreqs[[i]] <- pubFreqs[[i]][!(pubFreqs[[i]] %in% yearsToExclude)]; pubFreqs[[i]] # drop years to exclude
      pubFreqs[[i]] <- pubFreqs[[i]][!(is.na(pubFreqs[[i]]))]; pubFreqs[[i]] # drop NA
      # range of years
      years[[i]] <- sort(unique(unlist(pubFreqs[[i]]))); years[[i]]  # unqiue years
      yearsTmp <- years[[i]]; yearsTmp
      years[[i]] <- min(years[[i]]):max(years[[i]]) # full range of publication years
      # publication frequencies
      pubFreqs[[i]] <- table(pubFreqs[[i]]); pubFreqs[[i]]            # vector of pub freqs without years where no. pubs = 0 (corrected in next para)
      minPubYear <- as.numeric(names(pubFreqs[[i]])[1]); minPubYear
      tmp6 <- 1:length(years[[i]]); tmp6 # create vector of proper length
      tmp6[tmp6>0] <- 0; tmp6            # set all to 0
      names(tmp6) <- years[[i]]; tmp6    # assign names
      tmp6[names(tmp6) %in% names(pubFreqs[[i]])] <- pubFreqs[[i]]; tmp6 # replace 0 by pub freqs
      pubFreqs[[i]] <- tmp6
      # citations
      pubs[[i]]$year
      citeFreqs[[i]] <- pubs[[i]]$cites; citeFreqs[[i]]       # citations sort from current to first
      names(citeFreqs[[i]]) <- pubs[[i]]$year; citeFreqs[[i]] # assigne names (years)
      citeFreqs[[i]] <- citeFreqs[[i]][(!(is.na(names(citeFreqs[[i]]))))]; citeFreqs[[i]] # drop entries without year
      citeFreqs[[i]] <- citeFreqs[[i]][!(as.numeric(names(citeFreqs[[i]])) %in% yearsToExclude)]; citeFreqs[[i]] # drop years to exclude
      # vector of cite freqs
      tmp1 <- c()
      counter <- 0
      for (j in unique(names(citeFreqs[[i]]))) {
        counter <- counter + 1
        tmp1[counter] <- sum(citeFreqs[[i]][names(citeFreqs[[i]]) == j])
      }
      names(tmp1) <- unique(names(citeFreqs[[i]]))
      tmp1 <- rev(tmp1); tmp1
      tmp6[tmp6>0] <- 0; tmp6
      tmp6[names(tmp6) %in% names(tmp1)] <- tmp1; tmp6 # replace 0 by pub freqs
      citeFreqs[[i]] <- tmp6; citeFreqs[[i]]
      # cumulative Pub Frequencies
      tmpCumFreqs <- c();
      tmpCumFreqs <- pubFreqs[[i]][1]; tmpCumFreqs
      if (length(unlist(pubFreqs[[i]])) > 1) {
        for (j in 2:length(unlist(pubFreqs[[i]]))) tmpCumFreqs[j] <- tmpCumFreqs[j-1] + pubFreqs[[i]][j]
      }
      cumPubFreqs[[i]] <- tmpCumFreqs; cumPubFreqs[[i]]
      names(cumPubFreqs[[i]]) <- years[[i]]
      # cumulative Cite Frequencies
      tmpCumFreqs <- c();
      tmpCumFreqs <- citeFreqs[[i]][1]; tmpCumFreqs
      if (length(unlist(citeFreqs[[i]])) > 1) {
        for (j in 2:length(unlist(citeFreqs[[i]]))) tmpCumFreqs[j] <- tmpCumFreqs[j-1] + citeFreqs[[i]][j]
      }
      cumCiteFreqs[[i]] <- tmpCumFreqs; cumCiteFreqs[[i]]
      names(cumCiteFreqs[[i]]) <- years[[i]]
    } else {
      years[[i]] <-pubFreqs[[i]] <- cumPubFreqs[[i]]  <- citeFreqs[[i]] <- cumCiteFreqs[[i]] <- NA
    }
  }
  names(years) <- names(pubFreqs) <- names(cumPubFreqs)  <- names(citeFreqs) <- names(cumCiteFreqs)  <- researchers
  pubAnalysis <- list('inconsisten Names'=nameInconsitencies,
                      'Range of Publications Activities'=years,
                      'Publication Frequencies'=pubFreqs,
                      'Cum. Publication Frequencies'=cumPubFreqs,
                      'Citation Frequencies'=citeFreqs,
                      'Cum. Citation Frequencies'=cumCiteFreqs,
                      'All Publication Info'=pubs)

  invisible(pubAnalysis)
}
