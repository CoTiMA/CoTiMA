#' ctmaGetPub
#'
#' @description Retrieves publication and citation information from google scholar based on the supplied author names and their google ID (user)
#'
#' @param authorList list of authors and googe scholar addresses
#' @param yearsToExclude the years to be excluded (default = current year)
#' @param flush if TRUE, the cache will be cleared and the data reloaded from Google.
#'
#' @importFrom crayon red
#' @importFrom stringi stri_split_fixed
#' @importFrom scholar get_profile get_publications
#'
#' @export ctmaGetPub
#'
#' @note Set flush=TRUE only if retrieving is necessary (e.g., first retrieval on a day)
#'
#' @examples
#' \donttest{
#' pubList_8 <- ctmaGetPub(authorList = list( c("J; de Jonge",
#'               "https://scholar.google.de/citations?hl=de&user=0q27IckAAAAJ"),
#'               c("Arnold B.; Bakker", "user=FTl3bwUAAAAJ"),
#'               c("Evangelia; Demerouti", "user=9mj5LvMAAAAJ"),
#'               c("Joachim; Stoeber", "user=T9xdVusAAAAJ"),
#'               c("Claude; Fernet", "user=KwzjP4sAAAAJ"),
#'               c("Frederic; Guay", "user=99vnhX4AAAAJ"),
#'               c("Caroline; Senecal", "user=64ArFWQAAAAJ"),
#'               c("Stéphanie; Austin", "user=PPyTI7EAAAAJ")),
#'               flush=FALSE)
#' summary(pubList_8)
#' }
#'
#' @return list with (cumulative) frequencies and (cumulative) citations in google scholar
#'
ctmaGetPub <- function(authorList=NULL,
                       flush=FALSE, # if TRUE, the cache will be cleared and data is reloaded from Google
                       yearsToExclude=NULL) # exclude current year (default)

{
  researchers <- years <- pubs <- pubFreqs <- cumPubFreqs <- citeFreqs <- cumCiteFreqs <- list()
  nameInconsitencies <- list()
  researcherID <- c()

  if (is.null(yearsToExclude)) yearsToExclude <- 0
  if (is.na(yearsToExclude)) yearsToExclude <- 0

  if (length(authorList) == 0) {
    ErrorMsg <- "\nAt least one author name is required (e.g., list(\"Christian; Dormann\", NA))  \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (all(is.na(lapply(authorList, function(x) x[2])))) {
    Msg <- "No google scholar IDs (i.e., https) were provided. Not much to do. Processing continues anyway.   \n"
    message(Msg)
  }

  Msg <- "                       BE CAREFUL                                   \nGoogle Scholar might lock you out for a day or so if you retrieve \na particular author's information too frequently in short time.\n"
  message(Msg)

  counter <- 0
  for (i in 1:length(authorList)) {

    Msg <- paste0("Retrieving publication information of author no. ", i, " in the list provided. \n")
    message(Msg)

    researchers[[i]] <- (unlist(stringi::stri_split_fixed(authorList[[i]][1], "; "))); researchers[[i]]   # user provided name
    researcherID[i] <- unlist(stringi::stri_split_fixed(authorList[[i]][2], "user="))[2]; researcherID[i] # gs ID
    if (!(is.na(researcherID[i]))) {
      # use package scholar
      tmp1 <- suppressWarnings(unlist(stringi::stri_split_fixed(scholar::get_profile(researcherID[i])$name, " "))); tmp1 # get names (last = surname) from gs
      tmp2 <- grep(tolower(tmp1[length(tmp1)]), tolower(researchers[[i]][2])); tmp2  # check if gs name = user provided name
      if (length(tmp2) < 1) {
        counter <- counter +1; counter
        nameInconsitencies[[counter]] <- paste0(
          "Researcher's surname provided does not match surname extracted from google scholar.", "\n
          You provided the name", researchers[[i]][2], "with google scholar ID =", researcherID[i],
              " but I found the following name:", tmp1[length(tmp1)], ".",  "\n"
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
      years[[i]] <- min(years[[i]]):max(years[[i]]); years[[i]] # full range of publication years
      # publication frequencies
      pubFreqs[[i]] <- table(pubFreqs[[i]]); pubFreqs[[i]]            # vector of pub freqs without years where no. pubs = 0 (corrected in next para)
      minPubYear <- as.numeric(names(pubFreqs[[i]])[1]); minPubYear
      tmp6 <- 1:length(years[[i]]); tmp6 # create vector of proper length
      tmp6[tmp6>0] <- 0; tmp6            # set all to 0
      names(tmp6) <- years[[i]]; tmp6    # assign names
      tmp6[names(tmp6) %in% names(pubFreqs[[i]])] <- pubFreqs[[i]]; tmp6 # replace 0 by pub freqs
      pubFreqs[[i]] <- tmp6
      # citations
      #pubs[[i]]$year
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

  # if last authors do not have gs profile
  authors <- unlist(lapply(authorList, function(x) x[[1]])); authors
  if (length(years) < length(authors) ) {
    years[[length(authors)]] <- pubFreqs[[length(authors)]]  <- cumPubFreqs[[length(authors)]] <- NA
    citeFreqs[[length(authors)]]  <- cumCiteFreqs[[length(authors)]] <- NA
  }

  #names(years) <- names(pubFreqs) <- names(cumPubFreqs)  <- names(citeFreqs) <- names(cumCiteFreqs)  <- researchers
  names(years) <- names(pubFreqs) <- names(cumPubFreqs)  <- names(citeFreqs) <- names(cumCiteFreqs)  <- authors

  results <- list('inconsisten Names'=nameInconsitencies,
                  'Range of Publications Activities'=years,
                  'Publication Frequencies'=pubFreqs,
                  'Cum. Publication Frequencies'=cumPubFreqs,
                  'Citation Frequencies'=citeFreqs,
                  'Cum. Citation Frequencies'=cumCiteFreqs,
                  'All Publication Info'=pubs,
                  'Authors'= authors,
                  'AuthorList'=authorList,
                  summary=list('Publication Frequencies'=pubFreqs,
                               'Cum. Publication Frequencies'=cumPubFreqs,
                               'Citation Frequencies'=citeFreqs,
                               'Cum. Citation Frequencies'=cumCiteFreqs))

  class(results) <- "CoTiMAFit"

  invisible(results)
}
