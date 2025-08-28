#' ctmaExtract
#'
#' @description Extracts objects required for CoTiMA from list of results returned by \code{\link{ctmaGenData} to be used afterwards with \code{\link{ctmaPrep}
#'
#' @param activeDirectory defines active directory where raw data files are saved. No default.
#' @param ctmaGenDataList list returned by \code{\link{ctmaGenData}.
#' @param useRawData if FALSE (default) empcov objects (empcov1, empcov2, etc.) are created in the global environment.
#' If TRUE, raw data files are created in the activeDirectory and rawData objects (rawData1, rawData2, etc.) are created in the global environment.
#' Note that corresponding empcov or rawData objects are deleted in the global environment if the respective other object is requested.

#' @importFrom utils write.table
#' @importFrom ctsem ctLongToWide ctIntervalise

#' @export ctmaExtract

ctmaExtract <- function(
    activeDirectory = NULL,
    ctmaGenDataList = NULL,
    useRawData = FALSE
)
{
  if (is.null(activeDirectory)) {
    ErrorMsg <- "\nNo active directory has been specified! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (is.null(ctmaGenDataList)) {
    ErrorMsg <- "\nNo ctmaGenDataList object has been specified! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  for (i in 1:length(ctmaGenDataList)) {
    #i <-1
    data <- ctmaGenDataList[[i]]$data
    datawide <- invisible(
      suppressMessages(
        suppressWarnings(as.data.frame(ctsem::ctLongToWide(data, id="id", time="time", manifestNames = latentNames))
        )
      )
    )
    invisible(
      suppressMessages(
        suppressWarnings(datawide <- ctsem::ctIntervalise(datawide, Tpoints=length(ctmaGenDataList[[i]]$tpointTargets), n.manifest = 2, manifestNames = latentNames )
        )
      )
    )
    if (useRawData == FALSE) {
      targetVars <- grep("_", colnames(datawide)); targetVars
      assign(paste0("empcov", i), cor(datawide[ , targetVars]), envir = .GlobalEnv)
      obj <- paste0("rawData", i)
      if (exists(obj, envir = .GlobalEnv)) rm(list = obj, envir = .GlobalEnv)
    } else {
      utils::write.table(datawide, paste0(activeDirectory, "data", i, ".txt"))
      assign(paste0("rawData", i), list(fileName=paste0(activeDirectory, "data", i, ".txt"), studyNumbers=i,
                                        missingValues=c(-999999), standardize=TRUE, header=TRUE, dec=".", sep=" "), envir = .GlobalEnv)
      obj <- paste0("empcov", i)
      if (exists(obj, envir = .GlobalEnv)) rm(list = obj, envir = .GlobalEnv)

    }
    assign(paste0("delta_t", i), diff(ctmaGenDataList[[i]]$tpointTargets), envir = .GlobalEnv)
    assign(paste0("sampleSize", i), nrow(datawide), envir = .GlobalEnv)
    assign(paste0("moderator", i), modValues[[i]], envir = .GlobalEnv)
  }


}
