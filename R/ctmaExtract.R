#' ctmaExtract
#'
#' @description Extracts objects required for CoTiMA from list of results returned by \code{\link{ctmaGenData}} to be used afterwards with \code{\link{ctmaPrep}}
#'
#' @param activeDirectory defines active directory where raw data files are saved. No default.
#' @param ctmaGenDataList list returned by \code{\link{ctmaGenData}}.
#' @param envir environment where objects should be extracted too. NULL by default. Typically one would use the global environment (envir = globalenv()). Has to be specified if ctmaExtract is set to TRUE.
#' @param useRawData if FALSE (default) empcov objects (empcov1, empcov2, etc.) are created in the global environment.
#' If TRUE, raw data files are created in the activeDirectory and rawData objects (rawData1, rawData2, etc.) are created in the global environment.
#' Note that corresponding empcov or rawData objects are deleted in the global environment if the respective other object is requested.

#' @importFrom utils write.table
#' @importFrom ctsem ctLongToWide ctIntervalise
#' @importFrom stats cor

#' @export ctmaExtract

ctmaExtract <- function(
    activeDirectory = NULL,
    ctmaGenDataList = NULL,
    envir = NULL,
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

  if (is.null(envir)) {
    ErrorMsg <- "\nThe argument envir has to be set to an environment (e.g., envir = globalenv()). \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  latentNames <- ctmaGenDataList[[1]]$latentNames

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
        suppressWarnings(datawide <- ctsem::ctIntervalise(datawide, Tpoints=length(ctmaGenDataList[[i]]$tpointTargets), n.manifest = length(latentNames), manifestNames = latentNames )
        )
      )
    )
    if (useRawData == FALSE) {
      targetVars <- grep("_", colnames(datawide)); targetVars
      assign(paste0("empcov", i), stats::cor(datawide[ , targetVars]), envir = envir)
      obj <- paste0("rawData", i)
      if (exists(obj, envir = envir)) rm(list = obj, envir = envir)
    } else {
      utils::write.table(datawide, paste0(activeDirectory, "data", i, ".txt"))
      assign(paste0("rawData", i), list(fileName=paste0(activeDirectory, "data", i, ".txt"), studyNumbers=i,
                                        missingValues=c(-999999), standardize=TRUE, header=TRUE, dec=".", sep=" "), envir = envir)
      obj <- paste0("empcov", i)
      if (exists(obj, envir = envir)) rm(list = obj, envir = envir)

    }
    assign(paste0("delta_t", i), diff(ctmaGenDataList[[i]]$tpointTargets), envir = envir)
    assign(paste0("sampleSize", i), nrow(datawide), envir = envir)
    assign(paste0("moderator", i), ctmaGenDataList[[i]]$modValues, envir = envir)
  }


}
