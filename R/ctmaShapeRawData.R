#' ctmaShapeRawDataFiles
#'
#' @description Raw data objects are re-shaped (dealing with missing time points, wrong time intervals etc)
#'
#' @param dataFrame  an R object containing data
#' @param inputDataFrameFormat  "wide" or "long"
#' @param inputTimeFormat  "time" (default) or "delta"
#' @param missingValues  Missing value indicator, e.g., -999 or NA (default)
#' @param n.manifest  Number of process variables (e.g, 2 in a bivariate model)
#' @param Tpoints Number of time points in the data frame
#' @param allInputVariablesNames  Only required if the dataFrame (or the file header, see below) does not have column names
#' @param orderInputVariablesNames  "names" vs "time" (e.g., names: X1, X1, X3, Y1, Y2, X3 vs X1, Y1, X2, Y2, ... )
#' @param targetInputVariablesNames  The process variables in the dataFrame that should be used (in "names" or in "times" order)
#' @param targetInputTDpredNames  The actual TD names, e.g, 3, or 6, or 9, ... names if Tpoints = 3. Internally, each of the 3, 6, etc represents one TDpred.
#' @param targetInputTIpredNames  TIpred names in the dataFrame
#' @param targetTimeVariablesNames  The time variables name in the dataFrame
#' @param outputDataFrameFormat "long" (default) or "wide"
#' @param outputVariablesNames  "Y" (default), but can also be, e.g., c("X", "Y")
#' @param outputTDpredNames  Will become "TD" if not specified
#' @param outputTIpredNames  Will become "TI" if not specified
#' @param outputTimeVariablesNames "time" (default)
#' @param outputTimeFormat "time" (default) or "delta"
#' @param scaleTime  A scalar that is used to multiply the time variable. Typical use is to rescale primary study time to the time scale use in other primary studies
#' @param mininterval param supplied to ctIntervalise. Set to lower than any possible observed measurement interval, but above 0 - this is used for filling NA values where necessary and has no impact on estimates when set in the correct range.
#' @param minTolDelta Set, e.g. to 1/24, to delete variables from time points that are too close (1hr; or even before) after another time point.
#' @param maxTolDelta Set, e.g., to 7, to delete variables from time points that are too far after another time point (e.g., 7 days, if all cases should have responed within a week)
#' @param negTolDelta FALSE (default) or TRUE. Delete entire cases that have at least one negative delta ('unreliable responding'; use minTolDelta to delete certain variables only)
#' @param min.val.n.Vars default = 1. Minimum no. of valid variables per Tpoint. Default = 1 (retaines cases with only 1 valid variable), 0 would retain cases will all variables missing (not very useful).
#' @param min.val.Tpoints default = 1. Minimum no. of valid Tpoints (i.e. Tpoints where min.val.n.Var.per.Tpoint is met), default = 1 retains cases with full set of valid variables at least at one single Tpoint; 2 or more retains cases which provide longitudinal information.
#' @param experimental FALSE (default) or TRUE. Deprecated. (Was: Shift data left if all process variables are missing at a time point (even if time stamp is available)
#'
#' @examples
#' \dontrun{
#' tmpData <- data.frame(matrix(c(1,  2,  1, 2,  1, 2,  11, 26, 1,
#'                                NA, NA, 3, NA, 3, NA, 12, 27, 1,
#'                                1,  2,  1, 2,  1, 2,  NA, 24, 0 ),
#'                           nrow=3, byrow=TRUE))
#' colnames(tmpData) <- c("first_T0", "second_T0", "first_T1", "second_T1",
#'                          "TD1_0", "TD1_1",
#'                         "time1", "time2", "sex")
#' shapedData <- ctmaShapeRawData(dataFrame=tmpData,
#'                                inputDataFrameFormat="wide",
#'                                inputTimeFormat="time",
#'                                n.manifest=2,
#'                                Tpoints=2,
#'                                orderInputVariablesNames="time",
#'                                targetInputVariablesNames=c("first_T0", "second_T0",
#'                                                            "first_T1", "second_T1"),
#'                                targetInputTDpredNames=c("TD1_0", "TD1_1"),
#'                                targetInputTIpredNames="sex",
#'                                targetTimeVariablesNames=c("time1", "time2"),
#'                                scaleTime=1/12,
#'                                maxTolDelta=1.2)
#' head(shapedData)
#' }
#'
#' @importFrom  ctsem ctWideToLong ctDeintervalise
#' @importFrom  utils head
#'
#' @export ctmaShapeRawData
#'
#' @return A reshaped raw data file
#'
ctmaShapeRawData <- function(
    dataFrame=NULL,
    inputDataFrameFormat=NULL,
    inputTimeFormat="time",

    missingValues=NA,
    n.manifest=NULL,
    Tpoints=NULL,

    allInputVariablesNames=NULL,
    orderInputVariablesNames=NULL,
    targetInputVariablesNames=NULL,
    targetInputTDpredNames=NULL,
    targetInputTIpredNames=NULL,
    targetTimeVariablesNames=NULL,

    outputDataFrameFormat="long",
    outputVariablesNames="Y",
    outputTDpredNames=NULL,
    outputTIpredNames=NULL,
    outputTimeVariablesNames="time",
    outputTimeFormat="time",

    scaleTime=1,
    mininterval=0.0001,
    minTolDelta=NULL,
    maxTolDelta=NULL,
    negTolDelta=FALSE,

    min.val.n.Vars=1,
    min.val.Tpoints=1,

    experimental=FALSE
) {
  # some checks
  {
    if (is.null(n.manifest)) {
      ErrorMsg <- "\nThe number of manifest variables has to be specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (is.null(Tpoints)) {
      ErrorMsg <- "\nThe (maximum) number of time points has to be specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ( length(outputVariablesNames) > n.manifest) {
      ErrorMsg <- "\nYou provided more outputVariablesNames than you specified n.manifest! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ( !(orderInputVariablesNames) %in% c("names", "time")) {
      ErrorMsg <- "\nThe argument orderInputVariablesNames has to be either \"names\" or \"time\"! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ( !(inputTimeFormat) %in% c("time", "delta")) {
      ErrorMsg <- "\nThe argument inputTimeFormat has to be either \"time\" or \"delta\"! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ( !(outputTimeFormat) %in% c("time", "delta")) {
      ErrorMsg <- "\nThe argument inputTimeFormat has to be either \"time\" or \"delta\"! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.null(targetInputTDpredNames))) {
      if (length(targetInputTDpredNames) != Tpoints) {
        ErrorMsg <- "\nThe number of TD predictors names provided (\"targetInputTDpredNames\") should be equal to the number of time points (\"Tpoints\")! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (mininterval < .00001) {
      ErrorMsg <- "\nThe argument \"mininterval\" has been set to a value < .00001, which is currently not allowed! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (any(is.na(missingValues))) {
      Msg <- "Note: I assume that the missing values indicator in the dataFrame or dataFile is \"NA\" \n"
      message(Msg)
    }

    if (scaleTime == 1) {
      Msg <- "Note: Time is not scaled. \n"
      message(Msg)
    }

    if (!(is.null(minTolDelta))) {
      Msg <- paste0("Note: The shortest tolerated delta is ", minTolDelta, ". A subsequent time point closer to the preceeding one (afte possible time scaling) than ", minTolDelta," will be deleted. \n" )
      message(Msg)
    } else {
      minTolDelta = mininterval*2 # just slightly above the missing indicator
    }


    if (!(is.null(maxTolDelta))) {
      Msg <- paste0("Note: The longest tolerated Delta is ", maxTolDelta, ". All (!) subsequent time points (after possible time scaling) after T0 with an intervall larger than ", maxTolDelta," will be deleted. \n" )
      message(Msg)
    }

    if (is.null(maxTolDelta)) {
      Msg <- paste0("Note: All long deltas are specified to be acceptable (NULL). The shortest tolerate Delta is ", minTolDelta, ". \n" )
      message(Msg)
    }

    if ( !(is.null(allInputVariablesNames)) & (!(is.null(colnames(dataFrame)))) ) {
      if ( length(allInputVariablesNames) != length(colnames(dataFrame)) ) {
        ErrorMsg <- "\nThe argument \"allInputVariablesNames\" does not equal the no. of columns of the dataFrame provided! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      Msg <- "\nThe argument \"allInputVariablesNames\" has been provided, but the dataFrame provided has colnames, too. Take care you label variables correctly! \nGood luck for the next try!"
        message(ErrorMsg)
      }

    if (!(is.null(minTolDelta)) & !(is.null(maxTolDelta))) {
      if (minTolDelta > maxTolDelta) {
        ErrorMsg <- "\nThe argument minTolDelta has been set to a larger value than maxTolDelta  ! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (!(is.null(targetInputTDpredNames))) {
      tmp1 <- length(targetInputTDpredNames); tmp1
      if (tmp1/Tpoints != round(tmp1/Tpoints)) {
        ErrorMsg <- "\nThe number of variables specified in targetInputTDpredNames has to be a multifold of Tpoints! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (minTolDelta < mininterval) {
        ErrorMsg <- "\nThe argument minTolDelta has been set to a smaller value than mininterval (= indicator for missing)! \nGood luck for the next try!"
        stop(ErrorMsg)
    }


    tmp1 <- length(outputVariablesNames); tmp1
    if (tmp1 < n.manifest) {tmp1 <- rep(outputVariablesNames, n.manifest) } else {tmp1 <- outputVariablesNames}
    if (all(tmp1 == tmp1[1])) tmp1 <- paste0(tmp1[1], seq(1,length(tmp1),1))
    tmp2 <- rep("_T", 4); tmp2
    tmp3 <- paste0(tmp1, paste0(tmp2, c(0,0,1,1))); tmp3
    tmp3 <- paste(tmp3, collapse=" "); tmp3
    Msg <- paste0("Note: Output variable names will be ", tmp3, ", etc. \n")#)
    message(Msg)
  }


  ####################################################### Shape #######################################################
  ### Step 1 (Read raw data. Store in R-Object. Replacing missing value indicators with NA)
  tmpData <- data.frame(dataFrame)
  if (!(is.na(missingValues))) {
    tmp1 <- which(tmpData == missingValues, arr.ind = TRUE); tmp1
    tmpData[tmp1] <- NA
  }

  ### Step 2 - (Transpose data into wide format if they are in long format)
  {
  }

  ### Step 2b - (re-)label variables
  if ( !(is.null(allInputVariablesNames)) ) {
    colnames(dataFrame) <- allInputVariablesNames
  }

  # Step 3 (Select the desired "target variables" (at least X and Y and time) and kick out the remaining stuff.)
  #c(targetInputVariablesNames,  targetInputTDpredNames, targetTimeVariablesNames, targetInputTIpredNames)
  tmp1 <- c(targetInputVariablesNames,  targetInputTDpredNames, targetTimeVariablesNames, targetInputTIpredNames); tmp1
  tmpData <- tmpData[, tmp1]


  # Step 5 (Rename & re-arrange variables: X_T0, Y_T0, X_T1, Y_T1, ... time1, time2, ...)
  tmp1 <- length(outputVariablesNames); tmp1
  if (tmp1 < n.manifest) {tmp1 <- rep(outputVariablesNames, n.manifest) } else {tmp1 <- outputVariablesNames}
  if (all(tmp1 == tmp1[1])) tmp1 <- paste0(tmp1[1], seq(1,length(tmp1),1))
  newOutputVariablesNames <- tmp1; newOutputVariablesNames
  tmp2 <- sort(rep(seq(1, Tpoints, 1)-1, n.manifest)); tmp2
  allOutputVariablesNames <- paste0(tmp1, "_T", tmp2); allOutputVariablesNames
  # TD preds
  if (!(is.null(targetInputTDpredNames))) {
    n.TDpred <- length(targetInputTDpredNames)/Tpoints; n.TDpred
    if (is.null(outputTDpredNames)) {
      outputTDpredNames <- c()
      generalTDpredNames <- c()
      for (i in 1:n.TDpred) {
        generalTDpredNames <- c(generalTDpredNames, paste0("TD", i))
        for (j in 0:(Tpoints-1)) {
        outputTDpredNames <- c(outputTDpredNames, paste0("TD", i, "_T", j)); outputTDpredNames
        }
      }
    }
  } else {
    n.TDpred <- 0
    generalTDpredNames <- c()
    outputTDpredNames <- c()
  }
  # TI preds
  if (!(is.null(targetInputTIpredNames))) {
    n.TIpred <- length(targetInputTIpredNames); n.TIpred
    if (is.null(outputTIpredNames)) {
      outputTIpredNames <- paste0("TI", seq(1, length(targetInputTIpredNames), 1)); outputTIpredNames
    }
  } else {
    n.TIpred <- 0
    outputTIpredNames <- c()
  }
  # time
  allOutputTimeVariablesNames <- paste0("time", seq(0, (Tpoints-1), 1)); allOutputTimeVariablesNames

  # original order
  if (orderInputVariablesNames == "names") {
    variableOrder <- c()
    for (i in 1:Tpoints) {
      variableOrder <- c(variableOrder, seq(i, n.manifest*Tpoints, Tpoints)); variableOrder
    }
    targetInputVariablesNames <- targetInputVariablesNames[variableOrder]
  }

  tmpData <- tmpData[, c(targetInputVariablesNames, targetInputTDpredNames, targetTimeVariablesNames, targetInputTIpredNames)]
  if (inputTimeFormat == "delta") {
    dT0 <- data.frame(matrix(0, ncol=1, nrow=dim(tmpData)[1]))
    colnames(dT0) <- "dT0"
    tmpData <- cbind(tmpData[, c(targetInputVariablesNames, targetInputTDpredNames)],
                     dT0,
                     tmpData[, c(targetTimeVariablesNames, targetInputTIpredNames)])
  }
  colnames(tmpData) <- c(allOutputVariablesNames, outputTDpredNames, allOutputTimeVariablesNames, outputTIpredNames)

  #### Step 5b (make time out of delta if necessary)
  if (inputTimeFormat == "delta") {
    if (length(targetTimeVariablesNames) >= Tpoints) {
      ErrorMsg <- "\nYou specified time to be provided as time lags (deltas). The number of \"targetTimeVariablesNames\" provided exceeds the time lags in the data set! \nGood luck for the next try!"
      stop(ErrorMsg)
    }
    for (i in 1:(Tpoints-1)) {
      tmp1 <- tmpData[, paste0("time", i)] == mininterval
      tmpData[, paste0("time", i)] <- tmpData[ , paste0("time", i-1)] + tmpData[ , paste0("time", i)]
      tmpData[tmp1, paste0("time", i)] <- 0
      }
    allOutputTimeVariablesNames <- colnames(tmpData)[grep("time", colnames(tmpData))]; allOutputTimeVariablesNames
    tmp1 <- which(tmpData[, allOutputTimeVariablesNames[-1]] == 0, arr.ind = TRUE)
    tmpData[, allOutputTimeVariablesNames[-1]][tmp1] <- NA
  }

  # check
  if ( length(grep("time", colnames(tmpData)[c(allOutputVariablesNames, outputTDpredNames, outputTIpredNames)])) > 0) {
    if (minTolDelta > maxTolDelta) {
      ErrorMsg <- "\nThe name part \"time\" is only allowed in \"time\" or \"delta\" variables - not in latents, TIpreds, TDpreds! \nGood luck for the next try!"
      stop(ErrorMsg)
    }
  }

  ## at this stage, the variables should by in the order Y1_T0, Y2_T0, ..., Y1_T1, Y2_T1, ... TD1, TD2,... time1, time2,  ... TI1, TI2, ...

  # Step 6: Delete variables from time points for which no time stamp is available (without time information, ctsem is impossible)
  counter <- -1
  for (i in allOutputTimeVariablesNames) {
    counter <- counter + 1
    tmp1 <- which(is.na(tmpData[, i])); tmp1
    tmp2 <- grep(paste0("T", counter), allOutputVariablesNames); tmp2
    tmpData[tmp1, allOutputVariablesNames[tmp2]] <- NA
  }

  # Step 6b -  Scale time intervals
  tmpData[ , allOutputTimeVariablesNames] <- tmpData[ , allOutputTimeVariablesNames] * scaleTime

  # Step 6c - Delete all cases where all time stamps are missing
  if (inputTimeFormat == "time") { # if it is "delta" there should be at lease one time point
    tmp1 <- apply(tmpData[, allOutputTimeVariablesNames], 1, sum, na.rm=TRUE)
    tmp2 <- which(tmp1 == 0)
    if (length(tmp2) > 0) tmpData <- tmpData[-tmp2, ]
  }

  # Step 6d - Delete all cases where all process variables are missing
  tmp1 <- apply(tmpData[, allOutputVariablesNames], 1, sum, na.rm=TRUE)
  tmp2 <- which(tmp1 == 0)
  if (length(tmp2) > 0) tmpData <- tmpData[-tmp2, ]

  # Intermediate Step: delete cases for which conditions min.val.n.Vars and  min.val.Tpoints are not met
  # min.val.n.Vars
  #head(tmpData, 40)
  tmp1 <- apply(tmpData[, allOutputVariablesNames], 1, function(x) sum(!(is.na(x))))
  tmp2 <- which(tmp1 < min.val.n.Vars)
  if(length(tmp2) > 0 ) tmpData <- tmpData[-tmp2, ]
  # min.val.Tpoints
  validTpoints <- matrix(1, nrow=nrow(tmpData), ncol=Tpoints)
  for (i in 0:(Tpoints-1)) {
    tmp1 <- grep(paste0("T", i), colnames(tmpData))
    tmp2 <- apply(tmpData[, tmp1], 1, function(x) sum(!(is.na(x))))
    tmp3 <- which(tmp2 == 0)
    validTpoints[tmp3, i+1] <- 0
  }
  #head(validTpoints, 40)
  tmp1 <- apply(validTpoints, 1, function(x) sum(x))
  tmp2 <- which(tmp1 < min.val.Tpoints)
  if(length(tmp2) > 0 ) tmpData <- tmpData[-tmp2, ]


  # Step 6e - Shift data left if 1st time point is missing (otherwise lags will be not computed correctly later)
  tmpData2 <- tmpData
  n.TDpredPerWave <- length(targetInputTDpredNames)/Tpoints; n.TDpredPerWave
  for (t in 1:(Tpoints-1)) {
    # which T0 time stamp is missing
    tmp1 <- which(is.na(tmpData2[, allOutputTimeVariablesNames[1]])); tmp1
    # which substantive T0 variables are all missing
    tmp2 <- which(is.na(tmpData2[, allOutputVariablesNames[1:n.manifest]]), arr.ind = TRUE)
    tmp2 <- which(table(tmp2[, 1]) == n.manifest)
    tmp2 <- as.numeric(names(tmp2)); tmp2
    #combine
    tmp1 <- c(tmp1, tmp2); tmp1
    # shift substantive variables (allOutputVariablesNames)
    tmpData2[tmp1, allOutputVariablesNames[1:((Tpoints-1)*n.manifest)]] <-  tmpData2[tmp1, allOutputVariablesNames[(n.manifest+1):((Tpoints)*n.manifest)]]
    tmpData2[tmp1, allOutputVariablesNames[(n.manifest*(Tpoints-t)+1):((Tpoints+1-t)*n.manifest)]] <- NA
    # shift TDpreds (outputTDpredNames)
    tmpData2[tmp1, outputTDpredNames[1:((Tpoints-1)*n.TDpredPerWave)]] <-  tmpData2[tmp1, outputTDpredNames[(n.TDpredPerWave+1):((Tpoints)*n.TDpredPerWave)]]
    tmpData2[tmp1, outputTDpredNames[(n.TDpredPerWave*(Tpoints-t+1)):((Tpoints-t+1)*n.TDpredPerWave)]] <- NA
    # shift time variables
    tmpData2[tmp1, allOutputTimeVariablesNames[1:(Tpoints-t)]] <-  tmpData2[tmp1, allOutputTimeVariablesNames[(2):(Tpoints-t+1)]]
    tmpData2[tmp1, allOutputTimeVariablesNames[Tpoints+1-t]] <- NA
  }
  tmpData <- tmpData2

  # Step 6 Shift data left if all process variables are missing at a time point (even if time stamp is available)
  #if (experimental == TRUE) {
  if (Tpoints > 2) {
  for (tt in 2:(Tpoints-1)) {
      for (t in tt:(Tpoints-1)) {
        # which substantive T1 variables are all missing
        tmp2 <- which(is.na(tmpData2[, allOutputVariablesNames[((tt-1)*(n.manifest)+1):((tt-1)*(n.manifest)+n.manifest)]]), arr.ind = TRUE)
        tmp2 <- which(table(tmp2[, 1]) == n.manifest)
        tmp2 <- as.numeric(names(tmp2))
        # shift substantive variables (allOutputVariablesNames)
        tmpData2[tmp2, allOutputVariablesNames[((tt-1)*(n.manifest)+1):((Tpoints-1)*n.manifest)]] <- tmpData2[tmp2, allOutputVariablesNames[(tt*(n.manifest)+1):((Tpoints)*n.manifest)]]
        tmpData2[tmp2, allOutputVariablesNames[(n.manifest*(Tpoints-1)+1):(n.manifest*(Tpoints-1)+n.manifest)]] <- NA
        # shift TDpreds (outputTDpredNames)
        tmpData2[tmp2, outputTDpredNames[((tt-1)*(n.TDpredPerWave)+1):((Tpoints-1)*n.TDpredPerWave)]] <-  tmpData2[tmp2, outputTDpredNames[(tt*n.TDpredPerWave+1):((Tpoints)*n.TDpredPerWave)]]
        tmpData2[tmp2, outputTDpredNames[(n.TDpredPerWave*(Tpoints-1)+1):(n.TDpredPerWave*(Tpoints-1)+n.TDpredPerWave)]] <- NA
        # shift time variables
        tmpData2[tmp2, allOutputTimeVariablesNames[tt:(Tpoints-1)]] <- tmpData2[tmp2, allOutputTimeVariablesNames[(tt+1):(Tpoints)]]
        tmpData2[tmp2, allOutputTimeVariablesNames[(Tpoints)]] <- NA
        # delete last time stamp if process variables are missing at last time point
        tmp2 <- which(is.na(tmpData2[, allOutputVariablesNames[((Tpoints-1)*(n.manifest)+1):((Tpoints)*(n.manifest))]]), arr.ind = TRUE)
        tmp2 <- which(table(tmp2[, 1]) == n.manifest)
        tmp2 <- as.numeric(names(tmp2))
        tmpData2[tmp2, allOutputTimeVariablesNames[(Tpoints)]] <- NA
      }

      # delete time stamps and TDpreds if process variables are missing
      for (t in tt:(Tpoints-0)) {
        tmp2 <- which(is.na(tmpData2[, allOutputVariablesNames[((t-1)*n.manifest+1):(t*n.manifest)]]), arr.ind = TRUE)
        tmp2 <- which(table(tmp2[, 1]) == n.manifest)
        tmp2 <- as.numeric(names(tmp2))
        tmpData2[tmp2, allOutputTimeVariablesNames[(t)]] <- NA
      }
    }
    tmpData <- tmpData2
  }
  #head(tmpData)

  ### Step 6f - Determine possible lags that
  # - are longer than maxTolDelta
  # - are shorter than minTolDelta
  # and delete time points. Further, determine possible lags that
  # - are negative
  # and delete this cases (if negTolDelta is not set to TRUE)
  #
  # all possible lags (last value in name indicates the time point (0, 1, ... involved))
  tmp1 <- grep("time", colnames(tmpData)); tmp1
  timeMat <- tmpData[, tmp1]
  # test 1 wave lags first, then 2 wave lags, ... The first hit is the critical time point
  lagWidth <- 0
  for (j in 1:(Tpoints-1)) {
    lagWidth <- lagWidth + 1; lagWidth
    for (i in 1:(Tpoints-lagWidth)) {
      currentLags <- timeMat[,(i+lagWidth)]- timeMat[,i]; currentLags
      targetTimePoint <- i+lagWidth-1; targetTimePoint # 0, 1,
      timeVariableToDelete <- allOutputTimeVariablesNames[targetTimePoint+1]; timeVariableToDelete
      timePointsToDelete <- paste0("T", targetTimePoint); timePointsToDelete # just fro grepping the correct variable names
      variablesToDelete <- allOutputVariablesNames[c(grep(timePointsToDelete, allOutputVariablesNames))];variablesToDelete
      TDpredsToDelete <- outputTDpredNames[grep(timePointsToDelete, outputTDpredNames)]; TDpredsToDelete
      # delete variables involved in too short intervals
      targetCases <- which(currentLags < minTolDelta); targetCases
      #tmpData[targetCases, c(variablesToDelete, TDpredsToDelete, timeVariableToDelete)]
      if ( length(targetCases) > 0) tmpData[targetCases, c(variablesToDelete, TDpredsToDelete, timeVariableToDelete)] <- NA
      # delete variables involved in too long intervals
      targetCases <- which(currentLags > maxTolDelta); targetCases
      #tmpData[targetCases, c(variablesToDelete, TDpredsToDelete, timeVariableToDelete)]
      if ( length(targetCases) > 0) tmpData[targetCases, c(variablesToDelete, TDpredsToDelete, timeVariableToDelete)] <- NA
      # delete cases if a single delta is negative
      targetCases <- which(currentLags < 0); targetCases
      if ( (negTolDelta == FALSE) & (length(targetCases) > 0) ) {
        tmpData <- tmpData[-targetCases, ]
        timeMat <- timeMat[-targetCases, ]
      }
    }
  }

  # Step 7: ctIntervalise: Make time intervals out of time points if not already done.
  tmpData <- ctsem::ctIntervalise(tmpData, Tpoints = Tpoints, n.manifest = n.manifest,
                                  manifestNames =  newOutputVariablesNames,
                                  mininterval = mininterval,
                                  n.TDpred = n.TDpred, n.TIpred = n.TIpred,
                                  TDpredNames = generalTDpredNames, TIpredNames = outputTIpredNames)

  if ( (!(outputDataFrameFormat == "wide")) & (!(outputTimeFormat == "delta")) ) {
    # Step 8 (make long format again)
    tmpDataLong <- ctsem::ctWideToLong(tmpData, Tpoints = Tpoints, n.manifest = n.manifest,
                                       n.TDpred = n.TDpred, n.TIpred = n.TIpred,
                                       manifestNames =  newOutputVariablesNames,
                                       TDpredNames = generalTDpredNames, TIpredNames = outputTIpredNames)
    tmpDataLong <- data.frame(tmpDataLong)

    # Step 11 (ctsem::ctDeintervalise:)
    if (outputTimeFormat == "time") {
      tmpData <- ctsem::ctDeintervalise(tmpDataLong)
    }
  }


  ### make wide if required
  skip <- 0
  if (skip == 1) {
    if (outputDataFrameFormat == "wide") {
      tmpData <- ctsem::ctLongToWide(datalong = tmpData, id = "id", time = "time",
                                     manifestNames =  newOutputVariablesNames,
                                     TDpredNames = generalTDpredNames, TIpredNames = outputTIpredNames)
      if (outputTimeFormat == "delta") {
        tmpData <- ctIntervalise(datawide=tmpData,
                                 Tpoints=Tpoints,
                                 n.manifest=n.manifest,
                                 n.TDpred = n.TDpred,
                                 n.TIpred = n.TIpred,
                                 manifestNames = newOutputVariablesNames,
                                 TDpredNames = generalTDpredNames,
                                 TIpredNames = outputTIpredNames)
      }
    }
  }
  return(tmpData)
}
