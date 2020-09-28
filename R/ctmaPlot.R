#######################################################################################################################
################################################## Plotting ###########################################################
#######################################################################################################################
#' ctmaPlot
#'
#' @param ctmaFit ""
#' @param activeDirectory ""
#' @param saveFilePrefix ""
#' @param activateRPB ""
#' @param plotCrossEffects ""
#' @param plotAutoEffects ""
#' @param timeUnit ""
#' @param timeRange ""
#' @param yLimitsForEffects ""
#' @param mod.values ""
#' @param mod.num ""
#' @param aggregateLabel ""
#'
ctmaPlot <- function(
  # Primary Study Fits
  ctmaFit=list(),                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,
  #sourceDirectory= NULL,
  #resultsFilePrefix="ctmaPlot",           # the prefix for the result file that is created
  saveFilePrefix="ctmaPlot",

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,                      #set to TRUE to receive push messages with CoTiMA notifications on your phone

  # Figure Parameters
  plotCrossEffects=TRUE,                  # plot cross effects across a range of discrete intervals
  plotAutoEffects=TRUE,                   # plot auto effects across a range of discrete intervals
  timeUnit="timeUnit (not specified)",    # timelag unit to lable x-axis of plots
  timeRange=c(),                          # used for Plotting and Poweranalysis. Populated by 0 to 1.5*maxDelta (Steptwidth = 1) if not specified as c(min,max,stepWidth)
  yLimitsForEffects=c(),                  # used the y-axis of Drift-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
  mod.values=-2:2,
  mod.num=1,

  aggregateLabel="SUM"

)
{  # begin function definition (until end of file)

  # check if fit object is specified
  if (is.null(ctmaFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat("A fitted ctma object has to be supplied to plot something. \n")
    stop("Good luck for the next try!")
  }

  # some re-arrangements to cover all possibilities from a single object (list) to a list of list of objects (lists)
  {
    # extract general arguments for plot (unspecific for a particular fit) if a combined list is used
    plotParamsLong <- c("activeDirectory", "sourceDirectory", "resultsFilePrefix",
                        "saveFilePrefix", "activateRPB", "plotCrossEffects",
                        "plotAutoEffects", "timeUnit", "timeRange", "yLimitsForEffects")
    plotParams <- plotParamsLong[-c(1,2)]; plotParams


    # if fit object is not provided as list with name Fit make correction
    #names(ctmaFit)
    tmp <- which(names(ctmaFit) %in% c("activeDirectory", "plot.type", "model.type", "coresToUse",
                                       "n.studies", "n.latent", "studyList", "studyFitList",
                                       "emprawList", # in some fits
                                       "data",       # in some fits
                                       "statisticsList", "modelResults", "parameterNames", "summary")); tmp
    if (length(tmp) >= 10) {
      tmp2 <- ctmaFit
      ctmaFit <- list()
      ctmaFit$Fit <- tmp2
    }

    # identify position of possible fit objects
    fitPos <- which(names(ctmaFit) == "Fit"); fitPos
    # plot specifications (could be undone by more spceial ones=
    plotSpecs <- which(names(ctmaFit) %in% plotParamsLong); plotSpecs
    if (length(plotSpecs) > 0) for (i in rev(plotSpecs)) {
      assign(paste0(names(ctmaFit)[i]), ctmaFit[[i]])
      ctmaFit[[i]] <- NULL
    }

    # identify further plotting specs (and, if found, lift fit object one level up)
    for (i in fitPos) {
      #i <- fitPos[2]; i
      names(ctmaFit[[i]])
      plotSpecs <- which(names(ctmaFit[[i]]) %in% plotParams); plotSpecs
      if (length(plotSpecs) > 0) {
        # assign plotting specs
        for (j in rev(plotSpecs)) assign(paste0(names(ctmaFit[[i]])[j]), ctmaFit[[i]][[j]])
        tmp <- unlist(ctmaFit[[i]], recursive=FALSE)
        names(tmp)
        # lift fit object
        ctmaFit[[i]] <- tmp
      }
    }
    names(ctmaFit)

    n.fitted.obj <- length(ctmaFit); n.fitted.obj

    plot.type <- list()
    tmp <- activeDirectory
    if (is.null(tmp)) activeDirectory <- c()
    for (i in 1:n.fitted.obj) {
      plot.type[[i]] <- ctmaFit[[i]]$plot.type
      if (is.null(tmp)) {
        activeDirectory[i] <- ctmaFit[[i]]$activeDirectory
      } else {
        activeDirectory[i] <- tmp
      }
    }

    # check if fit object can be plotted
    if (length(plot.type)==0) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat("The fitted CoTiMA object provided cannot be plotted. \n")
      stop("Good luck for the next try!")
    }

    if (length(unique(plot.type)) > 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat("All fitted CoTiMA object to plot have to be of the same plot.type. \n")
      cat("The following ploty.type arguments were found: \n")
      cat(unique(plot.type), "\n")
      stop("Good luck for the next try!")
    }

    if (length(unique(activeDirectory)) > 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat("More than a single active directors was sepcified, which does not work. \n")
      cat("The following activeDirectory arguments were found: \n")
      cat(unique(activeDirectory), "\n")
      stop("Good luck for the next try!")
    }
    activeDirectory <- activeDirectory[[1]]; activeDirectory
  }


  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################
  {
    driftNames <- study.numbers <- list()
    n.latent <- n.studies <- c()
    DRIFTCoeff <- DRIFTSE <- sampleSize <- list()
    FixedEffect_Drift <- FixedEffect_DriftLow <- FixedEffect_DriftUp <- list()
    maxDelta <- minDelta <- meanDelta <- c()
    allDeltas <- list()
    requiredSampleSizes <- list()

    #n.fitted.obj
    for (i in 1:n.fitted.obj) {
      #i <- 1
      #if (plot.type[[i]] == "power") plot.type[[i]] <- c("power", "drift")
      if ("power" %in% plot.type[[i]]) plot.type[[i]] <- c("power", "drift")
      n.latent[i] <- unlist(ctmaFit[[i]]$n.latent); n.latent[i]
      driftNames[[i]] <- ctmaFit[[i]]$parameterNames$DRIFT; driftNames[[i]]
      n.studies[i] <- ctmaFit[[i]]$n.studies; n.studies[i]
      study.numbers[[i]] <- unlist(lapply(ctmaFit[[i]]$studyList, function(extract) extract$originalStudyNo)); study.numbers[[i]]
      study.numbers[[i]] <- study.numbers[[i]][-length(study.numbers[[i]])]; study.numbers[[i]]
      # hard correction
      if (length(study.numbers[[i]]) < n.studies[i]) study.numbers[[i]] <- 1:n.studies[i]
      #study.numbers[[i]]
      if (n.studies[i] == 1) {
        DRIFTCoeff[[i]] <- list(ctmaFit[[i]]$modelResults$DRIFT); DRIFTCoeff[[i]]
      } else {
        DRIFTCoeff[[i]] <- ctmaFit[[i]]$modelResults$DRIFT; DRIFTCoeff[[i]]
      }
      sampleSize[[i]] <- ctmaFit[[i]]$statisticsList$allSampleSizes; sampleSize[[i]]
      if ( ("funnel" %in% plot.type[[i]]) || ("forest" %in% plot.type[[i]]) ) {
        if (n.studies[i] == 1) {
          DRIFTSE[[i]] <- list(ctmaFit[[i]]$modelResults$DRIFTSE); DRIFTSE[[i]]
        } else {
          DRIFTSE[[i]] <- ctmaFit[[i]]$modelResults$DRIFTSE; DRIFTSE[[i]]
        }
        FixedEffect_Drift[[i]] <-  ctmaFit[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`[2,]; FixedEffect_Drift[[i]]
        FixedEffect_DriftLow[[i]] <-  ctmaFit[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftLowerLimit",]; FixedEffect_DriftLow[[i]]
        FixedEffect_DriftUp[[i]] <-  ctmaFit[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftUpperLimit",]; FixedEffect_DriftUp[[i]]
      }
      #if (("drift" %in% plot.type[[i]]) || ("power" %in% plot.type[[i]])) {
      if ("drift" %in% plot.type[[i]]) {
        allDeltas[[i]] <- ctmaFit[[i]]$statisticsList$allDeltas; allDeltas[[i]]
        maxDelta[i] <- max(allDeltas[[i]], na.rm=TRUE); maxDelta[i]
        minDelta[i] <- min(allDeltas[[i]], na.rm=TRUE); minDelta[i]
        meanDelta[i] <- mean(allDeltas[[i]], na.rm=TRUE); meanDelta[i]
      }
      if ("power" %in% plot.type[[i]]) {
        requiredSampleSizes[[i]] <- ctmaFit[[i]]$summary$estimates$`Required Sample Sizes`
        statisticalPower <- ctmaFit[[i]]$summary$estimates$`Requested Statistical Power`
      }
    }
  }
  nlatent <- n.latent[[1]]; n.latent


  #######################################################################################################################
  ################################################### funnel plots ######################################################
  #######################################################################################################################

    for (k in 1:n.fitted.obj) {
      #k <- 1
      if ("funnel" %in% unlist(plot.type[[k]])) {
        stretch <- 1.2
        figureName <- c()
        for (i in 1:(n.latent[k]^2)) {
          #i <- 1
          # Determine range of X- and Y-axes
          pairs <- cbind(DRIFTCoeff[[k]][,i], DRIFTSE[[k]][,i]); pairs
          minX <- min(DRIFTCoeff[[k]][,i]); minX
          maxX <- max(DRIFTCoeff[[k]][,i]); maxX
          maxAbsX <- max(abs(c(minX, maxX))); maxAbsX
          avgX <-FixedEffect_Drift[[k]][i]; avgX
          yMax2 <- stretch * max(DRIFTSE[[k]][,i]); yMax2
          # Determine pseudo confidence intervals (http://www.metafor-project.org/doku.php/plots:funnel_plot_variations)
          lowLeft <- FixedEffect_Drift[[k]][i] - 1.96 * yMax2; lowLeft
          lowRight <- FixedEffect_Drift[[k]][i] + 1.96 * yMax2; lowRight
          xMin2 <- 0 - max(abs(lowLeft), abs(lowRight), abs(DRIFTCoeff[[k]][,i])); xMin2
          xMax2 <- 0 + max(abs(lowLeft), abs(lowRight), abs(DRIFTCoeff[[k]][,i])); xMax2

          graphics::plot.new()
          plot(pairs, #xlab=colnames(DRIFTCoeff[[k]])[i],
               ylab="Standard Error", xlab="Continous Time Drift Coefficient",
               xlim = c(xMin2, xMax2), ylim = c(yMax2, 0))
          graphics::polygon(c(lowLeft, avgX, lowRight), c(stretch * yMax2, 0, stretch * yMax2), lty=3)
          graphics::abline(v=avgX, col="black", lwd=1.5, lty=1)

          figureContent <- colnames(DRIFTCoeff)[i]; figureContent
          if (is.null(figureContent)) figureContent <- colnames(DRIFTCoeff[[1]])[i]; figureContent
          figureName[i] <- paste0(activeDirectory, saveFilePrefix, "_",  "Funnel Plot for ", figureContent, ".png")
          grDevices::dev.copy(grDevices::png, quality = 100, figureName[i], width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (i in 1:(n.latent^2))
      } # END if ("funnel" %in% plot.type)
    } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  ################################################### FORREST plots #####################################################
  #######################################################################################################################
  for (k in 1:n.fitted.obj) {
    if ("forest" %in% plot.type[[k]]) {
    # extracting information
    {
      # see Lewis & Clarke 2001 BMJ
      autoNames <- diag(matrix(colnames(DRIFTCoeff[[k]]), n.latent)); autoNames
      crossNames <- colnames(DRIFTCoeff[[k]])[!(colnames(DRIFTCoeff[[k]]) %in% autoNames)]; crossNames
      # lower limits
      autoDRIFTCoeffLow <- DRIFTCoeff[[k]][,autoNames] - 1.96 * DRIFTSE[[k]][,autoNames]; autoDRIFTCoeffLow
      crossDRIFTCoeffLow <- DRIFTCoeff[[k]][,crossNames] - 1.96 * DRIFTSE[[k]][,crossNames]; crossDRIFTCoeffLow
      minAutoDRIFTCoeffLow <-  min(autoDRIFTCoeffLow); minAutoDRIFTCoeffLow
      minCrossDRIFTCoeffLow <-  min(crossDRIFTCoeffLow); minCrossDRIFTCoeffLow
      # upper limits
      autoDRIFTCoeffUp <- DRIFTCoeff[[k]][, autoNames] + 1.96 * DRIFTSE[[k]][, autoNames]; autoDRIFTCoeffUp
      crossDRIFTCoeffUp <- DRIFTCoeff[[k]][, crossNames] + 1.96 * DRIFTSE[[k]][, crossNames]; crossDRIFTCoeffUp
      maxAutoDRIFTCoeffUp <-  max(autoDRIFTCoeffUp); maxAutoDRIFTCoeffUp
      maxCrossDRIFTCoeffUp <-  max(crossDRIFTCoeffUp); maxCrossDRIFTCoeffUp
      # average effects
      # sample sizes (used for size of squares)
      precision <- unlist(sampleSize[k]); precision
      precision <- matrix((rep(precision, n.latent^2)), ncol=n.latent^2); precision
      minPrecision <- min(precision)-.1*min(precision); minPrecision
      maxPrecision <- max(precision)+.1*max(precision); maxPrecision
    }

    # figure size
    yMax <- 300 # 300 # arbitrarily set
    xMax <- 300 # 200 # arbitrarily set

    # adaptations based on figures size ('base' objects contained transformed coefficients and are used for plotting)
    {
      heigthPerStudy <- yMax/(n.studies+1); heigthPerStudy
      #maxSquareSize <- heigthPerStudy/2; maxSquareSize
      maxSquareSize <- heigthPerStudy; maxSquareSize
      squareSizeBase <- precision/maxPrecision * maxSquareSize; squareSizeBase # the size of the square to plot (based on 1/DRIFTSE^1)
      # some algebra to adapt raw coefficients to the scale (xMax, yMax) used for plotting
      betaAuto <-  -(xMax)/(minAutoDRIFTCoeffLow - maxAutoDRIFTCoeffUp); betaAuto
      betaCross <-  -(xMax)/(minCrossDRIFTCoeffLow - maxCrossDRIFTCoeffUp); betaCross
      constAuto <- -minAutoDRIFTCoeffLow * betaAuto; constAuto
      constCross <- -minCrossDRIFTCoeffLow * betaCross; constCross
      # drifts of primary studies
      autoDRIFTCoeffBase <- DRIFTCoeff[[k]][,autoNames] * betaAuto + constAuto; autoDRIFTCoeffBase
      crossDRIFTCoeffBase <- DRIFTCoeff[[k]][,crossNames] * betaCross + constCross; crossDRIFTCoeffBase
      #apply(DRIFTCoeff[[k]][,crossNames], 2, mean)
      #apply(crossDRIFTCoeffBase, 2, mean)
      # lower limits of primary studies
      autoDRIFTCoeffLowBase <- autoDRIFTCoeffLow * betaAuto + constAuto; autoDRIFTCoeffLowBase
      crossDRIFTCoeffLowBase <- crossDRIFTCoeffLow * betaCross + constCross; crossDRIFTCoeffLowBase
      # upper limits of primary studies
      autoDRIFTCoeffUpBase <- autoDRIFTCoeffUp * betaAuto + constAuto; autoDRIFTCoeffUpBase
      crossDRIFTCoeffUpBase <- crossDRIFTCoeffUp * betaCross + constCross; crossDRIFTCoeffUpBase
      # avg. drift
      autoFixedEffect_DriftBase <- FixedEffect_Drift[[k]][autoNames]  * betaAuto + constAuto; autoFixedEffect_DriftBase
      crossFixedEffect_DriftBase <- FixedEffect_Drift[[k]][crossNames]  * betaCross + constCross; crossFixedEffect_DriftBase
      # avg. drift lower limit
      autoFixedEffect_DriftLowBase <- FixedEffect_DriftLow[[k]][autoNames]  * betaAuto + constAuto; autoFixedEffect_DriftLowBase
      crossFixedEffect_DriftLowBase <- FixedEffect_DriftLow[[k]][crossNames]  * betaCross + constCross; crossFixedEffect_DriftLowBase
      # avg. drift upper limit
      autoFixedEffect_DriftUpBase <- FixedEffect_DriftUp[[k]][autoNames] * betaAuto + constAuto; autoFixedEffect_DriftUpBase
      crossFixedEffect_DriftUpBase <- FixedEffect_DriftUp[[k]][crossNames] * betaCross + constCross; crossFixedEffect_DriftUpBase
    }

    # plot
    graphics::plot.new()
    for (i in 1:(n.latent^2)) {
      #i <- 3
      # set frame by plotting invisible object
      plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(0, xMax), ylim = c(0, 300),
           xaxt='n', yaxt='n', ann=FALSE)
      graphics::par(new=F)
      for (j in 1:n.studies) {
        #j <- 1
        crossDRIFTCoeffBase
        # sample size and effect size conditional on effect size (identical x-scale for all cross and auto effects, respectively)
        if (colnames(DRIFTCoeff[[k]])[i] %in% autoNames) {
          currentXPos <- autoDRIFTCoeffBase[j, colnames(DRIFTCoeff[[k]])[i]]
        } else {
          currentXPos <- crossDRIFTCoeffBase[j, colnames(DRIFTCoeff[[k]])[i]]
        }
        currentYPos <- yMax - (j-1) * heigthPerStudy - heigthPerStudy/2; currentYPos
        xleft <- currentXPos - squareSizeBase[j, i]/2; xleft
        xright <- currentXPos + squareSizeBase[j, i]/2; xright
        ytop <- currentYPos + squareSizeBase[j, i]/2; ytop
        ybottom <- currentYPos - squareSizeBase[j, i]/2; ybottom
        graphics::rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col="grey")
        # error bars
        if (colnames(DRIFTCoeff[[k]])[i] %in% autoNames) {
          xleft <- autoDRIFTCoeffLowBase[j,colnames(DRIFTCoeff[[k]])[i]]
          xright <- autoDRIFTCoeffUpBase[j,colnames(DRIFTCoeff[[k]])[i]]
        } else {
          xleft <- crossDRIFTCoeffLowBase[j,colnames(DRIFTCoeff[[k]])[i]]
          xright <- crossDRIFTCoeffUpBase[j,colnames(DRIFTCoeff[[k]])[i]]
        }
        ytop <- currentYPos + 0; ytop
        ybottom <- currentYPos ; ybottom
        graphics::rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col="black")
      } # END for (j in 1:n.studies)

      # average effect
      if (colnames(DRIFTCoeff[[k]])[i] %in% autoNames) {
        xmiddle <- autoFixedEffect_DriftBase[colnames(DRIFTCoeff[[k]])[i]]
        xleft <- autoFixedEffect_DriftLowBase[colnames(DRIFTCoeff[[k]])[i]]
        xright <- autoFixedEffect_DriftUpBase[colnames(DRIFTCoeff[[k]])[i]]
      } else {
        xmiddle <- crossFixedEffect_DriftBase[colnames(DRIFTCoeff[[k]])[i]]
        xleft <- crossFixedEffect_DriftLowBase[colnames(DRIFTCoeff[[k]])[i]]
        xright <- crossFixedEffect_DriftUpBase[colnames(DRIFTCoeff[[k]])[i]]
      }
      xmiddle; xleft; xright
      ytop <- currentYPos - maxSquareSize/2 ; ytop
      ymiddle <- ytop - maxSquareSize/2; ymiddle
      ybottom <- ymiddle - maxSquareSize/2; ybottom
      #x <- c(xleft, xmiddle, xright, xmiddle, xleft); x # if based in 1/SE as precision/square size
      x <- c(xmiddle-maxSquareSize/2, xmiddle, xmiddle+maxSquareSize/2, xmiddle, xmiddle-maxSquareSize/2); x # if based on N
      y <- c(ymiddle, ytop, ymiddle, ybottom, ymiddle); y
      graphics::polygon(x=x, y=y, col="black")
      graphics::abline(v=xmiddle, lty=2)

      # y-axis (original study numbers)
      #atSeq <- seq((n.studies*heigthPerStudy- heigthPerStudy/2), 0, by = -heigthPerStudy); atSeq
      atSeq <- seq(((n.studies+1)*heigthPerStudy- heigthPerStudy/2), 0, by = -heigthPerStudy); atSeq
      #labelsSeq <- study.numbers; labelsSeq
      labelsSeq <- c(unlist(study.numbers), "SUM"); labelsSeq
      graphics::axis(side=2, at = atSeq, labels=labelsSeq, las=1)
      graphics::axis(side=2, at = atSeq, labels=labelsSeq, las=1)
      # x-axis (effect size)
      atSeq <- seq(0, xMax, by = 10); atSeq
      #minDRIFTCoeffLow <- min(DRIFTCoeff[[k]][ ,i]); minDRIFTCoeffLow
      #maxDRIFTCoeffUp <- max(DRIFTCoeff[[k]][ ,i]); maxDRIFTCoeffUp
      #xRange <- maxDRIFTCoeffUp - minDRIFTCoeffLow; xRange
      #currentMin <- min(crossDRIFTCoeffLow); currentMin
      #currentMax <- max(crossDRIFTCoeffHigh); currentMax
      #minCrossDRIFTCoeffLow
      #maxCrossDRIFTCoeffUp
      if (colnames(DRIFTCoeff[[k]])[i] %in% autoNames) {
        currentMin <- min(autoDRIFTCoeffLow); currentMin
        currentMax <- max(autoDRIFTCoeffUp); currentMax
      } else {
        currentMin <- min(crossDRIFTCoeffLow); currentMin
        currentMax <- max(crossDRIFTCoeffUp); currentMax
      }
      xRange <- currentMax - currentMin; xRange
      labelsSeq <- round(seq(currentMin, currentMax, by = (xRange/(length(atSeq)-1))), 2); labelsSeq
      graphics::axis(side=1, at = atSeq, labels=labelsSeq)

      # Add labels and title
      graphics::title(main = paste0("Forest Plot for the Effects of ", colnames(DRIFTCoeff[[k]])[i]), sub = NULL,
            ylab = "Study Number", xlab="Continous Time Drift Coefficient")
      figureFileName <- paste0(activeDirectory, saveFilePrefix," Forest Plot for ", colnames(DRIFTCoeff[[k]])[i], ".png"); figureFileName
      grDevices::dev.copy(grDevices::png, quality = 100, figureFileName, width = 8, height = 8, units = 'in', res = 300)
      grDevices::dev.off()
    } # END for (i in 1:(n.latent^2))
  } # END if ("forest" %in% plot.type)
  } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  ############################################## discrete time plots  ###################################################
  #######################################################################################################################

  if ("drift" %in% unlist(plot.type)) {

    # check time unit
    #if (timeUnit=="timeUnit") {
    #  if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    #  cat(crayon::red$bold("The default time interval (timeUnit) has been chosen.", "\n"))
    #  cat(red$bold("The default time interval (timeUnit) has been chosen.", "\n"))
    #  cat(blue("Press 'q' to quit and specify or 'c' to continue. Press ENTER afterwards ", "\n"))
    #  char <- readline(" ")
    #  while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
    #    cat((blue("Please press 'q' to quit and specify time interval or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
    #    char <- readline(" ")
    #  }
    #  if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    #}

    # Function to compute discrete parameters using drift parameters and time-scaling factors
    discreteDrift <-function(driftMatrix, timeScale, number) {
      discreteDriftValue <- expm::expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }


    #######################################################################################################################
    ############## Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  ######################
    #######################################################################################################################

    # can only plot (overlay) mutiple fitted objects if n.latent is identical
    if (length(unique(unlist(n.latent))) > 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
      cat(crayon::red$bold("A fitted CoTiMA object has to be supplied to plot something. \n"))
      stop("Good luck for the next try!")
    }


    #######################################################################################################################
    ############################ Specification of Parameters for Plotting, Optimal Lags  ##################################
    #######################################################################################################################
    {
      ## timeRange, stepWidth, & noOfSteps
      maxDelta <- max(unlist(maxDelta)); maxDelta
      minDelta <- min(unlist(minDelta)); minDelta
      if (length(timeRange) < 1) {
        stepWidth <- 1
        usedTimeRange <- seq(0, 1.5*maxDelta, stepWidth)
        # add empirical lags not yet included
        usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
        noOfSteps <- length(usedTimeRange); noOfSteps
      } else {
        stepWidth <- timeRange[3]
        usedTimeRange <- seq(timeRange[1], timeRange[2], stepWidth)
        # add empirical lags not yet included
        usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
        noOfSteps <- length(usedTimeRange); noOfSteps
      }

      ## yLimitsForEffects
      # discrete effects across time range
      discreteDriftCoeff <- list()
      for (g in 1:n.fitted.obj) {
        # check if moderator values should be plotted
        if (!(is.null(ctmaFit[[g]]$modelResults$MOD))) {
          n.studies[g] <- length(mod.values); n.studies[g]
          n.mod <- ctmaFit[[g]]$n.moderators; n.mod
          targetDriftNames <- rownames(ctmaFit[[g]]$modelResults$MOD); targetDriftNames
          targetDriftNames <- gsub(paste0(ctmaFit[[g]]$mod.names[mod.num], "_on_"), "", targetDriftNames); targetDriftNames
          targetDriftNames <- gsub("_", "", targetDriftNames); targetDriftNames
          targetDriftNames <- targetDriftNames[targetDriftNames %in% driftNames[[g]]]; targetDriftNames
          targetRow <- mod.num * length(targetDriftNames) - length(targetDriftNames) + 1; targetRow
          targetRow <- targetRow:(targetRow+length(targetDriftNames)-1); targetRow
          DRIFTCoeff[[g]] <- list()
          counter <- 1
          for (i in mod.values) {
            # adjust labels for plotting and used time Range according to the time label positiont
            ctmaFit[[g]]$studyList[[counter]]$originalStudyNo <- i # used for labeling in plot
            tmp1 <- stats::quantile(usedTimeRange, probs = seq(0, 1, 1/(n.studies[g]+1))); tmp1 # used for positioning of moderator value in plot
            usedTimeRange <- sort(c(tmp1, usedTimeRange)); usedTimeRange # correcting for added time points
            noOfSteps <- length(usedTimeRange); noOfSteps
            ctmaFit[[g]]$studyList[[counter]]$delta_t <- tmp1[counter+1]
            ### compute moderated drift matrices
            # main effects
            tmp1 <- ctmaFit[[g]]$summary$estimates[,3]; tmp1
            tmp1 <- tmp1[(names(tmp1) == "DRIFT")]; tmp1
            tmp1 <- matrix(tmp1, n.latent[[g]], byrow=TRUE); tmp1 # main effect
            # moderator effects
            tmp2 <- ctmaFit[[g]]$summary$mod.effects[,3]; tmp2
            tmp2 <- matrix(tmp2, n.latent[[g]], byrow=TRUE); tmp2 # moderator effect to be added to main effect

            DRIFTCoeff[[g]][[counter]] <- tmp1 + (i) * tmp2; DRIFTCoeff[[g]][[counter]]
            counter <- counter +1
          }
          names(DRIFTCoeff[[g]]) <- paste0("Moderator Value = ", mod.values, " SD from mean if standardized (default setting)")
          allDiags <- c()
          for (i in 1:length(DRIFTCoeff[[g]])) allDiags <- c(allDiags, diag(DRIFTCoeff[[g]][[i]]))

          if (!(all(allDiags < 0))) {
            if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
            cat(crayon::red$bold("Some of the moderated drift matrices have values > 0 in their diagonals. \n"))
            cat(crayon::red$bold("This is likely if the model used to create \"ctmaFit\" was not identified! \n"))
            cat(crayon::red$bold("You may want to try smaller moderator values (e.g., \"mod.values=c(-.5, 0., .5)\")! \n"))
            cat(crayon::red$bold("Some of the moderated drift matrices have values > 0 in their diagonals. \n"))
            cat(crayon::red$bold("This is likely if the model used to create \"ctmaFit\" was not identified! \n"))
            cat(crayon::red$bold("You may want to try smaller moderator values (e.g., \"mod.values=c(-.5, 0., .5)\")! \n"))
            stop("Good luck for the next try!")
          }
        }
        discreteDriftCoeff[[g]] <- array(dim=c(n.studies[[g]], noOfSteps, n.latent[[g]]^2))
        for (h in 1:n.studies[g]) {
          #h <- 1
          for (i in usedTimeRange[1]:noOfSteps){
            #i <- usedTimeRange[35]; i
            timeValue <- i * stepWidth; timeValue
            discreteDriftCoeff[[g]][h, i, 1:(n.latent[[g]]^2)] <- c(discreteDrift(matrix(unlist(DRIFTCoeff[[g]][[h]]), n.latent[[g]], n.latent[[g]]), timeValue))
          }
        }
      }

      # min, max, & yLimitsForEffects
      coeffSeq <- seq(1, n.latent[[1]]^2, 1)[!(seq(1, n.latent[[1]]^2, 1) %in% seq(1, n.latent[[1]]^2, (n.latent[[1]]+1)))]; coeffSeq
      maxYLimitsForEffects <- -99999999
      minYLimitsForEffects <- 99999999
      if (length(yLimitsForEffects) < 1) {
        for (g in 1:n.fitted.obj) {
          yLimitsForEffectsTmp1 <- round(min(discreteDriftCoeff[[g]][, , coeffSeq]) - .10, 1); yLimitsForEffectsTmp1
          yLimitsForEffectsTmp2 <- round(max(discreteDriftCoeff[[g]][, , coeffSeq]) + .10, 1); yLimitsForEffectsTmp2
          if (yLimitsForEffectsTmp1 < minYLimitsForEffects) {
            yLimitsForEffects[1] <- yLimitsForEffectsTmp1
            minYLimitsForEffects <- yLimitsForEffects[1]
          }
          if (yLimitsForEffectsTmp2 > maxYLimitsForEffects) {
            yLimitsForEffects[2] <- yLimitsForEffectsTmp2
            maxYLimitsForEffects <- yLimitsForEffects[2]
          }
        }
      } ## END if (length(yLimitsForEffects) < 1)
    } ### END Specification of Parameters for Plotting, Statistical Power, Optimal Lags ###

    if ( (plotCrossEffects == TRUE) || (plotAutoEffects == TRUE) ) {

      print(paste0("#################################################################################"))
      print(paste0("################################### Plotting ####################################"))
      print(paste0("#################################################################################"))

      ##################################### SELECT DRIFT MATRICES ########################################

      # Drift matrix used for plotting effects of primary studies (just renamed to adapt to older plotting procedures)
      DriftForPlot <- DRIFTCoeff; DriftForPlot

      ##################################### COMPUTE DOTS FOR PLOTTING ########################################
      plotPairs <- list()
      dotPlotPairs <- list()
      for (g in 1:n.fitted.obj) {
        plotPairs[[g]] <- array(dim=c(n.studies[[g]], noOfSteps, 2+n.latent[[g]]^2))
        dotPlotPairs[[g]] <- array(dim=c(n.studies[[g]], noOfSteps, 2+n.latent[[g]]^2))
        for (currentTimeScale in 1:(noOfSteps-1)){
          for (h in 1:n.studies[[g]]) {
            timeValue <- usedTimeRange[currentTimeScale+1]; timeValue
            plotPairs[[g]][h,currentTimeScale,1] <- currentTimeScale; plotPairs[[g]][h,currentTimeScale,1]
            plotPairs[[g]][h,currentTimeScale,2] <- timeValue; plotPairs[[g]][h,currentTimeScale,2]
            for (j in 1:(n.latent[[g]]^2)) {
              plotPairs[[g]][h,currentTimeScale,(2+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              if (n.studies[[g]] == 1) tmp <- round(meanDelta[[1]],0) else tmp <- ctmaFit[[g]]$studyList[[h]]$delta_t
              if (timeValue %in% tmp) {
                dotPlotPairs[[g]][h, currentTimeScale, 1] <- currentTimeScale
                dotPlotPairs[[g]][h, currentTimeScale, 2] <- timeValue
                dotPlotPairs[[g]][h, currentTimeScale, (2+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              }
            }
          } # END for (h in 1:n.studies[[g]])
        } # END for (currentTimeScale in 0:noOfSteps)
      } # END for (g in 1:n.fitted.obj)


      ##################################### PLOTTING PARAMETERS ##########################################
      yMin <- minYLimitsForEffects; yMin
      yMax <- maxYLimitsForEffects; yMax
      xMax <- max(usedTimeRange); xMax
      xMin <- usedTimeRange[1]; xMin
      targetRows <- max(usedTimeRange)/stepWidth; targetRows


      ############################################ PLOTTING ##############################################

      ## PLOT (auto effects)
      if (plotAutoEffects == TRUE) {
        graphics::plot.new()
        figureFileNameAuto <- list()
        counter <- 0
        nlatent <- n.latent[[1]]; n.latent
        for (j in seq(1, nlatent^2, (nlatent+1))) { # diagonal elements only
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            if (is.null(ctmaFit[[g]]$type)) plot..type <- "l" else plot..type <- ctmaFit[[g]]$type; plot..type
            if (is.null(ctmaFit[[g]]$col)) plot.col <- "grey" else plot.col <- ctmaFit[[g]]$col; plot.col
            if (is.null(ctmaFit[[g]]$lwd)) plot.lwd <- 1.5 else plot.lwd <- ctmaFit[[g]]$lwd; plot.lwd
            if (is.null(ctmaFit[[g]]$lty)) plot.lty <- 1 else plot.lty <- ctmaFit[[g]]$lty; plot.lty
            if (is.null(ctmaFit[[g]]$xMin)) plot.xMin <- xMin else plot.xMin <- ctmaFit[[g]]$xMin; plot.xMin
            if (is.null(ctmaFit[[g]]$xMax)) plot.xMax <- xMax else plot.xMax <- ctmaFit[[g]]$xMax; plot.xMax
            if (is.null(ctmaFit[[g]]$yMin)) plot.yMin <- yMin else plot.yMin <- ctmaFit[[g]]$yMin; yMin
            if (is.null(ctmaFit[[g]]$yMax)) plot.yMax <- yMax else plot.yMax <- ctmaFit[[g]]$yMax; yMax
            if (is.null(ctmaFit[[g]]$dot.type)) dot.plot.type <- "b" else dot.plot.type <- ctmaFit[[g]]$dot.type; dot.plot.type
            if (is.null(ctmaFit[[g]]$dot.col)) dot.plot.col <- "black" else dot.plot.col <- ctmaFit[[g]]$dot.col; dot.plot.col
            if (is.null(ctmaFit[[g]]$dot.lwd)) dot.plot.lwd <- .5 else dot.plot.lwd <- ctmaFit[[g]]$dot.lwd; dot.plot.lwd
            if (is.null(ctmaFit[[g]]$dot.lty)) dot.plot.lty <- 3 else dot.plot.lty <- ctmaFit[[g]]$dot.lty; dot.plot.lty
            if (is.null(ctmaFit[[g]]$dot.pch)) dot.plot.pch <- 16 else dot.plot.pch <- ctmaFit[[g]]$dot.pch; dot.plot.pch
            if (is.null(ctmaFit[[g]]$dot.cex)) dot.plot.cex <- 2 else dot.plot.cex <- ctmaFit[[g]]$dot.cex; dot.plot.cex
            for (h in 1:n.studies[[g]]) {
              plot(plotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax), ylim = c(plot.yMin, 1),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFit[[g]]$plotStudyNo)) || (ctmaFit[[g]]$plotStudyNo==TRUE) ) {
                # black circle
                plot(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     xlim = c(xMin, xMax), ylim = c(yMin, 1),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (n.studies[[g]] > 1) {
                  currentLabel <- ctmaFit[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  if (is.null(currentLabel)) currentLabel <- ctmaFit[[g]]$ctmaFit$studyList[[h]]$originalStudyNo; currentLabel
                  if (h < 10) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  if (h > 9) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                } else {
                  currentLabel=aggregateLabel
                  graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                }
                graphics::par(new=T)
              }
            } # END for (h in 1:n.studies[[g]])
          } # END for (g in 1:n.fitted.obj)

          # plot y-axis
          plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(xMin, xMax), ylim = c(yMin, 1), xaxt='n', ann=FALSE, las=1)
          # Correct differences in axis length
          atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
          labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows*stepWidth/12)); labelsSeq
          if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
          if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
          graphics::axis(side=1, at = atSeq*stepWidth, labels=labelsSeq, las=2)
          # add labels and title
          if (!(is.null(ctmaFit[[g]]$modelResults$MOD))) {
            graphics::title(main = paste0("Moderated Auto-regressive Effects of V", counter), sub = NULL,
                  xlab=paste0("Time Interval in ", timeUnit), ylab = "Auto-regressive Beta")
          } else {
            graphics::title(main = paste0("Auto-regressive Effects of V", counter), sub = NULL,
                  xlab=paste0("Time Interval in ", timeUnit), ylab = "Auto-regressive Beta")
          }

          # SAVE
          graphics::par(new=F)
          tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); tmp
          grDevices::dev.copy(grDevices::png, quality = 100, tmp, width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (j in seq(1, nlatent^2, (nlatent+1)))
      } ## END PLOT (auto effects)


      ## PLOT (cross effects)
      if (plotCrossEffects == TRUE & nlatent > 1) {
        graphics::plot.new()
        counter <- 0
        coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
        for (j in coeffSeq) {
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            if (is.null(ctmaFit[[g]]$type)) plot..type <- "l" else plot..type <- ctmaFit[[g]]$type; plot..type
            if (is.null(ctmaFit[[g]]$col)) plot.col <- "grey" else plot.col <- ctmaFit[[g]]$col; plot.col
            if (is.null(ctmaFit[[g]]$lwd)) plot.lwd <- 1.5 else plot.lwd <- ctmaFit[[g]]$lwd; plot.lwd
            if (is.null(ctmaFit[[g]]$lty)) plot.lty <- 1 else plot.lty <- ctmaFit[[g]]$lty; plot.lty
            if (is.null(ctmaFit[[g]]$xMin)) plot.xMin <- xMin else plot.xMin <- ctmaFit[[g]]$xMin; plot.xMin
            if (is.null(ctmaFit[[g]]$xMax)) plot.xMax <- xMax else plot.xMax <- ctmaFit[[g]]$xMax; plot.xMax
            if (is.null(ctmaFit[[g]]$yMin)) plot.yMin <- yMin else plot.yMin <- ctmaFit[[g]]$yMin; yMin
            if (is.null(ctmaFit[[g]]$yMax)) plot.yMax <- yMax else plot.yMax <- ctmaFit[[g]]$yMax; yMax
            if (is.null(ctmaFit[[g]]$dot.type)) dot.plot.type <- "b" else dot.plot.type <- ctmaFit[[g]]$dot.type
            if (is.null(ctmaFit[[g]]$dot.col)) dot.plot.col <- "black" else dot.plot.col <- ctmaFit[[g]]$dot.col
            if (is.null(ctmaFit[[g]]$dot.lwd)) dot.plot.lwd <- .5 else dot.plot.lwd <- ctmaFit[[g]]$dot.lwd
            if (is.null(ctmaFit[[g]]$dot.lty)) dot.plot.lty <- 3 else dot.plot.lty <- ctmaFit[[g]]$dot.lty
            if (is.null(ctmaFit[[g]]$dot.pch)) dot.plot.pch <- 16 else dot.plot.pch <- ctmaFit[[g]]$dot.pch
            if (is.null(ctmaFit[[g]]$dot.cex)) dot.plot.cex <- 2 else dot.plot.cex <- ctmaFit[[g]]$dot.cex
            for (h in 1:n.studies[[g]]) {
              plot(plotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax), ylim = c(plot.yMin, 1),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFit[[g]]$plotStudyNo)) || (ctmaFit[[g]]$plotStudyNo==TRUE) ) {
                plot(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     xlim = c(xMin, xMax), ylim = c(yMin, 1),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (n.studies[[g]] > 1) {
                  currentLabel <- ctmaFit[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  if (h < 10) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  if (h > 9) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                } else {
                  currentLabel=aggregateLabel
                  graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                }
                graphics::par(new=T)
              }
            }
          } # END for (g in 1:n.fitted.obj)

          # plot y-axis
          plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(xMin, xMax), ylim = c(yMin, 1), xaxt='n',ann=FALSE, las=1)

          # Correct differences in axis length
          atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
          labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows/12)*stepWidth); labelsSeq
          if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
          if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
          graphics::axis(side=1, at = atSeq*stepWidth, labels=labelsSeq, las=2)
          # Add labels and title
          if (!(is.null(ctmaFit[[g]]$modelResults$MOD))) {
            graphics::title(main = paste0("Moderated Cross-lagged Effects of ", driftNames[[g]][j]), sub = NULL,
                  xlab=paste0("Time Interval in ", timeUnit), ylab = "Cross-lagged Beta")
          } else {
            graphics::title(main = paste0("Cross-lagged Effects of ", driftNames[[g]][j]), sub = NULL,
                  xlab=paste0("Time Interval in ", timeUnit), ylab = "Cross-lagged Beta")
          }

          graphics::par(new=F)
          tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); tmp
          grDevices::dev.copy(grDevices::png, quality = 100, tmp, width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (j in coeffSeq)
      } ## END PLOT (cross effects)
    } ### END if (plotCrossEffects == TRUE | plotAutoEffects == TRUE)
  }  ## END if ("drift" %in% plot.type)


  #######################################################################################################################
  ########################################## required sample size plots  ################################################
  #######################################################################################################################

  if ("power" %in% unlist(plot.type)) {
    graphics::plot.new()
    g <- 1 # only a single power plot
    if (is.null(ctmaFit[[g]]$pow.type)) pow.plot.type <- "b" else pow.plot.type <- ctmaFit[[g]]$pow.type
    if (is.null(ctmaFit[[g]]$pow.col)) {
      pow.plot.col <- rep(c("black", "grey"), length(statisticalPower))
      } else {
        pow.plot.col <- ctmaFit[[g]]$pow.col
      }
    if (is.null(ctmaFit[[g]]$pow.lwd)) pow.plot.lwd <- .5 else pow.plot.lwd <- ctmaFit[[g]]$pow.lwd
    if (is.null(ctmaFit[[g]]$pow.lty)) pow.plot.lty <- 3 else pow.plot.lty <- ctmaFit[[g]]$pow.lty
    if (is.null(ctmaFit[[g]]$pow.yMin)) pow.plot.yMin <- 0 else pow.plot.lty <- ctmaFit[[g]]$pow.yMin
    if (is.null(ctmaFit[[g]]$pow.yMin)) pow.plot.yMax <- 2000 else pow.plot.lty <- ctmaFit[[g]]$pow.yMax

    tmp1 <- suppressWarnings(as.numeric(rownames(requiredSampleSizes[[g]])))
    tmp1 <- tmp1[!(is.na(tmp1))]; tmp1
    tmp1 <- round(tmp1, 0); tmp1
    tmp2 <- !(duplicated(tmp1)); tmp2
    currentRequiredSamleSizes <- requiredSampleSizes[[g]][tmp2,]; currentRequiredSamleSizes
    tmp4 <- nrow(currentRequiredSamleSizes); tmp4
    currentRequiredSamleSizes <- currentRequiredSamleSizes[-c((tmp4-2):tmp4),]; currentRequiredSamleSizes
    usedTimeRange <- min(tmp1):max(tmp1); usedTimeRange
    xMax <- round(max(usedTimeRange), 0); xMax
    xMin <- usedTimeRange[1]; xMin
    currentLWD <- c(3, 2, 1); currentLWD

    coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
    currentDriftNames <- driftNames[[1]][coeffSeq]; currentDriftNames
      for (j in 1:length(currentDriftNames)) {
        offset <- (j-1)*length(statisticalPower); offset
        graphics::par(new=F)
        for (h in 1:length(statisticalPower)) {
          forPlotting <- cbind(as.numeric(rownames(currentRequiredSamleSizes)),
                               currentRequiredSamleSizes[, h]); forPlotting
          c(xMin, xMax)
          plot(forPlotting,
               type="l", col=pow.plot.col[h], lwd=currentLWD[h],
               main=paste0("Required Sample Size For the Effect of ", currentDriftNames[j]),
               xlab=paste0("Time Interval in ", timeUnit),
               ylab="Required Sample Size",
               xlim = c(xMin, xMax),
               ylim = c(pow.plot.yMin, pow.plot.yMax),
               xaxt='n' #, yaxt='n', ann=FALSE
          )
          graphics::par(new=T)
        }

        # Correct differences in axis length
        atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
        atSeq <- xMin:xMax; atSeq
        mxNumOfLabels <- 25
        stepWidth <- round(mxNumOfLabels/(max(usedTimeRange)+1)); stepWidth
        labelsSeq <- seq(0, (max(usedTimeRange)+1), 1); labelsSeq
        if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]
        if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]
        graphics::axis(side=1, at = atSeq, labels=labelsSeq)
        graphics::legend('bottomright', legend=statisticalPower, lty=1, col=pow.plot.col, lwd=currentLWD, bty='n', cex=.75)

        tmp <- paste0(activeDirectory, saveFilePrefix," RequiredSampleSizesFor ", currentDriftNames[j], ".png"); tmp
        grDevices::dev.copy(grDevices::png, quality = 100, tmp, width = 8, height = 8, units = 'in', res = 300)
        grDevices::dev.off()
      }

  } ### END Plotting ###

} ### END function definition
