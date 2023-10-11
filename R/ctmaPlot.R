#' ctmaPlot
#'
#' @description Forest plot, funnel plots, plots of discrete time cross-lagged and autoregressive effect, and plots of required sample sizes
#'
#' @param ctmaFitObject 'CoTiMA' Fit object
#' @param activeDirectory  defines another active directory than the one used in ctmaInitFit
#' @param saveFilePrefix Prefix used for saved plots
#' @param activateRPB  set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param plotCrossEffects logical
#' @param plotAutoEffects logical
#' @param timeUnit label for x-axis when plotting discrete time plots
#' @param timeRange vector describing the time range for x-axis as sequence from/to/stepSize (e.g., c(1, 144, 1))
#' @param yLimitsForEffects range for y-axis
#' @param mod.values moderator values that should be used for plots
#' @param mod.number moderator number that should be used for plots
#' @param aggregateLabel label to indicate aggregated discrete time effects
#' @param xLabels labes used for x-axis
#' @param undoTimeScaling if TRUE, the original time scale is used (timeScale argument possibly used in \code{\link{ctmaInit}} is undone )
#' @param ... arguments passed through to plot()
#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red
#' @importFrom graphics plot plot.new polygon abline par rect axis title text legend
#' @importFrom grDevices dev.copy png dev.off
#' @importFrom OpenMx expm
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#' # cannot run without proper activeDirectory specified. Adapt!
#' CoTiMAFullFit_3$activeDirectory <- "/Users/tmp/" # adapt!
#' plot(ctmaFitList(CoTiMAInitFit_3, CoTiMAFullFit_3),
#'      timeUnit="Months", timeRange=c(1, 144, 1),
#'      plotAutoEffects=FALSE)
#' }
#'
#' @examples
#' \dontrun{
#' # cannot run without proper activeDirectory specified. Adapt!
#' CoTiMABiG_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' plot(CoTiMABiG_D_BO)
#' }
#'
#' @export ctmaPlot
#'
#' @return depending on the CoTiMA fit object supplied, generates funnel plots, forest plots, discrete time plots of
#' autoregressive and cross-lagged effects, plots of required samples sizes across a range of discrete time intervals
#' to achieve desired levels of statistical power, and post hoc power of primary studies. Plots are saved to disk.
#'
ctmaPlot <- function(
    ctmaFitObject=NULL,
    activeDirectory=NULL,
    saveFilePrefix="ctmaPlot",
    activateRPB=FALSE,
    plotCrossEffects=TRUE,
    plotAutoEffects=TRUE,
    timeUnit="timeUnit (not specified)",
    timeRange=c(),
    yLimitsForEffects=c(),
    mod.number=1,
    mod.values=-2:2,
    aggregateLabel="",
    xLabels=NULL,
    undoTimeScaling=TRUE,
    ...
)
{  # begin function definition (until end of file)

  par.original <- par("new"); par.original
  on.exit(par(new=par.original))


  { # some checks

    # check if fit object is specified
    if (is.null(ctmaFitObject)){
      if (activateRPB==TRUE) {
        RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ),
                            paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
      }
      ErrorMsg <- "A fitted ctma object has to be supplied to plot something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    # check #1 if object can be plotted
    #if (class(ctmaFitObject) == "list") testObject <- ctmaFitObject[[1]] else testObject <- ctmaFitObject
    if (is(ctmaFitObject) == "list") testObject <- ctmaFitObject[[1]] else testObject <- ctmaFitObject
    #if (class(testObject) != "CoTiMAFit")  {
    if (!(is(testObject, "CoTiMAFit"))) {
      if (activateRPB==TRUE) {
        RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ),
                            paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
      }
      ErrorMsg <- "This is nothing CoTiMA-related that I can plot. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    # some re-arrangements to cover all possibilities from a single object (list) to a list of list of objects (lists)
    if (!(is.null(names(ctmaFitObject)))) {
      tmp2 <- ctmaFitObject
      ctmaFitObject <- list()
      ctmaFitObject[[1]] <- tmp2
    }

    # check #2 if object can be plotted
    if ("none" %in% ctmaFitObject[[1]]$plot.type)  {
      if (activateRPB==TRUE) {
        RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ),
                            paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
      }
      ErrorMsg <- "This is nothing CoTiMA-related that I can plot. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (nchar(aggregateLabel) > 3) {
      Msg <- "The aggregate label has 4 or more characters. Plots will probably do not look nice.\n"
      message(Msg)
    }


    n.fitted.obj <- length(ctmaFitObject); n.fitted.obj # has to be done twice

    plot.type <- list() # has to be a list because a single fit could be used for different plots (e.g. "power")
    tmp <- activeDirectory; tmp
    if (n.fitted.obj == 1) {
      plot.type[[1]] <- ctmaFitObject[[1]]$plot.type; plot.type[[1]]
      if (is.null(tmp)) {
        activeDirectory <- ctmaFitObject[[1]]$activeDirectory
        if (is.null(activeDirectory)) activeDirectory <- ctmaFitObject$activeDirectory
        # CHD 27.2.23 added
        if (is.null(activeDirectory)) activeDirectory <- ctmaFitObject$argumentList$activeDirectory
      } else {
        activeDirectory <- tmp
      }
    }
    #activeDirectory

    if (n.fitted.obj != 1) {
      for (i in 1:n.fitted.obj) {
        if (is.null(ctmaFitObject[[1]]$plot.type)) {
          plot.type[[i]] <- "drift"
        } else {
          plot.type[[i]] <- ctmaFitObject[[i]]$plot.type
        }
      }
      if (!(is.null(tmp))) { # CHD 27.2.2023
        activeDirectory <- tmp
      }
    }

    if (is.null(activeDirectory)) activeDirectory <- ctmaFitObject[[1]]$activeDirectory
    if (is.null(activeDirectory)) activeDirectory <- ctmaFitObject[[1]]$argumentList$activeDirectory
    if (is.null(activeDirectory)) activeDirectory <- ctmaFitObject$activeDirectory
    if (is.null(activeDirectory)) activeDirectory <- ctmaFitObject$argumentList$activeDirectory


    # check if fit object can be plotted
    if (length(unlist(plot.type))==0) {
      if (activateRPB==TRUE) {
        RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ),
                            paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
      }
      ErrorMsg <- "The fitted CoTiMA object provided cannot be plotted. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (length(unique(unlist(plot.type))) > 1) {
      if (any(!(unique(unlist(plot.type)) %in% c("funnel", "forest")))) {
        if (activateRPB==TRUE) {
          RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ),
                              paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
        }
        ErrorMsg <- paste0("All fitted CoTiMA object to plot have to be of the same plot.type. \nThe following ploty.type arguments were found: \n", unique(unlist(plot.type)), "\nGood luck for the next try!")
        stop(ErrorMsg)
      }
    }

    if (length(unique(unlist(activeDirectory))) > 1) {
      if (activateRPB==TRUE) {
        RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ),
                            paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))
      }
      ErrorMsg <- paste0("More than a single active directors was sepcified, which does not work. \nThe following activeDirectory arguments were found: \n", unique(activeDirectory), "\nGood luck for the next try!")
      stop(ErrorMsg)
    }

    # detect study no > 4 nchar
    tmp1 <- c()
    for (i in 1:n.fitted.obj)  tmp1 <- c(tmp1, unlist(lapply(ctmaFitObject[[i]]$studyList, function(extract) extract$originalStudyNo)))
    if (!(is.null(tmp1))) {
      if (any(nchar(tmp1) > 4)) {
        Msg <- "Some orginal study numbers have 4 or more digits. Plots will probably do not look nice.\n"
        message(Msg)
      }
    }

  } # end some checks

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
    mod.values.backup <- mod.values; mod.values.backup
    mod.values <- list()
    n.primary.studies <- c()
    # CHD added Nov 2022
    allAvgDeltas <- c()


    # detect possible categorical moderator values
    for (i in 1:n.fitted.obj) {
      # CHD changes necessary because return list of ctmaFit changed in Sep 2022
      if (!(is.null(ctmaFitObject[[i]]$argumentList$mod.type))) ctmaFitObject[[i]]$mod.type <- ctmaFitObject[[i]]$argumentList$mod.type
      #
      if (!(is.null(ctmaFitObject[[i]]$mod.type))) {
        if (ctmaFitObject[[i]]$mod.type == "cat") {
          mod.values[[i]] <- c(-999, unique(as.numeric(substr(rownames(ctmaFitObject[[i]]$summary$mod.effects), 1,2))))
        }
      }
    }

    tmp <- length(mod.values); tmp
    for (i in 1:n.fitted.obj) {
      if (tmp != 0) {
        if (length(mod.values[[i]]) <= 1) mod.values[[i]] <- mod.values.backup; mod.values
      } else {
        mod.values[[i]] <- mod.values.backup; mod.values
      }
    }

    for (i in 1:n.fitted.obj) {
      #i <- 1
      if ("power" %in% plot.type[[i]]) plot.type[[i]] <- c("power", "drift")
      n.latent[i] <- unlist(ctmaFitObject[[i]]$n.latent); n.latent[i]
      driftNames[[i]] <- ctmaFitObject[[i]]$parameterNames$DRIFT; driftNames[[i]]
      n.studies[i] <- ctmaFitObject[[i]]$n.studies; n.studies[i]

      # CHD Jan 2023:
      #for (j in 1:n.fitted.obj)  tmp1 <- c(tmp1, unlist(lapply(ctmaFitObject[[i]]$studyList, function(extract) extract$originalStudyNo)))
      study.numbers[[i]] <- unlist(lapply(ctmaFitObject[[i]]$studyList, function(extract) extract$originalStudyNo)); study.numbers[[i]]
      tmp1 <- 0
      if (is.na(study.numbers[[i]][length(study.numbers[[i]])])) {
        study.numbers[[i]] <- study.numbers[[i]][1:(length(study.numbers[[i]])-1 )]
        tmp1 <- 1
      }

      n.primary.studies[i] <- length(ctmaFitObject[[i]]$studyList); n.primary.studies[i]
      if (tmp1 == 1) n.primary.studies[i] <- n.primary.studies[i] - 1

      if (n.studies[i] == 1) {
        DRIFTCoeff[[i]] <- list(ctmaFitObject[[i]]$modelResults$DRIFT); DRIFTCoeff[[i]]
        if (undoTimeScaling == TRUE) {
          if (!(is.null(ctmaFitObject[[i]]$modelResults$DRIFToriginal_time_scale))) {
            DRIFTCoeff[[i]] <- list(ctmaFitObject[[i]]$modelResults$DRIFToriginal_time_scale)
          }
        }
      } else {
        DRIFTCoeff[[i]] <- ctmaFitObject[[i]]$modelResults$DRIFT; DRIFTCoeff[[i]]
        if (undoTimeScaling == TRUE) {
          if (!(is.null(ctmaFitObject[[i]]$modelResults$DRIFToriginal_time_scale))) {
            DRIFTCoeff[[i]] <- ctmaFitObject[[i]]$modelResults$DRIFToriginal_time_scale
          }
        }
      }

      sampleSize[[i]] <- ctmaFitObject[[i]]$statisticsList$allSampleSizes; sampleSize[[i]]

      if ( ("funnel" %in% plot.type[[i]]) || ("forest" %in% plot.type[[i]]) ) {

        if (n.studies[i] == 1) {
          if (undoTimeScaling == TRUE) {
            DRIFTSE[[i]] <- list(ctmaFitObject[[i]]$modelResults$DRIFTSE); DRIFTSE[[i]]
          } else {
            DRIFTSE[[i]] <- list(ctmaFitObject[[i]]$modelResults$DRIFTSE_timeScaled); DRIFTSE[[i]]
          }
        } else {
          if (undoTimeScaling == TRUE) {
            DRIFTSE[[i]] <- ctmaFitObject[[i]]$modelResults$DRIFTSE; DRIFTSE[[i]]
          } else {
            DRIFTSE[[i]] <- ctmaFitObject[[i]]$modelResults$DRIFTSE_timeScaled; DRIFTSE[[i]]
          }
        }

        FixedEffect_Drift[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`[2,]; FixedEffect_Drift[[i]]
        FixedEffect_DriftLow[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftLowerLimit",]; FixedEffect_DriftLow[[i]]
        FixedEffect_DriftUp[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftUpperLimit",]; FixedEffect_DriftUp[[i]]
      }

      if ("drift" %in% plot.type[[i]]) {
        allDeltas[[i]] <- ctmaFitObject[[i]]$statisticsList$allDeltas; allDeltas[[i]]
        if (undoTimeScaling == FALSE) {
          if (!(is.null(ctmaFitObject[[i]]$summary$scaledTime)))  allDeltas[[i]] <- unlist(lapply(allDeltas[[i]], function(x) x * ctmaFitObject[[i]]$summary$scaledTime))
        }
        maxDelta[i] <- max(allDeltas[[i]], na.rm=TRUE); maxDelta[i]
        minDelta[i] <- min(allDeltas[[i]], na.rm=TRUE); minDelta[i]
        meanDelta[i] <- mean(allDeltas[[i]], na.rm=TRUE); meanDelta[i]
        # CHD added 4. Nov 2022
        if (n.studies[i] > 1) { # if init fit object
          #if (is.null(ctmaFitObject[[i]]$argumentList)) { # if init fit object
          allAvgDeltas[[i]] <- unlist(lapply(ctmaFitObject[[i]]$primaryStudyList$deltas, mean))
          #allAvgDeltas[[i]]
          tmp1 <- allAvgDeltas[[i]]
          # if deltas are not specified because raw data are provided
          deltaCounter <- 1
          for (j in 1:length(tmp1)) {
            tmp2 <- length(tmp1[[j]]); tmp2
            if (is.na(allAvgDeltas[[i]][j])) {
              allAvgDeltas[[i]][j] <- mean(allDeltas[[i]][deltaCounter:tmp2])
            }
            deltaCounter <- deltaCounter + tmp2
          }
        }
        #
        if (n.studies[i] == 1) { # if full fit object
          #if (!(is.null(ctmaFitObject[[i]]$argumentList))) { # if full fit object
          allAvgDeltas[[i]] <- mean(allDeltas[[i]])
        }
        #allAvgDeltas[[1]]
      }


      if ("power" %in% plot.type[[i]]) {
        requiredSampleSizes[[i]] <- ctmaFitObject[[i]]$summary$estimates$`Required Sample Sizes`
        statisticalPower <- ctmaFitObject[[i]]$summary$estimates$`Requested Statistical Power`
      }
    }
    nlatent <- unlist(n.latent[[1]]); nlatent  # nlatent used general specs; n.latent in special specs

  } ### END Extracting parameters

  #######################################################################################################################
  ################################################### funnel plots ######################################################
  #######################################################################################################################

  for (k in 1:n.fitted.obj) {
    if ("funnel" %in% unlist(plot.type[[k]])) {
      stretch <- 1.2
      figureName <- c()
      for (i in 1:(n.latent[k]^2)) {
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

        figureContent <- colnames(DRIFTCoeff)[i]
        if (is.null(figureContent)) figureContent <- colnames(DRIFTCoeff)[[i]]
        if (is.null(figureContent)) figureContent <- colnames(DRIFTCoeff[[1]])[i]; figureContent
        # Add labels and title
        graphics::title(main = paste0("Funnel Plot for the Effect of ", figureContent))
        figureName[i] <- paste0(activeDirectory, saveFilePrefix, "_",  "Funnel Plot for ", figureContent, ".png")
        grDevices::dev.copy(grDevices::png, figureName[i], width = 8, height = 8, units = 'in', res = 300)
        grDevices::dev.off()
      } # END for (i in 1:(n.latent^2))
    } # END if ("funnel" %in% plot.type)
  } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  #################################################### forest plots #####################################################
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

      # try
      yMin <- 0
      xMin <- 0
      if (is.null(ctmaFitObject[[k]]$xMin)) plot.xMin <- xMin else plot.xMin <- ctmaFitObject[[k]]$xMin; plot.xMin
      if (is.null(ctmaFitObject[[k]]$xMax)) plot.xMax <- xMax else plot.xMax <- ctmaFitObject[[k]]$xMax; plot.xMax
      if (is.null(ctmaFitObject[[k]]$yMin)) plot.yMin <- yMin else plot.yMin <- ctmaFitObject[[k]]$yMin; plot.yMin
      if (is.null(ctmaFitObject[[k]]$yMax)) plot.yMax <- yMax else plot.yMax <- ctmaFitObject[[k]]$yMax; plot.yMax


      # adaptations based on figures size ('base' objects contained transformed coefficients and are used for plotting)
      {
        heigthPerStudy <- yMax/(n.studies+1); heigthPerStudy
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
        # set frame by plotting invisible object
        plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(plot.xMin, plot.xMax), ylim = c(plot.yMin, plot.yMax),
             xaxt='n', yaxt='n', ann=FALSE)
        graphics::par(new=F)
        for (j in 1:n.studies) {
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
        ytop <- currentYPos - maxSquareSize/2 ; ytop
        ymiddle <- ytop - maxSquareSize/2; ymiddle
        ybottom <- ymiddle - maxSquareSize/2; ybottom
        x <- c(xmiddle-maxSquareSize/2, xmiddle, xmiddle+maxSquareSize/2, xmiddle, xmiddle-maxSquareSize/2); x # if based on N
        y <- c(ymiddle, ytop, ymiddle, ybottom, ymiddle); y
        graphics::polygon(x=x, y=y, col="black")
        graphics::abline(v=xmiddle, lty=2)

        # y-axis (original study numbers)
        atSeq <- seq(((n.studies+1)*heigthPerStudy- heigthPerStudy/2), 0, by = -heigthPerStudy); atSeq
        labelsSeq <- c(unlist(study.numbers), "SUM"); labelsSeq
        graphics::axis(side=2, at = atSeq, labels=labelsSeq, las=1)
        # x-axis (effect size)
        atSeq <- seq(plot.xMin, plot.xMax, by = 10); atSeq
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
        grDevices::dev.copy(grDevices::png, figureFileName, width = 8, height = 8, units = 'in', res = 300)
        grDevices::dev.off()
      } # END for (i in 1:(n.latent^2))
    } # END if ("forest" %in% plot.type)
  } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  ############################################## discrete time plots  ###################################################
  #######################################################################################################################

  if ("drift" %in% unlist(plot.type)) {

    if ( (plotCrossEffects == TRUE) || (plotAutoEffects == TRUE) ) {

      # Function to compute discrete parameters using drift parameters and time-scaling factors
      discreteDrift <-function(driftMatrix, timeScale, number) {
        discreteDriftValue <- OpenMx::expm(timeScale %x% driftMatrix)
        discreteDriftValue[number] }

      #######################################################################################################################
      ############## Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  ######################
      #######################################################################################################################

      # can only plot (overlay) multiple fitted objects if n.latent is identical
      if (length(unique(unlist(n.latent))) > 1) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        ErrorMsg <- "A fitted CoTiMA object has to be supplied to plot something. \nA fitted CoTiMA object has to be supplied to plot something. \nGood luck for the next try!"
        stop(ErrorMsg)
      }

      #######################################################################################################################
      ############################ Specification of Parameters for Plotting, Optimal Lags  ##################################
      #######################################################################################################################
      {
        ## timeRange, stepWidth, & noOfSteps
        maxDelta <- max(unlist(maxDelta)); maxDelta
        minDelta <- min(unlist(minDelta)); minDelta
        noDecims <- 0
        if (length(timeRange) < 1) {
          stepWidth <- 1
          usedTimeRange <- seq(1, 1.5*round(maxDelta + 1), stepWidth)
          # add empirical lags not yet included in requested timeRage
          # CHD changed 4 Nov 2022
          #usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
          usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas), unlist(allAvgDeltas)))); usedTimeRange
          noOfSteps <- length(usedTimeRange); noOfSteps
        }
        if (length(timeRange) > 0) {
          stepWidth <- timeRange[3]; stepWidth
          usedTimeRange <- seq(timeRange[1], timeRange[2], stepWidth); usedTimeRange
          # add empirical lags not yet included in requested timeRage
          # CHD changed 4 Nov 2022
          #tmp1 <- sort(unlist(allDeltas)/stepWidth); tmp1
          tmp1a <- sort(unlist(allDeltas)); tmp1a
          tmp1b <- sort(unlist(allAvgDeltas)); tmp1b
          tmp1c <- c()
          for (g in 1:n.fitted.obj) {
            tmp1c <- c(tmp1c, unlist(lapply(ctmaFitObject[[g]]$studyList, function(x) x$delta_t)))
            #CHD changed 10 Nov 2022
            if (length(ctmaFitObject[[g]]$studyList) > 1) {
              for (l in 1:length(ctmaFitObject[[g]]$studyList)) {
                tmp1c <- c(tmp1c, mean(ctmaFitObject[[g]]$studyList[[l]]$delta_t))
              }
            }
          }
          tmp1 <- sort(unique(c(tmp1a, tmp1b, tmp1c))); tmp1
          if (undoTimeScaling == FALSE) {
            tmp2 <- ctmaFitObject[[1]]$argumentList$scaleTime; tmp2
            if (is.null(tmp2)) tmp2 <- ctmaFitObject[[g]]$summary$scaledTime; tmp2
            tmp1 <- tmp1 * tmp2; tmp1
          }
          if (nchar(stepWidth) == 1) {
            noDecims <- 0
          } else {
            noDecims <- nchar((strsplit(as.character(stepWidth), "\\.")[[1]][2])); noDecims # number of decimal places
          }
          tmp1 <- round(tmp1, noDecims); tmp1
          #
          tmp2 <- which(tmp1 >= min(usedTimeRange) & tmp1 <= max(usedTimeRange) ); tmp2
          usedTimeRange <- sort(unique(c(usedTimeRange, tmp1[tmp2]))); usedTimeRange
          noOfSteps <- length(usedTimeRange); noOfSteps
        }

        # discrete effects across time range
        discreteDriftCoeff <- linearizedTIpredEffect <- DRIFThi <- DRIFTlo <- list()

        #g <- 1
        for (g in 1:n.fitted.obj) {
          toPlot <- n.studies[[g]]; toPlot

          ########################## start dealing with possible moderator values #############################################
          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {

            toPlot <- length(unlist(mod.values[[1]])); toPlot

            # augment usedTimeRange by time points (quantiles) where moderator values are plotted
            xValueForModValue <- stats::quantile(usedTimeRange, probs = seq(0, 1, 1/(toPlot+1))); xValueForModValue # used for positioning of moderator value in plot
            usedTimeRange <- unique(sort(c(xValueForModValue, usedTimeRange))); usedTimeRange # correcting for added time points
            noOfSteps <- length(usedTimeRange); noOfSteps

            DRIFTCoeff[[g]] <- linearizedTIpredEffect[[g]] <- DRIFThi[[g]] <- DRIFTlo[[g]] <- list()
            counter <- 1
            originalStudyNo <- delta_t <- c()

            for (i in unlist(mod.values[[g]]) ) {
              #i <- unlist(mod.values[[g]])[1]; i
              originalStudyNo[counter] <- i # used for labeling in plot
              delta_t[counter] <- xValueForModValue[counter+1]

              ### compute moderated drift matrices
              # main effects
              # CHD 11. Oct 203
              tmp1a <- which(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$matrix == "DRIFT"); tmp1a
              tmp1b <- which(!(is.na(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$param))); tmp1b
              tmp1c <- which(tmp1b %in% tmp1a); tmp1c
              #tmp1 <- ctmaFitObject[[g]]$studyFitList$stanfit$rawest[1:(n.latent^2)]; tmp1
              tmp1 <- ctmaFitObject[[g]]$studyFitList$stanfit$rawest[tmp1c]; tmp1
              tmp1 <- matrix(tmp1, n.latent[[g]], byrow=TRUE); tmp1

              # moderator effects (could be partial)
              if (n.primary.studies[[g]] > n.studies[[g]]) {
                n.TIpreds <- n.primary.studies[[g]]-1; n.TIpreds
              } else {
                n.TIpreds <- n.studies[[g]]-1; n.TIpreds
              }

              if (ctmaFitObject[[g]]$mod.type == "cont") {
                #ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT
                # could become necessary in newest ctsem version (look at ctsem 2022 workshop for getting moderated drift matrces)
                #tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT[,1:(n.latent^2),
                #                                                                                 ((n.TIpreds+1):(n.TIpreds+mod.number))]; tmp2
                # CHD 11. Oct 203
                tmp1a <- which(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$matrix == "DRIFT"); tmp1a
                tmp1b <- which(!(is.na(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$param))); tmp1b
                tmp1c <- which(tmp1b %in% tmp1a); tmp1c
                #tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,1:(n.latent^2),
                #                                                                                 ((n.TIpreds+1):(n.TIpreds+mod.number))]; tmp2
                tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,tmp1c,
                                                                                                 ((n.TIpreds+1):(n.TIpreds+mod.number))]; tmp2
                tmp3 <- rownames(ctmaFitObject[[g]]$modelResults$MOD); tmp3
                tmp4 <- c()
                for (l in 1:length(driftNames[[g]])) {
                  tmp5 <- grep(unlist(driftNames[[g]][l]), tmp3); tmp5
                  if (length(tmp5) == 0) tmp4 <- c(tmp4, NA) else tmp4 <- c(tmp4, tmp5)
                }
                tmp4[!(is.na(tmp4))] <- tmp2[!(is.na(tmp4))]
                tmp4[(is.na(tmp4))] <- 0

                tmp2 <- matrix(tmp4, n.latent[[g]], byrow=TRUE); tmp2 # raw moderator effect to be added to raw main effect (followed by tform)
                DRIFTCoeff[[g]][[counter]] <- tmp1 + unlist(mod.values[[g]])[counter] * tmp2; DRIFTCoeff[[g]][[counter]]
                names(DRIFTCoeff[[g]]) <- paste0("Moderator Value = ", mod.values, " SD from mean if standardized (default setting)")
                # compute matrices required  for linearizedTIpredEffect (JUST AS A CHECK)
                DRIFThi[[g]][[counter]] <- tmp1 + .01 * tmp2
                DRIFTlo[[g]][[counter]] <- tmp1 - .01 * tmp2
              }

              if (ctmaFitObject[[g]]$mod.type == "cat") {
                if (counter == 1) {
                  DRIFTCoeff[[g]][[counter]] <- tmp1 # copy main effects (= comparison group)
                  DRIFThi[[g]][[counter]] <- tmp1 #+ .01 * tmp2
                  DRIFTlo[[g]][[counter]] <- tmp1 #- .01 * tmp2
                } else {
                  n.mod.tmp <- length(mod.values[[g]])-1; n.mod.tmp
                  # added CHD 10 Nov 2022 to accound for possible modelling of means
                  tmp11 <- ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$matrix
                  tmp11 <- which(tmp11 == "MANIFESTMEANS"); tmp11
                  #is.na(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$param[tmp11])
                  tmp11 <- is.na(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$param[tmp11]); tmp11
                  offset <- length(which(tmp11 == FALSE)); offset
                  #if (!(is.null(ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull))) {
                  #  tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,1:(n.latent^2),
                  #                                                                                   n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  #} else {
                  #  tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT[,1:(n.latent^2),
                  #                                                                               n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  #}
                  # CHD changed 10 Nov 2022
                  if (!(is.null(ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull))) {
                    tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,(offset+1):(offset+n.latent^2),
                                                                                                     n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  } else {
                    tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT[,(offset+1):(offset+n.latent^2),
                                                                                                 n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  }

                  if (!(is.matrix(tmp2))) { # 29. Aug. 2022
                    tmp2 <- matrix(tmp2, n.latent, n.latent, byrow=TRUE); tmp2
                  } else {
                    tmp2 <- matrix(apply(tmp2, 2, mean), n.latent, n.latent, byrow=TRUE); tmp2 # new 18.7. 2022
                  }
                  DRIFTCoeff[[g]][[counter]] <- tmp1 + tmp2; DRIFTCoeff[[g]][[counter]]
                  # compute matrices required  for linearizedTIpredEffect (JUST AS A CHECK)
                  DRIFThi[[g]][[counter]] <- tmp1 + .01 * tmp2
                  DRIFTlo[[g]][[counter]] <- tmp1 - .01 * tmp2
                }
                names(DRIFTCoeff[[g]][[counter]]) <- paste0("Moderator Value = ", mod.values[[g]][counter])
              }
              counter <- counter +1
            } # END for (i in mod.values[[g]])

            ## apply tform to drift elements that should be tformed (extracted into tansforms)
            tmp1a <- ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars[, "transform"]; tmp1a
            tmp1b <- ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars[, "param"]; tmp1b
            transforms <- tmp1a[grep("toV", tmp1b)]; transforms

            for (k in 1:(length(DRIFTCoeff[[g]]))) {
              linearizedTIpredEffect[[g]][[k]] <- matrix(NA, n.latent, n.latent)
              counter <- 0
              for (l in 1:(n.latent)) {
                for (m in 1:(n.latent)) {
                  counter <- counter + 1
                  param <- DRIFTCoeff[[g]][[k]][l,m]; param
                  DRIFTCoeff[[g]][[k]][l,m] <- eval(parse(text=transforms[counter]))
                  # compute matrices required  for linearizedTIpredEffect (JUST AS A CHECK)
                  param <- DRIFThi[[g]][[k]][l,m]; param
                  DRIFThi[[g]][[k]][l,m] <- eval(parse(text=transforms[counter]))
                  param <- DRIFTlo[[g]][[k]][l,m]; param
                  DRIFTlo[[g]][[k]][l,m] <- eval(parse(text=transforms[counter]))
                  # undoTimScaling after tform
                  if ((undoTimeScaling == TRUE) & (!(is.null(ctmaFitObject[[g]]$summary$scaledTime)) ) ) {
                    DRIFTCoeff[[g]][[k]][l,m] <- DRIFTCoeff[[g]][[k]][l,m] * ctmaFitObject[[g]]$summary$scaledTime
                    DRIFThi[[g]][[k]][l,m] <- DRIFThi[[g]][[k]][l,m] * ctmaFitObject[[g]]$summary$scaledTime
                    DRIFTlo[[g]][[k]][l,m] <- DRIFTlo[[g]][[k]][l,m] * ctmaFitObject[[g]]$summary$scaledTime
                  }
                }
              }
              linearizedTIpredEffect[[g]][[k]] <- t((DRIFThi[[g]][[k]] - DRIFTlo[[g]][[k]])/.02) # transpose to make same order as in summary report of TIPREDeffects
            }
            # check all diags
            allDiags <- c()
            for (i in 1:length(DRIFTCoeff[[g]])) allDiags <- c(allDiags, diag(DRIFTCoeff[[g]][[i]]))
            if (!(all(allDiags < 0))) {
              if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
              ErrorMsg <- "Some of the moderated drift matrices have values > 0 in their diagonals. \nThis is likely if the model used to create \"ctmaFitObject\" was not identified! \nYou may want to try smaller moderator values (e.g., \"mod.values=c(-.5, 0., .5)\")! \nSome of the moderated drift matrices have values > 0 in their diagonals. \nThis is likely if the model used to create \"ctmaFitObject\" was not identified! \nYou may want to try smaller moderator values (e.g., \"mod.values=c(-.5, 0., .5)\")! \nGood luck for the next try!"
              stop(ErrorMsg)
            }
          } ########################## end dealing with possible moderator values #############################################

          if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) {
            discreteDriftCoeff[[g]] <- array(dim=c(n.studies[[g]], noOfSteps-1, n.latent[[g]]^2))
            for (h in 1:n.studies[g]) {
              for (i in usedTimeRange[1]:(noOfSteps-1)){
                timeValue <- i * stepWidth; timeValue
                discreteDriftCoeff[[g]][h, i, 1:(n.latent[[g]]^2)] <- c(discreteDrift(matrix(unlist(DRIFTCoeff[[g]][[h]]), n.latent[[g]], n.latent[[g]]), timeValue))
              }
            }
          }

          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
            discreteDriftCoeff[[g]] <- array(dim=c(length(mod.values[[g]]), noOfSteps-1, n.latent[[g]]^2))
            for (h in 1:length(mod.values[[g]])) {
              for (i in usedTimeRange[1]:(noOfSteps-1)) {
                timeValue <- i * stepWidth; timeValue
                discreteDriftCoeff[[g]][h, i, 1:(n.latent[[g]]^2)] <- c(discreteDrift(matrix(unlist(DRIFTCoeff[[g]][[h]]), n.latent[[g]], n.latent[[g]]), timeValue))
              } # end for (i in usedTimeRange[1]:(noOfSteps-1))
            } # end for (h in 1:length(mod.values[[g]])
          } # end if !(is.null(ctmaFitObject[[g]]$modelResults$MOD)))

        } # END for (g in 1:n.fitted.obj)

      } ### END Specification of Parameters for Plotting, Statistical Power, Optimal Lags ###


      Msg <- "      #################################################################################
      ################################### Plotting ####################################
      #################################################################################"
      message(Msg)

      ##################################### SELECT DRIFT MATRICES ########################################

      # Drift matrix used for plotting effects of primary studies (just renamed to adapt to older plotting procedures)
      DriftForPlot <- DRIFTCoeff; DriftForPlot

      ##################################### COMPUTE DOTS FOR PLOTTING ########################################
      plotPairs <- list()    # effect sizes and (scaled) time point
      dotPlotPairs <- list() # study symbol and (scaled) time point

      for (g in 1:n.fitted.obj) {
        #g <- 1
        if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[1]])
        plotPairs[[g]] <- array(dim=c(toPlot, length(usedTimeRange), 1+n.latent[[g]]^2))
        dotPlotPairs[[g]] <- array(dim=c(toPlot, length(usedTimeRange), 1+n.latent[[g]]^2))

        if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
          xValueForModValue2 <- xValueForModValue[-length(xValueForModValue)]; xValueForModValue2
          xValueForModValue2 <- xValueForModValue[-1]; xValueForModValue2
          # CHD 11. Oct 2023
          xValueForModValue2 <- xValueForModValue2[-length(xValueForModValue2)]; xValueForModValue2
        }

        for (h in 1:toPlot) {
          #h <- 5
          for (stepCounter in 1:length(usedTimeRange)){
            #stepCounter <- 64
            timeValue <- usedTimeRange[stepCounter]; timeValue
            plotPairs[[g]][h,stepCounter,1] <- timeValue; plotPairs[[g]][h,stepCounter,1]
            for (j in 1:(n.latent[[g]]^2)) {
              #j <- 2
              plotPairs[[g]][h,stepCounter,(1+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              #plotPairs[[g]][h,,]
              if (toPlot == 1) {
                tmp <- round(meanDelta[[1]],0)
              } else {
                if (exists("delta_t")) {
                  # CHD next line replaced by following two on 11 Oct 2022
                  #if (!(is.null(delta_t))) {
                  if (!(is.null(delta_t[h]))) {
                    if (!(is.na(delta_t[h])))  {
                      tmp <- delta_t[h]
                      #tmp
                    }
                  } else {
                    tmp <- mean(ctmaFitObject[[g]]$studyList[[h]]$delta_t)
                    if (undoTimeScaling == FALSE) {
                      if (!(is.null(ctmaFitObject[[g]]$summary$scaleTime))) {
                        tmp <- tmp * ctmaFitObject[[g]]$summary$scaleTime
                      }
                    }
                  }
                } else {
                  tmp <- mean(ctmaFitObject[[g]]$studyList[[h]]$delta_t); tmp
                  if (undoTimeScaling == FALSE) {
                    if (!(is.null(ctmaFitObject[[g]]$summary$scaleTime))) {
                      tmp <- tmp * ctmaFitObject[[g]]$summary$scaleTime
                    }
                  }
                }
              }
              # CHD added 5 Nov 2022
              tmp <- round(tmp, noDecims); tmp
              # CHD 11. Oct 2023
              #if (timeValue == tmp) { # plot only if the (used) time range includes the current study's mean time lag
              if ( (is.null(ctmaFitObject[[g]]$modelResults$MOD)) & (timeValue == tmp) )  { # set dots if NO moderator is plotted
                dotPlotPairs[[g]][h, stepCounter, 1] <- timeValue
                dotPlotPairs[[g]][h, stepCounter, (1+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              }
              if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD)))  { # set dots if moderator is plotted
                if (timeValue == xValueForModValue2[h]) {
                  tmp1 <- xValueForModValue[h+1]; tmp1
                  # CHD 11. Oct 2023
                  if (!(is.null(ctmaFitObject[[g]]$summary$scaleTime))) {
                    tmp1 <- tmp1 * ctmaFitObject[[g]]$summary$scaleTime
                  }
                  dotPlotPairs[[g]][h, stepCounter, 1] <- tmp1
                  dotPlotPairs[[g]][h, stepCounter, (1+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
                }
              }
            }
          } # END for (stepCounter in 0:noOfSteps)
          dotPlotPairs[[g]][h, ,]
        } # END for (h in 1:toPlot)
      } # END for (g in 1:n.fitted.obj)

      #dotPlotPairs[[1]][1,,]
      ##################################### PLOTTING PARAMETERS ##########################################
      {
        autoCols <- seq(1, nlatent^2, (nlatent+1)); autoCols
        crossCols <- (1:(nlatent^2))[!(1:(nlatent^2) %in% autoCols)]; crossCols
        yMinAuto <- yMinCross <-  999999
        yMaxAuto <- yMaxCross <- -999999
        for (g in 1:n.fitted.obj) {
          #g <- 1
          tmp1 <- dim(plotPairs[[g]])[3]; tmp1
          tmp2 <- plotPairs[[g]][, , -1, drop=FALSE]; tmp2 # array where in dim 3 there are n.latent dt effects sizes (do not drop if 1st dim=1)
          # y axis, auto
          yMinAutoTmp <- (min(tmp2[ , , autoCols])-.1); yMinAutoTmp
          if (yMinAutoTmp < yMinAuto) yMinAuto <- yMinAutoTmp; yMinAuto
          yMaxAutoTmp <- (max(tmp2[ , , autoCols])); yMaxAutoTmp
          if (yMaxAutoTmp > yMaxAuto) yMaxAuto <- yMaxAutoTmp; yMaxAuto
          # y axis, cross
          if (!(is.null(yLimitsForEffects))) {
            yMinCross <- yLimitsForEffects[1]
            yMaxCross <- yLimitsForEffects[2]
          } else {
            yMinCrossTmp <- round(min(tmp2[ , , crossCols]) - .1, 1); yMinCrossTmp
            if (yMinCrossTmp < yMinCross) yMinCross <- yMinCrossTmp; yMinCross
            yMaxCrossTmp <- round(max(tmp2[ , , crossCols]) + .1, 1); yMaxCrossTmp
            if (yMaxCrossTmp > yMaxCross) yMaxCross <- yMaxCrossTmp; yMaxCross
          }
        }
        # x axis,
        xMax <- max(usedTimeRange); xMax
        xMin <- usedTimeRange[1]; xMin
        #targetRows <- max(usedTimeRange)/stepWidth; targetRows
      }

      ############################################ PLOTTING ##############################################

      ## PLOT (auto effects)
      xLabelsBckup <- xLabels

      if (plotAutoEffects == TRUE) {
        graphics::plot.new()
        counter <- 0
        nlatent <- n.latent[[1]]; n.latent
        coeffSeq <- seq(1, nlatent^2, (nlatent+1)); coeffSeq

        for (j in coeffSeq) { # diagonal elements only
          #j <- 1
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            #g <- 1
            if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[1]])

            if (is.null(ctmaFitObject[[g]]$type)) plot..type <- "l" else plot..type <- ctmaFitObject[[g]]$type; plot..type
            if (is.null(ctmaFitObject[[g]]$col)) {
              plot.col <- "grey"
              if (toPlot == 1) plot.col <- "black"
            } else {
              plot.col <- ctmaFitObject[[g]]$col; plot.col
            }
            if (is.null(ctmaFitObject[[g]]$lwd)) {
              plot.lwd <- 1.5
              if (toPlot == 1) plot.lwd <- 2.5
            } else {
              plot.lwd <- ctmaFitObject[[g]]$lwd; plot.lwd
            }
            if (is.null(ctmaFitObject[[g]]$lty)) {
              plot.lty <- 1
              if (toPlot == 1) plot.lty <- 2
            } else {
              plot.lty <- ctmaFitObject[[g]]$lty; plot.lty
            }
            if (is.null(ctmaFitObject[[g]]$xMin)) plot.xMin <- xMin else plot.xMin <- ctmaFitObject[[g]]$xMin; plot.xMin
            if (is.null(ctmaFitObject[[g]]$xMax)) plot.xMax <- xMax else plot.xMax <- ctmaFitObject[[g]]$xMax; plot.xMax
            if (is.null(ctmaFitObject[[g]]$yMin)) plot.yMin <- yMinAuto else plot.yMin <- ctmaFitObject[[g]]$yMin; plot.yMin
            if (is.null(ctmaFitObject[[g]]$yMax)) plot.yMax <- yMaxAuto else plot.yMax <- ctmaFitObject[[g]]$yMax; plot.yMax
            if (is.null(ctmaFitObject[[g]]$dot.type)) dot.plot.type <- "b" else dot.plot.type <- ctmaFitObject[[g]]$dot.type; dot.plot.type
            if (is.null(ctmaFitObject[[g]]$dot.col)) dot.plot.col <- "black" else dot.plot.col <- ctmaFitObject[[g]]$dot.col; dot.plot.col
            if (is.null(ctmaFitObject[[g]]$dot.lwd)) dot.plot.lwd <- .5 else dot.plot.lwd <- ctmaFitObject[[g]]$dot.lwd; dot.plot.lwd
            if (is.null(ctmaFitObject[[g]]$dot.lty)) dot.plot.lty <- 3 else dot.plot.lty <- ctmaFitObject[[g]]$dot.lty; dot.plot.lty
            if (is.null(ctmaFitObject[[g]]$dot.pch)) dot.plot.pch <- 16 else dot.plot.pch <- ctmaFitObject[[g]]$dot.pch; dot.plot.pch
            if (is.null(ctmaFitObject[[g]]$dot.cex)) dot.plot.cex <- 3 else dot.plot.cex <- ctmaFitObject[[g]]$dot.cex; dot.plot.cex
            # CHD 5 Nov 2022
            if (is.null(ctmaFitObject[[g]]$text.cex)) text.plot.cex <- 2 else text.plot.cex <- ctmaFitObject[[g]]$text.cex; dot.plot.cex


            for (h in 1:toPlot) {
              #h <- 2
              currentPlotPair <- cbind(plotPairs[[g]][h, , 1], plotPairs[[g]][h, , 1+j])
              plot(currentPlotPair, type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax),
                   ylim = c(plot.yMin, plot.yMax),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFitObject[[g]]$plotStudyNo)) || (ctmaFitObject[[g]]$plotStudyNo==TRUE) ) {
                # black circle
                currentPlotPair <-cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                #currentPlotPair
                plot(currentPlotPair, type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     xlim = c(xMin, xMax), ylim = c(yMinAuto, yMaxAuto),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (toPlot > 1) {
                  if (exists("originalStudyNo")) {
                    if (!(is.null(originalStudyNo))) {
                      currentLabel <- originalStudyNo[h]
                    } else {
                      currentLabel <- ctmaFitObject[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                    }
                  } else {
                    currentLabel <- ctmaFitObject[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  }
                  if (currentLabel == -999) currentLabel <- "R"
                  if (is.null(currentLabel)) {
                    if (exists("originalStudyNo")) {
                      if (!(is.null(originalStudyNo))) {
                        currentLabel <- originalStudyNo[h]
                      } else {
                        currentLabel <- ctmaFitObject[[g]]$ctmaFitObject$studyList[[h]]$originalStudyNo; currentLabel
                      }
                    } else {
                      currentLabel <- ctmaFitObject[[g]]$ctmaFitObject$studyList[[h]]$originalStudyNo; currentLabel
                    }
                  }
                  if (nchar(currentLabel) == 1) graphics::text(currentPlotPair, labels=currentLabel, cex=3/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) == 2) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) > 2) graphics::text(currentPlotPair, labels=currentLabel, cex=1.5/5*text.plot.cex, col="white")
                } else {
                  currentLabel <- aggregateLabel
                  if (nchar(currentLabel) == 1) graphics::text(currentPlotPair, labels=currentLabel, cex=3/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) == 2) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) > 2) graphics::text(currentPlotPair, labels=currentLabel, cex=1.5/5*text.plot.cex, col="white")
                  currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                }
                graphics::par(new=T)
              }
            } # END for (h in 1:n.studies[[g]]) toPlot
          } # END for (g in 1:n.fitted.obj)

          # plot y-axis
          graphics::par(new=T)
          # CHD c(yMinCross, yMaxCross) changed to c(yMinAuto, yMaxAuto),
          plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(xMin, xMax), ylim = c(yMinAuto, yMaxAuto), xaxt='n', ann=FALSE, las=1)

          xLabels <- xLabelsBckup; xLabels
          if (is.null(xLabels)) xLabels <- round(seq(round(xMin,2), round((max(usedTimeRange+.4)),2), 1), 2); xLabels
          posForXLabel <- seq(xMin, xMax, abs(xMin-xMax)/(length(xLabels)-1)); posForXLabel
          if ( length(xLabels) != length(posForXLabel) ) posForXLabel <- posForXLabel[1:length(xLabels)]; xLabels; posForXLabel
          graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

          # add labels and title
          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
            graphics::title(main = paste0("Moderated Auto-regressive Effects of V", counter), sub = NULL,
                            xlab=paste0("Time Interval in ", timeUnit), ylab = "Auto-regressive Beta")
          } else {
            graphics::title(main = paste0("Auto-regressive Effects of V", counter), sub = NULL,
                            xlab=paste0("Time Interval in ", timeUnit), ylab = "Auto-regressive Beta")
          }

          # SAVE
          graphics::par(new=F)
          tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); tmp
          grDevices::dev.copy(grDevices::png, tmp, width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (j in coefSeq)
      } ## END PLOT (auto effects)


      ## PLOT (cross effects)
      if (plotCrossEffects == TRUE & nlatent > 1) {
        graphics::plot.new()
        counter <- 0
        nlatent <- n.latent[[1]]; n.latent
        coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
        for (j in coeffSeq) {
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[1]])
            #
            if (is.null(ctmaFitObject[[g]]$type)) plot..type <- "l" else plot..type <- ctmaFitObject[[g]]$type; plot..type
            if (is.null(ctmaFitObject[[g]]$col)) {
              plot.col <- "grey"
              if (n.studies[[g]] == 1) plot.col <- "black"
            } else {
              plot.col <- ctmaFitObject[[g]]$col; plot.col
            }
            if (is.null(ctmaFitObject[[g]]$lwd)) {
              plot.lwd <- 1.5
              if (n.studies[[g]] == 1) plot.lwd <- 2.5
            } else {
              plot.lwd <- ctmaFitObject[[g]]$lwd; plot.lwd
            }
            if (is.null(ctmaFitObject[[g]]$lty)) {
              plot.lty <- 1
              if (n.studies[[g]] == 1) plot.lty <- 2
            } else {
              plot.lty <- ctmaFitObject[[g]]$lty; plot.lty
            }
            if (is.null(ctmaFitObject[[g]]$xMin)) plot.xMin <- xMin else plot.xMin <- ctmaFitObject[[g]]$xMin; plot.xMin
            if (is.null(ctmaFitObject[[g]]$xMax)) plot.xMax <- xMax else plot.xMax <- ctmaFitObject[[g]]$xMax; plot.xMax
            if (is.null(ctmaFitObject[[g]]$yMin)) plot.yMin <- yMinCross else plot.yMin <- ctmaFitObject[[g]]$yMin; plot.yMin
            if (is.null(ctmaFitObject[[g]]$yMax)) plot.yMax <- yMaxCross else plot.yMax <- ctmaFitObject[[g]]$yMax; plot.yMax
            if (is.null(ctmaFitObject[[g]]$dot.type)) dot.plot.type <- "b" else dot.plot.type <- ctmaFitObject[[g]]$dot.type; dot.plot.type
            if (is.null(ctmaFitObject[[g]]$dot.col)) dot.plot.col <- "black" else dot.plot.col <- ctmaFitObject[[g]]$dot.col; dot.plot.col
            if (is.null(ctmaFitObject[[g]]$dot.lwd)) dot.plot.lwd <- .5 else dot.plot.lwd <- ctmaFitObject[[g]]$dot.lwd; dot.plot.lwd
            if (is.null(ctmaFitObject[[g]]$dot.lty)) dot.plot.lty <- 3 else dot.plot.lty <- ctmaFitObject[[g]]$dot.lty; dot.plot.lty
            if (is.null(ctmaFitObject[[g]]$dot.pch)) dot.plot.pch <- 16 else dot.plot.pch <- ctmaFitObject[[g]]$dot.pch; dot.plot.pch
            if (is.null(ctmaFitObject[[g]]$dot.cex)) dot.plot.cex <- 2 else dot.plot.cex <- ctmaFitObject[[g]]$dot.cex; dot.plot.cex
            # CHD 5 Nov 2022
            if (is.null(ctmaFitObject[[g]]$text.cex)) text.plot.cex <- 2 else text.plot.cex <- ctmaFitObject[[g]]$text.cex; dot.plot.cex

            #if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[1]])
            for (h in 1:toPlot) {
              currentPlotPair <- cbind(plotPairs[[g]][h, ,1], plotPairs[[g]][h, , 1+j])
              #currentPlotPair
              plot(currentPlotPair, type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax),
                   ylim = c(plot.yMin, plot.yMax),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFitObject[[g]]$plotStudyNo)) || (ctmaFitObject[[g]]$plotStudyNo==TRUE) ) {
                currentPlotPair <-cbind(dotPlotPairs[[g]][h, ,1], dotPlotPairs[[g]][h, ,1+j])
                plot(currentPlotPair, type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     xlim = c(plot.xMin, plot.xMax),
                     ylim = c(plot.yMin, plot.yMax),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (toPlot > 1) {
                  if (exists("originalStudyNo")) {
                    if (!(is.null(originalStudyNo))) {
                      currentLabel <- originalStudyNo[h]
                    } else {
                      currentLabel <- ctmaFitObject[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                    }
                  } else {
                    currentLabel <- ctmaFitObject[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  }
                  if (currentLabel == -999) currentLabel <- "R"
                  if (is.null(currentLabel)) {
                    if (exists("originalStudyNo")) {
                      if (!(is.null(originalStudyNo))) {
                        currentLabel <- originalStudyNo[h]
                      } else {
                        currentLabel <- ctmaFitObject[[g]]$ctmaFitObject$studyList[[h]]$originalStudyNo; currentLabel
                      }
                    } else {
                      currentLabel <- ctmaFitObject[[g]]$ctmaFitObject$studyList[[h]]$originalStudyNo; currentLabel
                    }
                  }

                  currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                  #if (h < 10) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  #if (h > 9) graphics::text(currentPlotPair, labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                  if (nchar(currentLabel) == 1) graphics::text(currentPlotPair, labels=currentLabel, cex=3/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) == 2) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) > 2) graphics::text(currentPlotPair, labels=currentLabel, cex=1.5/5*text.plot.cex, col="white")
                  #currentLabel <- aggregateLabel
                  #currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                  #graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                } else {
                  currentLabel <- aggregateLabel
                  if (nchar(currentLabel) == 1) graphics::text(currentPlotPair, labels=currentLabel, cex=3/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) == 2) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*text.plot.cex, col="white")
                  if (nchar(currentLabel) > 2) graphics::text(currentPlotPair, labels=currentLabel, cex=1.5/5*text.plot.cex, col="white")
                  currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                  #graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                }
                graphics::par(new=T)
              }
            }
          } # END for (g in 1:n.fitted.obj)

          # plot y-axis
          graphics::par(new=T)
          plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(xMin, xMax), ylim = c(yMinCross, yMaxCross), xaxt='n', ann=FALSE, las=1)

          xLabels <- xLabelsBckup; xLabels
          if (is.null(xLabels)) xLabels <- round(seq(round(xMin,2), round((max(usedTimeRange+.4)),2), 1), 2); xLabels
          posForXLabel <- seq(xMin, xMax, abs(xMin-xMax)/(length(xLabels)-1)); posForXLabel
          if ( length(xLabels) != length(posForXLabel) ) posForXLabel <- posForXLabel[1:length(xLabels)]; xLabels; posForXLabel
          graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

          # Add labels and title
          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
            driftNamesTmp <- c(t(matrix(driftNames[[1]], n.latent))); driftNamesTmp
          } else {
            driftNamesTmp <- driftNames[[1]]
          }

          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
            graphics::title(main = paste0("Moderated Cross-lagged Effects of ", driftNamesTmp[j]), sub = NULL,
                            xlab=paste0("Time Interval in ", timeUnit), ylab = "Cross-lagged Beta")
          } else {
            graphics::title(main = paste0("Cross-lagged Effects of ", driftNamesTmp[j]), sub = NULL,
                            xlab=paste0("Time Interval in ", timeUnit), ylab = "Cross-lagged Beta")
          }

          graphics::par(new=F)
          tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNamesTmp[j], ".png"); tmp
          grDevices::dev.copy(grDevices::png, tmp, width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (j in coeffSeq)

      } ## END PLOT (if (plotCrossEffects == TRUE & nlatent > 1))
    } ### END if (plotCrossEffects == TRUE | plotAutoEffects == TRUE)
  }  ## END if ("drift" %in% plot.type)


  #######################################################################################################################
  ########################################## required sample size plots  ################################################
  #######################################################################################################################

  if ("power" %in% unlist(plot.type)) {
    graphics::plot.new()
    g <- 1 # only a single power plot
    if (is.null(ctmaFitObject[[g]]$pow.type)) pow.plot.type <- "b" else pow.plot.type <- ctmaFitObject[[g]]$pow.type
    if (is.null(ctmaFitObject[[g]]$pow.col)) {
      pow.plot.col <- rep(c("black", "grey"), length(statisticalPower))
    } else {
      pow.plot.col <- ctmaFitObject[[g]]$pow.col
    }
    if (is.null(ctmaFitObject[[g]]$pow.lwd)) pow.plot.lwd <- .5 else pow.plot.lwd <- ctmaFitObject[[g]]$pow.lwd
    if (is.null(ctmaFitObject[[g]]$pow.lty)) pow.plot.lty <- 3 else pow.plot.lty <- ctmaFitObject[[g]]$pow.lty
    if (is.null(ctmaFitObject[[g]]$pow.yMin)) pow.plot.yMin <- 0 else pow.plot.yMin <- ctmaFitObject[[g]]$pow.yMin
    if (is.null(ctmaFitObject[[g]]$pow.yMax)) pow.plot.yMax <- 2000 else pow.plot.yMax <- ctmaFitObject[[g]]$pow.yMax

    tmp1 <- suppressWarnings(as.numeric(rownames(requiredSampleSizes[[g]]))); tmp1          # time lags
    tmp1 <- previouslyUsedTimeRange <- tmp1[!(is.na(tmp1))]; tmp1
    previousNoOfSteps <- length(previouslyUsedTimeRange); previousNoOfSteps
    tmp2 <- !(duplicated(tmp1)); tmp2
    currentRequiredSamleSizes <- requiredSampleSizes[[g]][tmp2,]; currentRequiredSamleSizes
    tmp4 <- nrow(currentRequiredSamleSizes); tmp4
    currentRequiredSamleSizes <- currentRequiredSamleSizes[-c((tmp4-2):tmp4),]; currentRequiredSamleSizes

    # adaptation of time range
    if (!(all(tmp1 %in% usedTimeRange) & all(usedTimeRange %in% tmp1))) {
      Msg <- paste0("      Note that required sample sizes can only be plotted for those time intervals
      that were provided in the \'timeRange\' argument of the ctmaPower function used before, and
      which was ", previouslyUsedTimeRange[1], ":", previouslyUsedTimeRange[2], ":", previouslyUsedTimeRange[3],
                    "..." , ":", previouslyUsedTimeRange[previousNoOfSteps-2], ":", previouslyUsedTimeRange[previousNoOfSteps-1],
                    ":", previouslyUsedTimeRange[previousNoOfSteps], ".
      => The \'timeRange\' argument was automatically adapted to the values used with ctmaPower <=
      You may need to re-run the ctmaPower function to suit your needs. \n")
      message(Msg)
    }

    if (min(tmp1) != min(usedTimeRange)) {
      if (min(tmp1) < min(usedTimeRange)) tmpText <- " longer " else tmpText <- " shorter "
      Msg <- paste0("      The shortest time interval provided in the \'timeRange\' argument
      of the current plot flunction is ", min(usedTimeRange), ", and it is", tmpText, "than the shortest time interval
      provided in \'timeRange\' argument of the ctmaPower function used before, which is ", min(tmp1), ".
      => The \'timeRange\' argument was automatically shortened <=
      You may need to re-run the ctmaPower function to suit your needs. \n")
      message(Msg)
    }
    if (max(tmp1) != max(usedTimeRange)) {
      if (max(tmp1) < max(usedTimeRange)) tmpText <- " longer " else tmpText <- " shorter "
      Msg <- paste0("      The longest time interval provided in the \'timeRange\' argument
      of the current plot flunction is ", max(usedTimeRange), ", and it is ", tmpText, "than the longest time interval
      provided in \'timeRange\' argument of the ctmaPower function used before, which was ", max(tmp1), ".
      => You might enlarge the \'timeRange\' argument of the current function, but this is no requirement <=
      You have to re-run the ctmaPower function if this does not suit your needs. \n")
      message(Msg)
    }
    tmp5 <- which(tmp1 >= min(usedTimeRange)); tmp5
    tmp6 <- which(tmp1 <= max(usedTimeRange)); tmp6
    usedTimeRange <- tmp1[tmp5 %in% tmp6]; usedTimeRange

    xMax <- max(usedTimeRange); xMax
    xMin <- usedTimeRange[1]; xMin
    tmp5 <- as.numeric(rownames(currentRequiredSamleSizes)); tmp5
    tmp6 <- which( (tmp5 >= xMin) & (tmp5 <=xMax) ); tmp6
    currentRequiredSamleSizes <- currentRequiredSamleSizes[tmp6 , ]; currentRequiredSamleSizes

    currentLWD <- c(3, 2, 1); currentLWD # line widths used later

    coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
    currentDriftNames <- driftNames[[1]][coeffSeq]; currentDriftNames
    for (j in 1:length(currentDriftNames)) {
      offset <- (j-1)*length(statisticalPower); offset
      graphics::par(new=F)
      for (h in 1:length(statisticalPower)) {
        forPlotting <- cbind(as.numeric(rownames(currentRequiredSamleSizes)),
                             currentRequiredSamleSizes[, h+offset]); forPlotting
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

      # plot y-axis (this is different compared to the discrete time plots)
      graphics::par(new=T)
      xLabels <- xLabelsBckup; xLabels
      if (is.null(xLabels)) xLabels <- sort(seq(max(usedTimeRange), min(usedTimeRange))); xLabels
      posForXLabel <- xLabels; posForXLabel
      if ( length(xLabels) != length(posForXLabel) ) posForXLabel <- posForXLabel[1:length(xLabels)]; xLabels; posForXLabel
      graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

      # legend
      graphics::legend('bottomright', legend=statisticalPower, lty=1, col=pow.plot.col, lwd=currentLWD, bty='n', cex=.75)

      tmp <- paste0(activeDirectory, saveFilePrefix," RequiredSampleSizesFor ", currentDriftNames[j], ".png"); tmp
      grDevices::dev.copy(grDevices::png, tmp, width = 8, height = 8, units = 'in', res = 300)
      grDevices::dev.off()
    }

  } ### END  ("power" %in% unlist(plot.type))

  if (!(is.null(ctmaFitObject[[1]]$modelResults$MOD))) {
    invisible(round(unlist(linearizedTIpredEffect)[-(1:(nlatent^2))], 4)) # just as a check
  }


} ### END function definition
