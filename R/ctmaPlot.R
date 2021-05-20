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
#' @param aggregateLabel label to indicate aggregated discrete time effects
#' @param xLabels labes used for x-axis
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
  mod.values=-2:2,
  aggregateLabel="",
  xLabels=NULL,
  ...
)
{  # begin function definition (until end of file)

  par.original <- par("new"); par.original
  on.exit(par(new=par.original))


  { # some checks

    # check if fit object is specified
    if (is.null(ctmaFitObject)){
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "A fitted ctma object has to be supplied to plot something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    # check #1 if object can be plotted
    if (class(ctmaFitObject) == "list") testObject <- ctmaFitObject[[1]] else testObject <- ctmaFitObject
    if (class(testObject) != "CoTiMAFit")  {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
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
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "This is nothing CoTiMA-related that I can plot. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    n.fitted.obj <- length(ctmaFitObject); n.fitted.obj # has to be done twice

    plot.type <- list() # has to be a list because a single fit could be used for different plots (e.g. "power")
    tmp <- activeDirectory; tmp
    if (n.fitted.obj == 1) {
      plot.type[[1]] <- ctmaFitObject[[1]]$plot.type; plot.type[[1]]
      if (is.null(tmp)) {
        activeDirectory <- ctmaFitObject[[1]]$activeDirectory
      } else {
        activeDirectory <- tmp
      }
    } else {
      for (i in 1:n.fitted.obj) {
        if (is.null(ctmaFitObject[[1]]$plot.type)) {
          plot.type[[i]] <- "drift"
        } else {
          plot.type[[i]] <- ctmaFitObject[[i]]$plot.type
        }
      }
      if (is.null(tmp)) activeDirectory <- ctmaFitObject$studyFitList[[1]]$activeDirectory else activeDirectory <- tmp
    }

    # check if fit object can be plotted
    if (length(unlist(plot.type))==0) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- "The fitted CoTiMA object provided cannot be plotted. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (length(unique(unlist(plot.type))) > 1) {
      if (any(!(unique(unlist(plot.type)) %in% c("funnel", "forest")))) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        ErrorMsg <- paste0("All fitted CoTiMA object to plot have to be of the same plot.type. \nThe following ploty.type arguments were found: \n", unique(unlist(plot.type)), "\nGood luck for the next try!")
        stop(ErrorMsg)
      }
    }

    if (length(unique(unlist(activeDirectory))) > 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      ErrorMsg <- paste0("More than a single active directors was sepcified, which does not work. \nThe following activeDirectory arguments were found: \n", unique(activeDirectory), "\nGood luck for the next try!")
      stop(ErrorMsg)
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

    # detect possible categorical moderator values
    for (i in 1:n.fitted.obj) {
      if (!(is.null(ctmaFitObject[[i]]$mod.type))) {
        if (ctmaFitObject[[i]]$mod.type == "cat") {
          mod.values[[i]] <- c(1, unique(as.numeric(substr(rownames(ctmaFitObject[[i]]$summary$mod.effects), 1,2))))
        }
      }
    }

    tmp <- length(mod.values); tmp
    for (i in 1:n.fitted.obj) {
      #i <- 2
      if (tmp != 0) {
        if (length(mod.values[[i]]) <= 1) mod.values[[i]] <- mod.values.backup; mod.values
      } else {
        mod.values[[i]] <- mod.values.backup; mod.values
      }
    }

    for (i in 1:n.fitted.obj) {
      if ("power" %in% plot.type[[i]]) plot.type[[i]] <- c("power", "drift")
      n.latent[i] <- unlist(ctmaFitObject[[i]]$n.latent); n.latent[i]
      driftNames[[i]] <- ctmaFitObject[[i]]$parameterNames$DRIFT; driftNames[[i]]
      n.studies[i] <- ctmaFitObject[[i]]$n.studies; n.studies[i]
      study.numbers[[i]] <- unlist(lapply(ctmaFitObject[[i]]$studyList, function(extract) extract$originalStudyNo)); study.numbers[[i]]
      if (n.studies[i] == 1) {
        DRIFTCoeff[[i]] <- list(ctmaFitObject[[i]]$modelResults$DRIFT); DRIFTCoeff[[i]]
      } else {
        DRIFTCoeff[[i]] <- ctmaFitObject[[i]]$modelResults$DRIFT; DRIFTCoeff[[i]]
      }
      sampleSize[[i]] <- ctmaFitObject[[i]]$statisticsList$allSampleSizes; sampleSize[[i]]

      if ( ("funnel" %in% plot.type[[i]]) || ("forest" %in% plot.type[[i]]) ) {
        if (n.studies[i] == 1) {
          DRIFTSE[[i]] <- list(ctmaFitObject[[i]]$modelResults$DRIFTSE); DRIFTSE[[i]]
        } else {
          DRIFTSE[[i]] <- ctmaFitObject[[i]]$modelResults$DRIFTSE; DRIFTSE[[i]]
        }
        FixedEffect_Drift[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`[2,]; FixedEffect_Drift[[i]]
        FixedEffect_DriftLow[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftLowerLimit",]; FixedEffect_DriftLow[[i]]
        FixedEffect_DriftUp[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftUpperLimit",]; FixedEffect_DriftUp[[i]]
      }
      if ("drift" %in% plot.type[[i]]) {
        allDeltas[[i]] <- ctmaFitObject[[i]]$statisticsList$allDeltas; allDeltas[[i]]
        maxDelta[i] <- max(allDeltas[[i]], na.rm=TRUE); maxDelta[i]
        minDelta[i] <- min(allDeltas[[i]], na.rm=TRUE); minDelta[i]
        meanDelta[i] <- mean(allDeltas[[i]], na.rm=TRUE); meanDelta[i]
      }
      if ("power" %in% plot.type[[i]]) {
        requiredSampleSizes[[i]] <- ctmaFitObject[[i]]$summary$estimates$`Required Sample Sizes`
        statisticalPower <- ctmaFitObject[[i]]$summary$estimates$`Requested Statistical Power`
      }
    }
    nlatent <- unlist(n.latent[[1]]); nlatent  # nlatent used general specs; n.latent in special specs
  }

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
        #i <- 3
        # set frame by plotting invisible object
        plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(plot.xMin, plot.xMax), ylim = c(plot.yMin, plot.yMax),
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
        x <- c(xmiddle-maxSquareSize/2, xmiddle, xmiddle+maxSquareSize/2, xmiddle, xmiddle-maxSquareSize/2); x # if based on N
        y <- c(ymiddle, ytop, ymiddle, ybottom, ymiddle); y
        graphics::polygon(x=x, y=y, col="black")
        graphics::abline(v=xmiddle, lty=2)

        # y-axis (original study numbers)
        atSeq <- seq(((n.studies+1)*heigthPerStudy- heigthPerStudy/2), 0, by = -heigthPerStudy); atSeq
        labelsSeq <- c(unlist(study.numbers), "SUM"); labelsSeq
        graphics::axis(side=2, at = atSeq, labels=labelsSeq, las=1)
        graphics::axis(side=2, at = atSeq, labels=labelsSeq, las=1)
        # x-axis (effect size)
        #atSeq <- seq(0, xMax, by = 10); atSeq
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
        if (length(timeRange) < 1) {
          stepWidth <- 1
          usedTimeRange <- seq(1, 1.5*maxDelta, stepWidth)
          # add empirical lags not yet included in requested timeRage
          usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
          noOfSteps <- length(usedTimeRange); noOfSteps
        } else {
          stepWidth <- timeRange[3]
          usedTimeRange <- seq(timeRange[1], timeRange[2], stepWidth)
          if (stepWidth != 1) timeUnit <- paste0(timeUnit, " * ", stepWidth); timeUnit
          # add empirical lags not yet included in requested timeRage
          tmp1 <- sort(unlist(allDeltas)/stepWidth); tmp1
          tmp2 <- which(tmp1 >= min(usedTimeRange) & tmp1 <= max(usedTimeRange) ); tmp2
          usedTimeRange <- sort(unique(c(usedTimeRange, tmp1[tmp2])))
          noOfSteps <- length(usedTimeRange); noOfSteps
        }

        # discrete effects across time range
        discreteDriftCoeff <- list()
        for (g in 1:n.fitted.obj) {
          #g <- 1
          ########################## start dealing with possible moderator values #############################################
          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
            if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(unlist(mod.values[[1]]))
            n.mod <- ctmaFitObject[[g]]$n.moderators; n.mod

            # augement usedTimeRange by time points (quantiles) where moderator values are plotted
            xValueForModValue <- stats::quantile(usedTimeRange, probs = seq(0, 1, 1/(toPlot+1))); xValueForModValue # used for positioning of moderator value in plot
            usedTimeRange <- unique(sort(c(xValueForModValue, usedTimeRange))); usedTimeRange # correcting for added time points
            noOfSteps <- length(usedTimeRange); noOfSteps

            DRIFTCoeff[[g]] <- list()
            counter <- 1

            for (i in mod.values[[g]]) {
              #i <- -2
              ctmaFitObject[[g]]$studyList[[counter]]$originalStudyNo <- i # used for labeling in plot
              ctmaFitObject[[g]]$studyList[[counter]]$delta_t <- xValueForModValue[counter+1]; xValueForModValue[counter+1]
              ### compute moderated drift matrices
              # main effects
              tmp1 <- ctmaFitObject[[g]]$modelResults$DRIFT; tmp1
              tmp1 <- matrix(tmp1, n.latent[[g]], byrow=TRUE); tmp1 # main effect
              # moderator effects (could be partial)
              tmp2 <- ctmaFitObject[[g]]$modelResults$MOD[,1]; tmp2
              tmp3 <- rownames(ctmaFitObject[[g]]$modelResults$MOD); tmp3
              tmp4 <- c()
              for (l in 1:length(driftNames[[g]])) {
                tmp5 <- grep(unlist(driftNames[[g]][l]), tmp3); tmp4
                if (length(tmp5) == 0) tmp4 <- c(tmp4, NA) else tmp4 <- c(tmp4, tmp5)
              }
              tmp4[!(is.na(tmp4))] <- tmp2
              tmp4[(is.na(tmp4))] <- 0
              tmp2 <- matrix(tmp4, n.latent[[g]], byrow=TRUE); tmp2 # moderator effect to be added to main effect

              if (ctmaFitObject[[1]]$mod.type == "cont") {
                DRIFTCoeff[[g]][[counter]] <- tmp1 + unlist(mod.values[[g]])[counter] * tmp2; DRIFTCoeff[[g]][[counter]]
                names(DRIFTCoeff[[g]]) <- paste0("Moderator Value = ", mod.values, " SD from mean if standardized (default setting)")
              }
              if (ctmaFitObject[[1]]$mod.type == "cat") {
                if (i == 1) {
                  DRIFTCoeff[[g]][[counter]] <- tmp1 # copy main effects (= comparison group)
                } else {
                  tmp2 <- ctmaFitObject[[g]]$summary$mod.effects[,1]; tmp2
                  #tmp2 <- tmp2[((i - 2) * n.latent^2 + 1): ((i - 2) * n.latent^2 + 0 + n.latent)]; tmp2
                  tmp2 <- tmp2[((i - 2) * n.latent^2 + 1): ((i - 2) * n.latent^2 + 0 + n.latent^2)]; tmp2
                  tmp2 <- matrix(tmp2, n.latent, n.latent, byrow=TRUE); tmp2
                  DRIFTCoeff[[g]][[counter]] <- tmp1 + (i-1) * tmp2; DRIFTCoeff[[g]][[counter]]
                }
                names(DRIFTCoeff[[g]][[counter]]) <- paste0("Moderator Value = ", mod.values[[g]][counter])
              }
              counter <- counter +1
            }
            #DRIFTCoeff
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
            } # end for (g in 1:n.fitted.obj)
          } else {
            discreteDriftCoeff[[g]] <- array(dim=c(length(mod.values[[g]]), noOfSteps-1, n.latent[[g]]^2))
            for (h in 1:length(mod.values[[g]])) {
              for (i in usedTimeRange[1]:(noOfSteps-1)){
                timeValue <- i * stepWidth; timeValue
                discreteDriftCoeff[[g]][h, i, 1:(n.latent[[g]]^2)] <- c(discreteDrift(matrix(unlist(DRIFTCoeff[[g]][[h]]), n.latent[[g]], n.latent[[g]]), timeValue))
              }
            } # end for (g in 1:n.fitted.obj)
          } # end if is.null(ctmaFitObject[[g]]$modelResults$MOD)) ... else ...

        } # end computing discrete effects across time range
      } ### END Specification of Parameters for Plotting, Statistical Power, Optimal Lags ###

      Msg <- "################################################################################# \n################################### Plotting #################################### \n#################################################################################"
      message(Msg)

      ##################################### SELECT DRIFT MATRICES ########################################

      # Drift matrix used for plotting effects of primary studies (just renamed to adapt to older plotting procedures)
      DriftForPlot <- DRIFTCoeff; DriftForPlot

      ##################################### COMPUTE DOTS FOR PLOTTING ########################################
      plotPairs <- list()    # effect sizes and (scaled) time point
      dotPlotPairs <- list() # study symbol and (scaled) time point

      for (g in 1:n.fitted.obj) {
        if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[1]])
        plotPairs[[g]] <- array(dim=c(toPlot, length(usedTimeRange), 1+n.latent[[g]]^2))
        dotPlotPairs[[g]] <- array(dim=c(toPlot, length(usedTimeRange), 1+n.latent[[g]]^2))

        for (h in 1:toPlot) {
          for (stepCounter in 1:length(usedTimeRange)){
            timeValue <- usedTimeRange[stepCounter]; timeValue
            plotPairs[[g]][h,stepCounter,1] <- timeValue; plotPairs[[g]][h,stepCounter,1]
            for (j in 1:(n.latent[[g]]^2)) {
              plotPairs[[g]][h,stepCounter,(1+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              plotPairs[[g]][h,stepCounter,(1+j)]
              if (toPlot == 1) tmp <- round(meanDelta[[1]],0) else tmp <- mean(ctmaFitObject[[g]]$studyList[[h]]$delta_t)
              tmp
              timeValue
              if (timeValue %in% (tmp / stepWidth)) { # plot only if the (used) time range includes the current study's mean time lag
                dotPlotPairs[[g]][h, stepCounter, 1] <- timeValue
                dotPlotPairs[[g]][h, stepCounter, (1+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              }
            }
          } # END for (stepCounter in 0:noOfSteps)
        } # END for (h in 1:toPlot)
      } # END for (g in 1:n.fitted.obj)


      ##################################### PLOTTING PARAMETERS ##########################################
      {
        autoCols <- seq(1, nlatent^2, (nlatent+1)); autoCols
        crossCols <- (1:(nlatent^2))[!(1:(nlatent^2) %in% autoCols)]; crossCols
        yMinAuto <- yMinCross <-  999999
        yMaxAuto <- yMaxCross <- -999999
        for (g in 1:n.fitted.obj) {
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
        targetRows <- max(usedTimeRange)/stepWidth; targetRows
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
            if (is.null(ctmaFitObject[[g]]$dot.cex)) dot.plot.cex <- 2 else dot.plot.cex <- ctmaFitObject[[g]]$dot.cex; dot.plot.cex

            for (h in 1:toPlot) {
              currentPlotPair <- cbind(plotPairs[[g]][h, , 1], plotPairs[[g]][h, , 1+j])
              plot(currentPlotPair, type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax),
                   ylim = c(plot.yMin, plot.yMax),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFitObject[[g]]$plotStudyNo)) || (ctmaFitObject[[g]]$plotStudyNo==TRUE) ) {
                # black circle
                currentPlotPair <-cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                plot(currentPlotPair, type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     xlim = c(xMin, xMax), ylim = c(yMinAuto, yMaxAuto),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (toPlot > 1) {
                  currentLabel <- ctmaFitObject[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  if (is.null(currentLabel)) currentLabel <- ctmaFitObject[[g]]$ctmaFitObject$studyList[[h]]$originalStudyNo; currentLabel
                  if (h < 10) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  if (h > 9) graphics::text(currentPlotPair, labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                } else {
                  currentLabel <- aggregateLabel
                  currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                  graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  currentPlotPair
                }
                graphics::par(new=T)
              }
            } # END for (h in 1:n.studies[[g]]) toPlot
          } # END for (g in 1:n.fitted.obj)

          # plot y-axis
          graphics::par(new=T)
          plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(xMin, xMax), ylim = c(yMinAuto, yMaxAuto), xaxt='n', ann=FALSE, las=1)

          xLabels <- xLabelsBckup; xLabels
          if (is.null(xLabels)) xLabels <- round(seq(round(xMin,2), round((max(usedTimeRange+.4)),2), 1), 2); xLabels
          posForXLabel <- (seq(1, noOfSteps, noOfSteps/length(xLabels))*stepWidth); posForXLabel
          if ( length(xLabels) < length(posForXLabel) ) xLabels <- xLabels[round(posForXLabel, 0)]; xLabels

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
          #j <- 3
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            #g <- 1
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

            #if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[1]])
            for (h in 1:toPlot) {
              #h <- 1
              currentPlotPair <- cbind(plotPairs[[g]][h, ,1], plotPairs[[g]][h, , 1+j])
              currentPlotPair
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
                  currentLabel <- ctmaFitObject[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  if (is.null(currentLabel)) currentLabel <- ctmaFitObject[[g]]$ctmaFitObject$studyList[[h]]$originalStudyNo; currentLabel
                  currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                  if (h < 10) graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  if (h > 9) graphics::text(currentPlotPair, labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                } else {
                  currentLabel <- aggregateLabel
                  currentPlotPair <- cbind(dotPlotPairs[[g]][h, ,1], plotPairs[[g]][h, ,1+j])
                  graphics::text(currentPlotPair, labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
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
          posForXLabel <- (seq(1, noOfSteps, noOfSteps/length(xLabels))*stepWidth); posForXLabel
          if ( length(xLabels) < length(posForXLabel) ) xLabels <- xLabels[round(posForXLabel, 0)]; xLabels

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
          #tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNamesTmp[[g]][j], ".png"); tmp
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
    if (is.null(ctmaFitObject[[g]]$pow.yMin)) pow.plot.yMin <- 0 else pow.plot.lty <- ctmaFitObject[[g]]$pow.yMin
    if (is.null(ctmaFitObject[[g]]$pow.yMin)) pow.plot.yMax <- 2000 else pow.plot.lty <- ctmaFitObject[[g]]$pow.yMax

    #requiredSampleSizes
    tmp1 <- suppressWarnings(as.numeric(rownames(requiredSampleSizes[[g]]))); tmp1
    tmp1 <- tmp1[!(is.na(tmp1))]; tmp1
    tmp1 <- round(tmp1, 0); tmp1
    tmp2 <- !(duplicated(tmp1)); tmp2
    currentRequiredSamleSizes <- requiredSampleSizes[[g]][tmp2,]; currentRequiredSamleSizes
    tmp4 <- nrow(currentRequiredSamleSizes); tmp4
    currentRequiredSamleSizes <- currentRequiredSamleSizes[-c((tmp4-2):tmp4),]; currentRequiredSamleSizes

    xMax <- max(usedTimeRange); xMax
    xMin <- usedTimeRange[1]; xMin
    currentRequiredSamleSizes <- currentRequiredSamleSizes[xMin:xMax , ]; currentRequiredSamleSizes
    #usedTimeRange <- min(tmp1):max(tmp1); usedTimeRange
    #xMax <- round(max(usedTimeRange), 0); xMax
    #xMin <- usedTimeRange[1]; xMin

    currentLWD <- c(3, 2, 1); currentLWD # line widths used later

    coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
    currentDriftNames <- driftNames[[1]][coeffSeq]; currentDriftNames
    for (j in 1:length(currentDriftNames)) {
      #j <- 1
      offset <- (j-1)*length(statisticalPower); offset
      graphics::par(new=F)
      for (h in 1:length(statisticalPower)) {
        #h <- 1
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

      # plot y-axis
      graphics::par(new=T)
      xLabels <- xLabelsBckup; xLabels
      if (is.null(xLabels)) xLabels <- round(seq(round(xMin,2), round((max(usedTimeRange)+1),2), 1), 2); xLabels
      posForXLabel <- (seq(1, noOfSteps, noOfSteps/length(xLabels))*stepWidth); posForXLabel
      if ( length(xLabels) < length(posForXLabel) ) xLabels <- xLabels[round(posForXLabel, 0)]; xLabels

      graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

      #if ( length(xLabels) > length(posForXLabel) ) {
      #  xLabels <- xLabels[1:length(posForXLabel)]
      #}

      graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

      # legend
      graphics::legend('bottomright', legend=statisticalPower, lty=1, col=pow.plot.col, lwd=currentLWD, bty='n', cex=.75)

      tmp <- paste0(activeDirectory, saveFilePrefix," RequiredSampleSizesFor ", currentDriftNames[j], ".png"); tmp
      grDevices::dev.copy(grDevices::png, tmp, width = 8, height = 8, units = 'in', res = 300)
      grDevices::dev.off()
    }

  } ### END Plotting ###

} ### END function definition
