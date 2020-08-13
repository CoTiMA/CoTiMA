#######################################################################################################################
################################################## Plotting ###########################################################
#######################################################################################################################
# debug <- 0
# if (debug == 1) {
#   # Primary Study Fits
#   ctmaFit = CoTiMAPowerFit #ctmaFit
#   ##### ENTER DEBUG INFORMATION BELOW THE FOLLOWING #######
#   # Directory names and file names
#   activeDirectory=NULL
#   #sourceDirectory= NULL
#   resultsFilePrefix="ctmaPlot"           # the prefix for the result file that is created
#   saveFilePrefix="ctmaPlot"
#   loadFilePrefix=NULL
#   # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
#   activateRPB=FALSE                      #set to TRUE to receive push messages with CoTiMA notifications on your phone
#   # Figure Parameters
#   plotCrossEffects=TRUE                  # plot cross effects across a range of discrete intervals
#   plotAutoEffects=TRUE                   # plot auto effects across a range of discrete intervals
#   timeUnit="timeUnit"                    # timelag unit to lable x-axis of plots
#   timeRange=c()                          # used for Plotting and Poweranalysis. Populated by 0 to 1.5*maxDelta (Steptwidth = 1) if not specified as c(min,max,stepWidth)
#   yLimitsForEffects=c()                  # used the y-axis of Drift-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
#   ##### ENTER DEBUG INFO HERE #######
#   #ctmaFit=ctmaFit
# }

#' ctmaPlot
#'
#' @param ctmaFit ?
#' @param activeDirectory ?
#' @param resultsFilePrefix ?
#' @param saveFilePrefix ?
#' @param activateRPB ?
#' @param plotCrossEffects ?
#' @param plotAutoEffects ?
#' @param timeUnit ?
#' @param timeRange ?
#' @param yLimitsForEffects ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaPlot <- function(
  # Primary Study Fits
  ctmaFit=list(),                    #list of lists: could be more than one fit object

  # Directory names and file names
  activeDirectory=NULL,
  #sourceDirectory= NULL,
  resultsFilePrefix="ctmaPlot",           # the prefix for the result file that is created
  saveFilePrefix="ctmaPlot",

  # Workflow (receive messages and request inspection checks to avoid proceeding with non admissible in-between results)
  activateRPB=FALSE,                      #set to TRUE to receive push messages with CoTiMA notifications on your phone

  # Figure Parameters
  plotCrossEffects=TRUE,                  # plot cross effects across a range of discrete intervals
  plotAutoEffects=TRUE,                   # plot auto effects across a range of discrete intervals
  timeUnit="timeUnit",                    # timelag unit to lable x-axis of plots
  timeRange=c(),                          # used for Plotting and Poweranalysis. Populated by 0 to 1.5*maxDelta (Steptwidth = 1) if not specified as c(min,max,stepWidth)
  yLimitsForEffects=c()                  # used the y-axis of Drift-Plots. Populated by c("round(min(effects)-.05, 1)", "round(max(effects)-.05, 1)")
)
{  # begin function definition (until end of file)

  #library('crayon')
  #library('ctsem')

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
    if (length(tmp) >= 13) {
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
    activeDirectory <- c()
    for (i in 1:n.fitted.obj) {
      plot.type[[i]] <- ctmaFit[[i]]$plot.type
      activeDirectory[i] <- ctmaFit[[i]]$activeDirectory
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

    for (i in 1:n.fitted.obj) {
      #i <- 1
      if (plot.type[[i]] == "power") plot.type[[i]] <- c("power", "drift")

      n.latent[i] <- unlist(ctmaFit[[i]]$n.latent); n.latent[i]
      driftNames[[i]] <- ctmaFit[[i]]$parameterNames$DRIFT; driftNames[[i]]
      n.studies[i] <- ctmaFit[[i]]$n.studies; n.studies[i]
      study.numbers[[i]] <- unlist(lapply(ctmaFit[[i]]$studyList, function(extract) extract$originalStudyNo)); study.numbers[[i]]
      study.numbers[[i]] <- study.numbers[[i]][-length(study.numbers[[i]])]; study.numbers[[i]]
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
        maxDelta[i] <- max(allDeltas[[i]]); maxDelta[i]
        minDelta[i] <- min(allDeltas[[i]]); minDelta[i]
        meanDelta[i] <- mean(allDeltas[[i]]); meanDelta[i]
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
          graphics::plot(pairs, xlab=colnames(DRIFTCoeff[[k]])[i], ylab="Standard Error", xlim = c(xMin2, xMax2), ylim = c(yMax2, 0))
          graphics::polygon(c(lowLeft, avgX, lowRight), c(stretch * yMax2, 0, stretch * yMax2), lty=3)
          graphics::abline(v=avgX, col="black", lwd=1.5, lty=1)

          figureContent <- colnames(DRIFTCoeff)[i]
          figureName[i] <- paste0(activeDirectory, saveFilePrefix, "_",  "Funnel Plot for ", figureContent, ".png")
          grDevices::dev.copy(png, figureName[i], width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (i in 1:(n.latent^2))
      } # END if ("funnel" %in% plot.type)
    } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  ################################################### FORREST plots #####################################################
  #######################################################################################################################
  for (k in 1:n.fitted.obj) {
    #k <- 1
    #("forest" %in% plot.type[[k]])
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
      #FixedEffect_Drift <- currentFit$summary$estimates$`Fixed Effects of Drift Coefficients`[2,]; FixedEffect_Drift
      #FixedEffect_DriftLow <- currentFit$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftLowerLimit",]; FixedEffect_DriftLow
      #FixedEffect_DriftUp <- currentFit$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_DriftUpperLimit",]; FixedEffect_DriftUp
      # precision (could be used for size of squares)
      #precision <- 1/(DRIFTSE^1); precision
      #minPrecision <- min(precision); minPrecision
      #maxPrecision <- max(precision); maxPrecision
      # sample sizes (could be used for size of squares)
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
      #minDRIFTCoeffLow * x + c = 0
      #maxDRIFTCoeffUp * x + c = xMax
      #(minDRIFTCoeffLow - maxDRIFTCoeffUp) * x = - xMax
      #(minDRIFTCoeffLow - maxDRIFTCoeffUp) * x = - xMax
      #(minDRIFTCoeffLow - maxDRIFTCoeffUp) / (-xMax) = 1/x
      # -(xMax)/(minDRIFTCoeffLow - maxDRIFTCoeffUp) = x
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
      #i <- 1
      # set frame by plotting invisible object
      graphics::plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(0, xMax), ylim = c(0, 300),
           xaxt='n', yaxt='n', ann=FALSE)
      graphics::par(new=F)
      for (j in 1:n.studies) {
        #j <- 1
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
      labelsSeq <- c(unlist(study.numbers), "∑"); labelsSeq
      graphics::axis(side=2, at = atSeq, labels=labelsSeq, las=1)
      # x-axis (effect size)
      atSeq <- seq(0, xMax, by = 10); atSeq
      minDRIFTCoeffLow <- min(DRIFTCoeff[[k]][ ,i]); minDRIFTCoeffLow
      maxDRIFTCoeffUp <- max(DRIFTCoeff[[k]][ ,i]); maxDRIFTCoeffUp
      xRange <- maxDRIFTCoeffUp - minDRIFTCoeffLow; xRange
      labelsSeq <- round(seq(minDRIFTCoeffLow, maxDRIFTCoeffUp, by = (xRange/(length(atSeq)-1))), 2); labelsSeq
      graphics::axis(side=1, at = atSeq, labels=labelsSeq)

      # Add labels and title
      graphics::title(main = paste0("Forest Plot fo the Effects of ", colnames(DRIFTCoeff[[k]])[i]), sub = NULL,
            ylab = "Original Study Number")

      figureFileName <- paste0(activeDirectory, saveFilePrefix," Forest Plot for ", colnames(DRIFTCoeff[[k]])[i], ".png"); figureFileName
      grDevices::dev.copy(png, figureFileName, width = 8, height = 8, units = 'in', res = 300)
      grDevices::dev.off()
    } # END for (i in 1:(n.latent^2))
  } # END if ("forest" %in% plot.type)
  } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  ############################################## discrete time plots  ###################################################
  #######################################################################################################################

  if ("drift" %in% unlist(plot.type)) {

    # check time unit
    if (timeUnit=="timeUnit") {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("The default time lag (timeUnit) has been chosen.", "\n"))
      cat(crayon::blue("Press 'q' to quit and specify or 'c' to continue. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and specify time lag or 'c' to continue without changes. Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    # Function to compute discrete parameters using drift parameters and time-scaling factors
    discreteDrift <-function(driftMatrix, timeScale, number) {
      discreteDriftValue <- OpenMx::expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }


    #######################################################################################################################
    ############## Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  ######################
    #######################################################################################################################

    # can only plot (overlay) mutiple fitted objects if n.latent is identical
    if (length(unique(unlist(n.latent))) > 1) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
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
        noOfSteps <- length(usedTimeRange)
      } else {
        stepWidth <- timeRange[3]
        usedTimeRange <- seq(timeRange[1], timeRange[2], stepWidth)
        # add empirical lags not yet included
        usedTimeRange <- sort(unique(c(usedTimeRange, unlist(allDeltas))))
        noOfSteps <- length(usedTimeRange)
      }

      ## yLimitsForEffects
      # discrete effects across time range
      discreteDriftCoeff <- list()
      for (g in 1:n.fitted.obj) {
        #g <- 1
        discreteDriftCoeff[[g]] <- array(dim=c(n.studies[[g]], noOfSteps, n.latent[[g]]^2))
        for (h in 1:n.studies[g]) {
          #h <- 5
          for (i in usedTimeRange[1]:noOfSteps){
            #i <- usedTimeRange[1]; i
            timeValue <- i * stepWidth; timeValue
            #c(discreteDrift(matrix(unlist(DRIFTCoeff[[g]][[h]]), n.latent, n.latent), timeValue))
            #discreteDriftCoeff[[g]]
            #discreteDriftCoeff[[g]][h, i, 1:(n.latent[[g]]^2)]
            discreteDriftCoeff[[g]][h, i, 1:(n.latent[[g]]^2)] <- c(discreteDrift(matrix(unlist(DRIFTCoeff[[g]][[h]]), n.latent[[g]], n.latent[[g]]), timeValue))
          }
        }
      }

      # min, max, & yLimitsForEffects
      coeffSeq <- seq(1, n.latent[[1]]^2, 1)[!(seq(1, n.latent[[1]]^2, 1) %in% seq(1, n.latent[[1]]^2, (n.latent[[1]]+1)))]; coeffSeq
      maxYLimitsForEffects <- -99999999
      minYLimitsForEffects <- 99999999
      #yLimitsForEffects <- c()
      if (length(yLimitsForEffects) < 1) {
        for (g in 1:n.fitted.obj) {
          yLimitsForEffectsTmp1 <- round(min(discreteDriftCoeff[[g]][, , coeffSeq]) - .10, 1); yLimitsForEffectsTmp1
          yLimitsForEffectsTmp2 <- round(max(discreteDriftCoeff[[g]][, , coeffSeq]) + .10, 1); yLimitsForEffectsTmp2
          if (yLimitsForEffectsTmp1 < minYLimitsForEffects) {
            yLimitsForEffects[1] <- yLimitsForEffectsTmp1
            minYLimitsForEffects <- yLimitsForEffects[1]
            #minYLimitsForEffects
          }
          if (yLimitsForEffectsTmp2 > maxYLimitsForEffects) {
            yLimitsForEffects[2] <- yLimitsForEffectsTmp2
            maxYLimitsForEffects <- yLimitsForEffects[2]
            #maxYLimitsForEffects
          }
        }
      } ## END if (length(yLimitsForEffects) < 1)
      minYLimitsForEffects
      maxYLimitsForEffects
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
        #g <- 2
        plotPairs[[g]] <- array(dim=c(n.studies[[g]], noOfSteps, 2+n.latent[[g]]^2))
        dotPlotPairs[[g]] <- array(dim=c(n.studies[[g]], noOfSteps, 2+n.latent[[g]]^2))
        for (currentTimeScale in 1:noOfSteps){
          #currentTimeScale <- 5
          for (h in 1:n.studies[[g]]) {
            #h <- 1
            timeValue <- usedTimeRange[currentTimeScale+1]; timeValue
            plotPairs[[g]][h,currentTimeScale,1] <- currentTimeScale; plotPairs[[g]][h,currentTimeScale,1]
            plotPairs[[g]][h,currentTimeScale,2] <- timeValue; plotPairs[[g]][h,currentTimeScale,2]
            for (j in 1:(n.latent[[g]]^2)) {
              #j <- 2
              #DriftForPlot[[g]]
              plotPairs[[g]][h,currentTimeScale,(2+j)] <- discreteDrift(matrix(unlist(DriftForPlot[[g]][h]), n.latent, n.latent), timeValue, j)
              #plotPairs[[g]][h,currentTimeScale,(2+j)]
              if (n.studies[[g]] == 1) tmp <- round(meanDelta[[1]],0) else tmp <- ctmaFit[[g]]$studyList[[h]]$delta_t
              tmp
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
          #j <- 1
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            #g <- 1
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
            #figureFileNameAuto[[g]] <- c(); figureFileNameAuto[[g]]
            for (h in 1:n.studies[[g]]) {
              #h <- 1
              graphics::plot(plotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax), ylim = c(plot.yMin, 1),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFit[[g]]$plotStudyNo)) || (ctmaFit[[g]]$plotStudyNo==TRUE) ) {
                # black circle
                graphics::plot(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     #pch=dot.plot.pch, lty=dot.plot.lty, cex=log(targetRows)/2*dot.plot.cex,
                     xlim = c(xMin, xMax), ylim = c(yMin, 1),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (n.studies[[g]] > 1) {
                  currentLabel <- ctmaFit[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  if (is.null(currentLabel)) currentLabel <- ctmaFit[[g]]$ctmaFit$studyList[[h]]$originalStudyNo; currentLabel
                  if (h < 10) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  if (h > 9) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                  #if (h < 10) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=log(targetRows)/5*dot.plot.cex, col="white")
                  #if (h > 9) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=log(targetRows)/5*3/5*dot.plot.cex, col="white")
                } else {
                  #currentLabel="∞"
                  currentLabel="∑"
                  graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  #graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=log(targetRows)/2*1/5*dot.plot.cex, col="white")
                }
                graphics::par(new=T)
              }
            }
          } # END for (g in 1:n.fitted.obj)

          # plot y-axis
          graphics::plot(c(0,0), type="l", col="white", lwd=1.5, xlim = c(xMin, xMax), ylim = c(yMin, 1), xaxt='n', ann=FALSE, las=1)
          # Correct differences in axis length
          atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
          labelsSeq <- seq(0, (max(usedTimeRange)+1), as.integer(targetRows*stepWidth/12)); labelsSeq
          if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]; atSeq
          if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]; labelsSeq
          graphics::axis(side=1, at = atSeq*stepWidth, labels=labelsSeq, las=2)
          # add labels and title
          graphics::title(main = paste0("Auto-regressive Effects of V", counter), sub = NULL,
                xlab=paste0("Time Lag in ", timeUnit), ylab = "Auto-regressive Beta")

          # SAVE
          graphics::par(new=F)
          #figureFileNameAuto[counter] <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); figureFileNameAuto[counter]
          tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); tmp
          grDevices::dev.copy(png, tmp, width = 8, height = 8, units = 'in', res = 300)
          grDevices::dev.off()
        } # END for (j in seq(1, nlatent^2, (nlatent+1)))
      } ## END PLOT (auto effects)


      ## PLOT (cross effects)
      if (plotCrossEffects == TRUE & nlatent > 1) {
        graphics::plot.new()
        #figureFileNameCross <- list()
        counter <- 0
        coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
        for (j in coeffSeq) {
          #j <- 2
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            #g <- 2
            if (is.null(ctmaFit[[g]]$type)) plot..type <- "l" else plot..type <- ctmaFit[[g]]$type; plot..type
            if (is.null(ctmaFit[[g]]$col)) plot.col <- "grey" else plot.col <- ctmaFit[[g]]$col; plot.col
            if (is.null(ctmaFit[[g]]$lwd)) plot.lwd <- 1.5 else plot.lwd <- ctmaFit[[g]]$lwd; plot.lwd
            if (is.null(ctmaFit[[g]]$lty)) plot.lty <- 1 else plot.lty <- ctmaFit[[g]]$lty; plot.lty
            if (is.null(ctmaFit[[g]]$xMin)) plot.xMin <- xMin else plot.xMin <- ctmaFit[[g]]$xMin; plot.xMin
            if (is.null(ctmaFit[[g]]$xMax)) plot.xMax <- xMax else plot.xMax <- ctmaFit[[g]]$xMax; plot.xMax
            if (is.null(ctmaFit[[g]]$yMin)) plot.yMin <- yMin else plot.yMin <- ctmaFit[[g]]$yMin; yMin
            if (is.null(ctmaFit[[g]]$dot.type)) dot.plot.type <- "b" else dot.plot.type <- ctmaFit[[g]]$dot.type
            if (is.null(ctmaFit[[g]]$dot.col)) dot.plot.col <- "black" else dot.plot.col <- ctmaFit[[g]]$dot.col
            if (is.null(ctmaFit[[g]]$dot.lwd)) dot.plot.lwd <- .5 else dot.plot.lwd <- ctmaFit[[g]]$dot.lwd
            if (is.null(ctmaFit[[g]]$dot.lty)) dot.plot.lty <- 3 else dot.plot.lty <- ctmaFit[[g]]$dot.lty
            if (is.null(ctmaFit[[g]]$dot.pch)) dot.plot.pch <- 16 else dot.plot.pch <- ctmaFit[[g]]$dot.pch
            if (is.null(ctmaFit[[g]]$dot.cex)) dot.plot.cex <- 2 else dot.plot.cex <- ctmaFit[[g]]$dot.cex
            for (h in 1:n.studies[[g]]) {
              #h <- 1
              plot(plotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
                   xlim = c(plot.xMin, plot.xMax), ylim = c(plot.yMin, 1),
                   xaxt='n', yaxt='n', ann=FALSE)
              graphics::par(new=T)
              if ( (is.null(ctmaFit[[g]]$plotStudyNo)) || (ctmaFit[[g]]$plotStudyNo==TRUE) ) {
                plot(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
                     pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
                     #pch=dot.plot.pch, lty=dot.plot.lty, cex=log(targetRows)/2*dot.plot.cex,
                     xlim = c(xMin, xMax), ylim = c(yMin, 1),
                     xaxt='n', yaxt='n', ann=FALSE)
                graphics::par(new=T)
                if (n.studies[[g]] > 1) {
                  currentLabel <- ctmaFit[[g]]$studyList[[h]]$originalStudyNo; currentLabel
                  if (h < 10) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  if (h > 9) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=1/5*dot.plot.cex, col="white")
                  #if (h < 10) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=log(targetRows)/5*dot.plot.cex, col="white")
                  #if (h > 9) graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=log(targetRows)/5*3/5*dot.plot.cex, col="white")
                } else {
                  #currentLabel="∞"
                  currentLabel="∑"
                  #tmp1 <- which(round(plotPairs[[g]][h, ,2], 0) == round(meanDelta[[1]],0)); tmp1
                  graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=2/5*dot.plot.cex, col="white")
                  #graphics::text(dotPlotPairs[[g]][h, ,2], plotPairs[[g]][h, ,2+j], labels=currentLabel, cex=log(targetRows)/2*1/5*dot.plot.cex, col="white")
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
          graphics::title(main = paste0("Cross-lagged Effects of ", driftNames[[g]][j]), sub = NULL,
                xlab=paste0("Time Lag in ", timeUnit), ylab = "Cross-lagged Beta")

          graphics::par(new=F)
          #figureFileNameCross[counter] <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); figureFileNameCross[counter]
          tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNames[[g]][j], ".png"); tmp
          grDevices::dev.copy(png, tmp, width = 8, height = 8, units = 'in', res = 300)
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

      #yMin <-  yLimitsForPower[1]; yMin
      #yMin <-  0; yMin
      #yMax <-  yLimitsForPower[2]; yMax
      #yMax <-  1000; yMax
      #yStep <- round(abs(yMin-yMax)/100, 0); yStep
    tmp <- suppressWarnings(as.numeric(rownames(requiredSampleSizes[[g]])))
    usedTimeRange <- tmp[!(is.na(tmp))]; usedTimeRange
      #targetRows <- max(usedTimeRange)/stepWidth; targetRows
    targetRows <- length(tmp[!(is.na(tmp))]); targetRows
      xMax <- targetRows; xMax
      xMin <- usedTimeRange[1]; xMin
      currentLWD <- c(3, 2, 1); currentLWD
      #pow.plot.col <- rep(c("black", "grey"), length(statisticalPower)); pow.plot.col

      coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
      currentDriftNames <- driftNames[[1]][coeffSeq]; currentDriftNames
      for (j in 1:length(currentDriftNames)) {
        #j <- 1
        offset <- (j-1)*length(statisticalPower); offset
        graphics::par(new=F)
        for (h in 1:length(statisticalPower)) {
          #h <- 1
          plot(requiredSampleSizes[[g]][1:targetRows, (h+offset)],
               type="l", col=pow.plot.col[h], lwd=currentLWD[h],
               main=paste0("Required Sample Size For the Effect of ", currentDriftNames[j]),
               xlab=paste0("Time Lag in ", timeUnit),
               ylab="Required Sample Size",
               xlim = c(xMin, xMax),
               ylim = c(pow.plot.yMin, pow.plot.yMax),
               xaxt='n' #, yaxt='n', ann=FALSE
          )
          graphics::par(new=T)
        }

        # Correct differences in axis length
        atSeq <- seq(0, targetRows, by = as.integer(targetRows/12)); atSeq
        mxNumOfLabels <- 25
        stepWidth <- round(mxNumOfLabels/(max(usedTimeRange)+1)); stepWidth
        labelsSeq <- seq(0, (max(usedTimeRange)+1), 1); labelsSeq
        if(length(atSeq) > length(labelsSeq)) atSeq <- atSeq[0: length(labelsSeq)]
        if(length(atSeq) < length(labelsSeq)) labelsSeq <- labelsSeq[0: length(atSeq)]
        graphics::axis(side=1, at = atSeq, labels=labelsSeq)
        graphics::legend('bottomright', legend=statisticalPower, lty=1, col=pow.plot.col, lwd=currentLWD, bty='n', cex=.75)

        #tmp <- paste0(saveFilePrefix," RequiredSampleSizesFor ", currentDriftNames[j], ".png"); tmp
        tmp <- paste0(activeDirectory, saveFilePrefix," RequiredSampleSizesFor ", currentDriftNames[j], ".png"); tmp
        grDevices::dev.copy(png, tmp, width = 8, height = 8, units = 'in', res = 300)
        grDevices::dev.off()
      }

  } ### END Plotting ###

} ### END function definition
