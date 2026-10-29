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
#' @param timeRange Numeric start/end/positive-step vector (e.g., c(1, 144, 1)); start must be nonnegative and end at least start. Power plots retain all available time lags within these bounds, including fractional lags.
#' @param yLimitsForEffects range for y-axis
#' @param mod.values moderator values that should be used for plots
#' @param mod.number moderator number that should be used for plots
#' @param aggregateLabel label to indicate aggregated discrete time effects
#' @param xLabels Optional finite numeric tick positions on the time axis; labels show these actual time values.
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
  # Non-interactive devices inhibit recording by default; dev.copy needs it.
  grDevices::dev.control(displaylist = "enable")
  graphics::par(graphics::par(c("mar", "oma", "cex", "ps", "family", "font",
    "font.main", "font.axis", "font.lab", "las", "mgp")))


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

    # A single fit and a (possibly named) list of fits are different shapes.
    single <- inherits(ctmaFitObject, "CoTiMAFit") &&
      any(c("plot.type", "n.latent", "modelResults") %in% names(ctmaFitObject))
    if (single) ctmaFitObject <- list(ctmaFitObject)
    if (!is.list(ctmaFitObject) || !length(ctmaFitObject) ||
        !all(vapply(ctmaFitObject, function(x) inherits(x, "CoTiMAFit"), logical(1)))) stop(
      "Supply a CoTiMAFit object or a nonempty list of CoTiMAFit objects.", call. = FALSE)
    ctmaFitObject <- lapply(ctmaFitObject, .ctmaPlotCompatible)
    if (length(timeRange)) .ctmaPlotTimeRange(timeRange, numeric())
    if (length(yLimitsForEffects) && (length(yLimitsForEffects) != 2L ||
        any(!is.finite(yLimitsForEffects)) || yLimitsForEffects[1] >= yLimitsForEffects[2])) stop(
      "yLimitsForEffects must contain two finite increasing limits.", call. = FALSE)
    if (!length(yLimitsForEffects)) yLimitsForEffects <- NULL

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
        if (is.null(ctmaFitObject[[i]]$plot.type)) {
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

    if (!is.character(activeDirectory) || length(activeDirectory) != 1L ||
        is.na(activeDirectory) || !dir.exists(activeDirectory)) stop(
      "activeDirectory must name an existing directory for saved plots.", call. = FALSE)
    activeDirectory <- file.path(activeDirectory, "")

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
    mod.values <- rep(list(mod.values.backup), length(ctmaFitObject))
    n.primary.studies <- c()
    # CHD added Nov 2022
    allAvgDeltas <- list()


    # detect possible categorical moderator values
    for (i in 1:n.fitted.obj) {
      # CHD changes necessary because return list of ctmaFit changed in Sep 2022
      if (!(is.null(ctmaFitObject[[i]]$argumentList$mod.type))) ctmaFitObject[[i]]$mod.type <- ctmaFitObject[[i]]$argumentList$mod.type
      #
      if (!(is.null(ctmaFitObject[[i]]$mod.type))) {
        if (ctmaFitObject[[i]]$mod.type == "cat") {
          mod.values[[i]] <- c(-999, .ctmaPlotCategories(rownames(ctmaFitObject[[i]]$summary$mod.effects)))
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
        DRIFTCoeff[[i]] <- .ctmaPlotStudyMatrix(DRIFTCoeff[[i]], n.studies[i], n.latent[i], driftNames[[i]])


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

        DRIFTSE[[i]] <- .ctmaPlotStudyMatrix(DRIFTSE[[i]], n.studies[i], n.latent[i], driftNames[[i]])
        FixedEffect_Drift[[i]] <-  ctmaFitObject[[i]]$summary$estimates$`Fixed Effects of Drift Coefficients`["FixedEffect_Drift",]; FixedEffect_Drift[[i]]
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
          allAvgDeltas[[i]] <- .ctmaPlotAverageLags(ctmaFitObject[[i]], n.studies[i])
          if (!undoTimeScaling) {
            scale <- ctmaFitObject[[i]]$summary$scaledTime
            if (is.null(scale)) scale <- ctmaFitObject[[i]]$argumentList$scaleTime
            if (!is.null(scale)) allAvgDeltas[[i]] <- allAvgDeltas[[i]] * scale
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
        graphics::polygon(c(lowLeft, avgX, lowRight), c(yMax2, 0, yMax2), lty=3)
        graphics::abline(v=avgX, col="black", lwd=1.5, lty=1)

        figureContent <- colnames(DRIFTCoeff[[k]])[i]
        if (is.null(figureContent) || is.na(figureContent)) figureContent <- driftNames[[k]][i]
        # Add labels and title
        graphics::title(main = paste0("Funnel Plot for the Effect of ", figureContent))
        figureName[i] <- file.path(activeDirectory, paste0(saveFilePrefix,
          if (n.fitted.obj > 1L) paste0("_fit", k) else "", "_Funnel Plot for ", figureContent, ".png"))
        .ctmaPlotSave(figureName[i])
      } # END for (i in 1:(n.latent^2))
    } # END if ("funnel" %in% plot.type)
  } # END for (k in 1:n.fitted.obj)


  #######################################################################################################################
  #################################################### forest plots #####################################################
  #######################################################################################################################
  for (k in seq_len(n.fitted.obj)) {
    if ("forest" %in% plot.type[[k]]) {
      nl <- n.latent[k]; ns <- n.studies[k]
      estimates <- .ctmaPlotStudyMatrix(DRIFTCoeff[[k]], ns, nl, driftNames[[k]])
      errors <- .ctmaPlotStudyMatrix(DRIFTSE[[k]], ns, nl, driftNames[[k]])
      if (any(!is.finite(estimates)) || any(!is.finite(errors)) || any(errors < 0)) stop(
        "Forest plots require finite effect estimates and non-negative finite SEs.", call. = FALSE)
      lower <- estimates - 1.96 * errors; upper <- estimates + 1.96 * errors
      avg <- as.numeric(FixedEffect_Drift[[k]])
      low <- as.numeric(FixedEffect_DriftLow[[k]]); up <- as.numeric(FixedEffect_DriftUp[[k]])
      if (any(!is.finite(c(avg, low, up))) || any(low > up)) stop(
        "Forest plots require finite ordered aggregate confidence limits.", call. = FALSE)
      auto <- seq.int(1, nl^2, by = nl + 1L)
      groups <- list(auto, setdiff(seq_len(nl^2), auto))
      ticks <- if (!is.null(xLabels)) .ctmaPlotAxis(xLabels, c(0, 1))$at else NULL
      ranges <- lapply(groups, function(cols) {
        if (!length(cols)) return(NULL)
        r <- range(c(lower[, cols], upper[, cols], low[cols], up[cols], ticks))
        if (diff(r) == 0) r <- r + c(-1, 1) * max(1, abs(r[1])) * .1
        r
      })
      precision <- as.numeric(sampleSize[[k]])
      if (length(precision) != ns || any(!is.finite(precision)) || any(precision <= 0)) stop(
        "Forest plots require one positive sample size per study.", call. = FALSE)
      yMax <- 300; xMax <- 300; height <- yMax / (ns + 1L)
      y <- yMax - (seq_len(ns) - .5) * height
      pooled_y <- yMax - (ns + .5) * height
      for (i in seq_len(nl^2)) {
        r <- ranges[[if (i %in% auto) 1L else 2L]]
        scale <- function(value) (value - r[1]) / diff(r) * xMax
        plot(c(0, 0), type = "n",
          xlim = c(if (is.null(ctmaFitObject[[k]]$xMin)) 0 else ctmaFitObject[[k]]$xMin,
            if (is.null(ctmaFitObject[[k]]$xMax)) xMax else ctmaFitObject[[k]]$xMax),
          ylim = c(if (is.null(ctmaFitObject[[k]]$yMin)) 0 else ctmaFitObject[[k]]$yMin,
            if (is.null(ctmaFitObject[[k]]$yMax)) yMax else ctmaFitObject[[k]]$yMax),
          xaxt = "n", yaxt = "n", ann = FALSE)
        size <- precision / (1.1 * max(precision)) * height
        graphics::rect(scale(estimates[, i]) - size/2, y - size/2,
          scale(estimates[, i]) + size/2, y + size/2, col = "grey")
        graphics::segments(scale(lower[, i]), y, scale(upper[, i]), y)
        graphics::polygon(x = scale(c(low[i], avg[i], up[i], avg[i], low[i])),
          y = c(pooled_y, pooled_y + height/3, pooled_y, pooled_y - height/3, pooled_y), col = "black")
        graphics::abline(v = scale(avg[i]), lty = 2)
        graphics::axis(2, at = c(y, pooled_y), labels = c(study.numbers[[k]], "SUM"), las = 1)
        axis <- .ctmaPlotAxis(xLabels, r)
        graphics::axis(1, at = scale(axis$at), labels = axis$labels)
        graphics::title(main = paste0("Forest Plot for the Effects of ", driftNames[[k]][i]),
          ylab = "Study Number", xlab = "Continuous Time Drift Coefficient")
        tag <- if (n.fitted.obj > 1L) paste0("_fit", k) else ""
        .ctmaPlotSave(file.path(activeDirectory, paste0(saveFilePrefix, tag,
          " Forest Plot for ", driftNames[[k]][i], ".png")))
      }
    }
  }


  #######################################################################################################################
  ############################################## discrete time plots  ###################################################
  #######################################################################################################################

  xLabelsBckup <- xLabels
  if ("drift" %in% unlist(plot.type)) {
    usedTimeRange <- .ctmaPlotTimeRange(timeRange, c(allDeltas, allAvgDeltas))
    stepWidth <- if (length(timeRange)) timeRange[3] else 1
    noDecims <- 10L
    moderatorPositions <- vector("list", n.fitted.obj)


    if ( (plotCrossEffects == TRUE) || (plotAutoEffects == TRUE) ) {

      # Function to compute discrete parameters using drift parameters and time-scaling factors
      discreteDrift <-function(driftMatrix, timeScale, number) {
        discreteDriftValue <- OpenMx::expm(timeScale %x% driftMatrix)
        discreteDriftValue[number] }

      #######################################################################################################################
      ############## Extracting Parameters from Fitted Primary Studies created with ctmaInit Function  ######################
      #######################################################################################################################

      # can only plot (overlay) multiple fitted objects if nlatent is identical
      if (length(unique(unlist(n.latent))) > 1) {
        if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
        ErrorMsg <- "A fitted CoTiMA object has to be supplied to plot something. \nA fitted CoTiMA object has to be supplied to plot something. \nGood luck for the next try!"
        stop(ErrorMsg)
      }

      #######################################################################################################################
      ############################ Specification of Parameters for Plotting, Optimal Lags  ##################################
      #######################################################################################################################
      {
        noOfSteps <- length(usedTimeRange)

        # discrete effects across time range
        discreteDriftCoeff <- linearizedTIpredEffect <- DRIFThi <- DRIFTlo <- list()

        #g <- 1
        for (g in 1:n.fitted.obj) {
          toPlot <- n.studies[[g]]; toPlot

          ########################## start dealing with possible moderator values #############################################
          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {

            toPlot <- length(unlist(mod.values[[g]])); toPlot

            # augment usedTimeRange by time points (quantiles) where moderator values are plotted
            xValueForModValue <- stats::quantile(usedTimeRange, probs = seq(0, 1, 1/(toPlot+1))); xValueForModValue # used for positioning of moderator value in plot
            moderatorPositions[[g]] <- xValueForModValue
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
              #tmp1 <- ctmaFitObject[[g]]$studyFitList$stanfit$rawest[1:(nlatent^2)]; tmp1
              tmp1 <- ctmaFitObject[[g]]$studyFitList$stanfit$rawest[tmp1c]; tmp1
              if (length(tmp1) != nlatent^2) stop(
                "This moderated fit lacks a complete retained DRIFT parameter mapping.", call. = FALSE)
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
                #tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT[,1:(nlatent^2),
                #                                                                                 ((n.TIpreds+1):(n.TIpreds+mod.number))]; tmp2
                # CHD 11. Oct 203
                tmp1a <- which(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$matrix == "DRIFT"); tmp1a
                tmp1b <- which(!(is.na(ctmaFitObject[[g]]$studyFitList$ctstanmodelbase$pars$param))); tmp1b
                tmp1c <- which(tmp1b %in% tmp1a); tmp1c
                #tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,1:(nlatent^2),
                #                                                                                 ((n.TIpreds+1):(n.TIpreds+mod.number))]; tmp2
                effects <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT
                if (is.null(effects)) effects <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT
                if (length(dim(effects)) != 3L || n.TIpreds + mod.number > dim(effects)[3]) stop(
                  "The requested moderator has no retained TIPREDEFFECT draws.", call. = FALSE)
                tmp2 <- apply(effects[, tmp1c, n.TIpreds + mod.number, drop = FALSE], 2L, mean)
                tmp3 <- rownames(ctmaFitObject[[g]]$modelResults$MOD); tmp3
                tmp4 <- c()
                for (l in 1:length(driftNames[[g]])) {
                  tmp5 <- grep(unlist(driftNames[[g]][l]), tmp3); tmp5
                  if (length(tmp5) == 0) tmp4 <- c(tmp4, NA) else tmp4 <- c(tmp4, tmp5)
                }
                if (length(tmp2) != length(tmp4)) stop(
                  "Moderator coefficients do not match the DRIFT parameter mapping.", call. = FALSE)
                tmp4[!(is.na(tmp4))] <- tmp2[!(is.na(tmp4))]
                tmp4[(is.na(tmp4))] <- 0

                tmp2 <- matrix(tmp4, n.latent[[g]], byrow=TRUE); tmp2 # raw moderator effect to be added to raw main effect (followed by tform)
                DRIFTCoeff[[g]][[counter]] <- tmp1 + unlist(mod.values[[g]])[counter] * tmp2; DRIFTCoeff[[g]][[counter]]
                names(DRIFTCoeff[[g]])[counter] <- paste0("Moderator Value = ", mod.values[[g]][counter], " SD from mean if standardized (default setting)")
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
                  #  tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,1:(nlatent^2),
                  #                                                                                   n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  #} else {
                  #  tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT[,1:(nlatent^2),
                  #                                                                               n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  #}
                  # CHD changed 10 Nov 2022
                  if (!(is.null(ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull))) {
                    tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedparsfull$TIPREDEFFECT[,(offset+1):(offset+nlatent^2),
                                                                                                     n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  } else {
                    tmp2 <- ctmaFitObject[[g]]$studyFitList$stanfit$transformedpars$TIPREDEFFECT[,(offset+1):(offset+nlatent^2),
                                                                                                 n.TIpreds+(n.mod.tmp*(mod.number-1)+counter)-1]; tmp2
                  }

                  if (!(is.matrix(tmp2))) { # 29. Aug. 2022
                    tmp2 <- matrix(tmp2, nlatent, nlatent, byrow=TRUE); tmp2
                  } else {
                    tmp2 <- matrix(apply(tmp2, 2, mean), nlatent, nlatent, byrow=TRUE); tmp2 # new 18.7. 2022
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
              linearizedTIpredEffect[[g]][[k]] <- matrix(NA, nlatent, nlatent)
              counter <- 0
              for (l in 1:(nlatent)) {
                for (m in 1:(nlatent)) {
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
        if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[g]])
        plotPairs[[g]] <- array(dim=c(toPlot, length(usedTimeRange), 1+n.latent[[g]]^2))
        dotPlotPairs[[g]] <- array(dim=c(toPlot, length(usedTimeRange), 1+n.latent[[g]]^2))

        if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
          xValueForModValue <- moderatorPositions[[g]]
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
              plotPairs[[g]][h,stepCounter,(1+j)] <- discreteDrift(.ctmaPlotDriftAt(DriftForPlot[[g]], h, n.latent[[g]]), timeValue, j)
              #plotPairs[[g]][h,,]
              if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) {
                tmp <- if (toPlot == 1L) meanDelta[g] else allAvgDeltas[[g]][h]
              } else tmp <- moderatorPositions[[g]][h + 1L]
              # CHD added 5 Nov 2022
              tmp <- round(tmp, noDecims); tmp
              # CHD 11. Oct 2023
              #if (isTRUE(all.equal(timeValue, tmp, tolerance = 1e-8))) { # plot only if the (used) time range includes the current study's mean time lag
              if ( (is.null(ctmaFitObject[[g]]$modelResults$MOD)) & (isTRUE(all.equal(timeValue, tmp, tolerance = 1e-8))) )  { # set dots if NO moderator is plotted
                dotPlotPairs[[g]][h, stepCounter, 1] <- timeValue
                dotPlotPairs[[g]][h, stepCounter, (1+j)] <- discreteDrift(.ctmaPlotDriftAt(DriftForPlot[[g]], h, n.latent[[g]]), timeValue, j)
              }
              if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD)))  { # set dots if moderator is plotted
                if (timeValue == xValueForModValue2[h]) {
                  tmp1 <- xValueForModValue[h+1]; tmp1
                  # CHD 11. Oct 2023
                  if (!(is.null(ctmaFitObject[[g]]$summary$scaleTime))) {
                    tmp1 <- tmp1 * ctmaFitObject[[g]]$summary$scaleTime
                  }
                  dotPlotPairs[[g]][h, stepCounter, 1] <- tmp1
                  dotPlotPairs[[g]][h, stepCounter, (1+j)] <- discreteDrift(.ctmaPlotDriftAt(DriftForPlot[[g]], h, n.latent[[g]]), timeValue, j)
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
          tmp2 <- plotPairs[[g]][, , -1, drop=FALSE]; tmp2 # array where in dim 3 there are nlatent dt effects sizes (do not drop if 1st dim=1)
          # y axis, auto
          yMinAutoTmp <- (min(tmp2[ , , autoCols])-.1); yMinAutoTmp
          if (yMinAutoTmp < yMinAuto) yMinAuto <- yMinAutoTmp; yMinAuto
          yMaxAutoTmp <- (max(tmp2[ , , autoCols])); yMaxAutoTmp
          if (yMaxAutoTmp > yMaxAuto) yMaxAuto <- yMaxAutoTmp; yMaxAuto
          # y axis, cross
          if (!length(crossCols)) {
            yMinCross <- -.1; yMaxCross <- .1
          } else if (!(is.null(yLimitsForEffects))) {
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
        nlatent <- n.latent[[1]]; nlatent
        coeffSeq <- seq(1, nlatent^2, (nlatent+1)); coeffSeq

        for (j in coeffSeq) { # diagonal elements only
          #j <- 1
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            #g <- 1
            if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[g]])

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


            originalStudyNo <- if (!is.null(ctmaFitObject[[g]]$modelResults$MOD)) mod.values[[g]] else study.numbers[[g]]
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
          axis <- .ctmaPlotAxis(xLabels, c(xMin, xMax))
          xLabels <- axis$labels; posForXLabel <- axis$at
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
          tmp <- file.path(activeDirectory, paste0(saveFilePrefix," ", driftNames[[1]][j], ".png")); tmp
          .ctmaPlotSave(tmp)
        } # END for (j in coefSeq)
      } ## END PLOT (auto effects)


      ## PLOT (cross effects)
      if (plotCrossEffects == TRUE & nlatent > 1) {
        graphics::plot.new()
        counter <- 0
        nlatent <- n.latent[[1]]; nlatent
        coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
        for (j in coeffSeq) {
          counter <- counter + 1
          for (g in 1:n.fitted.obj) {
            if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[g]])
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

            #if (is.null(ctmaFitObject[[g]]$modelResults$MOD)) toPlot <- n.studies[[g]] else toPlot <- length(mod.values[[g]])
            originalStudyNo <- if (!is.null(ctmaFitObject[[g]]$modelResults$MOD)) mod.values[[g]] else study.numbers[[g]]
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
          axis <- .ctmaPlotAxis(xLabels, c(xMin, xMax))
          xLabels <- axis$labels; posForXLabel <- axis$at
          graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

          # Add labels and title
          if (!(is.null(ctmaFitObject[[g]]$modelResults$MOD))) {
            driftNamesTmp <- c(t(matrix(driftNames[[1]], nlatent))); driftNamesTmp
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
          .ctmaPlotSave(tmp)
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
    statisticalPower <- ctmaFitObject[[g]]$summary$estimates$`Requested Statistical Power`
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

    times <- suppressWarnings(as.numeric(rownames(requiredSampleSizes[[g]])))
    valid <- is.finite(times) & times > 0 & !duplicated(times)
    bounds <- if (length(timeRange)) timeRange[1:2] else range(times[valid])
    valid <- valid & times >= bounds[1] & times <= bounds[2]
    if (!any(valid)) stop("No stored power time points fall inside the requested timeRange.", call. = FALSE)
    currentRequiredSamleSizes <- requiredSampleSizes[[g]][valid, , drop = FALSE]
    order <- order(times[valid]); currentRequiredSamleSizes <- currentRequiredSamleSizes[order, , drop = FALSE]
    usedTimeRange <- sort(times[valid])
    xMin <- min(usedTimeRange); xMax <- max(usedTimeRange)
    if (xMin == xMax) { xMin <- xMin - .5; xMax <- xMax + .5 }

    currentLWD <- if (is.null(ctmaFitObject[[g]]$pow.lwd)) rev(seq_along(statisticalPower)) else rep_len(pow.plot.lwd, length(statisticalPower))

    coeffSeq <- seq(1, nlatent^2, 1)[!(seq(1, nlatent^2, 1) %in% seq(1, nlatent^2, (nlatent+1)))]; coeffSeq
    currentDriftNames <- driftNames[[1]][coeffSeq]; currentDriftNames
    for (j in seq_along(currentDriftNames)) {
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
      axis <- .ctmaPlotAxis(xLabels, c(xMin, xMax))
      xLabels <- axis$labels; posForXLabel <- axis$at
      graphics::axis(side=1, at = posForXLabel, labels=xLabels, las=2)

      # legend
      graphics::legend('bottomright', legend=statisticalPower, lty=1, col=pow.plot.col, lwd=currentLWD, bty='n', cex=.75)

      tmp <- paste0(activeDirectory, saveFilePrefix," RequiredSampleSizesFor ", currentDriftNames[j], ".png"); tmp
      .ctmaPlotSave(tmp)
    }

  } ### END  ("power" %in% unlist(plot.type))

  if (!(is.null(ctmaFitObject[[1]]$modelResults$MOD))) {
    invisible(round(unlist(linearizedTIpredEffect)[-(1:(nlatent^2))], 4)) # just as a check
  }


} ### END function definition


.ctmaPlotTimeRange <- function(timeRange, observed = numeric()) {
  observed <- sort(unique(as.numeric(unlist(observed))))
  observed <- observed[is.finite(observed) & observed > 0]
  if (!length(timeRange)) {
    if (!length(observed)) stop("No finite positive time intervals are available for plotting.", call. = FALSE)
    return(sort(unique(c(seq(1, max(1, 1.5 * round(max(observed) + 1)), by = 1), observed))))
  }
  if (!is.numeric(timeRange) || length(timeRange) != 3L ||
      any(!is.finite(timeRange)) || timeRange[1] < 0 ||
      timeRange[2] < timeRange[1] || timeRange[3] <= 0) stop(
    "timeRange must contain a finite non-negative start, end >= start, and a positive step.", call. = FALSE)
  sort(unique(c(seq(timeRange[1], timeRange[2], by = timeRange[3]),
    observed[observed >= timeRange[1] & observed <= timeRange[2]])))
}

.ctmaPlotAxis <- function(xLabels, limits) {
  if (is.null(xLabels)) {
    ticks <- pretty(limits)
    ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
    if (!length(ticks)) ticks <- unique(limits)
  } else {
    if (!is.numeric(xLabels) || !length(xLabels) || any(!is.finite(xLabels))) stop(
      "xLabels must be a nonempty vector of finite numeric tick positions.", call. = FALSE)
    ticks <- xLabels
  }
  list(at = ticks, labels = ticks)
}

.ctmaPlotCategories <- function(labels) {
  token <- sub("^\\s*([+-]?[0-9]+(?:\\.[0-9]+)?).*", "\\1", labels, perl = TRUE)
  value <- suppressWarnings(as.numeric(token))
  if (any(!is.finite(value))) stop("Cannot read numeric category values from moderator labels.", call. = FALSE)
  unique(value)
}

.ctmaPlotAverageLags <- function(x, n) {
  values <- .ctmaPlotIntervals(x, n)
  vapply(values, function(v) mean(v[is.finite(v) & v > 0]), numeric(1))
}

.ctmaPlotIntervals <- function(x, count) {
  sizes <- as.integer(x$statisticsList$allTpoints) - 1L
  flat <- as.numeric(unlist(x$statisticsList$allDeltas))
  mapped <- NULL
  if (length(sizes) == count && all(is.finite(sizes)) && all(sizes >= 1) && sum(sizes) == length(flat)) {
    ends <- cumsum(sizes); starts <- c(1L, head(ends, -1L) + 1L)
    mapped <- lapply(seq_len(count), function(i) flat[seq.int(starts[i], ends[i])])
  }
  lapply(seq_len(count), function(i) {
    candidates <- list()
    if (length(x$studyList) >= i) candidates <- c(candidates, list(x$studyList[[i]]$delta_t))
    if (length(x$primaryStudyList$deltas) >= i) candidates <- c(candidates, list(x$primaryStudyList$deltas[[i]]))
    if (length(x$emprawList) >= i) {
      raw <- x$emprawList[[i]]
      if (is.data.frame(raw) || is.matrix(raw)) {
        cols <- grep("^dT[0-9]+$", colnames(raw), value = TRUE)
        if (length(cols)) candidates <- c(candidates, list(as.numeric(as.matrix(raw[, cols, drop = FALSE]))))
      }
    }
    if (!is.null(mapped)) candidates <- c(candidates, list(mapped[[i]]))
    if (count == 1L) candidates <- c(candidates, list(flat))
    for (v in candidates) {
      v <- as.numeric(v); v <- v[is.finite(v) & v > 0]
      if (length(v)) return(v)
    }
    stop("Study-specific time intervals cannot be reconstructed unambiguously from this fit.", call. = FALSE)
  })
}

.ctmaPlotCompatible <- function(x) {
  n <- as.integer(x$n.studies)[1]; nl <- as.integer(x$n.latent)[1]
  if (!is.finite(n) || n < 1 || !is.finite(nl) || nl < 1) stop(
    "The fitted object must supply positive n.studies and n.latent.", call. = FALSE)
  if (is.null(x$studyList) || !length(x$studyList)) {
    ids <- x$statisticsList$originalStudyNumbers
    if (!length(ids)) ids <- x$primaryStudyList$studyNumbers
    count <- max(n, length(ids))
    if (!length(ids)) ids <- seq_len(count)
    if (length(ids) != count) stop("Study identifiers do not match the stored study count.", call. = FALSE)
    # Forest/funnel plots need labels but do not require per-study time intervals.
    deltas <- tryCatch(.ctmaPlotIntervals(x, count), error = function(e) {
      if (n == 1L) rep(list(as.numeric(unlist(x$statisticsList$allDeltas))), count)
      else rep(list(NA_real_), count)
    })
    x$studyList <- lapply(seq_len(count), function(i) list(originalStudyNo = ids[i], delta_t = deltas[[i]]))
  }
  # Recover a single population/reference matrix only when every cell is explicit.
  if (is.null(x$modelResults$DRIFT)) {
    table <- x$summary$estimates
    if (!is.data.frame(table) || n != 1L || !all(c("row", "col", "Mean") %in% names(table))) stop(
      "This legacy fit does not contain recoverable DRIFT estimates.", call. = FALSE)
    rows <- grep("^DRIFT(?:[ _]|$)", rownames(table), perl = TRUE)
    cells <- (table$col[rows] - 1L) * nl + table$row[rows]
    if (length(rows) != nl^2 || anyDuplicated(cells) || !setequal(cells, seq_len(nl^2))) stop(
      "Legacy summary does not uniquely identify every DRIFT matrix cell.", call. = FALSE)
    recover <- function(tab) {
      v <- rep(NA_real_, nl^2); v[cells] <- as.numeric(tab$Mean[rows]); v
    }
    x$modelResults$DRIFT <- recover(table)
    original <- x$summary$estimates_original_time_scale
    if (is.data.frame(original) && identical(rownames(original), rownames(table))) {
      x$modelResults$DRIFToriginal_time_scale <- recover(original)
    }
    labels <- rep(NA_character_, nl^2)
    labels[cells] <- sub("^DRIFT[ _]*", "", rownames(table)[rows])
    x$parameterNames$DRIFT <- labels
    message("Using population/reference DRIFT point estimates recovered from the legacy summary; these are not study-specific curves.")
  }
  if (length(x$parameterNames$DRIFT) != nl^2) {
    labels <- colnames(x$modelResults$DRIFT)
    if (length(labels) != nl^2) labels <- names(x$modelResults$DRIFT)
    if (length(labels) != nl^2) labels <- paste0("DRIFT", seq_len(nl^2))
    x$parameterNames$DRIFT <- labels
  }
  x
}

.ctmaPlotStudyMatrix <- function(x, n, nl, labels) {
  if (is.list(x)) x <- do.call(rbind, lapply(x, as.numeric))
  if (length(x) != n * nl^2) stop("The DRIFT/SE dimensions do not match the stored study count.", call. = FALSE)
  if (n == 1L) x <- matrix(as.numeric(x), nrow = 1L) else x <- as.matrix(x)
  if (!identical(dim(x), c(as.integer(n), as.integer(nl^2)))) stop(
    "The DRIFT/SE matrix must have one row per study and one column per effect.", call. = FALSE)
  colnames(x) <- labels
  x
}

.ctmaPlotDriftAt <- function(x, h, nl) {
  value <- if (is.matrix(x)) x[h, , drop = TRUE] else x[[h]]
  if (length(value) != nl^2) stop("A plotted DRIFT matrix has the wrong number of cells.", call. = FALSE)
  matrix(as.numeric(value), nrow = nl, ncol = nl)
}

.ctmaPlotSave <- function(path) {
  grDevices::dev.copy(grDevices::png, path, width = 8, height = 8, units = "in", res = 300)
  grDevices::dev.off()
  invisible(path)
}
