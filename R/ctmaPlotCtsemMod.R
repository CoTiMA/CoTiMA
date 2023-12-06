#' ctmaPlotCtsemMod
#'
#' @description Plots moderator models using \code{\link{ctsem}} fit objects
#'
#' @param ctStanFitObject The fit object with moderator (TIpred) effects to be plotted
#' @param Tipred.pos the Tipred that represent the moderator. Could be more than one in case of categorical moderators (e.g., Tipred.pos = c(3,4))
#' @param scaleTime factor to increase or decrease the time scale (e.g., 1/12 if estimates where based on yearly intervals and figure should show monthly intervals)
#' @param activeDirectory defines the active directory (where to save plots)
#' @param saveFilePrefix Prefix used for saved plots
#' @param fitSummary Mainl ofr debugging purpose. Saves computation time if provided in addition to the fit object
#' @param mod.sd.to.plot The standard deviation vlaues (default -1, 0, +1) for which the drift effects are plotted
#' @param timeUnit Label for the x-axis
#' @param timeRange time range across which drift effects are plotted
#' @param mod.type Could be either "cont" or "cat"
#' @param no.mod.cats Need to be specified if type = "cat". The number of categories should usually be equal the number of dummy variables used to represent the categorical moderator + 1.
#' @param n.x.labels How many values to be used for indicating time points on the x-axis (0 is automatically added and should not be counted)
#' @param plot.xMin default = 0
#' @param plot.xMax default = NULL
#' @param plot.yMin default = -1
#' @param plot.yMax default = 1
#' @param plot..type default = "l", # 2 dots .. are correct
#' @param plot.lty default = 1
#' @param plot.col default = "grey"
#' @param plot.lwd default =  1.5
#' @param dot.plot.type default =  "b" for the dots indicating the moderator values
#' @param dot.plot.col default ="black" for the dots indicating the moderator values
#' @param dot.plot.lwd default =  .5 for the dots indicating the moderator values
#' @param dot.plot.lty default = 3 for the dots indicating the moderator values
#' @param dot.plot.pch default = 16 for the dots indicating the moderator values
#' @param dot.plot.cex default = 3 for the dots indicating the moderator values
#'
#' @importFrom OpenMx expm
#' @importFrom stats quantile
#' @importFrom graphics plot plot.new polygon abline par rect axis title text legend
#' @importFrom grDevices dev.copy png dev.off
#'
#' @export ctmaPlotCtsemMod
#'
#' @examples
#' #Plot a categorical moderator
#' \dontrun{
#' ctmaPlotCtsemMod(ctStanFitObject = ctsemFit,
#'                  activeDirectory=NULL,
#'                  mod.sd.to.plot = NULL,
#'                  timeUnit = "Months",
#'                  timeRange = c(0, 12, .1),
#'                  mod.type = "cat",
#'                  no.mod.cats = NULL
#' }
#'
#' @return writes png figures to disc using the path specified in the activeDirectory arguments.
#'

ctmaPlotCtsemMod <- function(ctStanFitObject = NULL,
                             fitSummary = NULL,
                             activeDirectory = NULL,
                             Tipred.pos=1,
                             saveFilePrefix="Moderator Plot ",
                             scaleTime=1,
                             mod.sd.to.plot = -1:1,
                             timeUnit = "not specified",
                             timeRange = NULL,
                             mod.type = "cont",
                             no.mod.cats = NULL,
                             n.x.labels = NULL,
                             #
                             plot.xMin = 0,
                             plot.xMax = NULL,
                             plot.yMin = -1,
                             plot.yMax = 1,
                             plot..type = "l", # 2 dots .. are correct
                             plot.lty = 1,
                             plot.col = "grey",
                             plot.lwd = 1.5,
                             dot.plot.type = "b",
                             dot.plot.col = "black",
                             dot.plot.lwd = .5,
                             dot.plot.lty = 3,
                             dot.plot.pch = 16,
                             dot.plot.cex = 3
) {
  # arguments to be added
  {
    #if (is.null(mod.sd.to.plot)) mod.sd.to.plot <- -2:2

    if (mod.type == "cont") no.mod.cats <- NULL

    if (is.null(ctStanFitObject)) {
      ErrorMsg <- "\nNo ctsem fit object has been specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (is(ctStanFitObject) != "ctStanFit") {
      ErrorMsg <- "\nNo The fit object provided is not a ctStanFit fit object! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (is.null(activeDirectory)) {
      ErrorMsg <- "\nNo active directory has been specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ((mod.type == "cat") & (is.null(no.mod.cats)))  {
      ErrorMsg <- "\nThe moderator was specified to be categorical (cat), the number of different categories (no.mod.cats) was not. I need that! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if ((mod.type == "cont") & (!(is.null(no.mod.cats)))  ) {
      Msg <- "Note: The moderator was specified to be continous (cont), and you also provided a value for the number of different categories (no.mod.cats). \n The no.mod.cats argument is ignorded!"
      message(Msg)
    }

    if (is.null(timeRange))  {
      tmp <- max(ctStanFitObject$standata$time); tmp
      timeRange <- seq(0, 1.5 * tmp, .1); timeRange
      Msg <- paste0("Note: Time range was not specified. I from 0 to 1.5 times the largest interval (", round(tmp, 2), ") in the data in .1 steps!")
      message(Msg)
    } else { # CHD Aug 2023
      timeRange <- seq(timeRange[1], timeRange[2], timeRange[3])
    }

    # CHD Aug 2023
    #if (is.null(n.x.labels)) n.x.labels <- max(timeRange)/3  # the 0 is added automatically later, so that e.g., 12 will becomes 13, and 4 X labels will be printed
    if (is.null(n.x.labels)) n.x.labels <- length(timeRange)/3  # the 0 is added automatically later, so that e.g., 12 will becomes 13, and 4 X labels will be printed
    #n.x.labels

    # CHD Aug 2023
    #if (n.x.labels < 1)  {
    if (max(timeRange) < 1)  {
      ErrorMsg <- "\nThe largest time interval is < 1. Cannot plot. Please specify larger time range, e.g., timeRange=seq(0, 10, .1). \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    mod.no.to.plot <- Tipred.pos # only a single continuous moderator can be plotted in a single plot, e.g., 1st, 3rd ...
    #                             note that each categorical moderator counts as k, with k = number of moderator categories - 1. (i.e., dummies)
  }


  # some derivations from input
  {
    print(paste0("#################################################################################"))
    print(paste0("####################### Computing model fit summary. ############################"))
    print(paste0("#################################################################################"))
    if (is.null(fitSummary)) ctStanFitObject_summary <- summary(ctStanFitObject) else ctStanFitObject_summary <- fitSummary
    n.latent <- ctStanFitObject$ctstanmodelbase$n.latent; n.latent
    tmp1 <- which(!(is.na(ctStanFitObject$ctstanmodelbase$pars$param)))
    tmpPars <- ctStanFitObject$ctstanmodelbase$pars[tmp1,]; tmpPars
    driftPos <- which(tmpPars$matrix == "DRIFT"); driftPos
    driftNames <- tmpPars$param[driftPos]; driftNames # this is row wise
    if (mod.type == "cat") modPos <- mod.no.to.plot:(mod.no.to.plot - 1 + no.mod.cats - 1)
    if (mod.type == "cont") modPos <- mod.no.to.plot
    TIpred.values <- matrix(ctStanFitObject$data$tipredsdata[, modPos], ncol=length(modPos)); head(TIpred.values)
    m.TIpred<- apply(TIpred.values, 2, mean); m.TIpred
    sd.TIpred <- apply(TIpred.values, 2, sd); sd.TIpred
    mod.values.to.plot <- c()
    if (mod.type == "cat") {
      mod.values.to.plot <- apply(TIpred.values, 2, unique); mod.values.to.plot
      tmp1 <- which(mod.values.to.plot == 0); tmp1
      tmp2 <- which(mod.values.to.plot != 0); tmp2
      if (length(tmp1) > 0 ) mod.values.to.plot <- c(0, tmp2) else mod.values.to.plot <- c(tmp1, tmp2)
    } else{
      for (k in mod.sd.to.plot) mod.values.to.plot <- c(mod.values.to.plot, (m.TIpred + (k * sd.TIpred)))
    }
    n.mod.values.to.plot <- toPlot <- length(mod.values.to.plot); n.mod.values.to.plot

    # warning that selected moderator values may be outliers
    if (mod.type == "cont") {
      tmp1 <- which(TIpred.values < min(mod.values.to.plot)); tmp1
      tmp2 <- length(tmp1)/length(TIpred.values); tmp2
      tmp3 <- min(mod.sd.to.plot); tmp3
      tmp4 <- which(TIpred.values > max(mod.values.to.plot)); tmp4
      tmp5 <- length(tmp4)/length(TIpred.values); tmp5
      tmp6 <- max(mod.sd.to.plot); tmp6

      print(paste0("Note: In your sample ", length(tmp1), " ( = ", round(tmp2*100, 4), "%) people have smaller moderator values than ", tmp3, "SD below the sample mean"))
      print(paste0("Note: In your sample ", length(tmp4), " ( = ", round(tmp5*100, 4), "%) people have larger moderator values than ", tmp6, "SD above the sample mean"))
    }

    if ((mod.type == "cat") & (!(is.null(no.mod.cats))) ) {
      if (no.mod.cats == 1) {
        ErrorMsg <- "\nYou want to plot effects of a categorical moderator. The argument \"no.mod.cats\" should be > 1 (The number of dummy variables to represent the categories + 1)
     \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }
  }

  {
    # poke moderator values that should be plotted into time range (for placing moderator labels in the plot)
    xValueForModValue <- stats::quantile(timeRange, probs = seq(0, 1, 1/(toPlot+1))); xValueForModValue # used for positioning of moderator value in plot
    usedTimeRange <- unique(sort(c(xValueForModValue, timeRange))); usedTimeRange # correcting for added time points
    noOfSteps <- length(usedTimeRange); noOfSteps # can be placed later when generation plotPairs

    DRIFTCoeff <- list()
    #
    # raw main effects (correct order in matrix)
    tmp1 <- ctStanFitObject$stanfit$rawest[driftPos]; tmp1
    rawDrift <- matrix(tmp1, n.latent, byrow=TRUE); rawDrift
    #
    counter <- 0
    for (i in 1:n.mod.values.to.plot ) {
      #i <- 1
      counter <- counter + 1 ; counter
      # horizontal position to plot the moderator labels
      currentXValueForModValue <- xValueForModValue[counter+1]; currentXValueForModValue
      ### compute moderated drift matrices
      # raw moderator effects (could be partial - therefore the complex code)
      if (mod.type == "cont") {
        tmp2 <- ctStanFitObject$stanfit$transformedparsfull$TIPREDEFFECT[,driftPos, modPos]; tmp2 # (rowwise extracted)
        #tmp3 <- rownames(ctStanFitObject_summary$tipreds); tmp3
        # CHD 10. Aug 2023
        tmp3 <- which(colnames(ctStanFitObject$ctstanmodelbase$pars) == "sdscale") + modPos ; tmp3
        modName <- colnames(ctStanFitObject$ctstanmodelbase$pars)[tmp3]; modName
        modName <- gsub("_effect", "", modName); modName
        tmp3 <- rownames(ctStanFitObject_summary$tipreds)[grep(modName, rownames(ctStanFitObject_summary$tipreds))]; tmp3
        #
        if (length(driftNames) != (n.latent^2)) { # if some drift effects were specified as NA; CHD added 28.4.23
          tmp4 <- c()
          for (l in 1:length(driftNames)) {
            tmp5 <- grep(unlist(driftNames[l]), tmp3); tmp5
            if (length(tmp5) == 0) tmp4 <- c(tmp4, NA) else tmp4 <- c(tmp4, tmp5)
          }
          tmp4 <- unique(tmp4); tmp4
          tmp4[!(is.na(tmp4))] <- tmp2[!(is.na(tmp4))]
          tmp4[(is.na(tmp4))] <- 0
        } else {
          tmp4 <- tmp2
        }
        rawMod <- matrix(tmp4, n.latent, byrow=TRUE); rawMod # raw moderator effect to be added to raw main effect (followed by tform) (correct order in matrix)
        tmpNames <- paste0("Raw Drift for Moderator Value = ", mod.sd.to.plot[counter], " SD from mean of moderator"); tmpNames
        DRIFTCoeff[[tmpNames]] <- rawDrift + mod.values.to.plot[counter] * rawMod; DRIFTCoeff[[counter]]
      }
      #
      if (mod.type == "cat") {
        if (counter == 1) {
          tmpNames <- paste0("Raw Drift for Moderator Category No ", counter, ". (= raw Drift)"); tmpNames
          DRIFTCoeff[[counter]] <- matrix(tmp1, n.latent, n.latent); DRIFTCoeff[[counter]] # copy main effects (= comparison group)
        } else {
          tmp2 <- ctStanFitObject$stanfit$transformedparsfull$TIPREDEFFECT[,driftPos, modPos[counter-1]]; tmp2
          rawMod <- matrix(tmp2, n.latent, n.latent, byrow=TRUE); rawMod
          tmpNames <- paste0("Raw Drift for Moderator Category No ", counter, "."); tmpNames
          DRIFTCoeff[[tmpNames]] <- rawDrift + rawMod; DRIFTCoeff[[counter]]
        }
      }
    }
    #DRIFTCoeff
    #DRIFTCoeffBackup <- DRIFTCoeff
    #lapply(DRIFTCoeff, function(x) x * 12)

    ### apply tform to drift elements that should be tformed (extracted into tforms) (tforms are rowwise) (effects in correct order)
    # get tforms for drift (rowwise extraction)
    tmp1a <- ctStanFitObject$ctstanmodelbase$pars[, "transform"]; tmp1a
    tmp1b <- ctStanFitObject$ctstanmodelbase$pars[, "param"]; tmp1b
    tmp1c <- tmp1b %in% driftNames; tmp1c
    transforms <- tmp1a[tmp1c]; transforms
    #
    # compute tformed drift effects
    for (k in 1:(length(DRIFTCoeff))) {
      counter <- 0
      for (l in 1:(n.latent)) {
        for (m in 1:(n.latent)) {
          counter <- counter + 1
          param <- DRIFTCoeff[[k]][l,m]; param
          DRIFTCoeff[[k]][l,m] <- eval(parse(text=transforms[counter])); DRIFTCoeff[[k]][l,m]
          DRIFTCoeff[[k]][l,m] <- DRIFTCoeff[[k]][l,m] * scaleTime
        }
      }
    }
    #DRIFTCoeff
    #DRIFTCoeffBackup

    # compute dt drift coefficients (extracted rowumnwise; in the order of driftnames, not yet in ctmaPlot)
    # Function to compute discrete parameters using drift parameters and time-scaling factors
    discreteDrift <- function(driftMatrix, timeScale, number) {
      discreteDriftValue <- OpenMx::expm(timeScale %x% driftMatrix)
      discreteDriftValue[number] }
    # values where to plot symbol for moderator value/cat
    xValueForModValue2 <- xValueForModValue[-length(xValueForModValue)]; xValueForModValue2
    xValueForModValue2 <- xValueForModValue[-1]; xValueForModValue2
    #
    discreteDriftCoeff <- array(dim=c(n.mod.values.to.plot, length(usedTimeRange), n.latent^2))
    for (h in 1:n.mod.values.to.plot) {
      counter <- 0
      for (i in usedTimeRange) {
        counter <- counter + 1
        discreteDriftCoeff[h, counter, 1:(n.latent^2)] <- c(t(discreteDrift(DRIFTCoeff[[h]], i)))
      }
    }

    ## plotting
    toPlot <- n.mod.values.to.plot; toPlot

    plotPairs <- discreteDriftCoeff
    dotPlotPairs <- discreteDriftCoeff
    # eliminate effect sizes if time point should not be printed (all but one for each moderator value)
    for (h in 1:n.mod.values.to.plot) {
      tmp1 <- which(!(usedTimeRange %in% xValueForModValue2[h])); tmp1
      dotPlotPairs[h, tmp1, 1:n.latent^2] <- NA
    }

    ## plot
    autoEffects <- seq(1, n.latent^2, (n.latent+1)); autoEffects
    crossEffects <- seq(1, n.latent^2, 1)[!(seq(1, n.latent^2, 1) %in% seq(1, n.latent^2, (n.latent+1)))]; crossEffects
    plot.xMax <- length(plotPairs[1, ,1]); plot.xMax

    for (i in 1:n.latent^2) {
      #i <- 2
      graphics::plot.new()
      graphics::par(new=F)
      if (i %in% autoEffects) {
        plot.yMin <- min(plotPairs[, ,autoEffects])[1]; plot.yMin
        plot.yMax <- max(plotPairs[, ,autoEffects])[1]; plot.yMax
      } else {
        plot.yMin <- min(plotPairs[, ,crossEffects])[1]; plot.yMin
        plot.yMax <- max(plotPairs[, ,crossEffects])[1]; plot.yMax
      }
      for (h in 1:toPlot) {
        #h <- 1
        # slopes
        currentPlotPair <- plotPairs[h, ,i]
        currentPlotPair <- cbind(1:length(usedTimeRange), currentPlotPair)
        graphics::par(new=T)
        if (h == 1) {
          plot(currentPlotPair, type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
               xlim = c(0, plot.xMax),
               ylim = c(plot.yMin, plot.yMax),
               xaxt='n',
               ann=FALSE)
        } else  {
          plot(currentPlotPair, type=plot..type, col=plot.col, lwd=plot.lwd, lty=plot.lty,
               xlim = c(0, plot.xMax),
               ylim = c(plot.yMin, plot.yMax),
               xaxt='n',
               yaxt='n',
               ann=FALSE)
        }
        # dots
        currentPlotPair <- dotPlotPairs[h, ,i]
        currentPlotPair <- cbind(1:length(usedTimeRange), currentPlotPair)
        tmp1 <- which(!(is.na(currentPlotPair[,2]))); tmp1 #retain only first and the dot position
        tmp1 <- c(1, tmp1); tmp1
        currentPlotPair <- currentPlotPair[tmp1,]
        graphics::par(new=T)
        plot(currentPlotPair[,], type=dot.plot.type, col=dot.plot.col, lwd=dot.plot.lwd,
             pch=dot.plot.pch, lty=dot.plot.lty, cex=dot.plot.cex,
             xlim = c(0, plot.xMax), ylim = c(plot.yMin, plot.yMax),
             xaxt='n', yaxt='n', ann=FALSE)
        #
        currentLabel <- ""
        if (mod.type == "cont") currentLabel <- mod.sd.to.plot[h]; currentLabel
        if (mod.type == "cat") currentLabel <- h-1; currentLabel
        #graphics::par(new=T)
        if (nchar(currentLabel) == 1) graphics::text(currentPlotPair, labels=currentLabel, cex=1.2, col="white")
        if (nchar(currentLabel) == 2) graphics::text(currentPlotPair, labels=currentLabel, cex=.8, col="white")
        if (nchar(currentLabel) == 3) graphics::text(currentPlotPair, labels=currentLabel, cex=.6, col="white")
        if (nchar(currentLabel) == 4) graphics::text(currentPlotPair, labels=currentLabel, cex=.6, col="white")
        if (nchar(currentLabel) > 4) graphics::text(currentPlotPair, labels=currentLabel, cex=.4, col="white")
        #graphics::par(new=T)
      }
      # axis
      x.labels <- seq(0, max(timeRange),(max(timeRange)/ n.x.labels))[-1]; x.labels # without 0
      stepSize <- (plot.xMax/length(x.labels)); stepSize
      x.pos <- seq(0, plot.xMax, stepSize); x.pos
      x.labels <- round(seq(0, max(timeRange),(max(timeRange)/ n.x.labels)), 2); x.labels # now with 0
      graphics::par(new=T)
      axis(1, labels=x.labels,
           at = x.pos, las=2)

      # titel and axis labels
      if (i %in% autoEffects) tmp1 <- "Autoregressive" else tmp1 <- "Cross-lagged"
      graphics::title(main = paste0("Moderated ", tmp1, " Effects of ", driftNames[i]), sub = NULL,
                      xlab=paste0("Time Interval in ", timeUnit), ylab = paste0(tmp1, " Beta"))

      #tmp <- paste0(activeDirectory, "Moderator Plot "," ", driftNames[i], ".png"); tmp
      tmp <- paste0(activeDirectory, saveFilePrefix," ", driftNames[i], ".png"); tmp

      grDevices::dev.copy(grDevices::png, tmp, width = 8, height = 8, units = 'in', res = 300)
      grDevices::dev.off()
    }
  }
  graphics::par(new=F)
  return(DRIFTCoeff)
}

