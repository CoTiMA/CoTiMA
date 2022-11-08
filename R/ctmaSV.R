#' ctmaSV
#'
#' @description derives start values by average discrete time SEM effects, converting them to continuous time, and inversely apply transformations used by 'ctsem'
#'
#' @param ctmaInitFit object to which all single 'ctsem' fits of primary studies has been assigned to (i.e., what has been returned by \code{\link{ctmaInit}})
#' @param activeDirectory defines another active directory than the one used in \code{\link{ctmaInit}}
#' @param coresToUse if negative, the value is subtracted from available cores, else value = cores to use
#' @param primaryStudies if ctmaInitFit does not contain the primaryStudies object created with  \code{\link{ctmaPrep}} it could be added
#' @param replaceSV if TRUE replaces startValues in primaryStudies, else it saves them as list element inits
#'
#'
#' @importFrom crayon red blue
#' @importFrom parallel detectCores
#' @importFrom ctsem ctModel ctWideToLong ctDeintervalise
#' @importFrom OpenMx vech2full expm logm
#' @importFrom lavaan sem inspect
#' @importFrom stats uniroot
#'
#' @examples \dontrun{
#' newPrimaryStudyList <- ctmaSV(ctmaInitFit=CoTiMAInitFit_6)
#' }
#'
#' @export ctmaSV
#'
#' @return returns a modified list of primary studies with starting values added or replaced

#'
ctmaSV <- function(
  ctmaInitFit=NULL,
  activeDirectory=NULL,
  primaryStudies=NULL,
  coresToUse=1,
  replaceSV=TRUE)

{  # begin function definition (until end of file)

  { ### CHECKS

    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory

    # check if fit object is specified
    if (is.null(ctmaInitFit)){
      ErrorMsg <- "\nA fitted CoTiMA object (\"ctmaInitFit\") has to be supplied to analyse something. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.null(primaryStudies))) ctmaInitFit$primaryStudyList <- primaryStudies

  } # END Checks


  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  {
    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
    }

    if (coresToUse >= parallel::detectCores()) {
      coresToUse <- parallel::detectCores() - 1
      Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash. \n"
      message(Msg)
    }
  }

  #######################################################################################################################
  ############## Extracting Parameters from Fitted Primary Studies created with ctmaINIT Function  ######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  {

    if (is.null(ctmaInitFit$primaryStudyList)) ctmaInitFit$primaryStudyList <- ctmaInitFit$studyList # possible compatibility with earlier version

    n.latent <- ctmaInitFit$n.latent; n.latent
    n.manifest <- ctmaInitFit$n.manifest; n.manifest
    if (is.null(n.manifest)) n.manifest <- 0
    if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
    n.studies <- ctmaInitFit$n.studies; n.studies
    manifestNames <- ctmaInitFit$studyFitList[[1]]$ctstanmodel$manifestNames; manifestNames
    driftNames <- driftNamesBackup <- ctmaInitFit$parameterNames$DRIFT; driftNames
    # check if drift effects were fixed to 0
    tmp1 <- rownames(ctmaInitFit$studyFitList[[1]]$resultsSummary$popmeans); tmp1
    tmp2 <- which(!(driftNames %in% tmp1)); tmp2
    if (length(tmp2) != 0) driftNames[tmp2] <- "0"

    driftNames <- c(matrix(driftNames, n.latent, byrow=FALSE));     driftNames

    diffusionNames <- diffusionNamesBackup <- ctmaInitFit$parameterNames$DIFFUSION; diffusionNames
    T0varNames <- T0varNamesBackup <- ctmaInitFit$parameterNames$T0VAR; T0varNames

    allSampleSizes <- ctmaInitFit$statisticsList$allSampleSizes; allSampleSizes
  }


  Msg <- "#################################################################################\n############ Use estimates from lavaan model to compute start values ############ \n################################################################################# \n################# Set up required discrete time lavaan models ################### \n#################################################################################"
  message(Msg)
  {
    {    # lavaan model setup
      modelText <- c()
      # check for drift effects fixed to 0 (reverse order here)
      driftNamesTmp <- driftNamesBackup; driftNamesTmp
      tmp1 <- rownames(ctmaInitFit$studyFitList[[1]]$resultsSummary$popmeans); tmp1
      tmp2 <- which(!(driftNamesTmp %in% tmp1)); tmp2
      if (length(tmp2) != 0) driftNamesTmp[tmp2] <- "0"
      driftNamesTmp <- c(matrix(driftNamesTmp, n.latent, byrow=FALSE)); driftNamesTmp

      counter <- 0
      for (i in 1:n.latent) {
        #i <- 1
        counter <- counter + 1
        modelText[counter] <- ""
        for (j in 1:n.latent) {
          tmp1 <- paste0("V", j, "toV", i); tmp1
          if(tmp1 %in% driftNamesTmp) {  # if drift effect is not fixed to 0
            if (j == 1) modelText[counter] <- paste0(modelText[counter], paste0("V", i, "T1 ~ V", j, "T0"))
            if (j != 1) modelText[counter] <- paste0(modelText[counter], paste0(" + V", j, "T0"))
          }
          modelText
        }
      }
      counter <- n.latent; counter
      for (i in 1:n.latent) {
        for (j in i:n.latent) {
          if (i != j) {
            counter <- counter + 1; counter
            modelText[counter] <- paste0("V", i, "T0 ~~ V", j, "T0")
            counter <- counter + 1; counter
            modelText[counter] <- paste0("V", i, "T1 ~~ V", j, "T1")
          }
        }
      }
    } # END # lavaan model setup

    tmp1 <- paste(modelText, sep="", collapse="\n "); tmp1
    model.2w <- paste0("\n ", tmp1, "\n"); model.2w

    model.full <- list()
    for (k in 1:n.studies) {
      #k <- 16
      model.full[[k]] <- model.2w
      if (length(ctmaInitFit$primaryStudyList$deltas[[k]]) > 1 ) {
        for (l in length(ctmaInitFit$primaryStudyList$deltas[[k]]):2) {
          #l <- 2
          tmp2 <- tmp1; tmp2
          seek1 <- "T1"; seek1
          repl1 <- paste0("T", l); repl1
          tmp2 <- gsub(seek1, repl1, tmp2); tmp2
          seek2 <- "T0"; seek2
          repl2 <- paste0("T", l-1); repl2
          tmp2 <- gsub(seek2, repl2, tmp2); tmp2
          model.full[[k]] <- paste0(model.full[[k]], " ", tmp2, "\n")
        }
      }
    }

    #model.full.fit <- list()
    driftSV <- diffSV <- array(NA, dim=c(n.studies, n.latent, n.latent) ); driftSV
    model.full.fit.summary <- list()
    for (k in 1:n.studies) {
      #k <- 16
      dataTmp <- ctmaInitFit$emprawList[[k]]
      colnames(dataTmp) <- gsub("_", "", colnames(dataTmp))
      model.full.fit <- suppressWarnings(lavaan::sem(model.full[[k]], # its okay to "duplicated elements in model syntax have been ignored":" - it was easier to leave them in
                                                     data = dataTmp,
                                                     sample.nobs = ctmaInitFit$primaryStudyList$sampleSizes[[k]]))
      tmp1 <- model.full.fit@ParTable; tmp1
      # model.full.fit.summary <- summary(model.full.fit) # THIS DOES NOT WORK, weird code required, cf ctmaPower
      model.full.fit.summary$PE <- as.data.frame(cbind(tmp1$lhs, tmp1$op, tmp1$rhs, tmp1$est)); model.full.fit.summary$PE
      colnames(model.full.fit.summary$PE) <- c("lhs", "op", "rhs", "est"); model.full.fit.summary$PE

      # get effects
      for (j in 1:n.latent) {
        for (h in 1:n.latent) {
          #j <- 1
          #h <- 1
          DV <- paste0("V", j); DV
          IV <- paste0("V", h); IV
          # drift
          tmp1a <- grep(DV, model.full.fit.summary$PE$lhs); tmp1a
          tmp1b <- grep(IV, model.full.fit.summary$PE$rhs); tmp1b
          tmp1 <- tmp1a[(tmp1a %in% tmp1b)]; tmp1

          tmp3a <- grep("~", model.full.fit.summary$PE$op); tmp3a
          tmp3b <- grep("~~", model.full.fit.summary$PE$op); tmp3b
          tmp3 <- tmp3a[!(tmp3a %in% tmp3b)]; tmp3

          tmp4 <- tmp1[tmp1 %in% tmp3]; tmp4
          #model.full.fit.summary$PE
          driftSV[k, h,j] <- mean(as.numeric(model.full.fit.summary$PE$est[tmp4])); driftSV[k,h,j]
          # diffusion
          tmp2 <- grep("T0", model.full.fit.summary$PE$lhs); tmp2
          tmp2 <- c(tmp2, grep("T0", model.full.fit.summary$PE$rhs)); tmp2
          tmp3 <- tmp3b[!(tmp3b %in% tmp2)]; tmp3
          tmp4 <- tmp1[tmp1 %in% tmp3]; tmp4
          #tmp5 <- mean(model.full.fit.summary$PE$est[tmp4]); tmp5 #diffSV[k, h, j]
          #tmp5[upper.tri(tmp5)] <- tmp5[lower.tri(tmp5)]; tmp5
          #diffSV[k, , ] <- tmp5; diffSV[k, , ]
          diffSV[k, h, j]  <- mean(as.numeric(model.full.fit.summary$PE$est[tmp4])) #diffSV[k, h, j]
        }
      }
      tmp5 <- diffSV[k, , ]
      tmp5[upper.tri(tmp5)] <- tmp5[lower.tri(tmp5)]; tmp5
      diffSV[k, , ] <- tmp5; diffSV[k, , ]
    }
  }

  Msg <- "################################################################################# \n############ Transform lavaan estimates into continuous time effects ############ \n#################################################################################"
  message(Msg)
  {

    # OLD: functions to invert the tform function of ctsem
    #tformAutoInv <- function(x) {for (i in 1:length(x)) x[i] <- -.5 * log(exp(-.5 * x[i] - 5 * 10^(-7)) - 1); return(x)}
    #tformDiffInv <- function(x) {for (i in 1:length(x)) x[i] <- .5 * log(exp(.1 * x[i] - 10^(-6)) - 1); return(x)}
    #tformT0VarInv <- function(x) {for (i in 1:length(x)) x[i] <- .5 * log(exp(.2 * x[i] - 2 * 10^(-7)) - 1); return(x)}

    # NEW: peek tform functions from model setup
    modStr <- ctmaInitFit$studyFitList[[1]]$ctstanmodelbase$pars
    tmp1 <- which(modStr$matrix == "DRIFT")
    tmp2 <- which(modStr$matrix == "DIFFUSION")
    tmp3 <- which(modStr$matrix == "T0VAR")
    tformDrift <- modStr[tmp1, "transform"] # not yet the inverse
    tformDiff <- modStr[tmp2, "transform"]  # not yet the inverse
    tformT0Var <- modStr[tmp3, "transform"] # not yet the inverse
    # see: https://stackoverflow.com/questions/10081479/solving-for-the-inverse-of-a-function-in-r
    inverse <- function (f, lower = -100, upper = 100) {
      function (y) stats::uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
    } # will be used further below (helps finding inverse of tform functions)

    newInit <- list()
    for (k in 1:n.studies) {
      driftTmp  <- OpenMx::logm(driftSV[k,,]) / mean(ctmaInitFit$studyList[[k]]$delta_t); driftTmp
      driftTmpKron <- driftTmp %x% diag(1, n.latent, n.latent) + diag(1, n.latent, n.latent) %x% driftTmp; driftTmpKron
      driftTmpKronSolv <- solve(driftTmpKron); driftTmpKronSolv
      diffTmp <- driftTmpKronSolv %*% (OpenMx::expm(driftTmpKron * 1) - diag(1, dim(driftTmpKron)[1], dim(driftTmpKron)[2])) %*% c(diffSV[k, , ]); diffTmp
      diffTmp <- solve(driftTmpKronSolv %*% ( OpenMx::expm(driftTmpKron * mean(ctmaInitFit$studyList[[k]]$delta_t)) -
                                                diag(1, dim(driftTmpKron)[1], dim(driftTmpKron)[2]) )) %*% c(diffSV[k, , ]); diffTmp
      diffTmp <- matrix(diffTmp, n.latent, n.latent); diffTmp
      if (!(is.na(ctmaInitFit$primaryStudyList$empcovs[[k]][1]))) {
        T0varTmp <- ctmaInitFit$primaryStudyList$empcovs[[k]][1:n.latent, 1:n.latent]; T0varTmp
      } else {
        if ( length(ctmaInitFit$modelResults$T0VAR[[k]]) == n.latent^2) {
          T0varTmp <- matrix(ctmaInitFit$modelResults$T0VAR[[k]], n.latent, n.latent)
        } else {
          T0varTmp <- OpenMx::vech2full(ctmaInitFit$modelResults$T0VAR[[k]])
        }
      }
      ## tform inv
      #
      T0varTmp <- t(chol(T0varTmp)); T0varTmp
      #diag(T0varTmp) <- tformT0VarInv(diag(T0varTmp)); T0varTmp  # OLD
      counter <- 0
      for (i in 1:n.latent) {
        for (j in 1:n.latent) {
          counter <- counter + 1
          tform_inverse <- inverse(function (param) eval(parse(text=tformT0Var[counter])), -1000, 1000)
          if (!(is.na(tformT0Var[counter]))) T0varTmp[i,j] <- tform_inverse(T0varTmp[i,j])$root
        }
      }
      #
      #diag(driftTmp) <- tformAutoInv(diag(driftTmp)); driftTmp # OLD
      counter <- 0
      for (i in 1:n.latent) {
        for (j in 1:n.latent) {
          counter <- counter + 1
          tform_inverse = inverse(function (param) eval(parse(text=tformDrift[counter])), -1000, 1000)
          if (!(is.na(tformDrift[counter]))) driftTmp[i,j] <- tform_inverse(driftTmp[i,j])$root
        }
      }
      #
      diffTmp <- t(chol(diffTmp)); diffTmp
      #diag(diffTmp) <- tformDiffInv(diag(diffTmp)); diffTmp  # OLD
      counter <- 0
      for (i in 1:n.latent) {
        for (j in 1:n.latent) {
          counter <- counter + 1
          tform_inverse = inverse(function (param) eval(parse(text=tformDiff[counter])), -1000, 1000)
          if (!(is.na(tformDiff[counter]))) diffTmp[i,j] <- tform_inverse(diffTmp[i,j])$root
        }
      }
      #
      newInit[[k]] <- c(c(driftTmp), c(diffTmp[lower.tri(diffTmp, diag=T)]), c(T0varTmp[lower.tri(T0varTmp, diag=T)])); newInit
    }
  }

  newPrimaryStudyList <- ctmaInitFit$primaryStudyList
  if (replaceSV == TRUE) newPrimaryStudyList$startValues <- newInit else newPrimaryStudyList$inits <- newInit
  newPrimaryStudyList$emprawList <- ctmaInitFit$emprawList

  results <- newPrimaryStudyList
  class(results) <- "CoTiMAFit"

  invisible(results)

} ### END function definition
