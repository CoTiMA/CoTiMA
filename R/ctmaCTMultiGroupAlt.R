##### Alternative Function to replace ctMultigroupFit (of ctsem) and mxFitFunctionMultigroup (of OpenMx) #####

#' ctmaCTMultiGroupAlt
#'
#' @param dat ?
#' @param groupings ?
#' @param retryattempts ?
#' @param ctmodelobj ?
#' @param startValues ?
#' @param compSVMethod ?
#' @param fixedmodel ?
#' @param moderators ?
#' @param moderatorLabels ?
#' @param modmodel ?
#' @param extraTries ?
#' @param fastIterativeCoTiMA ?
#' @param tryHard ?
#' @param wtgcsv ?
#' @param finetuneGradient ?
#' @param checkHess ?
#' @param method ?
#' @param activateRPB ?
#'
#' @return
#' @export
#'
#' @examples
#'
ctmaCTMultiGroupAlt <- function(dat=NULL,
                               groupings=NULL,
                               retryattempts=20,
                               ctmodelobj=NULL,  # model to be tested
                               startValues=NULL, # start values as provided by CoTiMAgetStartValues.R or a single study fits (then start values are computed below)
                               compSVMethod = "all",
                               fixedmodel=NULL, # if not specified than take the ctmodelobj
                               moderators = NULL, # vector of moderator values of primary studies (N = 1)
                               moderatorLabels = NULL,
                               modmodel = NULL,
                               extraTries=30,
                               fastIterativeCoTiMA=FALSE,
                               tryHard=FALSE,
                               wtgcsv=c("prev","best","initial"), # where to get start values from
                               #initialGradientIterations=mxAutoOptionValue('Gradient iterations'),
                               finetuneGradient=TRUE,
                               checkHess=TRUE,
                               method=2,     # 1 = use DRIFTbig, 2 = use discreteDRIFTbig
                               activateRPB=FALSE #MH "added 20200302"
)
{
  #######################################################################################################################
  ############################################ CoTiMA (as MODERATED ctsem) ##############################################
  #######################################################################################################################

  # debug <- 0
  # if (debug == 1) {
  #   dat=datawide_all
  #   groupings = groups
  #   retryattempts = retryattempts
  #   startValues = homDRIFTallStartValues
  #   ctmodelobj = homDRIFTallCTmodelobj
  #   fixedmodel = homDRIFTallFixedModel
  #   method = method
  #   moderators = NULL # vector of moderator values of primary studies (N = 1)
  #   moderatorLabels = NULL
  #   modmodel = NULL
  #   extraTries=30
  #   fastIterativeCoTiMA=FALSE
  #   tryHard=FALSE
  #   wtgcsv=c("prev","best","initial") # where to get start values from
  #   #initialGradientIterations=mxAutoOptionValue('Gradient iterations'),
  #   finetuneGradient=TRUE
  #   checkHess=TRUE
  #   method=2
  # }
  # debug <- 0

  { # Model Preparation

    # number of latent variables in models
    nlatents <- ctmodelobj$n.latent; nlatents
    #if (is.null(nlatents)) nlatents <- ctmodelobj$nlatent$values

    # make fixedmodel = ctmodelobj if not spefied otherwise (all effects moderated by groups)
    if (is.null(fixedmodel)) fixedmodel <- ctmodelobj

    # check start Values
    if ( is.null(startValues) ) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("At (",Sys.time(),")" ),
                                     paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("Start values have to be provided for ctMultigroupFitAlt !", sep="\n"))
      cat(crayon::red$bold("\n"))
      stop("Good luck for the next try!")
    }
    if ( is.list(startValues) ) {
      startValues <- ctmaCompSV(singleStudyFits = startValues,
                                        ctModelObj = ctmodelobj, fixedModel = fixedmodel,
                                        modModel = modmodel,
                                        compSVMethod=compSVMethod)
    }


    # make modmodel a model that does not allow for moderation
    if (is.null(modmodel)) {
      modmodel <- ctmodelobj
      modmodel$DRIFT <- matrix("groupfixed", nlatents, nlatents)
      for (i in 1:nlatents) {
        for (j in i:nlatents) {
          modmodel$DIFFUSION[j, i] <- "groupfixed"
          modmodel$T0VAR[j, i] <- "groupfixed"
        }
      }
      modmodel$CINT <- matrix("groupfixed", nlatents, 1)
    }

    # check
    if (any(which(modmodel$DRIFT == "moderated") %in% which(fixedmodel$DRIFT == "groupfixed")) == TRUE) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("At (",Sys.time(),")" ),
                                     paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold("At least one drift effect was specified as fixed accross groups and at the same time moderated!", sep="\n"))
      cat(crayon::red$bold("\n"))
      cat(crayon::red$bold("A \"groufixed\" coefficients cannot be \"moderated\" and vice versa!", sep="\n"))
      cat(crayon::red$bold("\n"))
      stop("Good luck for the next try!")
    }

    groupLabels <- names(groupings)[1]; groupLabels
    if (is.null(groupLabels)) groupLabels <- "Group"; groupLabels
    if (is.null(moderators)) moderators <- 1; moderators
    if (is.null(moderatorLabels)) moderatorLabels <- "Moderator_"; moderatorLabels

    # a function
    stringToMxAlgebra <- function(algString, name=NA, dimnames=NA) {
      eval(substitute(mxAlgebra(tExp, name=name, dimnames=dimnames), list(tExp = parse(text=algString)[[1]])))
    }

    # create OpenMx model to be modfied below

    #utils::head(dat)
    #ctmodelobj
    nullModelFit <- ctsem::ctFit(dat=dat, dataform = 'wide', ctmodelobj=ctmodelobj, retryattempts=retryattempts, nofit=TRUE)
    OpenMxModel <- nullModelFit$mxobj

    # add group (as moderator) and moderator to ObenMx Model
    targetNames <- colnames(OpenMxModel$data@observed); targetNames
    # group
    newGroupingID <- match(groupings, sort(unique(groupings))); newGroupingID # consecutive numbering
    OpenMxModel$data@observed <- cbind(OpenMxModel$data@observed, as.numeric(newGroupingID)); utils::head(OpenMxModel$data@observed)
    colnames(OpenMxModel$data@observed) <- c(targetNames, paste0("groupID_T", 0:(length(groupLabels) - 1)))

    # moderator & moderator x group
    newModeratorID <- match(moderators, sort(unique(moderators))); newModeratorID # consecutive numbering
    grpXmodID <- newGroupingID + length(unique(newGroupingID)) * newModeratorID - length(unique(newGroupingID)); grpXmodID
    OpenMxModel$data@observed <- cbind(OpenMxModel$data@observed,
                                       as.numeric(newModeratorID),
                                       as.numeric(grpXmodID))
    colnames(OpenMxModel$data@observed) <- c(targetNames,
                                             paste0("groupID_T", 0:(length(groupLabels) - 1)),
                                             paste0("moderatorID_T", 0:(length(moderatorLabels) - 1)),
                                             paste0("grpXmodID_T", 0:(length(moderatorLabels) - 1)))

    # change model
    newData <- OpenMx::mxData(observed=OpenMxModel$data@observed, type="raw")
    OpenMxModel <- OpenMx::mxModel(OpenMxModel, newData)

    # define required objects
    intervalNames <- colnames(OpenMxModel$data@observed)[grep("dT", colnames(OpenMxModel$data@observed))]; intervalNames
    uniqueIntervals <- sort(unique(c(OpenMxModel$data@observed[, intervalNames]))); uniqueIntervals
    intervalNamesID <- colnames(OpenMxModel$data@observed)[grep("intervalID_", colnames(OpenMxModel$data@observed))]; intervalNamesID
    uniqueIntervalsID <- unique(c(OpenMxModel$data@observed[, intervalNamesID])); uniqueIntervalsID
    uniqueGroups <- unique(groupings); uniqueGroups
    #uniqueModerators <- sort(unique(moderators)); uniqueModerators
    uniqueModerators <- sort(unique(newModeratorID)); uniqueModerators
    maxTpoints <- ctmodelobj$Tpoints; maxTpoints

    # matrix with defintion variables for Groups (Person Index for picking the correct drift matrix later)
    groupID_T <- list()
    groupID_T <- OpenMxModel$matrices$intervalID_T; groupID_T
    groupID_T <- groupID_T[, 1:length(groupLabels)]; groupID_T
    groupID_T$labels <- sub("interval", "group", groupID_T$labels); groupID_T
    groupID_T$labels <- paste0("data.groupID_T", rep(0:(maxTpoints-2))); groupID_T
    groupID_T$values <- 0; groupID_T
    groupID_T$name <- paste0("groupID_T")

    # matrix with defintion variables for Moderators (Person Index for picking the correct drift matrix later)
    moderatorID_T <- list()
    moderatorID_T <- OpenMxModel$matrices$intervalID_T; moderatorID_T
    moderatorID_T <- moderatorID_T[, 1:length(groupLabels)]; moderatorID_T
    moderatorID_T$labels <- sub("interval", "moderator", moderatorID_T$labels); moderatorID_T
    moderatorID_T$labels <- paste0("data.moderatorID_T", rep(0:(maxTpoints-2))); moderatorID_T
    moderatorID_T$values <- 0; moderatorID_T
    moderatorID_T$name <- paste0("moderatorID_T")

    # matrix with defintion variables for Grous x Moderators (Person Index for picking the correct drift matrix later)
    grpXmodID_T <- list()
    grpXmodID_T <- OpenMxModel$matrices$intervalID_T; grpXmodID_T
    grpXmodID_T <- grpXmodID_T[, 1:length(moderatorLabels)]; grpXmodID_T
    grpXmodID_T$labels <- sub("interval", "grpXmod", grpXmodID_T$labels); grpXmodID_T
    grpXmodID_T$labels <- paste0("data.grpXmodID_T", rep(0:(maxTpoints-2))); grpXmodID_T
    grpXmodID_T$values <- 0; grpXmodID_T
    grpXmodID_T$name <- paste0("grpXmodID_T")

    ###### SET OF START VALUES FOR ALL GROUPS (to be changed if matrices are excluded from group-based moderation) #####
    ### Single Studies
    driftValuesForGroup <- diffusionValuesForGroup <- T0varValuesForGroup <- cintValuesForGroup <- list()
    for (i in 1:length(uniqueGroups)) {
      # DRIFT
      currentValues <- c()
      for (j in 1:nlatents) {
        for (k in 1:nlatents) {
          currentValues <- c(currentValues, startValues[paste0("V",j,"toV",k,"_G",i)])
        }
      }
      driftValuesForGroup[[i]] <- matrix(currentValues, nlatents, nlatents)

      # DIFFUSION
      diffusionValuesForGroup[[i]] <- matrix(0, nlatents, nlatents)
      counter <- 0
      for (j in 1:nlatents) {
        for (k in j:nlatents) {
          counter <- counter +1
          diffusionValuesForGroup[[i]][k,j] <- startValues[paste0("diffusion_eta",k,"_eta",j,"_G",i)]
        }
      }

      # T0Var
      T0varValuesForGroup[[i]] <- matrix(0, nlatents, nlatents)
      counter <- 0
      for (j in 1:nlatents) {
        for (k in j:nlatents) {
          counter <- counter +1
          startValues[paste0("T0var_eta",k,"_eta",j,"_G",i)]
          T0varValuesForGroup[[i]][k,j] <- startValues[paste0("T0var_eta",k,"_eta",j,"_G",i)]
        }
      }

      # CINT
      currentValues <- c()
      for (j in 1:nlatents) currentValues <- c(currentValues, startValues[paste0("cint",j,"_G",i)])
      cintValuesForGroup[[i]] <- matrix(currentValues, nlatents, 1); cintValuesForGroup[[i]]
    }


    ### Moderator Categories
    driftValuesForModerator <- diffusionValuesForModerator <- T0varValuesForModerator <- cintValuesForModerator <- list()
    for (i in 1:length(uniqueModerators)) {
      # DRIFT
      currentValues <- c()
      for (j in 1:nlatents) {
        for (k in 1:nlatents) {
          currentValues <- c(currentValues, startValues[paste0("V",j,"toV",k,"_M",i)])
        }
      }
      driftValuesForModerator[[i]] <- matrix(currentValues, nlatents, nlatents)

      # DIFFUSION
      diffusionValuesForModerator[[i]] <- matrix(0, nlatents, nlatents)
      counter <- 0
      for (j in 1:nlatents) {
        for (k in j:nlatents) {
          counter <- counter +1
          diffusionValuesForModerator[[i]][k,j] <- startValues[paste0("diffusion_eta",k,"_eta",j,"_M",i)]
        }
      }

      # T0Var
      T0varValuesForModerator[[i]] <- matrix(0, nlatents, nlatents)
      counter <- 0
      for (j in 1:nlatents) {
        for (k in j:nlatents) {
          counter <- counter +1
          T0varValuesForModerator[[i]][k,j] <- startValues[paste0("T0var_eta",k,"_eta",j,"_M",i)]
        }
      }

      # CINT
      currentValues <- c()
      for (j in 1:nlatents) currentValues <- c(currentValues, startValues[paste0("cint",j,"_M",i)])
      cintValuesForModerator[[i]] <- matrix(currentValues, nlatents, 1); cintValuesForModerator[[i]]
    }

    ### Across all Groups
    driftValuesForAll <- diffusionValuesForAll <- T0varValuesForAll <- cintValuesForALL <- matrix(NA, nlatents, nlatents)
    # DRIFT
    currentValues <- c()
    for (j in 1:nlatents) {
      for (k in 1:nlatents) {
        currentValues <- c(currentValues, startValues[paste0("V",j,"toV",k)])
      }
    }
    driftValuesForAll <- matrix(currentValues, nlatents, nlatents); driftValuesForAll
    tmp1 <- startValues[grep("to", names(startValues))]; tmp1
    tmp2 <- startValues[grep("_G", names(startValues))]; tmp2
    tmp1 <- tmp1[!(tmp1 %in% tmp2)]; tmp1
    driftValuesForAll <- tmp1; driftValuesForAll

    # DIFFUSION
    counter <- 0
    for (j in 1:nlatents) {
      for (k in j:nlatents) {
        counter <- counter +1
        diffusionValuesForAll[k,j] <- startValues[paste0("diffusion_eta",k,"_eta",j)]
      }
    }
    diffusionValuesForAll
    tmp1 <- startValues[grep("diffusion", names(startValues))]; tmp1
    tmp2 <- startValues[grep("_G", names(startValues))]; tmp2
    tmp1 <- tmp1[!(tmp1 %in% tmp2)]; tmp1
    diffusionValuesForAll <- tmp1; diffusionValuesForAll

    # T0Var
    counter <- 0
    for (j in 1:nlatents) {
      for (k in j:nlatents) {
        counter <- counter +1
        T0varValuesForAll[k,j] <- startValues[paste0("T0var_eta",k,"_eta",j)]
      }
    }
    T0varValuesForAll
    tmp1 <- startValues[grep("T0var", names(startValues))]; tmp1
    tmp2 <- startValues[grep("_G", names(startValues))]; tmp2
    tmp1 <- tmp1[!(tmp1 %in% tmp2)]; tmp1
    T0varValuesForAll <- tmp1; T0varValuesForAll

    # CINT
    currentValues <- c()
    for (j in 1:nlatents) currentValues <- c(currentValues, startValues[paste0("cint",j)])
    cintValuesForAll <- matrix(currentValues, nlatents, 1); cintValuesForAll
    tmp1 <- startValues[grep("cint", names(startValues))]; tmp1
    tmp2 <- startValues[grep("_G", names(startValues))]; tmp2
    tmp1 <- tmp1[!(tmp1 %in% tmp2)]; tmp1
    cintValuesForAll <- tmp1; cintValuesForAll

  } # END Model Preparation


  #######################################################################################################################
  ################################################  CHANGE DRIFTs #######################################################
  #######################################################################################################################

  {
    algebrasString <- unlist(lapply(OpenMxModel$algebras, function(extract) extract$name)); algebrasString

    ### create new DRIFTs, invDRIFTs, DRIFTHATCHs
    oldDrift <- OpenMxModel$matrices$DRIFT; oldDrift
    oldDriftHatch <- OpenMxModel$algebras$DRIFTHATCH; oldDriftHatch
    oldInvDrift <- OpenMxModel$algebras$invDRIFT; oldInvDrift
    oldDiscreteDrift_i <-  OpenMxModel$algebras[algebrasString[grep("discreteDRIFT_i", algebrasString)]]; oldDiscreteDrift_i
    oldDiscreteDrift_T <-  OpenMxModel$algebras[algebrasString[grep("discreteDRIFT_T", algebrasString)]]; oldDiscreteDrift_T

    ### DRIFT_G
    newDrift <- list()
    newDrift[[1]] <- oldDrift # make list to handle both single and mutliple drift matrices
    if (!(any(is.na(driftValuesForAll)))) {
      tmp1 <- which(newDrift[[1]]$labels %in% names(driftValuesForAll)); tmp1
      tmp2 <- which(names(driftValuesForAll) %in% newDrift[[1]]$labels); tmp2
      #newDrift[[1]]$values[tmp1] <- driftValuesForAll; newDrift[[1]]$values
      newDrift[[1]]$values[tmp1] <- driftValuesForAll[tmp2]; newDrift[[1]]$values
      if (!(is.null(fixedmodel))) {
        tmp2 <- suppressWarnings(as.numeric(fixedmodel$DRIFT)); tmp2
        tmp3 <- which(!(is.na(tmp2))); tmp3
        tmp3 <- tmp3[!is.na(tmp3)]; tmp3
        newDrift[[1]]$free[tmp3] <- FALSE
      }
    }
    newDrift

    newInvDrift <- oldInvDrift
    newDriftHatch <- oldDriftHatch
    discreteDrift_im <- oldDiscreteDrift_i
    newDiscreteDrift_T <- oldDiscreteDrift_T
    newDriftTemplate <- newDrift[[1]]

    if (any(!(is.na(driftValuesForGroup[[1]])))) {   # if start values for different groups/studies exists (some drift are free across groups/studies)
      newDrift <- list()
      newInvDrift <- list()
      newDriftHatch <- list()
      discreteDrift_im <- list()
      newDiscreteDrift_T <- list()
      for (i in 1:length(uniqueGroups)) {
        # newDrift
        newDrift[[i]] <- newDriftTemplate
        # labels for parameters
        tmp1 <- which(newDrift[[i]]$free); tmp1
        if (length(tmp1) > 0) {
          newLabels <- paste0(oldDrift$labels[tmp1], "_G", i); newLabels
          newDrift[[i]]$labels[tmp1] <- newLabels; newDrift[[i]]$labels
        } else {
          newLabels <- paste0(oldDrift$labels, "_G", i); newLabels
          newDrift[[i]]$labels <- newLabels; newDrift[[i]]$labels
        }
        # change labels back to unmoderated model if effect is not moderated by group:
        newValues <- driftValuesForGroup[[i]]; newValues
        tmp1 <- which(!is.na(newValues)); tmp1
        newDrift[[i]]$values[tmp1] <- newValues[tmp1]; newDrift[[i]]$values
        newDrift[[i]]$labels[fixedmodel$DRIFT == "groupfixed"] <- oldDrift$labels[fixedmodel$DRIFT == "groupfixed"]
        # label for matrix
        newLabels <- paste0("DRIFT", "_", i)
        newDrift[[i]]$name <- newLabels; newDrift[[i]]$name
      }
    } else {  # dublicate Drift - but identical labels for parameters
      for (i in 1:length(uniqueGroups)) {
        newDrift[[i]] <- newDrift[[1]]
        # label for matrix
        newLabels <- paste0("DRIFT", "_", i)
        newDrift[[i]]$name <- newLabels; newDrift[[i]]$name
      }
    }
    newDrift

    if (length(uniqueModerators > 1)) {
      newDrift <- rep(newDrift, length(uniqueModerators)) # multiply list of drift matrices
      counter <- 0
      #for (i in 1:(length(newDrift)/length(uniqueModerators))) { # assign new labels to drift coefficients. (for each drift matrix)
      #  for (j in 1:length(uniqueModerators)) {
      for (j in 1:length(uniqueModerators)) {
       for (i in 1:(length(newDrift)/length(uniqueModerators))) { # assign new labels to drift coefficients. (for each drift matrix)
      counter <- counter +1
          #newDrift[[counter]]$labels[modmodel$DRIFT == "moderated"] <- paste0(oldDrift$labels[modmodel$DRIFT == "moderated"], "_M", uniqueModerators[j])
          newDrift[[counter]]$labels[modmodel$DRIFT == "moderated"] <- paste0(oldDrift$labels[modmodel$DRIFT == "moderated"], "_M", unique(newModeratorID)[j])
          newLabels <- paste0("DRIFT", "_", counter)
          newDrift[[counter]]$name <- newLabels
          newDrift[[counter]]$values[modmodel$DRIFT == "moderated"] <- driftValuesForModerator[[j]][modmodel$DRIFT == "moderated"]
        }
      }
    }
    newDrift

    # get correct start values for drift and set ubound for drift
    if (fastIterativeCoTiMA == TRUE) {
      for (k in 1:length(newDrift)) {
        newDrift[[k]]$ubound <- 1
        newDrift[[k]]$lbound <- -1
        diag(newDrift[[k]]$ubound) <- 0
        newDrift[[k]]$values <- startValues[grep("toV", names(startValues))]
      }
    }
    newDrift

    if (method == 1) {
      ### DRIFTbig (of ct drift matrices)
      bigString <-unlist(lapply(newDrift, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "DRIFTbig"
      driftBig <-  stringToMxAlgebra(x1, name=x2); driftBig

      ### DRIFT (select from ct DRIFTbig)
      x1 <- paste0("DRIFTbig[ ((grpXmodID_T[1, ", 1, "] - 1) * nlatent + 1):(grpXmodID_T[1, ", 1, "] * nlatent), 1:nlatent]"); x1
      x2 <- paste0("DRIFT")
      drift <- stringToMxAlgebra(x1, name=x2); drift

      ### invDRIFT_G
      newInvDrift <- oldInvDrift; newInvDrift

      ### DRIFTHATCH_G
      newDriftHatch <- oldDriftHatch; newDriftHatch

      ### discreteDRIFT_im
      discreteDrift_im <- oldDiscreteDrift_i; discreteDrift_im

      ### discreteDRIFTbig
      bigString <-unlist(lapply(discreteDrift_im, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "discreteDRIFTbig"
      newDiscreteDriftBig <-  stringToMxAlgebra(x1, name=x2); newDiscreteDriftBig

      ### discreteDRIFT_T
      newDiscreteDrift_T <- oldDiscreteDrift_T; newDiscreteDrift_T
    } # END if (method == 1)

    if (method ==2 ) {
      # invDRIFT_M
      newInvDrift <- list()
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("solve(DRIFT", "_", i, ")"); x1
        x2 <- paste0("invDRIFT_", i); x2
        newInvDrift[[i]] <-  stringToMxAlgebra(x1, name=x2); newInvDrift[[i]]
      }
      newInvDrift

      # DRIFTHACT_M
      newDriftHatch <- list()
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("DRIFT_", i, " %x% II + II %x% DRIFT_", i); x1
        newDriftHatch[[i]] <- stringToMxAlgebra(x1, name=paste0("DRIFTHATCH_", i)); newDriftHatch[[i]]
      }
      newDriftHatch

      # New algebras for discreteDRIFT_im
      discreteDrift_im <- list()
      counter <- 0
      for (i in 1:length(uniqueIntervals)) {
        for (j in 1:(length(uniqueModerators) * length(uniqueGroups))) {
          counter <- counter +1
          x1 <- paste0("expm(DRIFT", "_", j, " %x% ", uniqueIntervals[i], ")"); x1
          x2 <- paste0("discreteDRIFT_im", (i-1)*(length(uniqueModerators) * length(uniqueGroups)) +j ); x2
          discreteDrift_im[[counter]] <-  stringToMxAlgebra(x1, name=x2); discreteDrift_im[[counter]]
        }
      }
      discreteDrift_im

      # New algebras for discreteDRIFTbig
      bigString <-unlist(lapply(discreteDrift_im, function(extract) extract$name)); bigString
      tmp1  <- paste(as.character(bigString), sep="' '", collapse=", "); tmp1
      x1 <- paste0("rbind(", tmp1, ")")
      x2 <- "discreteDRIFTbig"
      newDiscreteDriftBig <-  stringToMxAlgebra(x1, name=x2); newDiscreteDriftBig

      # New algebra for discreteDRIFT_T1, 2, etc
      newDiscreteDrift_T <- list()
      #tmp <- length(unique(grpXmodID)); tmp
      tmp <- length(uniqueModerators) * length(uniqueGroups); tmp
      for (i in 1:(maxTpoints-1)) {
        x1 <- paste0("discreteDRIFTbig[
                     ( ( (", tmp, " * intervalID_T[1, ", i, "] ) - ( ", tmp, " - (grpXmodID_T[1, 1] ) ) )  * nlatent - nlatent + 1):
                     (( (", tmp, " * intervalID_T[1, ", i, "] ) - ( ", tmp, " - (grpXmodID_T[1, 1] ) ) )  * nlatent - nlatent + nlatent),
                     1:nlatent]")
        x2 <- paste0("discreteDRIFT_T", i)
        newDiscreteDrift_T[[i]] <- stringToMxAlgebra(x1, name=x2);
      }
      newDiscreteDrift_T
    } # END if (method == 2)

  }

  #######################################################################################################################
  ################################################  CHANGE DIFFUSIONs ###################################################
  #######################################################################################################################
  {
    oldDiffusionBase <- OpenMxModel$matrices$DIFFUSIONbase; oldDiffusionBase
    oldDiffusionChol <- OpenMxModel$algebras$DIFFUSIONchol; oldDiffusionChol
    oldDiffusion <- OpenMxModel$algebras$DIFFUSION ; oldDiffusion
    oldAsymDiffusion <- OpenMxModel$matrices$asymDIFFUSION; oldAsymDiffusion
    oldAsymDiffusionAlg <- OpenMxModel$algebras$asymDIFFUSIONalg; oldAsymDiffusionAlg
    oldDiscreteDiffusion_i <- OpenMxModel$algebras[algebrasString[grep("discreteDIFFUSION_i", algebrasString)]]; oldDiscreteDiffusion_i
    oldDiscreteDiffusion_T <- OpenMxModel$algebras[algebrasString[grep("discreteDIFFUSION_T", algebrasString)]]; oldDiscreteDiffusion_T

    newDiffusionBase <- list() # make list to handle both single and multiple diffusion matrices
    newDiffusionBase[[1]] <- oldDiffusionBase # make list to handle both single and mutliple drift matrices
    if (!(any(is.na(diffusionValuesForAll)))) {
      tmp1 <- which(newDiffusionBase[[1]]$labels %in% names(diffusionValuesForAll)); tmp1
      tmp2 <- which(names(diffusionValuesForAll) %in% newDiffusionBase[[1]]$labels); tmp2
      #newDiffusionBase[[1]]$values[tmp1] <- diffusionValuesForAll; newDiffusionBase[[1]]$values
      newDiffusionBase[[1]]$values[tmp1] <- diffusionValuesForAll[tmp2]; newDiffusionBase[[1]]$values
      if (!(is.null(fixedmodel))) {
        tmp2 <- suppressWarnings(as.numeric(fixedmodel$DIFFUSION)); tmp2
        tmp3 <- which(!(is.na(tmp2))); tmp3
        tmp3 <- tmp3[!is.na(tmp3)]; tmp3
        newDiffusionBase[[1]]$free[tmp3] <- FALSE
      }
    }
    newDiffusionBase
    newDiffusionBaseTemplate <- newDiffusionBase; newDiffusionBaseTemplate

    newDiffusionChol <- oldDiffusionChol
    newDiffusion <- oldDiffusion
    newAsymDiffusion <- oldAsymDiffusion
    newAsymDiffusionAlg <- oldAsymDiffusionAlg
    discreteDiffusion_im <- oldDiscreteDiffusion_i
    newDiscreteDiffusion_T <- oldDiscreteDiffusion_T

    if (any(!(is.na(diffusionValuesForGroup[[1]])))) {
      newDiffusionChol <- list()
      newDiffusion <- list()
      newAsymDiffusion <- list()
      newAsymDiffusionAlg <- list()
      discreteDiffusion_im <- list()
      newDiscreteDiffusion_T <- list()
      for (i in 1:length(uniqueGroups)) {
        ### DIFFUSIONbase_G
        newDiffusionBase[[i]] <- newDiffusionBaseTemplate[[1]]; newDiffusionBase[[i]]
        # labels for parameters
        tmp <- newDiffusionBase[[i]]
        tmp1 <- which(newDiffusionBase[[i]]$free); tmp1
        if (length(tmp1) > 0) {
          newLabels <- paste0(oldDiffusionBase$labels[tmp1], "_G", i); newLabels
          newDiffusionBase[[i]]$labels[tmp1] <- newLabels; newDiffusionBase[[i]]$labels
        } else {
          newLabels <- paste0(oldDiffusionBase$labels, "_G", i); newLabels
          newDiffusionBase[[i]]$labels <- newLabels; newDiffusionBase[[i]]$labels
        }
        # change labels back to unmoderated model if effect is not moderated by groupd:
        newValues <- diffusionValuesForGroup[[i]]; newValues
        tmp1 <- which(!is.na(newValues)); tmp1
        newDiffusionBase[[i]]$values[tmp1] <- newValues[tmp1]; newDiffusionBase[[i]]$values
        newDiffusionBase[[i]]$labels[fixedmodel$DIFFUSION == "groupfixed"] <- oldDrift$labels[fixedmodel$DIFFUSION == "groupfixed"]
        # label for matrix
        newLabels <- paste0("DIFFUSIONbase", "_", i)
        newDiffusionBase[[i]]$name <- newLabels; newDiffusionBase[[i]]$name
      }
    } else {
      for (i in 1:length(uniqueGroups)) {
        newDiffusionBase[[i]] <- newDiffusionBase[[1]]
        # label for matrix
        newLabels <- paste0("DIFFUSIONbase", "_", i)
        newDiffusionBase[[i]]$name <- newLabels; newDiffusionBase[[i]]$name
      }
    }
    newDiffusionBase

    if (length(uniqueModerators > 1)) {
      newDiffusionBase <- rep(newDiffusionBase, length(uniqueModerators)) # multiply list of drift matrices
      newDiffusionBase
      counter <- 0
      #for (i in 1:(length(newDiffusionBase)/length(uniqueModerators))) { # assign new labels to drift coefficients.
      #for (j in 1:length(uniqueModerators)) {
        for (j in 1:length(uniqueModerators)) {
          for (i in 1:(length(newDiffusionBase)/length(uniqueModerators))) { # assign new labels to drift coefficients.
            counter <- counter +1
          #newDiffusionBase[[counter]]$labels[modmodel$DIFFUSION == "moderated"] <- paste0(oldDiffusionBase$labels[modmodel$DIFFUSION == "moderated"], "_", uniqueModerators[j])
          newDiffusionBase[[counter]]$labels[modmodel$DIFFUSION == "moderated"] <- paste0(oldDiffusionBase$labels[modmodel$DIFFUSION == "moderated"], "_", unique(newModeratorID)[j])
          newLabels <- newDiffusionBase[[counter]]$labels
          newLabels[grep("NA", newLabels)] <- NA
          newDiffusionBase[[counter]]$labels <- newLabels
          newLabels <- paste0("DIFFUSIONbase", "_", counter)
          newDiffusionBase[[counter]]$name <- newLabels
          newDiffusionBase[[counter]]$values[modmodel$DIFFUSION == "moderated"] <- diffusionValuesForModerator[[j]][modmodel$DIFFUSION == "moderated"]
        }
      }
    }
    newDiffusionBase

    # get correct start values for diffusion base
    if (fastIterativeCoTiMA == TRUE) {
      for (k in 1:length(newDiffusionBase)) {
        counter <- 1
        for (l in 1:(nlatents)) {
          for (m in l:(nlatents)) {
            tmp1 <- startValues[grep(paste0("diffusion"), names(startValues))]; tmp1
            tmp2 <- tmp1[grep(paste0("_G", k), names(tmp1))]; tmp2
            newDiffusionBase[[k]]$values[m,l] <- tmp2[counter]
            counter <- counter + 1
          }
        }
      }
    }
    newDiffusionBase

    if (method == 1) {
      ### DIFFUSIONbaseBig
      bigString <-unlist(lapply(newDiffusionBase, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "DIFFUSIONbaseBig"
      diffusionBaseBig <-  stringToMxAlgebra(x1, name=x2); diffusionBaseBig

      ### DIFFUSIONbase
      x1 <- paste0("DIFFUSIONbaseBig[ ((grpXmodID_T[1, ", 1, "] - 1) * nlatent + 1):(grpXmodID_T[1, ", 1, "] * nlatent), 1:nlatent]"); x1
      x2 <- paste0("DIFFUSIONbase")
      diffusionBase <- stringToMxAlgebra(x1, name=x2); diffusionBase

      ### DIFFUSIONchol
      newDiffusionChol <- oldDiffusionChol; newDiffusionChol

      ### DIFFUSION
      newDiffusion <- oldDiffusion; newDiffusion

      ### asymDIFFUSION
      newAsymDiffusion <- oldAsymDiffusion; newAsymDiffusion

      ### asymDIFFUSIONalg_G
      newAsymDiffusionAlg <- oldAsymDiffusionAlg; newAsymDiffusionAlg

      ### discreteDIFFUSION_im
      discreteDiffusion_im <- oldDiscreteDiffusion_i; discreteDiffusion_im

      ### discreteDIFFUSIONbig
      bigString <- unlist(lapply(discreteDiffusion_im, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "discreteDIFFUSIONbig"
      newDiscreteDiffusionBig <-  stringToMxAlgebra(x1, name=x2); newDiscreteDiffusionBig

      ### discreteDIFFUSION_T
      newDiscreteDiffusion_T <- oldDiscreteDiffusion_T; newDiscreteDiffusion_T
    } #END if (method == 1)

    if (method == 2) {
      # DIFFUSIONchol
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("(vec2diag(exp(diag2vec(DIFFUSIONbase", "_", i, "))) + DIFFUSIONbase", "_", i, " - vec2diag(diag2vec(DIFFUSIONbase", "_", i, ")))"); x1 # test
        x2 <- paste0("DIFFUSIONchol", "_", i); x2
        newDiffusionChol[[i]] <- stringToMxAlgebra(x1, name=x2)
      }
      newDiffusionChol

      # DIFFUSION
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("(DIFFUSIONchol", "_", i, " %*% t(DIFFUSIONchol", "_", i, "))"); x1
        x2 <- paste0("DIFFUSION", "_", i); x2
        newDiffusion[[i]] <- stringToMxAlgebra(x1, name=x2)
      }
      newDiffusion

      # asymDIFFUSION
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        newAsymDiffusion[[i]] <- oldAsymDiffusion
        newLabels <- paste0(oldAsymDiffusion$labels)
        newLabels <- sub("alg", paste0("alg_", i),  newLabels); newLabels
        newAsymDiffusion[[i]]$labels <- newLabels; newAsymDiffusion[[i]]$labels
        newLabels <- paste0("asymDIFFUSION", "_", i)
        newAsymDiffusion[[i]]$name <- newLabels; newAsymDiffusion[[i]]$name
      }
      newAsymDiffusion

      # asymDIFFUSIONalg_M
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("((-solve(DRIFTHATCH_", i, ")) %*% cvectorize(DIFFUSION_", i, "))"); x1
        x2 <- paste0("asymDIFFUSIONalg_", i); x2
        newAsymDiffusionAlg[[i]] <- stringToMxAlgebra(x1, name=x2)
      }
      newAsymDiffusionAlg

      # New algebras for discreteDIFFUSION_im
      counter <- 0
      for (i in 1:length(uniqueIntervals)) {
        for (j in 1:(length(uniqueModerators) * length(uniqueGroups))) {
          counter <- counter +1
          x1 <- paste0("(asymDIFFUSION_", j, " - (discreteDRIFT_im", counter, " %&% asymDIFFUSION_", j,"))"); x1
          x2 <- paste0("discreteDIFFUSION_im", (i-1)*(length(uniqueModerators) * length(uniqueGroups)) +j ); x2
          discreteDiffusion_im[[counter]] <-  stringToMxAlgebra(x1, name=x2)
        }
      }
      discreteDiffusion_im

      # New algebras for discreteDIFFUSIONbig
      bigString <- unlist(lapply(discreteDiffusion_im, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "discreteDIFFUSIONbig"
      newDiscreteDiffusionBig <-  stringToMxAlgebra(x1, name=x2)
      newDiscreteDiffusionBig

      # New algebra for discreteDIFFUSION_T1, 2, etc
      #tmp <- length(unique(grpXmodID)); tmp
      tmp <- length(uniqueModerators) * length(uniqueGroups); tmp
      for (i in 1:(maxTpoints-1)) {
        x1 <- paste0("discreteDIFFUSIONbig[
                     ( ( (", tmp, " * intervalID_T[1, ", i, "] ) - ( ", tmp, " - (grpXmodID_T[1, 1] ) ) )  * nlatent - nlatent + 1):
                     (( (", tmp, " * intervalID_T[1, ", i, "] ) - ( ", tmp, " - (grpXmodID_T[1, 1] ) ) )  * nlatent - nlatent + nlatent),
                     1:nlatent]")
        x2 <- paste0("discreteDIFFUSION_T", i)
        newDiscreteDiffusion_T[[i]] <- stringToMxAlgebra(x1, name=x2);
      }
      newDiscreteDiffusion_T
    } # END if (method == 2)
  }


  ######################################################################################################################
  ################################################  CHANGE CINTs #######################################################
  ######################################################################################################################
  {
    ### create new CINT, discreteCINT_T1
    oldCint <- OpenMxModel$matrices$CINT; oldCint
    oldDiscreteCint_i <- OpenMxModel$algebras[algebrasString[grep("discreteCINT_i", algebrasString)]]; oldDiscreteCint_i
    oldDiscreteCint_T <- OpenMxModel$algebras[algebrasString[grep("discreteCINT_T", algebrasString)]]; oldDiscreteCint_T

    newCint <- list()  # make list to handle both single and multiple cint matrices
    newCint[[1]] <- oldCint
    if (!(any(is.na(cintValuesForAll))))       {
      tmp1 <- which(newCint[[1]]$labels %in% names(cintValuesForAll)); tmp1
      tmp2 <- which(names(cintValuesForAll) %in% newCint[[1]]$labels); tmp2
      #newCint[[1]]$values[tmp1] <- cintValuesForAll; newCint[[1]]$values
      newCint[[1]]$values[tmp1] <- cintValuesForAll[tmp2]; newCint[[1]]$values
      if (!(is.null(fixedmodel))) {
        tmp2 <- suppressWarnings(as.numeric(fixedmodel$CINT)); tmp2
        tmp3 <- which(!(is.na(tmp2))); tmp3
        tmp3 <- tmp3[!is.na(tmp3)]; tmp3
        newCint[[1]]$free[tmp3] <- FALSE
      }
    }
    newCint
    newCintTemplate <- newCint; newCintTemplate

    discreteCint_im <- oldDiscreteCint_i
    newDiscreteCint_T <- oldDiscreteCint_T

    ### CINT_G
    if (any(!(is.na(cintValuesForGroup[[1]])))) {
      discreteCint_im <- list()
      newDiscreteCint_T <- list()
      for (i in uniqueGroups) {
        newCint[[i]] <- newCintTemplate[[1]]
        # labels for parameters
        tmp1 <- which(newCint[[i]]$free); tmp1
        if (length(tmp1) > 0) {
          newLabels <- paste0(oldCint$labels[tmp1], "_G", i); newLabels
          newCint[[i]]$labels[tmp1] <- newLabels; newCint[[i]]$labels
        } else {
          newLabels <- paste0(oldCint$labels, "_G", i); newLabels
          newCint[[i]]$labels <- newLabels; newCint[[i]]$labels
        }
        # Change labels back to unmoderated model if effect is not moderated: not yet set up
        newValues <- cintValuesForGroup[[i]]; newValues
        tmp1 <- which(!is.na(newValues)); tmp1
        newCint[[i]]$values[tmp1] <- newValues[tmp1]; newCint[[i]]$values
        newCint[[i]]$labels[fixedmodel$CINT == "groupfixed"] <- oldDrift$labels[fixedmodel$CINT == "groupfixed"]
        # label for matrix
        newLabels <- paste0("CINT", "_", i)
        newCint[[i]]$name <- newLabels; newCint[[i]]$name
      }
    } else {
      for (i in 1:length(uniqueGroups)) {
        newCint[[i]] <- newCint[[1]]
        # label for matrix
        newLabels <- paste0("CINT", "_", i)
        newCint[[i]]$name <- newLabels; newCint[[i]]$name
      }
    }
    newCint

    if (length(uniqueModerators > 1)) {
      newCint <- rep(newCint, length(uniqueModerators)) # multiply list of drift matrices
      counter <- 0
      #for (i in 1:(length(newCint)/length(uniqueModerators))) { # assign new labels to drift coefficients.
      #for (j in 1:length(uniqueModerators)) {
        for (j in 1:length(uniqueModerators)) {
          for (i in 1:(length(newCint)/length(uniqueModerators))) { # assign new labels to drift coefficients.
            counter <- counter +1
          #newCint[[counter]]$labels[modmodel$CINT == "moderated"] <- paste0(oldCint$labels[modmodel$CINT == "moderated"], "_", uniqueModerators[j])
          newCint[[counter]]$labels[modmodel$CINT == "moderated"] <- paste0(oldCint$labels[modmodel$CINT == "moderated"], "_", unique(newModeratorID)[j])
          newLabels <- paste0("CINT", "_", counter)
          newCint[[counter]]$name <- newLabels
          newCint[[counter]]$values[modmodel$CINT == "moderated"] <- cintValuesForModerator[[j]][modmodel$CINT == "moderated"]
        }
      }
    }
    newCint

    # get correct start values for cint
    if (fastIterativeCoTiMA == TRUE) {
      for (k in 1:length(newCint)) {
        tmp1 <- startValues[grep(paste0("_cint"), names(startValues))]; tmp1
        tmp2 <- tmp1[grep(paste0("_G", k), names(tmp1))]; tmp2
        if (length(tmp2) > 0) newCint[[k]]$values <- tmp2
      }
    }
    newCint

    if (method == 1) {
      ### CINTbig
      bigString <-unlist(lapply(newCint, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "CINTbig"
      cintBig <-  stringToMxAlgebra(x1, name=x2); cintBig

      ### CINT
      x1 <- paste0("CINTbig[ ((grpXmodID_T[1, ", 1, "] - 1) * nlatent + 1):(grpXmodID_T[1, ", 1, "] * nlatent), 1]"); x1
      x2 <- paste0("CINT")
      cint <- stringToMxAlgebra(x1, name=x2); cint

      ### discreteCINT_IM
      discreteCint_im <- oldDiscreteCint_i; discreteCint_im

      ### discreteCINTbig
      bigString <-unlist(lapply(discreteCint_im, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "discreteCINTbig"
      newDiscreteCintBig <-stringToMxAlgebra(x1, name=x2); newDiscreteCintBig

      ### discreteCINT_T
      newDiscreteCint_T <- oldDiscreteCint_T; newDiscreteCint_T
    }

    if (method == 2) {
      ### discreteCINT_IM
      discreteCint_im <- list()
      counter <- 0
      for (i in 1:length(uniqueIntervals)) {
        for (j in 1:(length(uniqueModerators) * length(uniqueGroups))) {
          counter <- counter +1
          x1 <- paste0("solve(DRIFT", "_", j, ") %*% (expm(DRIFT", "_", j, " %x% ", uniqueIntervals[i], ") - II) %*% CINT_", j); x1
          x2 <- paste0("discreteCINT_im", (i-1)*(length(uniqueModerators) * length(uniqueGroups)) +j ); x2
          discreteCint_im[[counter]] <-  stringToMxAlgebra(x1, name=x2); discreteCint_im[[counter]]
        }
      }
      discreteCint_im

      bigString <-unlist(lapply(discreteCint_im, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "discreteCINTbig"
      newDiscreteCintBig <-stringToMxAlgebra(x1, name=x2); newDiscreteCintBig

      newDiscreteCint_T <- list()
      tmp <- length(unique(grpXmodID)); tmp
      tmp <- length(uniqueGroups) * length(uniqueModerators); tmp
      for (i in 1:(maxTpoints-1)) {
        x1 <- paste0("discreteCINTbig[
                     ( ( (", tmp, " * intervalID_T[1, ", i, "] ) - ( ", tmp, " - (grpXmodID_T[1, 1] ) ) )  * nlatent - nlatent + 1):
                     (( (", tmp, " * intervalID_T[1, ", i, "] ) - ( ", tmp, " - (grpXmodID_T[1, 1] ) ) )  * nlatent - nlatent + nlatent),
                     1]")
        x2 <- paste0("discreteCINT_T", i)
        newDiscreteCint_T[[i]] <- stringToMxAlgebra(x1, name=x2);
      }
      newDiscreteCint_T

    } # END method == 2
  }

  ######################################################################################################################
  ################################################  CHANGE T0VARs ######################################################
  ######################################################################################################################
  {
    oldT0VarBase <- OpenMxModel$matrices$T0VARbase; oldT0VarBase
    oldT0VarChol <- OpenMxModel$algebras$T0VARchol; oldT0VarChol
    oldT0Var <- OpenMxModel$algebras$T0VAR; oldT0Var

    newT0VarBase <- list() # make list to handle both single and multiple T0var matrices
    newT0VarBase[[1]] <- oldT0VarBase
    if (!(any(is.na(T0varValuesForAll)))) {
      tmp1 <- which(newT0VarBase[[1]]$labels %in% names(T0varValuesForAll)); tmp1
      tmp2 <- which(names(T0varValuesForAll) %in% newT0VarBase[[1]]$labels); tmp2
      #newT0VarBase[[1]]$values[tmp1] <- T0varValuesForAll; newT0VarBase[[1]]$values
      newT0VarBase[[1]]$values[tmp1] <- T0varValuesForAll[tmp2]; newT0VarBase[[1]]$values
      if (!(is.null(fixedmodel))) {
        tmp2 <- suppressWarnings(as.numeric(fixedmodel$T0VAR)); tmp2
        tmp3 <- which(!(is.na(tmp2))); tmp3
        tmp3 <- tmp3[!is.na(tmp3)]; tmp3
        newT0VarBase[[1]]$free[tmp3] <- FALSE
      }
    }
    newT0VarBase
    newT0VarBaseTemplate <- newT0VarBase

    newT0VarChol <- oldT0VarChol
    newT0Var <- oldT0Var

    if (any(!(is.na(T0varValuesForGroup[[1]])))) {
      newT0VarChol <- list()
      newT0Var <- list()
      for (i in 1:length(uniqueGroups)) {
        ### T0VARbase
        newT0VarBase[[i]] <- newT0VarBaseTemplate[[1]]
        tmp1 <- which(newT0VarBase[[i]]$free); tmp1
        # labels for parameters
        if (length(tmp1) > 0) {
          newLabels <- paste0(oldT0VarBase$labels[tmp1], "_G", i); newLabels
          newT0VarBase[[i]]$labels[tmp1] <- newLabels; newT0VarBase[[i]]$labels
        } else {
          newLabels <- paste0(oldT0VarBase$labels, "_G", i); newLabels
          newT0VarBase[[i]]$labels <- newLabels; newT0VarBase[[i]]$labels
        }
        # change labels back to unmoderated model if effect is not moderated by groupd:
        newValues <- T0varValuesForGroup[[i]]; newValues
        tmp1 <- which(!is.na(newValues)); tmp1
        newT0VarBase[[i]]$values[tmp1] <- newValues[tmp1]; newT0VarBase[[i]]$values
        newT0VarBase[[i]]$labels[fixedmodel$T0VAR == "groupfixed"] <- oldT0VarBase$labels[fixedmodel$T0VAR == "groupfixed"]
        # label for matrix
        newLabels <- paste0("T0VARbase", "_", i)
        newT0VarBase[[i]]$name <- newLabels; newT0VarBase[[i]]$name
      }
    } else {
      for (i in 1:length(uniqueGroups)) {
        newT0VarBase[[i]] <- newT0VarBase[[1]]
        # label for matrix
        newLabels <- paste0("T0VARbase", "_", i)
        newT0VarBase[[i]]$name <- newLabels; newT0VarBase[[i]]$name
      }
    }
    newT0VarBase

    if (length(uniqueModerators > 1)) {
      newT0VarBase <- rep(newT0VarBase, length(uniqueModerators)) # multiply list of drift matrices
      counter <- 0
      #for (i in 1:(length(newT0VarBase)/length(uniqueModerators))) { # assign new labels to drift coefficients.
      #for (j in 1:length(uniqueModerators)) {
        for (j in 1:length(uniqueModerators)) {
          for (i in 1:(length(newT0VarBase)/length(uniqueModerators))) { # assign new labels to drift coefficients.
            counter <- counter +1
          #newT0VarBase[[counter]]$labels[modmodel$T0VAR == "moderated"] <- paste0(oldT0VarBase$labels[modmodel$T0VAR == "moderated"], "_", uniqueModerators[j])
          newT0VarBase[[counter]]$labels[modmodel$T0VAR == "moderated"] <- paste0(oldT0VarBase$labels[modmodel$T0VAR == "moderated"], "_", unique(newModeratorID)[j])
          newLabels <- newT0VarBase[[counter]]$labels
          newLabels[grep("NA", newLabels)] <- NA
          newT0VarBase[[counter]]$labels <- newLabels
          newLabels <- paste0("T0VARbase", "_", counter)
          newT0VarBase[[counter]]$name <- newLabels
          newT0VarBase[[counter]]$values[modmodel$T0VAR == "moderated"] <- T0varValuesForModerator[[j]][modmodel$T0VAR == "moderated"]
        }
      }
    }
    newT0VarBase

    # get correct start values for diffusion base
    if (fastIterativeCoTiMA == TRUE) {
      for (k in 1:length(newT0VarBase)) {
        counter <- 1
        for (l in 1:(nlatents)) {
          for (m in l:(nlatents)) {
            tmp1 <- startValues[grep(paste0("T0var"), names(startValues))]; tmp1
            tmp2 <- tmp1[grep(paste0("_G", k), names(tmp1))]; tmp2
            newT0VarBase[[k]]$values[m,l] <- tmp2[counter]
            counter <- counter + 1
          }
        }
      }
    }

    if (method == 1) {
      ### T0VARbaseBig
      bigString <-unlist(lapply(newT0VarBase, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "T0VARbaseBig"
      T0VarBaseBig <-  stringToMxAlgebra(x1, name=x2); T0VarBaseBig

      ### T0VAR
      x1 <- paste0("T0VARbaseBig[ ( (grpXmodID_T[1, ", 1, "] - 1) * nlatent + 1):(grpXmodID_T[1, ", 1, "] * nlatent), 1:nlatent]"); x1
      #x2 <- paste0("T0VAR")
      x2 <- paste0("T0VARbase")
      T0Var <- stringToMxAlgebra(x1, name=x2); T0Var

      ### T0VARchol
      newT0VarChol <- oldT0VarChol; newT0VarChol

      ### T0VAR
      newT0Var <- oldT0Var; newT0Var

      ### discreteT0VARbig
      ### discreteT0VAR
    }

    if (method == 2) {
      #T0VARchol
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("(vec2diag(exp(diag2vec(T0VARbase", "_", i, "))) + T0VARbase", "_", i, " - vec2diag(diag2vec(T0VARbase", "_", i, ")))"); x1 # test
        x2 <- paste0("T0VARchol", "_", i); x2
        newT0VarChol[[i]] <- stringToMxAlgebra(x1, name=x2)
      }
      newT0VarChol

      # T0VAR
      for (i in 1:(length(uniqueModerators) * length(uniqueGroups))) {
        x1 <- paste0("(T0VARchol", "_", i, " %*% t(T0VARchol", "_", i, "))"); x1
        x2 <- paste0("T0VAR", "_", i); x2
        newT0Var[[i]] <- stringToMxAlgebra(x1, name=x2)
      }
      newT0Var

      # New algebras for discreteT0VARbig
      bigString <- unlist(lapply(newT0Var, function(extract) extract$name)); bigString
      helper  <- paste(as.character(bigString), sep="' '", collapse=", "); helper
      x1 <- paste0("rbind(", helper, ")")
      x2 <- "discreteT0VARbig"
      newDiscreteT0VarBig <-  stringToMxAlgebra(x1, name=x2); newDiscreteT0VarBig

      # New algebra for discreteT0VAR (required fro matrix S)
      x1 <- paste0("discreteT0VARbig[ (grpXmodID_T[1, ", 1, "] * nlatent + 1 - nlatent):(grpXmodID_T[1, ", 1, "] * nlatent), 1:nlatent]"); x1
      x2 <- paste0("T0VAR")
      newDiscreteT0Var <- stringToMxAlgebra(x1, name=x2); newDiscreteT0Var

    }
  }


  ######################################################################################################################
  #################################################### Create Model ####################################################
  ######################################################################################################################

  OpenMxModel2 <- OpenMxModel
  #OpenMxModel2$matrices$DIFFUSIONbase
  matricesToDelete <- c("DRIFT", "DIFFUSIONbase", "asymDIFFUSION", "T0VARbase", "CINT")
  algebrasString <-unlist(lapply(OpenMxModel2$algebras, function(extract) extract$name)); algebrasString
  algebrasToDelete <- algebrasString[grep("DIFFUSION", algebrasString)]; algebrasToDelete
  algebrasToDelete <- c(algebrasToDelete, algebrasString[grep("discrete", algebrasString)]); algebrasToDelete
  algebrasToDelete <- c(algebrasToDelete, algebrasString[grep("CINT", algebrasString)]); algebrasToDelete
  algebrasToDelete <- c(algebrasToDelete, algebrasString[grep("T0VAR", algebrasString)]); algebrasToDelete
  algebrasToDelete <- c(algebrasToDelete, algebrasString[grep("DRIFT", algebrasString)]); algebrasToDelete

  OpenMxModel2 <- OpenMx::mxModel(OpenMxModel2, matricesToDelete, remove=TRUE) # remove old matrices
  OpenMxModel2 <- OpenMx::mxModel(OpenMxModel2, algebrasToDelete, remove=TRUE) # remove old algebras
  if (method == 1) {
    elementsToAdd <- list(grpXmodID_T,
                          newDrift,
                          driftBig,
                          drift,
                          newInvDrift,
                          newDriftHatch,
                          discreteDrift_im,
                          newDiscreteDriftBig,
                          newDiscreteDrift_T,
                          newDiffusionBase,
                          diffusionBaseBig,
                          diffusionBase,
                          newDiffusionChol,
                          newDiffusion,
                          newAsymDiffusion,
                          newAsymDiffusionAlg,
                          discreteDiffusion_im,
                          newDiscreteDiffusionBig,
                          newDiscreteDiffusion_T,
                          newCint,
                          cintBig,
                          cint,
                          discreteCint_im,
                          newDiscreteCintBig,
                          newDiscreteCint_T,
                          newT0VarBase,
                          T0VarBaseBig,
                          T0Var,
                          newT0VarChol,
                          newT0Var)
  }

  if (method == 2) {
    elementsToAdd <- list(grpXmodID_T,
                          newDrift,
                          #driftBig,
                          #drift,
                          newInvDrift,
                          newDriftHatch,
                          discreteDrift_im,
                          newDiscreteDriftBig,
                          newDiscreteDrift_T,
                          newDiffusionBase,
                          #diffusionBaseBig,
                          #diffusionBase,
                          newDiffusionChol,
                          newDiffusion,
                          newAsymDiffusion,
                          newAsymDiffusionAlg,
                          discreteDiffusion_im,
                          newDiscreteDiffusionBig,
                          newDiscreteDiffusion_T,
                          newCint,
                          #cintBig,
                          #cint,
                          discreteCint_im,
                          newDiscreteCintBig,
                          newDiscreteCint_T,
                          newT0VarBase,
                          #T0VarBaseBig,
                          newDiscreteT0VarBig,
                          newDiscreteT0Var,
                          #T0Var,
                          newT0VarChol,
                          newT0Var)
  }

  OpenMxModel3 <- OpenMx::mxModel(OpenMxModel2, elementsToAdd)
  if (fastIterativeCoTiMA == TRUE) {
    OpenMxModel3 <- OpenMx::mxOption(OpenMxModel3, "Calculate Hessian", "No")
    OpenMxModel3 <- OpenMx::mxOption(OpenMxModel3, "Standard Errors"  , "No")
  }
  #if (fastIterativeCoTiMA == TRUE) {
  #  tryCatch(OpenMxModel3Fit <- OpenMx::mxRun(OpenMxModel3), error = function(e) e)
  #  if (!(exists("OpenMxModel3Fit"))) OpenMxModel3Fit <- OpenMx::mxTryHard(OpenMxModel3, extraTries = extraTries)
  #} else {
  if (tryHard == TRUE ) {
    OpenMxModel3Fit <- OpenMx::mxTryHard(OpenMxModel3, extraTries = extraTries,
                                 checkHess = checkHess) #,
                                 #wtgcsv=wtgcsv, #initialGradientIterations=initialGradientIterations,
                                 #finetuneGradient=finetuneGradient)
  } else {
    #OpenMxModel3 <- OpenMx::mxOption(OpenMxModel3, "Nudge zero starts"  , FALSE)
    OpenMxModel3Fit <- OpenMx::mxRun(OpenMxModel3)
    if (OpenMxModel3Fit$output$conditionNumber %in% c(5, 6)) {
      cat("The solution was either not convex or had nonzero gradient. Try again")
      cat(" ", "", sep="\n")
      cat("using the current estimates as starting values with mxTryHard")
    }
  }
  #}
  return(OpenMxModel3Fit)
}
