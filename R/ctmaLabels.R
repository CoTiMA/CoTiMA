#' ctmaLabels
#'
#' @description used for consistent labeling of names and parameters
#'
#' @param n.latent n.latent
#' @param n.manifest n.manifest
#' @param lambda lambda
#' @param manifestVars manifestVar
#' @param drift drift
#' @param diff diffusion
#' @param invariantDrift invariantDrift
#' @param moderatedDrift moderatedDrift
#' @param equalDrift equalDrift
#' @param T0means T0means
#' @param manifestMeans manifestMeans
#'
#' @return returns consistently named parameters (e.g., "V1toV2") as well es their symbolic values, which are used to fix or free
#' parameters when fitting a 'CoTiMA' model
#'
ctmaLabels <- function(
    n.latent=NULL,
    n.manifest=0,
    lambda=NULL,
    manifestVars=NULL,
    drift=NULL,
    diff=NULL,
    invariantDrift=NULL,
    moderatedDrift=NULL,
    equalDrift=NULL,
    T0means=0,
    manifestMeans=0)
{

  n.var <- max(c(n.manifest, n.latent)); n.var

  driftNames <- diffNames <- c()
  driftParams <- diffParams <- c()
  for (i in 1:(n.latent)) {
    for (j in 1:(n.latent)) {
      driftNames <- c(driftNames, paste0("V",i,"toV", j))
      if (i != j) diffNames <- c(diffNames, paste0("diff_eta", j, "_eta", i))
      if (i == j) diffNames <- c(diffNames, paste0("diff_eta", j))
    }
  }

  #driftNames <- c(t(matrix(driftNames, n.latent))); driftNames
  driftParams <- driftNames; driftParams
  diffParams <- diffNames; diffParams
  driftFullNames <- driftNames; driftFullNames
  diffFullNames <- diffNames; diffFullNames
  driftFullNames <- c(t(matrix(driftFullNames, n.latent))); driftFullNames
  diffFullNames <- c(t(matrix(diffFullNames, n.latent))); diffFullNames

  if (!(is.null(drift))) {
    tmp1 <- which(!(driftParams %in% drift)); tmp1
    # changed 11. Aug. 2022
    #driftParams[tmp1] <- "0"
    # undo change on 12.7.23
    #driftParams[tmp1] <- drift; driftParams
    tmp2 <- which(driftParams %in% drift); tmp2
    driftParams[tmp1] <- drift[tmp1]; driftParams
    #tmp1 <- which(drift == 0); tmp1
    tmp1 <- suppressWarnings(which(!(is.na(as.numeric(drift))))); tmp1
    if (length(tmp1) > 0) driftNames <- driftNames[-tmp1]
  }

  # added 11. Aug. 2022
  if (!(is.null(diff))) {
    tmp1 <- which(!(diffParams %in% drift)); tmp1
    diffParamsParams[tmp1] <- diff; driftParams
    tmp1 <- which(diff == 0); tmp1
    if (length(tmp1) > 0) diffNames <- diffNames[-tmp1]
  }

  # backup full names for labeling output later
  tmp1 <- tmp2 <- c()
  if (!(is.null(drift))) {
    # check validity of user-provided drift names
    # taken out on 11. Aug 2022
    #tmp1 <- which(c(driftParams) %in% driftNames); tmp1
    #tmp2 <- which(c(driftParams) == "0"); tmp2
    #if ( (length(tmp1)+length(tmp2)) != length(driftParams) ) {
    #ErrorMsg <- "\nDrift names provided by user do not match requirements.\nThey should be of the type V1toV2 or just 0. \nGood luck for the next try!"
    #stop(ErrorMsg)
    #}
    driftNames <- drift # replace
  }

  if (is.null(invariantDrift)) invariantDrift <- driftNames
  invariantDriftParams <- invariantDriftNames <- invariantDrift; invariantDriftParams
  tmp1 <- which(driftNames %in% invariantDriftNames); tmp1
  driftNames[tmp1] <- paste0(driftNames[tmp1], " (invariant)"); driftNames

  if (!(is.null(moderatedDrift))) {
    moderatedDriftNames <- moderatedDrift
    if (length(moderatedDriftNames) < 2)  {
      if (moderatedDriftNames == "all") moderatedDriftNames <- driftFullNames
    }
  } else {
    moderatedDriftNames <- NULL
  }

  equalDriftParams <- equalDriftNames <- equalDrift; equalDriftParams
  if (!(is.null(equalDrift))) {
    #tmp1 <- which(driftNames %in% equalDriftNames); tmp1
    #driftNames[tmp1] <- paste0(driftNames[tmp1], " (equal)"); driftNames
    tmp1 <- which(driftParams %in% equalDriftNames); tmp1
    tmpNames <- driftParams[tmp1[1]]; tmpNames
    for (t in 2:length(tmp1)) tmpNames <- c(tmpNames, "_eq_", driftParams[tmp1[t]])
    tmpNames <- paste(tmpNames,collapse=""); tmpNames
    driftParams[tmp1] <- tmpNames; driftParams
    tmpNames <- paste(tmpNames, "(invariant)", collapse=""); tmpNames
    driftNames[tmp1] <- tmpNames; driftNames
    # ensure that equal params are invarant params, too.
    invariantDriftNames <- unique(c(invariantDrift, equalDriftNames))
    tmp1 <- c()
    for (t in 1:length(equalDrift)) tmp1 <- c(tmp1,  grep(equalDrift[t], invariantDriftNames))
    invariantDriftParams[tmp1] <- invariantDriftNames[tmp1] <- driftParams[tmp1]
    invariantDriftNames <- paste0(invariantDriftNames, " (invariant)"); invariantDriftNames
  }
  driftNames; driftParams; invariantDriftNames
  invariantDriftParams


  tmp0 <- matrix(diffParams, n.latent); tmp0
  tmp0[upper.tri(tmp0, diag=FALSE)] <- 0; tmp0
  diffParams <- c(tmp0); diffParams

  # Adaptations if latent variables are measured with multiple indicators
  # loadings
  if (n.manifest > n.latent) {
    LAMBDA <- lambda
  } else {
    LAMBDA=diag(n.latent)
  }

  # error variances
  if(!(is.null(manifestVars))) manifestVarsParams <- manifestVars else manifestVarsParams <- 0

  # T0 variance
  T0VAR <- "auto"
  skip <- 1
  if (skip == 1) {
    tmp1 <- which(LAMBDA == "0")
    tmp2 <- which(LAMBDA == "1")
    if ( (length(tmp1) + length(tmp2)) < n.var * n.latent ) {
      tmp3 <- suppressWarnings(matrix(as.numeric(LAMBDA), nrow=nrow(LAMBDA))); tmp3
      tmp3 <- as.data.frame(tmp3)
      targetVar <- which(is.na(colSums(tmp3))); targetVar
      T0VAR <- matrix(0, n.latent, n.latent); T0VAR
      for (k in 1:n.latent) {
        for (m in k:n.latent) {
          T0VAR[k, m] <- paste0("T0VAR", k, m)
          if ( (k == targetVar) & (k == m) ) T0VAR[k, m] <- 1
        }
      }
      T0VAR <- t(T0VAR)
    } else {
      T0VAR <- "auto"
    }
  }

  # manifest means
  # CHD 9.6.2023
  if (manifestMeans == 'auto') {
    MANIFESTMEANS <- 'auto'
    print(MANIFESTMEANS)
  } else {
    MANIFESTMEANS <- 0
    if (!(is.null(invariantDrift))) {
      if ( (length(tmp1) + length(tmp2)) < n.var * n.latent ) {
        # added 16. Aug 2022 (if else)
        if (manifestMeans == 0) {
          MANIFESTMEANS <- rep("0", n.manifest); MANIFESTMEANS
          targetVar <- which(is.na(rowSums(tmp3))); targetVar
          MANIFESTMEANS[targetVar] <- paste0("mean_", targetVar); MANIFESTMEANS
        } #else {
        #  MANIFESTMEANS <- manifestMeans
        #}
      } else {
        # changed 16. Aug 2022
        #MANIFESTMEANS <- 0
        MANIFESTMEANS <- manifestMeans
        #}
      }
    }
    print(MANIFESTMEANS)
  }

  #T0Means
  if (!(is.null(invariantDrift))) {
    T0meansParams <- T0means
  }

  results <- list(driftNames=driftNames,
                  driftFullNames=driftFullNames,
                  driftParams=driftParams,
                  diffNames=diffNames,
                  diffFullNames=diffFullNames,
                  diffParams=diffParams,
                  invariantDriftNames=invariantDriftNames,
                  invariantDriftParams=invariantDriftParams,
                  moderatedDriftNames=moderatedDriftNames,
                  equalDriftNames=equalDriftNames,
                  equalDriftParams=equalDriftParams,
                  lambdaParams=LAMBDA,
                  T0VARParams=T0VAR,
                  manifestMeansParams=MANIFESTMEANS,
                  T0meansParams=T0meansParams,
                  manifestVarsParams=manifestVarsParams)

  invisible(results)
}
