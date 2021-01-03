#' ctmaLabels
#'
#' @param n.latent n.latent
#' @param n.manifest n.manifest
#' @param lambda lambda
#' @param manifestVar manifestVar
#' @param drift drift
#' @param invariantDrift invariantDrift
#' @param moderatedDrift moderatedDrift
#' @param equalDrift equalDrift
#'
#' @export ctmaLabels
#'
ctmaLabels <- function(
  n.latent=NULL,
  n.manifest=0,
  lambda=NULL,
  manifestVar=NULL,
  drift=NULL,
  invariantDrift=NULL,
  moderatedDrift=NULL,
  equalDrift=NULL
) {

  #n.latent=2
  #n.manifest=0
  #lambda=NULL
  #manifestVar=NULL
  #drift=c("V1toV1", "V1toV2", "V2toV2")
  #invariantDrift=NULL
  #moderatedDrift=NULL
  #equalDrift=NULL

  n.var <- max(c(n.manifest, n.latent)); n.var

  driftNames <- diffNames <- c()
  driftParams <- diffParams <- c()
  for (i in 1:(n.latent)) {
    for (j in 1:(n.latent)) {
      driftNames <- c(driftNames, paste0("V",i,"toV", j))
      if (i != j) diffNames <- c(diffNames, paste0("diff_eta", j, "_eta", i)) else diffNames <- c(diffNames, paste0("diff_eta", j))
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
    driftParams[tmp1] <- "0"
    driftNames <- driftNames[-tmp1]
  }

  # backup full names for labelling output later
  if (!(is.null(drift))) {
    # check validity of user-provided drift names
    tmp1 <- which(c(driftParams) %in% driftNames); tmp1
    tmp2 <- which(c(driftParams) == "0"); tmp2
    if ( (length(tmp1)+length(tmp2)) != length(driftParams) ) {
      cat(crayon::red$bold("Drift names provided by user do not match requirements.", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::blue("They should be of the type V1toV2 or just 0.", "\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      stop("Good luck for the next try!")
    }
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
  tmp1 <- which(driftNames %in% equalDriftNames); tmp1
  driftNames[tmp1] <- paste0(driftNames[tmp1], " (equal)"); driftNames
  #driftParams <- c(t(matrix(driftParams, n.latent))); driftParams
  #driftNames <- c(t(matrix(driftNames, n.latent))); driftNames

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
  if(!(is.null(manifestVar))) manifestVarParams <- manifestVar else manifestVarParams <- 0

  # T0 variance
  T0VAR <- "auto"
  skip <- 0
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

  # manifest means (as interindividually varying params, which replaces error auto-correlations)
  MANIFESTMEANS <- 0
  skip <- 0
  if (skip == 1) {
    if ( (length(tmp1) + length(tmp2)) < n.var * n.latent ) {
      MANIFESTMEANS <- rep("0", n.manifest); MANIFESTMEANS
      targetVar <- which(is.na(rowSums(tmp3))); targetVar
      MANIFESTMEANS[targetVar] <- paste0("mean_", targetVar); MANIFESTMEANS
    } else {
      MANIFESTMEANS <- 0
    }
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
                 manifestVarParams=manifestVarParams)

  invisible(results)
}

