#' ctmaGenData
#'
#' @description Generates data from lists of parameters (drift, diffusion etc). Experimental!!
#'
#' @param activeDirectory defines active directory where files are saved. No default.
#' @param burnin vector of initial time points to be deleted (default = 0)
#' @param cint list of cint matrices. By default (NULL), cint matrices will be used that create a steady-state (i.e., means at all time points = T0means)
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use.
#' @param diff list of diffusion matrices. By default (NULL), diffusion matrices will be used that create a steady-state (i.e., covariance at all time points = T0var)
#' @param digits number of digits used for rounding (in outputs).
#' @param doPar parallel generating of data. if TRUE, data are generated in coresToUse parallel loops during which no output is generated (screen remains silent).
#' @param drift list of drift matrices. No default (all = NULL).
#' @param empirical whether (default) or not generated data that should be independent as generally assumed (e.g., diffusions at different time points, T0var, etc)  are truly independent. Allow for exact estimation of parameters (no random variance; not useful for MC simulations). May require large sample sizes.
#' @param envir environment where objects should be extracted too. NULL by default. Typically one would use the global environment (envir = globalenv()). Has to be specified if ctmaExtract is set to TRUE.
#' @param sampleSizes vector of sample sizes. Default = 100.
#' @param lambda list of matrices. By default all are diagonal matrices with 1 in the diagonal.
#' @param latentNames names for latent variables (default = NULL using generic names)
#' @param manifestMeans list of manifest mean matrices. By default all are = 0.
#' @param manifestVars list of manifest error (co-)variances matrices. By default all are = 0.
#' @param manifestNames  names for manifest variables (default = NULL using generic names)
#' @param modValues list of moderator values (possible used to generate the list of drift matrices provided). By default all are = 0.
#' @param missings proportion of missings (default = 0, which does not delete any value)
#' @param n.latent number of latent variables of the model No default (all = NULL).
#' @param n.manifest number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).
#' @param randomIntercepts list of (co-)variances matrices of TIpreds (traits). By default all are = 0.
#' @param T0means list of Time 0 mean levels. By default all are = 0.
#' @param T0var list of Time 0 (co-)variance matrices. No default (all = NULL).
#' @param TIpreds list of time-independent predictors, e.g., the moderators that were used to create the drift matrices. No default (all = NULL).
#' @param tpoints vector of number of tpoints to be generated (default = 10).
#' @param tpointTargets list of vectors of tpoints to be selectd (default = burnin:tpoints).
#' @param useRawData if TRUE (or FALSE) and ctmaExtract is also TRUE, creates rawData (or empcov) objects required for CoTiMA into the environment specified using the argument envir.
#' @param ctmaExtract if TRUE (default = FALSE) uses ctmaExtract to extract objects required for CoTiMA into the environment specified using the argument envir. Requires the argument useRawData to be set to TRUE or FALSE.
#'
#' @importFrom parallel detectCores makeCluster
#' @importFrom ctsem ctDeintervalise ctLongToWide ctIntervalise ctWideToLong ctModel ctStanFit ctExtract ctCollapse
#' @importFrom doParallel registerDoParallel
#' @importFrom OpenMx expm
#' @importFrom Matrix bdiag
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @export ctmaGenData
#'
#' @examples
#' # Fit a ctsem model to all three primary studies summarized in
#' # CoTiMAstudyList_3 and save the three fitted models
#' \dontrun{
#' CoTiMAInitFit_3 <- ctmaInit(primaryStudies=CoTiMAstudyList_3,
#'                             n.latent=2,
#'                             checkSingleStudyResults=FALSE,
#'                             activeDirectory="/Users/tmp/") # adapt!
#' summary(CoTiMAInitFit_3)
#' }
#'
#' @return ctmaGenData returns a list containing ...
#'
ctmaGenData <- function(
    activeDirectory = NULL,
    burnin = 0,
    cint = NULL,
    coresToUse = 2,
    ctmaExtract=FALSE,
    diff = NULL,
    digits = 4,
    doPar = FALSE,
    drift = NULL,
    empirical = TRUE,
    envir = NULL,
    sampleSizes = 100,
    lambda = NULL,
    latentNames = NULL,
    manifestMeans = NULL,
    manifestVars = NULL,
    manifestNames = NULL,
    missings = NULL,
    modValues = 0,
    n.latent = NULL,
    n.manifest = NULL,
    randomIntercepts = NULL,
    T0means = 0,
    T0var = NULL,
    TIpreds = NULL,
    tpoints = 10,
    tpointTargets = NULL,
    useRawData=NULL
)
{

  original.options <- getOption("scipen"); original.options
  options(scipen = 999); options("scipen") # turn scientific notation off.
  on.exit(options(scipen = original.options))  # scientific notation as user's original

  # Check if not yet working arguments are set #####
  {
    if (!(is.null(TIpreds))) {
      ErrorMsg <- "\n Argument TIpreds not yet possible"
      stop(ErrorMsg)
    }


    if (empirical == FALSE) {
      if (burnin < 100) {
        Msg <- "The argument empirical == FALSE was used and burnin < 100. Longer burnin phases are recommended to achieve a steady state.\n"
        message(Msg)
      }
    }
  }

  # Check if arguments are properly used #####
  { # start checks
    if (utils::packageDescription("ctsem")$Version > "3.10.2") type <- "ct" else type <- "stanct"

    if (is.null(n.latent)) {
      Msg <- "n.latent not specified. I infer n.latent from the dimensions of the drift matrix.\n"
      message(Msg)
      n.latent <- ncol(drift[[1]]); n.latent
    }

    if (is.null(drift)) {
      ErrorMsg <- "\n No drift matrix (drift) specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.list(drift))) {
      ErrorMsg <- "\nThe drift argument has to be a list. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (length(unique(unlist(lapply(drift, function(x) ncol(x))))) != 1) {
      ErrorMsg <- "\nThe drift matrices supplied have different dimensions. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.null(manifestVars))) {
      if (!(is.list(manifestVars))) {
        ErrorMsg <- "\nThe manifestVars argument has to be a list (of matrices). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(manifestVars, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe manifestVars matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( ncol(manifestVars[[1]]) != nrow(lambda[[1]]) ) {
        ErrorMsg <- "\nThe manifestVars matrices supplied have a different dimensions than the rows of the lambda matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (!all(unlist((lapply(manifestVars, function(x) isSymmetric(x, tol = 1e-8)))))) {
        ErrorMsg <- "\nAt least one of the manifestVars matrices supplied is not symmetric. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (!(is.null(n.manifest))) {
      if (is.null(lambda)) {
        ErrorMsg <- "\n You specified n.manifest. I also need a list of lambda matrices."
        stop(ErrorMsg)
      }
    }

    if (!(is.null(randomIntercepts))) {
      if (!(is.list(randomIntercepts))) {
        ErrorMsg <- "\nThe randomIntercepts argument has to be a list (of matrices). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(randomIntercepts, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe randomIntercepts matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( ncol(randomIntercepts[[1]]) != nrow(drift[[1]]) ) {
        ErrorMsg <- "\nThe randomIntercepts matrices supplied have a different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (!all(unlist((lapply(randomIntercepts, function(x) isSymmetric(x, tol = 1e-8)))))) {
        ErrorMsg <- "\nAt least one of the randomIntercept matrices supplied is not symmetric. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (is.null(randomIntercepts)) {
      randomIntercepts <- rep(list(matrix(0, n.latent, n.latent)), length(drift))
    }

    if (is.null(tpointTargets)) {
      tpointTargets <- rep(list(burnin:(tpoints-1)), length(drift))
    }

    if (!(is.null(tpointTargets))) {
      if (!(is.list(tpointTargets))) {
        ErrorMsg <- "\nThe tpointTargets argument has to be a list (of vectors). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( length(tpointTargets) != length(drift) ) {
        ErrorMsg <- "\nThere have to be as many tpointTargets in the list as drift matrices supplied. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (any(sapply(tpointTargets, function(v) any(v %% 1 != 0)))) {
      Msg <- "some time points specified in tpointTargets are not interger values. I will round them to 0 digits.\n"
      message(Msg)
      tmp <- which(sapply(tpointTargets, function(v) any(v %% 1 != 0)))
      for (i in tmp) tpointTargets[[i]] <- round(tpointTargets[[i]], 0)
    }

    if (is.null(T0var)) {
      Msg <- "No T0var matrix specified. I set all T0var to Identity matrices! \n"
      T0var <- rep(list(diag(1, n.latent, n.latent)), length(drift))
      message(Msg)
    }

    if (!(is.list(T0var))) {
      ErrorMsg <- "\nThe T0var argument has to be a list. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (length(T0var) != length(drift)) {
      ErrorMsg <- "\nThe number of T0var matrices provided (T0var) does not match the number of drift matrices provided (drift). \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (length(unique(unlist(lapply(T0var, function(x) ncol(x))))) != 1) {
      ErrorMsg <- "\nThe T0var matrices supplied have different dimensions. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!all(unlist((lapply(T0var, function(x) isSymmetric(x, tol = 1e-8)))))) {
      ErrorMsg <- "\nAt least one of the T0var matrices supplied is not symmetric. \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    # trait means
    TRAITMEANS <- matrix(0, nrow=n.latent, ncol=1)

    if (is.null(n.manifest)) {
      if (is.null(lambda)) {
        Msg <- "n.manifest not specified. I assume n.manifest is equal to n.latent.\n"
        message(Msg)
        n.manifest <- n.latent
        Msg <- "lambda not specified. I assume all lambdas are diagonal matrices with 1 in the diagonals.\n"
        message(Msg)
        lambda <- rep(list(diag(1, n.latent, n.latent)), length(drift))
      }
    }

    if (is.null(n.manifest)) {
      if (!(is.null(lambda))) {
        Msg <- "n.manifest not specified but lambda. I infer n.manifest from the number of rows of lambda.\n"
        message(Msg)
        n.manifest <- nrow(lambda[[1]])
      }
    }


    if (!(is.null(lambda))) {
      if (!(is.list(lambda))) {
        ErrorMsg <- "\nThe lambda argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(lambda, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe lambda matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(lambda, function(x) nrow(x))))) != 1) {
        ErrorMsg <- "\nThe lambda matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(lambda, function(x) ncol(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe lambda matrices supplied assume more latents than included in the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (any(diff(unique(unlist(lapply(lambda, function(x) nrow(x))))) != 0) ) {
        ErrorMsg <- "\nThe lambda matrices supplied assume different numbers of manifest variables. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (!(is.null(n.manifest))) {
        n.manifest.tmp <- n.manifest
        n.manifest <- unique(unlist(lapply(lambda, function(x) nrow(x))) )[1]
        if (n.manifest.tmp != n.manifest) {
          ErrorMsg <- "n.manifest specified differs from the number of rows of lambda. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
      if (is.null(manifestVars)) {
        Msg <- "lambda matrices supplied but no manifestVars matrices. I assume all manifest error variances are 0.\n"
        message(Msg)
        manifestVars <- rep(list(diag(0, n.manifest, n.manifest)), length(drift))
      }
    }

    if (length(T0means[[1]]) != 1) {
      if (!(is.list(T0means))) {
        ErrorMsg <- "\nThe T0means argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(T0means, function(x) nrow(x))))) != 1) {
        ErrorMsg <- "\nThe T0means matrices supplied have different dimensions or were not in matrix form. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(T0means, function(x) nrow(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe T0means matrices supplied have different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    } else {
      if (T0means[[1]] == 0) {
        T0means <- rep(list(rep(0, n.latent)), length(drift))
      } else {
        print("Problem with T0means")
      }
    }

    if (!(is.null(manifestMeans))) {
      if (!(is.list(manifestMeans))) {
        ErrorMsg <- "\nThe manifestMeans argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(manifestMeans, function(x) length(x))))) != 1) {
        ErrorMsg <- "\nThe manifestMeans vectors supplied have different lengths \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(manifestMeans, function(x) length(x)))) != n.manifest) {
        ErrorMsg <- "\nThe lengths of the manifestMeans vectors supplied do not match n.manifest. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (is.null(manifestMeans)) {
      if (!(is.null(manifestVars))) {
        Msg <- "manifestVars specified but no manifestMeans. I assume all manifest means are 0.\n"
        message(Msg)
        manifestMeans <- rep(list(rep(0, n.manifest)), length(drift))
      }
    }


    if (!(is.null(diff))) {
      if (!(is.list(diff))) {
        ErrorMsg <- "\nThe diff argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(diff, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe diff matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(diff, function(x) ncol(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe diff matrices supplied have different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (!all(unlist((lapply(diff, function(x) isSymmetric(x, tol = 1e-8)))))) {
        ErrorMsg <- "\nAt least one of the diff matrices supplied is not symmetric. \nGood luck for the next try!"
        stop(ErrorMsg)
      }

    }

    if (is.null(diff)) {
      Msg <- "diff not specified. I will set all diffusion matrices to values leading to a steady-state of the (co-)variances among latents.\n"
      message(Msg)
    }


    if (!(is.null(cint))) {
      if (!(is.list(cint))) {
        ErrorMsg <- "\nThe cint argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      (length(unique(unlist(lapply(cint, function(x) ncol(x))))) != 1)
      if (length(unique(unlist(lapply(cint, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe cint matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(cint, function(x) nrow(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe cint matrices supplied have different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (is.null(cint)) {
      Msg <- "cint not specified. I will set all cint matrices to values leading to a steady-state of the means of latents.\n"
      message(Msg)
    }

    if (!(is.null(TIpreds))) {
      if (!(is.list(TIpreds))) {
        ErrorMsg <- "\nThe TIpreds argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(TIpreds, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe TIpreds supplied have different numbers of elements. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }


    if (is.null(activeDirectory)) {
      ErrorMsg <- "\nNo active directory has been specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if  (length(coresToUse) > 0) {
      if (coresToUse < 1)  coresToUse <- parallel::detectCores() + coresToUse
    }

    if (coresToUse >= parallel::detectCores()) {
      coresToUse <- parallel::detectCores() - 1
      Msg <- "No of coresToUsed was set to >= all cores available. Reduced to max. no. of cores - 1 to prevent crash.\n"
      message(Msg)
    }


    skip <- FALSE
    if (!skip) {
      if (doPar == TRUE) {
        myCluster <- parallel::makeCluster(coresToUse)
      } else {
        myCluster <- parallel::makeCluster(1)
      }
      on.exit(parallel::stopCluster(myCluster))
      doParallel::registerDoParallel(myCluster)
      '%dopar%' <- foreach::'%dopar%'
    }

    if (n.manifest > n.latent ) {
      if (is.null(lambda)) {
        ErrorMsg <- "\nManifest variables specified, but not matrix of loadings (lambda) specified, which is required. \nGood luck for the next try!"
        stop(ErrorMsg)
      } else {
        if ( (dim(lambda[[1]])[1] != n.manifest) | (dim(lambda[[1]])[2] != n.latent) ) {
          ErrorMsg <- "\nDimensions of loadings matrix (lambda) do not match n.latent and n.manifest. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
    }

    if (length(sampleSizes) == 1) {
      if (length(sampleSizes[[1]]) == 1) {
        sampleSizes <- as.list(rep(sampleSizes[[1]], length(drift)))
      }
    } else {
      if (!(is.list(sampleSizes))) {
        ErrorMsg <- "\nThe sampleSizes argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (length(randomIntercepts) == 1) {
      if (!(is.list(randomIntercepts))) {
        if (randomIntercepts == 0) {
          randomIntercepts <- replicate(length(drift), matrix(0, n.latent, n.latent), simplify = FALSE)
        }
      }
    } else {
      if (!(is.list(randomIntercepts))) {
        ErrorMsg <- "\nThe randomIntercepts argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(randomIntercepts, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe randomIntercepts matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(randomIntercepts, function(x) ncol(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe randomIntercepts matrices supplied have different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if ((n.manifest <= n.latent ) &  (is.null(lambda))) {lambda <- diag(n.latent)}

    if (!(is.null(diff))) {
      if ( length(diff) != length(drift) ) {
        ErrorMsg <- "\nThe number of diffusion matrices provided (diff) does not match the number of drift matrices provided (drift). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (is.list(diff) == FALSE) {
        ErrorMsg <- "\nThe diff argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (ctmaExtract == TRUE) {
      if (is.null(useRawData)) {
        ErrorMsg <- "\nSince ctmaExtract is set to true, the argument useRawData has to be set to TRUE or FALSE. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ((useRawData != TRUE) & (useRawData != FALSE)) {
        ErrorMsg <- "\nSince ctmaExtract is set to true, useRawData is possibly misspelled. It has to be set to TRUE or FALSE. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (is.null(envir)) {
        ErrorMsg <- "\nSince ctmaExtract is set to true, the argument envir has to be set to an environment (e.g., envir = globalenv()). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (length(modValues) == 1) {
      modValues <- replicate(length(drift), modValues, simplify = FALSE)
    }

    if (!(is.list(modValues))) {
      ErrorMsg <- "\nThe modValues argument has to be a list. \nGood luck for the next try!"
      stop(ErrorMsg)
    }
    if ( length(modValues) != length(drift) ) {
      ErrorMsg <- "\nThe number of modValues provided does not match the number of drift matrices provided (drift). \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (!(is.null(missings))) {
      if (!(is.list(missings))) {
        ErrorMsg <- "\nThe missings argument has to be a list (of vectors). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( length(missings) != length(drift) ) {
        ErrorMsg <- "\nThere have to be as many missings in the list as drift matrices supplied. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      #Msg <- "Missingness proportions were specified. Note that this does not yet work with manifests, only for latents.\n"
      #message(Msg)
      if (any(unlist(missings) < 0) | any(unlist(missings) > 1)) {
        ErrorMsg <- "\nAt least one of the missings proportions supplied is < 0 or > 1. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (is.null(missings)) {
      Msg <- "missings not specified. I will not generate missing data.\n"
      message(Msg)
      missings <- rep(list(0), length(drift))
    }

    # compute diffusion covariance (used later)
    resvar <- list()
    for (i in 1:length(drift)) {
      T1cov_impl <- expm(drift[[i]]) %*% ( T0var[[i]] ) %*% t(expm(drift[[i]])) + randomIntercepts[[i]] ; T1cov_impl
      resvar[[i]] <- as.matrix(T0var[[i]] + randomIntercepts[[i]] - T1cov_impl); resvar[[i]]
    }


    ## Compute diffusions to achieve steady state  ####
    if (is.null(diff)) {
      diff <- list()
      for (i in 1:length(drift)) {
        #i <- 1
        T1cov_impl <- expm(drift[[i]]) %*% ( T0var[[i]] ) %*% t(expm(drift[[i]])) + randomIntercepts[[i]] ; T1cov_impl
        #resvar <- as.matrix(T0var[[i]] + randomIntercepts[[i]] - T1cov_impl); resvar
        if (any(diag(resvar[[i]]) < 0)) {
          ErrorMsg <- paste0("\nCannot generate difussions for achieving steady states because negative diffusion variances are implied for Study ", i, ",",
                             "\nwhich may be due too large T0vars in combination with large randomIntercepts (variances).",
                             "\nT0var = ", unlist(T0var[[i]]), ",",
                             "\nrandomIntercepts = ", unlist(randomIntercepts[[i]]), ",",
                             "\nT1cov_impl = ", unlist(T1cov_impl), ",",
                             "\nresvar = ", unlist(resvar[[i]]), ",",
                             "\nGood luck for the next try!")
          stop(ErrorMsg)
        }
        #length(unique(round(abs(c(resvar)), 5)))
        if (n.latent > 1) {
          if (length(unique(round(abs(c(resvar[[i]])), 5))) == 1) {
            ErrorMsg <- paste0("\nCannot generate data because of singularity issues with Study ", i, ",",
                               "\nwhich may be due all drift elements having identical magnitudes (e.g., all -.1 or .1).",
                               "\nYour provided drift = ", paste0(c(drift[[i]]), collapse = " "),
                               "\nGood luck for the next try!")
            stop(ErrorMsg)
          }
        }
        DIAG <- diag(1, nrow(drift[[i]]), ncol(drift[[i]])); DIAG
        DRIFT_hatch <- drift[[i]] %x% DIAG + DIAG %x% drift[[i]]; DRIFT_hatch
        DIAG_hatch <- diag(1, nrow(DRIFT_hatch), ncol(DRIFT_hatch)); DIAG_hatch
        Q <- solve((expm(DRIFT_hatch * 1) - DIAG_hatch)) %*% DRIFT_hatch %*% c(resvar[[i]]); Q
        diff[[i]] <- matrix(Q, n.latent, n.latent); diff[[i]]
        #diff_dt[[i]] <- matrix(resvar, n.latent, n.latent); diff_dt[[i]] # done later

        if (any(eigen(diff[[i]])$values < 0)) {
          ErrorMsg <- paste0("\nCannot generate data because of negative eigenvalues in the dt diffusion matrix of Study ", i, ",",
                             "\nwhich may be due too large cross or auto effects.",
                             "\nDrift = ", unlist(drift[[i]]), ",",
                             "\nDriff = ", unlist(diff[[i]]), ",",
                             "\nGood luck for the next try!")
          stop(ErrorMsg)
        }
      }
    }

    if (!(is.null(cint))) {
      if (length(cint) != length(cint)) {
        ErrorMsg <- "\nThe number ofcint matrices provided (cint) does not match the number of drift matrices provided (drift). \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (!(is.list(cint)))  {
        ErrorMsg <- "\nThe cint argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    ## Compute cints to achieve steady state #####
    if (is.null(cint)) {
      cint <- list()
      for (i in 1:length(drift)) {
        asymCINT <- (T0means[[i]] + TRAITMEANS); asymCINT
        cint[[i]] <- -drift[[i]] %*% asymCINT; cint[[i]]
      }
    }

    ## check if all provided arguments assume the same number of studies ###ä
    # Objects to check
    objs <- list(
      diff            = diff,
      drift           = drift,
      sampleSizes     = sampleSizes,
      lambda          = lambda,
      manifestMeans   = manifestMeans,
      manifestVars    = manifestVars,
      missings        = missings,
      modValues       = modValues,
      #n.latent        = n.latent,
      #n.manifest      = n.manifest,
      randomIntercepts= randomIntercepts,
      T0means         = T0means,
      T0var           = T0var,
      #TIpreds         = TIpreds,
      #tpoints         = tpoints,
      tpointTargets   = tpointTargets
    )

    # 1. Check that all objects are lists
    are_lists <- sapply(objs, is.list)
    if (!all(are_lists)) {
      bad <- names(are_lists)[!are_lists]
      stop(sprintf("The following objects are not lists: %s", paste(bad, collapse = ", ")))
    }

    # 2. Check that all lists have identical lengths
    lens <- sapply(objs, length)
    if (length(unique(lens)) != 1) {
      stop(sprintf(
        "Lists do not all have the same length. Lengths are: %s",
        paste(names(lens), lens, sep = "=", collapse = ", ")
      ))
    }

    #message("All objects are lists and all have identical length = ", unique(lens))

    ## Define DT matrices and check for standardization #####
    drift_dt <- diff_dt <- cint_dt <- list()
    for (i in 1:length(drift))  {
      #i <- 1
      drift_dt[[i]] <- expm(drift[[i]]); drift_dt[[i]]

      DIAG <- diag(1, nrow(drift[[i]]), ncol(drift[[i]])); DIAG
      DRIFT_hatch <- drift[[i]] %x% DIAG + DIAG %x% drift[[i]]; DRIFT_hatch
      diff_dt[[i]] <- matrix((solve(DRIFT_hatch) %*% (expm(DRIFT_hatch * 1) - diag(n.latent^2)) %*% c(diff[[i]])), n.latent, n.latent); diff_dt[[i]]

      cint_dt[[i]] <- solve(drift[[i]]) %*% (expm(drift[[i]]) - diag(1, n.latent, n.latent)) %*% cint[[i]]; cint_dt[[i]]

      T1cov_impl <- expm(drift[[i]]) %*% ( T0var[[i]] ) %*% t(expm(drift[[i]])) + randomIntercepts[[i]] ; T1cov_impl

      if (!(all(round(diag(T1cov_impl + resvar[[i]]), 2) == 1))) {
        Msg <- paste0("\nThe steady state-variances rounded to 2 digits of Study ", i, " are not 1.0,",
                      "\nwhich implies so that the data are unstandardized and not useful for CoTiMA.",
                      "\nThey are = ",
                      paste0(diag(T1cov_impl + resvar[[i]]), collapse = " "),
                      "\nGood luck for the next try!")
        message(Msg)
      }
    }

    if (!(is.null(latentNames))) {
      if (length(latentNames) != ncol(drift[[1]]) ) {
        ErrorMsg <- "\nThe number of latentNames provided does not match the dimensions of the drift matrix. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (is.null(latentNames)) {
      latentNames <- LETTERS[1:n.latent]; latentNames
    }

    if (is.null(manifestNames)) {
      manifestNames <- LETTERS[1:n.manifest]; manifestNames
    }

    if (!(is.null(manifestNames))) {
      if (length(manifestNames) != nrow(lambda[[1]]) ) {
        ErrorMsg <- "\nThe number of manifestNames provided does not match the dimensions (rows) of the lambda matrix. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }


  } # end checks

  # Generate data ####
  #studies <- list()
  #for (i in 1:length(drift)) {
  studies <- foreach::foreach(i = 1:length(drift), .packages = 'ctsem', .export = c("ctmaExtract")) %dopar% {
    #i <- 1
    if (empirical == TRUE) {
      ## T0 data, diffusion, and traitvar all independent ####
      ### diffusions #####
      diff.var_tmp <- list(diff_dt[[i]]); diff.var_tmp
      diff.var <- rep(diff.var_tmp, tpoints); diff.var
      diff.var <- as.matrix(Matrix::bdiag(diff.var))
      rows1 <- nrow(diff.var); rows1
      ### T0var #####
      diffT0.var <- cbind(diff.var, matrix(0, ncol=n.latent, nrow=rows1))# diffvar & T0var
      cols1 <- ncol(diffT0.var); cols1
      diffT0.var <- rbind(diffT0.var, matrix(0, ncol=cols1, nrow=n.latent))
      diffT0.var[(cols1-n.latent+1):cols1, (cols1-n.latent+1):cols1] <- T0var[[i]]
      rows2 <- nrow(diffT0.var); rows2
      ### Traitvar ####
      diffT0Trait.var <- cbind(diffT0.var, matrix(0, ncol=n.latent, nrow=rows2)) # diffvar & T0var & traitvar
      cols2 <- ncol(diffT0Trait.var); cols2
      diffT0Trait.var <- rbind(diffT0Trait.var, matrix(0, ncol=cols2, nrow=n.latent))
      diffT0Trait.var[(cols2-n.latent+1):cols2, (cols2-n.latent+1):cols2] <- randomIntercepts[[i]]
      rows3 <- nrow(diffT0Trait.var); rows3
      ### manifestVar (measurement error) ####
      err.var_tmp <- list(manifestVars[[i]]); err.var_tmp
      err.var <- rep(err.var_tmp, tpoints); err.var
      err.var <- as.matrix(Matrix::bdiag(err.var))
      cols3 <- ncol(err.var); cols3
      rows4 <- nrow(err.var); rows4
      diffT0TraitERR.var <- matrix(0, nrow = rows3 + rows4, ncol = cols2 + cols3)
      diffT0TraitERR.var[1:rows3, 1:cols2] <- diffT0Trait.var
      diffT0TraitERR.var[(rows3+1):(rows3 + rows4), (cols2+1):(cols2 + cols3)] <- err.var
      cols3 <- ncol(diffT0TraitERR.var); cols3
      ## mvrnorm to create independent data for diffusions, T0 variables, random intercepts (traits), and measurement error
      #tmp <- length(c(rep(0, n.latent*tpoints), T0means[[i]], TRAITMEANS, rep(0, n.latent*tpoints)))
      tmp <- length(c(rep(0, n.latent*tpoints), T0means[[i]], TRAITMEANS, rep(0, n.manifest*tpoints)))
      if ( tmp > sampleSizes[[i]]) {
        ErrorMsg <- paste0("\n Requested sample size too small. It should be at least: sampleSizes =", tmp, " (or try empirical = FALSE).")
        stop(ErrorMsg)
      }
      allInit.dat <- MASS::mvrnorm(n=sampleSizes[[i]],
                                   mu=c(rep(0, n.latent*(tpoints)), T0means[[i]], TRAITMEANS, rep(0, n.manifest*tpoints)),
                                   Sigma = diffT0TraitERR.var, empirical=empirical)
      diff.dat <- allInit.dat[, (1:(n.latent*tpoints))]; dim(diff.dat)
      T0.dat <- allInit.dat[, (cols1 -n.latent+1):(cols1)]; dim(T0.dat)
      trait.dat <- allInit.dat[, (cols2 -n.latent+1):(cols2)]; dim(trait.dat)
      err.dat <- allInit.dat[, (cols2+1):(cols3)]; dim(err.dat)
      #
      data <- T0.dat
    } else {
      data  <- MASS::mvrnorm(n=sampleSizes[[i]], mu=T0means[[i]], Sigma = T0var[[i]], empirical=empirical)
      trait.dat <- MASS::mvrnorm(n=sampleSizes[[i]], mu=TRAITMEANS, Sigma = randomIntercepts[[i]], empirical=empirical)
      err.dat <- MASS::mvrnorm(n=sampleSizes[[i]], mu=rep(0, n.manifest*tpoints), Sigma = as.matrix(Matrix::bdiag( rep(list(manifestVars[[i]]), tpoints) )), empirical=empirical)
    }
    if (n.latent == 1) data <- as.matrix(data)

    # data that do not vary among cases and tpoints
    cint.dat <- matrix(t(as.matrix(cint_dt[[i]])),  nrow=sampleSizes[[i]], ncol=n.latent, byrow = T)
    manifestMeans.dat <- matrix(t(as.matrix(manifestMeans[[i]])),  nrow=sampleSizes[[i]], ncol=n.manifest, byrow = T)

    ## manifests at T0
    dataMM <- data + trait.dat
    dataMM <- dataMM %*% t(lambda[[i]])
    if (empirical == TRUE) {
      dataMM <- dataMM + err.dat[, (1:n.manifest)]
    } else {
      dataMM <- dataMM + MASS::mvrnorm(n=sampleSizes[[i]], mu=rep(0, n.manifest), Sigma = manifestVars[[i]], empirical=FALSE)
    }
    #round(cov(dataMM), 3)

    #### T1, T2, ... all subsequent Tpoints ####
    for (t in 1:(tpoints-1)) {
      #t <- 1
      tmpData <- data[, ((t-1)*n.latent+1):((t-1)*n.latent+n.latent)]
      if (n.latent == 1) tmpData <- as.matrix(tmpData)
      # apply drift
      tmp <- t(apply(tmpData, 1, function(x) as.matrix(drift_dt[[i]] %*% x)))
      if (n.latent == 1) tmp <- t(tmp)
      # add cint
      tmp <- tmp + cint.dat
      # add diffusion
      if (empirical == TRUE) {
        tmp <- tmp + diff.dat[, (n.latent*(t-1)+1):(n.latent*(t-1)+n.latent)]
      } else {
        tmp <- tmp + MASS::mvrnorm(n=sampleSizes[[i]], mu=rep(0, n.latent), Sigma = diff_dt[[i]], empirical=FALSE)
      }
      # add trait to latents (done later, otherwise it would be a cint rather than a trait)
      # tmp <- tmp + trait.dat

      # combine latents
      data <- cbind(data, tmp)

      # make manifests
      tmpMM <- tmp + trait.dat
      tmpMM <- tmpMM %*% t(lambda[[i]]) # manifest
      # add measurement error to manifests
      if (empirical == TRUE) {
        tmpMM <- tmpMM + err.dat[, (((t-1)*n.manifest+1):((t-1)*n.manifest+n.manifest))]
      } else {
        tmpMM <- tmpMM + MASS::mvrnorm(n=sampleSizes[[i]], mu=rep(0, n.manifest), Sigma = manifestVars[[i]], empirical=FALSE)
      }
      # add constant manifest means
      tmpMM <- tmpMM + manifestMeans.dat
      # combine manifests
      dataMM <- cbind(dataMM, tmpMM)
    }

    #### Add trait.dat to latents ####
    for (t in seq(1, ncol(data), n.latent)) data[, c(t:(t+(n.latent-1)))] <- data[, (t:(t+(n.latent-1)))] + trait.dat

    #### Add manifestMeans.dat to latents (in data without manifests) ####
    #for (t in seq(1, ncol(data), n.latent)) data[, c(t:(t+(n.latent-1)))] <- data[, c(t:(t+(n.latent-1)))] + manifestMeans.dat

    # label latents
    colnames(data) <- paste0(paste0(latentNames, "_T"), sort(rep(seq(0,(tpoints-1),1), n.latent)))
    data <- cbind(data, matrix(seq(0, tpoints-1, 1), nrow=nrow(data), ncol=tpoints, byrow=T))
    colnames(data)[(ncol(data)-tpoints+1):ncol(data)] <- paste0("T", sort(rep(seq(0,(tpoints-1),1))))
    #head(data); dim(data)

    # label manifests
    colnames(dataMM) <- paste0(paste0(manifestNames, "_T"), sort(rep(seq(0,(tpoints-1),1), n.manifest)))
    dataMM <- cbind(dataMM, matrix(seq(0, tpoints-1, 1), nrow=nrow(dataMM), ncol=tpoints, byrow=T))
    colnames(dataMM)[(ncol(dataMM)-tpoints+1):ncol(dataMM)] <- paste0("T", sort(rep(seq(0,(tpoints-1),1))))


    ### Make long data ####

    Msg <- "creating data. May take a while.\n"
    message(Msg)

    #### latent data
    datawide <- invisible(
      suppressMessages(
        suppressWarnings(
          ctIntervalise(data, Tpoints=tpoints, manifestNames=latentNames, n.manifest=n.latent)
        )
      )
    )
    datalong <- invisible(
      suppressMessages(
        suppressWarnings(
          ctWideToLong(datawide, n.manifest=n.latent, Tpoints=tpoints, manifestNames=latentNames )
        )
      )
    )
    datalong <- invisible(
      suppressMessages(
        suppressWarnings(as.data.frame(ctsem::ctDeintervalise(datalong))
        )
      )
    )
    datalong <- datalong[datalong$time >= burnin,]
    datalong$time <- datalong$time-(burnin)

    datalong <- datalong[datalong$time %in% tpointTargets[[i]], ]

    #### manifest data
    datawideMM <- invisible(
      suppressMessages(
        suppressWarnings(
          ctIntervalise(dataMM, Tpoints=tpoints, manifestNames=manifestNames, n.manifest=n.manifest)
        )
      )
    )
    datalongMM <- invisible(
      suppressMessages(
        suppressWarnings(
          ctWideToLong(datawideMM, n.manifest=n.manifest, Tpoints=tpoints, manifestNames=manifestNames)
        )
      )
    )
    datalongMM <- invisible(
      suppressMessages(
        suppressWarnings(as.data.frame(ctsem::ctDeintervalise(datalongMM))
        )
      )
    )
    datalongMM <- datalongMM[datalongMM$time >= burnin,]
    datalongMM$time <- datalongMM$time-(burnin)

    # tpointTargets
    datalong <- datalong[datalong$time %in% tpointTargets[[i]], ]
    datalongMM <- datalongMM[datalongMM$time %in% tpointTargets[[i]], ]
    #head(datalong)


    #studies[[i]] <- list() # required for normal loop (not dopar)
    #studies[[i]]$data <- datalong
    #studies[[i]]$tpointTargets <- tpointTargets[[i]]
    #studies[[i]]$drift <- drift[[i]]
    #studies[[i]]$drift_dt <- drift_dt[[i]]
    #studies[[i]]$diff <- diff[[i]]
    #studies[[i]]$diff_dt <- diff_dt[[i]]
    #studies[[i]]$cint <- cint[[i]]
    #studies[[i]]$cint_dt <- cint_dt[[i]]
    #studies[[i]]$T0means <- T0means[[i]]
    #studies[[i]]$T0var <- T0var[[i]]
    #studies[[i]]$randomIntercepts <- randomIntercepts[[i]]
    #studies[[i]]$manifestVars <- manifestVars[[i]]
    #studies[[i]]$manifestMeans <- manifestMeans[[i]]
    #studies[[i]]$TIpreds <- TIpreds[[i]]
    #head(studies[[i]]$data)

    tmp <- (list(data = datalong, # this (last) computation is automatically returned by doPar
                 dataMM = datalongMM,
                 tpointTargets = tpointTargets[[i]],
                 drift = drift[[i]],
                 drift_dt = as.matrix(drift_dt[[i]]),
                 diff = diff[[i]],
                 diff_dt = diff_dt[[i]],
                 cint = cint[[i]],
                 cint_dt = as.matrix(cint_dt[[i]]),
                 T0means = T0means[[i]],
                 T0var = T0var[[i]],
                 randomIntercepts = randomIntercepts[[i]],
                 manifestVars = manifestVars[[i]],
                 manifestMeans = manifestMeans[[i]],
                 TIpreds = TIpreds[[i]],
                 latentNames = latentNames,
                 modValues = modValues[[i]]))
  }

  # generate missings
  if (any(unlist(missings) != 0)) {
    for (i in 1:length(drift)) {
      # latents
      missingData <- studies[[i]]$data
      k <- round(missings[[i]] * nrow(missingData)); k
      for (c in latentNames) {
        random_values <- sort(sample(1:(sampleSizes[[i]] * length(tpointTargets[[i]])), k))
        missingData[random_values, c] <- NA
      }
      studies[[i]]$data <- missingData
      # manifests
      missingDataMM <- studies[[i]]$dataMM
      k <- round(missings[[i]] * nrow(missingDataMM)); k
      for (c in manifestNames) {
        random_values <- sort(sample(1:(sampleSizes[[i]] * length(tpointTargets[[i]])), k))
        missingDataMM[random_values, c] <- NA
      }
      studies[[i]]$dataMM <- missingDataMM
    }
  }

  # ctmaExtract ####
  if(ctmaExtract == TRUE) {
    if (useRawData == TRUE) tmp <- "rawData objects (rawDat1, rawData2, etc)"
    if (useRawData == FALSE) tmp <- "empcov objects (empcov1, empcov2, etc)"
    Msg <- paste0("\n\nctmaExtract was set to TRUE. Creating required CoTiMA objects incl. ", tmp, " in the environment specified in the argument envir.\n")
    message(Msg)

    #str(studies)
    ctmaExtract(activeDirectory = activeDirectory,
                ctmaGenDataList = studies,
                envir = envir,
                useRawData = useRawData)}

  # Return ####
  return(studies)
}
