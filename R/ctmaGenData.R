#' ctmaGenData
#'
#' @description Generates data from lists of parameters (drift, diffusion etc). Experimental!!
#'
#' @param activeDirectory defines active directory where files are saved. No default.
#' @param burnin vector of inititial time points to be deleted (default = 0)
#' @param cint list of cint matrices. By default (NULL), cint matrices will be used that create a steady-state (i.e., means at all time points = T0means)
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use.
#' @param diff list of diffusion matrices. By default (NULL), diffusion matrices will be used that create a steady-state (i.e., covariance at all time points = T0var)
#' @param digits number of digits used for rounding (in outputs).
#' @param doPar parallel generating of data. A value > 1 will generate data for doPar studies (default = 1) in parallel mode during which no output is generated (screen remains silent).
#' @param drift list of drift matrices. No default (all = NULL).
#' @param empirical whether (default) or not generated data should allow for exact estimation of parameters (no random variance; not useful for MC simulations). This applies only if all Tpoints are used (see argument missings)
#' @param sampleSizes vector of sample sizes. Default = 100.
#' @param lambda list of matrices. By default all are diagonal matrices with 1 in the diagonal.
#' @param latentNames names for latent variables (default = NULL using generic names)
#' @param manifestMeans list of manifest mean matrices. By default all are = 0.
#' @param manifestVars list of manifest error (co-)variances matrices. By default all are = 0.
#' @param missings proportion of missings (default = 0, which does not delete any value)
#' @param n.latent number of latent variables of the model No default (all = NULL).
#' @param n.manifest number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).
#' @param randomIntercepts list of (co-)variances matrices of TIpreds (traits). By default all are = 0.
#' @param T0means list of Time 0 mean levels. By default all are = 0.
#' @param T0var list of Time 0 (co-)variance matrices. No default (all = NULL).
#' @param TIpreds list of time-independent predictors, e.g., the moderators that were used to create the drift matrices. No default (all = NULL).
#' @param tpoints vector of number of tpoints to be generated (default = 10).
#' @param tpointTargets list of vectors of tpoints to be selectd (default = burnin:tpoints).

#' @importFrom parallel detectCores makeCluster
#' @importFrom ctsem ctDeintervalise ctLongToWide ctIntervalise ctWideToLong ctModel ctStanFit ctExtract ctCollapse
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom OpenMx expm
#' @importFrom Matrix bdiag
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
    cint  = NULL,
    coresToUse = 2,
    diff = NULL,
    digits = 4,
    doPar = 1,
    drift = NULL,
    empirical = TRUE,
    sampleSizes = 100,
    lambda = NULL,
    latentNames = NULL,
    manifestMeans = 0,
    manifestVars = 0,
    missings = 0,
    n.latent = NULL,
    n.manifest = 0,
    randomIntercepts = 0,
    T0means = 0,
    T0var = NULL,
    TIpreds = NULL,
    tpoints = 10,
    tpointTargets = NULL
)
{

  original.options <- getOption("scipen"); original.options
  options(scipen = 999); options("scipen") # turn scientific notation off.
  on.exit(options(scipen = original.options))  # scientific notation as user's original

  # Check if not yet working arguments are set #####
  {
    if (manifestVars != 0) {
      ErrorMsg <- "\n Argument manifestVars not yet possible"
      stop(ErrorMsg)
    }

    if (n.manifest != 0) {
      ErrorMsg <- "\n Argument n.manifest not yet possible"
      stop(ErrorMsg)
    }
    if (!(is.null(TIpreds))) {
      ErrorMsg <- "\n Argument TIpreds not yet possible"
      stop(ErrorMsg)
    }
    if (!(is.null(lambda))) {
      ErrorMsg <- "\n Argument lambda (i.e, measurement models) not yet possible"
      stop(ErrorMsg)
    }
    #if (empirical == FALSE) {
    #  ErrorMsg <- "\n Argument empirical = FALSE not yet possible"
    #  stop(ErrorMsg)
    #}
  }

  # Check if arguments are properly used #####
  { # start checks
    if (utils::packageDescription("ctsem")$Version > "3.10.2") type <- "ct" else type <- "stanct"

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

    if (is.null(tpointTargets)) {
      tpointTargets <- rep(list(burnin:tpoints), length(drift))
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

    if (is.null(T0var)) {
      ErrorMsg <- "\n No T0var matrix (drift) specified! \nGood luck for the next try!"
      stop(ErrorMsg)
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

    if (is.null(n.latent)) {
      Msg <- "n.latent not specified. I infer n.latent from the dimensions of the drift matrix.\n"
      message(Msg)
      n.latent <- ncol(drift[[1]]); n.latent
    }

    ## trait means ####
    TRAITMEANS <- matrix(0, nrow=n.latent, ncol=1)

    if (n.manifest == 0) {
      Msg <- "n.manifest not specified. I assume n.manifest is equal to n.latent.\n"
      message(Msg)
      n.manifest <- n.latent
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

    if (length(manifestMeans[[1]]) != 1) {
      if (!(is.list(manifestMeans))) {
        ErrorMsg <- "\nThe manifestMeans argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(manifestMeans, function(x) nrow(x))))) != 1) {
        ErrorMsg <- "\nThe manifestMeans matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(manifestMeans, function(x) nrow(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe manifestMeans matrices supplied have different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    } else {
      if (manifestMeans[[1]] == 0) {
        manifestMeans <- rep(list(rep(0, n.latent)), length(drift))
      } else {
        print("Problem with manifestMeans")
      }
    }

    #if (manifestMeans == 0) {manifestMeans <- rep(list(matrix(rep(0, n.latent), nrow=n.latent, ncol=1)), length(drift))}

    if (manifestVars != 0) {
      if (!(is.list(manifestVars))) {
        ErrorMsg <- "\nThe manifestVars argument has to be a list. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (length(unique(unlist(lapply(manifestVars, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe manifestVars matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(manifestVars, function(x) ncol(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
        ErrorMsg <- "\nThe manifestVars matrices supplied have different dimensions than the drift matrices. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    if (manifestVars == 0) {manifestVars <- rep(list(matrix(0, n.latent, n.latent)), length(drift))}

    if (is.null(lambda)) {
      Msg <- "lambda not specified. I assume all lambdas are diagonal matrices with 1 in the diagonals.\n"
      message(Msg)
      lambda <- diag(1, n.latent, n.latent)
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
      if (length(unique(unlist(lapply(cint, function(x) ncol(x))))) != 1) {
        ErrorMsg <- "\nThe cint matrices supplied have different dimensions. \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if ( unique(unlist(lapply(cint, function(x) ncol(x)))) != unique(unlist(lapply(drift, function(x) ncol(x))))) {
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

    if (doPar > 1) {
      myCluster <- parallel::makeCluster(coresToUse)
      on.exit(parallel::stopCluster(myCluster))
      doParallel::registerDoParallel(myCluster)
      '%dopar%' <- foreach::'%dopar%'
    }

    if (n.manifest > n.latent ) {
      if (is.null(lambda)) {
        ErrorMsg <- "\nManifest variables specified, but not matrix of loadings (lambda) specified, which is required. \nGood luck for the next try!"
        stop(ErrorMsg)
      } else {
        if ( (dim(lambda)[1] != n.manifest) | (dim(lambda)[2] != n.latent) ) {
          ErrorMsg <- "\nDimensions of loadings matrix (lambda) do not match n.latent and n.manifest. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
    }

    if (length(sampleSizes) == 1) {
      if (length(sampleSizes[[1]]) == 1) {
        #if (sampleSizes == 100) sampleSizes <- as.list(rep(sampleSizes, length(drift)))
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

    ## Compute diffusions to achieve steady state  ####
    if (is.null(diff)) {
      diff <- list()
      for (i in 1:length(drift)) {
        #T1cov_impl <- expm(drift[[i]]) %*% T0var[[i]] %*% t(expm(drift[[i]])); T1cov_impl
        T1cov_impl <- expm(drift[[i]]) %*% (T0var[[i]] + randomIntercepts[[i]]) %*% t(expm(drift[[i]])); T1cov_impl
        #resvar <- T0var[[i]] - T1cov_impl; resvar
        resvar <- (T0var[[i]] + randomIntercepts[[i]]) - T1cov_impl; resvar
        if (length(unique(round(abs(c(resvar)), 5))) == 1) {
          ErrorMsg <- paste0("\nCannot generate data because of singularity issues with Study ", i, ",",
                             "\nwhich may be due all drift elements having identical magnitudes (e.g., all -.1 or .1).",
                             "\nYour provided drift = ", paste0(c(drift[[i]]), collapse = " "),
                             "\nGood luck for the next try!")
          stop(ErrorMsg)
        }
        DIAG <- diag(1, nrow(drift[[i]]), ncol(drift[[i]])); DIAG
        DRIFT_hatch <- drift[[i]] %x% DIAG + DIAG %x% drift[[i]]; DRIFT_hatch
        DIAG_hatch <- diag(1, nrow(DRIFT_hatch), ncol(DRIFT_hatch)); DIAG_hatch
        #solve(expm(DRIFT_hatch * 1) - DIAG_hatch)
        Q <- solve((expm(DRIFT_hatch * 1) - DIAG_hatch)) %*% DRIFT_hatch %*% c(resvar); Q
        diff[[i]] <- matrix(Q, n.latent, n.latent); diff[[i]]

        if (any(eigen(diff[[i]])$values < 0)) {
          ErrorMsg <- paste0("\nCannot generate data because of negative eigenvalues in the diffusion matrix of Study ", i, ",",
                             "\nwhich may be due too large cross or auto effects.",
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

    ## compute T0var if it was not provided as the asymDiffcov (does nto work because Q cannot be computed without T0var) ####
    # Solve A S + S t(A) + Q = 0 for S
    #if (missT0var == 1) { # set to 1 in the beginning if missing
    #  T0var <- list()
    #  for (i in 1:length(drift)) {
    #
    #    solve_lyapunov <- function(A, Q) {
    #      n <- nrow(A)
    #      K <- kronecker(diag(n), A) + kronecker(A, diag(n))  # I⊗A + A⊗I
    #      s <- solve(K, -as.vector(Q))
    #      matrix(s, n, n)
    #    }
    #    T0var[[i]] <- solve_lyapunov(A, Q)
    #    S
    #  }
    #}

    ## Compute cints to achieve steady state #####
    if (is.null(cint)) {
      cint <- list()
      for (i in 1:length(drift)) {
        #asymCINT <- T0means[[i]]; asymCINT
        asymCINT <- (T0means[[i]] + TRAITMEANS); asymCINT
        cint[[i]] <- -drift[[i]] %*% asymCINT; cint[[i]]
      }
    }

    ## Define DT matrices and check for standardization #####
    drift_dt <- diff_dt <- cint_dt <- list()
    for (i in 1:length(drift))  {
      drift_dt[[i]] <- expm(drift[[i]]); drift_dt[[i]]

      DIAG <- diag(1, nrow(drift[[i]]), ncol(drift[[i]])); DIAG
      DRIFT_hatch <- drift[[i]] %x% DIAG + DIAG %x% drift[[i]]; DRIFT_hatch
      diff_dt[[i]] <- matrix((solve(DRIFT_hatch) %*% (expm(DRIFT_hatch * 1) - diag(n.latent^2)) %*% c(diff[[i]])), n.latent, n.latent); diff_dt[[i]]

      cint_dt[[i]] <- solve(drift[[i]]) %*% (expm(drift[[i]]) - diag(1, n.latent, n.latent)) %*% cint[[i]]; cint_dt[[i]]

      T1cov_impl <- expm(drift[[i]]) %*% T0var[[i]] %*% t(expm(drift[[i]])); T1cov_impl
      resvar <- T0var[[i]] - T1cov_impl; resvar
      if (!(all(round(diag(T1cov_impl + resvar), 2) == 1))) {
        Msg <- paste0("\nThe steady state-variances rounded to 2 digits of Study ", i, " are not 1.0,",
                      "\nwhich implies that the data are unstandardized and not useful for CoTiMA.",
                      "\nThey are = ",
                      paste0(diag(T1cov_impl + resvar), collapse = " "),
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

  } # end checks

  # Generate data ####
  studies <- list()
  for (i in 1:length(drift)) {
    #i <- 1

    ## T0 data, diffusion, and traitvar all independent ####
    ### diffusions #####
    diff.var_tmp <- list(diff_dt[[i]]); diff.var_tmp
    diff.var <- rep(diff.var_tmp, tpoints); diff.var
    diff.var <- as.matrix(Matrix::bdiag(diff.var))
    rows1 <- nrow(diff.var)
    #round(diff.var, 2)
    ### T0var #####
    diffT0.var <- cbind(diff.var, matrix(0, ncol=n.manifest, nrow=rows1)); diffT0.var # diffvar & T0var
    cols1 <- ncol(diffT0.var); cols1
    diffT0.var <- rbind(diffT0.var, matrix(0, ncol=cols1, nrow=n.manifest))
    diffT0.var[(cols1-n.manifest+1):cols1, (cols1-n.manifest+1):cols1] <- T0var[[i]]
    rows2 <- nrow(diffT0.var); rows2
    #round(diffT0.var, 2)
    ### Traitvar ####
    diffT0Trait.var <- cbind(diffT0.var, matrix(0, ncol=n.manifest, nrow=rows2)) # diffvar & T0var & traitvar
    cols2 <- ncol(diffT0Trait.var); cols2
    diffT0Trait.var <- rbind(diffT0Trait.var, matrix(0, ncol=cols2, nrow=n.manifest))
    diffT0Trait.var[(cols2-n.manifest+1):cols2, (cols2-n.manifest+1):cols2] <- randomIntercepts[[i]]
    rows3 <- nrow(diffT0Trait.var); rows3
    #round(diffT0Trait.var, 2)
    ### manifestVar (measurement error) ####
    err.var_tmp <- list(manifestVars[[i]]); err.var_tmp
    err.var <- rep(err.var_tmp, tpoints); err.var
    err.var <- as.matrix(Matrix::bdiag(err.var))
    a_rows <- nrow(diffT0Trait.var)
    a_cols <- ncol(diffT0Trait.var)
    b_rows <- nrow(err.var)
    b_cols <- ncol(err.var)
    diffT0TraitERR.var <- matrix(0, nrow = a_rows + b_rows, ncol = a_cols + b_cols)
    diffT0TraitERR.var[1:a_rows, 1:a_cols] <- diffT0Trait.var
    diffT0TraitERR.var[(a_rows+1):(a_rows+b_rows), (a_cols+1):(a_cols+b_cols)] <- err.var
    cols3 <- ncol(diffT0TraitERR.var); cols3
    #dim(diffT0TraitERR.var)
    #round(diffT0TraitERR.var, 2)
    ## mvrnorm to create independent data for diffusions, T0 variables, random intercepts (traits), and measurement error
    #allInit.dat <- MASS::mvrnorm(n=sampleSizes[[i]], mu=c(rep(manifestMeans[[i]], n.manifest*tpoints), T0means[[i]], TRAITMEANS), Sigma = diffT0TraitERR.var, empirical=empirical)
    tmp <- length(c(rep(0, n.latent*tpoints), T0means[[i]], TRAITMEANS, rep(0, n.latent*tpoints)))
    if ( tmp > sampleSizes[[i]]) {
      ErrorMsg <- paste0("\n Requested sample size too small. It should at least: sampleSizes=", tmp, ".")
      stop(ErrorMsg)
    }
    allInit.dat <- MASS::mvrnorm(n=sampleSizes[[i]], mu=c(rep(0, n.latent*tpoints), T0means[[i]], TRAITMEANS, rep(0, n.latent*tpoints)), Sigma = diffT0TraitERR.var, empirical=empirical)
    #round(cov(allInit.dat)[1:22, 1:22], 2)
    #round(cov(allInit.dat)[23:44, 23:44], 2)
    diff.dat <- allInit.dat[, (1:(n.manifest*tpoints))]; dim(diff.dat)
    #round(cov(diff.dat), 2)
    T0.dat <- allInit.dat[, (cols1 -n.manifest+1):(cols1)]; dim(T0.dat)
    #round(cov(T0.dat), 2)
    trait.dat <- allInit.dat[, (cols2 -n.manifest+1):(cols2)]; dim(trait.dat)
    #round(cov(trait.dat), 2)
    #round(cov(cbind(T0.dat, trait.dat)), 2)
    #apply(cbind(T0.dat, trait.dat), 2, mean)
    err.dat <- allInit.dat[, (a_cols+1):(cols3)]; dim(err.dat)
    #round(cov(err.dat), 2)
    # data that do not vary among cases and tpoints
    cint.dat <- matrix(t(cint_dt[[i]]),  nrow=sampleSizes[[i]], ncol=n.latent, byrow = T)
    #head(cint.dat)
    manifestMeans.dat <- matrix(t(manifestMeans[[i]]),  nrow=sampleSizes[[i]], ncol=n.latent, byrow = T)
    #head(manifestMeans.dat)
    #
    data <- T0.dat
    #data <- data + trait.dat
    #head(data)
    #data <- data + matrix(t(T0means[[i]]), ncol=n.latent, nrow=nrow(data), byrow=T)
    #head(matrix(t(T0means[[i]]), ncol=n.latent, nrow=nrow(data), byrow=T))
    #head(data)

    #### T1, T2, ... all subsequent Tpoints ####
    for (t in 1:(tpoints-1)) {
      tmpData <- data[, ((t-1)*n.latent+1):((t-1)*n.latent+n.latent)]
      tmp <- t(apply(tmpData, 1, function(x) drift_dt[[i]] %*% x))
      # add diffusion
      if (empirical == TRUE) {
        tmp <- tmp + diff.dat[, (2*(t-1)+1):(2*(t-1)+n.manifest)]
      } else {
        tmp <- tmp + MASS::mvrnorm(n=sampleSizes[[i]], mu=rep(0, n.latent), Sigma = diff_dt[[i]], empirical=FALSE)
        #print("Problem with empirical argument!!")
      }
      # add cint
      tmp <- tmp + cint.dat
      # add trait (done later)
      # tmp <- tmp + trait.dat
      # combine
      data <- cbind(data, tmp)
    }
    #head(data)

    #### Add constant trait.dat ####
    #head(trait.dat)
    for (t in seq(1, ncol(data), n.latent)) data[, c(t, t+1)] <- data[, c(t, t+1)] + trait.dat
    #for (t in seq(n.latent+1, ncol(data), n.latent)) data[, c(t, t+1)] <- data[, c(t, t+1)] + trait.dat
    #head(data)

    #### Add constant manifestMeans.dat ####
    for (t in seq(1, ncol(data), n.latent)) data[, c(t, t+1)] <- data[, c(t, t+1)] + manifestMeans.dat
    #head(data)

    # Make measurement model
    # to be done (and tbd further below not here)
    # add manifestVar
    # to be done  (and tbd further below not here)

    colnames(data) <- paste0(paste0(latentNames, "_T"), sort(rep(seq(0,(tpoints-1),1), 2)))
    data <- cbind(data, matrix(seq(0, tpoints-1, 1), nrow=nrow(data), ncol=tpoints, byrow=T))
    colnames(data)[(ncol(data)-tpoints+1):ncol(data)] <- paste0("T", sort(rep(seq(0,(tpoints-1),1))))
    #head(data); dim(data)

    ### Make long data ####
    datawide <- invisible(
      suppressMessages(
        suppressWarnings(
          ctIntervalise(data, Tpoints=tpoints, manifestNames=latentNames, n.manifest=n.manifest)
        )
      )
    )
    datalong <- invisible(
      suppressMessages(
        suppressWarnings(
          ctWideToLong(datawide, n.manifest=n.manifest, Tpoints=tpoints, manifestNames=latentNames )
        )
      )
    )
    datalong <- as.data.frame(ctsem::ctDeintervalise(datalong))
    datalong <- datalong[datalong$time >= burnin,]
    datalong$time <- datalong$time-(burnin)
    #head(datalong, 10)

    # tpointTargets
    datalong <- datalong[datalong$time %in% tpointTargets[[i]], ]

    studies[[i]] <- list()
    studies[[i]]$data <- datalong
    studies[[i]]$tpointTargets <- tpointTargets[[i]]
    studies[[i]]$drift <- drift[[i]]
    studies[[i]]$drift_dt <- drift_dt[[i]]
    studies[[i]]$diff <- diff[[i]]
    studies[[i]]$diff_dt <- diff_dt[[i]]
    studies[[i]]$cint <- cint[[i]]
    studies[[i]]$cint_dt <- cint_dt[[i]]
    studies[[i]]$T0means <- T0means[[i]]
    studies[[i]]$T0var <- T0var[[i]]
    studies[[i]]$randomIntercepts <- randomIntercepts[[i]]
    studies[[i]]$manifestVars <- manifestVars[[i]]
    studies[[i]]$manifestMeans <- manifestMeans[[i]]
    studies[[i]]$TIpreds <- TIpreds[[i]]
    #head(studies[[i]]$data)
  }

  #dat1 <- studyData[[1]]
  #cov(dat1[dat1$time == 0 ,])
  #cov(dat1[dat1$time == 9 ,])
  #apply(dat1[dat1$time == 0 ,], 2, mean)
  #apply(dat1[dat1$time == 9 ,], 2, mean)

  #dat2 <- studyData[[2]]
  #cov(dat2[dat2$time == 0 ,])
  #cov(dat2[dat2$time == 9 ,])
  #apply(dat2[dat2$time == 0 ,], 2, mean)
  #apply(dat2[dat2$time == 9 ,], 2, mean)

  #saveRDS(datalong, file=paste0(activeDirectory, "Study_", i, "_data_all.rds))

  # generate missings
  if (missings != 0) {
    for (i in 1:length(drift)) {
      missingData <- studies[[i]]$data
      k <- round(missings * nrow(missingData)); k
      for (c in latentNames) {
        random_values <- sort(sample(1:(sampleSizes[[i]] * (tpoints-burnin)), k))
        missingData[random_values, c] <- NA
      }
      missingData$miss <- apply(missingData[,latentNames], 1, sum, na.rm=T)
      missingData <- missingData[missingData$miss !=0,]
      missingData$miss <- NULL
      studies[[i]]$data <- missingData
    }
  }

  return(studies)
}

