#' ctmaPower
#'
#' @description Fits or loads an invariant continuous-time model, then computes plug-in Gaussian discrete-time regression power and required sample sizes across positive time lags. Verified unrestricted regressions use covariance calculations; supported restrictions use lavaan regression fits.
#'
#' @param ctmaInitFit object to which all single 'ctsem' fits of primary studies has been assigned to (i.e., what has been returned by \code{\link{ctmaInit}})
#' @param activeDirectory defines another active directory than the one used in \code{\link{ctmaInit}}
#' @param statisticalPower Distinct requested power probabilities strictly between zero and one; if empty, uses c(.9, .8).
#' @param failSafeN  sample size used to determine across which time intervals effects become non-significant
#' @param failSafeP  p-value used to determine across which time intervals effects become non-significant
#' @param timeRange Nonnegative start/end/positive-step triple. Only positive lags within the bounds are analysed; observed study lags inside those bounds are included.
#' @param useMBESS Use MBESS's legacy upper-tail probability with a two-sided critical value. The default uses the complete two-sided probability.
#' @param coresToUse if negative, the value is subtracted from available cores, else value = cores to use
#' @param indVarying Allows continuous time intercepts to vary at the individual level (random effects model, accounts for unobserved heterogeneity)
#' @param activateRPB set to TRUE to receive push messages with 'CoTiMA' notifications on your phone
#' @param silentOverwrite overwrite old files without asking
#' @param digits number of digits used for rounding (in outputs)
#' @param loadAllInvFit load the fit of fully constrained 'CoTiMA' model
#' @param saveAllInvFit save the fit of fully constrained 'CoTiMA' model
#' @param loadAllInvWOSingFit Load reduced lavaan fits with a filename prefix. Forces the lavaan route; cache covariance/model metadata must match. Legacy bare-fit caches must be regenerated.
#' @param saveAllInvWOSingFit Save reduced lavaan fits with a filename prefix. Forces the lavaan route; cache covariance/model metadata must match. Legacy bare-fit caches must be regenerated.
#' @param skipScaling Must be TRUE; rescaling must be performed explicitly before invariant fitting.
#' @param useSampleFraction Must be NULL; subset the input studies explicitly instead.
#' @param optimize TRUE requests optimization of the invariant CT model; FALSE requests posterior sampling with the selected backend. Power calculations use plug-in matrices in either case.
#' @param priors if FALSE, any priors are disabled – sometimes desirable for optimization
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param iter number of iterations (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#' @param customPar logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting

#'
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red blue
#' @importFrom parallel detectCores
#' @importFrom ctsem ctModel ctWideToLong ctDeintervalise
#' @importFrom OpenMx vech2full expm
#' @importFrom stats cov2cor pt qt na.omit median
#' @importFrom lavaan sem inspect
#' @importFrom MBESS ss.power.reg.coef
#' @importFrom rootSolve uniroot.all
#'
#' @examples \dontrun{
#' CoTiMAInitFit_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAPower_D_BO <- ctmaPower(ctmaInitFit=CoTiMAInitFit_D_BO,
#'                               statisticalPower = c(.50, .80, .95),
#'                               finishsamples = 10000)
#' summary(CoTiMAPower_D_BO)
#' }
#'
#' @param powerMethod Regression calculation route. \code{"auto"} uses direct
#'   covariance calculations only for a verified unrestricted Gaussian regression;
#'   otherwise supported Gaussian restrictions use lavaan. \code{"covariance"}
#'   requires eligibility; \code{"lavaan"} explicitly requests regression fits.
#' @param regressionModel Optional lavaan syntax for the discrete-time Gaussian
#'   regression on variables \code{V1T0,...,VpT0,V1T1,...,VpT1}. NULL constructs
#'   a complete unrestricted regression. Fixed/omitted target slopes have NA power.
#'   Shared labels, equality constraints, latent measurement models and non-Gaussian
#'   indicators are outside the one-coefficient power approximation.
#' @details The fast default includes both tails of the noncentral-t test and
#'   finds the minimum integer N by bisection. \code{useMBESS=TRUE} retains MBESS's
#'   legacy upper-tail probability. Calculations concern plug-in discrete-time
#'   regression effects, not the power of the pooled CT drift test. Missing post-hoc
#'   powers remain NA. The legacy Opt. Lag row selects the first actual minimizing
#'   lag; all ties are returned separately. Console and summary notes identify
#'   whether covariance calculations or lavaan fits were used.
#'   The state projection requires directly observed variables, identity loadings
#'   and zero measurement error; recorded latent measurement/augmented-state
#'   models require simulation-based power analysis.
#'
#' @export ctmaPower
#'
#' @return A \code{"CoTiMAFit"} object containing the invariant model, matrix
#'   estimates, regression eligibility/details, calculation/interpretation notes,
#'   required sample sizes and post-hoc power tables. \code{power_method} records
#'   the route. The legacy \code{Opt. Lag} row selects the first tied minimum;
#'   \code{summary$estimates[["All optimal lags"]]} retains every tied lag.
#'   Infinite N denotes an unattainable target or a requirement beyond the
#'   exact-integer search range; unavailable powers remain NA.
#'
ctmaPower <- function(
    ctmaInitFit=NULL,
    activeDirectory=NULL,
    statisticalPower=c(),
    failSafeN =NULL,
    failSafeP=NULL,
    timeRange=NULL,
    useMBESS=FALSE,
    coresToUse=1,
    digits=4,
    indVarying=FALSE,
    activateRPB=FALSE,
    silentOverwrite=FALSE,
    loadAllInvFit=c(),
    saveAllInvFit=c(),
    loadAllInvWOSingFit=c(),
    saveAllInvWOSingFit=c(),
    skipScaling=TRUE,
    useSampleFraction=NULL,
    optimize=TRUE,
    priors=FALSE,
    finishsamples=NULL,
    iter=NULL,
    chains=NULL,
    verbose=NULL,
    customPar=FALSE,
    scaleTime=NULL,
    powerMethod=c("auto", "covariance", "lavaan"),
    regressionModel=NULL
) {
  backend <- "stan"; backendArgs <- list()
  stan_args <- CoTiMAStanctArgs
  if (is.null(scaleTime)) scaleTime <- 1
  stan_args$scaleTime <- scaleTime
  stan_args$optimize <- optimize; stan_args$priors <- priors
  if (!is.null(finishsamples)) stan_args$optimcontrol$finishsamples <- finishsamples
  if (!is.null(iter)) stan_args$iter <- iter
  if (!is.null(chains)) stan_args$chains <- chains
  if (!is.null(verbose)) stan_args$verbose <- verbose
  available <- parallel::detectCores()
  if (!is.finite(available)) available <- 1L
  if (length(coresToUse)!=1L || !is.finite(coresToUse)) stop("coresToUse must be finite.", call.=FALSE)
  coresToUse <- max(1, min(available, if (coresToUse<1) available+coresToUse else coresToUse))
  resources <- list(resolved=coresToUse)
  # 1. Validate the design and build the evaluated lag grid.
  design <- .ctmaPowerDesign(ctmaInitFit, timeRange, statisticalPower, failSafeN,
    failSafeP, activeDirectory, skipScaling, useSampleFraction, digits)
  powerMethod <- match.arg(powerMethod)
  if (!is.logical(useMBESS) || length(useMBESS)!=1L || is.na(useMBESS)) stop("useMBESS must be TRUE or FALSE.", call.=FALSE)
  if (!is.finite(scaleTime) || length(scaleTime)!=1L || scaleTime<=0) stop("scaleTime must be positive.", call.=FALSE)
  for (prefix in list(loadAllInvFit,saveAllInvFit,loadAllInvWOSingFit,saveAllInvWOSingFit))
    if (length(prefix)) .ctmaPowerPath(design$directory,prefix)
  # Reject known non-Gaussian indicators before a potentially expensive fit.
  .ctmaPowerGaussian(ctmaInitFit,list())

  # 2. Obtain the invariant fit. Data generation and backend fitting belong to
  # ctmaAllInvFit; ctmaPower no longer duplicates and discards that preparation.
  if (length(loadAllInvFit)) {
    .ctmaPowerMessage("Loading the cached invariant model.")
    pooled <- readRDS(.ctmaPowerPath(design$directory,loadAllInvFit))
  } else {
    .ctmaPowerMessage("Fitting the invariant model used for plug-in power calculations.")
    fit_args <- list(ctmaInitFit=ctmaInitFit,activeDirectory=file.path(design$directory,""),
      activateRPB=activateRPB,digits=digits,coresToUse=coresToUse,
      n.manifest=if (is.null(ctmaInitFit$n.manifest)) 0 else ctmaInitFit$n.manifest,
      indVarying=indVarying,scaleTime=scaleTime,optimize=optimize,priors=priors,
      finishsamples=finishsamples,iter=iter,chains=chains,verbose=stan_args$verbose,
      customPar=customPar)
    fit_args$CoTiMAStanctArgs <- stan_args
    pooled <- do.call(ctmaAllInvFit,fit_args)
  }
  if (length(saveAllInvFit)) .ctmaPowerSave(pooled,
    .ctmaPowerPath(design$directory,saveAllInvFit),silentOverwrite)

  # 3. Extract matrices, choose the verified regression route, evaluate each
  # lag once, and calculate probabilities/sample sizes for all requested powers.
  analysis <- .ctmaPowerAnalyse(ctmaInitFit,pooled,design,useMBESS,powerMethod,
    regressionModel,loadAllInvWOSingFit,saveAllInvWOSingFit,silentOverwrite,digits)
  controls <- list(backend=backend,backendArgs=backendArgs,resource_resolution=resources,
    cores=coresToUse,coresToUse=coresToUse,CoTiMAStanctArgs=if (backend=="stan") stan_args else NULL,
    scaleTime=scaleTime,optimize=optimize,priors=priors,
    powerMethod=powerMethod,regressionModel=regressionModel,useMBESS=useMBESS,
    finishsamples=finishsamples,iter=iter,chains=chains)
  # 4. Assemble the established plotting/summary tables and structured diagnostics.
  invisible(.ctmaPowerResult(ctmaInitFit,pooled,design,analysis,controls,digits))
}
