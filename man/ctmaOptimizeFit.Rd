% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctmaOptimizeFit.R
\name{ctmaOptimizeFit}
\alias{ctmaOptimizeFit}
\title{ctmaOptimizeFit}
\usage{
ctmaOptimizeFit(
  activateRPB = FALSE,
  activeDirectory = NULL,
  checkSingleStudyResults = FALSE,
  coresToUse = c(2),
  CoTiMAStanctArgs = NULL,
  ctmaFitFit = NULL,
  ctmaInitFit = NULL,
  customPar = FALSE,
  finishsamples = NULL,
  indVarying = FALSE,
  lambda = NULL,
  manifestMeans = 0,
  manifestVars = NULL,
  n.latent = NULL,
  posLL = TRUE,
  primaryStudies = NULL,
  problemStudy = NULL,
  randomPar = FALSE,
  randomScaleTI = FALSE,
  randomScaleTime = c(1, 1),
  saveModelFits = FALSE,
  shuffleStudyList = FALSE,
  reFits = NULL,
  scaleMod = NULL,
  scaleTime = NULL,
  scaleTI = TRUE,
  T0means = 0,
  transfMod = NULL,
  parallel = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{activateRPB}{set to TRUE to receive push messages with 'CoTiMA' notifications on your phone}

\item{activeDirectory}{activeDirectory}

\item{checkSingleStudyResults}{displays estimates from single study 'ctsem' models and waits for user input to continue.}

\item{coresToUse}{if neg., the value is subtracted from available cores, else value = cores to use}

\item{CoTiMAStanctArgs}{parameters that can be set to improve model fitting of the \code{\link{ctStanFit}} Function}

\item{ctmaFitFit}{a object fitted with \code{\link{ctmaFit}}}

\item{ctmaInitFit}{the ctmaInitFit object that was used to create the ctmaFitFit object with \code{\link{ctmaFit}}}

\item{customPar}{logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)}

\item{finishsamples}{number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).}

\item{indVarying}{control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means}

\item{lambda}{R-type matrix with pattern of fixed (=1) or free (any string) loadings.}

\item{manifestMeans}{Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.}

\item{manifestVars}{define the error variances of the manifests within a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest. Useful to check estimates before they are saved.}

\item{n.latent}{number of latent variables of the model (hast to be specified)!}

\item{posLL}{logical. Allows (default = TRUE) of positive loglik (neg -2ll) values}

\item{primaryStudies}{list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}}

\item{problemStudy}{number (position in list) where the problem study in primaryStudies is found}

\item{randomPar}{logical (default = FALSE). Overrides arguments used for customPar and randomly sets customPar either TRUE or FALSE}

\item{randomScaleTI}{logical (default = FALSE). Overrides arguments used for scaleTI and randomly sets scaleTI either TRUE or FALSE}

\item{randomScaleTime}{lower and upper limit (default = c(1,1)) of uniform distribution from which timeScale argument for ctmaInit is uniformly shuffled (integer)}

\item{saveModelFits}{save the fit of each Fit attempt (default = FALSE).}

\item{shuffleStudyList}{(default = FALSE) randomly re-arranges studies in primaryStudyList. We encountered a few cases where this mattered, even though it should not. Only works if ctmaFit is optimized.}

\item{reFits}{how many reFits should be done}

\item{scaleMod}{scale moderator variables - TRUE (default) recommended for continuous and categorical moderators, to separate withing and betwen efeccts}

\item{scaleTime}{scale time (interval) - sometimes desirable to improve fitting}

\item{scaleTI}{scale TI predictors - not recommended until version 0.5.3.1. Does not change aggregated results anyways, just interpretation of effects for dummies representing primary studies.}

\item{T0means}{Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.}

\item{transfMod}{more general option to change moderator values. A vector as long as number of moderators analyzed (e.g., c("mean(x)", "x - median(x)"))}

\item{parallel}{(default = FALSE). When set to trUe parallel fitting on clusters is enabled (could save some time when many refits are done)}

\item{verbose}{integer from 0 to 2. Higher values print more information during model fit – for debugging}
}
\value{
returns a list with bestFit (= the best fit achieved), all_minus2ll (= all -2ll values for all fitted models), and summary, which
is printed if the summary function is applied to the returned object, and which shows the summary information of the ctsem model with the
best fit.

%dopar%' <- foreach::'%dopar%' deprecated
}
\description{
Replaces deprecated \code{\link{ctmaOptimizeInit}}, which was limited to initial fitting
(i.e., applies \code{\link{ctmaInit}}) of a primary study reFits times to capitalize on chance for obtaining
a hard-to-find optimal fit.
Now, optimizing a CoTiMA model generated with \code{\link{ctmaFit}} can also be done.
Using \code{\link{ctmaOptimizeFit}} could be helpful if a model yields out-of-range estimates, which could happen if the fitting
algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
\code{\link{ctmaOptimizeFit}} is like gambling, hoping that at least one set of starting values (the number it tries is specified in the reFits argument)
enables finding the global optimal fit. On unix-like machines (e.g. MacOS), this could be done in parallel mode if coresToUse > 1.
}
\note{
All but one of multiple cores are used on unix-type machines for parallel fitting

During fitting, not output is generated. Be patient.
}
\examples{
\dontrun{
optimFit313 <- ctmaOptimizeFit(primaryStudies=CoTiMAstudyList_3,
                                activeDirectory="/Users/tmp/",  # adapt!
                                problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
                                reFits=10,
                                n.latent=2)
summary(optimFit313)
}

}
