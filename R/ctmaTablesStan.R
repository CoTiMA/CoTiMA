# Regular ctsem stores Stan arrays differently from the native Julia backend.
# Use the fitted Stan program for raw-to-model transformations; never guess
# transforms from parameter names or treat rawest as a distribution of draws.
.ctmaTablesStanSampled <- function(fit) {
  sf<-fit$stanfit$stanfit
  isS4(sf)&&'sim' %in% methods::slotNames(sf)&&length(methods::slot(sf,'sim'))>0L
}
.ctmaTablesJointArray <- function(x) {
  is.array(x)&&length(dim(x))==3L&&is.numeric(x)&&dim(x)[1L]>=2L&&
    dim(x)[2L]>=1L&&dim(x)[2L]==dim(x)[3L]&&all(is.finite(x))
}
.ctmaTablesStanDriftDraws <- function(fit) {
  if(!inherits(fit,'ctStanFit'))return(NULL)
  a<-fit$stanfit[["transformedpars"]]$pop_DRIFT
  if(.ctmaTablesJointArray(a))return(a)
  if(.ctmaTablesStanSampled(fit)) {
    sf<-fit$stanfit$stanfit;sim<-methods::slot(sf,'sim')
    if(!is.null(sim$pars_oi)&&!'pop_DRIFT' %in% sim$pars_oi)return(NULL)
    if(!requireNamespace('rstan',quietly=TRUE)) .ctmaTablesStop('Reading saved Stan draws requires package rstan.')
    a<-tryCatch(rstan::extract(sf,pars='pop_DRIFT',permuted=TRUE,inc_warmup=FALSE)$pop_DRIFT,
      error=function(e)NULL)
    if(.ctmaTablesJointArray(a))return(a)
  }
  NULL
}
.ctmaTablesStanReady <- function(fit,label) {
  if(!inherits(fit,'ctStanFit')||.ctmaTablesStanSampled(fit))return(fit)
  cached<-fit$stanfit[["transformedpars"]]$pop_DRIFT
  if(.ctmaTablesJointArray(cached))return(fit)
  raw<-fit$stanfit[["rawposterior"]]
  if(!is.matrix(raw)||!is.numeric(raw)||nrow(raw)<2L)
    .ctmaTablesStop('The saved Stan fit for ',label,' lacks joint uncertainty draws. rawest is a point estimate. Supply driftDraws or use dtPointOnly=TRUE for DT point estimates.')
  if(any(!is.finite(raw))) .ctmaTablesStop('The retained raw Stan draws for ',label,' contain nonfinite values.')
  if(is.null(fit$stanmodel)||!is.list(fit$standata))
    .ctmaTablesStop('Transforming raw Stan draws for ',label,' requires the saved stanmodel and standata.')
  if(!requireNamespace('ctsem',quietly=TRUE)) .ctmaTablesStop('Transforming saved Stan draws requires ctsem.')
  ns<-asNamespace('ctsem')
  if(!exists('stan_constrainsamples',envir=ns,inherits=FALSE))
    .ctmaTablesStop('This ctsem version cannot reconstruct the missing Stan transformation cache.')
  transform<-get('stan_constrainsamples',envir=ns,inherits=FALSE)
  message('CoTiMA tables: transforming ',nrow(raw),' saved joint raw Stan draws for ',label,'; no model refitting.')
  a<-tryCatch(transform(sm=fit$stanmodel,standata=fit$standata,samples=raw,cores=1L,
      savescores=FALSE,savesubjectmatrices=FALSE,dokalman=FALSE,quiet=TRUE),
    error=function(e).ctmaTablesStop('Could not transform saved Stan draws for ',label,': ',conditionMessage(e)))
  if(!.ctmaTablesJointArray(a$pop_DRIFT)||dim(a$pop_DRIFT)[1L]!=nrow(raw))
    .ctmaTablesStop('Stan transformation did not preserve all joint DRIFT draws for ',label,'. No draws were silently discarded.')
  fit$stanfit$transformedpars<-a
  fit
}
