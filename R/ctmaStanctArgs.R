#' This are preset arguments
#'
#' @name CoTiMAStanctArgs
#' @export CoTiMAStanctArgs
#'
CoTiMAStanctArgs<-list(test=TRUE,
                       scaleTI=FALSE, scaleMod=FALSE, scaleClus=TRUE,
                       scaleTD=FALSE, scaleLongData=FALSE, scaleTime=1/1,
                       scaleTimeAuto=NULL,
                       savesubjectmatrices=FALSE, verbose=1,
                       datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
                       iter=1000, intoverstates=TRUE,
                       binomial=FALSE, fit=TRUE,
                       intoverpop="auto", stationary=FALSE,
                       plot=FALSE, derrind="all",
                       optimize=TRUE, optimcontrol=list(is=FALSE, stochastic=FALSE, finishsamples=1000),
                       nlcontrol=list(),
                       nopriors=TRUE,
                       chains=2,
                       cores=1,
                       warmup=500,
                       inits=NULL, forcerecompile=FALSE,
                       savescores=FALSE, gendata=FALSE,
                       control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001))
