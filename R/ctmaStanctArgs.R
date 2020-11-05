#' This are preset arguments
#'
#' @name CoTiMAStanctArgs
#' @docType data
#' @author My Name \email{blahblah@@roxygen.org}
#' @references \url{data_blah.com}
#' @keywords data
"CoTiMAStanctArgs"

CoTiMAStanctArgs<-list(test=TRUE,
                       scaleTI=TRUE, scaleMod=TRUE, scaleClus=FALSE, scaleLongData=FALSE,
                       scaleTime=1/1,
                       savesubjectmatrices=FALSE, verbose=1,
                       datalong=NA, ctstanmodel=NA, stanmodeltext = NA,
                       iter=1000, intoverstates=TRUE,
                       binomial=FALSE, fit=TRUE,
                       intoverpop='auto', stationary=FALSE,
                       plot=FALSE, derrind="all",
                       optimize=TRUE, optimcontrol=list(is=F, stochastic=FALSE),
                       nlcontrol=list(),
                       nopriors=TRUE,
                       chains=2,
                       cores=1,
                       inits=NULL, forcerecompile=FALSE,
                       savescores=FALSE, gendata=FALSE,
                       control=list(adapt_delta = .8, adapt_window=2, max_treedepth=10, adapt_init_buffer=2, stepsize = .001),
                       verbose=0,
                       scaleTD=FALSE)
