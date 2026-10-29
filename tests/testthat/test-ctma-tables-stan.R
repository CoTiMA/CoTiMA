library(testthat)
stan_tables_fixture <- function(mode='ML',cache=TRUE) {
 a<-array(seq_len(16)/100,c(4,2,2))
 structure(list(args=list(optimize=mode!='Bayes',priors=mode!='ML'),
   stanmodel='saved-model',standata=list(priors=as.integer(mode!='ML')),
   stanfit=list(rawest=1:3,rawposterior=matrix(seq_len(12),4,3),
     transformedpars=if(cache)list(pop_DRIFT=a) else NULL)),class='ctStanFit')
}
test_that('regular Stan direct metadata identifies ML MAP and Bayes',{
 for(mode in c('ML','MAP','Bayes')) {
   f<-stan_tables_fixture(mode)
   expect_identical(.ctmaTablesMethod(f),if(mode=='Bayes')'Bayesian' else mode)
   wrapped<-structure(list(studyFitList=f),class='CoTiMAFit')
   expect_identical(.ctmaTablesMethod(wrapped),if(mode=='Bayes')'Bayesian' else mode)
 }
})
test_that('saved Stan arrays are already joint transformed draws',{
 f<-stan_tables_fixture()
 expect_identical(.ctmaTablesDriftDraws(f),f$stanfit$transformedpars$pop_DRIFT)
 f$stanfit$transformedpars$pop_DRIFT<-array(1:4,c(1,2,2))
 expect_null(.ctmaTablesDriftDraws(f))
 f$stanfit$rawposterior<-NULL
 expect_error(.ctmaTablesStanReady(f,'primary 7'),'joint uncertainty draws')
})
test_that('missing Stan cache uses fitted program without Kalman filtering or refitting',{
 f<-stan_tables_fixture(cache=FALSE);original<-f;calls<-0L
 expected<-array(seq_len(16)/100,c(4,2,2))
 local_mocked_bindings(stan_constrainsamples=function(sm,standata,samples,cores,
     savescores,savesubjectmatrices,dokalman,quiet) {
   calls<<-calls+1L
   expect_identical(sm,'saved-model');expect_identical(samples,f$stanfit$rawposterior)
   expect_false(dokalman);expect_false(savescores);expect_false(savesubjectmatrices)
   expect_identical(cores,1L)
   list(pop_DRIFT=expected)
 },.package='ctsem')
 z<-.ctmaTablesStanReady(f,'primary 7')
 expect_identical(z$stanfit$transformedpars$pop_DRIFT,expected)
 expect_identical(f,original);expect_identical(calls,1L)
 z2<-.ctmaTablesStanReady(z,'primary 7')
 expect_identical(z2,z);expect_identical(calls,1L)
})
test_that('inadmissible or discarded Stan draws fail explicitly',{
 f<-stan_tables_fixture(cache=FALSE)
 local_mocked_bindings(stan_constrainsamples=function(...)list(pop_DRIFT=array(1,c(3,2,2))),.package='ctsem')
 expect_error(.ctmaTablesStanReady(f,'primary 7'),'preserve all joint DRIFT draws')
 f$stanfit$rawposterior[1,1]<-NA_real_
 expect_error(.ctmaTablesStanReady(f,'primary 7'),'nonfinite')
})
test_that('raw transformations are reused across requested summary lags',{
 f<-stan_tables_fixture(cache=FALSE);count<-0L;lags<-numeric()
 local_mocked_bindings(.ctmaTablesStanReady=function(fit,label) {count<<-count+1L;fit},
   .ctmaTablesRunSummary=function(target,lag,digits,label) {lags<<-c(lags,lag);list(parmatrices=data.frame())})
 .ctmaTablesGenerateSummaries(f,c(1,5),4,'primary 7')
 expect_identical(count,1L);expect_identical(lags,c(1,5))
})

test_that('legacy cache names are never accepted through partial matching',{
 f<-stan_tables_fixture();f$stanfit$transformedpars_old<-f$stanfit$transformedpars
 f$stanfit$transformedpars<-NULL
 expect_null(.ctmaTablesStanDriftDraws(f))
 f$stanfit$rawposterior_old<-f$stanfit$rawposterior;f$stanfit$rawposterior<-NULL
 expect_error(.ctmaTablesStanReady(f,'primary 7'),'joint uncertainty draws')
})
test_that('fresh native summaries retain precision for the four-decimal display',{
 f<-stan_tables_fixture()
 summary.stan_tables_precision<-function(object,timeinterval=1,digits=3,parmatrices=TRUE,priorcheck=TRUE,residualcov=TRUE,...) {
   stopifnot(!priorcheck,!residualcov)
   list(parmatrices=data.frame(matrix='DRIFT',row=1,col=1,Mean=round(1.23456789,digits)),npars=3)
 }
 registerS3method('summary','stan_tables_precision',summary.stan_tables_precision,envir=asNamespace('base'))
 class(f)<-c('stan_tables_precision','ctStanFit')
 x<-ctmaTables(f,parameters=drift,statistics='mean',format='none')
 expect_equal(x$data$Mean,1.23456789);expect_identical(x$display$M,'1.2346')
 old<-summary(f,priorcheck=FALSE,residualcov=FALSE)
 y<-ctmaTables(f,summary=old,parameters=drift,statistics='mean',format='none')
 expect_equal(y$data$Mean,1.235)
})

test_that('native discrete-time matrices are not exponentiated as CT drift',{
 f<-stan_tables_fixture();f$ctstanmodel<-list(continuoustime=FALSE)
 result<-list(data=data.frame(matrix='DRIFT',row=1,col=1,parameter='a',Mean=.8),display=data.frame(Parameter='a',Mean='0.8'))
 expect_error(.ctmaTablesAddDT(result,f,1,NULL,FALSE,4,'error'),'continuous-time model')
})
test_that('an augmented Stan cache uses full native DT summaries without replacing CT rows',{
 f<-stan_tables_fixture();f$stanfit$transformedpars$pop_DRIFT<-array(.1,c(4,3,3))
 d<-data.frame(matrix='DRIFT',row=c(1,1,2,2),col=c(1,2,1,2),Mean=c(-.4,.1,-.1,-.5))
 f$summary<-list(parmatrices=d,npars=3)
 calls<-numeric()
 local_mocked_bindings(.ctmaTablesRunSummary=function(target,lag,digits,label) {
   calls<<-c(calls,lag);dt<-d;dt$matrix<-'dtDRIFT';dt$Mean<-lag+seq_len(4)/100
   ct<-d;ct$Mean<-999
   list(parmatrices=rbind(ct,dt),npars=3)
 })
 x<-ctmaTables(f,parameters=drift,statistics='mean',timeIntervals=5,recomputeSummary='yes',format='none')
 expect_equal(calls,5)
 expect_equal(sort(x$data$Mean[x$data$matrix=='DRIFT']),sort(d$Mean))
 expect_equal(sort(x$data$Mean[x$data$matrix=='dtDRIFT@5']),5+seq_len(4)/100)
})
