library(testthat)
refresh_base <- function(bayes=FALSE) {
  d <- data.frame(row=c(1,1,2,2,1),col=c(1,2,1,2,1),
    Mean=c(-.3,.2,-.1,-.4,.5),sd=rep(.02,5),
    check.names=FALSE)
  d[['2.5%']] <- d$Mean-.04;d[['50%']] <- d$Mean+.001
  d[['97.5%']] <- d$Mean+.04
  rownames(d) <- c('DRIFT V1toV1 (invariant)','DRIFT V2toV1 (invariant)',
                  'DRIFT V1toV2 (invariant)','DRIFT V2toV2 (invariant)','asymCINT_1_1')
  if (bayes) { d$Rhat <- rep(1.01,5);d$n_eff <- rep(350,5) }
  structure(list(summary=list(estimates=d,minus2ll=123.456,n.parameters=5,
      uncertaintyNote=if(bayes)'Posterior draws' else 'Transformed Hessian draws'),
      argumentList=list(optimize=!bayes,priors=bayes)),class='CoTiMAFit')
}
refresh_init <- function(n=2L,ns=3L) {
 k<-n*n;ids<-seq_len(ns)*3L
 labels<-as.vector(t(outer(seq_len(n),seq_len(n),function(r,c)paste0('V',c,'toV',r))))
 vals<-matrix(seq_len(ns*k)/100,ns,k,byrow=TRUE)
 wide<-matrix('',ns,1L+2L*k);wide[,1L]<-paste('Reference',ids)
 wide[,seq.int(2L,2L*k,by=2L)]<-vals;wide[,seq.int(3L,2L*k+1L,by=2L)]<-.01
 colnames(wide)<-c('Reference',as.vector(rbind(labels,rep('SE',k))))
 rownames(wide)<-paste('Study No',ids)
 diff<-wide;diff[,seq.int(2L,2L*k,by=2L)]<-vals+1
 ci<-matrix(0,ns,6L*k)
 for(b in 0:2) {
  ci[,b*2L*k+seq.int(1L,2L*k,by=2L)]<-vals+b-.02
  ci[,b*2L*k+seq.int(2L,2L*k,by=2L)]<-vals+b+.02
 }
 structure(list(n.latent=n,n.studies=ns,statisticsList=list(originalStudyNumbers=ids),
   parameterNames=list(DRIFT=labels,T0VAR=paste0('t0_',seq_len(k))),
   argumentList=list(optimize=TRUE,priors=FALSE),
   modelResults=list(T0VAR=lapply(seq_len(ns),function(i)vals[i,]+2)),
   summary=list(estimates=wide,diff_estimates_original_time_scale=diff,
     confidenceIntervals=ci,scaleTime=1,minus2ll=123,n.parameters=ns*k)),class='CoTiMAFit')
}
library(testthat)
refresh_fixture <- function() {
 f<-refresh_base();f$summary$estimates[['50%']]<-NULL
 class(f)<-'ctStanFit';f
}
refresh_summary <- function(target,lag,digits,label) {
 d<-.ctmaTablesParameters(refresh_base()$summary$estimates,list())
 d[['50%']]<-d$Mean+.012345
 dt<-d[tolower(d$matrix)=='drift',];dt$matrix<-'dtDRIFT'
 dt$Mean<-lag+dt$row/10+dt$col/100
 dt[['50%']]<-dt$Mean+.001;dt[['2.5%']]<-dt$Mean-.123;dt[['97.5%']]<-dt$Mean+.456
 list(parmatrices=rbind(d,dt),npars=5,loglik=-61.728)
}
test_that('four decimals is the default and does not round returned numeric data',{
 expect_identical(formals(ctmaTables)$digits,4)
 f<-refresh_base();f$summary$estimates$Mean[1]<--.31234567
 x<-ctmaTables(f,parameters='drift[1,1]',statistics='estimate',format='none')
 expect_equal(x$data$Mean,-.31234567);expect_identical(x$display$M,'-0.3123')
})
test_that('missing statistics and DT uncertainty use native summaries at requested lags',{
 calls<-numeric()
 local_mocked_bindings(.ctmaTablesStanReady=function(fit,label)fit,.ctmaTablesRunSummary=function(target,lag,digits,label) {
   calls<<-c(calls,lag);refresh_summary(target,lag,digits,label)
 })
 f<-refresh_fixture();original<-f
 expect_error(ctmaTables(f,parameters=drift,statistics='median',format='none',interactive=FALSE),'recomputeSummary')
 x<-ctmaTables(f,parameters=drift,statistics=c('estimate','median','2.5','97.5'),timeIntervals=c(1,5),recomputeSummary='yes',format='none',interactive=FALSE)
 expect_equal(calls,c(1,5));expect_equal(nrow(x$data),12L)
 dt<-x$data[x$data$matrix=='dtDRIFT@5'&x$data$row==1&x$data$col==2,]
 expect_equal(dt$Mean,5.12);expect_equal(dt[['2.5%']],4.997);expect_equal(dt[['97.5%']],5.576)
 expect_identical(f,original);expect_length(x$generated_summaries,2L)
 bad<-f$summary;bad$estimates$Mean[1]<-999
 expect_error(ctmaTables(f,summary=bad,parameters=drift,statistics='median',recomputeSummary='yes',format='none'),'does not agree')
 expect_equal(length(calls),2L)
})
test_that('no calculations are made for complete summaries or when declined',{
 local_mocked_bindings(.ctmaTablesStanReady=function(fit,label)fit,.ctmaTablesRunSummary=function(...)stop('Unexpected summary'))
 f<-refresh_fixture()
 expect_s3_class(ctmaTables(f,parameters=drift,statistics='estimate',recomputeSummary='yes',format='none'),'ctmaTables')
 x<-ctmaTables(f,parameters=drift,statistics=c('estimate','median'),recomputeSummary='no',missing='note',format='none')
 expect_true(all(is.na(x$data$median)))
 # A raw point vector must never be mistaken for joint draws.
 f$stanfit<-list(rawest=c(-.3,.2,-.1,-.4))
 expect_null(.ctmaTablesDriftDraws(f))
})
test_that('all saved primary fits are summarized once per lag after a single consent',{
 f<-refresh_init(ns=2L);f$studyFitList<-list(`6`=refresh_fixture(),`3`=refresh_fixture())
 f$studyFitList[['6']]$studyID<-6;f$studyFitList[['3']]$studyID<-3
 calls<-list();consents<-0L
 local_mocked_bindings(.ctmaTablesStanReady=function(fit,label)fit,.ctmaTablesRunSummary=function(target,lag,digits,label) {
   calls[[length(calls)+1L]]<<-c(target$studyID,lag)
   refresh_summary(target,lag,digits,label)
 },.ctmaTablesSummaryConsent=function(...) {consents<<-consents+1L;TRUE})
 x<-ctmaTables(f,parameters=drift,statistics=c('mean','median'),timeIntervals=c(1,5),format='none')
 expect_equal(consents,1L);expect_equal(calls,list(c(3,1),c(3,5),c(6,1),c(6,5)))
 expect_equal(nrow(x$data),24L);expect_identical(unique(x$data$Study),c('Study No 3','Study No 6'))
 expect_length(x$generated_summaries,2L)
})
test_that('raw-draw transformation is delegated to the actual ctsem summary method',{
 local_mocked_bindings(.ctmaTablesStanReady=function(fit,label)fit)
 f<-refresh_fixture();f$stanfit<-list(rawest=c(1,2),rawposterior=matrix(1:20,10,2))
 summary.saved_tables_test<-function(object,timeinterval=1,digits=3,parmatrices=TRUE,priorcheck=TRUE,residualcov=TRUE,...) {
   stopifnot(identical(object$stanfit$rawposterior,f$stanfit$rawposterior),!priorcheck,!residualcov,parmatrices)
   refresh_summary(object,timeinterval,digits,'saved')
 }
 registerS3method('summary','saved_tables_test',summary.saved_tables_test,envir=asNamespace('base'))
 class(f)<-c('saved_tables_test','ctStanFit')
 x<-ctmaTables(f,parameters=drift,statistics='mean',timeIntervals=5,recomputeSummary='yes',format='none')
 expect_equal(nrow(x$data),8L);expect_true(any(grepl('summary(saved_fit',x$notes,fixed=TRUE)))
})
