library(testthat)
tables_update_fixture <- function(n=2L,ns=3L) {
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
test_that('bare case-insensitive matrix families select all entries',{
 f<-tables_update_fixture()
 a<-ctmaTables(f,parameters=c(drift),statistics='estimate',format='none',interactive=FALSE)
 expect_equal(nrow(a$data),12L);expect_true(all(a$data$matrix=='DRIFT'))
 b<-ctmaTables(f,parameters=c(DRIft,DIFFuSION),statistics='estimate',format='none',interactive=FALSE)
 expect_equal(nrow(b$data),24L)
 expect_setequal(unique(b$data$matrix),c('DRIFT','DIFFUSIONcov'))
 quoted<-ctmaTables(f,parameters=c('drift','diffusion'),statistics='estimate',format='none',interactive=FALSE)
 expect_identical(b$data,quoted$data)
 drift<-999;diffusion<-function()stop('must not run')
 expect_equal(nrow(ctmaTables(f,parameters=c(drift,diffusion),statistics='estimate',format='none')$data),24L)
 selected<-c('drift','diffusion')
 expect_identical(ctmaTables(f,parameters=selected,statistics='estimate',format='none')$data,b$data)
 expect_identical(ctmaTables(f,parameters=tolower(c('DRIFT','DIFFUSION')),statistics='estimate',format='none')$data,b$data)
 mixed<-ctmaTables(f,parameters=c(drift,'T0VAR[1,1]'),statistics='estimate',format='none')
 expect_equal(nrow(mixed$data),15L)
 calls<-0L
 ctmaTables(f,parameters={calls<-calls+1L;selected},statistics='estimate',format='none')
 expect_equal(calls,1L)
})
test_that('matrix families also expand in conventional fitted models',{
 f<-structure(list(argumentList=list(optimize=TRUE,priors=FALSE),summary=list(estimates=data.frame(
  matrix=rep(c('DRIFT','DIFFUSIONcov'),each=4),row=rep(c(1,1,2,2),2),col=rep(c(1,2,1,2),2),Mean=1:8,
  parameter=c('drift',paste0('p',2:8))))),class='CoTiMAFit')
 a<-ctmaTables(f,parameters=c(drift,diffusion),statistics='estimate',format='none')
 expect_equal(nrow(a$data),8L)
 expect_equal(nrow(ctmaTables(f,parameters=Drift,statistics='estimate',format='none')$data),4L)
})
test_that('initial fits preserve every primary and matrix coordinate',{
 for(n in c(1L,2L,3L)) {
  f<-tables_update_fixture(n,ns=2L)
  x<-ctmaTables(f,parameters=c(drift,diffusion,T0VAR),statistics='estimate',format='none',defaultOrder=FALSE)
  expect_equal(nrow(x$data),2L*3L*n*n)
  expect_identical(unique(x$data$Study),c('Study No 3','Study No 6'))
  expect_true(all(c('Study','Parameter','M') %in% names(x$display)))
  expect_equal(x$metadata$minus2ll,123)
  expect_true(any(grepl('totals across',x$notes,fixed=TRUE)))
 }
 f<-tables_update_fixture(n=3L,ns=1L)
 x<-ctmaTables(f,parameters='drift[2,3]',statistics=c('estimate','2.5','97.5'),format='none')
 expect_equal(x$data$Mean,.06);expect_equal(x$data[['2.5%']],.04);expect_equal(x$data[['97.5%']],.08)
})
test_that('initial summaries validate provenance and do not manufacture missing statistics',{
 f<-tables_update_fixture();s<-f$summary
 expect_s3_class(ctmaTables(f,summary=s,parameters=drift,statistics='estimate',format='none'),'ctmaTables')
 s$estimates[1,2]<-'999'
 expect_error(ctmaTables(f,summary=s,parameters=drift,statistics='estimate',format='none'),'does not agree')
 expect_error(ctmaTables(f,parameters=drift,statistics='median',format='none'),'Unavailable')
 x<-ctmaTables(f,parameters=drift,statistics=c('estimate','median'),format='none',missing='note')
 expect_true(all(is.na(x$data$median)));expect_true(any(grepl('do not retain medians',x$notes,fixed=TRUE)))
 expect_error(ctmaTables(f,parameters=T0VAR,statistics='sd',format='none'),'not retained')
 x<-ctmaTables(f,parameters=T0VAR,statistics='sd',format='none',missing='note')
 expect_true(all(is.na(x$data$sd)))
 f$summary$estimates<-f$summary$estimates[3:1,,drop=FALSE]
 expect_error(ctmaTables(f,parameters=drift,statistics='estimate',format='none'),'different orders')
})
test_that('initial diffusion units are reconciled without changing drift units',{
 f<-tables_update_fixture(ns=1L);f$summary$scaleTime<-.1
 cols<-seq.int(2L,8L,by=2L);sdcols<-cols+1L
 f$summary$diff_estimates_original_time_scale[,cols]<-as.numeric(f$summary$diff_estimates_original_time_scale[,cols])*.1
 f$summary$diff_estimates_original_time_scale[,sdcols]<-.001
 x<-ctmaTables(f,parameters=c(drift,diffusion),statistics=c('estimate','sd'),format='none')
 expect_equal(x$data$Mean[x$data$matrix=='DIFFUSIONcov'],c(1.01,1.02,1.03,1.04))
 expect_equal(x$data$sd[x$data$matrix=='DIFFUSIONcov'],rep(.01,4))
 expect_equal(sort(x$data$Mean[x$data$matrix=='DRIFT']),c(.01,.02,.03,.04))
 expect_true(any(grepl('recorded multiplier 0.1',x$notes,fixed=TRUE)))
})
test_that('packaged ctmaInit summary can be used without retained primary fit objects',{
 pkg<-getNamespaceName(environment(ctmaTables));e<-new.env()
 utils::data('CoTiMAInitFit_3',package=pkg,envir=e);f<-e$CoTiMAInitFit_3
 x<-ctmaTables(f,parameters=c(drift,diffusion),statistics=c('estimate','2.5','97.5'),format='none')
 expect_equal(nrow(x$data),24L)
 expect_identical(unique(x$data$Study),c('Study No 1','Study No 4','Study No 313'))
 z<-x$data$Study=='Study No 4'&x$data$matrix=='DRIFT'&x$data$row==2&x$data$col==1
 expect_equal(x$data$Mean[z],as.numeric(f$summary$estimates[2,6]))
 z<-x$data$Study=='Study No 313'&x$data$matrix=='DRIFT'&x$data$row==1&x$data$col==2
 expect_equal(x$data[['2.5%']][z],unname(f$summary$confidenceIntervals[3,3]))
})

test_that('legacy Bayesian primary summaries preserve their method and study roster',{
 pkg<-getNamespaceName(environment(ctmaTables));e<-new.env()
 utils::data('CoTiMAInitFit_6_NUTS',package=pkg,envir=e)
 x<-ctmaTables(e$CoTiMAInitFit_6_NUTS,parameters=DRIFT,statistics='estimate',format='none')
 expect_equal(nrow(x$data),24L)
 expect_identical(x$metadata$method,'Bayesian')
 expect_length(unique(x$data$Study),6L)
})
