library(testthat)
tables_fixture <- function() {
 d<-data.frame(row=c(1,1,2,2),col=c(1,2,1,2),Mean=c(-.3,.2,-.1,-.4),sd=.02,check.names=FALSE)
 d[['2.5%']]<-d$Mean-.04;d[['25%']]<-d$Mean-.01;d[['50%']]<-d$Mean;d[['97.5%']]<-d$Mean+.04
 rownames(d)<-c('DRIFT V1toV1 (invariant)','DRIFT V2toV1 (invariant)','DRIFT V1toV2 (invariant)','DRIFT V2toV2 (invariant)')
 structure(list(summary=list(estimates=d,minus2ll=123,n.parameters=4),argumentList=list(optimize=TRUE,priors=FALSE)),class='CoTiMAFit')
}
ct <- function(f=tables_fixture(),...)ctmaTables(f,parameters='drift',statistics=c('estimate','sd','2.5','median','97.5'),format='none',interactive=FALSE,...)
test_that('distinct quantiles and statistical methods retain their definitions',{
 f<-tables_fixture();x<-ctmaTables(f,parameters='drift',statistics=c('2.5','25%'),format='none',interactive=FALSE)
 expect_true(all(c('2.5%','25%') %in% names(x$data)))
 expect_false(identical(x$data[['2.5%']],x$data[['25%']]))
 expect_identical(ct(f)$metadata$method,'ML')
 f$argumentList$priors<-TRUE;expect_identical(ct(f)$metadata$method,'MAP')
 f$argumentList$optimize<-FALSE;expect_identical(ct(f)$metadata$method,'Bayesian')
 f$argumentList<-NULL;expect_identical(ct(f)$metadata$method,'unspecified')
 f$summary$minus2ll<-NULL;f$summary$loglik<-10;expect_equal(ct(f)$metadata$minus2ll,-20)
 f$summary$loglik<-NULL;f$summary$logposterior<-10;expect_true(is.na(ct(f)$metadata$minus2ll))
})
test_that('DT effects transform joint draws rather than endpoints',{
 f<-tables_fixture();draws<-array(0,c(3,2,2));A<-matrix(f$summary$estimates$Mean,2,2,byrow=TRUE)
 for(i in 1:3) draws[i,,]<-A+diag(c(-.02,.03)*(i-2))
 x<-ct(f,timeIntervals=c(1,2),driftDraws=draws)
 expect_equal(nrow(x$data),12)
 z<-which(x$data$matrix=='dtDRIFT@2'&x$data$row==2&x$data$col==1)
 direct<-vapply(1:3,function(i)expm::expm(draws[i,,]*2)[2,1],numeric(1))
 expect_equal(x$data$Mean[z],mean(direct));expect_equal(x$data$sd[z],sd(direct))
 expect_equal(x$data[['2.5%']][z],unname(quantile(direct,.025)))
 expect_equal(x$data[['50%']][z],median(direct))
 labelled<-ct(f,timeIntervals=1,driftDraws=draws,labels=c(V1toV2='X to Y'))
 expect_true(any(labelled$display$Parameter=='DT X to Y (lag 1)'))
 expect_error(ct(f,timeIntervals=1),'joint transformed drift draws')
 y<-ct(f,timeIntervals=1,dtPointOnly=TRUE)
 expect_true(all(is.na(y$data$sd[y$data$matrix=='dtDRIFT@1'])))
 for(v in list(c(0,1),c(1,1),NA_real_,Inf))expect_error(ct(f,timeIntervals=v,driftDraws=draws))
 expect_error(ct(f,timeIntervals=1,driftDraws=array(1,c(3,3,3))),'dimensions')
 expect_error(ct(f,timeIntervals=1,driftDraws=array(NA_real_,c(3,2,2))))
})
test_that('model comparisons align coordinates and preserve missing cells',{
 f<-tables_fixture();g<-f;g$summary$estimates$Mean<-g$summary$estimates$Mean*2
 g$summary$estimates<-g$summary$estimates[c(4,3,1),];names(g$summary$estimates)[names(g$summary$estimates)=='sd']<-'Sd'
 x<-ctmaTables(list(A=f,B=g),parameters='drift',statistics=c('estimate','sd'),format='none',interactive=FALSE)
 expect_equal(nrow(x$display),4);expect_equal(nrow(x$data),7)
 expect_true(all(c('A: M','B: M','A: SD','B: SD') %in% names(x$display)))
 expect_equal(sum(x$display[['B: M']]==''),1)
 expect_equal(x$metadata$A$minus2ll,123)
 expect_error(ctmaTables(list(f,g),parameters='drift',statistics='estimate',modelNames=c('a','a'),format='none'))
 expect_error(ctmaTables(list(f,g),parameters='drift',statistics='estimate',summary=list(f$summary),format='none'))
 expect_error(ctmaTables(list(f,g),parameters=NULL,statistics='estimate',format='none'))
})
test_that('export validation occurs before publishing any files',{
 directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
 x<-ct();d<-x$display;notes<-x$notes
 for(v in c('../escape','a/b','a..b','/absolute'))expect_error(.ctmaTablesWrite(d,notes,'tex',directory,v,'Test'))
 expect_error(.ctmaTablesWrite(d,notes,c('tex','bad'),directory,'valid','Test'))
 expect_length(list.files(directory,all.files=FALSE),0)
 p<-.ctmaTablesWrite(d,notes,'tex',directory,'valid','Title & percent % _ \\ { } $ # ~ ^')
 expect_true(file.exists(p));expect_error(.ctmaTablesWrite(d,notes,'tex',directory,'valid','Test'))
 before<-readLines(p);expect_true(any(grepl('\\textbackslash{}',before,fixed=TRUE)))
 .ctmaTablesWrite(d,notes,'tex',directory,'valid','Updated',overwrite=TRUE)
 expect_false(identical(before,readLines(p)))
 if(requireNamespace('xml2',quietly=TRUE)) {
  d$Parameter[1]<-'A & B < C "D"';p<-.ctmaTablesWrite(d,'Note & check','xml',directory,'xml','Title < & >')
  xml<-xml2::read_xml(p);expect_true(inherits(xml,'xml_document'))
  expect_true(grepl('A & B < C "D"',xml2::xml_text(xml),fixed=TRUE))
 }
})
test_that('failed publication restores all previously existing files',{
 directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
 oldtex<-file.path(directory,'restore.tex');oldxml<-file.path(directory,'restore.xml')
 writeLines('ORIGINAL TEX',oldtex);writeLines('ORIGINAL XML',oldxml)
 writer<-.ctmaTablesWrite;e<-new.env(parent=environment(writer))
 e$file.rename<-function(from,to) {
   if(basename(to)=='restore.xml'&&!startsWith(basename(from),'previous-'))return(FALSE)
   base::file.rename(from,to)
 }
 environment(writer)<-e
 expect_error(writer(ct()$display,'Test',c('tex','xml'),directory,'restore','Test',overwrite=TRUE),'publish')
 expect_true(file.exists(oldtex));expect_true(file.exists(oldxml))
 if(file.exists(oldtex))expect_identical(readLines(oldtex),'ORIGINAL TEX')
 if(file.exists(oldxml))expect_identical(readLines(oldxml),'ORIGINAL XML')
})
test_that('backups remain recoverable if the filesystem refuses restoration',{
 directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
 writeLines('ORIGINAL XML',file.path(directory,'retain.xml'))
 writer<-.ctmaTablesWrite;e<-new.env(parent=environment(writer))
 e$file.rename<-function(from,to) {
   if(basename(to)=='retain.xml')return(FALSE)
   base::file.rename(from,to)
 }
 environment(writer)<-e
 expect_warning(expect_error(writer(ct()$display,'Test','xml',directory,'retain','Test',overwrite=TRUE),'publish'),'backups are retained')
 paths<-list.files(directory,all.files=TRUE,full.names=TRUE,no..=TRUE)
 backups<-file.path(paths[startsWith(basename(paths),'.ctmatables-')],'previous-retain.xml')
 expect_length(backups,1)
 expect_identical(readLines(backups),'ORIGINAL XML')
})
