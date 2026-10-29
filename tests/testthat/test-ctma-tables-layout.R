library(testthat)
layout_fit <- function() {
 d<-data.frame(matrix=c(rep('DRIFT',4),'CINT'),row=c(1,1,2,2,1),col=c(1,2,1,2,1),
   parameter=c('a11','a12','a21','a22','intercept'),Mean=c(-.3,.2,-.1,-.4,.5),sd=rep(.02,5),check.names=FALSE)
 structure(list(summary=list(estimates=d,minus2ll=100,n.parameters=5),argumentList=list(optimize=TRUE,priors=FALSE)),class='CoTiMAFit')
}
layout_init <- function() {
 f<-layout_fit();wide<-matrix(c('Study A',-.3,.02,.2,.02,-.1,.02,-.4,.02,
   'Study B',-.6,.03,.4,.03,-.2,.03,-.8,.03),2,9,byrow=TRUE)
 colnames(wide)<-c('Reference','a11','SE','a12','SE','a21','SE','a22','SE');rownames(wide)<-c('Study No 7','Study No 9')
 structure(list(n.latent=2,n.studies=2,statisticsList=list(originalStudyNumbers=c(7,9)),
   summary=list(estimates=wide,scaleTime=1),argumentList=f$argumentList),class='CoTiMAFit')
}
layout_table <- function(f=layout_fit(),...)ctmaTables(f,parameters=c('drift','CINT'),statistics=c('mean','sd'),format='none',interactive=FALSE,...)
test_that('page orientation and row arrangement are independent',{
 x<-layout_table();y<-layout_table(layout='wide')
 expect_identical(x$data,y$data);expect_identical(x$display,y$display)
 expect_identical(x$layout$orientation,'portrait');expect_identical(y$layout$orientation,'landscape')
 expect_equal(nrow(x$display),5L)
 expect_error(layout_table(layout='automatic'),'arg')
 expect_error(layout_table(rowLayout='unknown'),'arg')
})
test_that('initial rows sort by study or by parameter without changing values',{
 call<-function(...)ctmaTables(layout_init(),parameters=drift,statistics='mean',format='none',defaultOrder=FALSE,...)
 study<-call();parameter<-call(sortBy='parameter')
 expect_identical(study$data$Study,rep(c('Study No 7','Study No 9'),each=4L))
 expect_identical(parameter$data$Study,rep(c('Study No 7','Study No 9'),4L))
 expect_equal(parameter$data$Mean,c(-.3,-.6,.2,.4,-.1,-.2,-.4,-.8))
 expect_equal(sort(study$data$Mean),sort(parameter$data$Mean))
 expect_identical(layout_table(sortBy='parameter')$data,layout_table()$data)
})
test_that('compact rows retain every selected parameter and statistic',{
 base<-layout_table(defaultOrder=FALSE)
 all<-layout_table(rowLayout='all',defaultOrder=FALSE)
 drift<-layout_table(rowLayout='drift',defaultOrder=FALSE)
 expect_identical(all$data,base$data);expect_identical(drift$data,base$data)
 expect_equal(nrow(all$display),1L);expect_equal(nrow(drift$display),2L)
 expect_equal(ncol(all$display),11L)
 expect_identical(all$display[['a12 | M']],'0.2000')
 expect_identical(all$display[['intercept | SD']],'0.0200')
 expect_identical(drift$display$Parameter,c('DRIFT','intercept'))
 expect_identical(drift$display[['a12 | M']],c('0.2000',''))
 x<-ctmaTables(layout_init(),parameters=drift,statistics=c('mean','sd'),rowLayout='all',format='none')
 expect_equal(nrow(x$display),2L);expect_identical(x$display$Study,c('Study No 7','Study No 9'))
})
test_that('DT drifts get separate compact rows for each lag',{
 draws<-array(rep(matrix(c(-.3,-.1,.2,-.4),2),2),c(2,2,2))
 # Use explicit draws with identical draw slices to avoid a native engine dependency.
 draws[1,,]<-draws[2,,]<-matrix(c(-.3,-.1,.2,-.4),2)
 x<-ctmaTables(layout_fit(),parameters=drift,statistics='mean',timeIntervals=c(1,5),driftDraws=draws,rowLayout='drift',format='none')
 expect_equal(nrow(x$display),3L)
 expect_identical(x$display$Parameter,c('DRIFT','DT DRIFT (lag 1)','DT DRIFT (lag 5)'))
 expect_equal(nrow(x$data),12L)
})
test_that('column panels cover all statistic columns exactly once',{
 x<-layout_table(rowLayout='all',layout='long');y<-layout_table(rowLayout='all',layout='wide')
 expect_gt(length(x$panels),1L);expect_lte(length(y$panels),length(x$panels))
 cols<-unlist(lapply(x$panels,function(p)setdiff(names(p),'Parameter')),use.names=FALSE)
 expect_identical(cols,setdiff(names(x$display),'Parameter'));expect_false(anyDuplicated(cols)>0)
 for(p in x$panels)expect_true(all(.ctmaTablesWidths(p,'long')$widths>=.55))
 expect_true(any(grepl('no selected parameter or statistic was dropped',x$notes,fixed=TRUE)))
})
test_that('all exporters honor orientation and grouped headers',{
 directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
 for(layout in c('long','wide')) {
  x<-layout_table(rowLayout='all',layout=layout)
  paths<-.ctmaTablesWrite(x$display,x$notes,c('tex','xml'),directory,layout,'Compact',layout=layout)
  tex<-readLines(paths[['tex']]);xml<-readLines(paths[['xml']])
  expect_equal(any(grepl(',landscape',tex,fixed=TRUE)),layout=='wide')
  expect_equal(any(grepl('w:orient="landscape"',xml,fixed=TRUE)),layout=='wide')
  expect_true(any(grepl('multicolumn',tex,fixed=TRUE)))
  expect_true(any(grepl('gridSpan',xml,fixed=TRUE)))
  if(requireNamespace('officer',quietly=TRUE)&&requireNamespace('flextable',quietly=TRUE)&&requireNamespace('xml2',quietly=TRUE)) {
   path<-.ctmaTablesWrite(x$display,x$notes,'docx',directory,paste0(layout,'-word'),'Compact',layout=layout)
   unpack<-file.path(directory,paste0(layout,'-unpacked'));utils::unzip(path,exdir=unpack)
   doc<-xml2::read_xml(file.path(unpack,'word/document.xml'));ns<-xml2::xml_ns(doc)
   orientation<-xml2::xml_attr(xml2::xml_find_first(doc,'.//w:sectPr/w:pgSz',ns),'orient')
   expect_identical(orientation,if(layout=='wide')'landscape' else 'portrait')
   expect_true(length(xml2::xml_find_all(doc,'.//w:gridSpan',ns))>0)
  }
 }
})
test_that('negative and large numeric values fit without breaking their signs',{
 f<-layout_fit();f$summary$estimates$Mean[1]<--1234567.8912
 x<-layout_table(f,rowLayout='all')
 for(p in x$panels) {
   w<-.ctmaTablesWidths(p,'long')
   expect_true(all(w$widths[w$numeric]+1e-9>=w$required))
 }
 tex<-.ctmaTablesRenderTex(x$panels,x$notes,'Signed values','long')
 expect_true(any(grepl('\\mbox{-1234567.8912}',tex,fixed=TRUE)))
 initial<-ctmaTables(layout_init(),parameters=drift,statistics=c('mean','sd'),rowLayout='drift',layout='wide',format='none')
 expect_false('Parameter' %in% names(initial$display))
 expect_length(initial$panels,1L)
})
test_that('compact comparison rows retain both models and disambiguate custom labels',{
 f<-layout_fit();g<-layout_fit();g$summary$estimates$Mean<-g$summary$estimates$Mean+1
 x<-ctmaTables(list(First=f,Second=g),parameters=drift,statistics='mean',rowLayout='all',format='none')
 expect_equal(nrow(x$display),1L)
 expect_equal(x$display[['a12 | First: M']],'0.2000')
 expect_equal(x$display[['a12 | Second: M']],'1.2000')
 y<-layout_table(rowLayout='all',labels=c(a11='Same label',a12='Same label'))
 expect_false(anyDuplicated(names(y$display))>0)
 expect_true(any(grepl('DRIFT:1:1',names(y$display),fixed=TRUE)))
})
