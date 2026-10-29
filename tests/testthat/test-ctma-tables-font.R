library(testthat)
font_fit <- function() structure(list(summary=list(estimates=data.frame(
 matrix='DRIFT',row=1,col=1,parameter='a11',Mean=-.3,sd=.02),minus2ll=100,n.parameters=1),
 argumentList=list(optimize=TRUE,priors=FALSE)),class='CoTiMAFit')
test_that('font size validates and is retained without changing estimates',{
 call<-function(...)ctmaTables(font_fit(),parameters='drift',statistics='mean',format='none',...)
 a<-call();b<-call(fontSize=13)
 expect_equal(a$layout$fontSize,10);expect_equal(b$layout$fontSize,13)
 expect_identical(a$data,b$data);expect_identical(a$display,b$display)
 for(bad in list(NULL,NA,NaN,Inf,0,-1,'10',c(10,12)))expect_error(call(fontSize=bad),'fontSize')
})
test_that('font size reaches all formatted exporters',{
 d<-tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE))
 for(size in c(10,13)) {
  formats<-c('tex','xml')
  word<-all(vapply(c('officer','flextable','xml2'),requireNamespace,logical(1),quietly=TRUE))
  if(word)formats<-c(formats,'docx')
  z<-ctmaTables(font_fit(),parameters='drift',statistics='mean',fontSize=size,
    activeDirectory=d,file=paste0('font',size),format=formats)
  expect_true(any(grepl(paste0('\\fontsize{',size,'}{'),readLines(z$files[['tex']]),fixed=TRUE)))
  expect_true(any(grepl(paste0('<w:sz w:val="',2*size,'"'),readLines(z$files[['xml']]),fixed=TRUE)))
  if(word) {
   unpack<-file.path(d,paste0('unpack',size));utils::unzip(z$files[['docx']],exdir=unpack)
   doc<-xml2::read_xml(file.path(unpack,'word/document.xml'));ns<-xml2::xml_ns(doc)
   sizes<-xml2::xml_attr(xml2::xml_find_all(doc,'.//w:rPr/w:sz',ns),'val')
   expect_gt(length(sizes),0);expect_true(all(sizes==as.character(2*size)))
  }
 }
})
