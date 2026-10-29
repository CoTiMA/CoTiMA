library(testthat)
tables_fixture <- function(bayes=FALSE) {
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
tables_call <- function(f=tables_fixture(),p='drift',s=c('estimate','median'),...)
  ctmaTables(f,parameters=p,statistics=s,format='none',interactive=FALSE,...)

test_that('selectors and statistics retain exact source values and ordering', {
  x <- tables_call();expect_s3_class(x,'ctmaTables');expect_equal(nrow(x$data),4)
  a <- tables_call(p='V1toV2');b <- tables_call(p='DRiFt[2,1]')
  expect_equal(a$data$Mean,-.1);expect_equal(b$data$Mean,-.1)
  expect_equal(a$data[['50%']],-.099)
  x <- tables_call(p=c('V1toV1','V1toV2'),defaultOrder=FALSE)
  expect_equal(x$data$Mean,c(-.3,-.1))
  x <- tables_call(p=c('drift[1,1]','DRIFT[1,1]'),defaultOrder=FALSE)
  expect_equal(nrow(x$data),1)
  x <- tables_call(p='asymCINT');expect_equal(x$data$Mean,.5)
  expect_true(any(grepl('equilibrium positions',x$notes,fixed=TRUE)))
})
test_that('invalid selections and ambiguous names fail safely', {
  for(p in c('drift[0,1]','drift[9,1]','doesnotexist','drift[1,1];stop(1)'))
    expect_error(tables_call(p=p))
  expect_error(ctmaTables(tables_fixture(),statistics='estimate',format='none',interactive=FALSE))
  expect_error(ctmaTables(tables_fixture(),parameters='drift',format='none',interactive=FALSE))
  expect_error(tables_call(p=character()))
  f <- tables_fixture();rownames(f$summary$estimates)[2] <- 'DRIFT V1toV1'
  expect_error(tables_call(f,p='V1toV1'),'[Aa]mbigu')
})
test_that('supplied summary is compared to the fitted values', {
  f <- tables_fixture();expect_s3_class(tables_call(f,summary=f$summary),'ctmaTables')
  s <- f$summary;s$estimates$Mean[1] <- 900
  expect_error(tables_call(f,summary=s))
  s <- f$summary;s$minus2ll <- 999
  expect_error(tables_call(f,summary=s))
  valid <- tables_call(f)$summary
  valid$estimates$Mean[1] <- 900
  expect_error(tables_call(f,summary=valid))
})
test_that('Bayesian and unavailable statistics are not relabelled', {
  f <- tables_fixture(TRUE)
  x <- tables_call(f,s=c('estimate','sd','Rhat','ESS'))
  expect_true('sd' %in% names(x$data));expect_true('Rhat' %in% names(x$data))
  expect_error(tables_call(f,s='SE'))
  expect_error(tables_call(f,s='bulk_ess'))
  x <- tables_call(f,s='bulk_ess',missing='note')
  expect_true(any(grepl('bulk',x$notes,ignore.case=TRUE)))
  expect_false(any(grepl('confidence interval',x$notes,ignore.case=TRUE)))
  expect_equal(x$metadata$minus2ll,123.456)
  f$summary$minus2ll <- NULL;f$summary$logposterior <- -99
  expect_true(is.null(tables_call(f)$metadata$minus2ll) || is.na(tables_call(f)$metadata$minus2ll))
})
