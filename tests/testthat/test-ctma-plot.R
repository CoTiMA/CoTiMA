.plot_fixture <- function(type="drift", studies=1L, nl=2L) {
  labels <- paste0("effect", seq_len(nl^2))
  a <- diag(-seq(.2,.4,length.out=nl),nl); a[row(a)!=col(a)] <- .01
  values <- matrix(rep(as.numeric(a),each=studies),studies,nl^2,dimnames=list(NULL,labels))
  avg <- colMeans(values)
  fixed <- rbind(MeanOfDriftValues=avg,FixedEffect_Drift=avg,
    FixedEffect_DriftLowerLimit=avg-.05,FixedEffect_DriftUpperLimit=avg+.05)
  structure(list(plot.type=type,n.latent=nl,n.studies=studies,activeDirectory=tempdir(),
    studyList=lapply(seq_len(studies),function(i) list(originalStudyNo=i,delta_t=i+1)),
    parameterNames=list(DRIFT=labels),modelResults=list(DRIFT=if(studies==1) avg else values,
      DRIFTSE=matrix(.03,studies,nl^2,dimnames=list(NULL,labels))),
    statisticsList=list(allDeltas=c(1,2,3),allSampleSizes=rep(100,studies)),
    summary=list(estimates=list(`Fixed Effects of Drift Coefficients`=fixed))),class="CoTiMAFit")
}

.plot_run <- function(object, ...) {
  state <- new.env(); state$polygons <- state$axes <- state$plots <- list(); state$files <- character()
  f <- ctmaPlot
  env <- new.env(parent=environment(f))
  env$.ctmaPlotSave <- function(path) { state$files <- c(state$files,path); invisible(path) }
  env$.capture_polygon <- function(x,y=NULL,...) {
    state$polygons[[length(state$polygons)+1L]] <- list(x=x,y=y)
    graphics::polygon(x,y,...)
  }
  env$.capture_axis <- function(side,at=NULL,labels=TRUE,...) {
    state$axes[[length(state$axes)+1L]] <- list(side=side,at=at,labels=labels)
    graphics::axis(side,at=at,labels=labels,...)
  }
  env$plot <- function(x,...) {
    state$plots[[length(state$plots)+1L]] <- list(x=x,args=list(...))
    graphics::plot(x,...)
  }
  rewrite <- function(x) {
    if(missing(x)) return(quote(expr=))
    if(!is.call(x)) return(x)
    if(identical(x[[1]],quote(graphics::polygon))) x[[1]] <- as.name(".capture_polygon")
    if(identical(x[[1]],quote(graphics::axis))) x[[1]] <- as.name(".capture_axis")
    as.call(lapply(as.list(x),rewrite))
  }
  body(f) <- rewrite(body(f)); environment(f) <- env
  body(f) <- as.call(c(as.list(body(f)),quote(invisible(environment()))))
  grDevices::pdf(file=NULL); device <- grDevices::dev.cur()
  on.exit(if(device %in% grDevices::dev.list()) grDevices::dev.off(device),add=TRUE)
  state$local <- suppressMessages(do.call(f,c(list(ctmaFitObject=object),list(...))))
  state
}

test_that("forest diamonds encode aggregate confidence limits", {
  x <- .plot_fixture("forest",2L)
  for(width in c(.03,.3)) {
    table <- x$summary$estimates[[1]]
    table["FixedEffect_DriftLowerLimit",] <- table["FixedEffect_Drift",]-width
    table["FixedEffect_DriftUpperLimit",] <- table["FixedEffect_Drift",]+width
    x$summary$estimates[[1]] <- table
    r <- .plot_run(x)
    axes <- Filter(function(a) a$side==1,r$axes)
    expect_length(r$polygons,4L)
    for(i in 1:4) {
      ticks <- axes[[i]]
      raw <- function(v) ticks$labels[1]+(v-ticks$at[1])*diff(ticks$labels[1:2])/diff(ticks$at[1:2])
      expect_equal(raw(r$polygons[[i]]$x[c(1,3)]),unname(table[c("FixedEffect_DriftLowerLimit","FixedEffect_DriftUpperLimit"),i]),tolerance=1e-10)
    }
  }
})

test_that("funnel confidence boundaries have the correct SE slope", {
  r <- .plot_run(.plot_fixture("funnel",2L))
  for(p in r$polygons) expect_equal(unname(abs(p$x[1]-p$x[2])/p$y[1]),1.96)
})

test_that("time coordinates do not serve as array indices", {
  x <- .plot_fixture()
  for(range in list(c(10,12,1),c(.1,1,.1),c(1,144,10),c(0,3,1))) {
    r <- .plot_run(x,timeRange=range)
    p <- r$local$plotPairs[[1]]
    a <- matrix(x$modelResults$DRIFT,2)
    for(j in seq_len(dim(p)[2])) expect_equal(as.numeric(p[1,j,-1]),
      as.numeric(OpenMx::expm(a*p[1,j,1])),tolerance=1e-10)
  }
})

test_that("tick positions agree with the numeric labels", {
  x <- .plot_fixture(); x$statisticsList$allDeltas <- c(.5,1,2,3)
  r <- .plot_run(x)
  for(a in Filter(function(a) a$side==1,r$axes)) expect_equal(a$at,a$labels)
  r <- .plot_run(x,timeRange=c(.5,6,.5),xLabels=c(.5,2.5,6))
  for(a in Filter(function(a) a$side==1,r$axes)) expect_equal(a$at,c(.5,2.5,6))
})

test_that("power-only plots select all stored times inside the requested window", {
  x <- .plot_fixture("power")
  times <- c(1:10,13.5,15:20)
  samples <- matrix(seq_len((length(times)+3L)*4L)+100,ncol=4,
    dimnames=list(c(as.character(times),"summaryA","summaryB","summaryC"),NULL))
  x$summary$estimates <- list(`Requested Statistical Power`=c(.5,.8),`Required Sample Sizes`=samples)
  for(range in list(c(3,6,1),c(5,20,1))) {
    r <- .plot_run(x,plotAutoEffects=FALSE,plotCrossEffects=FALSE,timeRange=range)
    expected <- times[times>=range[1]&times<=range[2]]
    expect_equal(r$local$usedTimeRange,expected)
    expect_equal(as.numeric(rownames(r$local$currentRequiredSamleSizes)),expected)
    expect_equal(unname(r$plots[[1]]$x[,1]),expected)
    expect_true(all(dirname(r$files)==tempdir()))
  }
  expect_error(.plot_run(x,plotAutoEffects=FALSE,plotCrossEffects=FALSE,timeRange=c(30,40,1)),"No stored power")
})

test_that("headless devices produce PNG files without manual display-list setup", {
  directory <- withr::local_tempdir()
  draw <- function() {
    grDevices::pdf(file=NULL); device <- grDevices::dev.cur()
    on.exit(if(device %in% grDevices::dev.list()) grDevices::dev.off(device))
    ctmaPlot(.plot_fixture("forest",studies=1,nl=1),activeDirectory=directory)
  }
  expect_no_error(draw())
  files <- list.files(directory,pattern="[.]png$",full.names=TRUE)
  expect_length(files,1L)
  expect_true(all(file.info(files)$size>1000))
})

test_that("multiple forest and funnel fits keep dimensions labels and files separate", {
  a <- .plot_fixture("forest",2L); b <- .plot_fixture("forest",3L,nl=3L)
  r <- .plot_run(structure(list(first=a,second=b),class="CoTiMAFit"))
  expect_length(r$files,13L); expect_equal(anyDuplicated(r$files),0L)
  a$plot.type <- b$plot.type <- "funnel"
  b$parameterNames$DRIFT <- paste0("other",1:9)
  r <- .plot_run(structure(list(a,b),class="CoTiMAFit"))
  expect_length(r$files,13L); expect_equal(anyDuplicated(r$files),0L)
  expect_true(any(grepl("fit2_Funnel Plot for other",r$files,fixed=TRUE)))
})

test_that("study matrices are not recycled to a single coefficient in curves", {
  x <- .plot_fixture(studies=2L)
  x$modelResults$DRIFT[2,] <- c(-.8,.05,-.02,-.3)
  r <- .plot_run(x,timeRange=c(1,3,1))
  for(h in 1:2) expect_equal(as.numeric(r$local$plotPairs[[1]][h,1,-1]),
    as.numeric(OpenMx::expm(matrix(x$modelResults$DRIFT[h,],2))),tolerance=1e-10)
})

test_that("invalid ranges are rejected and one-variable plots work", {
  for(range in list(c(1,2),c(1,2,0),c(3,1,1),c(NA,2,1))) {
    expect_error(.plot_run(.plot_fixture(),timeRange=range),"timeRange")
  }
  expect_error(.plot_run(.plot_fixture(),yLimitsForEffects=c(1,1)),"yLimitsForEffects")
  expect_no_error(.plot_run(.plot_fixture(nl=1),timeRange=c(1,3,1)))
  expect_no_error(.plot_run(.plot_fixture("forest",studies=1,nl=1)))
})

test_that("raw-data lag fallback preserves the study-to-interval mapping", {
  x <- .plot_fixture(studies=3L)
  x$studyList <- lapply(1:3,function(i) list(originalStudyNo=i,delta_t=NA_real_))
  x$statisticsList$allTpoints <- c(3,2,4)
  x$statisticsList$allDeltas <- c(1,3,10,20,30,40)
  r <- .plot_run(x,timeRange=c(1,40,2))
  expect_equal(r$local$allAvgDeltas[[1]],c(2,10,30))
  x$statisticsList$allTpoints <- NULL
  expect_error(.plot_run(x),"unambiguously")
})

test_that("categorical labels preserve signed and multi-digit values", {
  parser <- get(".ctmaPlotCategories",environment(ctmaPlot))
  expect_equal(parser(c("10 (category value)","100 (category value)","-10 (category value)")),c(10,100,-10))
  expect_error(parser("unknown category"),"category")
})

test_that("packaged legacy examples remain plottable", {
  package <- environmentName(environment(ctmaPlot))
  if(!package %in% c("CoTiMA","CoTiMA2")) package <- "CoTiMA2"
  for(name in c("CoTiMABiG_D_BO","CoTiMAFullFit_3","CoTiMAFullInv23Fit_6","CoTiMAPower_D_BO")) {
    path <- system.file("example-data",paste0(name,".rda"),package=package)
    if(!nzchar(path)) path <- system.file("data",paste0(name,".rda"),package=package)
    if(!nzchar(path)) next
    e <- new.env(); load(path,e); x <- e[[name]]
    expect_no_error(.plot_run(x,activeDirectory=tempdir(),timeRange=c(1,6,1)))
  }
})

test_that("moderator overlays use each fit's own categories and dot positions", {
  make <- function(categories) {
    x <- .plot_fixture()
    labels <- c("V1toV1","V2toV1","V1toV2","V2toV2")
    x$parameterNames$DRIFT <- labels
    x$mod.type <- "cat"
    x$modelResults$MOD <- matrix(.01,length(categories),4)
    x$summary$mod.effects <- matrix(0,length(categories),1,
      dimnames=list(paste(categories,"(category value)"),NULL))
    effects <- array(0,c(2,4,length(categories)))
    for(j in seq_along(categories)) effects[,,j] <- matrix(rep(c(-.05,.01,.02,-.05)*j,each=2),2)
    x$studyFitList <- list(ctstanmodelbase=list(pars=data.frame(matrix="DRIFT",param=labels,
      transform=c("-exp(param)","param","param","-exp(param)"))),
      stanfit=list(rawest=c(-1,0,0,-1),transformedparsfull=list(TIPREDEFFECT=effects)))
    x
  }
  a <- make(c(100,-10)); b <- make(10)
  r <- .plot_run(structure(list(a,b),class="CoTiMAFit"),timeRange=c(1,4,1))
  expect_equal(dim(r$local$plotPairs[[1]])[1],3L)
  expect_equal(dim(r$local$plotPairs[[2]])[1],2L)
  expect_equal(r$local$mod.values[[1]],c(-999,100,-10))
  expect_equal(r$local$mod.values[[2]],c(-999,10))
  expect_true(all(r$local$moderatorPositions[[1]][2:4] %in% r$local$usedTimeRange))
})

test_that("continuous moderator curves average draws and select only mod.number", {
  x <- .plot_fixture()
  labels <- c("V1toV1","V2toV1","V1toV2","V2toV2")
  x$parameterNames$DRIFT <- labels; x$mod.type <- "cont"
  x$modelResults$MOD <- matrix(0,4,1,dimnames=list(labels,NULL))
  effects <- array(0,c(2,4,2))
  effects[,,1] <- rbind(c(-.1,.01,.02,-.1),c(-.3,.03,.06,-.3))
  effects[,,2] <- effects[,,1]*2
  x$studyFitList <- list(ctstanmodelbase=list(pars=data.frame(matrix="DRIFT",param=labels,
    transform=c("-exp(param)","param","param","-exp(param)"))),
    stanfit=list(rawest=c(-1,0,0,-1),transformedparsfull=list(TIPREDEFFECT=effects)))
  for(number in 1:2) {
    expect_no_warning(r <- .plot_run(x,timeRange=c(1,4,1),mod.values=c(-1,0,1),mod.number=number))
    a <- matrix(c(-1,0,0,-1)+colMeans(effects[,,number]),2,byrow=TRUE)
    diag(a) <- -exp(diag(a))
    expect_equal(as.numeric(r$local$plotPairs[[1]][3,1,-1]),as.numeric(OpenMx::expm(a)),tolerance=1e-10)
  }
})
