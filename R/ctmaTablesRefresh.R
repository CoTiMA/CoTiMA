# Optional reconstruction of missing table quantities from saved fitted objects.
# Only summary() is called here; no estimation/optimization function is invoked.
.ctmaTablesSummaryTarget <- function(fit) {
  target<-.ctmaTablesNative(fit)
  if(inherits(target,'ctmaJuliaFit')) target<-target$ctsem_fit
  if(inherits(target,c('ctStanFit','ctJuliaFit'))) target else NULL
}
.ctmaTablesSummaryConsent <- function(policy,interactive,description) {
  if(policy=='yes')return(TRUE)
  if(policy=='no')return(FALSE)
  if(!isTRUE(interactive)) .ctmaTablesStop('Additional summaries are needed for ',description,
    '. Set recomputeSummary="yes" to permit summary calculations, or "no" to keep existing results. No model refitting is involved.')
  repeat {
    answer<-tolower(trimws(readline(paste0('Compute summary() from saved fits for ',description,'? No models will be refitted. [yes/no] '))))
    if(answer %in% c('yes','y'))return(TRUE)
    if(answer %in% c('no','n'))return(FALSE)
    if(!nzchar(answer)) .ctmaTablesStop('Summary calculation cancelled.')
    message('Please answer yes or no.')
  }
}
.ctmaTablesRunSummary <- function(target,lag,digits,label) {
  method<-NULL
  for(cl in class(target)) {
    method<-utils::getS3method('summary',cl,optional=TRUE)
    if(!is.null(method))break
  }
  if(is.null(method)||!'timeinterval' %in% names(formals(method)))
    .ctmaTablesStop('The saved fit for ',label,' has no supported summary(timeinterval=...) method.')
  if(inherits(target,'ctStanFit'))target<-.ctmaTablesStanReady(target,label)
  args<-list(object=target,timeinterval=lag,digits=max(8L,digits))
  fields<-names(formals(method))
  if('parmatrices' %in% fields)args$parmatrices<-TRUE
  if('priorcheck' %in% fields)args$priorcheck<-FALSE
  if('residualcov' %in% fields)args$residualcov<-FALSE
  message('CoTiMA tables: computing summary for ',label,' (timeinterval = ',lag,').')
  tryCatch({invisible(utils::capture.output(s<-do.call(base::summary,args)));.ctmaTablesTable(s);s},
    error=function(e).ctmaTablesStop('Could not summarise ',label,' at timeinterval=',lag,': ',conditionMessage(e)))
}
.ctmaTablesGenerateSummaries <- function(target,lags,digits,label) {
  if(is.null(lags))lags<-1
  if(!is.numeric(lags)||!length(lags)||anyNA(lags)||any(!is.finite(lags))||any(lags<=0)||anyDuplicated(lags))
    .ctmaTablesStop('timeIntervals must contain distinct finite positive lags in fitted model time units.')
  if(inherits(target,'ctStanFit'))target<-.ctmaTablesStanReady(target,label)
  out<-lapply(lags,function(lag).ctmaTablesRunSummary(target,lag,digits,label))
  stats::setNames(out,vapply(lags,function(x)format(x,digits=15),character(1)))
}
.ctmaTablesPrepare <- function(fit,provided,parameters,statistics,lags,draws,pointOnly,policy,interactive,digits) {
  stored<-.ctmaTablesSummary(fit,digits=if(is.null(provided))max(8L,digits) else 3L)
  if(!is.null(provided)&&!isTRUE(all.equal(.ctmaTablesSignature(provided),.ctmaTablesSignature(stored),check.attributes=TRUE)))
    .ctmaTablesStop('The supplied summary does not agree with the fit.')
  fit$summary<-stored
  d<-.ctmaTablesParameters(.ctmaTablesTable(stored),fit)
  available<-setdiff(names(d)[vapply(d,is.numeric,logical(1))],c('row','col'))
  requested<-if(is.null(statistics))character() else .ctmaTablesStats(statistics,available)
  selected<-if(is.null(parameters))seq_len(nrow(d)) else .ctmaTablesSelect(d,parameters)
  if(!is.null(lags)&&!any(tolower(d$matrix[selected])=='drift'))
    .ctmaTablesStop('Select DRIFT entries when requesting DT effects.')
  needCT<-length(setdiff(requested,available))>0L
  cached<-if(is.null(lags)||!is.null(draws)||isTRUE(pointOnly))NULL else .ctmaTablesDriftDraws(fit)
  driftRows<-d[tolower(d$matrix)=='drift',,drop=FALSE]
  expected<-if(nrow(driftRows))max(driftRows$row,driftRows$col) else 0L
  needDT<-!is.null(lags)&&is.null(draws)&&!isTRUE(pointOnly)&&
    (is.null(cached)||any(dim(cached)[2:3]!=expected))
  if(!needCT&&!needDT)return(fit)
  target<-.ctmaTablesSummaryTarget(fit)
  if(is.null(target))return(fit)
  description<-paste(c(if(needCT)'missing requested statistics',if(needDT)'DT effects at the requested lags'),collapse=' and ')
  if(!.ctmaTablesSummaryConsent(policy,interactive,description))return(fit)
  generated<-.ctmaTablesGenerateSummaries(target,if(needDT)lags else NULL,digits,'the fitted model')
  fit$summary<-if(needCT)generated[[1L]] else stored
  if(needDT)fit$ctma_tables_dt_summaries<-generated
  fit$ctma_tables_generated_summaries<-generated
  fit
}
.ctmaTablesDTFromSummary <- function(result,s,lag,digits,missing) {
  source<-.ctmaTablesParameters(.ctmaTablesTable(s),list())
  source<-source[tolower(source$matrix)=='dtdrift',,drop=FALSE]
  if(!nrow(source)) .ctmaTablesStop('summary(..., timeinterval=',lag,') did not return dtDRIFT estimates.')
  chosen<-which(tolower(result$data$matrix)=='drift')
  if(!length(chosen)) .ctmaTablesStop('Select DRIFT entries when requesting DT effects.')
  d<-result$data[chosen,,drop=FALSE]
  z<-match(paste(d$row,d$col),paste(source$row,source$col))
  if(anyNA(z)) .ctmaTablesStop('A requested DT matrix entry is absent from summary(..., timeinterval=',lag,').')
  stats<-setdiff(names(d),c('matrix','row','col','parameter'))
  available<-setdiff(names(source)[vapply(source,is.numeric,logical(1))],c('row','col'))
  map<-.ctmaTablesStats(stats,available);bad<-!map %in% available
  if(any(bad)&&missing=='error') .ctmaTablesStop('Requested statistics unavailable in the DT summary: ',paste(stats[bad],collapse=', '))
  d$matrix<-paste0('dtDRIFT@',format(lag,digits=15))
  d$parameter<-paste0('DT ',d$parameter,' (lag ',format(lag,digits=15),')')
  for(j in seq_along(stats))d[[stats[j]]]<-if(bad[j])rep(NA_real_,nrow(d)) else source[[map[j]]][z]
  display<-data.frame(Parameter=paste0('DT ',result$display$Parameter[chosen],' (lag ',format(lag,digits=15),')'),stringsAsFactors=FALSE)
  for(stat in stats)display[[stat]]<-ifelse(is.finite(d[[stat]]),formatC(d[[stat]],format='f',digits=digits),'')
  result$data<-rbind(result$data,d);result$display<-rbind(result$display,display)
  result$notes<-unique(c(result$notes,paste0('DT = discrete time. DT effects and their available uncertainty were obtained from summary(saved_fit, timeinterval=',lag,') in fitted model time units.'),
    if(any(bad))paste0('Unavailable DT statistics are blank: ',paste(stats[bad],collapse=', '),'.')))
  result$metadata$timeIntervals<-c(result$metadata$timeIntervals,lag)
  result
}

.ctmaTablesRefreshInit <- function(fit,input,parameters,statistics,lags,draws,pointOnly,policy,interactive,digits) {
  ns<-length(input$studies)
  out<-list(input=input,targets=vector('list',ns),generated=vector('list',ns),recomputed=FALSE)
  saved<-fit$studyFitList
  if(is.list(saved)&&length(saved)==ns) {
    ids<-as.character(fit$statisticsList$originalStudyNumbers)
    if(!is.null(names(saved))&&length(ids)==ns&&all(ids %in% names(saved)))saved<-saved[ids]
    out$targets<-lapply(saved,.ctmaTablesSummaryTarget)
  }
  needCT<-any(vapply(input$studies,function(x) {
    selected<-.ctmaTablesSelect(x$table,parameters)
    requested<-.ctmaTablesStats(statistics,names(x$available))
    any(vapply(requested,function(stat)is.null(x$available[[stat]])||any(!x$available[[stat]][selected]),logical(1)))
  },logical(1)))
  needDT<-!is.null(lags)&&is.null(draws)&&!isTRUE(pointOnly)&&
    any(vapply(out$targets,function(x) {
      a<-if(is.null(x))NULL else .ctmaTablesDriftDraws(x)
      is.null(a)||any(dim(a)[2:3]!=fit$n.latent)
    },logical(1)))
  if(!needCT&&!needDT)return(out)
  if(any(vapply(out$targets,is.null,logical(1)))) {
    if(policy=='yes') .ctmaTablesStop('Additional summaries require one saved native fit per primary study in studyFitList. No model will be refitted.')
    return(out)
  }
  description<-paste0('all ',ns,' saved primary-study fits (',paste(c(if(needCT)'missing requested statistics',if(needDT)'DT effects at the requested lags'),collapse=' and '),')')
  if(!.ctmaTablesSummaryConsent(policy,interactive,description))return(out)
  n<-fit$n.latent
  for(i in seq_len(ns)) {
    generated<-.ctmaTablesGenerateSummaries(out$targets[[i]],lags,digits,input$labels[i])
    d<-.ctmaTablesParameters(.ctmaTablesTable(generated[[1L]]),out$targets[[i]])
    d<-d[d$row<=n&d$col<=n,,drop=FALSE]
    old<-input$studies[[i]]$table
    z<-match(paste(tolower(d$matrix),d$row,d$col),paste(tolower(old$matrix),old$row,old$col))
    ok<-!is.na(z);d$parameter[ok]<-old$parameter[z[ok]]
    stats<-setdiff(names(d)[vapply(d,is.numeric,logical(1))],c('row','col'))
    availability<-as.data.frame(stats::setNames(lapply(stats,function(x)rep(TRUE,nrow(d))),stats))
    out$input$studies[[i]]<-list(table=d,available=availability)
    out$generated[[i]]<-generated # Keep each generated summary as an audit record.
  }
  out$recomputed<-TRUE
  out
}
