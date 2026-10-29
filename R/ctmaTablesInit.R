# Adapter for the aggregate summary retained by ctmaInit(), including old fits
# whose individual ctsem objects were not saved. No fitting is performed.
.ctmaTablesIsInit <- function(fit) {
  if(!inherits(fit,'CoTiMAFit')||!is.list(fit$summary)) return(FALSE)
  d<-fit$summary$estimates;n<-fit$n.latent
  is.numeric(n)&&length(n)==1L&&is.finite(n)&&n>=1&&n==floor(n)&&
    (is.matrix(d)||is.data.frame(d))&&ncol(d)==1L+2L*n*n&&
    !all(c('row','col') %in% colnames(d))
}
.ctmaTablesInitNumber <- function(x) {
  v<-suppressWarnings(as.numeric(x))
  bad<-is.na(v)&!is.na(x)&!tolower(trimws(as.character(x))) %in% c('','na','nan')
  if(any(bad)) .ctmaTablesStop('A ctmaInit estimate/uncertainty cell is not numeric.')
  v
}
.ctmaTablesInitTables <- function(fit,s) {
  n<-as.integer(fit$n.latent);k<-n*n;wide<-as.matrix(s$estimates);ns<-nrow(wide)
  if(!ns) .ctmaTablesStop('The ctmaInit summary has no primary-study rows.')
  if(!is.null(fit$n.studies)&&!identical(as.numeric(fit$n.studies),as.numeric(ns)))
    .ctmaTablesStop('Study count and ctmaInit summary rows do not agree.')
  ids<-rownames(wide)
  if(is.null(ids)||identical(ids,as.character(seq_len(ns)))) {
    nums<-fit$statisticsList$originalStudyNumbers
    if(length(nums)!=ns) nums<-seq_len(ns)
    ids<-paste('Study No',nums)
  }
  if(anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids)) .ctmaTablesStop('Primary-study row labels must be nonempty and unique.')
  nums<-fit$statisticsList$originalStudyNumbers
  if(length(nums)==ns&&all(startsWith(ids,'Study No '))&&
      !identical(ids,paste('Study No',nums)))
    .ctmaTablesStop('Primary-study summary rows and the retained study metadata have different orders.')
  scale<-s$scaleTime
  if(is.null(scale)) scale<-fit$argumentList$scaleTime
  if(is.null(scale)) scale<-1
  if(!is.numeric(scale)||length(scale)!=1L||!is.finite(scale)||scale<=0) .ctmaTablesStop('Invalid ctmaInit time multiplier.')
  fullCI<-s$confidenceIntervals
  if(!is.null(fullCI)) {
    fullCI<-as.matrix(fullCI)
    if(nrow(fullCI)!=ns||ncol(fullCI)!=6L*k) .ctmaTablesStop('Unsupported ctmaInit confidenceIntervals layout; expected paired DRIFT, DIFFUSION and T0VAR limits.')
  }
  ciBlock<-function(i,block,dedicated=NULL) {
    if(!is.null(fullCI)) return(.ctmaTablesInitNumber(fullCI[i,block*2L*k+seq_len(2L*k)]))
    if(!is.null(dedicated)) {
      z<-as.matrix(dedicated)
      if(nrow(z)!=ns||ncol(z)!=2L*k) .ctmaTablesStop('A dedicated ctmaInit interval table has incompatible dimensions.')
      return(.ctmaTablesInitNumber(z[i,])/scale)
    }
    NULL
  }
  namesFor<-function(w,field) {
    labels<-colnames(w)[seq.int(2L,2L*k,by=2L)]
    if(length(labels)!=k||anyNA(labels)||any(!nzchar(labels))) labels<-fit$parameterNames[[field]]
    if(length(labels)!=k||anyNA(labels)||any(!nzchar(labels))) labels<-rep(NA_character_,k)
    labels
  }
  add<-function(matrix,means,sd=NULL,ci=NULL,labels=NULL) {
    z<-data.frame(matrix=matrix,row=rep(seq_len(n),each=n),col=rep(seq_len(n),times=n),
      Mean=.ctmaTablesInitNumber(means),sd=if(is.null(sd))rep(NA_real_,k) else .ctmaTablesInitNumber(sd),check.names=FALSE)
    z[['2.5%']]<-if(is.null(ci))NA_real_ else ci[seq.int(1L,2L*k,by=2L)]
    z[['97.5%']]<-if(is.null(ci))NA_real_ else ci[seq.int(2L,2L*k,by=2L)]
    z$parameter<-if(is.null(labels)||length(labels)!=k) paste0(matrix,'[',z$row,',',z$col,']') else labels
    bad<-is.na(z$parameter)|!nzchar(z$parameter)
    z$parameter[bad]<-paste0(matrix,'[',z$row[bad],',',z$col[bad],']')
    attr(z,'available')<-data.frame(Mean=rep(TRUE,k),sd=rep(!is.null(sd),k),
      `2.5%`=rep(!is.null(ci),k),`97.5%`=rep(!is.null(ci),k),check.names=FALSE)
    z
  }
  out<-lapply(seq_len(ns),function(i) {
    drift<-add('DRIFT',wide[i,seq.int(2L,2L*k,by=2L)],wide[i,seq.int(3L,2L*k+1L,by=2L)],
      ciBlock(i,0L,s$drift_CI_original_time_scale),namesFor(wide,'DRIFT'))
    parts<-list(drift)
    diff<-s$diff_estimates_original_time_scale
    if(!is.null(diff)) {
      diff<-as.matrix(diff)
      if(nrow(diff)!=ns||ncol(diff)!=1L+2L*k) .ctmaTablesStop('The ctmaInit diffusion estimate table has incompatible dimensions.')
      means<-.ctmaTablesInitNumber(diff[i,seq.int(2L,2L*k,by=2L)])/scale
      raw<-fit$modelResults$DIFFUSION
      if(length(raw)) {
        if(!is.list(raw)||length(raw)!=ns||length(raw[[i]])!=k) .ctmaTablesStop('Retained diffusion coefficients do not match the primary-study layout.')
        means<-.ctmaTablesInitNumber(raw[[i]])
      }
      parts[[length(parts)+1L]]<-add('DIFFUSIONcov',means,
        .ctmaTablesInitNumber(diff[i,seq.int(3L,2L*k+1L,by=2L)])/scale,
        ciBlock(i,1L,s$diff_CI_original_time_scale),namesFor(diff,'DIFFUSION'))
    }
    t0<-fit$modelResults$T0VAR
    if(length(t0)) {
      if(!is.list(t0)||length(t0)!=ns||length(t0[[i]])!=k) .ctmaTablesStop('The retained T0VAR estimates do not match the primary-study layout.')
      parts[[length(parts)+1L]]<-add('T0cov',t0[[i]],ci=ciBlock(i,2L),labels=fit$parameterNames$T0VAR)
    }
    availability<-do.call(rbind,lapply(parts,attr,'available'))
    table<-do.call(rbind,parts);rownames(table)<-NULL;attr(table,'available')<-NULL
    list(table=table,available=availability)
  })
  list(studies=out,labels=ids,sources=as.character(wide[,1L]),scale=scale)
}
.ctmaTablesInit <- function(fit,provided,parameters,statistics,defaultOrder,digits,labels,
    missing,timeUnit,tableNote,interactive,timeIntervals,driftDraws,dtPointOnly,studyNames,recomputeSummary) {
  s<-fit$summary
  if(!is.null(provided)) {
    candidate<-provided;attr(candidate,'ctmaTables_fingerprint')<-NULL
    original<-s;attr(original,'ctmaTables_fingerprint')<-NULL
    if(!isTRUE(all.equal(candidate,original,check.attributes=TRUE))) .ctmaTablesStop('The supplied ctmaInit summary does not agree with the fitted object.')
  }
  input<-.ctmaTablesInitTables(fit,s);ns<-length(input$studies)
  if(!is.null(studyNames)) {
    if(!is.character(studyNames)||length(studyNames)!=ns||anyNA(studyNames)||any(!nzchar(studyNames))||anyDuplicated(studyNames)) .ctmaTablesStop('modelNames must uniquely label every primary study.')
    input$labels<-studyNames
  }
  d<-input$studies[[1L]]$table
  availableIDs<-paste0(d$matrix,'[',d$row,',',d$col,']')
  if(is.null(parameters)) parameters<-.ctmaTablesAsk(availableIDs,'parameters',interactive)
  selected<-.ctmaTablesSelect(d,parameters)
  availability<-input$studies[[1L]]$available
  common<-names(availability)[vapply(availability,function(x)all(x[selected]),logical(1))]
  if(is.null(statistics)) statistics<-.ctmaTablesAsk(common,'statistics',interactive)
  if(!is.null(driftDraws)&&(!is.list(driftDraws)||length(driftDraws)!=ns)) .ctmaTablesStop('For ctmaInit, driftDraws must be a list containing one joint drift array per primary study.')
  refreshed<-.ctmaTablesRefreshInit(fit,input,parameters,statistics,timeIntervals,driftDraws,dtPointOnly,recomputeSummary,interactive,digits)
  input<-refreshed$input
  if(isTRUE(refreshed$recomputed))common<-Reduce(intersect,lapply(input$studies,function(x) {
    z<-.ctmaTablesSelect(x$table,parameters)
    names(x$available)[vapply(x$available,function(v)all(v[z]),logical(1))]
  }))
  tables<-lapply(seq_len(ns),function(i) {
    primary<-structure(list(summary=list(estimates=input$studies[[i]]$table),argumentList=fit$argumentList,ctModel=fit$ctModel),class='CoTiMAFit')
    if(!is.null(refreshed$targets[[i]]))primary$studyFitList<-list(refreshed$targets[[i]])
    if(!is.null(refreshed$generated[[i]]))primary$ctma_tables_dt_summaries<-refreshed$generated[[i]]
    t<-.ctmaTablesSingle(primary,parameters=parameters,statistics=statistics,defaultOrder=defaultOrder,
      digits=digits,labels=labels,missing=missing,interactive=FALSE,format='none')
    sourceTable<-input$studies[[i]]$table
    z<-match(paste(t$data$matrix,t$data$row,t$data$col),paste(sourceTable$matrix,sourceTable$row,sourceTable$col))
    for(stat in setdiff(names(t$data),c('matrix','row','col','parameter'))) {
      retained<-input$studies[[i]]$available[[stat]]
      if(!is.null(retained)&&any(!retained[z])&&missing=='error') .ctmaTablesStop('Statistic ',stat,' was not retained for some selected ctmaInit parameters. Use missing="note" for blank cells.')
    }
    if(!is.null(timeIntervals)) t<-.ctmaTablesAddDT(t,primary,timeIntervals,driftDraws[[i]],dtPointOnly,digits,missing)
    t
  })
  canonical<-setdiff(names(tables[[1L]]$data),c('matrix','row','col','parameter'))
  for(i in seq_len(ns)) {
    old<-setdiff(names(tables[[i]]$data),c('matrix','row','col','parameter'))
    if(length(old)!=length(canonical)) .ctmaTablesStop('Primary-study statistics have inconsistent mappings.')
    names(tables[[i]]$data)[match(old,names(tables[[i]]$data))]<-canonical
    names(tables[[i]]$display)[-1L]<-canonical
  }
  combined<-lapply(seq_len(ns),function(i) {
    z<-tables[[i]]$data;z$Study<-input$labels[i];z$Source<-input$sources[i];z
  })
  data<-do.call(rbind,combined);rownames(data)<-NULL
  display<-do.call(rbind,lapply(seq_len(ns),function(i)cbind(Study=input$labels[i],tables[[i]]$display,stringsAsFactors=FALSE)))
  rownames(display)<-NULL
  noted<-s
  noted$uncertaintyNote<-if(isTRUE(refreshed$recomputed))'Requested missing statistics were obtained from saved native primary-study fits using summary(); their uncertainty interpretation is that of the native fits.' else 'Primary-study uncertainty columns labelled SE by ctmaInit originate from sd/Sd summaries and are reported here as SD. Aggregate summaries do not retain medians or Bayesian ESS/Rhat.'
  stats<-setdiff(names(tables[[1]]$data),c('matrix','row','col','parameter'))
  notes<-.ctmaTablesNotes(fit,noted,data,stats,timeUnit,tableNote)
  notes$notes<-c(notes$notes,if(isTRUE(refreshed$recomputed))'Primary-study estimates were regenerated from saved native fits; no models were refitted.' else 'Primary-study estimates come from the retained ctmaInit summary and saved coefficient lists; no models were refitted.',
    'Reported -2LL and parameter counts are totals across the separate primary-study fits, not statistics for a single pooled fit.')
  if(input$scale!=1&&!isTRUE(refreshed$recomputed)) notes$notes<-c(notes$notes,paste0('Original-time diffusion summary entries were mapped back to fitted units using the recorded multiplier ',input$scale,
    '. Retained fitted coefficients are used when available; recovered uncertainty is limited by the rounding of the stored summary.'))
  if(isTRUE(refreshed$recomputed))for(i in seq_len(ns)) {
    note<-refreshed$generated[[i]][[1L]]$uncertaintyNote
    if(length(note))notes$notes<-c(notes$notes,paste0(input$labels[i],': ',note))
  }
  if(any(!is.finite(as.matrix(data[,stats,drop=FALSE])))) notes$notes<-c(notes$notes,'Statistics not retained or not finite are shown as blank; no median or uncertainty statistic was reconstructed from a point estimate.')
  refs<-!is.na(input$sources)&nzchar(input$sources)&!input$sources %in% c('Reference not provided','NA')
  if(any(refs)) notes$notes<-c(notes$notes,paste0(input$labels[refs],': ',input$sources[refs]))
  if(!is.null(s$note)) notes$notes<-c(notes$notes,as.character(s$note))
  problems<-unlist(s$estimationProblems,use.names=TRUE)
  if(length(problems)) notes$notes<-c(notes$notes,paste0('Source estimation-problem notes: ',paste(problems,collapse='; ')))
  if(!is.null(timeIntervals)) notes$notes<-unique(c(notes$notes,tables[[1L]]$notes[grepl('DT |drift representation',tables[[1L]]$notes)]))
  validated<-if(is.null(provided))s else provided
  attr(validated,'ctmaTables_fingerprint')<-.ctmaTablesFingerprint(list(summary=s,means=fit$modelResults$T0VAR))
  notes$metadata$primary_studies<-input$labels;notes$metadata$study_sources<-input$sources
  notes$metadata$summary_validation<-'The full retained ctmaInit summary agrees with the supplied summary, when provided.'
  structure(list(data=data,display=display,notes=unique(notes$notes),files=character(),summary=validated,
    metadata=notes$metadata,generated_summaries=refreshed$generated,tables=stats::setNames(tables,input$labels),available_parameters=availableIDs,
    available_statistics=common),class='ctmaTables')
}
