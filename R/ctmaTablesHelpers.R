# Independently corrected prototype based on DeepSeek's initial core draft.
# Extraction, provenance, selection, uncertainty and DT calculation helpers.
.ctmaTablesStop <- function(...) stop(..., call. = FALSE)
.ctmaTablesParameterInput <- function(expression, caller) {
  families <- c('drift','diffusion','diffusioncov','t0var','t0cov','t0means',
    'cint','asymcint','asymdiffusioncov','lambda','manifestmeans','manifestvar',
    'manifestcov','pars','tdpredeffect','dtdrift')
  resolve <- function(x) {
    if(is.symbol(x) && tolower(as.character(x)) %in% families) return(as.character(x))
    if(is.call(x) && (identical(x[[1L]],quote(c)) || identical(x[[1L]],quote(base::c))))
      return(do.call(base::c,lapply(as.list(x)[-1L],resolve)))
    eval(x,envir=caller)
  }
  resolve(expression)
}
.ctmaTablesTable <- function(s) {
  if (!is.list(s)) .ctmaTablesStop('A summary must be a list with estimates or parmatrices.')
  d <- if (is.data.frame(s$parmatrices)) s$parmatrices else s$estimates
  if (!is.data.frame(d) || !all(c('row','col') %in% names(d)))
    .ctmaTablesStop('Unsupported summary: a matrix-coordinate parameter data frame is required; select a fitted model, not a power-analysis or multi-study summary.')
  if (!nrow(d)) .ctmaTablesStop('The summary contains no parameter estimates.')
  d
}
.ctmaTablesSignature <- function(s) {
  d <- .ctmaTablesTable(s)
  list(table=d, statistics=s[intersect(c('minus2ll','loglik','n.parameters','npars'),names(s))])
}
.ctmaTablesFingerprint <- function(x) {
  p <- tempfile();on.exit(unlink(p));saveRDS(x,p,version=2)
  unname(tools::md5sum(p))
}
.ctmaTablesNative <- function(fit) {
  if (inherits(fit,c('ctJuliaFit','ctStanFit','ctmaJuliaFit'))) return(fit)
  x <- fit$studyFitList
  if (inherits(x,c('ctJuliaFit','ctStanFit','ctmaJuliaFit'))) return(x)
  if (is.list(x) && length(x)==1L && inherits(x[[1]],c('ctJuliaFit','ctStanFit','ctmaJuliaFit'))) return(x[[1]])
  NULL
}
.ctmaTablesSummary <- function(fit,digits=8L) {
  if (is.list(fit$summary) && (is.data.frame(fit$summary$estimates) || is.data.frame(fit$summary$parmatrices))) return(fit$summary)
  if (is.list(fit$ctsem_summary) && is.data.frame(fit$ctsem_summary$parmatrices)) return(fit$ctsem_summary)
  if (is.list(fit$cotima_summary) && is.data.frame(fit$cotima_summary$estimates)) return(fit$cotima_summary)
  native <- .ctmaTablesNative(fit)
  if (is.null(native)) .ctmaTablesStop('No supported parameter summary is stored. Select one native ctsem fit explicitly for separately fitted primary studies.')
  args <- list(object=native)
  if (inherits(native,c('ctJuliaFit','ctStanFit'))) {
    method<-utils::getS3method('summary',class(native)[1L],optional=TRUE)
    fields<-if(is.null(method))character() else names(formals(method))
    if('digits' %in% fields)args$digits<-digits
    if('priorcheck' %in% fields)args$priorcheck<-FALSE
    if('residualcov' %in% fields)args$residualcov<-FALSE
    if(inherits(native,'ctStanFit'))native<-.ctmaTablesStanReady(native,'the fitted model')
    args$object<-native
  }
  invisible(utils::capture.output(s <- do.call(base::summary,args)))
  .ctmaTablesTable(s);s
}
.ctmaTablesParameters <- function(d,fit) {
  rn <- rownames(d)
  if (!'matrix' %in% names(d)) {
    if (is.null(rn) || all(rn==as.character(seq_len(nrow(d))))) .ctmaTablesStop('The summary lacks matrix identities.')
    d$matrix <- sub('[ _].*$','',rn)
  }
  d$matrix <- as.character(d$matrix)
  if (anyNA(d[,c('matrix','row','col')]) || !is.numeric(d$row) || !is.numeric(d$col) ||
      any(d$row<1 | d$col<1 | d$row!=floor(d$row) | d$col!=floor(d$col))) .ctmaTablesStop('Invalid parameter matrix coordinates.')
  if (anyDuplicated(paste(d$matrix,d$row,d$col))) .ctmaTablesStop('Ambiguous duplicate matrix coordinates in summary.')
  labels <- if ('parameter' %in% names(d)) as.character(d$parameter) else rep(NA_character_,nrow(d))
  pars <- fit$ctstanmodelbase$pars
  if (is.null(pars)) pars <- fit$modelbase$pars
  if (is.null(pars)) pars <- fit$ctModel$pars
  if (is.data.frame(pars) && all(c('matrix','row','col','param') %in% names(pars))) {
    key <- paste(tolower(pars$matrix),pars$row,pars$col)
    z <- match(paste(tolower(d$matrix),d$row,d$col),key)
    v <- as.character(pars$param[z]);ok <- !is.na(v)&nzchar(v)
    labels[ok & is.na(labels)] <- v[ok & is.na(labels)]
  }
  if (!is.null(rn)) {
    v <- sub('^[^ ]+ +','',rn)
    v <- sub(' +\\(invariant\\)$','',v)
    ok <- is.na(labels) & grepl(' ',rn,fixed=TRUE)
    labels[ok] <- v[ok]
  }
  fallback <- paste0(d$matrix,'[',d$row,',',d$col,']')
  labels[is.na(labels)|!nzchar(labels)] <- fallback[is.na(labels)|!nzchar(labels)]
  d$parameter <- labels;d
}
.ctmaTablesSelect <- function(d,parameters) {
  if (!is.character(parameters)||!length(parameters)||anyNA(parameters)||any(!nzchar(parameters))) .ctmaTablesStop('Select at least one parameter.')
  selected <- integer()
  family <- function(x) {
    alias <- c(t0var='t0cov',diffusion='diffusioncov',manifestvar='manifestcov')
    x <- tolower(x)
    if (x %in% names(alias) && alias[[x]] %in% tolower(d$matrix)) x <- alias[[x]]
    x
  }
  for (p in parameters) {
    group <- which(tolower(d$matrix)==family(p))
    if(length(group)) {selected<-c(selected,group);next}
    exact <- which(d$parameter==p)
    if (length(exact)>1L) .ctmaTablesStop('Ambiguous parameter label: ',p,'; use matrix coordinates.')
    if (length(exact)) {selected <- c(selected,exact);next}
    parts <- regmatches(p,regexec('^([A-Za-z][A-Za-z0-9_]*)\\[ *([1-9][0-9]*) *, *([1-9][0-9]*) *\\]$',p))[[1]]
    if (length(parts)) {
      m <- family(parts[2]);r <- as.numeric(parts[3]);cc <- as.numeric(parts[4])
      z <- which(tolower(d$matrix)==m & d$row==r & d$col==cc)
      if (!length(z) && m %in% c('t0cov','diffusioncov','manifestcov','asymdiffusioncov'))
        z <- which(tolower(d$matrix)==m & d$row==cc & d$col==r)
    } else z <- which(tolower(d$matrix)==family(p))
    if (!length(z)) .ctmaTablesStop('No match for parameter selector: ',p)
    selected <- c(selected,z)
  }
  unique(selected)
}
.ctmaTablesAsk <- function(items,kind,interactive) {
  if (!isTRUE(interactive)) .ctmaTablesStop('Non-interactive use requires explicit ',kind,'.')
  keep <- vapply(items,function(x) repeat {
    ans <- tolower(trimws(readline(paste0('Include ',kind,' ',x,'? [yes/no] '))))
    if (ans %in% c('y','yes')) return(TRUE)
    if (ans %in% c('n','no')) return(FALSE)
    if (!nzchar(ans)) .ctmaTablesStop('Selection cancelled; provide explicit ',kind,'.')
    message('Please answer yes or no.')
  },logical(1))
  items[keep]
}
.ctmaTablesStats <- function(requested,available) {
  if (!is.character(requested)||!length(requested)||anyNA(requested)||any(!nzchar(requested))) .ctmaTablesStop('Select at least one statistic.')
  norm <- function(x) sub('^n[.]eff$','neff',tolower(gsub('[%_ ]','',x)))
  groups <- list(mean=c('estimate','mean'),sd=c('sd'),se=c('se','stderr','standarderror'),
    `50`=c('median','50'),`2.5`=c('2.5'),`97.5`=c('97.5'),rhat='rhat',
    neff=c('ess','neff'),bulkess=c('bulkess','essbulk'),tailess=c('tailess','esstail'))
  # Quantiles are matched before decimal punctuation is removed.
  canonical <- function(x) {
    x <- tolower(gsub('%','',x,fixed=TRUE))
    for(g in names(groups)) if(x %in% groups[[g]] || norm(x) %in% norm(groups[[g]])) return(g)
    norm(x)
  }
  vapply(requested,function(x) {
    z <- which(vapply(available,canonical,character(1))==canonical(x))
    if(length(z)>1L) .ctmaTablesStop('Ambiguous statistic: ',x)
    if(length(z)) available[z] else x
  },character(1),USE.NAMES=FALSE)
}
.ctmaTablesMethod <- function(fit) {
  a <- fit$argumentList
  if (is.null(a$optimize)) a <- fit$args$resolved
  if (is.null(a$optimize)) a <- fit$args
  if (is.null(a$optimize)) a <- .ctmaTablesNative(fit)$args$resolved
  if (is.null(a$optimize)) a <- .ctmaTablesNative(fit)$args
  if (isFALSE(a$optimize)) return('Bayesian')
  if (isTRUE(a$optimize) && isTRUE(a$priors)) return('MAP')
  if (isTRUE(a$optimize) && isFALSE(a$priors)) return('ML')
  'unspecified'
}
.ctmaTablesScalar <- function(...) {
  for(x in list(...)) if(is.numeric(x)&&length(x)==1L&&is.finite(x)) return(as.numeric(x))
  NA_real_
}
.ctmaTablesNotes <- function(fit,s,d,stats,timeUnit,tableNote) {
  method <- .ctmaTablesMethod(fit)
  ll <- .ctmaTablesScalar(s$minus2ll,fit$minus2ll)
  if(is.na(ll)) {v <- .ctmaTablesScalar(s$loglik,fit$loglik);if(!is.na(v)) ll <- -2*v}
  np <- .ctmaTablesScalar(s$n.parameters,s$npars,fit$n.parameters,fit$npars)
  notes <- paste0('Estimation: ',method,'.')
  if(is.finite(np)) notes <- c(notes,paste0('Number of estimated parameters = ',np,'.'))
  if(is.finite(ll)) notes <- c(notes,paste0('-2LL = ',format(ll,digits=12),' (minus twice the reported log likelihood).'))
  else notes <- c(notes,'-2LL is unavailable; no posterior objective was substituted for the log likelihood.')
  interval <- any(grepl('%|median|^2[.]5$|^97[.]5$|^50$',stats,ignore.case=TRUE))
  if(interval) notes <- c(notes,if(method=='Bayesian') 'CrI = credible interval; quantiles summarize posterior draws.' else if(method=='ML') 'CI = confidence interval; quantiles summarize the reported uncertainty distribution.' else if(method=='MAP') 'MAP quantiles summarize the stored uncertainty approximation; no frequentist confidence-interval interpretation is assumed.' else 'Interval interpretation is unspecified by the fit metadata.')
  if(any(grepl('ess|n_eff',stats,ignore.case=TRUE))) notes <- c(notes,'ESS = effective sample size. Legacy n_eff, bulk ESS and tail ESS retain their source definitions.')
  if(any(grepl('rhat',stats,ignore.case=TRUE))) notes <- c(notes,'Rhat = the source-reported between-chain convergence diagnostic.')
  if(any(tolower(d$matrix)=='drift')) notes <- c(notes,'CT = continuous time. Matrix row denotes the outcome; column denotes the predictor.')
  glossary <- c(t0means='initial latent means',t0cov='initial latent covariance',
    cint='continuous-time intercepts',asymcint='long-term equilibrium positions',
    lambda='measurement loadings',manifestmeans='manifest intercepts',
    manifestcov='manifest measurement-error covariance',
    diffusioncov='continuous-time process-noise covariance',
    asymdiffusioncov='asymptotic latent covariance')
  for(m in unique(d$matrix)) if(tolower(m) %in% names(glossary))
    notes <- c(notes,paste0(m,' = ',glossary[[tolower(m)]],'.'))
  u <- s$uncertaintyNote
  if(is.character(u)&&length(u)) notes <- c(notes,u) else notes <- c(notes,'Uncertainty source: unspecified by the stored summary; SD was not relabelled as an analytic SE.')
  notes <- c(notes,if(is.null(timeUnit)) 'Time unit: fitted model units; no automatic rescaling was performed.' else paste0('Time unit label: ',timeUnit,'. No automatic rescaling was performed.'))
  if(is.numeric(s$scaledTime)&&length(s$scaledTime)==1L) notes <- c(notes,paste0('Stored model time multiplier = ',s$scaledTime,'.'))
  if(isFALSE(fit$optim$converged)||isFALSE(fit$optimization$converged)) notes <- c(notes,'WARNING: the fit explicitly reports non-convergence.')
  if(isFALSE(.ctmaTablesNative(fit)$optim$converged)) notes <- c(notes,'WARNING: the native fit explicitly reports non-convergence.')
  list(notes=unique(c(notes,tableNote)),metadata=list(method=method,minus2ll=ll,nparameters=np))
}
.ctmaTablesOrder <- function(d) {
  m <- tolower(d$matrix)
  groups <- list(c('lambda','manifestcov'),c('t0means','t0cov'),c('manifestmeans','cint','asymcint'),
                 'drift','diffusioncov',c('asymdiffusioncov'))
  g <- rep(99L,nrow(d));for(i in seq_along(groups)) g[m %in% groups[[i]]] <- i
  order(g,ifelse(m=='drift',d$row==d$col,FALSE),d$matrix,d$row,d$col)
}
.ctmaTablesSingle <- function(fit, summary=NULL, parameters=NULL, statistics=NULL,
    activeDirectory=NULL, defaultOrder=TRUE, format='docx', file='ctma_table',
    title='Model parameter estimates', digits=3, labels=NULL, missing=c('error','note'),
    overwrite=FALSE, timeUnit=NULL, tableNote=NULL, interactive=base::interactive()) {
  missing <- match.arg(missing)
  if(!is.list(fit)) .ctmaTablesStop('fit must be a CoTiMA or ctsem fitted object.')
  if(!is.logical(defaultOrder)||length(defaultOrder)!=1L||is.na(defaultOrder)) .ctmaTablesStop('defaultOrder must be TRUE or FALSE.')
  if(!is.numeric(digits)||length(digits)!=1L||is.na(digits)||digits<0||digits>12||digits!=floor(digits)) .ctmaTablesStop('digits must be an integer from 0 to 12.')
  if(!is.character(format)||!length(format)||anyNA(format)||any(!format %in% c('none','docx','xml','pdf','tex'))||('none' %in% format && length(format)>1L)) .ctmaTablesStop('Unsupported output format.')
  stored <- .ctmaTablesSummary(fit)
  if(!is.null(summary) && !isTRUE(all.equal(.ctmaTablesSignature(summary),.ctmaTablesSignature(stored),check.attributes=TRUE)))
    .ctmaTablesStop('The supplied summary does not agree with the fit. Use its stored summary or regenerate it with matching summary options.')
  validated <- if(is.null(summary)) stored else summary
  attr(validated,'ctmaTables_fingerprint') <- .ctmaTablesFingerprint(.ctmaTablesSignature(stored))
  native <- .ctmaTablesNative(fit)
  d <- .ctmaTablesParameters(.ctmaTablesTable(stored),if(is.null(native))fit else native)
  available <- setdiff(names(d)[vapply(d,is.numeric,logical(1))],c('row','col'))
  ids <- paste0(d$matrix,'[',d$row,',',d$col,']')
  if(is.null(parameters)) parameters <- .ctmaTablesAsk(ids,'parameters',interactive)
  z <- .ctmaTablesSelect(d,parameters);d <- d[z,,drop=FALSE]
  if(defaultOrder) d <- d[.ctmaTablesOrder(d),,drop=FALSE]
  if(is.null(statistics)) statistics <- .ctmaTablesAsk(available,'statistics',interactive)
  stats <- unique(.ctmaTablesStats(statistics,available));bad <- setdiff(stats,available)
  if(length(bad)&&missing=='error') .ctmaTablesStop('Unavailable statistics: ',paste(bad,collapse=', '))
  out <- d[,c('matrix','row','col','parameter',intersect(stats,available)),drop=FALSE]
  for(n in bad) out[[n]] <- NA_real_
  out <- out[,c('matrix','row','col','parameter',stats),drop=FALSE]
  notes <- .ctmaTablesNotes(fit,stored,d,stats,timeUnit,tableNote)
  if(length(bad)) notes$notes <- c(notes$notes,paste0('Unavailable statistics shown as blank: ',paste(bad,collapse=', '),'.'))
  display <- data.frame(Parameter=out$parameter,stringsAsFactors=FALSE)
  if(!is.null(labels)) {
    if(!is.character(labels)||is.null(names(labels))||anyDuplicated(names(labels))||anyNA(labels)) .ctmaTablesStop('labels must be a uniquely named character vector.')
    id <- paste0(out$matrix,'[',out$row,',',out$col,']')
    for(n in names(labels)) {j<-which(out$parameter==n|tolower(id)==tolower(n));if(!length(j)) .ctmaTablesStop('No selected parameter matches label: ',n);display$Parameter[j]<-labels[[n]]}
  }
  for(n in stats) display[[n]] <- ifelse(!is.finite(out[[n]]),'',formatC(out[[n]],format='f',digits=digits))
  if(any(!is.finite(as.matrix(out[,stats,drop=FALSE])))) notes$notes<-c(notes$notes,'Nonfinite or unavailable statistics are shown as blank; numeric source values are retained in the returned data.')
  notes$metadata$summary_validation <- 'Full parameter table and explicit fit statistics agree with the stored or regenerated summary; historical origin of an untagged legacy summary is not provable.'
  files <- character()
  if(!identical(format,'none')) {
    if(is.null(activeDirectory)) activeDirectory <- if(is.null(fit$activeDirectory))getwd() else fit$activeDirectory
    files <- .ctmaTablesWrite(display,notes$notes,format,activeDirectory,file,title,overwrite)
  }
  structure(list(data=out,display=display,notes=notes$notes,files=files,summary=validated,
    metadata=notes$metadata,available_parameters=ids,available_statistics=available),class='ctmaTables')
}

.ctmaTablesDriftDraws <- function(fit,draws=NULL) {
  if(!is.null(draws)) {
    if(!is.array(draws)||length(dim(draws))!=3L||dim(draws)[1]<2L||dim(draws)[2]!=dim(draws)[3]||!is.numeric(draws)||any(!is.finite(draws))) .ctmaTablesStop('driftDraws must be finite joint transformed drift draws with dimensions draws x outcomes x predictors.')
    return(draws)
  }
  native <- .ctmaTablesNative(fit)
  stan <- .ctmaTablesStanDriftDraws(native)
  if(!is.null(stan))return(stan)
  tp <- native$transformedpars
  if(is.list(tp$layout)&&is.matrix(tp$flat)) {
    k <- which(tp$layout$matrix=='DRIFT')
    if(length(k)==1L) {
      nr <- tp$layout$nrow[k];nc <- tp$layout$ncol[k];offset <- tp$layout$offset[k]
      if(nr==nc&&ncol(tp$flat)>=2L&&offset+nr*nc<=nrow(tp$flat))
        return(array(t(tp$flat[offset+seq_len(nr*nc),,drop=FALSE]),dim=c(ncol(tp$flat),nr,nc)))
    }
  }
  NULL
}
.ctmaTablesAddDT <- function(result,fit,lags,draws,pointOnly,digits,missing) {
  if(!is.numeric(lags)||!length(lags)||anyNA(lags)||any(!is.finite(lags))||any(lags<=0)||anyDuplicated(lags)) .ctmaTablesStop('timeIntervals must contain distinct finite positive lags in fitted model time units.')
  native <- .ctmaTablesNative(fit)
  if(isFALSE(native$ctstanmodel$continuoustime)||isFALSE(native$model$continuoustime))
    .ctmaTablesStop('DT interval conversion by a matrix exponential requires a continuous-time model; native discrete-time transitions need a different conversion.')
  pars <- native$model$pars
  if(is.null(pars)) pars <- fit$ctModel$pars
  if(is.data.frame(pars)&&all(c('matrix','stated') %in% names(pars))&&any(pars$matrix=='DRIFT'&pars$stated %in% TRUE))
    .ctmaTablesStop('DT conversion requires a constant drift matrix; state-dependent drift needs model-specific calculations.')
  cached <- fit$ctma_tables_dt_summaries
  if(!is.null(cached)&&is.null(draws)&&!isTRUE(pointOnly)) {
    for(lag in lags) {
      key<-format(lag,digits=15)
      if(is.null(cached[[key]])) .ctmaTablesStop('No regenerated DT summary for lag ',lag,'.')
      result<-.ctmaTablesDTFromSummary(result,cached[[key]],lag,digits,missing)
    }
    return(result)
  }
  if(!requireNamespace('expm',quietly=TRUE)) .ctmaTablesStop('DT effects require package expm; install.packages("expm").')
  draws <- if(isTRUE(pointOnly))NULL else .ctmaTablesDriftDraws(fit,draws)
  if(is.null(draws)&&!isTRUE(pointOnly)) .ctmaTablesStop('DT uncertainty needs joint transformed drift draws or a native summary with uncertainty. If saved native fits are available, set recomputeSummary="yes" to transform their draws through summary(timeinterval=...). Otherwise supply driftDraws (draws x outcomes x predictors), or explicitly request dtPointOnly=TRUE. A raw point estimate alone is insufficient.')
  d <- result$data;sel <- which(tolower(d$matrix)=='drift')
  if(!length(sel)) .ctmaTablesStop('Select at least one DRIFT entry when requesting timeIntervals.')
  full <- .ctmaTablesParameters(.ctmaTablesTable(.ctmaTablesSummary(fit)),fit)
  a <- full[tolower(full$matrix)=='drift',,drop=FALSE]
  meanCol <- .ctmaTablesStats('estimate',names(a))
  if(!meanCol %in% names(a)) .ctmaTablesStop('No drift mean is available for DT point estimates.')
  n <- max(c(a$row,a$col))
  if(nrow(a)!=n*n||anyDuplicated(paste(a$row,a$col))) .ctmaTablesStop('A complete constant drift matrix is required for DT conversion.')
  A <- matrix(NA_real_,n,n);A[cbind(a$row,a$col)]<-a[[meanCol]]
  if(any(!is.finite(A))) .ctmaTablesStop('Nonfinite drift estimates cannot be converted to DT effects.')
  if(!is.null(draws)&&any(dim(draws)[2:3]!=n)) .ctmaTablesStop('driftDraws dimensions do not match the fitted drift matrix.')
  stats <- setdiff(names(d),c('matrix','row','col','parameter'))
  added <- list();note <- character()
  for(lag in lags) {
    joint <- if(!is.null(draws)) vapply(seq_len(dim(draws)[1]),function(i)as.vector(expm::expm(matrix(draws[i,,],n,n)*lag)),numeric(n*n)) else NULL
    point <- as.vector(expm::expm(A*lag))
    dt <- d[sel,,drop=FALSE];dt$matrix <- paste0('dtDRIFT@',format(lag,digits=15))
    dt$parameter <- paste0('DT ',d$parameter[sel],' (lag ',format(lag,digits=15),')')
    for(s in stats) {
      values <- rep(NA_real_,nrow(dt));kind <- tolower(gsub('%','',s,fixed=TRUE))
      for(j in seq_len(nrow(dt))) {
        k <- (dt$col[j]-1L)*n+dt$row[j]
        if(isTRUE(pointOnly)) {if(kind %in% c('mean','estimate')) values[j]<-point[k]}
        else if(kind %in% c('mean','estimate')) values[j]<-mean(joint[k,])
        else if(kind=='sd') values[j]<-stats::sd(joint[k,])
        else if(kind=='median') values[j]<-stats::median(joint[k,])
        else if(grepl('^[0-9]+([.][0-9]+)?$',kind)&&as.numeric(kind)>0&&as.numeric(kind)<100)
          values[j]<-unname(stats::quantile(joint[k,],as.numeric(kind)/100))
      }
      if(all(is.na(values))) {
        if(missing=='error'&&!isTRUE(pointOnly)) .ctmaTablesStop('Statistic ',s,' is unavailable for DT effects; choose missing="note" or omit this statistic.')
        note <- c(note,paste0(s,' is unavailable for DT effects and is shown blank.'))
      }
      dt[[s]] <- values
    }
    added[[length(added)+1L]] <- dt
  }
  extra <- do.call(rbind,added);result$data <- rbind(result$data,extra)
  dtLabels <- unlist(lapply(lags,function(lag) paste0('DT ',result$display$Parameter[sel],
    ' (lag ',format(lag,digits=15),')')),use.names=FALSE)
  display <- data.frame(Parameter=dtLabels,stringsAsFactors=FALSE)
  for(s in stats) display[[s]]<-ifelse(is.na(extra[[s]]),'',formatC(extra[[s]],format='f',digits=digits))
  result$display<-rbind(result$display,display)
  result$notes<-unique(c(result$notes,'DT = discrete time. Reported lags are in fitted model time units; no time rescaling was performed.',
    if(isTRUE(pointOnly)) 'DT estimates are plug-in matrix exponentials of the mean drift; DT intervals are not available.' else 'DT summaries transform each joint drift draw by the matrix exponential, then summarize the resulting distribution.',
    'DT effects describe the selected drift representation; they are not averages over subject- or study-specific trajectories.',note))
  result$metadata$timeIntervals<-lags;result
}
