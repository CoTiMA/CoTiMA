# Presentation only: keep numeric values in their original long schema.
.ctmaTablesLayout <- function(result,isInit,multi,layout,sortBy,rowLayout) {
  d<-result$data;display<-result$display
  key<-paste(d$matrix,d$row,d$col,sep=':')
  if(multi) {
    keep<-!duplicated(key);catalog<-d[keep,,drop=FALSE]
  } else catalog<-d
  if(nrow(catalog)!=nrow(display)) .ctmaTablesStop('Display rows do not match parameter identities.')
  catalog$key<-paste(catalog$matrix,catalog$row,catalog$col,sep=':')
  if(isInit&&sortBy=='parameter') {
    parameters<-unique(catalog$key);studies<-unique(catalog$Study)
    order<-order(match(catalog$key,parameters),match(catalog$Study,studies))
    catalog<-catalog[order,,drop=FALSE];display<-display[order,,drop=FALSE]
    result$data<-d[order,,drop=FALSE];rownames(result$data)<-NULL
  }
  if(rowLayout!='parameter') {
    statCols<-setdiff(names(display),c('Study','Source','Parameter'))
    if(!length(statCols)) .ctmaTablesStop('A compact row layout requires at least one statistic.')
    parameterKeys<-unique(catalog$key)
    first<-match(parameterKeys,catalog$key)
    parameterLabels<-as.character(display$Parameter[first])
    duplicatedLabels<-duplicated(parameterLabels)|duplicated(parameterLabels,fromLast=TRUE)
    parameterLabels[duplicatedLabels]<-paste0(parameterLabels[duplicatedLabels],' [',parameterKeys[duplicatedLabels],']')
    group<-catalog$key
    if(rowLayout=='all')group[]<-'all'
    else {
      drift<-tolower(catalog$matrix) %in% c('drift','dtdrift')|grepl('^dtdrift@',tolower(catalog$matrix))
      group[drift]<-catalog$matrix[drift]
    }
    studies<-if(isInit)as.character(catalog$Study) else rep('',nrow(catalog))
    rows<-unique(data.frame(study=studies,group=group,stringsAsFactors=FALSE))
    titles<-vapply(seq_len(nrow(rows)),function(i) {
      z<-which(studies==rows$study[i]&group==rows$group[i])[1L]
      if(rowLayout=='all')'All selected parameters'
      else if(tolower(rows$group[i])=='drift')'DRIFT'
      else if(tolower(rows$group[i])=='dtdrift')'DT DRIFT'
      else if(grepl('^dtdrift@',tolower(rows$group[i])))paste0('DT DRIFT (lag ',sub('^[^@]+@','',rows$group[i]),')')
      else as.character(display$Parameter[z])
    },character(1))
    out<-data.frame(Parameter=titles,stringsAsFactors=FALSE,check.names=FALSE)
    if(isInit)out<-cbind(Study=rows$study,out,stringsAsFactors=FALSE)
    columns<-data.frame(key=names(out),group=names(out),statistic=names(out),stringsAsFactors=FALSE)
    for(p in seq_along(parameterKeys)) for(stat in statCols) {
      col<-paste(parameterLabels[p],stat,sep=' | ')
      if(col %in% names(out)) .ctmaTablesStop('Ambiguous compact-layout column: ',col,'. Use unique parameter labels.')
      out[[col]]<-rep('',nrow(out))
      z<-which(catalog$key==parameterKeys[p])
      dest<-match(paste(studies[z],group[z],sep='\r'),paste(rows$study,rows$group,sep='\r'))
      if(anyDuplicated(dest)) .ctmaTablesStop('Repeated parameter cells cannot be combined into one row.')
      out[[col]][dest]<-display[[stat]][z]
      columns<-rbind(columns,data.frame(key=col,group=parameterLabels[p],statistic=stat,stringsAsFactors=FALSE))
    }
    if(isInit&&length(unique(titles))==1L) {
      out$Parameter<-NULL
      columns<-columns[columns$key!='Parameter',,drop=FALSE]
    }
    display<-out
    attr(display,'ctmaTables_columns')<-columns
  }
  rownames(display)<-NULL
  result$display<-display
  result$layout<-list(page=layout,orientation=if(layout=='wide')'landscape' else 'portrait',
    sortBy=if(isInit)sortBy else NULL,rowLayout=rowLayout)
  result$panels<-.ctmaTablesPanels(display,layout)
  if(length(result$panels)>1L)result$notes<-unique(c(result$notes,
    paste0('The selected columns continue across ',length(result$panels),' labelled panels; no selected parameter or statistic was dropped.')))
  result
}
.ctmaTablesWidths <- function(display,layout) {
  text<-which(names(display) %in% c('Study','Source','Parameter'))
  if(!length(text))text<-1L
  numeric<-setdiff(seq_len(ncol(display)),text)
  total<-if(layout=='wide')10.1 else 6.2
  widths<-numeric(ncol(display))
  widths[text]<-if(length(text)==1L)1.9 else ifelse(names(display)[text]=='Study',1.1,1.8)
  available<-total-sum(widths[text])
  required<-if(length(numeric))vapply(numeric,function(j)max(.65,.078*max(nchar(display[[j]]),0L)+.10),numeric(1)) else numeric()
  if(length(numeric))widths[numeric]<-required+(available-sum(required))/length(numeric)
  else widths[text]<-total/length(text)
  list(text=text,numeric=numeric,widths=widths,total=total,required=required,available=available)
}
.ctmaTablesPanels <- function(display,layout) {
  widths<-.ctmaTablesWidths(display,layout)
  numeric<-widths$numeric
  required<-stats::setNames(widths$required,as.character(numeric))
  available<-widths$available
  if(any(required>available+1e-9)) .ctmaTablesStop('A numeric cell is too wide for this page. Use layout="wide" or fewer display digits.')
  if(sum(required)<=available+1e-9)return(list(display))
  mapping<-attr(display,'ctmaTables_columns')
  units<-if(is.null(mapping))as.list(numeric) else {
    group<-mapping$group[match(names(display)[numeric],mapping$key)]
    split(numeric,cumsum(c(TRUE,head(group,-1L)!=utils::tail(group,-1L))))
  }
  chunks<-list();current<-integer()
  cost<-function(cols)sum(required[as.character(cols)])
  for(unit in units) {
    if(cost(unit)>available+1e-9) {
      if(length(current)){chunks[[length(chunks)+1L]]<-current;current<-integer()}
      for(col in unit) {
        if(cost(c(current,col))>available+1e-9){chunks[[length(chunks)+1L]]<-current;current<-integer()}
        current<-c(current,col)
      }
      if(length(current)){chunks[[length(chunks)+1L]]<-current;current<-integer()}
    } else {
      if(cost(c(current,unit))>available+1e-9){chunks[[length(chunks)+1L]]<-current;current<-integer()}
      current<-c(current,unit)
    }
  }
  if(length(current))chunks[[length(chunks)+1L]]<-current

  lapply(chunks,function(cols) {
    x<-display[,sort(c(widths$text,cols)),drop=FALSE]
    if(!is.null(mapping))attr(x,'ctmaTables_columns')<-mapping[match(names(x),mapping$key),,drop=FALSE]
    x
  })
}
