.ctmaTablesXml <- function(x) {
  x <- gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE)
  x<-gsub('>','&gt;',x,fixed=TRUE);x<-gsub('"','&quot;',x,fixed=TRUE)
  gsub("'",'&apos;',x,fixed=TRUE)
}
.ctmaTablesTex <- function(x) {
  replacements <- c('\\'='\\textbackslash{}','&'='\\&','%'='\\%','_'='\\_',
    '#'='\\#','$'='\\$','{'='\\{','}'='\\}','~'='\\textasciitilde{}','^'='\\textasciicircum{}')
  vapply(x,function(s) paste0(vapply(strsplit(s,'',fixed=TRUE)[[1]],function(ch)
    if(ch %in% names(replacements)) replacements[[ch]] else ch,character(1)),collapse=''),character(1),USE.NAMES=FALSE)
}
.ctmaTablesWrite <- function(display,notes,format,activeDirectory,file,title,overwrite=FALSE,layout="long",fontSize=10) {
  fail <- .ctmaTablesStop
  if(!is.data.frame(display)||!nrow(display)||!ncol(display)||anyNA(display)) fail('A nonempty formatted table without NA cells is required.')
  if(!is.character(format)||!length(format)||anyNA(format)||any(!format %in% c('docx','xml','tex','pdf'))) fail('Unsupported output format.')
  layout <- match.arg(layout,c("long","wide"))
  panels <- .ctmaTablesPanels(display,layout)
  format <- unique(format)
  if(!is.character(file)||length(file)!=1L||is.na(file)||!grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',file)||grepl('..',file,fixed=TRUE)) fail('file must be a safe filename base without a path or double dots.')
  if(!is.character(activeDirectory)||length(activeDirectory)!=1L||is.na(activeDirectory)||!dir.exists(activeDirectory)) fail('activeDirectory must be an existing directory.')
  if(!is.logical(overwrite)||length(overwrite)!=1L||is.na(overwrite)) fail('overwrite must be TRUE or FALSE.')
  text <- c(title,notes,unlist(display,use.names=FALSE),names(display))
  if(!is.character(title)||length(title)!=1L||is.na(title)||anyNA(text)||any(grepl('[\001-\010\013\014\016-\037]',text))) fail('Invalid title, note, or table text.')
  destination <- normalizePath(activeDirectory,mustWork=TRUE)
  paths <- stats::setNames(file.path(destination,paste0(file,'.',format)),format)
  if(any(dir.exists(paths))) fail('An output target is a directory.')
  if(!overwrite&&any(file.exists(paths))) fail('Output exists; choose another file or set overwrite=TRUE.')
  if('docx' %in% format) for(p in c('officer','flextable')) if(!requireNamespace(p,quietly=TRUE))
    fail('Word output requires package ',p,'. Install it with install.packages("',p,'").')
  compiler <- Sys.which('xelatex')
  if('pdf' %in% format && !nzchar(compiler)) fail('PDF export requires XeLaTeX. Install a TeX distribution (e.g. TinyTeX with XeLaTeX), or export tex/docx.')
  stage <- tempfile('.ctmatables-',tmpdir=destination);if(!dir.create(stage)) fail('Cannot create output staging directory.')
  keepStage <- FALSE
  on.exit(if(!keepStage)unlink(stage,recursive=TRUE),add=TRUE)
  staged <- stats::setNames(file.path(stage,paste0(file,'.',format)),format)
  if(any(c('tex','pdf') %in% format)) {
    tex <- .ctmaTablesRenderTex(panels,notes,title,layout,fontSize=fontSize)
    texfile <- file.path(stage,paste0(file,'.tex'));writeLines(tex,texfile,useBytes=TRUE)
    if('pdf' %in% format) {
      old <- setwd(stage);on.exit(setwd(old),add=TRUE)
      log <- file.path(stage,'compile.log')
      status <- system2(compiler,c('-no-shell-escape','-interaction=nonstopmode','-halt-on-error',shQuote(basename(texfile))),stdout=log,stderr=log)
      setwd(old)
      if(status!=0L||!file.exists(staged[['pdf']])) fail('PDF compilation failed. ',paste(utils::tail(readLines(log,warn=FALSE),8),collapse=' '))
    }
  }
  if('xml' %in% format) {
    xml <- .ctmaTablesRenderXml(panels,notes,title,layout,fontSize=fontSize)
    writeLines(xml,staged[['xml']],useBytes=TRUE)
  }
  if('docx' %in% format) {
    .ctmaTablesRenderDocx(panels,notes,title,layout,staged[['docx']],fontSize=fontSize)
  }
  # Publish only after all requested renderers succeed; restore replaced files on failure.
  committed <- character();backups <- character();success <- FALSE
  on.exit(if(!success){
    unlink(committed)
    failed <- character()
    for(p in names(backups)) if(!suppressWarnings(file.rename(backups[[p]],p))) failed<-c(failed,p)
    if(length(failed)) {
      keepStage<-TRUE
      warning('Could not restore prior outputs; backups are retained in ',stage,call.=FALSE)
    }
  },add=TRUE,after=FALSE)
  for(p in paths) if(file.exists(p)) {
    b <- file.path(stage,paste0('previous-',basename(p)))
    if(!file.rename(p,b)) fail('Cannot preserve existing output: ',p)
    backups[p] <- b
  }
  for(n in names(paths)) {
    if(!file.rename(staged[[n]],paths[[n]])) fail('Cannot publish output: ',paths[[n]])
    committed <- c(committed,paths[[n]])
  }
  success <- TRUE;paths
}

