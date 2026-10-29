.ctmaTablesTexWrap <- function(x) {
  vapply(x,function(text) {
    words<-strsplit(text,' ',fixed=TRUE)[[1L]]
    paste(vapply(words,function(word) {
      if(nchar(word)<=18L)return(.ctmaTablesTex(word))
      starts<-seq.int(1L,nchar(word),by=10L)
      paste(.ctmaTablesTex(substring(word,starts,pmin(starts+9L,nchar(word)))),collapse='\\allowbreak{}')
    },character(1)),collapse=' ')
  },character(1),USE.NAMES=FALSE)
}
# Shared page width, grouped headers and column panels across output formats.
.ctmaTablesPanelTitle <- function(title,i,n) {
  if(n==1L)title else paste0(title,' (column panel ',i,' of ',n,')')
}
.ctmaTablesRenderTex <- function(panels,notes,title,layout,fontSize=10) {
  text<-c('\\documentclass[12pt,a4paper]{article}',
    paste0('\\usepackage[margin=0.7in',if(layout=='wide')',landscape' else '',']{geometry}'),
    '\\usepackage{fontspec,booktabs,longtable,array}',
    '\\IfFontExistsTF{Times New Roman}{\\setmainfont{Times New Roman}}{}',
    '\\setlength{\\tabcolsep}{3pt}','\\begin{document}',
    paste0('\\fontsize{',fontSize,'}{',1.2*fontSize,'}\\selectfont'))
  for(i in seq_along(panels)) {
    d<-panels[[i]];w<-.ctmaTablesWidths(d,layout);mapping<-attr(d,'ctmaTables_columns')
    # Include intercolumn padding in the width budget rather than overflowing it.
    widths<-w$widths*(w$total-2*ncol(d)*3/72.27)/w$total
    spec<-vapply(seq_len(ncol(d)),function(j)paste0('>{',if(j %in% w$text)'\\raggedright' else '\\raggedleft',
      '\\arraybackslash}p{',formatC(widths[j],format='f',digits=4),'in}'),character(1))
    if(i>1L)text<-c(text,'\\clearpage')
    text<-c(text,'\\noindent\\textbf{Table 1}\\par',paste0('\\noindent\\textit{',
      .ctmaTablesTex(.ctmaTablesPanelTitle(title,i,length(panels))),'}\\par'),
      paste0('\\begin{longtable}{',paste0(spec,collapse=''),'}'),'\\toprule')
    if(!is.null(mapping)) {
      groups<-mapping$group;starts<-c(1L,which(head(groups,-1L)!=tail(groups,-1L))+1L)
      ends<-c(starts[-1L]-1L,length(groups))
      cells<-vapply(seq_along(starts),function(j) {
        span<-ends[j]-starts[j]+1L
        width<-sum(widths[starts[j]:ends[j]])+2*(span-1L)*3/72.27
        paste0('\\multicolumn{',span,'}{>{\\centering\\arraybackslash}p{',formatC(width,format='f',digits=4),
          'in}}{',.ctmaTablesTexWrap(groups[starts[j]]),'}')
      },character(1))
      text<-c(text,paste0(paste(cells,collapse=' & '),' \\\\'))
      header<-mapping$statistic;header[w$text]<-''
    } else header<-names(d)
    rows<-apply(d,1,function(x) {
      cells<-.ctmaTablesTexWrap(x)
      cells[w$numeric]<-paste0('\\mbox{',.ctmaTablesTex(x[w$numeric]),'}')
      paste(cells,collapse=' & ')
    })
    text<-c(text,paste0(paste(.ctmaTablesTexWrap(header),collapse=' & '),' \\\\'),
      '\\midrule','\\endhead',paste0(rows,' \\\\'),'\\bottomrule','\\end{longtable}')
  }
  c(text,paste0('\\noindent\\textit{Note.} ',paste(.ctmaTablesTex(notes),collapse=' ')),'\\end{document}')
}
.ctmaTablesRenderXml <- function(panels,notes,title,layout,fontSize=10) {
  paragraph<-function(x)paste0('<w:p><w:r><w:rPr><w:sz w:val="',round(2*fontSize),'"/><w:szCs w:val="',round(2*fontSize),'"/></w:rPr><w:t xml:space="preserve">',.ctmaTablesXml(x),'</w:t></w:r></w:p>')
  cell<-function(x,width,span=1L)paste0('<w:tc><w:tcPr><w:tcW w:w="',round(width*1440),
    '" w:type="dxa"/>',if(span>1L)paste0('<w:gridSpan w:val="',span,'"/>') else '',
    '</w:tcPr>',paragraph(x),'</w:tc>')
  row<-function(x,widths,header=FALSE)paste0('<w:tr>',if(header)'<w:trPr><w:tblHeader/></w:trPr>' else '',
    paste0(vapply(seq_along(x),function(j)cell(x[j],widths[j]),character(1)),collapse=''),'</w:tr>')
  text<-c('<?xml version="1.0" encoding="UTF-8"?>','<?mso-application progid="Word.Document"?>',
    '<w:wordDocument xmlns:w="http://schemas.microsoft.com/office/word/2003/wordml"><w:body>')
  for(i in seq_along(panels)) {
    d<-panels[[i]];w<-.ctmaTablesWidths(d,layout);mapping<-attr(d,'ctmaTables_columns')
    if(i>1L)text<-c(text,'<w:p><w:r><w:br w:type="page"/></w:r></w:p>')
    text<-c(text,paragraph('Table 1'),paragraph(.ctmaTablesPanelTitle(title,i,length(panels))),
      '<w:tbl><w:tblPr><w:tblBorders><w:top w:val="single" w:sz="6"/><w:bottom w:val="single" w:sz="6"/><w:left w:val="nil"/><w:right w:val="nil"/><w:insideH w:val="nil"/><w:insideV w:val="nil"/></w:tblBorders></w:tblPr>',
      paste0('<w:tblGrid>',paste0('<w:gridCol w:w="',round(w$widths*1440),'"/>',collapse=''),'</w:tblGrid>'))
    if(!is.null(mapping)) {
      groups<-mapping$group;starts<-c(1L,which(head(groups,-1L)!=tail(groups,-1L))+1L);ends<-c(starts[-1L]-1L,length(groups))
      cells<-vapply(seq_along(starts),function(j)cell(groups[starts[j]],sum(w$widths[starts[j]:ends[j]]),ends[j]-starts[j]+1L),character(1))
      text<-c(text,paste0('<w:tr><w:trPr><w:tblHeader/></w:trPr>',paste0(cells,collapse=''),'</w:tr>'))
      header<-mapping$statistic;header[w$text]<-''
    } else header<-names(d)
    text<-c(text,row(header,w$widths,TRUE),apply(d,1,row,widths=w$widths),'</w:tbl>')
  }
  size<-if(layout=='wide')'<w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>' else '<w:pgSz w:w="11906" w:h="16838"/>'
  c(text,paragraph(paste('Note.',paste(notes,collapse=' '))),paste0('<w:sectPr>',size,
    '<w:pgMar w:top="1008" w:right="1008" w:bottom="1008" w:left="1008"/></w:sectPr>'),'</w:body></w:wordDocument>')
}
.ctmaTablesRenderDocx <- function(panels,notes,title,layout,path,fontSize=10) {
  doc<-officer::read_docx()
  doc<-officer::body_set_default_section(doc,officer::prop_section(
    page_size=officer::page_size(orient=if(layout=='wide')'landscape' else 'portrait'),
    page_margins=officer::page_mar(top=.7,bottom=.7,left=.7,right=.7)))
  font<-officer::fp_text(font.family='Times New Roman',font.size=fontSize,color='black')
  doc<-officer::docx_set_paragraph_style(doc,style_id='Title',style_name='Title',fp_t=font)
  for(i in seq_along(panels)) {
    d<-panels[[i]];w<-.ctmaTablesWidths(d,layout);mapping<-attr(d,'ctmaTables_columns')
    if(i>1L)doc<-officer::body_add_break(doc)
    doc<-officer::body_add_fpar(doc,officer::fpar(officer::ftext('Table 1',officer::fp_text(font.family='Times New Roman',font.size=fontSize,bold=TRUE))),style='Title')
    doc<-officer::body_add_fpar(doc,officer::fpar(officer::ftext(.ctmaTablesPanelTitle(title,i,length(panels)),officer::fp_text(font.family='Times New Roman',font.size=fontSize,italic=TRUE))))
    ft<-flextable::flextable(d)
    if(!is.null(mapping)) {
      headers<-data.frame(col_keys=mapping$key,Parameter=mapping$group,Statistic=mapping$statistic,stringsAsFactors=FALSE,check.names=FALSE)
      ft<-flextable::set_header_df(ft,mapping=headers,key='col_keys')
      ft<-flextable::merge_h(ft,i=1,part='header')
      ft<-flextable::merge_v(ft,j=w$text,part='header')
    }
    ft<-flextable::theme_booktabs(ft)
    ft<-flextable::font(ft,fontname='Times New Roman',part='all')
    ft<-flextable::fontsize(ft,size=fontSize,part='all')
    ft<-flextable::align(ft,j=w$text,align='left',part='all')
    if(length(w$numeric))ft<-flextable::align(ft,j=w$numeric,align='right',part='body')
    if(length(w$numeric))ft<-flextable::align(ft,j=w$numeric,align='center',part='header')
    ft<-flextable::width(ft,j=seq_len(ncol(d)),width=w$widths)
    ft<-flextable::padding(ft,padding.top=2,padding.bottom=2,padding.left=3,padding.right=3,part='all')
    ft<-flextable::line_spacing(ft,space=1,part='all')
    ft<-flextable::set_table_properties(ft,layout='fixed',width=1,align='left',opts_word=list(split=FALSE))
    ft<-flextable::paginate(ft,init=FALSE,hdr_ftr=TRUE)
    doc<-flextable::body_add_flextable(doc,ft)
  }
  doc<-officer::body_add_fpar(doc,officer::fpar(officer::ftext(paste('Note.',paste(notes,collapse=' ')),font)))
  print(doc,target=path)
}
