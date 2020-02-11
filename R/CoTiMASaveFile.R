##### CoTiMA fitting procedure - singleStudyModelFit
CoTiMASaveFile <- function(activateRPB,
                           workingDirectory=Null,
                           SaveObject,
                           FileName,
                           Directory)
{
  if (!(dir.exists(Directory))){dir.create(Directory);print("Directory added!")}
  setwd(Directory)
  wdFileList <- list.files()        
  if (FileName %in% wdFileList) {
    if (activateRPB==TRUE) {pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    cat(red$bold(paste("The filename", FileName, "does already exist.", "\n")))
    cat(blue("Do you wish to override the existing file?", "\n"))
    cat(blue("Press 'y' to OVERRIDE the existing file or 'n' to continue WITHOUT SAVING.", "\n"))
    cat(blue("If you wish to CHANGE the filename please press 'c' or press 'q' to QUIT", "\n"))
    cat(blue("Press ENTER afterwards to confirm your choice."))
    char <- readline(" ")
    while (!(char == 'y') & !(char == 'Y') & !(char == 'n') & !(char == 'N') & !(char == 'c') & !(char == 'C') 
           & !(char == 'q') & !(char == 'Q')) {
      cat(blue("Please press 'y' to override, or 'n' to continue without saving. Press 'c' to change the filename or 'q to quit. Confirm your choice by pressing ENTER afterwards.", "\n"))
      char <- readline(" ")
    }
    if (char == 'y' | char == 'Y') {saveRDS(SaveObject, file=FileName)
    } else if (char == 'n' | char == 'N'){char
    } else if (char == 'c' | char == 'C'){
      cat(blue("Please enter new filename. Press ENTER afterwards to confirm your choice.", "\n"))
      chars <- readline(" ")
      FileName <- paste0(chars, ".rds"); FileName
      saveRDS(SaveObject, file=FileName)
    } else if (char == 'q' | char == 'Q') {stop("Good luck for the next try!")}
  } else {saveRDS(SaveObject, file=FileName)}
  setwd(workingDirectory)
}

