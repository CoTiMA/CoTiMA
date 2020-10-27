##### CoTiMA fitting procedure - singleStudyModelFit

#' ctmaSaveFile
#'
#' @param activateRPB "set to TRUE to receive push messages with CoTiMA notifications on your phone"
#' @param activeDirectory "Directory name"
#' @param SaveObject "Object to save"
#' @param FileName "Filename"
#' @param Directory "Directory to save file in"
#' @param silentOverwrite "override old files without asking"
#'
<<<<<<< HEAD
#' @importFrom RPushbullet pbPost
#' @importFrom crayon red blue
=======
#' @import RPushbullet
#' @import crayon
>>>>>>> 04891f85cf4de7f719709e5405bff9d1d6a51695
#'
#' @export ctmaSaveFile
#'
ctmaSaveFile <- function(activateRPB,
                           #workingDirectory=Null,
                           activeDirectory=activeDirectory,
                           SaveObject,
                           FileName,
                           Directory,
                           silentOverwrite=FALSE)
{
  x1 <- paste0(activeDirectory, "/", Directory); x1
  #if (!(dir.exists(Directory))) {dir.create(Directory);print("Directory added!")}
  if (!(dir.exists(x1))) {dir.create(x1);print("Directory added!")}
  #setwd(Directory)
  wdFileList <- list.files(x1)
  if (FileName %in% wdFileList) {
    if (silentOverwrite==FALSE) {
      if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
      cat(crayon::red$bold(paste("The filename", FileName, "does already exist.", "\n")))
      cat(crayon::blue("Do you wish to override the existing file?", "\n"))
      cat(crayon::blue("Press 'y' to OVERRIDE the existing file or 'n' to continue WITHOUT SAVING.", "\n"))
      cat(crayon::blue("If you wish to CHANGE the filename please press 'c' or press 'q' to QUIT", "\n"))
      cat(crayon::blue("Press ENTER afterwards to confirm your choice."))
      char <- readline(" ")
      while (!(char == 'y') & !(char == 'Y') & !(char == 'n') & !(char == 'N') & !(char == 'c') & !(char == 'C')
             & !(char == 'q') & !(char == 'Q')) {
        cat(crayon::blue("Please press 'y' to override, or 'n' to continue without saving. Press 'c' to change the filename or 'q to quit. Confirm your choice by pressing ENTER afterwards.", "\n"))
        char <- readline(" ")
      }
      if (char == 'y' | char == 'Y') {
        FileName <- paste0(x1, FileName) ######
        saveRDS(SaveObject, file=FileName)
      } else if (char == 'n' | char == 'N'){char
      } else if (char == 'c' | char == 'C'){
        cat(crayon::blue("Please enter new filename. Press ENTER afterwards to confirm your choice.", "\n"))
        chars <- readline(" ")
        FileName <- paste0(chars, ".rds"); FileName
        FileName <- paste0(x1, FileName)
        saveRDS(SaveObject, file=FileName)
      } else if (char == 'q' | char == 'Q') {stop("Good luck for the next try!")}
    } else {
      FileName <- paste0(x1, FileName)
      saveRDS(SaveObject, file=FileName)
    }
  } else {
    FileName <- paste0(x1, FileName)
    saveRDS(SaveObject, file=FileName)
    }
  #setwd(workingDirectory)
}
