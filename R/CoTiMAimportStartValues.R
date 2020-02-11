importStartValues <- function(filePrefix=NULL, modelName=NULL, modelNumber=NULL, directory=NULL) {
  if ("OpenMx" %in% (.packages())) detach("package:OpenMx", force=TRUE) #, unload=TRUE)
  if ("ctsem" %in% (.packages())) detach("package:ctsem", force=TRUE) #, unload=TRUE)
  Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
  library(ctsem)
  mxOption(key='Number of Threads', value=parallel::detectCores()) #now
  
  modelFit <- readRDS(paste0(directory, "/", filePrefix, " ", modelName, modelNumber, ".rds"))
  #print(modelFit$mxobj$output$estimate)
  return(modelFit$mxobj@output$estimate) # @ instead of $ introduced 
}
