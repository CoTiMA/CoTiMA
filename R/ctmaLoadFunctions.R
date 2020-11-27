{

  print("You will need to install the following packages using \"install.packages(packagename)\" ")
  print("install.packages(ctsem)")
  print("install.packages(ctsemOMX)")
  print("install.packages(OpenMx)")
  print("install.packages(MASS)")
  print("install.packages(MBESS)")
  print("install.packages(rootSolve)")
  print("install.packages(doParallel)")
  print("install.packages(crayon)")
  print("install.packages(psych)")
  print("install.packages(RPushbullet)")

  mainPath <- "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/GitHub/CoTiMA/R/"

  source(paste0(mainPath, "ctmaInit.R"))
  source(paste0(mainPath, "ctmaBias.R"))
  source(paste0(mainPath, "ctmaPlot.R"))
  source(paste0(mainPath, "ctmaFit.R"))
  source(paste0(mainPath, "ctmaEqual.R"))
  source(paste0(mainPath, "ctmaSummary.R"))                     # called by generic "summary" function
  source(paste0(mainPath, "ctmaCompFit.R"))

  source(paste0(mainPath, "ctmaEmpCov.R"))                      # called by User (Reorganizes Matrix reported in Articles)
  source(paste0(mainPath, "ctmaPrep.R"))                        # called by User
  source(paste0(mainPath, "ctmaFittingSSMF.R"))                 # called by various functions
  source(paste0(mainPath, "ctmaPRaw.R"))                        # called by ctmaPrep.R
  source(paste0(mainPath, "ctmaCombPRaw.R"))                    # called by ctmaPrep.R
  source(paste0(mainPath, "ctmaSaveFile.R"))                    # called by various functions
  source(paste0(mainPath, "ctmaPower.R"))
  source(paste0(mainPath, "ctmaModFull.R"))
  source(paste0(mainPath, "ctmaModMultiple.R"))
  source(paste0(mainPath, "ctmaStanResample.R"))
  source(paste0(mainPath, "ctmaStanctArgs.R"))

  print("If you see warnings, check if the following packages are installed:")
  print("ctsem, MASS, MBESS, rootSolve, doParallel, crayon, psych, RPushbullet.")

}
