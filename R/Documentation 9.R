#library('ctsem')
#source("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/CoTiMA 2.1.0.7.R")
#source("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/CoTiMAcompileListOfPrimaryStudies.R")
#source("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/CoTiMAimportStartValues.R")
#source("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/CoTiMAfitUserSpecifiedModel.R")
#source("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/CoTiMAcombinePseudoRawData.R")
#source("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/CoTiMAcomputeStartValues.R")
#setwd("/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/")

################################ Primary Study 1 ###############################
empcov1 <- matrix(c(1.00, 0.47, 0.61, 0.34,
                    0.47, 1.00, 0.34, 0.69,
                    0.61, 0.34, 1.00, 0.44,
                    0.34, 0.69, 0.44, 1.00), nrow=4, ncol=4)
delta_t1 <-18
sampleSize1 <- 1378
moderator1 <- c(1, 1, 2)

################################ Primary Study 2 ###############################
delta_t2 <- c(12, 11)
sampleSize2 <- 668
empcov2 <- matrix(c(1.00, 0.40, 0.56,   NA, 0.53, 0.30,
                    0.40, 1.00, 0.31,   NA, 0.46, 0.53,
                    0.56, 0.31, 1.00,   NA, 0.62, 0.36,
                    NA,   NA,   NA,   NA,   NA,   NA,
                    0.53, 0.46, 0.62,   NA, 1.00, 0.46,
                    0.30, 0.53, 0.36,   NA, 0.46, 1.00), nrow=6, ncol=6)
empcov2 <- matrix(c(1.00, 0.41, 0.62,   NA, 0.58, 0.36,
                    0.41, 1.00, 0.33,   NA, 0.32, 0.58,
                    0.62, 0.33, 1.00,   NA, 0.58, 0.37,
                      NA,   NA,   NA,   NA,   NA,   NA,
                    0.58, 0.32, 0.58,   NA, 1.00, 0.46,
                    0.36, 0.58, 0.37,   NA, 0.46, 1.00), nrow=6, ncol=6)
moderator2 <- c(1, 2, 2)

################################ Primary Study 3 ###############################
empcov3 <- matrix(c(1.00, 0.44, 0.74, 0.36, 0.71, 0.32,
                    0.44, 1.00, 0.35, 0.66, 0.38, 0.65,
                    0.74, 0.35, 1.00, 0.43, 0.83, 0.35,
                    0.36, 0.66, 0.43, 1.00, 0.41, 0.71,
                    0.71, 0.38, 0.83, 0.41, 1.00, 0.44,
                    0.32, 0.65, 0.35, 0.71, 0.44, 1.00), nrow=6, ncol=6)
delta_t3 <- c(12,12)
moderator3 <- c(2, 2, 2)
pairwiseN3 <- matrix(c(658, 650,   0,  631,  599, 579,
                       650, 658,   0,  630,  588, 578,
                       0,   0,     0,    0,   0,    0,
                       631, 630,   0,  638, 623,  600,
                       599, 588,   0,  623, 630,  560,
                       579, 578,   0,  600, 560,  610), nrow=6, ncol=6)
isSymmetric(pairwiseN3)

################################ Primary Study 4 ###############################
delta_t4 <- 96
sampleSize4 <-556
empcov4 <- matrix(c(1.00, 0.39, 0.41, 0.19,
                    0.39, 1.00, 0.23, 0.44,
                    0.41, 0.23, 1.00, 0.26,
                    0.19, 0.44, 0.26, 1.00), nrow=4, ncol=4)
moderator4 <- c(2, 1, 1)
startValues4 <- c(-0.05, -0.02, 0.01,
                  -0.03, -1.19, 0.13, -1.48,
                  0.00, 0.45, -0.12)

############################### Primary Study 5 ################################
rawData5 <- list(fileName="data5.dat", missingValues=-99, standardize=TRUE,  header=TRUE, dec=".", sep=" ")
moderator5 <- c(2, 1, 1)

###################### COMPILE list() OF PRIMARY STUDIES #######################

firstStudyList <- compileListOfPrimaryStudies(selectedStudies=1:5)

################### Conduct CoTiMA - Study 2 will fit poorly ###################
CoTiMA(
  # Directory names and file names
  workingDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/",
  sourceDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/",
  resultsFileName="CoTiMA 1.txt",
  filePrefix="CoTiMA 1",
  # Primary Study Information
  primaryStudies=firstStudyList,
  nlatents=2,
  # Figure Parameters
  timeUnit="Months"
  )

####### COMPILE list() of studies including only primary study 2 and ... #######
###### ... do a CoTiMA again with smaller number of refits (and new name). #####
listOfPrimaryStudies <- compileListOfPrimaryStudies(selectedStudies=2)
CoTiMA(
  # Directory names and file names
  workingDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/",
  sourceDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/",
  resultsFileName="tmp.txt",
  filePrefix="tmp",
  # Primary Study Information
  primaryStudies=listOfPrimaryStudies,
  nlatents=2,
  # Figure Parameters
  timeUnit="Months",
  # Fitting Parameters
  refits=200,
  saveSingleStudyModelFit=c("tmp", 2)
)

#### Use importStartValues Function to assign them to the object startValues2 ###
startValues2 <- importStartValues(filePrefix="tmp", modelName="studyFit", modelNumber="1",
                         directory="/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/tmp singleStudyFits")
round(startValues2, 4)

### COMPILE list() of primary studies again (now start values for primary ... ##
# ... study 2 will be included) and re-do first CoTiMA ( and save model fits). #
listOfPrimaryStudies <- compileListOfPrimaryStudies(selectedStudies=1:5)
CoTiMA2Fit <- CoTiMA(
  # Directory names and file names
  workingDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/",
  sourceDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/",
  resultsFileName="CoTiMA 2.txt",
  filePrefix="CoTiMA 2",
  # Primary Study information
  primaryStudies=listOfPrimaryStudies,
  nlatents=2,
  # Figure Parameters
  timeUnit="Months",                    # timelag unit to lable x-axis of plots
  # Fitting Parameters
  refits=1,
  confidenceIntervals=TRUE,
  saveSingleStudyModelFit=c("CoTiMA 2", 1:5),
  saveDRIFTAllModelFit=c("CoTiMA 2")
)

### Moderator Model ###
CoTiMA3Fit <- CoTiMA(
  # Directory names and file names
  workingDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/",
  sourceDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/",
  resultsFileName="CoTiMA 3.txt",
  filePrefix="CoTiMA 3",
  # Primary Study information
  primaryStudies=listOfPrimaryStudies,
  nlatents=2,
  # Workflow
  checkSingleStudyResults=FALSE,
  # Figure Parameters
  timeUnit="Months",
  # Specific Model Tests
  testModeratorModel=TRUE,
  moderatorNumber = 1,
  # Fitting Parameters
  refits=1,
  confidenceIntervals=TRUE,
  loadSingleStudyModelFit=c("CoTiMA 2", 1:5),
  loadDRIFTAllModelFit=c("CoTiMA 2"),
  saveModeratorModelFit=c("CoTiMA 3")
)
CoTiMA3Fit

## FIT USER SPECIFIED MODEL
moderatorValues1 <- c(moderator1[1], moderator2[1], moderator3[1], moderator4[1], moderator5[1]); moderatorValues1

fixedUserModel1T0VAR <- matrix(c("groupfixed", "0",
                                 "groupfixed", "groupfixed"), 2, 2, byrow=TRUE)
fixedUserModel1DRIFT <- matrix(c("groupfixed", .070,
                                 "V1toV2", "V2toV2"), 2, 2, byrow=TRUE)
moderatedUserModel1DRIFT <- matrix(c("V1toV1", "moderated", # see text
                                     "moderated", "V2toV2"), 2, 2, byrow=TRUE)

listOfFixedModelMatrices1 <- list(DRIFT=fixedUserModel1DRIFT, T0VAR=fixedUserModel1T0VAR)
listOfModeratorModelMatrices1 <- list(DRIFT=moderatedUserModel1DRIFT)

CoTiMA4Fit <- CoTiMA(
  # Directory names and file names
  testUserSpecifiedModel=list(singleStudyModelFits = 1:5,
                          moderatorValues = moderatorValues1,
                          listOfFixedModelMatrices = listOfFixedModelMatrices1,
                          listOfModeratorModelMatrices = listOfModeratorModelMatrices1),
  workingDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/DOCUMENTATION/",
  sourceDirectory= "/Users/cdormann/SynologyDrive/Drive/CHRISTIAN/TEXTE/METHODEN/R/SYNTAXSAMMLUNG/CoTiMA/CURRENT VERSION/",
  resultsFileName="CoTiMA 4.txt",
  filePrefix="CoTiMA 4",
  # Primary Study information
  primaryStudies=listOfPrimaryStudies,
  nlatents=2,
  # Workflow
  checkSingleStudyResults=FALSE,
  # Figure Parameters
  timeUnit="Months",
  # Specific Model Tests
  #testModeratorModel=FALSE,
  #moderatorNumber = 1,
  # Fitting Parameters
  refits=1,
  confidenceIntervals=TRUE,
  testDRIFTallModel=FALSE,
  saveUserSpecifiedModel=c("CoTiMA 4"),
  loadSingleStudyModelFit=c("CoTiMA 2", 1:5)
)

1-pchisq(36528.48-35773.7, 50-26)


