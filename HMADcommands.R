##
## R code for the HMAD project (Head movements Automatic Detection)
## author : S. Rauzy, LPL
## date : 18/02/2022

##
## Define the working directory where you have installed the HMAD-master directory
##
## WARNING: You have to change the next line in order to match your own instal path
HMAD_DIRECTORY <- "C:/Users/rauzy/Desktop/HMAD-master"; 
setwd(HMAD_DIRECTORY);

##
## Load the sources and packages
##
source("sourceR/loadSourcesAndPackages.R");

##
## Define the path of directory where OpenFace executable files are installed
## WARNING: You have to change the next line in order to match your own 
## instal path, in particular check the version number of OpenFace which 
## is updated regularly 
##
OPENFACE_DIRECTORY <- "C:/Users/rauzy/Desktop/OpenFace_2.0.5_win_x86";
##
## Define the software used to track the face
##
TRACKING_SOFTWARE <- "OPENFACE";

##
## List the HMAD projects
##
lpn <- printProjectNames();
##
## Define the project and create it
##
## The project name, please use CAPITAL letters
##
PROJECT_NAME <- "MY_NEW_PROJECT";
##
## Create the project
##
lpl.R.dev.hmad.createNewProject(TRACKING_SOFTWARE, PROJECT_NAME);

##
## Create the OPENFACE output
##
lpl.R.dev.hmad.createOpenFaceOutput(PROJECT_NAME, OPENFACE_DIRECTORY);

##
## Create Head Model and the internal facial landmarks residuals
## 
df <- lpl.R.dev.hmad.createHeadModelAndResiduals(TRACKING_SOFTWARE, PROJECT_NAME);
##
## Create the Action Units table
## 
audf <- lpl.R.dev.hmad.createActionUnitTable(TRACKING_SOFTWARE, PROJECT_NAME);
##
## Command to load the OpenFace csv file (if required) 
##
ofcsvdf <- lpl.R.dev.hmad.loadCsvData(TRACKING_SOFTWARE, PROJECT_NAME);

##
## EBMAD 
##
## EyeBrows Movements Automatic Detection
##
## Create the files for EBMAD
##
lpl.R.dev.ebmad.computeEBMAD(TRACKING_SOFTWARE, PROJECT_NAME, df);
##
## Create the Elan annotation file
##
## Define the parameters and file names
##
MEDIA_FILE_NAME <- lpl.R.dev.ebmad.getVideoFileName(PROJECT_NAME);
VIDEO_MIME_TYPE <- lpl.R.dev.ebmad.getElanMimeType(PROJECT_NAME);
type <- "RAISE_AND_FROWN";
EAF_ANNOTATION_FILE_NAME <- lpl.R.dev.ebmad.getDefaultAnnotationFileName(type);
##
## Create the Elan eaf annotation file name
##
lpl.R.dev.ebmad.createElanEAFFile(type, TRACKING_SOFTWARE, PROJECT_NAME, MEDIA_FILE_NAME, EAF_ANNOTATION_FILE_NAME, VIDEO_MIME_TYPE);

##
## SMAD
##
## Smile Movements Automatic Detection
##
## Create the files for SMAD
##
smad <- lpl.R.dev.smad.computeSMAD(TRACKING_SOFTWARE, PROJECT_NAME, df, audf);
##
## Create the Elan annotation file
##
## Define the parameters and file names
##
MEDIA_FILE_NAME <- lpl.R.dev.ebmad.getVideoFileNameWithRegex(PROJECT_NAME, "mp4");
VIDEO_MIME_TYPE <- lpl.R.dev.ebmad.getElanMimeTypeWithRegex(PROJECT_NAME, "mp4");
EAF_ANNOTATION_FILE_NAME <- lpl.R.dev.smad.getDefaultAnnotationFileName();
##
## Create the Elan eaf annotation file name
##
lpl.R.dev.smad.createElanEAFFile(TRACKING_SOFTWARE, PROJECT_NAME, MEDIA_FILE_NAME, EAF_ANNOTATION_FILE_NAME, VIDEO_MIME_TYPE);

##
## BMAD
##
## Blink Movements Automatic Detection
##
## Create the files for BMAD
##
bmad <- lpl.R.dev.bmad.computeBMAD(TRACKING_SOFTWARE, PROJECT_NAME, audf, 1);
##
## Create the Elan annotation file
##
## Define the parameters and file names
##
MEDIA_FILE_NAME <- lpl.R.dev.ebmad.getVideoFileNameWithRegex(PROJECT_NAME, "mp4");
VIDEO_MIME_TYPE <- lpl.R.dev.ebmad.getElanMimeTypeWithRegex(PROJECT_NAME, "mp4");
##
## Create the Elan eaf annotation file name
##
EAF_ANNOTATION_FILE_NAME <- "bmad.eaf";
lpl.R.dev.bmad.createElanEAFFile(TRACKING_SOFTWARE, PROJECT_NAME, MEDIA_FILE_NAME, EAF_ANNOTATION_FILE_NAME, VIDEO_MIME_TYPE, "bmad.txt");

##
## Plot the head model
##
library(ggplot2)
source("sourceR/lpl.R.dev.hmplot.R")
##
lpl.R.dev.hmplot.plotOpenFaceHeadModel(PROJECT_NAME);

