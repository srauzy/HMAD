##
## R code for the HMAD project (Head movements Automatic Detection)
## author : S. Rauzy, LPL
## date : 27/11/2017

##
## Define the working directory where you have installed the HMAD-master directory
##
HMAD_DIRECTORY <- "C:/Users/rauzy/Desktop/HMAD-master"; 
setwd(HMAD_DIRECTORY);

##
## Load the sources and packages
##
source("sourceR/loadSourcesAndPackages.R");

##
## Define the path of directory where OpenFace executable files are installed
##
OPENFACE_DIRECTORY <- "C:/Users/rauzy/Desktop/OpenFace_2.0.5_win_x86"
##
## Define the software used to track the face
##
TRACKING_SOFTWARE <- "OPENFACE"

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
df <- lpl.R.dev.hmad.createHeadModelAndResiduals(TRACKING_SOFTWARE, PROJECT_NAME)
##
## Create the Action Units table
## 
audf <- lpl.R.dev.hmad.createActionUnitTable(TRACKING_SOFTWARE, PROJECT_NAME)
