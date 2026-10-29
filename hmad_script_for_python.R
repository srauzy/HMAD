# -------------------------------------------------------------------------
# Copyright (C) 2018-2046 Stéphane Rauzy
# Laboratoire Parole et Langage, Aix-en-Provence, France
#
# Use of this software is governed by the GNU Public License, version 3.
#
# HMAD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HMAD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HMAD If not, see <http://www.gnu.org/licenses/>.
#
# This banner notice must not be removed.
# -------------------------------------------------------------------------

##
## R code allowing to call HMAD from python
## author : S. Rauzy, LPL
## date : 28/10/2026

## Retrieve the arguments of the python call 
args <- commandArgs(trailingOnly = TRUE)

## Define the working directory where the HMAD-master directory have been installed
HMAD_DIRECTORY <- args[1];
cat("HMAD working directory", HMAD_DIRECTORY, "\n")
flush.console() 

## Change the working directory of R for the hmad directory
setwd(HMAD_DIRECTORY);

## =========================================================================
## R VERSION CHECK & DEPENDENCY MANAGEMENT (renv)
## =========================================================================

## Define the minimum required R version for your code
MIN_R_VERSION <- "4.1.0"

if (getRversion() < MIN_R_VERSION) {
	stop(paste("FATAL ERROR: Your R version is", getRversion(), ". HMAD requires R version", MIN_R_VERSION, "or higher.", "Please update R from https://cran.r-project.org/"))
}

## Ensure 'renv' is available; if not, try to bootstrap it or install it locally
if (!requireNamespace("renv", quietly = TRUE)) {
	message("-> 'renv' package not found. Attempting to install it...")
	install.packages("renv", repos = "https://cloud.r-project.org/")
}

## Automatically restore the exact project packages using renv.lock without prompting the user
if (file.exists("renv.lock")) {
	message("-> Restoring R environment dependencies via renv...")
  	renv::restore(prompt = FALSE)
} else {
  	warning("-> Warning: 'renv.lock' not found. Packages will load from current library.")
}

## =========================================================================
## Load the sources and packages
## =========================================================================
source("sourceR/loadSourcesAndPackages.R");
source("sourceObjectR/loadObjectSourcesAndPackages.R");

## The tracking software
TRACKING_SOFTWARE <- "OPENFACE"

## The project name, please use CAPITAL letters
PROJECT_NAME <- args[2];
cat("Project name : ", PROJECT_NAME, "\n")
flush.console() 

## The location of the OpenFace directory if you want to launch OpenFace from R.
## Not recommanded, it works only on Windows system, prefer the python OpenFace integrated in Panamic distribution
OPENFACE_DIRECTORY <- args[3];
cat("OpenFace directory : ", OPENFACE_DIRECTORY, "\n")
flush.console() 

## Flag indicating if one needs to lauch OpenFace from R or if it has already been done previously,
## for example by a python call to OpenFace
FLAG_COMPUTE_OPENFACE_OUTPUT_WITH_PYTHON <- args[4];
cat("Flag compute OpenFace output with R :", FLAG_COMPUTE_OPENFACE_OUTPUT_WITH_PYTHON, "\n")
flush.console() 

## The video file path to treat
VIDEO_FILE_PATH <- args[5];
cat("Video file path : ", VIDEO_FILE_PATH, "\n")
flush.console() 

## Compute OpenFace output with the Windows executable version of OPenFace
if (FLAG_COMPUTE_OPENFACE_OUTPUT_WITH_PYTHON == "FALSE") {
	lpl.R.dev.hmad.createOpenFaceOutputFromVideoFilePath(PROJECT_NAME, VIDEO_FILE_PATH, OPENFACE_DIRECTORY)
}

## Create or load Head Model and the internal facial landmarks residuals
df <- lpl.R.dev.hmad.createHeadModelAndResiduals(TRACKING_SOFTWARE, PROJECT_NAME);

## Create or laod the Action Units table
audf <- lpl.R.dev.hmad.createActionUnitTable(TRACKING_SOFTWARE, PROJECT_NAME);

## Compute ebmad flag
COMPUTE_FLAG <- args[6];
if (COMPUTE_FLAG == "TRUE") {
	## EBMAD 
	## EyeBrows Movements Automatic Detection
	## Create the files for EBMAD
	lpl.R.dev.ebmad.computeEBMAD(TRACKING_SOFTWARE, PROJECT_NAME, df);
	## Create the Elan annotation file
	## Define the parameters and file names
	MEDIA_FILE_NAME <- lpl.R.dev.ebmad.getVideoFileName(PROJECT_NAME);
	VIDEO_MIME_TYPE <- lpl.R.dev.ebmad.getElanMimeType(PROJECT_NAME);
	type <- "RAISE_AND_FROWN";
	EAF_ANNOTATION_FILE_NAME <- lpl.R.dev.ebmad.getDefaultAnnotationFileName(type);
	## Create the Elan eaf annotation file name
	lpl.R.dev.ebmad.createElanEAFFile(type, TRACKING_SOFTWARE, PROJECT_NAME, MEDIA_FILE_NAME, EAF_ANNOTATION_FILE_NAME, VIDEO_MIME_TYPE);
}

## Compute smad flag
COMPUTE_FLAG <- args[7];
if (COMPUTE_FLAG == "TRUE") {
	## SMAD 
	## Smile Movements Automatic Detection
	## Create the files for SMAD
	smad <- lpl.R.dev.smad.computeSMAD(TRACKING_SOFTWARE, PROJECT_NAME, df, audf);
	## Create the Elan annotation file
	## Define the parameters and file names
	VIDEO_MIME_TYPE <- lpl.R.dev.smad.getElanMimeType(VIDEO_FILE_PATH);
	EAF_ANNOTATION_FILE_NAME <- lpl.R.dev.smad.getDefaultAnnotationFileName();
	filename <- basename(VIDEO_FILE_PATH)
	lpl.R.dev.smad.createElanEAFFileFromFilePath(TRACKING_SOFTWARE, PROJECT_NAME, VIDEO_FILE_PATH, filename, EAF_ANNOTATION_FILE_NAME, VIDEO_MIME_TYPE);
}

## Compute bagmad flag
COMPUTE_FLAG <- args[8];
if (COMPUTE_FLAG == "TRUE") {
	## BAGMAD
	## Blink And Gaze Movements Automatic Detection
	## Print comments in the console
	verbose <- TRUE;
	## Create the HMAD input from the OpenFaceOutput if it not already done
	hmad_object <- lpl.R.dev.hmad.hmadObject.createDefaultHmadObject(HMAD_DIRECTORY, PROJECT_NAME, verbose);
	## Define the project directory where the HMAD input files are stored
	PROJECT_DIRECTORY <- paste(HMAD_DIRECTORY, "projects", PROJECT_NAME, sep="/");
	## The video file and its extension for creating a link for Elan annotation file
	media_file <- VIDEO_FILE_PATH
	## Create the blink and gaze automatic annotation
	bagmad_object <- lpl.R.dev.hmad.bagmadObject.createDefaultBagmadObject(HMAD_DIRECTORY, PROJECT_NAME, hmad_object, media_file, verbose); 
}

##
## Plot the head model
#library(ggplot2)
#source("sourceR/lpl.R.dev.hmplot.R")
#lpl.R.dev.hmplot.plotOpenFaceHeadModel(PROJECT_NAME);
