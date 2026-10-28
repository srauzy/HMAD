##
## lpl.R.dev.hmad.hmadObject.R
## The functions defining an hmad object for head mouvements automatic detection
##

##
## Create the hmad object and return it
## return the hmad object
##
lpl.R.dev.hmad.hmadObject.newHmadObject <- function() {

	object <- vector(mode = "list", length = 9);

	object[[1]] <- "_HMAD_OBJECT";

	names(object) <- c("type", "openface_object", "landmark_file", "df_landmark", "au_file", "df_au", "gaze_file", "df_gaze", "FPS");  

	return (object);
}

##
## Set the OpenFace object
## object : The hmad object
## openface_object : The OpenFace object
## return the hmad object completed with the HMAD data frames
##
lpl.R.dev.hmad.hmadObject.setOpenfaceObject <- function(object, openface_object) {

	## Instanciate the file variable
	object[[2]] <- openface_object;

	return (object);
}

##
## Create the default hmad object and return it
## hmad_directory : The HMAD directory containing the "projects" directory 
## project_name : The project name
## verbose : put this flag to TRUE if you want message in the console
## return the hmad object completed with the 3 HMAD data frames and file names
##
lpl.R.dev.hmad.hmadObject.createDefaultHmadObject <- function(hmad_directory, project_name, verbose) {

	project_directory <- paste(hmad_directory, "projects", project_name, sep="/");
       	if (!file.exists(project_directory)) stop("Project does not exist, Please first create a project!\n");
	
	create_hmad_input <- lpl.R.dev.hmad.hmadObject.isHMADInputMissing(project_directory);

	TRACKING_SOFTWARE <- "OPENFACE";
	PROCESSED_OPENFACE_DIRECTORY <- paste(project_directory, TRACKING_SOFTWARE, "processed", sep="/");
	openface_csv_name <- paste(project_name, "ofo.csv", sep="_");
	openface_csv_file <-  paste(PROCESSED_OPENFACE_DIRECTORY, openface_csv_name, sep="/");

	openface_object <- lpl.R.dev.hmad.openfaceObject.newOpenfaceObject();
	openface_object <- lpl.R.dev.hmad.openfaceObject.setOpenFaceCsvFile(openface_object, openface_csv_file);
	
	if (create_hmad_input) {
		## Read the the OpenFace output and store it in an openface object
		openface_object <- lpl.R.dev.hmad.openfaceObject.readCsvFile(openface_object, verbose); 
	}

	hmad_object <- lpl.R.dev.hmad.hmadObject.newHmadObject();
	## Set the OpenFace output in the hmad object
	hmad_object <- lpl.R.dev.hmad.hmadObject.setOpenfaceObject(hmad_object, openface_object);
	## Release the openface object
	openface_object <- NULL;

	hmad_object <- lpl.R.dev.hmad.hmadObject.createHMADInput(hmad_object, project_directory, verbose);

	return (hmad_object); 
}

##
## Check if an HMAD input is missing and thus have to be created
## project_directory : The project directory 
## return TRUE if an HMAD input is missing
##
lpl.R.dev.hmad.hmadObject.isHMADInputMissing <- function(project_directory) {

	file <- paste(project_directory, "OPENFACE", "tables", "itdcomplete.txt", sep="/");
	if (!file.exists(file)) return (TRUE);

	file <- paste(project_directory, "OPENFACE", "tables", "audf.txt", sep="/");
	if (!file.exists(file)) return (TRUE);

	file <- paste(project_directory, "OPENFACE", "tables", "egdf.txt", sep="/");
	if (!file.exists(file)) return (TRUE);

	return (FALSE);
}

##
## Create the 3 HMAD data frames (landmark, au and gaze) from the csv OpenFace output data
## object : The hmad object
## project_directory : The directory where to write the 3 files (in project_directory/OPENFACE/tables directory) 
## verbose : The directory where to write the 3 files (in project_directory/OPENFACE/tables directory) 
## return the hmad object completed with the 3 HMAD data frames and file names
##
lpl.R.dev.hmad.hmadObject.createHMADInput <- function(object, project_directory, verbose) {

	openface_subdirectory <- paste(project_directory, "OPENFACE", sep="/");
	if (!file.exists(openface_subdirectory)) {
		dir.create(openface_subdirectory, recursive=TRUE);
	}	
	object[[3]] <- paste(project_directory, "OPENFACE", "tables", "itdcomplete.txt", sep="/");
	if (!file.exists(object[[3]])) {
		dir.create(paste(project_directory, "OPENFACE", "tables", sep="/"), recursive = TRUE);
		df_landmark <- lpl.R.dev.openFaceOutputAnalysis.createHeadModelAndResiduals(paste(project_directory, "OPENFACE", sep="/"), object$openface_object$df_openface);
	} else {
		if (verbose) cat("Read the HMAD landmarks data in", object[[3]], "...\n");
		df_landmark <- read.delim(object[[3]], fileEncoding="UTF-8");
	}
	object[[4]] <- df_landmark;
	## The function save the data frame in the file

	object[[5]] <- paste(project_directory, "OPENFACE", "tables", "audf.txt", sep="/");
	if (!file.exists(object[[5]])) {
		df_au <- lpl.R.dev.openFaceOutputAnalysis.createActionUnitTable(object$openface_object$df_openface);
		## Save the data frame in the file
		write.table(df_au, object$au_file, sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8");
		if (verbose) cat("Save AU in", object$au_file, "\n"); 
	} else {
		if (verbose) cat("Read the HMAD AU data in", object[[5]], "...\n");
		df_au <- read.delim(object[[5]], fileEncoding="UTF-8");
	}
	object[[6]] <- df_au;

	object[[7]] <- paste(project_directory, "OPENFACE", "tables", "egdf.txt", sep="/");
	if (!file.exists(object[[7]])) {
		df_gaze <- lpl.R.dev.openFaceOutputAnalysis.createEyesGazeTable(object$openface_object$df_openface);
		## Save the data frame in the file
		write.table(df_gaze, object$gaze_file, sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8");
		if (verbose) cat("Save gaze data in", object$gaze_file, "\n"); 
	} else {
		if (verbose) cat("Read the HMAD gaze data in", object[[7]], "...\n");
		df_gaze <- read.delim(object[[7]], fileEncoding="UTF-8");
	}
	object[[8]] <- df_gaze;
	
	## Instantiate the FPS variable too	
	object <- lpl.R.dev.hmad.hmadObject.computeFrameRatePerSecond(object);

	return (object);
}

##
## Load the default 3 HMAD data frames (landmark, au and gaze) and the 3 associated files
## object : The hmad object
## project_directory : The directory containing the 3 files (in project_directory/OPENFACE/tables directory)
## verbose : put this flag to TRUE if you want message in the console
## return the hmad object completed with the 3 HMAD data frames and file names
##
lpl.R.dev.hmad.hmadObject.loadDefaultHMADInput <- function(object, project_directory, verbose) {

	object[[3]] <- paste(project_directory, "OPENFACE", "tables", "itdcomplete.txt", sep="/");
	object[[5]] <- paste(project_directory, "OPENFACE", "tables", "audf.txt", sep="/");
	object[[7]] <- paste(project_directory, "OPENFACE", "tables", "egdf.txt", sep="/");
	
	object <- lpl.R.dev.hmad.hmadObject.readHMADInput(object, verbose); 

	## Instantiate the FPS variable too	
	object <- lpl.R.dev.hmad.hmadObject.computeFrameRatePerSecond(object);

	return (object);
}

##
## Set the 3 files to read the HMAD input in
## object : The hmad object
## landmark_file : The landmark file
## au_file : The AU file
## gaze_file : The gaze file
## return the hmad object completed with the 3 HMAD input files
##
lpl.R.dev.hmad.hmadObject.setHMADInputFiles <- function(object, landmark_file, au_file, gaze_file) {

	## Instanciate the file variable
	object[[3]] <- landmark_file;
	object[[5]] <- au_file;
	object[[7]] <- gaze_file;

	return (object);
}

##
## Read the 3 HMAD input file
## object : The hmad object
## verbose : put this flag to TRUE if you want message in the console
## return the hmad object completed with the 3 HMAD data frame input instantiated
##
lpl.R.dev.hmad.hmadObject.readHMADInput <- function(object, verbose) {

	## read the files
	if (verbose) cat("Read the HMAD landmarks data in", object$landmark_file, "...\n");
	object[[4]] <- read.delim(object$landmark_file, fileEncoding="UTF-8");
	if (verbose) cat("Read the HMAD AU data in", object$au_file, "...\n");
	object[[6]] <- read.delim(object$au_file, fileEncoding="UTF-8");
	if (verbose) cat("Read the HMAD gaze data in", object$gaze_file, "...\n");
	object[[8]] <- read.delim(object$gaze_file, fileEncoding="UTF-8");

	## Instantiate the FPS variable too	
	object <- lpl.R.dev.hmad.hmadObject.computeFrameRatePerSecond(object);

	return (object);
}

##
## Compute the frame rate per second of the video (FPS,  i.e. number of frames per second for data frame, or 1/timestep) 
## object : The hmad object
## return the hmad object completed with the FPS instantiated
##
lpl.R.dev.hmad.hmadObject.computeFrameRatePerSecond <- function(object) {

	if (!is.null(object$df_landmark)) {
		object[[9]] <- round(1/(object$df_landmark$time[2] - object$df_landmark$time[1]), 0);
	}
	return (object);
}

##
## Get the name of the file (not directory) in the project directory (if all is right, it is the video file) and return it
## project_directory : The project directory
## extension : The extension of the video file
##
## return the video file name or null in case of problem
##
lpl.R.dev.hmad.hmadObject.getVideoFileName <- function(project_directory, extension) {

	regex <- paste(".*\\.", extension, sep="");
	listfiles <- list.files(project_directory, pattern=regex);
	if (length(listfiles) < 1) {
		return (cat(paste("No video file found in project directory", project_directory, "...")));
	} else if (length(listfiles) > 1) {
		return (cat(paste("More than one candidate for the video file in the project directory", project_directory, "...")));
	}
	return (listfiles[1]);
}

##
## Get the video file (not directory) in the project directory (if all is right, it is the video file) and return it
## project_directory : The project directory
## extension : The extension of the video file
##
## return the video file or null in case of problem
##
lpl.R.dev.hmad.hmadObject.getVideoFile <- function(project_directory, extension) {

	regex <- paste(".*\\.", extension, sep="");
	listfiles <- list.files(project_directory, pattern=regex);
	if (length(listfiles) < 1) {
		return (cat(paste("No video file found in project directory", project_directory, "...")));
	} else if (length(listfiles) > 1) {
		return (cat(paste("More than one candidate for the video file in the project directory", project_directory, "...")));
	}
	return (paste(project_directory, listfiles[1], sep="/"));
}

##
## Get the video file (not directory) in the project directory (if all is right, it is the video file) and return it
## project_directory : The project directory
## extension : The extension of the video file
##
## return the video file or null in case of problem
##
lpl.R.dev.hmad.hmadObject.getVideoFileOld <- function(project_directory) {

	df_video <- read.delim(paste(project_directory, "OPENFACE", "default_video_file.txt", sep="/"), fileEncoding="UTF-8", header=F);	
	return (df_video[1,1]);
}
