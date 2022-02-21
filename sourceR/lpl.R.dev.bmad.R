##
## lpl.R.dev.bmad.R, creator S. Rauzy, LPL, 18/02/2022
##

##
## software : The software used to created the output
## projectname : The name of the project
## df : The output data frame created by HMAD from the parameters along time
## audf : The output data frame created by HMAD from the OpenFace AU output
##
lpl.R.dev.bmad.computeBMAD <- function(software, projectname, df, audf) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdir <- paste("projects/", projectname, sep="");
	projectdir <- paste(projectdir, "/", software, sep="");

	## Define the folders	
	FOLDER_MODEL <- paste(projectdir, "models", sep="/");

	blinkdir <- paste(projectdir, "/", "blink", sep="");
	if (!file.exists(blinkdir)) {
		dir.create(blinkdir);
	}
	FOLDER_TABLES <- paste(blinkdir, "tables", sep="/");
	if (!file.exists(FOLDER_TABLES)) {
		dir.create(FOLDER_TABLES);
	}

	filename <-  "bmad.txt";
	au_45_threshold <- 1;
	pitch_threshold <- -25;

	## If the result file for automatic annotation of head nods has already been created, just load it
	if (fileExists(FOLDER_TABLES, filename)) {
		cat("Automatic annotations of blinks previously done, loading the file", filename, "...\n");
		drf <- loadInternalDataFrame(FOLDER_TABLES, filename);
	} else {
		cat("Blink threshold =", au_45_threshold, "\n");
		cat("Compute the blink intervals annotation from AU45 values...\n");
		## Create the automatic annotations by using the AU45 variable 
		drf <- lpl.R.dev.bmad.createBlinkAnnotations(df, audf, au_45_threshold, pitch_threshold);
		
		## Save the result in a file
		saveInternalDataFrame(drf, FOLDER_TABLES, filename);
	}

	return (drf);
}

##
## Create the automatic blinks annotations predicted from the audf output data frame created by HMAD
## df : The output data frame created by HMAD from the parameters along time
## audf : The output data frame created by HMAD from the OpenFace AU output
## au_45_threshold : The threshold on the AU45_r value (Blink Action Unit of OpenFace, from 0 to 5).
## Return a data frame including tmin tmax time intervals with blink annotations (column blink)
##
lpl.R.dev.bmad.createBlinkAnnotations <- function(df, audf, au_45_threshold, pitch_threshold) {

	time <- audf$time;
	binary_values <- audf$AU45_c;
	values <- audf$AU45_r;
	pitch_values <- df$pitch;

	tmin <- 0;
	tmax <- 0;

	result <- NULL

	blink_area <- FALSE;
	mean_pitch <- 0;
	n_frames <- 0;

	for (i in c(1:nrow(audf))) {

		##if (!blink_area & binary_values[i] == 1) {
		if (!blink_area & values[i] > au_45_threshold) {
			tmin <- time[i];
			blink_area <- TRUE;
			mean_pitch <- 0;
			n_frames <- 0;
			
		}

		mean_pitch <- mean_pitch + pitch_values[i];
		n_frames <- n_frames + 1;

		##if (blink_area & binary_values[i] == 0) {
		if (blink_area  & values[i] < au_45_threshold) {
			tmax <- time[i];
			blink_area <- FALSE;
			mean_pitch <- mean_pitch/n_frames;
			if (mean_pitch < pitch_threshold) {
				result <- rbind(result, c(tmin, tmax, "P"));
			} else {
				result <- rbind(result, c(tmin, tmax, "B"));
			}
		}	
	}

	dfr <- data.frame(result);
	names(dfr) <- c("tmin", "tmax", "blink");

	return (dfr);
}

##
## Return the name of the eaf output file name function of the type
##
lpl.R.dev.bmad.getDefaultAnnotationFileName<- function() {

	return ("bmad.eaf");
}

##
## Create the ELAN eaf output file format for editing the automatic annotation of theblinks (in directory
## in projects/projectname/software/blink/result)
##
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.bmad.createElanEAFFile<- function(software, projectname, mediafilename, videomimetype) {

	bdatafile <- "bmad.txt";
	eaffilename <- "bmad.eaf";

	return (lpl.R.dev.bmad.createElanEAFFileWithName(software, projectname, mediafilename, eaffilename, videomimetype, bdatafile));
}

##
## Create the ELAN eaf output file format for editing the automatic annotation of theblinks (in directory
## in projects/projectname/software/blink/result)
##
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/blink/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
## bdatafile : The name of the  file name in projects/projectname/software/blink/tables containing the table of blink annotations
##
lpl.R.dev.bmad.createElanEAFFileWithName <- function(software, projectname, mediafilename, eaffilename, videomimetype, bdatafile) {

	## Define the folders
	projectdir <- paste("projects/", projectname, sep="");

	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	hnmaddir <- paste(projectdir, "/", "blink", sep="");
	
	## The folder for the tables
	FOLDER_TABLES <- paste(hnmaddir, "tables", sep="/");

	if (!fileExists(FOLDER_TABLES, bdatafile)) {
		cat("Please first run the lpl.R.dev.bmad.computeBMAD routine...");
	} else {
		bmad <- loadInternalDataFrame(FOLDER_TABLES, bdatafile);
		lpl.R.dev.bmad.cawriteElanEAFFile(bmad, software, projectname, mediafilename, eaffilename, videomimetype); 
	}
}

##
## Create and write the ELAN eaf output file format for editing the automatic annotation of the blinks (in directory
## in projects/projectname/software/blink/result)
##
## bmad : The data frame containing the data
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.bmad.cawriteElanEAFFile <- function(bmad, software, projectname, mediafilename, eaffilename, videomimetype) {

	projectdir <- paste("projects/", projectname, sep="");
	videodir <- projectdir;

	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	hnmaddir <- paste(projectdir, "/", "blink", sep="");

	FOLDER_RESULT <- paste(hnmaddir, "result", sep="/"); 
	if (!file.exists(FOLDER_RESULT)) {
		dir.create(FOLDER_RESULT);
	}
	WORKING_DIRECTORY <- getwd(); 

	MEDIA_URL <- "file://"; 
	MEDIA_URL <- paste(MEDIA_URL, WORKING_DIRECTORY, sep="/");
	MEDIA_URL <- paste(MEDIA_URL, "projects", sep="/");
	MEDIA_URL <- paste(MEDIA_URL, projectname, sep="/");
	MEDIA_URL <- paste(MEDIA_URL, mediafilename, sep="/");

	RELATIVE_MEDIA_URL <- "..";
	RELATIVE_MEDIA_URL <- paste(RELATIVE_MEDIA_URL, mediafilename, sep="/");

	cat(paste("Creating the Elan eaf annotation file in", FOLDER_RESULT, "...\n"));
	df <- lpl.R.dev.bmad.createElanEAFTable(bmad, MEDIA_URL, RELATIVE_MEDIA_URL, videomimetype);

	eaffile <- paste(FOLDER_RESULT, eaffilename, sep="/");

	write.table(df, eaffile, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");

	cat(paste("You can now open the annotation file", paste(FOLDER_RESULT, eaffilename, sep="/"), "with the ELAN Software.\n"));
}

##
## Create a one column data frame mimicking the xml ouput of the ELAN eaf format 
##
## da : The data frame containing the automatic annotation head nods automatic annotations to transform 
## mediaurl : the directory of the media file name
## relativemediaurl : the directory of the relative media file name
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
## return the xml document as a data frame
##
lpl.R.dev.bmad.createElanEAFTable <- function(da, mediaurl, relativemediaurl, videomimetype) {

	df <- NULL;

	## Create the xml header
	df <- addXmlHeader(df);

	df <- openAnnotationDocument(df, "BMAD", Sys.time());

	cat("Add head nods automatic annotations...\n");
	current_annotation_index = 0;
	cat("Add tier blink...\n");

       	tiersblink <- createAnnotationTier(current_annotation_index, da, "blink", "blink");

	current_annotation_index <- tiersblink[[2]];
	lastUsedAnnotationId <- current_annotation_index;

	df <- addHeaderHeader(df, mediaurl, relativemediaurl,  lastUsedAnnotationId, videomimetype); 
	cat("Add time order...\n");
	df <- addTimeOrder(df, da);

	cat("Add tiers...\n");
	df <- rbind(df, tiersblink[[1]]);

	df <- addConstraints(df);

	## Close the xml header
	df <- closeAnnotationDocument(df);
	colnames(df) <- NULL;
	
       return (df);	
}
