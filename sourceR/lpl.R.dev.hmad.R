## R code for computing head model and residuals motions, part of the HMAD project
## author : S. Rauzy, LPL
## date : 27/02/2017

##
## Compute the internal facial movements from the csv output file
##
## software : The software used for extracting the video output 
## projectname : The project name
##
## return the data frame containing the residuals
##
lpl.R.dev.hmad.createHeadModelAndResiduals <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	## OpenFace output
	if (software == "OPENFACE") {
	
		FOLDER_TABLES <- paste(projectdir, "tables", sep="/");

		if (fileExists(FOLDER_TABLES, "itdcomplete.txt")) {
			cat("STEP 1 previously done, load the data with residuals...\n");
			d <- loadInternalDataFrame(FOLDER_TABLES, "itdcomplete.txt");	
		} else {
			filename <- paste(projectname, "_ofo.csv", sep="");
			csvprojectdir <- paste(projectdir, "/processed/", sep="");
			cat(paste("LOAD CSV FILE", filename, "in directory", csvprojectdir, "...\n"));
			csv <- lpl.R.dev.openFaceOutputAnalysis.loadOpenFaceSingleFaceOutput(csvprojectdir, filename);

			d <- lpl.R.dev.openFaceOutputAnalysis.createHeadModelAndResiduals(projectdir, csv);
		}
	## IntraFace output	
	} else if (software == "INTRAFACE") {
		d <- lpl.R.dev.hmad.createHeadModelSFactorAndResiduals(software, projectname); 
	}

	FOLDER_RESULT <- paste(projectdir, "result", sep="/"); 
	if (!file.exists(FOLDER_RESULT)) {
		dir.create(FOLDER_RESULT);
	}
	##if (!fileExists(FOLDER_RESULT, "description.html")) {
		lpl.R.dev.hmad.createHtmlDescription(d, software, projectname);
	##}

	return (d);
}

##
## Create the Action Units table from the csv output file
##
## software : The software used for extracting the video output 
## projectname : The project name
##
## return the data frame containing the AU measurements
##
lpl.R.dev.hmad.createActionUnitTable <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	## OpenFace output
	if (software == "OPENFACE") {
	
		FOLDER_TABLES <- paste(projectdir, "tables", sep="/");

		if (fileExists(FOLDER_TABLES, "audf.txt")) {
			cat("STEP 1 previously done, load the data with action units...\n");
			d <- loadInternalDataFrame(FOLDER_TABLES, "audf.txt");	
		} else {
			filename <- paste(projectname, "_ofo.csv", sep="");
			csvprojectdir <- paste(projectdir, "/processed/", sep="");
			cat(paste("LOAD CSV FILE", filename, "in directory", csvprojectdir, "...\n"));
			csv <- lpl.R.dev.openFaceOutputAnalysis.loadOpenFaceSingleFaceOutput(csvprojectdir, filename);

			d <- lpl.R.dev.openFaceOutputAnalysis.createActionUnitTable(csv);

			if (!file.exists(FOLDER_TABLES)) {
				dir.create(FOLDER_TABLES);
			}
			cat("Save OpenFace action units file in tables/audf.txt ...\n");
			saveInternalDataFrame(d, FOLDER_TABLES, "audf.txt");
		}
	}

	return (d);
}

##
## Load the csv output file of the tracking software
##
## software : The software used for extracting the video output 
## projectname : The project name
##
## return the data frame containing the csv
##
lpl.R.dev.hmad.loadCsvData <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	## OpenFace output
	if (software == "OPENFACE") {
		filename <- paste(projectname, "_ofo.csv", sep="");
		csvprojectdir <- paste(projectdir, "/processed/", sep="");
		cat(paste("LOAD CSV FILE", filename, "in directory", csvprojectdir, "...\n"));
		csv <- lpl.R.dev.openFaceOutputAnalysis.loadOpenFaceSingleFaceOutput(csvprojectdir, filename);
	}

	return (csv);
}

##
## cvs output file
##
## software : The software used for extracting the video output 
## projectname : The project name
##
## return the data frame containing the AU measurements
##
lpl.R.dev.hmad.createEyebrowActionUnitTable <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	## OpenFace output
	if (software == "OPENFACE") {
	
		FOLDER_TABLES <- paste(projectdir, "tables", sep="/");

		filename <- paste(projectname, "_ofo.csv", sep="");
		csvprojectdir <- paste(projectdir, "/processed/", sep="");
		cat(paste("LOAD CSV FILE", filename, "in directory", csvprojectdir, "...\n"));
		csv <- lpl.R.dev.openFaceOutputAnalysis.loadOpenFaceSingleFaceOutput(csvprojectdir, filename);

		d <- lpl.R.dev.openFaceOutputAnalysis.createEyebrowActionUnitTable(csv);
	}

	return (d);
}

##
## Create an html document containing the summary of the analysis
##
## data : The data table as a data frame
## software : The software used for extracting the video output 
## projectname : The project name
##
lpl.R.dev.hmad.createHtmlDescription <- function(data, software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	topprojectdir <- projectdir;
	projectdir <- paste(projectdir, "/", software, sep="");

	FOLDER_RESULT <- paste(projectdir, "result", sep="/"); 

	WORKING_DIRECTORY <- getwd();

	df <- NULL;
	df <- lpl.R.dev.htmldf.openHtmlTag(df);
	df <- lpl.R.dev.htmldf.openHtmlHeader(df, paste(projectname, "project"));
	df <- lpl.R.dev.htmldf.closeHtmlHeader(df);

	df <- lpl.R.dev.htmldf.openHtmlBody(df);

	line <- paste("Project name:", projectname, sep=" ");
	
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("R working directory:", WORKING_DIRECTORY, sep=" ");
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);

	line <- paste("Project directory:", paste(WORKING_DIRECTORY, topprojectdir, sep="/"));
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	lof <- getVideoFile(projectname);
	line <- paste("In the project directory, media file:",  lof[1], sep=" ");
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);

	line <- paste("Total number of frames:",  (nrow(data)), sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	line <- paste("Total duration:",  (nrow(data)/25), "seconds or", (lpl.R.dev.faceOutputAnalysis.stomnsString(nrow(data)/25)), sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	## OpenFace output
	if (software == "OPENFACE") {
		csvprojectdir <- paste(projectdir, "/processed/", sep="");
		line <- paste("OpenFace software ouput in directory", csvprojectdir, sep=" ");
		line <- lpl.R.dev.htmldf.boldHtml(line);
		df <- lpl.R.dev.htmldf.addLine(df, line);
		lof <- lpl.R.dev.hmad.getCsvFiles(csvprojectdir);
		line <- paste("Csv output:",  lof[1], sep=" ");
		df <- lpl.R.dev.htmldf.addLine(df, line);

		nlc <- nrow(subset(data, data$label == "N"));
		nt <- nrow(data);
		line <- paste("Number of frames with low confidence in head tracking (c < 0.8):",  nlc, "on", nt, sep=" ");
		df <- lpl.R.dev.htmldf.addLine(df, line);	
	}

	df <- lpl.R.dev.htmldf.addLine(df, "");
	line <- lpl.R.dev.htmldf.boldHtml("Characteristics of the video, head model and residual landmark movements:");
	df <- lpl.R.dev.htmldf.addLine(df, line);

	globaldf <- paste(projectdir, "tables", "itdcomplete.txt", sep="/");
	line <- paste("Global data file with residuals:", globaldf, sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.openHtmlTable(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Total number of frames");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(data));
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of frames far from the frontal head position");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(subset(data, data$label == "L")));
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of frames far from the central XY position");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(subset(data, data$label == "P")));
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of frames far from the mean depth position (Z-axis)");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(subset(data, data$label == "S")));
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);
	
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of frames with low confidence in head tracking");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(subset(data, data$label == "N")));
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of frames remaining for head model computation");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(subset(data, data$label == "Y")));
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.closeHtmlTable(df);

	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("Head model file:", paste(projectdir, "model", "dhmt.txt", sep="/"), sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	line <- paste("Measurement of the dimension of the head in pixels:", round(lpl.R.dev.openFaceOutputAnalysis.meanHeadDimensionInPixels(software, projectname), 2) , sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	
 
	line <- paste("Dispersion of residual landmark movements file:", paste(projectdir, "model", "rrsdt.txt", sep="/"), sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	line <- paste("Mean dispersion of residuals in pixels:", round(lpl.R.dev.faceOutputAnalysis.meanDispersionInPixels(software, projectname), 2) , sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	cat(line);

	df <- lpl.R.dev.htmldf.closeHtmlBody(df);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("Eyebrows raising and frowning actions:", "<A HREF=\"ebmaddescription.html\">Result of the automatic annotation</A>", sep=" ");
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);


	df <- lpl.R.dev.htmldf.closeHtmlTag(df);

	descriptionfile <- paste(FOLDER_RESULT, "description.html", sep="/");
	cat("\n");
	cat(paste("Description of the video can be found in the file description.html (folder ", FOLDER_RESULT, ")", sep=""));
	cat("\n");

	write.table(df, descriptionfile, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");
}

##
## Compute the internal facial movements from the cvs output file
##
## software : The software used for extracting the video output 
## projectname : The project name 
##
lpl.R.dev.hmad.createHeadModelSFactorAndResiduals <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	## OpenFace output
	if (software == "OPENFACE") {

		FOLDER_TABLES <- paste(projectdir, "tables", sep="/");

		if (fileExists(FOLDER_TABLES, "itdcomplete.txt")) {
		
			##
			## Transform the OpenFace csv file in a R data frame and save it as an internal file
			##
			filename <- paste(projectname, "_ofo.csv", sep="");
		
			if (fileExists(FOLDER_TABLES, "iodf.txt")) {
				cat("STEP 1 previously done, load the R data frame of the OpenFace ouput...\n");
				iodf <- loadInternalDataFrame(FOLDER_TABLES, "iodf.txt");
			} else {
				cat("---------- STEP 1 ------------------\n");
				cat("Transform the OpenFace csv output file in an R data frame...\n");
				iodf  <- lpl.R.dev.openFaceOutputAnalysis.loadOpenFaceSingleFaceOutput(projectdir, filename);
				if (!file.exists(FOLDER_TABLES)) {
					dir.create(FOLDER_TABLES);
				}
				## Save the table in an internal file
				saveInternalDataFrame(iodf, FOLDER_TABLES, "iodf.txt");
			}

			###
			## Second step : Transform the output with R data frame format,
			## and prepare them (data cleaning) in order to create head model.
			## 

			## Transform the OpenFace csv output file in an R data frame
			if (fileExists(FOLDER_TABLES, "iodfcleaned.txt")) {
				cat("STEP 2 previously done, load the cleaned data...\n");
				iodfc <- loadInternalDataFrame(FOLDER_TABLES, "iodfcleaned.txt");
			} else {
				cat("---------- STEP 2 ------------------\n");
				cat("Prepare the data (data cleaning) in order to create head model...\n");
				iodfc <- lpl.R.dev.openFaceOutputAnalysis.prepareDataFrameOpenFaceOutput(iodf);
				## Save the table in an internal file
				saveInternalDataFrame(iodfc, FOLDER_TABLES, "iodfcleaned.txt");
			}
		}

	} else if (software == "INTRAFACE") {
		
		## Transform the csv file in a R data frame an save it as an internal file.
		FOLDER_TABLES <- paste(projectdir, "tables", sep="/");

		## Transform the Intraface csv output file in an R data frame
		if (fileExists(FOLDER_TABLES, "iodf.txt")) {
			cat("STEP 1 previously done, load the R data frame of the Intraface ouput...\n");
			iodf <- loadInternalDataFrame(FOLDER_TABLES, "iodf.txt");
		} else {
			## Get the Intraface csv output file
			csvfiles <- lpl.R.dev.hmad.getCsvFiles(projectdir);

			## If you have more than csv file in your folder, please select
			## the index of the csv file you want to work on, else let it to 1
			CSV_FILE_INDEX = 1
			csvfile <- csvfiles[CSV_FILE_INDEX];
			
			cat("---------- STEP 1 ------------------\n");
			cat("Transform the Intraface csv output file in an R data frame...\n");
			iodf  <- lpl.R.dev.intraFaceOutputAnalysis.loadIntrafaceOutput12Points(projectdir, csvfile);
			if (!file.exists(FOLDER_TABLES)) {
				dir.create(FOLDER_TABLES);
			}
			## Save the table in an internal file
			saveInternalDataFrame(iodf, FOLDER_TABLES, "iodf.txt");
		}
		
		###
		## Second step : Transform the output with R data frame format,
		## and prepare them (data cleaning) in order to create head model.
		## 

		FOLDER_TABLES <- paste(projectdir, "tables", sep="/");
		## Transform the Intraface csv output file in an R data frame
		if (fileExists(FOLDER_TABLES, "iodfcleaned.txt")) {
			cat("STEP 2 previously done, load the cleaned data...\n");
			iodfc <- loadInternalDataFrame(FOLDER_TABLES, "iodfcleaned.txt");
		} else {
			cat("---------- STEP 2 ------------------\n");
			cat("Prepare the data (data cleaning) in order to create head model...\n");
			iodfc <- lpl.R.dev.intraFaceOutputAnalysis.prepareDataFrameIntrafaceOutput(iodf);
			## Save the table in an internal file
			saveInternalDataFrame(iodfc, FOLDER_TABLES, "iodfcleaned.txt");
		}

	}

	FOLDER_MODEL <- paste(projectdir, "model", sep="/");
	## If the model has already been computed, do nothing
	if (fileExists(FOLDER_TABLES, "itdcomplete.txt")) {
		cat("STEP 3 previously done, load the data with residuals...\n");
		df <- loadInternalDataFrame(FOLDER_TABLES, "itdcomplete.txt");
	} else {
		cat("---------- STEP 3 ------------------\n");
		cat("Create the head model from the filtered data (data cleaning) and compute residuals...\n");
		d <- lpl.R.dev.faceOutputAnalysis.createHeadAndSFactorModel(iodfc);
		## Compute the head model with the las version of the head factor
		drc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(d, "direct", "XYZ")
		dhmt <- lpl.R.dev.faceOutputAnalysis.createDirectHeadModel(d, drc, FALSE);

		irc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(d, "inverse", "XYZ")
		r <- lpl.R.dev.faceOutputAnalysis.createRelativeProjectedResiduals(d, irc, dhmt);
		rrsdt <- lpl.R.dev.faceOutputAnalysis.computeProjectedResidualsStandardDispersion(r);

		if (!file.exists(FOLDER_MODEL)) {
			dir.create(FOLDER_MODEL);
		}
		## Save the head model data frame
		write.table(dhmt, paste(FOLDER_MODEL, "dhmt.txt", sep="/"), sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE);
		
		df <- lpl.R.dev.faceOutputAnalysis.addResidualsAndErrorsAmplitudes(iodfc, dhmt, rrsdt);
		rsdt <- lpl.R.dev.faceOutputAnalysis.computeProjectedResidualsStandardDispersion(df);
	       	## Save the head model data frame
		write.table(rsdt, paste(FOLDER_MODEL, "rrsdt.txt", sep="/"), sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE);	
		saveInternalDataFrame(df, FOLDER_TABLES, "itdcomplete.txt");
	}

	return (df);
}

##
## Return true if the file exists
## directory : The directory
## filename : The fine name
##
fileExists <- function(directory, filename) {

	return (file.exists(paste(directory, filename, sep="/")));
}

##
## Print project names on the console
##
printProjectNames <- function() {

	if (!file.exists("projects")) {
		cat("No project found!\n");
		cat("Please use the command 'lpl.R.dev.hmad.createNewProject(TRACKING_SOFTWARE, PROJECT_NAME);' to create a project...\n");
	} else {	
		listfiles <- list.files("projects");
	
		cat(paste("List of existing projects :\n")); 
		for (i in c(1:length(listfiles))) {
				cat(paste("Project", i, ":", listfiles[i], "\n")); 
		}
		return (listfiles);
	}
}

##
## Create a new project (a directory named projectname in the directory  "HMAD/projects")
##
## software : The software used for extracting the video output 
## projectname : The project name 
##
lpl.R.dev.hmad.createNewProject <- function(software, projectname) {

		## Create the directory "projects" at the root of the "HMAD" folder if it does not exist
	if (!file.exists("projects")) {
		dir.create("projects");	
	}

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	softprojectdir <- paste(projectdir, "/", software , sep="");

	if (!file.exists(projectdir) | !file.exists(softprojectdir)) {
		
		if (!file.exists(projectdir)) {
			dir.create(projectdir);
		}
		
		dir.create(softprojectdir);
		cat(paste("Project", projectname, "has been created.\n")); 
		if (software == "INTRAFACE") {
			cat(paste("Please drop or copy your IntraFace csv output file in the folder", softprojectdir, "\n")); 	
		} else if (software == "OPENFACE") {
			listfiles <- list.files(softprojectdir);
			if (length(listfiles) != 1) {
				cat(paste("Please drop or copy the video file you want to analyse in the folder", projectdir, "\n", "Afterwards, please run the command:\nlpl.R.dev.hmad.createOpenFaceOutput(PROJECT_NAME, OPENFACE_DIRECTORY)\n", "in order to create the OpenFace csv output file.", "\n")); 
			}	
		}
	} else {
		cat(paste("Project", projectname, "already exists !\n"));
		cat("------------------\n");
		printProjectNames();
		cat("------------------\n");
		warning(paste("Please choose another project name.", "\n", "Or if you wish to replace the existing project, please first delete manually the folder", projectdir)); 	
	}
}

##
## Get the list of csv files in the directory "projecdir"
##
## projectdir : The project directory to be scanned
##
lpl.R.dev.hmad.getCsvFiles <- function(projectdir) {

	
	listfiles <- list.files(projectdir, pattern=".*csv");
	cat(paste("In directory", projectdir, ":\n")); 
	for (i in c(1:length(listfiles))) {
		cat(paste("Csv file", i, ":", listfiles[i], "\n")); 
	}
	if (length(listfiles) < 1) {
		warning(paste("WARNING : Please copy the Intraface csv output file in your project folder", projectdir, "!\n")); 
	} else if ((length(listfiles) > 1)) {
	 	warning("WARNING : You have more than one csv file in your project folder! Please select the index of the csv file you want to work on.\n"); 
	}
	return (listfiles);
}

getVideoFile <- function(projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	listfiles <- list.files(projectdir, pattern=".*\\.mp.*|.*\\.avi|.*\\.MXF|.*\\.mxf|.*\\.wmv");
	cat(paste("In directory", projectdir, ":\n")); 
	for (i in c(1:length(listfiles))) {
		cat(paste("File", i, ":", listfiles[i], "\n")); 
	}
	if (length(listfiles) < 1) {
		warning(paste("WARNING : Please drop the video file in your project folder", projectdir, "!\n"));
	} else if ((length(listfiles) > 1)) {
	 	warning("WARNING : You have more than one file in your folder project ! Please select the index of the file you want to work on.\n");
	}
	return (listfiles);	
}

##
## Get and return the list of files (not directory) in the directory "projecdir"
## Print a message if there is not a single file 
##
## projectname : The project name to be scanned
##
lpl.R.dev.hmad.getListOfFiles <- function(projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	listfiles <- setdiff(list.files(projectdir), list.dirs(projectdir, recursive = FALSE, full.names = FALSE));
	
	cat(paste("In directory", projectdir, ":\n")); 
	for (i in c(1:length(listfiles))) {
		cat(paste("File", i, ":", listfiles[i], "\n")); 
	}
	if (length(listfiles) < 1) {
		warning(paste("WARNING : Please drop the file in your project folder", projectdir, "!\n"));
	} else if ((length(listfiles) > 1)) {
	 	warning("WARNING : You have more than one file in your project folder! Please select the index of the file you want to work on.\n");
	}
	return (listfiles);	
}

##
## Create the cvs output file of the OpenFace software
##
## projectname : The name of the project in the directory "projects"
## openfacedirectory : The absolute path where the OpenFace executable file "FeatureExtraction.exe" is installed
##
## see https://github.com/TadasBaltrusaitis/OpenFace/wiki/Command-line-arguments
## and https://github.com/TadasBaltrusaitis/OpenFace/wiki/Output-Format#featureextraction
##
lpl.R.dev.hmad.createOpenFaceOutput <- function(projectname, openfacedirectory) {

	projectdir <- paste("projects/", projectname, sep="");
	
	listfiles <- lpl.R.dev.hmad.getListOfFiles(projectname);

	## Store the HMAD working directory 
	hmaddirectory <- getwd();

	if (length(listfiles) == 1) {
		videofile <- paste(hmaddirectory, "/", projectdir, "/", listfiles[1], sep="");
		ofname =  paste(hmaddirectory, "/", projectdir, "/OPENFACE/processed/", projectname, "_ofo.csv", sep="");
		## Check whether the video has not been already processed (csv file already there) 
		if (file.exists(ofname)) {
			warning(paste("WARNING : It looks like the video file has been already processed!\nIf you wish to re-process the video, please delete manually the directory ", paste(hmaddirectory, "/", projectdir, "/OPENFACE/processed", sep=""), " and launch the command createOpenFaceOutput(PROJECT_NAME, OPENFACE_DIRECTORY) again!\n"));
		} else {
			## Change the working directory to OpenFace directory executable
			setwd(openfacedirectory);
			inputname = paste("FeatureExtraction.exe", "-f", videofile, "-of", ofname);
			inputname = paste(inputname, "-verbose");
			system(command = "cmd.exe", input = inputname);
		        ## Back to the HMAD working directory  	
			setwd(hmaddirectory);
		}
	}
}

##
## Create the cvs output file of the OpenFace software (old version)
##
## projectname : The name of the project in the directory "projects"
## openfacedirectory : The absolute path where the OpenFace executable file "FeatureExtraction.exe" is installed 
## saveTrackedVideo : Put this flag to TRUE if you want to save the decorated video output 
##
lpl.R.dev.hmad.createOpenFaceOutput1 <- function(projectname, openfacedirectory, saveTrackedVideo) {

	projectdir <- paste("projects/", projectname, sep="");
	
	listfiles <- getVideoFile(projectname);

	hmaddirectory <- getwd();
	setwd(openfacedirectory);

	if (length(listfiles) == 1) {
		videofile <- paste(hmaddirectory, "/", projectdir, "/", listfiles[1], sep="");
		ofname =  paste(hmaddirectory, "/", projectdir, "/OPENFACE/", projectname, "_ofo.csv", sep="");
		inputname = paste("FeatureExtraction.exe", "-f", videofile, "-of", ofname);
		inputname = paste(inputname, "-verbose");
		if (saveTrackedVideo) {
			ovname =  paste(hmaddirectory, "/", projectdir, "/OPENFACE/", projectname, "_ofo.avi", sep="");
			inputname = paste(inputname, "-ov", ovname); 
		}
		system(command = "cmd.exe", input = inputname); 
	}
	setwd(hmaddirectory);
}

##
## Get the list of csv files in the directory "projecdir"
##
## projectdir : The project directory to be scanned
##
getAviFiles <- function(projectdir) {

	listfiles <- list.files(projectdir, pattern=".*avi");
	cat(paste("In directory", projectdir, ":\n")); 
	for (i in c(1:length(listfiles))) {
		cat(paste("Avi file", i, ":", listfiles[i], "\n")); 
	}
	if (length(listfiles) < 1) {
		cat(paste("WARNING : Please drop the avi media file in your project folder", projectdir, "!\n")); 
	} else if ((length(listfiles) > 1)) {
	 	print("WARNING : You have more than one avi file in your project folder ! Please select the index of the avi file you want to work on.\n"); 
	}
	return (listfiles);
}

##
## Print avi files on the console
##
## projectname The directory project name
##
printAviFiles <- function(projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	listfiles <- getAviFiles(projectdir);

	return (listfiles);
}

##
## Print files with the specific extension on the console
##
## projectname The directory project name
## extension The extension of the files
##
printFiles <- function(projectname, extension) {

	projectdir <- paste("projects/", projectname, sep="");
	listfiles <- list.files(projectdir, pattern=extension);

	return (listfiles);
}
