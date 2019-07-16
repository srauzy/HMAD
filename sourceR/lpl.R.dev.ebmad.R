## lpl.R.dev.ebmad.R, creator S. Rauzy, LPL, 12/02/2017

##
## Create the automatic annotation of eyebrow raise and frown motion from the residual motion of the landmarks
## software : The software used for extracting the landmark residuals 
## projectname : The project name
## df : The data frame containing the residual facial motions
##
lpl.R.dev.ebmad.runEBMADetection <- function(software, projectname, df) {

	## Create step by step the automatic annotation of eyebrow raise and frown motion
	ebmad <- lpl.R.dev.ebmad.computeEBMAD(software, projectname, df);

	return (ebmad);
}

##
## Create error function from the error amplitudes in X and Y, the coefficients and their standard dispersions
## df : The data frame containing residual facial movements
## dfc : The data frame containing the coefficients for each residual on X and Y axis
## dfsd : The data frame containing standard deviation for each residual on X and Y axis
##
errorFunctionColumn <- function(df, dfc, dfsd) {

	wt <- 0;
	func <- numeric(nrow(df));
	for (i in c(1:nrow(dfc))) {
		landmark <- as.character(dfc$landmark[i]);
		dfsd_landmark_row_index <- which(dfsd$point == landmark);
		if (dfc$coefx[i] != 0) {
			sd <- dfsd$sdx[dfsd_landmark_row_index];
			wt <- wt + 1/sd;
			func <- func + dfc$coefx[i]*df$eaX;
		}
		if (dfc$coefy[i] != 0) {
			sd <- dfsd$sdy[dfsd_landmark_row_index];
			wt <- wt + 1/sd;
			func <- func + dfc$coefy[i]/sd*df$eaY;
		}
	}	
	func <- func/wt;
	
	return (func);
}


##
## Create the raise function from the residuals and their standard dispersion
## df : The data frame containing residual motion of face points
## sdrx : The data frame containing standard deviation of residuals on the x axis
## sdry : The data frame containing standard deviation of residuals on the y axis
##
raiseAdhocFunctionColumn <- function(df, sdrx, sdry) {

	wt <- 0.5/sdrx$sr[1] + 0.5/sdrx$sr[4] + 0.5/sdrx$sr[7] + 0.5/sdrx$sr[8]+ 0.5/sdrx$sr[9]+ 0.5/sdrx$sr[10]+ 0.5/sdry$sr[1] + 1/sdry$sr[2] + 1/sdry$sr[3] + 0.5/sdry$sr[4] + 0.5/sdry$sr[5]+ 0.5/sdry$sr[6]+ 1/sdry$sr[8]+ 1/sdry$sr[9];
	raiseFunction <- ((0.5*df$rP01x/sdrx$sr[1]-0.5*df$rP04x/sdrx$sr[4]+0.5*df$rP07x/sdrx$sr[7]+0.5*df$rP08x/sdrx$sr[8]-0.5*df$rP09x/sdrx$sr[9]-0.5*df$rP10x/sdrx$sr[10]) + (0.5*df$rP01y/sdry$sr[1]+df$rP02y/sdry$sr[2]+df$rP03y/sdry$sr[3]+0.5*df$rP04y/sdry$sr[4])-0.5*df$rP05y/sdry$sr[5]-0.5*df$rP06y/sdry$sr[6]-df$rP08y/sdry$sr[8]-df$rP09y/sdry$sr[9])/wt; 
	
	return (raiseFunction);
}

##
## Create the frown function from the residuals and their standard dispersion
## df : The data frame containing residual motion of face points
## sdrx : The data frame containing standard deviation of residuals on the x axis
## sdry : The data frame containing standard deviation of residuals on the y axis
##
frownAdhocFunctionColumn <- function(df, sdrx, sdry) {

	wt <- 0.5/sdrx$sr[2] + 0.5/sdrx$sr[3] + 0.5/sdrx$sr[7] + 0.5/sdrx$sr[8]+ 0.5/sdrx$sr[9]+ 0.5/sdrx$sr[10]+ 0.5/sdry$sr[2] + 0.5/sdry$sr[3] + 0.5/sdry$sr[5]+ 0.5/sdry$sr[6]+ 1/sdry$sr[8]+ 1/sdry$sr[9]+ 0.5/sdry$sr[11]+ 0.5/sdry$sr[12];
	frownFunction <- ((0.5*df$rP02x/sdrx$sr[2]-0.5*df$rP03x/sdrx$sr[3]-0.5*df$rP07x/sdrx$sr[7]-0.5*df$rP08x/sdrx$sr[8]+0.5*df$rP09x/sdrx$sr[9]+0.5*df$rP10x/sdrx$sr[10]) + (-0.5*df$rP02y/sdry$sr[2]-0.5*df$rP03y/sdry$sr[3]+0.5*df$rP05y/sdry$sr[5]+0.5*df$rP06y/sdry$sr[6]+df$rP08y/sdry$sr[8]+df$rP09y/sdry$sr[9]-0.5*df$rP11y/sdry$sr[11]-0.5*df$rP12y/sdry$sr[12]))/wt; 
	return (frownFunction);
}

##
## Create the frown error function from the errors amplitude and the standard dispersion of residuals
## eax : The error amplitude on x
## eax : The error amplitude on y
## sdrx : The data frame containing standard deviation of residuals on the x axis
## sdry : The data frame containing standard deviation of residuals on the y axis
##
frownErrorFunctionColumn <- function(eaX, eaY, sdrx, sdry, residualtype) {

	wt <- 1/sdrx$sr[2] + 1/sdrx$sr[3] + 1/sdry$sr[2] + 1/sdry$sr[3] + 1/sdry$sr[5];
	frownErrorFunction <- (eaX/sdrx$sr[2] + eaX/sdrx$sr[3] + eaY/sdry$sr[2] + eaY/sdry$sr[3] + eaY/sdry$sr[5])/wt; 
	return (frownErrorFunction);
}

##
## Create the raise error function from the errors amplitude and the standard dispersion of residuals
## eax : The error amplitude on x
## eax : The error amplitude on y
## sdrx : The data frame containing standard deviation of residuals on the x axis
## sdry : The data frame containing standard deviation of residuals on the y axis
##
raiseErrorFunctionColumn <- function(eaX, eaY, sdrx, sdry, residualtype) {

	wt <- 1/sdrx$sr[1] + 1/sdrx$sr[4] + 1/sdry$sr[1] + 1/sdry$sr[2] + 1/sdry$sr[3] + 1/sdry$sr[4] + 1/sdry$sr[5];
	raiseErrorFunction <- ((eaX/sdrx$sr[1]+eaX/sdrx$sr[4]) + (eaY/sdry$sr[1]+eaY/sdry$sr[2]+eaY/sdry$sr[3]+eaY/sdry$sr[4]+eaY/sdry$sr[5]))/wt; 
	return (raiseErrorFunction);
}

##
## Create step by step (first the eyebrow raise automatic annotation and second eyebrow frown one) the automatic
## annotations. Intermediate files are saved, so if the program is stopped when running, the computation will
## start for the next run after the last file created.
##
## software : The software used to created the output
## projectname : The name of the project
## df : The data frame containing residual motion of the facial landmarks
##
lpl.R.dev.ebmad.computeEBMAD <- function(software, projectname, df) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdir <- paste("projects/", projectname, sep="");
	projectdir <- paste(projectdir, "/", software, sep="");

	## Define the folders	
	FOLDER_MODEL <- paste(projectdir, "model", sep="/");

	ebmaddir <- paste(projectdir, "/", "EBMAD", sep="");
	if (!file.exists(ebmaddir)) {
		dir.create(ebmaddir);
	}
	FOLDER_TABLES <- paste(ebmaddir, "tables", sep="/");
	if (!file.exists(FOLDER_TABLES)) {
		dir.create(FOLDER_TABLES);
	}

	qualitythreshold = 0.3;

	## If the result file for automatic annotation of eyebrow raise has already been created, just load it
	if (fileExists(FOLDER_TABLES, "ebrmad.txt")) {
		cat("Automatic annotations of eyebrow raise motion previously done...\n");
		ebrmad <- loadInternalDataFrame(FOLDER_TABLES, "ebrmad.txt");
	} else {
		if (fileExists(FOLDER_TABLES, "dfwcebr.txt")) {
			df1Dqf <- loadInternalDataFrame(FOLDER_TABLES, "df1Dqfebr.txt");
			dfwc <- loadInternalDataFrame(FOLDER_TABLES, "dfwcebr.txt");
		} else {
			cat("Create wavelet 2D coefficients for raise function...\n");

			## Load the raise function coefficients model	
			dfrf <- lpl.R.dev.hmad.loadFunctionModel(software, "ebraise_model.txt");
			## Load the residual standard dispersion table
			dfrsd <- loadInternalDataFrame(FOLDER_MODEL, "rrsdt.txt");
			## Create the filter to discard problematic frames
			if (software == "INTRAFACE") {
				filter <- lpl.dev.intraFaceOutputAnalysis.createFilterForWaveletAnalysis(df);
			} else if (software == "OPENFACE") {
				filter <- lpl.dev.openFaceOutputAnalysis.createFilterForWaveletAnalysis(df);
			}
			## Create the raise function column from residuals, coefficients, and residual standard dispersions
			raiseFunction <- lpl.R.dev.hmad.functionColumn(df, dfrf, dfrsd);

			minimalscale <- 0.16;

			cat(paste("Minimal duration parameter for eyebrow motion detection =", minimalscale, "seconds", "\n"));
			dfwc <- lpl.R.dev.ebmotion.createWaveletTable(df, raiseFunction, filter, qualitythreshold, minimalscale, 2.4, 0.08, 1);

			df1Dqf <- lpl.R.dev.ebmotion.create1DQualityFunction(df, filter, minimalscale, 1);

			saveInternalDataFrame(df1Dqf, FOLDER_TABLES, "df1Dqfebr.txt");
			saveInternalDataFrame(dfwc, FOLDER_TABLES, "dfwcebr.txt");
		}

		ebrmad <- lpl.R.dev.ebmotion.automaticAnnotationEBM(dfwc, df1Dqf, qualitythreshold, "Raise");
		## Save the result in a file
		saveInternalDataFrame(ebrmad, FOLDER_TABLES, "ebrmad.txt");
	}

	## If the result file for automatic annotation of eyebrow frown has already been created, just load it
	if (fileExists(FOLDER_TABLES, "ebfmad.txt")) {
		cat("Automatic annotations of eyebrow frown motion previously done...\n");
		ebfmad <- loadInternalDataFrame(FOLDER_TABLES, "ebfmad.txt");
	} else {
		if (fileExists(FOLDER_TABLES, "dfwcebf.txt")) {
			df1Dqf <- loadInternalDataFrame(FOLDER_TABLES, "df1Dqfebf.txt");
			dfwc <- loadInternalDataFrame(FOLDER_TABLES, "dfwcebf.txt");
		} else {
			cat("Create wavelet 2D coefficients for frown function...\n");

			## Load the raise function coefficients model	
			dfrf <- lpl.R.dev.hmad.loadFunctionModel(software, "ebfrown_model.txt");
			## Load the residual standard dispersion table
			dfrsd <- loadInternalDataFrame(FOLDER_MODEL, "rrsdt.txt");
			## Create the filter to discard problematic frames
			if (software == "INTRAFACE") {
				filter <- lpl.dev.intraFaceOutputAnalysis.createFilterForWaveletAnalysis(df);
			} else if (software == "OPENFACE") {
				filter <- lpl.dev.openFaceOutputAnalysis.createFilterForWaveletAnalysis(df);
			}
			## Create the raise function column from residuals, coefficients, and residual standard dispersions
			frownFunction <- lpl.R.dev.hmad.functionColumn(df, dfrf, dfrsd);

			minimalscale <- 0.32;

			cat(paste("Minimal duration parameter for eyebrow motion detection =", minimalscale, "seconds", "\n"));
			dfwc <- lpl.R.dev.ebmotion.createWaveletTable(df, frownFunction, filter, qualitythreshold, minimalscale, 3.6, 0.08, 1);

			df1Dqf <- lpl.R.dev.ebmotion.create1DQualityFunction(df, filter, minimalscale, 1);

			saveInternalDataFrame(df1Dqf, FOLDER_TABLES, "df1Dqfebf.txt");
			saveInternalDataFrame(dfwc, FOLDER_TABLES, "dfwcebf.txt");
		}

		ebfmad <- lpl.R.dev.ebmotion.automaticAnnotationEBM(dfwc, df1Dqf, qualitythreshold, "Frown");
		## Save the result in a file
		saveInternalDataFrame(ebfmad, FOLDER_TABLES, "ebfmad.txt");
	}

	if (!fileExists(FOLDER_TABLES, "ebmad.txt")) {
		cat("Merge annotation of eyebrow raise and frown...\n");
		lmr <- subset(ebrmad, ebrmad$class != "N" &  ebrmad$class != "X");
		lmf <- subset(ebfmad, ebfmad$class != "N" &  ebfmad$class != "X");
		ebmad <- mergeLocalMaximaByZscore(lmr, lmf);
		cat("----------------\n");
		cat("Add bad quality intervals (class X)...\n");
		qualitythreshold = 0.3;
		df1Dqf <- loadInternalDataFrame(FOLDER_TABLES, "df1Dqfebr.txt");
		ebmad <- lpl.R.dev.ebmotion.addBadQualityAreaWithZscore(ebmad, df1Dqf, qualitythreshold, "NA");
		cat("Add intervals presumed without eyebrow motion (class N)...\n");
		ebmad <- lpl.R.dev.ebmotion.addNullDetectionAreaWithZscore(ebmad, df1Dqf, "NA");
		saveInternalDataFrame(ebmad, FOLDER_TABLES, "ebmad.txt");
	} else {
		cat("Automatic annotations of eyebrow raise and frown motions previously done...\n");
		ebmad <- loadInternalDataFrame(FOLDER_TABLES, "ebmad.txt");	
	}

	lpl.R.dev.ebmad.createHtmlDescription(ebmad, software, projectname);

	return (ebmad);
}

##
## Create an html document containing the summary of the analysis
## data : The data table as a data frame
## software : The software used for extracting the video output 
## projectname : The project name
##
lpl.R.dev.ebmad.createHtmlDescription <- function(data, software, projectname) {

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
	df <- lpl.R.dev.htmldf.openHtmlHeader(df, paste(projectname, "ebm result"));
	df <- lpl.R.dev.htmldf.closeHtmlHeader(df);

	df <- lpl.R.dev.htmldf.openHtmlBody(df);

	line <- paste("Project name:", projectname, sep=" ");
	
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("Automatic prediction of eyebrows raising and frowning actions from", software, "output:", sep=" ");
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	globaldf <- paste(projectdir, "EBMAD/tables", sep="/");
	line <- paste("Table of results in",  globaldf, ": ebfmad.txt, ebrmad.txt and ebmad.txt", sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);

	globaldf <- paste(projectdir, "result", sep="/");

	line <- paste("Elan editable output format in",  globaldf, ": ebmadraf.eaf", sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("Total duration:",  (nrow(data)/25), "seconds or", (lpl.R.dev.faceOutputAnalysis.stomnsString(nrow(data)/25)), sep=" ");

	df <- lpl.R.dev.htmldf.openHtmlTable(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Class");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of<BR>intervals");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Duration<BR>(in seconds)");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Ratio of<BR>Raise");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- data
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "W");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "X");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "X");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class != "X");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "M");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "A");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "B");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "B");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "C");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "C");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "D");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "D");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
		if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "E");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "E");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "N");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "N");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "AB");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B" | data$class == "C");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "ABC");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B" | data$class == "C" | data$class == "D");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "ABCD");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B" | data$class == "C" | data$class == "D" | data$class == "E");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "ABCDE");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.closeHtmlTag(df);

	descriptionfile <- paste(FOLDER_RESULT, "ebmaddescription.html", sep="/");
	cat(paste("Description of the automatic annotation of eyebrows actions can be found in the file ebmaddescription.html (folder ", FOLDER_RESULT, ")", sep=""));
	cat("\n");

	write.table(df, descriptionfile, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");
}

##
## Create from the ebmad table a table containing for each interval of time the quality
## of the automatic annotation
## ebmad : the ebmad table
##
createTierEBMADQualityClassTable <- function(ebmad) {

	tmin <- numeric(nrow(ebmad));
	tmax <- numeric(nrow(ebmad));
	annotation <- character(nrow(ebmad));

	tmin <- round(ebmad$tmin*1000, 0);
	tmax <- round(ebmad$tmax*1000, 0);
	annotation <- ebmad$class;

	return (data.frame(tmin, tmax, annotation));
}

createTierEBMADTypeTable <- function(ebmad) {

	tmin <- numeric(nrow(ebmad));
	tmax <- numeric(nrow(ebmad));
	annotation <- character(nrow(ebmad));

	tmin <- round(ebmad$tmin*1000, 0);
	tmax <- round(ebmad$tmax*1000, 0);
	for (i in c(1:nrow(ebmad))) {
		if (ebmad$class[i] != "X" & ebmad$class[i] != "N") {
			annotation[i] <- "Raise";
		} else {
			annotation[i] <- "";
		}
	}

	return (data.frame(tmin, tmax, annotation));
}

##
## Get the name of the file (not directory) in the project directory (if all is right, it is the video file) and return it
##
## projectname : The name of the project in the directory "projects"
##
## return the video file name or null in case of problem
##
lpl.R.dev.ebmad.getVideoFileName<- function(projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	listfiles <- list.files(projectdir, pattern=".*\\..*");
	if (length(listfiles) < 1) {
		return (cat(paste("No video file found in project directory", projectdir, "...")));
	} else if (length(listfiles) > 1) {
		return (cat(paste("More than one candidate for the video file in the project directory", projectdir, "...")));
	}
	return (listfiles[1]);
}

##
## Get the name of the file (not directory) in the project directory (if all is right, it is the video file) and return it
##
## projectname : The name of the project in the directory "projects"
##
## return the video file name or null in case of problem
##
lpl.R.dev.ebmad.getVideoFileNameWithRegex <- function(projectname, regex) {

	projectdir <- paste("projects/", projectname, sep="");
	listfiles <- list.files(projectdir, pattern=regex);
	if (length(listfiles) < 1) {
		return (cat(paste("No video file found in project directory", projectdir, "...")));
	} else if (length(listfiles) > 1) {
		return (cat(paste("More than one candidate for the video file in the project directory", projectdir, "...")));
	}
	return (listfiles[1]);
}


##
## Return the ELAN video Mime type in function of the extension of the video file or null if not found
##
## projectname : The name of the project in the directory "projects"
##
lpl.R.dev.ebmad.getElanMimeType<- function(projectname) {

	vfn <- lpl.R.dev.ebmad.getVideoFileName(projectname);
	if (endsWith(vfn, ".mp4")) {
		return ("video/mp4");
	} else if (endsWith(vfn, ".avi")) {
	       return ("video/*");
	} else if (endsWith(vfn, ".mov")) {
	       return ("video/quicktime");
	} else if (endsWith(vfn, ".mpg") || endsWith(vfn, ".mpeg")) {
		return ("video/mpeg");
	}

	cat("Mime Type not found! Please consult the Elan instruction guideline (https://www.mpi.nl/corpus/html/lamus/apa.html) and instantiate manually the VIDEO_MIME_TYPE R variable...");
}	

##
## Return the ELAN video Mime type in function of the extension of the video file or null if not found
##
## projectname : The name of the project in the directory "projects"
##
lpl.R.dev.ebmad.getElanMimeTypeWithRegex<- function(projectname, regex) {

	vfn <- lpl.R.dev.ebmad.getVideoFileNameWithRegex(projectname, regex);
	if (endsWith(vfn, ".mp4")) {
		return ("video/mp4");
	} else if (endsWith(vfn, ".avi")) {
	       return ("video/*");
	} else if (endsWith(vfn, ".mov")) {
	       return ("video/quicktime");
	} else if (endsWith(vfn, ".mpg") || endsWith(vfn, ".mpeg")) {
		return ("video/mpeg");
	}

	cat("Mime Type not found! Please consult the Elan instruction guideline (https://www.mpi.nl/corpus/html/lamus/apa.html) and instantiate manually the VIDEO_MIME_TYPE R variable...");
}	
	
##
## Return the name of the eaf output file name function of the type
##
## type : The type of annotations ("RAISE", "FROWN" or "RAISE_AND_FROWN")
##
lpl.R.dev.ebmad.getDefaultAnnotationFileName<- function(type) {
	
	if (type == "RAISE") {
		return ("ebmadr.eaf");
	} else if (type == "FROWN") {
		return ("ebmadf.eaf");
	} else {
		return ("ebmadraf.eaf");
	}
}

##
## Create the ELAN eaf output file format for editing the automatic annotation of the eyebrows action (in directory
## in projects/projectname/software/result)
##
## ebmadtype : The type of annotations : "RAISE", "FROWN" or "RAISE_AND_FROWN"
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.ebmad.createElanEAFFile<- function(ebmadtype, software, projectname, mediafilename, eaffilename, videomimetype) {

	## Define the folders
	projectdir <- paste("projects/", projectname, sep="");

	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	ebmaddir <- paste(projectdir, "/", "EBMAD", sep="");
	
	## The folder for the tables
	FOLDER_TABLES <- paste(ebmaddir, "tables", sep="/");

	if (ebmadtype == "RAISE") {
		if (!fileExists(FOLDER_TABLES, "ebrmad.txt")) {
			cat("Please first run the lpl.R.dev.ebmad.computeEBMAD routine...");
		} else {
			ebmad <- loadInternalDataFrame(FOLDER_TABLES, "ebrmad.txt");
			lpl.R.dev.ebmad.cawriteElanEAFFile(ebmad, software, projectname, mediafilename, eaffilename, videomimetype); 
		}
	} else if (ebmadtype == "FROWN") {
		if (!fileExists(FOLDER_TABLES, "ebfmad.txt")) {
			cat("Please first run the lpl.R.dev.ebmad.computeEBMAD routine...");
		} else {
			ebmad <- loadInternalDataFrame(FOLDER_TABLES, "ebfmad.txt");
			lpl.R.dev.ebmad.cawriteElanEAFFile(ebmad, software, projectname, mediafilename, eaffilename, videomimetype); 
		}
	} else {
		if (!fileExists(FOLDER_TABLES, "ebmad.txt")) {
			cat("Please first run the lpl.R.dev.ebmad.computeEBMAD routine...");
		} else {
			ebmad <- loadInternalDataFrame(FOLDER_TABLES, "ebmad.txt");
			lpl.R.dev.ebmad.cawriteElanEAFFile(ebmad, software, projectname, mediafilename, eaffilename, videomimetype); 
		}
	}
}

##
## Create and write the ELAN eaf output file format for editing the automatic annotation of the eyebrows action (in directory
## in projects/projectname/software/result)
##
## ebmadtype : The type of annotations : "RAISE", "FROWN" or "RAISE_AND_FROWN"
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.ebmad.cawriteElanEAFFile <- function(ebmad, software, projectname, mediafilename, eaffilename, videomimetype) {

	projectdir <- paste("projects/", projectname, sep="");
	videodir <- projectdir;

	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	FOLDER_RESULT <- paste(projectdir, "result", sep="/"); 
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
	df <- lpl.R.dev.elan.createElanEAFTable(ebmad, MEDIA_URL, RELATIVE_MEDIA_URL, videomimetype);

	eaffile <- paste(FOLDER_RESULT, eaffilename, sep="/");

	write.table(df, eaffile, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");

	cat(paste("You can now open the annotation file", paste(FOLDER_RESULT, eaffilename, sep="/"), "with the ELAN Software.\n"));
}
