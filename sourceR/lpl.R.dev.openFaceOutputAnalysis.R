## R code for computing head model and residuals motions from OpenFace toolkit output 
## author : S. Rauzy, LPL
## date : 15/09/2017

##
## Extract from the cvs output file of OpenFace toolkit the table containing the pertinent information for our analysis
##
## projectdir : The directory of the project where the cvs output file is
## filename : The file name of cvs OpenFacee output file
##
## return the table containing the pertinent information for our analysis
##
lpl.R.dev.openFaceOutputAnalysis.loadOpenFaceSingleFaceOutput <- function(projectdir, filename) {

	filename <- paste(projectdir, filename, sep="/");

	## Read the cvs output file
	d <- read.table(filename, h=TRUE, sep=",");

	return (d);
}

lpl.R.dev.openFaceOutputAnalysis.meanHeadDimensionInPixels <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	topprojectdir <- projectdir;
	projectdir <- paste(projectdir, "/", software, sep="");

	FOLDER_MODEL <- paste(projectdir, "model", sep="/");
        dhmt <- loadInternalDataFrame(FOLDER_MODEL, "dhmt.txt");

	xn <- dhmt[,2] - dhmt[28, 2];
	yn <- dhmt[,4] - dhmt[28, 4];

	return ((sd(xn) + sd(yn))/2);
}

lpl.R.dev.openFaceOutputAnalysis.meanHead2DimensionInPixels <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	topprojectdir <- projectdir;
	projectdir <- paste(projectdir, "/", software, sep="");

	FOLDER_MODEL <- paste(projectdir, "model", sep="/");
        dhmt <- loadInternalDataFrame(FOLDER_MODEL, "dhmt.txt");

	xmin = min(dhmt[,2]);
	xmax = max(dhmt[,2]);
        ymin = min(dhmt[,4]);	
        ymax = max(dhmt[,4]);

	xymean <- numeric(2);
	xymean[1] = xmax - xmin;
	xymean[2] = ymax - ymin;

	xn <- dhmt[,2] - dhmt[28, 2];
	yn <- dhmt[,4] - dhmt[28, 4];

	return (xymean);
	return ((sd(xn) + sd(yn))/2);
}

##
## Create the head model and the residuals internal facial movements for each landmarks,
## and save the head model and the table of residual dispersions of each landmark in files
##
## projectdir : The directory of the project where the cvs output file is
## d : The data frame containing the pertinent information for our analysis
##
## return the data frame containing residuals, angles, ...
##
lpl.R.dev.openFaceOutputAnalysis.createHeadModelAndResiduals <- function(projectdir, d) {

	dd <- lpl.R.dev.openFaceOutputAnalysis.formPrimaryData(projectdir, d);

	cat("Create the head model from the filtered data (data cleaning)...\n");

	irc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(dd, "inverse", "XYZ");

	d3D <- lpl.R.dev.openFaceOutputAnalysis.createRotationCorrected3DVector(d, irc);

	d3Ds <- subset(d3D, dd$label == "Y");
	dhmt <- lpl.R.dev.openFaceOutputAnalysis.createHeadModel(d3Ds, FALSE);

	FOLDER_MODEL <- paste(projectdir, "model", sep="/");
	if (!file.exists(FOLDER_MODEL)) {
			dir.create(FOLDER_MODEL);
	}
	## Save the head model data frame
	write.table(dhmt, paste(FOLDER_MODEL, "dhmt.txt", sep="/"), sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE);

	dd <- lpl.R.dev.openFaceOutputAnalysis.addProjectedResidualsEstimate(dd, d, irc, dhmt);

	dd <- lpl.R.dev.faceOutputAnalysis.addProjectedErrorAmplitudes(dd, irc);

	FOLDER_TABLES <- paste(projectdir, "tables", sep="/");
	if (!file.exists(FOLDER_TABLES)) {
		dir.create(FOLDER_TABLES);
	}
	cat("Save the internal facial movements file in tables/itdcomplete.txt ...\n");
	saveInternalDataFrame(dd, FOLDER_TABLES, "itdcomplete.txt");

	cat("Compute the standard dispersions for each landmarks...\n");
	ds <- subset(dd, dd$label == "Y");
	dfsd <- lpl.R.dev.faceOutputAnalysis.computeProjectedResidualsStandardDispersion(ds);
	saveInternalDataFrame(dfsd, FOLDER_MODEL, "rrsdt.txt");

	return (dd);
}

##
## Create the projected residuals
##
## data : The data frame containing the S factor value and the angles
## irc : The inverse rotation coefficients 
## dhmt : The direct head model table
##
## return the original data augmented from the residual estimates
##
lpl.R.dev.openFaceOutputAnalysis.addProjectedResidualsEstimate <- function(data, d, irc, dhmt) {

	number_of_points <- 68;
	
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);
	
	cxx <- irc$cxx;
	cxy <- irc$cxy;
	cxz <- irc$cxz;

	cyx <- irc$cyx;
	cyy <- irc$cyy;
	cyz <- irc$cyz;

	czx <- irc$czx;
	czy <- irc$czy;
	czz <- irc$czz;
	
	for (j in c(1:number_of_points)) {
		
		listpointX[[j]] <- numeric(nrow(d));
		listpointY[[j]] <- numeric(nrow(d));

		ix <- which(colnames(d) == paste("X_", (j-1), sep=""));
		xobs <- d[ ,ix] - d$pose_Tx;
		
		iy <- which(colnames(d) == paste("Y_", (j-1), sep=""));
		yobs <- -d[ ,iy] + d$pose_Ty;
	
		Xobs = (cxx - cxz*czx/czz)*xobs + (cxy - cxz*czy/czz)*yobs;
		Yobs = (cyx - cyz*czx/czz)*xobs + (cyy - cyz*czy/czz)*yobs;

		listpointX[[j]] = round(Xobs - dhmt$X[j] + cxz/czz*dhmt$Z[j], 4);
		listpointY[[j]] = round(Yobs - dhmt$Y[j] + cyz/czz*dhmt$Z[j], 4);

		if (j == 1) {
			dd <- data.frame(listpointX[[j]], listpointY[[j]]);
		} else {
			dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		}
		colnames(dd)[ncol(dd)-1] <- paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
	}

	dd <- data.frame(data, dd);

	return (dd);
}

##
## Estimate the parameters of the head model for all the landmarks
##
## d3D : The data frame containing the rotation corrected residuals on X, Y and Z axis (rP..x, rP..y, rP..z)
## removeOutliers : Put this flag to TRUE to remove outliers
##
## return the data frame containing the head model parameters
##
lpl.R.dev.openFaceOutputAnalysis.createHeadModel <- function(d3D, removeOutliers) {

	number_of_points <- 68;
	cat(paste(number_of_points, "landmarks:\n"));   
	df <- NULL;
	for (i in c(1:number_of_points)) {
		df <- rbind(df, lpl.R.dev.openFaceOutputAnalysis.estimateDirectXYZ(d3D, i, removeOutliers));
	}
	cat("\n");   
	
	return (df);
}

##
## Estimate the parameter of the head model for a given landmark
##
## d3D : The data frame containing the rotation corrected vectors on X, Y and Z axis (vP..x, vP..y, vP..z)
## i : The index of the landmark
## removeOutliers : Put this flag to TRUE to remove outliers
##
## return the data frame with direct head model parameters
##
lpl.R.dev.openFaceOutputAnalysis.estimateDirectXYZ <- function(d3D, i, removeOutliers) {

	column_index_x <- which(colnames(d3D) == paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(i), sep=""), "x", sep=""));
	df <- d3D;
	if (removeOutliers) {
		df <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(df, column_index_x, 3);
	}
	lm <- lm(df[, column_index_x] ~ 1);
	X <- lm$coefficients[1];
	seX <- coef(summary(lm))[, 2][1];

	column_index_y <- which(colnames(d3D) == paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(i), sep=""), "y", sep=""));
	df <- d3D;
	if (removeOutliers) {
		df <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(df, column_index_y, 3);
	}
	lm <- lm(df[, column_index_y] ~ 1);
	Y <- lm$coefficients[1];
	seY <- coef(summary(lm))[, 2][1];

	column_index_z <- which(colnames(d3D) == paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(i), sep=""), "z", sep=""));
	df <- d3D;
	if (removeOutliers) {
		df <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(df, column_index_z, 3);
	}
	lm <- lm(df[, column_index_z] ~ 1);
	Z <- lm$coefficients[1];
	seZ <- coef(summary(lm))[, 2][1];

	pointname <- paste("P", lpl.R.dev.faceOutputAnalysis.getStringNumber(i), sep="");
	cat(".");     
	df <- data.frame(pointname, X, seX, Y, seY, Z, seZ);
	names(df) <- c("point", "X", "seX", "Y", "seY", "Z", "seZ");
	rownames(df) <- NULL;

	return (df);
}

lpl.R.dev.openFaceOutputAnalysis.createHeadModel1 <- function(d, dd) {

	line_number <- nrow(d);

	dd <- lpl.R.dev.openFaceOutputAnalysis.formPrimaryData(projectname, d);	

	number_of_points <- 68;
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);
	listpointZ <- vector(mode="list", length=number_of_points);

	## Remove the mean X, Y, Z (head model) from the residuals
	for (j in c(1:number_of_points)) {
		ix <- which(colnames(dd) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
		dd[ , ix] <-  round((dd[ , ix] - dhmt$X[j]), 4);
		
		iy <- which(colnames(dd) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));
		dd[ , iy] <- round((dd[ , iy] - dhmt$Y[j]), 4);

		iz <- which(colnames(dd) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "z", sep=""));
		dd[ , iz] <- round((dd[ , iz] - dhmt$Z[j]), 4);
		ds <- lpl.R.dev.openFaceOutputAnalysis.filterOnConfidenceLevel(dd, 0.8);
	}
	FOLDER_TABLES <- paste(projectdir, "tables", sep="/");
	if (!file.exists(FOLDER_TABLES)) {
		dir.create(FOLDER_TABLES);
	}
	cat("Save the internal facial movements file in table/itdcomplete.txt ...\n");
	saveInternalDataFrame(dd, FOLDER_TABLES, "itdcomplete.txt");

	cat("Compute the standard dispersions for each landmarks...\n");
	ds <- lpl.R.dev.openFaceOutputAnalysis.filterOnConfidenceLevel(dd, 0.8);
	dfsd <- lpl.R.dev.faceOutputAnalysis.computeProjectedResidualsStandardDispersion(ds);
	saveInternalDataFrame(dfsd, FOLDER_MODEL, "rrsdt.txt");

	return (dd);
}


lpl.R.dev.openFaceOutputAnalysis.createRotationCorrected3DVector <- function(d, irc) {

	number_of_points <- 68;
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);
	listpointZ <- vector(mode="list", length=number_of_points);

	d3D <- NULL;

	line_number <- nrow(d);

	## The X and Y and Z vectors for the 68 points
	for (j in c(1:number_of_points)) {
		listpointX[[j]] <- numeric(line_number);
		listpointY[[j]] <- numeric(line_number);
		listpointZ[[j]] <- numeric(line_number);

		ix <- which(colnames(d) == paste("X_", (j-1), sep=""));
		difx <- d[ ,ix] - d$pose_Tx;
		
		iy <- which(colnames(d) == paste("Y_", (j-1), sep=""));
		dify <- -d[ ,iy] + d$pose_Ty;
		
	       	iz <- which(colnames(d) == paste("Z_", (j-1), sep=""));
		difz <- -d[ ,iz] + d$pose_Tz;

		##listpointX[[j]] <- round((irc$cxx*difx + irc$cxy*dify + irc$cxz*difz)/dd$S, 4);
	       	##listpointY[[j]] <- round((irc$cyx*difx + irc$cyy*dify + irc$cyz*difz)/dd$S, 4);	
		listpointX[[j]] <- round((irc$cxx*difx + irc$cxy*dify + irc$cxz*difz), 4);
	       	listpointY[[j]] <- round((irc$cyx*difx + irc$cyy*dify + irc$cyz*difz), 4);	
		listpointZ[[j]] <- round((irc$czx*difx + irc$czy*dify + irc$czz*difz), 4);

		if (is.null(d3D)) {
			d3D <- data.frame(listpointX[[j]], listpointY[[j]], listpointZ[[j]]);
		} else {
			d3D <- data.frame(d3D, listpointX[[j]], listpointY[[j]], listpointZ[[j]]);
		}
		colnames(d3D)[ncol(d3D)-2] <-  paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(d3D)[ncol(d3D)-1] <-  paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
		colnames(d3D)[ncol(d3D)] <-  paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "z", sep="");	
	}

	return (d3D);
}

lpl.R.dev.openFaceOutputAnalysis.checkLine <- function(line) {

	if (is.na(line[11]) || is.na(line[12])) {
		return (FALSE);
	}
	return (TRUE);
}

lpl.R.dev.openFaceOutputAnalysis.changeLine <- function(line) {

	for (i in c(3:length(line))) {
		line[i] <- 0;
	}
	return (line);
}

lpl.R.dev.openFaceOutputAnalysis.changeProblematicLine <- function(d) {

	for (i in c(1:nrow(d))) {
		print(i);
		if (!lpl.R.dev.openFaceOutputAnalysis.checkLine(d[i,])) {
			line <- lpl.R.dev.openFaceOutputAnalysis.changeLine(d[i,]);
			
			for (j in c(1:length(line))) {
				d[i,j] <- line[j];
			}
		}
	}
	return (d);
}

lpl.R.dev.openFaceOutputAnalysis.computeFocalLengthAndOpticalCenter <- function(d) {

	v_0 <- d$X_0/d$Z_0;
	lm <- lm(d$x_0 ~ v_0);
	cx <- round(lm$coefficients[1], 4);
	fx <- round(lm$coefficients[2], 4);

	v_0 <- d$Y_0/d$Z_0;
	lm <- lm(d$y_0 ~ v_0);
	cy <- round(lm$coefficients[1], 4);
	fy <- round(lm$coefficients[2], 4);

	# Create a filter for 3-sigma outliers on the TZ variable
	cis = which(colnames(d) == "pose_Tz");
	ds <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(d, cis, 3);
	## Remove measurements with low confidence level
	dc <- lpl.R.dev.openFaceOutputAnalysis.filterOnConfidenceLevel(ds, 0.8);
	
	## Compute the mean pose_Tz
	meanTz <- mean(dc$pose_Tz);

	dd <- data.frame(fx, fy, cx, cy, meanTz);
	colnames(dd) <-  c("fx", "fy", "cx", "cy", "mZ");

	return (dd);
}

lpl.R.dev.openFaceOutputAnalysis.computeCameraProjection <- function(dfloc, X, Y, Z) {

	lcm <- vector(mode="list", length=2);
	lcm[[1]] <- dfloc$fx*X/Z + dfloc$cx;
	lcm[[2]] <- dfloc$fy*Y/Z + dfloc$cy;

	return (lcm);
}

lpl.R.dev.openFaceOutputAnalysis.convertInPixel <- function(dfloc, values, axis) {

	nv <- numeric(length(values));
	if (axis == "x") {
		nv <- dfloc$fx*values/dfloc$mZ;
	} else {
		nv <- dfloc$fy*values/dfloc$mZ;
	}
	return (nv);
}

lpl.R.dev.openFaceOutputAnalysis.formPrimaryData <- function(projectdir, d) {

	line_number <- nrow(d);

	cat(paste("TREATING THE CSV FILE...\n"));
	cat(paste("Number of frames :", line_number, ", video duration :", lpl.R.dev.faceOutputAnalysis.stomnsString(line_number/25), "...\n"));

	frame <- numeric(line_number);
	time <- numeric(line_number);
	timem <- numeric(line_number);
	times <- numeric(line_number);
	confidence <- numeric(line_number);
	pitch <- numeric(line_number);
	yaw <- numeric(line_number);
	roll <- numeric(line_number);

	time <- round(as.numeric(d$timestamp), 4);

	for (i in c(1:line_number)) {
		frame[i] = i;
		## time in seconds 
		timem[i] <-  round(as.numeric(floor(time[i]/60)), 4);
		## time in minutes and seconds
		times[i] <- round(as.numeric(time[i] - timem[i]*60), 4);
	}

	## Degree to radian conversion coefficient
	alpha = 2*pi/360;

	confidence <-  round(as.numeric(as.character(d$confidence)), 4);
	pitch <- round(as.numeric(as.character(-d$pose_Rx/alpha)), 4);
	yaw <- round(as.numeric(as.character(d$pose_Ry/alpha)), 4);
	roll <- round(as.numeric(as.character(d$pose_Rz/alpha)), 4);

	dfloc <- lpl.R.dev.openFaceOutputAnalysis.computeFocalLengthAndOpticalCenter(d);
	mcp <- lpl.R.dev.openFaceOutputAnalysis.computeCameraProjection(dfloc, d$pose_Tx, d$pose_Ty, d$pose_Tz);

	MOx <- round(mcp[[1]], 4);
	MOy <- round(mcp[[2]], 4);

	dd <- data.frame(frame, time, timem, times, confidence, pitch, yaw, roll, MOx, MOy);

	## Add the vM.. variables
	number_of_points <- 68;
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);

	for (j in c(1:number_of_points)) {

		listpointX[[j]] <- numeric(line_number);
		listpointY[[j]] <- numeric(line_number);

		ix <- which(colnames(d) == paste("x_", (j-1), sep=""));
		x <- d[ ,ix] - MOx;
		
		iy <- which(colnames(d) == paste("y_", (j-1), sep=""));
		y <- d[ ,iy] - MOy;
		listpointX[[j]] <- round(x, 4);
	       	listpointY[[j]] <- round(y, 4);	
		dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		colnames(dd)[ncol(dd)-1] <-  paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-  paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");	
	}

	# Create a filter for 3-sigma outliers on the MOz variable
	cis = which(colnames(d) == "pose_Tz");
	ds <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(d, cis, 3);
	## Remove measurements with low confidence level
	dc <- lpl.R.dev.openFaceOutputAnalysis.filterOnConfidenceLevel(ds, 0.8);
	
	## Compute the mean pose_Tz
	meanTz <- mean(dc$pose_Tz);
	## The scale factor is defined like that
	S <- round(as.numeric(as.character(meanTz/d$pose_Tz)), 4);

	## Form the data frame
	dd <- data.frame(dd, S);

	## Add the clean data flag
	dd <- lpl.R.dev.openFaceOutputAnalysis.addOpenFaceLabelToCleanData(dd);

	return (dd);
}

##
## Create a data frame with the AU (Action Unit) measurements (intensity and detection column) specific to eyebrows 
## d : The data frame containing the pertinent information for our analysis
##
lpl.R.dev.openFaceOutputAnalysis.createEyebrowActionUnitTable <- function(d) {

	line_number <- nrow(d);

	cat(paste("TREATING THE CSV FILE...\n"));
	cat(paste("Number of frames :", line_number, ", video duration :", lpl.R.dev.faceOutputAnalysis.stomnsString(line_number/25), "...\n"));

	AU01 <- numeric(line_number);
	AU01b <- numeric(line_number);
	AU02 <- numeric(line_number);
	AU02b <- numeric(line_number);
	AU04 <- numeric(line_number);
	AU04b <- numeric(line_number);

	i <- which(colnames(d) == "AU01_r");
	AU01 <- d[ ,i];
	i <- which(colnames(d) == "AU02_r");
	AU02 <- d[ ,i];
	i <- which(colnames(d) == "AU04_r");
	AU04 <- d[ ,i];
	i <- which(colnames(d) == "AU01_c");
	AU01b <- d[ ,i];
	i <- which(colnames(d) == "AU02_c");
	AU02b <- d[ ,i];
	i <- which(colnames(d) == "AU04_c");
	AU04b <- d[ ,i];

	dd <- data.frame(AU01, AU02, AU04, AU01b, AU02b, AU04b);

	return (dd);
}

##
## Create a data frame with the AU (Action Unit) measurements (intensity and detection column) 
## d : The data frame containing the pertinent information for our analysis
##
lpl.R.dev.openFaceOutputAnalysis.createActionUnitTable <- function(d) {

	## Retrieve the list of column indexes of the AU
	list_column_indexes_AU <- grep('AU.*', colnames(d), value = FALSE);
	## The number of columns
	number_of_columns <- length(list_column_indexes_AU);

	j = 1;
	df <- data.frame(d[, list_column_indexes_AU[j]]);
	colnames(df)[j] <-  colnames(d)[list_column_indexes_AU[j]];

	for (j in c(2:number_of_columns)) {
		df <- data.frame(df, d[, list_column_indexes_AU[j]]);
		colnames(df)[j] <-  colnames(d)[list_column_indexes_AU[j]];
	}

	list_column_indexes_AU_r <- grep('AU.*_r', colnames(d), value = FALSE);
	cat(paste("Number of AU with intensity measurements :", length(list_column_indexes_AU_r), "\n"));
	list_column_indexes_AU_c <- grep('AU.*_c', colnames(d), value = FALSE);
	cat(paste("Number of AU with detection measurements :", length(list_column_indexes_AU_c), "\n"));

	return (df);
}

lpl.R.dev.openFaceOutputAnalysis.checkData <- function(projectname, d) {

	line_number <- nrow(d);

	filename <- paste(projectname, "_ofo.csv", sep="");
	cat(paste("TREATING CSV FILE", filename, "...\n"));
	cat(paste("Number of frames :", line_number, ", video duration :", lpl.R.dev.faceOutputAnalysis.stomnsString(line_number/25), "...\n"));

	frame <- numeric(line_number);
	time <- numeric(line_number);
	timem <- numeric(line_number);
	times <- numeric(line_number);
	confidence <- numeric(line_number);
	pitch <- numeric(line_number);
	yaw <- numeric(line_number);
	roll <- numeric(line_number);

	time <- round(as.numeric(d$timestamp), 4);

	for (i in c(1:line_number)) {
		frame[i] = i;
		## time in seconds 
		timem[i] <-  round(as.numeric(floor(time[i]/60)), 4);
		## time in minutes and seconds
		times[i] <- round(as.numeric(time[i] - timem[i]*60), 4);
	}

	## Degree to radian conversion coefficient
	alpha = 2*pi/360;

	confidence <-  round(as.numeric(as.character(d$confidence)), 4);
	pitch <- round(as.numeric(as.character(-d$pose_Rx/alpha)), 4);
	yaw <- round(as.numeric(as.character(d$pose_Ry/alpha)), 4);
	roll <- round(as.numeric(as.character(d$pose_Rz/alpha)), 4);

	MOx <- round(as.numeric(as.character(d$pose_Tx)), 4);
	MOy <- round(as.numeric(as.character(-d$pose_Ty)), 4);
	MOz <- round(as.numeric(as.character(-d$pose_Tz)), 4);

	dd <- data.frame(frame, time, timem, times, confidence, pitch, yaw, roll, MOx, MOy, MOz);
	# Create a filter for 3-sigma outliers on the MOz variable
	cis = which(colnames(dd) == "MOz");
	ds <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(dd, cis, 3);
	## Remove measurements with low confidence level
	dc <- lpl.R.dev.openFaceOutputAnalysis.filterOnConfidenceLevel(ds, 0.8);
	
	## Compute the mean MOz
	meanMOz <- mean(dc$MOz);
	## The scale factor is defined like that
	S <- round(as.numeric(as.character(meanMOz/MOz)), 4);

	dd <- data.frame(dd, S);

	## Add the vM.. variables
	number_of_points <- 68;
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);

	r <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(dd, "direct", "XYZ");

	for (j in c(1:number_of_points)) {

		listpointX[[j]] <- numeric(line_number);
		listpointY[[j]] <- numeric(line_number);

		ix <- which(colnames(d) == paste("X_", (j-1), sep=""));
		difx <- d[ ,ix] - MOx;
		
		iy <- which(colnames(d) == paste("Y_", (j-1), sep=""));
		dify <- -d[ ,iy] - MOy;

		iz <- which(colnames(d) == paste("Z_", (j-1), sep=""));
		difz <- -d[ ,iz] - MOz;

		X <- round((r$cxx*difx + r$cxy*dify + r$cxz*difz) + MOx, 4);
	       	Y <- round((r$cyx*difx + r$cyy*dify + r$cyz*difz) + MOy, 4);	
		Z <- round((r$czx*difx + r$czy*dify + r$czz*difz) + MOz, 4);

		dd <- data.frame(dd, X, Y, Z);
		colnames(dd)[ncol(dd)-2] <-  paste(paste("M", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(dd)[ncol(dd)-1] <-  paste(paste("M", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
		colnames(dd)[ncol(dd)] <-  paste(paste("M", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "z", sep="");	
	}

	return (dd);
}

##
## Add a column named label, with various code (N for NA data, S for too large or small face 
## dimension, P for too large positions for P05, L for too large angles, and Y for times which
## will be used to calibrate the head model
## data : The data frame
##
lpl.R.dev.openFaceOutputAnalysis.addOpenFaceLabelToCleanData <- function(data) {

	cat(paste("Total number of frames :", nrow(data), "\n"));

	cat("Annotate frame with low or medium confidence (c < 0.8) value with label N...\n");
	d1 <- lpl.R.dev.openFaceOutputAnalysis.addLabelForLowConfidenceFrame(data);	
	cat(paste("Number of N frames :", sum(d1$label == "N"), "\n"));

	cat("Annotate frame with too large or small face size (S factor) with label S...\n");
	d2 <- lpl.R.dev.faceOutputAnalysis.addLabelSFactorOutlier(d1, d1$S, 2);
	cat(paste("Number of S frames :", sum(d2$label == "S"), "\n"));

	cat("Annotate frame with too large or small MO X and Y coordinates with label P...\n");
	d3 <- lpl.R.dev.faceOutputAnalysis.addLabelOutlierHeadPosition(d2, 3);
	cat(paste("Number of P frames :", sum(d3$label == "P"), "\n"));

	maximal_angle = 15;
	cat(paste("Annotate with too large or small angle value around the mean (", maximal_angle, "degree) with label L...\n"));
	d4 <- lpl.R.dev.faceOutputAnalysis.addLabelForLargeAngle(d3, -maximal_angle, maximal_angle, -maximal_angle, maximal_angle, -maximal_angle, maximal_angle);	
	cat(paste("Number of L frames :", sum(d4$label == "L"), "\n"));

	cat(paste("Number of clean frames (annotated with Y) :", sum(d4$label == "Y"), "\n"));

	## Add the factor information to the label column
	d4$label <- as.factor(d4$label);

	return (d4);
}

##
## Create the label column with N label for data with low confidence values and Y for other points
## data : The data frame
##
lpl.R.dev.openFaceOutputAnalysis.addLabelForLowConfidenceFrame <- function(data) {

	label <- character(nrow(data));
	for (i in c(1:nrow(data))) {
		if (data$confidence[i] < 0.8) {
			label[i]<- "N";
		} else {
			label[i] <- "Y";
		}
	}

	return (data.frame(data, label, stringsAsFactors=FALSE));
}

##
## Filter the data of measurements below the confidence level threshold
##  data :The data frame
##  confidencelevelthreshold :The confidence level threshold
##
lpl.R.dev.openFaceOutputAnalysis.filterOnConfidenceLevel <- function(data, confidencelevelthreshold) {

	return (subset(data, data$confidence > confidencelevelthreshold));
}

##
##  Create a filter for data of measurements below the confidence level threshold
##  data :The data frame
##  confidencelevelthreshold : The confidence level threshold
##
##  The filter, column of LOGIGAL with TRUE value for row to be filtered and FALSE otherwise
##
lpl.R.dev.openFaceOutputAnalysis.createConfidenceLevelFilter <- function(data, confidencelevelthreshold) {

	f <- logical(nrow(data));
	f <- (data$confidence <= confidencelevelthreshold)
	
	return (f);
}

##
## Prepare the OpenFace output (transformed in data frame) by adding a label
## column allowing to clean the data before applying the head model operation
## iodf : The OpenFace output (transformed in data frame)
##
lpl.R.dev.openFaceOutputAnalysis.prepareDataFrameOpenFaceOutput <- function(iodf) {

	iodfc <- lpl.R.dev.faceOutputAnalysis.centerOn(28, iodf);
	iodfc <- lpl.R.dev.faceOutputAnalysis.invertYAxis(iodfc);
	iodfc <- lpl.R.dev.openFaceOutputAnalysis.addOpenFaceLabelToCleanData(iodfc);

	return (iodfc);
}

##
## Create the data filter for the wavelet analysis.
## df : The data frame containing residuals
##
## return the the filter (a LOGICAL column with TRUE value for rows to discard and FALSE for rows to conserve) 
##
lpl.dev.openFaceOutputAnalysis.createFilterForWaveletAnalysis <- function(df) {

	fc <- lpl.R.dev.openFaceOutputAnalysis.createConfidenceLevelFilter(df, 0.8);
	## Create a filter for 3-sigma outliers on the S factor 
	cis = which(colnames(df) == "S");
	fs <- lpl.R.dev.faceOutputAnalysis.createFilterOutliersForColumn(df, cis, 3);
	## Merge the two filters
	f1 <- lpl.R.dev.faceOutputAnalysis.mergeFilter(fs, fc);
	## Create a filter for 3-sigma outliers on the projection amplitude error variable on y
	cieay = which(colnames(df) == "eaY");
	feay <- lpl.R.dev.faceOutputAnalysis.createFilterOutliersForColumn(df, cieay, 3);
	## Merge the filters
	f <- lpl.R.dev.faceOutputAnalysis.mergeFilter(f1, feay);

	return (f);
}





lpl.R.dev.openFaceOutputAnalysis.addProjectedResidualsEstimate1 <- function(dd, d3D, dhmt) {

	number_of_points <- 68;
	line_number <- nrow(dd);
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);
	listpointZ <- vector(mode="list", length=number_of_points);

	## Remove the mean X, Y, Z (head model) from the residuals
	for (j in c(1:number_of_points)) {
		
		listpointX[[j]] <- numeric(line_number);
		listpointY[[j]] <- numeric(line_number);
		listpointZ[[j]] <- numeric(line_number);

		ix <- which(colnames(d3D) == paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
		listpointX[[j]]  <-  round((d3D[ , ix] - dhmt$X[j]), 4);
		
		iy <- which(colnames(d3D) == paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));
		listpointY[[j]] <- round((d3D[ , iy] - dhmt$Y[j]), 4);

		iz <- which(colnames(d3D) == paste(paste("vP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "z", sep=""));
		listpointZ[[j]] <- round((d3D[ , iz] - dhmt$Z[j]), 4);

		dd <- data.frame(dd, listpointX[[j]], listpointY[[j]], listpointZ[[j]]);
		colnames(dd)[ncol(dd)-2] <-  paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(dd)[ncol(dd)-1] <-  paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
		colnames(dd)[ncol(dd)] <-  paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "z", sep="");	
	}

	return (dd);
}	

##
## Create the projected residuals  
## data : The data frame containing the S factor value and the angles
## irc : The inverse rotation coefficients 
## dhmt : The dirext head model table
##
lpl.R.dev.openFaceOutputAnalysis.addProjectedResidualsEstimate2 <- function(data, irc, dhmt) {

	number_of_points <- 68;
	
	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);
	
	cxx <- irc$cxx;
	cxy <- irc$cxy;
	cxz <- irc$cxz;

	cyx <- irc$cyx;
	cyy <- irc$cyy;
	cyz <- irc$cyz;

	czx <- irc$czx;
	czy <- irc$czy;
	czz <- irc$czz;
	
	for (j in c(1:number_of_points)) {
		
		listpointX[[j]] <- numeric(nrow(data));
		listpointY[[j]] <- numeric(nrow(data));

		ix <- which(colnames(data) == paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
		xobs <- data[ ,ix]/data$S;
		
		iy <- which(colnames(d) ==paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));
		yobs <- -data[ ,iy]/data$S;
	
		Xobs = (cxx - cxz*czx/czz)*xobs + (cxy - cxz*czy/czz)*yobs;
		Yobs = (cyx - cyz*czx/czz)*xobs + (cyy - cyz*czy/czz)*yobs;

		listpointX[[j]] = round(Xobs - dhmt$X[j] + cxz/czz*dhmt$Z[j], 4);
		listpointY[[j]] = round(Yobs - dhmt$Y[j] + cyz/czz*dhmt$Z[j], 4);

		if (j == 1) {
			dd <- data.frame(listpointX[[j]], listpointY[[j]]);
		} else {
			dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		}
		colnames(dd)[ncol(dd)-1] <- paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
	}

	dd <- data.frame(data, dd);

	return (dd);
}

