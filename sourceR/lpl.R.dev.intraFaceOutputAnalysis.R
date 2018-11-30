## R code for computing head model and residuals motions from Intrace software output 
## author : S. Rauzy, LPL
## date : 15/12/2016

##
## Extract from the cvs output file of Intraface software the table containing the 
## pertinent information for our analysis
##
## projectdir : The directory of the project where the cvs output file is
## filename : The file name of cvs Intraface output file
##
## return the data frame containing the angles, landmark positions, etc
##
lpl.R.dev.intraFaceOutputAnalysis.loadIntrafaceOutput12Points <- function(projectdir, filename) {

	filename <- paste(projectdir, filename, sep="/");

	## Read the cvs output file
	d <- read.table(filename, h=TRUE, sep=";");
	if (ncol(d) == 1) {
		d <- read.table(filename, h=TRUE, sep=",");
		d <- d[-c(1),]
	}

	line_number <- nrow(d);
	frame <- numeric(line_number);
	time <- numeric(line_number);
	timem <- numeric(line_number);
	times <- numeric(line_number);
	confidence <- numeric(line_number);
	pitch <- numeric(line_number);
	yaw <- numeric(line_number);
	roll <- numeric(line_number);

	number_of_points <- 12;
	listpointX <- vector(mode="list", length=12);
	listpointY <- vector(mode="list", length=12);

	## The X and Y vector for the 12 points
	for (j in c(1:number_of_points)) {
		listpointX[[j]] <- numeric(line_number);
		listpointY[[j]] <- numeric(line_number);
	}

	cat(paste("TREATING CSV FILE", filename, "...\n"));
	cat(paste("Number of frames :", line_number, ", video duration :", lpl.R.dev.faceOutputAnalysis.stomnsString(line_number/25), "...\n"));

	starting_time = getTimeInHundredthSecond();

	## Loop on the frames
	for (i in c(1:line_number)) {
		
		frame[i] = i;
		## time in seconds (25 frames per second)
		time[i] = (i-1)*0.04;
		## time in minutes and seconds
		timem[i] <- floor(time[i]/60);
		times[i] <- time[i] - timem[i]*60;

		## X and Y positions of the 12 points
		positions = toString(d[i, 4]);
		## Remove the surrounding brackets 
		positions = substr(positions, 2, nchar(positions)-1);

		## Case with no signal, everything is put to NA
		if (nchar(positions) == 0) {

			for (j in c(1:number_of_points)) {
				listpointX[[j]][i] <- NA;
		       		listpointY[[j]][i] <- NA; 
			}

			## Confidence level for the measurements and head position angles
			confidence[i] <- d[i, 2];
			pitch[i] <- NA;
			yaw[i] <- NA;
			roll[i] <- NA;

		## Case with signal
		} else {
			## get X and Y positions separated by a semicolumn
			positionsXY <-  strsplit(positions, ";");
			## For X and Y get the values separated by a coma
			positionsX <- strsplit(positionsXY[[1]][1], ",");
			positionsY <- strsplit(positionsXY[[1]][2], ",");

			## Loop on the points
			for (j in c(1:number_of_points)) {
				listpointX[[j]][i] <- positionsX[[1]][j];
		       		listpointY[[j]][i] <- positionsY[[1]][j]; 
			}

			## Confidence level for the measurements and head position angles
			confidence[i] <- round(as.numeric(as.character(d[i, 2])), 4);
			pitch[i] <- round(as.numeric(as.character(d[i, 5])), 4);
			yaw[i] <- round(as.numeric(as.character(d[i, 6])), 4);
			roll[i] <- round(as.numeric(as.character(d[i, 7])), 4);
		}
		if (i == 25*60-1) {
			ed <- estimateDurationInSecond(durationInSecond(starting_time, getTimeInHundredthSecond()), (line_number/25/60));
			cat(paste("PREDICTED COMPUTING TIME :", lpl.R.dev.faceOutputAnalysis.stomnsString(ed), "\n"));
		}
		if (floor(i/25/60) == (i/25/60)) {
			cat(paste(round((100*i/line_number), 0), "% "));
		}
				   
	}
	cat("\n");

	
	## Form the data frame
	dd <- data.frame(frame, time, timem, times, confidence, pitch, yaw, roll);

	## Add the X and Y measurements for the 12 points
	for (j in c(1:number_of_points)) {
		dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		colnames(dd)[ncol(dd)-1] <-  paste(paste("fl", j, sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-  paste(paste("fl", j, sep=""), "y", sep="");
	}
	return (dd);
}

##
## Modify the table to present the ouput as the point 5, and the 11 remaining
## points as vectors relative to point 5
##
## d : The data frame containing the Intraface output transformed in data frame
##
## return the modified table
##
lpl.R.dev.intraFaceOutputAnalysis.centerOnP5 <- function(d) {

	## Copy the first 8 columns
	dd <- data.frame(d[, 1:8]);

	dx5 <- as.numeric(as.character(d[, (8+2*5-1)]));
	dy5 <- as.numeric(as.character(d[, (8+2*5)]));

	## Add the X and Y of point 5
	dd <- data.frame(dd, dx5, dy5);
	colnames(dd)[ncol(dd)-1] <-  paste(paste("M", "05", sep=""), "x", sep="");
	colnames(dd)[ncol(dd)] <-  paste(paste("M", "05", sep=""), "y", sep="");

	for (j in c(1:12)) {
		difx <- numeric(nrow(d));
		dify <- numeric(nrow(d));
		if (j != 5) {
			dx <- as.numeric(as.character(d[, (8+2*j-1)]));
			dy <- as.numeric(as.character(d[, (8+2*j)]));
			
			difx <- dx - dx5;
			dify <- dy - dy5;
			dd <- data.frame(dd, difx, dify);
			colnames(dd)[ncol(dd)-1] <-  paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
			colnames(dd)[ncol(dd)] <-  paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
		}
	}
	return (dd);
}

##
## Prepare the Intraface output (transformed in data frame) by adding a label
## column allowing to clean the data before applying the head model operation
##
## iodf : The Intraface output (transformed in data frame)
##
## return the prepared data frame
##
lpl.R.dev.intraFaceOutputAnalysis.prepareDataFrameIntrafaceOutput <- function(iodf) {

	iodfc <- lpl.R.dev.intraFaceOutputAnalysis.centerOnP5(iodf);
	iodfc <- lpl.R.dev.faceOutputAnalysis.invertYAxis(iodfc);
	df <- lpl.R.dev.intraFaceOutputAnalysis.addLabelToCleanData(iodfc);

	return (df);
}

##
## Compute the S factor for each frame without correction from the head rotation effect
##
## data : The data frame containing the relative coordinates of the landmarks
## normalize : The normalization flag (TRUE for a S factor equals 1 in average) 
##
## return the raw S factor column
##
lpl.R.dev.intraFaceOutputAnalysis.computeRawSFactor <- function(data, normalize) {

	S <- numeric(nrow(data));

	jc = 1;
	for (j in c(1:12)) {
		if (j != 5) {
		
			column_index_x = which(colnames(data) == paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
			column_index_y = which(colnames(data) == paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));

			xobs <- abs(data[, column_index_x]);
			yobs <- abs(data[, column_index_y]);

			S <- S + xobs + yobs;

			jc <- jc+1;
		}
	}

	if (normalize) {
		mS <- mean(S, na.rm = T);
		S <- S/mS;
	}
	return (S);
}

##
## Create the data filter for the wavelet analysis
##
## df : The data frame containing residuals
##
## return the the filter (a LOGICAL column with TRUE value for rows to discard and FALSE for rows to conserve) 
##
lpl.dev.intraFaceOutputAnalysis.createFilterForWaveletAnalysis <- function(df) {

	## Create a filter for 3-sigma outliers on the S factor 
	cis = which(colnames(df) == "S");
	fs <- lpl.R.dev.faceOutputAnalysis.createFilterOutliersForColumn(df, cis, 3);
	## Create a filter for 3-sigma outliers on the projection amplitude error variable on y
	cieay = which(colnames(df) == "eaY");
	feay <- lpl.R.dev.faceOutputAnalysis.createFilterOutliersForColumn(df, cieay, 3);
	## Merge the two filters
	f <- lpl.R.dev.faceOutputAnalysis.mergeFilter(fs, feay);

	return (f);
}

##
## Add a column named label, with various code (N for NA data, S for too large or small face 
## dimension, P for too large positions for P05, L for too large angles, and Y for times which
## will be used to calibrate the head model
##
## data : The data frame
##
## return the data frame augmented of the label column
##
lpl.R.dev.intraFaceOutputAnalysis.addLabelToCleanData <- function(data) {

	cat(paste("Total number of frames :", nrow(data), "\n"));

	cat("Annotate frame with NA detection value with label N...\n");
	d1 <- lpl.R.dev.intraFaceOutputAnalysis.addLabelForNAFrame(data);	
	cat(paste("Number of N frames :", sum(d1$label == "N"), "\n"));

	cat("Annotate frame with too large or small face size (raw S factor) with label S...\n");
	S <- lpl.R.dev.intraFaceOutputAnalysis.computeRawSFactor(d1, TRUE);
	d2 <- lpl.R.dev.faceOutputAnalysis.addLabelSFactorOutlier(d1, S, 2);
	cat(paste("Number of S frames :", sum(d2$label == "S"), "\n"));

	cat("Annotate frame with too large or small P05 X and Y coordinates with label P...\n");
	d3 <- lpl.R.dev.intraFaceOutputAnalysis.addLabelOutlierP05(d2, 3);
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
## Create the label column with N label for data with NA values and Y for other points
##
## data : The data frame
##
## return the data frame with the label column
##
lpl.R.dev.intraFaceOutputAnalysis.addLabelForNAFrame <- function(data) {

	label <- character(nrow(data));
	for (i in c(1:nrow(data))) {
		if (is.na(data$pitch[i])) {
			label[i]<- "N";
		} else {
			label[i] <- "Y";
		}
	}

	return (data.frame(data, label, stringsAsFactors=FALSE));
}

##
## Add a label P to P05 x and y coordinates too large
##
## data : The data frame
## threshold : The threshold in sigma for removing outliers
##
## return the data frame with the modified label column
##
lpl.R.dev.intraFaceOutputAnalysis.addLabelOutlierP05 <- function(data, threshold) {

	dx <- as.numeric(as.character(data$P05x));
	xmean <- mean(dx, na.rm=T);
	xsd <- sd(dx, na.rm=T);

	dy <- as.numeric(as.character(data$P05y));
	ymean <- mean(dy, na.rm=T);
	ysd <- sd(dy, na.rm=T);

	for (i in c(1:nrow(data))) {
		if (as.character(data$label[i]) == "Y") {
			if (!is.na(dx[i])) {
				if (((dx[i]-xmean)/xsd) < -threshold | (threshold < (dx[i]- xmean)/xsd)) {
					data$label[i]	<- "P";
				}
				if (((dy[i]-ymean)/ysd) < -threshold | (threshold < (dy[i]- ymean)/ysd)) {
					data$label[i]	<- "P";
				}
			}
		}
	}

	return (data);
}
