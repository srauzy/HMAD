# -------------------------------------------------------------------------
# Copyright (C) 2018-2023  Stéphane Rauzy
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

## R code for plotting various graphics related to head motion
## author : S. Rauzy, LPL
## date : 27/02/2017

## 
## Return a plot X-Y of the head model
##
## projectname : The project name
##
## return the plot object
lpl.R.dev.hmplot.plotOpenFaceHeadModel <- function(projectname) {

	return (lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlag(projectname, 0, 0, 0, FALSE));
}

##
## Plot and return the head model in the X-Y plane with a given angle
##
## projectname : The project name
## pitchangle: The pitch angle of rotation
## rollangle: The roll angle of rotation
## yawangle: The yaw angle of rotation
## angleflag: put to true to write the angle on the plot
##
## return the plot object
##
lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlag <- function(projectname, pitch, roll, yaw, angleflag) {

	FOLDER_PROJECT <- paste("projects", projectname, sep="/");
	directory <- paste(FOLDER_PROJECT, "OPENFACE/model", sep="/");
	filename <- "dhmt.txt";

	return (lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlagDF(directory, filename, pitch, roll, yaw, angleflag)); 
}

##
## Plot and return the head model in the X-Y plane with a given angle
##
## directory : The directory of the head model file
## filename : The file of the head model 
## pitchangle: The pitch angle of rotation
## rollangle: The roll angle of rotation
## yawangle: The yaw angle of rotation
## angleflag: put to true to write the angle on the plot
##
## return the plot object
##
lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlagDF <- function(directory, filename, pitch, roll, yaw, angleflag) {

	dhmt <- loadInternalDataFrame(directory, filename);
	index <- numeric(nrow(dhmt));
	if (nrow(dhmt) == 67) {
		jc <- 1;
		 for (j in c(1:(length(index)+1))) {
			if (j != 28) {
				index[jc] <- j;
				jc <- jc + 1;
			}
		}
		df <- data.frame(dhmt, index);

		dff <- data.frame("P28", 0, 0, 0, 0, 0, 0, 28);
		names(dff) <- c("point", "X", "seX", "Y", "seY", "Z", "seZ", "index");
		rownames(dff) <- NULL;
		df <- rbind(df, dff);
	} else {
		index <- numeric(nrow(dhmt));
		for (j in c(1:nrow(dhmt))) {
			index[j] <- j;
		}
		df <- data.frame(dhmt, index);
	}
	
	dfr <- rotate(df, pitch, roll, yaw);

	xmin <- min(dfr$X);
	xmax <- max(dfr$X);	
	ymin <- min(dfr$Y);
	ymax <- max(dfr$Y);
	xlength <- xmax-xmin;
	ylength <- ymax-ymin;

	qp <- qplot(data=dfr, x=X, y=Y, label=index) + geom_point(size=I(5), shape=I(19), colour="blue");
	qp <- qp + geom_text(size=I(3), vjust = I(0), nudge_y = I(-1.2), colour="white");
	qp <- qp + scale_x_continuous(limits = c(xmin - xlength/10, xmax + xlength/10));
	yminshift <- ylength/10;
	if (angleflag) {
		textangle <- paste("pitch=", pitch, "°, roll=", roll, "°, yaw=", yaw, "°", sep="");
        	qp <- qp +  annotate("text", label = textangle, x = xmin - xlength/10, y = ymin - 1.8*ylength/10, hjust = 0);
		yminshift <- 2*ylength/10;
	}	
       	qp <- qp + scale_y_continuous(limits = c(ymin - yminshift, ymax + ylength/10)) + coord_fixed(ratio = 1, expand = TRUE);
		##qp <- qp + ggtitle("Head model");	
		##qp <- qp + ggtitle(paste(projectname, "head model"));
		qp <- qp + ylab("Y (in millimeters)");
		qp <- qp + xlab("X (in millimeters)");
		qp <- qp + theme(text = element_text(size=8));

	return (qp);
}

##
## Rotate the data 
##
## data : The data with the genuine 3D coordinates X, Y and Z
## pitchangle: The pitch angle of rotation
## rollangle: The roll angle of rotation
## yawangle: The yaw angle of rotation
##
## return a data frame with the new 3D coordinates
##
rotate <- function(data, pitchangle, rollangle, yawangle) {

	pitch <- numeric(1);
	roll <- numeric(1);
	yaw <- numeric(1);

	pitch[1] <- pitchangle;
	roll[1] <- rollangle;
	yaw[1] <-yawangle;

	dangle <- data.frame(pitch, roll, yaw);
	drc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(dangle, "direct", "XYZ");

	X <- numeric(nrow(data));
	Y <- numeric(nrow(data));
	index <- numeric(nrow(data));

	X <- drc$cxx*data$X + drc$cxy*data$Y + drc$cxz*data$Z;
	Y <- drc$cyx*data$X + drc$cyy*data$Y + drc$cyz*data$Z;
	index  <- data$index;

	dd <- data.frame(index, X, Y);

	return (dd);
}


plotSupportQuality <- function(projectname, ebmtype) {

	FOLDER_PROJECT <- paste("project", projectname, sep="/");
	FOLDER_TABLES <- paste(FOLDER_PROJECT, "tables", sep="/");

	if (ebmtype == "Raise") {
		w2Ds <- loadInternalDataFrame(FOLDER_TABLES, "w2Dsr.txt");
	} else {
		w2Ds <- loadInternalDataFrame(FOLDER_TABLES, "w2Dsf.txt");
	}

	cat("Compute standard deviations (in the wavelet space) in function of the quality...\n");  
	sdwcq <- createWaveletCoefficientStandardDeviationByQualityAndScale(w2Ds);

	return(ggplot(aes(x = quality, y = sdwc, color = duration, group = duration), data = sdwcq) +  geom_point() + geom_line());
}

plotHistogramPredictionObservation <- function(projectname, ebmtype) {

	FOLDER_PROJECT <- paste("project", projectname, sep="/");
	FOLDER_EVAL <- paste(FOLDER_PROJECT, "eval", sep="/");

	if (ebmtype == "Raise") {
		eval <- loadInternalDataFrame(FOLDER_EVAL, "evalr.txt");
	} else {
		eval <- loadInternalDataFrame(FOLDER_EVAL, "evalf.txt");
	}

	e <- subset(eval, eval$igs != -1 & eval$class != "X" & eval$class != "N");
	e$dtmin <- e$tmin - e$gstmin;
	e$dtmax <- e$tmax - e$gstmax;
	e$dduration <- e$duration - e$gsduration;
	e$dt <- (e$dtmin +  e$dtmax)/2;

	dtminhist <- qplot() + geom_histogram(aes(e$dtmin), fill = "black", binwidth=0.1, origin=0)+ xlab("tmin") 

	cat("Compute standard deviations (in the wavelet space) in function of the quality...\n");  
	sdwcq <- createWaveletCoefficientStandardDeviationByQualityAndScale(w2Ds);

	return(ggplot(aes(x = quality, y = sdwc, color = duration, group = duration), data = sdwcq) +  geom_point() + geom_line());
}

##
## Plot and return the head model in the X-Y plane with a given angle
##
## projectname : The project name
## pitchangle: The pitch angle of rotation
## rollangle: The roll angle of rotation
## yawangle: The yaw angle of rotation
## angleflag: put to true to write the angle on the plot
##
## return the plot object
##
lpl.R.dev.hmplot.plotLandmarkProjectedResidualDispersion <- function(projectname, axis) {

	FOLDER_PROJECT <- paste("projects", projectname, sep="/");
	directory <- paste(FOLDER_PROJECT, "OPENFACE/model", sep="/");
	filename <- "rrsdt.txt";

	rrsdt <- loadInternalDataFrame(directory, filename);
	if (axis=="y") {
		values <- rrsdt$sdy;
	} else {
		values <- rrsdt$sdx;
	}

	ymax <- max(values);
	qp <- ggplot(aes(x = landmark, y = values), data = rrsdt) +  geom_point()
	qp <- qp + scale_y_continuous(limits = c(0, ymax + ymax/10));
	qp <- qp + ylab(paste("Residuals dispersion on" , axis, "(in millimeters)"));
	qp <- qp + theme(text = element_text(size=10));

	return (qp);
}

plotStandardDispersionResidual <- function(projectname, sdt, axis) {

	sdm <- melt(rrsdt, c=("landmark"));
	ymax <- max(sdm$value);

	qp <- ggplot(aes(x = point, y = value, color = variable, group = variable), data = sdm) +  geom_point() + geom_line();
	qp <- qp + scale_y_continuous(limits = c(0, ymax + ymax/10));
	qp <- qp + ylab(paste("Residuals dispersion on" , axis, "(in millimeters)"));
	qp <- qp + ggtitle(paste(projectname, "head model"));
	qp <- qp + theme(text = element_text(size=10));

	return (qp);
}

plotErrorAmplitudeFunction <- function(axis) {

	roll <- numeric(101*101);
	pitch <- numeric(101*101);
	yaw <- numeric(101*101);
	S <- numeric(101*101);

	for (i in c(1:101)) {
		for (j in c(1:101)) {
			pitch[(i-1)*101+j] <- i - 51;
			yaw[(i-1)*101+j] <- j - 51;
			S[(i-1)*101+j] <- 1;
		}
	}

	d <- data.frame(pitch, yaw, roll, S);
	irc <- createRotationCoefficientsTable(d, "inverse", "XYZ");
	d <- addInverseProjectedErrorAmplitudes(d, irc);

	if (axis == "x") {
		p <- ggplot(data = d, aes(x = yaw, y = pitch, fill = eaX, z =
eaX)) + geom_tile() + stat_contour(bins=10) + scale_fill_gradientn(limits =
c(1,max(d$eaX)), colours=c("navyblue", "darkmagenta", "darkorange1",
"yellow"), na.value="navyblue") + ggtitle("Error amplitude on the X axis");
	} else {
		p <- ggplot(data = d, aes(x = yaw, y = pitch, fill = eaY, z =
eaY)) + geom_tile() + stat_contour(bins=10) + scale_fill_gradientn(limits =
c(1,max(d$eaY)), colours=c("navyblue", "darkmagenta", "darkorange1",
"yellow"), na.value="navyblue") + ggtitle("Error amplitude on the Y axis");
	}

	return (p);
}
