## R code for plotting various graphics related to head motion
## author : S. Rauzy, LPL
## date : 27/02/2017

lpl.R.dev.hmplot.plotOpenFaceHeadModel <- function(projectname) {

	return (lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlag(projectname, 0, 0, 0, FALSE));
}


plotHeadModel <- function(projectname) {

	return (plotHeadModelWithAngleAndFlag(projectname, 0, 0, 0, FALSE));
}

plotHeadModelWithAngle <- function(projectname, pitch, roll, yaw) {

	return (plotHeadModelWithAngleAndFlag(projectname, pitch, roll, yaw, TRUE));
}

plotHeadModelWithAngleAndFlag <- function(projectname, pitch, roll, yaw, angleflag) {

	FOLDER_PROJECT <- paste("projects", projectname, sep="/");
	FOLDER_MODEL <- paste(FOLDER_PROJECT, "model", sep="/");

	dhmt <- loadInternalDataFrame(FOLDER_MODEL, "dhmt.txt");
	index <- numeric(nrow(dhmt));
	jc <- 1;
	 for (j in c(1:(length(index)+1))) {
		if (j != 5) {
			index[jc] <- j;
			jc <- jc + 1;
		}
	}
	df <- data.frame(dhmt, index);

	dff <- data.frame("P05", 0, 0, 0, 0, 0, 0, 5);
	names(dff) <- c("point", "X", "seX", "Y", "seY", "Z", "seZ", "index");
	rownames(dff) <- NULL;
	df <- rbind(df, dff);
	
	dfr <- rotate(df, pitch, roll, yaw);

	xmin <- min(dfr$X);
	xmax <- max(dfr$X);	
	ymin <- min(dfr$Y);
	ymax <- max(dfr$Y);
	xlength <- xmax-xmin;
	ylength <- ymax-ymin;

	qp <- qplot(data=dfr, x=X, y=Y, label=index) + geom_point(size=I(7), shape=I(19), colour="blue");
	qp <- qp + geom_text(size=I(4), vjust = I(0), nudge_y = I(-1.2), colour="white");
	qp <- qp + scale_x_continuous(limits = c(xmin - xlength/10, xmax + xlength/10));
	yminshift <- ylength/10;
	if (angleflag) {
		textangle <- paste("pitch=", pitch, "°, roll=", roll, "°, yaw=", yaw, "°", sep="");
        	qp <- qp +  annotate("text", label = textangle, x = xmin - xlength/10, y = ymin - 1.8*ylength/10, hjust = 0);
		yminshift <- 2*ylength/10;
	}	
       	qp <- qp + scale_y_continuous(limits = c(ymin - yminshift, ymax + ylength/10)) + coord_fixed(ratio = 1, expand = TRUE);
        if (is.null(projectname)) {
		
	} else {
		##qp <- qp + ggtitle("Head model");	
		##qp <- qp + ggtitle(paste(projectname, "head model"));
		qp <- qp + ylab("Y (pixels)");
		qp <- qp + xlab("X (pixels)");
		qp <- qp + theme(text = element_text(size=8));
	}

	return (qp);

}

lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlag <- function(projectname, pitch, roll, yaw, angleflag) {

	FOLDER_PROJECT <- paste("projects", projectname, sep="/");
	directory <- paste(FOLDER_PROJECT, "OPENFACE/model", sep="/");
	filename <- "dhmt.txt";

	return (lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlag1(directory, filename, pitch, roll, yaw, angleflag)); 
}

plotOpenFaceHeadModelDifferenceWithAngleAndFlag1 <- function(directory, filename, pitch, roll, yaw, angleflag) {

	dhmt <- loadInternalDataFrame(directory, filename);
	index <- numeric(nrow(dhmt));
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
	
	##dfr <- rotate(df, pitch, roll, yaw);
	dfr <- df;

	dfr$Y <- dfr$Y - df$Y;
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
        if (is.null(projectname)) {
		
	} else {
		##qp <- qp + ggtitle("Head model");	
		##qp <- qp + ggtitle(paste(projectname, "head model"));
		qp <- qp + ylab("Y (pixels)");
		qp <- qp + xlab("X (pixels)");
		qp <- qp + theme(text = element_text(size=8));
	}

	return (qp);

}

lpl.R.dev.hmplot.plotOpenFaceHeadModelWithAngleAndFlag1 <- function(directory, filename, pitch, roll, yaw, angleflag) {

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
	
	##dfr <- rotate(df, pitch, roll, yaw);
	dfr <- df;

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
		qp <- qp + ylab("Y (pixels)");
		qp <- qp + xlab("X (pixels)");
		qp <- qp + theme(text = element_text(size=8));

	return (qp);

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



rotate <- function(data, pitchangle, rollangle, yawangle) {

	pitch <- numeric(1);
	roll <- numeric(1);
	yaw <- numeric(1);

	pitch[1] <- pitchangle;
	roll[1] <- rollangle;
	yaw[1] <-yawangle;

	dangle <- data.frame(pitch, roll, yaw);
	drc <- createRotationCoefficientsTable(dangle, "direct", "XYZ");

	X <- numeric(nrow(data));
	Y <- numeric(nrow(data));
	index <- numeric(nrow(data));

	X <- drc$cxx*data$X + drc$cxy*data$Y + drc$cxz*data$Z;
	Y <- drc$cyx*data$X + drc$cyy*data$Y + drc$cyz*data$Z;
	index  <- data$index;

	dd <- data.frame(index, X, Y);

	return (dd);
}

plotStandardDispersionResidual <- function(projectname, sdt, axis) {

	sdm <- melt(sdt, c=("point"));
	ymax <- max(sdm$value);

	qp <- ggplot(aes(x = point, y = value, color = variable, group = variable), data = sdm) +  geom_point() + geom_line();
	qp <- qp + scale_y_continuous(limits = c(0, ymax + ymax/10));
	qp <- qp + ylab(paste("Residuals dispersion on" , axis, "(pixels)"));
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
