## R code for computing head model and residuals motions from Intrace software output 
## author : S. Rauzy, LPL
## date : 15/12/2016

##
## Modify the table to present the ouput as the point center_index, and the remaining
## points as vectors relative to point center_index

## center_index : The index of the center landmark
## d : The data frame containing the OpenFace output transformed in data frame
##
## return the modified data frame
##
lpl.R.dev.faceOutputAnalysis.centerOn <- function(center_index, d) {

	## Copy the first 8 columns
	dd <- data.frame(d[, 1:8]);

	dxc <- as.numeric(as.character(d[, (8+2*center_index-1)]));
	dyc <- as.numeric(as.character(d[, (8+2*center_index)]));

	## Add the X and Y of point center_index
	dd <- data.frame(dd, dxc, dyc);
	colnames(dd)[ncol(dd)-1] <-  paste(paste("M", lpl.R.dev.faceOutputAnalysis.getStringNumber(center_index), sep=""), "x", sep="");
	colnames(dd)[ncol(dd)] <-  paste(paste("M",  lpl.R.dev.faceOutputAnalysis.getStringNumber(center_index), sep=""), "y", sep="");

	number_of_points <- length(grep("fl.*", colnames(d), perl=TRUE, value=FALSE))/2;

	for (j in c(1:number_of_points)) {
		difx <- numeric(nrow(d));
		dify <- numeric(nrow(d));
		if (j != center_index) {
			dx <- as.numeric(as.character(d[, (8+2*j-1)]));
			dy <- as.numeric(as.character(d[, (8+2*j)]));
			
			difx <- round(dx - dxc, 4);
			dify <- round(dy - dyc, 4);
			dd <- data.frame(dd, difx, dify);
			colnames(dd)[ncol(dd)-1] <-  paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
			colnames(dd)[ncol(dd)] <-  paste(paste("vM", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
		}
	}
	return (dd);
}

##
## Invert the y-axis (all the y become negative)
##
## d : The data frame containing centered data
##
## return the modified data frame
##
lpl.R.dev.faceOutputAnalysis.invertYAxis <- function(d) {

	iy <- grep("M..y", colnames(d), perl=TRUE, value=FALSE)[1];
	d[, iy] <- -d[, iy];

	li <- grep("vM..y", colnames(d), perl=TRUE, value=FALSE);
	number_of_points <- length(li);
	for (j in c(1:number_of_points)) {
		iy <- li[j];
		d[, iy] <- -d[, iy];
	}
	return (d);
}

##
## Convert a number by adding 0 at its left when it is lower than 10
##
## number : The number
##
## return the Sting resul of the operation
##
lpl.R.dev.faceOutputAnalysis.getStringNumber <- function(number) {

	if (number < 10) {
		return (paste("0", number, sep=""));
	} else {
		return (as.character(number));
	}
}

##
## Create the matrix of rotation coefficients for a 3 given rotation angles
##
## pitch : the pitch angle
## yaw : the yaw angle
## roll : the roll angle
## direction : "direct" for the direct rotation, "inverse" for the inverse rotation (XYZ convention)
##
## return the 3*3 rotation matrix
##
lpl.R.dev.faceOutputAnalysis.createRotationCoefficients <- function(pitch, yaw, roll, direction) {
	
	## Degree to radian conversion coefficient
	alpha = 2*pi/360;

	if (direction == "inverse") {
		pitch = -pitch;
		roll = -roll;
		yaw = -yaw;
	}

	sx = sin(alpha*pitch);
	cx = cos(alpha*pitch);
	sz = sin(alpha*roll);
	cz = cos(alpha*roll);
	sy = sin(alpha*yaw);
	cy = cos(alpha*yaw);

	rc <- matrix(nrow=3, ncol=3);
	## Compute the sinus and cosinus of rotation angles, if direction is "inverse", 
	## we compute this quantities for the inverse rotation
	if (direction == "direct") {
		## Compute the 9 rotation coefficients for the Rx*Ry*Rz rotation
		rc[1, 1] <- cz*cy;
		rc[1, 2] <- sz*cx + cz*sy*sx;
		rc[1, 3] <- sz*sx - cz*sy*cx;
		rc[2, 1] <- -sz*cy;
		rc[2, 2] <- cz*cx - sz*sy*sx;
		rc[2, 3] <- cz*sx + sz*sy*cx;
		rc[3, 1] <- sy;
		rc[3, 2] <- -cy*sx;
		rc[3, 3] <- cy*cx;
	} else if (direction == "inverse") {
		## Compute the 9 rotation coefficients for the Rz*Ry*Rx rotation
		rc[1, 1] <- cz*cy;
		rc[1, 2] <- sz*cy;
		rc[1, 3] <- -sy;
		rc[2, 1] <- -sz*cx + cz*sy*sx;
		rc[2, 2] <- cz*cx + sz*sy*sx;
		rc[2, 3] <- cy*sx;
		rc[3, 1] <- sz*sx + cz*sy*cx;
		rc[3, 2] <- -cz*sx + sz*sy*cx;
		rc[3, 3] <- cy*cx;
	}
	return (rc);
}

##
## Create the rotation coefficients
##
## d : The data frame containing pitch, yaw and roll angles
## direction : "direct" for the direct rotation, "inverse" for the inverse rotation
## type : The six possibilities for the rotation (here, we use "XYZ" convention)
##
## return a data frame containing the 9 columns of rotation coefficients
##
lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable <- function(d, direction, type) {
	
	## Degree to radian conversion coefficient
	alpha = 2*pi/360;

	## Pitch is the rotation around x-axis, counterclockwise
	## Roll is the rotation around z-axis, counterclockwise, it is inverted because 
	## we keep the y-axis pointing downward
	## Yaw is the rotation around the y-axis, counterclockwise, it is inverted because 
	## we keep the y-axis pointing downward
	## Compute the sinus and cosinus of rotation angles, if direction is "inverse", 
	## we compute this quantities for the inverse rotation
	if (direction == "direct") {
		pitch = d$pitch;
		roll = d$roll;
		yaw = d$yaw;
	} else if (direction == "inverse") {
		pitch = -d$pitch;
		roll = -d$roll;
		yaw = -d$yaw;
	}

	sx = sin(alpha*pitch);
	cx = cos(alpha*pitch);
	sz = sin(alpha*roll);
	cz = cos(alpha*roll);
	sy = sin(alpha*yaw);
	cy = cos(alpha*yaw);

	if (type == "XYZ") {
		## Compute the sinus and cosinus of rotation angles, if direction is "inverse", 
		## we compute this quantities for the inverse rotation
		if (direction == "direct") {
			## Compute the 9 rotation coefficients for the Rx*Ry*Rz rotation
			cxx <- cz*cy;
			cxy <- sz*cx + cz*sy*sx;
			cxz <- sz*sx - cz*sy*cx;
			cyx <- -sz*cy;
			cyy <- cz*cx - sz*sy*sx;
			cyz <- cz*sx + sz*sy*cx;
			czx <- sy;
			czy <- -cy*sx;
			czz <- cy*cx;
		} else if (direction == "inverse") {
			## Compute the 9 rotation coefficients for the Rz*Ry*Rx rotation
			cxx <- cz*cy;
			cxy <- sz*cy;
			cxz <- -sy;
			cyx <- -sz*cx + cz*sy*sx;
			cyy <- cz*cx + sz*sy*sx;
			cyz <- cy*sx;
			czx <- sz*sx + cz*sy*cx;
			czy <- -cz*sx + sz*sy*cx;
			czz <- cy*cx;
		}
	} else if (type == "ZYX") {
		if (direction == "direct") {
			## Compute the 9 rotation coefficients for the Rz*Ry*Rx rotation
			cxx <- cz*cy;
			cxy <- sz*cy;
			cxz <- -sy;
			cyx <- -sz*cx + cz*sy*sx;
			cyy <- cz*cx + sz*sy*sx;
			cyz <- cy*sx;
			czx <- sz*sx + cz*sy*cx;
			czy <- -cz*sx + sz*sy*cx;
			czz <- cy*cx;	
		} else if (direction == "inverse") {
			## Compute the 9 rotation coefficients for the Rx*Ry*Rz rotation
			cxx <- cz*cy;
			cxy <- sz*cx + cz*sy*sx;
			cxz <- sz*sx - cz*sy*cx;
			cyx <- -sz*cy;
			cyy <- cz*cx - sz*sy*sx;
			cyz <- cz*sx + sz*sy*cx;
			czx <- sy;
			czy <- -cy*sx;
			czz <- cy*cx;
		}
	} else if (type == "XZY") {
		if (direction == "direct") {
			## Compute the 9 rotation coefficients for the Rx*Rz*Ry rotation
			cxx <- cz*cy;
			cxy <- cx*sz*cy + sx*sy;
			cxz <- sx*sz*cy - cx*sy;
			cyx <- -sz;
			cyy <- cx*cz;
			cyz <- sx*cz;
			czx <- cz*sy;
			czy <- cx*sz*sy - sx*cy;
			czz <- sx*sz*sy + cx*cy;	
		} else if (direction == "inverse") {
			## Compute the 9 rotation coefficients for the Ry*Rz*Rx rotation
			cxx <- cy*cz;
			cxy <- sz;
			cxz <- -sy*cz;
			cyx <- -cy*sz*cx + sy*sx;
			cyy <- cz*cx;
			cyz <- sy*sz*cx + cy*sx;
			czx <- cy*sz*sx + sy*cx;
			czy <- -cz*sx;
			czz <- -sy*sz*sx + cy*cx;
		}
	} else if (type == "YZX") {
		if (direction == "direct") {
			## Compute the 9 rotation coefficients for the Ry*Rz*Rx rotation
			cxx <- cy*cz;
			cxy <- sz;
			cxz <- -sy*cz;
			cyx <- -cy*sz*cx + sy*sx;
			cyy <- cz*cx;
			cyz <- sy*sz*cx + cy*sx;
			czx <- cy*sz*sx + sy*cx;
			czy <- -cz*sx;
			czz <- -sy*sz*sx + cy*cx;	
		} else if (direction == "inverse") {
			## Compute the 9 rotation coefficients for the Rx*Rz*Ry rotation
			cxx <- cz*cy;
			cxy <- cx*sz*cy + sx*sy;
			cxz <- sx*sz*cy - cx*sy;
			cyx <- -sz;
			cyy <- cx*cz;
			cyz <- sx*cz;
			czx <- cz*sy;
			czy <- cx*sz*sy - sx*cy;
			czz <- sx*sz*sy + cx*cy;
		}
	} else if (type == "YXZ") {
		if (direction == "direct") {
			## Compute the 9 rotation coefficients for the Ry*Rx*Rz rotation
			cxx <- cy*cz + sy*sx*sz;
			cxy <- cx*sz;
			cxz <- -sy*cx+cy*sx*sz;
			cyx <- -cy*sz + sy*sx*cz;
			cyy <- cx*cz;
			cyz <- sy*sz + cy*sx*cz;
			czx <- sy*cx;
			czy <- -sx;
			czz <- cy*cx;	
		} else if (direction == "inverse") {
			## Compute the 9 rotation coefficients for the Rz*Rx*Ry rotation
			cxx <- cz*cy - sz*sx*sy;
			cxy <- sz*cy + cz*sx*sy;
			cxz <- -cx*sy;
			cyx <- -sz*cx;
			cyy <- cz*cx;
			cyz <- sx;
			czx <- cz*sy + sz*sx*cy;
			czy <- sz*sy - cx*sx*cy;
			czz <- cx*cy;
		}
	} else if (type == "ZXY") {
		if (direction == "direct") {
			## Compute the 9 rotation coefficients for the Rz*Rx*Ry rotation
			cxx <- cz*cy - sz*sx*sy;
			cxy <- sz*cy + cz*sx*sy;
			cxz <- -cx*sy;
			cyx <- -sz*cx;
			cyy <- cz*cx;
			cyz <- sx;
			czx <- cz*sy + sz*sx*cy;
			czy <- sz*sy - cx*sx*cy;
			czz <- cx*cy;	
		} else if (direction == "inverse") {
			## Compute the 9 rotation coefficients for the Ry*Rx*Rz rotation
			cxx <- cy*cz + sy*sx*sz;
			cxy <- cx*sz;
			cxz <- -sy*cx+cy*sx*sz;
			cyx <- -cy*sz + sy*sx*cz;
			cyy <- cx*cz;
			cyz <- sy*sz + cy*sx*cz;
			czx <- sy*cx;
			czy <- -sx;
			czz <- cy*cx;	
		}
	}

	return (data.frame(cxx, cxy, cxz, cyx, cyy, cyz, czx, czy, czz));
}

##
## Create the head model and the S factor from the data with clean label
##
## iodfc : The Intraface output (transformed in data frame) with clean annotation
##
lpl.R.dev.faceOutputAnalysis.createHeadAndSFactorModel <- function(iodfc) {

	d <- subset(iodfc, iodfc$label == "Y");
	hasfm <- lpl.R.dev.faceOutputAnalysis.iterativeHeadAndFactorSModel(d, 0.05);

	return (hasfm);
}


##
## Load a data frame written in file filename in directory directory
##
## directory : The directory
## filename : The fine name
##
## return the data frame
##
loadInternalDataFrame <- function(directory, filename) {

	df <- read.table(paste(directory, filename, sep="/"), h=TRUE);

	return (df);
}

loadInternalDataFrameWithoutFactor <- function(directory, filename) {

	df <- read.table(paste(directory, filename, sep="/"), h=TRUE, stringsAsFactors=FALSE);

	return (df);
}

##
## Save a data frame in file filename in directory directory
##
## data : The data frame
## directory : The directory
## filename : The fine name
##
saveInternalDataFrame <- function(data, directory, filename) {

	if (!file.exists(directory)) {
		dir.create(directory);
	}
	write.table(data,  paste(directory, filename, sep="/"), sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE);	
}

searchForBinIndex <- function(value, binfrontier) {

	nbrbin <- length(binfrontier)+1;
	for (k in c(1:nbrbin)) {
		if (k == 1) {
			if (value <= binfrontier[k]) {
				return (k);
			}
		} else if (k == nbrbin) {
			if (value > binfrontier[k-1]) {
				return (k);
			}
		} else {
			if (value <= binfrontier[k] & value > binfrontier[k-1]) {
				return (k);
			}
		}
	}
}

##
## Add the projected residuals computed from the head model
##
## data : The data frame containing the S factor value and the angles
## irc : The inverse rotation coefficients 
## dhmt : The head model table
## rrsdt : The table containing the relative residual standard dispersion
##
## return the data frame with residuals added
##
lpl.R.dev.faceOutputAnalysis.addProjectedResiduals <- function(data, irc, dhmt, rrsdt) {

	d <- lpl.R.dev.faceOutputAnalysis.createProjectedResiduals(data, irc, dhmt, rrsdt);
	return (data.frame(data, d));
}

##
## Create the projected residuals rPi - rP05
##
## data : The data frame containing the S factor value and the angles
## irc : The inverse rotation coefficients 
## dhmt : The dirext head model table
##
## return data frame containing the relative projected residuals
##
lpl.R.dev.faceOutputAnalysis.createRelativeProjectedResiduals <- function(data, irc, dhmt) {

	lix <- grep("vM..x", colnames(data), perl=TRUE, value=FALSE);
	liy <- grep("vM..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);
	
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
		pointindex <- substring(colnames(data)[lix[j]], 3,4);
		listpointX[[j]] <- numeric(nrow(data));
		listpointY[[j]] <- numeric(nrow(data));
		column_index_x = lix[j];
		column_index_y = liy[j];

		xobs <- data[, column_index_x]/data$S;
		yobs <- data[, column_index_y]/data$S;
	
		Xobs = (cxx - cxz*czx/czz)*xobs + (cxy - cxz*czy/czz)*yobs;
		Yobs = (cyx - cyz*czx/czz)*xobs + (cyy - cyz*czy/czz)*yobs;

		listpointX[[j]] = round(Xobs - dhmt$X[j] + cxz/czz*dhmt$Z[j], 4);
		listpointY[[j]] = round(Yobs - dhmt$Y[j] + cyz/czz*dhmt$Z[j], 4);

		if (j == 1) {
			dd <- data.frame(listpointX[[j]], listpointY[[j]]);
		} else {
			dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		}
		colnames(dd)[ncol(dd)-1] <- paste(paste("rP", pointindex, sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-paste(paste("rP", pointindex, sep=""), "y", sep="");
	}
	
	return (dd);
}

#
## Create the projected residuals
##
## data : The data frame containing the S factor value and the angles
## irc : The inverse rotation coefficients 
## dhmt : The dirext head model table
## rrsdt : The table containing the relative residual standard dispersion
##
## return the data frame containing the residuals
##
lpl.R.dev.faceOutputAnalysis.createProjectedResiduals <- function(data, irc, dhmt, rrsdt) {

	lix <- grep("vM..x", colnames(data), perl=TRUE, value=FALSE);
	liy <- grep("vM..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	r <- lpl.R.dev.faceOutputAnalysis.createRelativeProjectedResiduals(data, irc, dhmt);
	rP05x <- numeric(nrow(r));
	rP05y <- numeric(nrow(r));

	## The center residuals are minus the weighted bulk flow of the other landmarks
	for (j in c(1:number_of_points)) {
		column_index_x = 2*(j-1)+1;
		column_index_y = 2*(j-1)+2;
		rP05x <- rP05x + r[, column_index_x]/rrsdt$sdx[j];
		rP05y <- rP05y + r[, column_index_y]/rrsdt$sdy[j];
	}
	rP05x <- -rP05x/sum(1/rrsdt$sdx);
	rP05y <- -rP05y/sum(1/rrsdt$sdy);

	listpointX <- vector(mode="list", length=number_of_points+1);
	listpointY <- vector(mode="list", length=number_of_points+1);

	for (j in c(1:length(listpointX))) {
		column_index_x = which(colnames(r) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
		column_index_y = which(colnames(r) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));
		if (length(column_index_x) != 0) {
			listpointX[[j]] <- round(r[, column_index_x] + rP05x, 4);
			listpointY[[j]] <- round(r[, column_index_y] + rP05y, 4);
		
		} else {
			listpointX[[j]] <- round(rP05x, 4);
			listpointY[[j]] <- round(rP05y, 4);
		}
		if (j == 1) {
			dd <- data.frame(listpointX[[j]], listpointY[[j]]);
		} else {
			dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		}
			
		colnames(dd)[ncol(dd)-1] <- paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep="");
	}
	
	return (dd);
}

##
## Add the projected error amplitude an X and Y which depends on the rotation angles
##
## data : The data frame containing the angles
## irc : The inverse rotation coefficients 
##
## return the data frame with the two columns added
##
lpl.R.dev.faceOutputAnalysis.addProjectedErrorAmplitudes <- function(data, irc) {

	return (data.frame(data, lpl.R.dev.faceOutputAnalysis.createProjectedErrorAmplitudes(data, irc)));
}

##
## Create the projected error amplitude on X and Y which depends on the rotation angles
##
## data : The data frame containing the angles
## irc : The inverse rotation coefficients 
##
## return a data frame with the two columns errors amplitudes
##
lpl.R.dev.faceOutputAnalysis.createProjectedErrorAmplitudes <- function(data, irc) {

	cxx <- irc$cxx;
	cxy <- irc$cxy;
	cxz <- irc$cxz;

	cyx <- irc$cyx;
	cyy <- irc$cyy;
	cyz <- irc$cyz;

	czx <- irc$czx;
	czy <- irc$czy;
	czz <- irc$czz;

	eaX <- numeric(nrow(data));
	eaY <- numeric(nrow(data));
	
	vXobs = (cxx - cxz*czx/czz)^2 + (cxy - cxz*czy/czz)^2;
	vYobs = (cyx - cyz*czx/czz)^2 + (cyy - cyz*czy/czz)^2;

	eaX = round(sqrt(vXobs/(data$S^2) + (cxz/czz)^2), 4);
	eaY = round(sqrt(vYobs/(data$S^2) + (cyz/czz)^2), 4);

	return (data.frame(eaX, eaY));
}

##
## Compute the S factor for each frame corrected from the head rotation effect
##
## data : The data frame containing the relative coordinates of the landmarks
## drc : The direct rotation coefficients for each frame
## dhmt : The data frame containing the (direct) head model parameters
## normalize : The normalization flag (TRUE for a S factor equals 1 in average) 
##
## return the corrected S factor
##
lpl.R.dev.faceOutputAnalysis.computeRotationCorrectedSFactor <- function(data, drc, dhmt, normalize) {

	S <- numeric(nrow(data));
	w <- numeric(nrow(data));

	cxx <- drc$cxx;
	cxy <- drc$cxy;
	cxz <- drc$cxz;

	cyx <- drc$cyx;
	cyy <- drc$cyy;
	cyz <- drc$cyz;

	lix <- grep("vM..x", colnames(data), perl=TRUE, value=FALSE);
	liy <- grep("vM..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	for (j in c(1:number_of_points)) {
		
		column_index_x = lix[j];
		column_index_y = liy[j];

		xden <-abs(cxx*dhmt$X[j] + cxy*dhmt$Y[j] + cxz*dhmt$Z[j]);
		yden <- abs(cyx*dhmt$X[j] + cyy*dhmt$Y[j] + cyz*dhmt$Z[j]);
		
		w <- w + xden + yden;

		xobs <- abs(data[, column_index_x]);
		yobs <- abs(data[, column_index_y]);	

		S <- S + xobs + yobs;
	}

	S <- round(S/w, 4);
	
	if (normalize) {
		mS <- mean(S, na.rm = T);
		S <- round(S/mS, 4);
	}
	return (S);
}

lpl.R.dev.faceOutputAnalysis.computeCenterResidualData <- function(data, irc, ihmt) {

	lix <- grep("vM..x", colnames(data), perl=TRUE, value=FALSE);
	liy <- grep("vM..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	listpointX <- vector(mode="list", length=number_of_points);
	listpointY <- vector(mode="list", length=number_of_points);

	dd <- data;

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

		column_index_x = lix[j];
		column_index_y = liy[j];

		xobs <- data[, column_index_x]/data$S;
		yobs <- data[, column_index_y]/data$S;
	
		Xobs = (cxx - cxz*czx/czz)*xobs + (cxy - cxz*czy/czz)*yobs;
		Yobs = (cyx - cyz*czx/czz)*xobs + (cyy - cyz*czy/czz)*yobs;

		listpointX[[j]] <- round(Xobs-ihmt$X[j] + cxz/czz*ihmt$Z[j], 4);
		listpointY[[j]] <- round(Yobs-ihmt$Y[j] + cyz/czz*ihmt$Z[j], 4);

		if (j == 1) {
			dd <- data.frame(listpointX[[j]], listpointY[[j]]);
		} else {
			dd <- data.frame(dd, listpointX[[j]], listpointY[[j]]);
		}
		pointindex <- substring(colnames(data)[lix[j]], 3,4);
		colnames(dd)[ncol(dd)-1] <- paste(paste("rsP", pointindex, sep=""), "x", sep="");
		colnames(dd)[ncol(dd)] <-paste(paste("rsP", pointindex, sep=""), "y", sep="");
	}
	return (dd);
}

##
## Compute for each point the weight representing the contribution of each point to the calculation of the S factor
##
## fsd : The data frame containing for each point and each time the individual value of the S factor
##
lpl.R.dev.faceOutputAnalysis.computeCenterResidualModel <- function(rp5d, axis) {

	lix <- grep("rsP..x", colnames(rp5d), perl=TRUE, value=FALSE);
	liy <- grep("rsP..y", colnames(rp5d), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	point <- character(number_of_points);
	w <- numeric(number_of_points);
	sw <- 0;

	for (j in c(1:number_of_points)) {
		pointindex <- substring(colnames(rp5d)[lix[j]], 4,5);
		column_index = which(colnames(rp5d) == paste(paste("rsP", pointindex, sep=""), axis, sep=""));
		s <- rp5d[, column_index];
	
		w[j] <- 1/sd(s, na.rm = T);
		sw <- sw + w[j];
	}

	w <- round(w/sw, 4);
	return (data.frame(point, w));
}

##
## Compute the S factor given the position of points relative to P05, the
## inverse rotation coefficients, the head model the individual S factor and the S model
##
## fsd : The indivual values of the S factor for each point over X and Y axis
## fsm : The S model containing the weights for each point for the S factor calculation
##
lpl.R.dev.faceOutputAnalysis.computeCenterResidual <- function(rp5d, rp5m, axis) {

	V <- numeric(nrow(rp5d));
	
	lix <- grep("rsP..x", colnames(rp5d), perl=TRUE, value=FALSE);
	liy <- grep("rsP..y", colnames(rp5d), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	for (j in c(1:number_of_points)) {
		pointindex <- substring(colnames(rp5d)[lix[j]], 4,5);
		column_index = which(colnames(rp5d) == paste(paste("rsP", pointindex, sep=""), axis, sep=""));
		p5r <- rp5d[, column_index];
		p5r <- p5r + rp5m$w[j]*p5r;
	}

	return (p5r);
}

##
## Add to the data frame the S factor column initialized to 1
##
## data : The data frame
##
## return the resulting data frame
##
pl.R.dev.faceOutputAnalysis.addFiducialSFactor <- function(data) {

	data$S  <- 1;
	return (data);
}

##
## Create iteratively the head model and S factor
##
## data : The data frame containing the position of points relative to P05
## precision : The threshold precision under which the iteration is stopped
## (the quantity estimated is the difference between two consecutive head models)
##
## return the data frame with removed outliers and with the S factor column
##
lpl.R.dev.faceOutputAnalysis.iterativeHeadAndFactorSModel <- function(data, precision) {

	data <- pl.R.dev.faceOutputAnalysis.addFiducialSFactor(data);
	drc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(data, "direct", "XYZ");

	lix <- grep("vM..x", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix) + 1;
	## First pass, obtain S1(t)
	cat(paste("Compute the", number_of_points, "landmark positions...\n"));
	dhmt <- lpl.R.dev.faceOutputAnalysis.createDirectHeadModel(data, drc, FALSE);
	irc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(data, "inverse", "XYZ");
	r <- lpl.R.dev.faceOutputAnalysis.createRelativeProjectedResiduals(data, irc, dhmt);	
	sdr <- lpl.R.dev.faceOutputAnalysis.computeMeanProjectedResidualsStandardDispersion(r);
	cat("--------------------\n");
	cat(paste("First standard deviation for residuals :", sdr, "\n"));
	cat("First calculation for the S factor...\n");
	S <- lpl.R.dev.faceOutputAnalysis.computeRotationCorrectedSFactor(data, drc, dhmt, TRUE);
	data <- lpl.R.dev.faceOutputAnalysis.replaceFactorS(data, S);

	rdP5 <- lpl.R.dev.faceOutputAnalysis.computeCenterResidualData(data, irc, dhmt);
	rdP5xm <- lpl.R.dev.faceOutputAnalysis.computeCenterResidualModel(rdP5, "x");
	rdP5ym <- lpl.R.dev.faceOutputAnalysis.computeCenterResidualModel(rdP5, "y");
	rP5x <- lpl.R.dev.faceOutputAnalysis.computeCenterResidual(rdP5, rdP5xm, "x");
	rP5y <- lpl.R.dev.faceOutputAnalysis.computeCenterResidual(rdP5, rdP5ym, "y");

	cis = which(colnames(data) == "S");
	df <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(data, cis, 3);
	irc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(df, "inverse", "XYZ");
	drc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(df, "direct", "XYZ");
	
	count = 0;

	stop = FALSE;

	while (!stop) {
		
		dhmtsafe <- dhmt;
		dhmt <- lpl.R.dev.faceOutputAnalysis.createDirectHeadModel(df, drc, FALSE);
		hmdif <- lpl.R.dev.faceOutputAnalysis.computeHeadModelsDifference(dhmtsafe, dhmt);
		cat("--------------------\n");
		cat(paste("Iteration",(count+1), "in the calculation of the S factor and head model...\n"));
		cat(paste("Head models difference :", hmdif, "\n"));
		r <- lpl.R.dev.faceOutputAnalysis.createRelativeProjectedResiduals(df, irc, dhmt);	
		sdr <- lpl.R.dev.faceOutputAnalysis.computeMeanProjectedResidualsStandardDispersion(r);
		cat(paste("Standard deviation for residuals :", sdr, "\n"));
		S <- lpl.R.dev.faceOutputAnalysis.computeRotationCorrectedSFactor(df, drc, dhmt, TRUE);
		c <- cor(S, df$S);
		count = count+1;
		cat("Correlation with previous S =", c, "\n");
		df <- lpl.R.dev.faceOutputAnalysis.replaceFactorS(df, S);
		if (hmdif < precision | count > 10) {
			stop <- TRUE;
		}
	}
	cat("--------------------\n");

	return (df);
}

##
## Compute the difference between two head models
##
## hmt1 : The head model table 1 
## hmt2 : The head model table 2 
##
## return the value of the difference
##
lpl.R.dev.faceOutputAnalysis.computeHeadModelsDifference <- function(hmt1, hmt2) {

	varX <- 0;
	varY <- 0;
	varZ <- 0;
	for (j in c(1:nrow(hmt1))) {
			varX <- varX + (hmt1$X[j] - hmt2$X[j])^2;
			varY <- varY + (hmt1$Y[j] - hmt2$Y[j])^2;
			varZ <- varZ + (hmt1$Z[j] - hmt2$Z[j])^2;
	}
	return (sqrt((varX+varY+varZ)/3));
}

##
## Add to the data frame the residuals on X and Y and the errors amplitudes
##
## df : The data frame containing the position of the points relative to P05
## dhmt : The head model (direct solution)
## rrsdt : The table containing the relative residual standard dispersion
##
## return the data frame with residuals and amplitude errors columns added
##
lpl.R.dev.faceOutputAnalysis.addResidualsAndErrorsAmplitudes <- function(df, dhmt, rrsdt) {

	cat("Create residuals from the head and S factor model...\n");
	drc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(df, "direct", "XYZ");
	S <- lpl.R.dev.faceOutputAnalysis.computeRotationCorrectedSFactor(df, drc, dhmt, FALSE);
	df <- lpl.R.dev.faceOutputAnalysis.replaceFactorS(df, S);

	irc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficientsTable(df, "inverse", "XYZ");
	df <- lpl.R.dev.faceOutputAnalysis.addProjectedResiduals(df, irc, dhmt, rrsdt);
	df <- lpl.R.dev.faceOutputAnalysis.addProjectedErrorAmplitudes(df, irc);

	return(df);
}

##
## Compute the mean standard dispersion of residuals over the values (X and Y)
##
## data : The data frame containing the values of the residuals (data filtering has to be done before this step)
##
## return a data frame containing for each residual (first column), the residuals on X, Y and potentially Z column if the residuals exits
##
lpl.R.dev.faceOutputAnalysis.computeProjectedResidualsStandardDispersion <- function(data) {

	tix <- grep("rP..x", colnames(data), perl=TRUE, value=FALSE);
	tiy <- grep("rP..y", colnames(data), perl=TRUE, value=FALSE);
	tiz <- grep("rP..z", colnames(data), perl=TRUE, value=FALSE);

	sdx <- numeric(length(tix));
	sdy <- numeric(length(tix));
	landmark <- character(length(tix));
		
	for (j in c(1:length(sdx))) {
		sdx[j] <- round(sd(data[, tix[j]], na.rm = T), 4);
		sdy[j] <- round(sd(data[, tiy[j]], na.rm = T), 4);
		landmark[j] <- paste("P", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep="");
	}

	d <- data.frame(landmark, sdx, sdy);

	if (length(tiz) != 0) {
		sdz <- numeric(length(tiz));
		for (j in c(1:length(sdz))) {
			sdz[j] <- round(sd(data[, tiz[j]], na.rm = T), 4);
		}
		d <- data.frame(d, sdz);
	}
	
	return (d);
}

##
## Create or replace in the data frame data the column data$S by S
##
## data : The data frame
## S : The new values for S
##
## return the modified data frame
##
lpl.R.dev.faceOutputAnalysis.replaceFactorS <- function(data, S) {

		cis = which(colnames(data) == "S");
		if (length(cis) == 0) {
			data <- data.frame(data, S);
		} else {
			data[, cis] <- S;
		}

		return (data);
}

##
## Compute the mean standard dispersion of residuals over the X and Y values
##
## data : The data frame containing the values of the residuals
##
## return the mean standard dispersion
##
lpl.R.dev.faceOutputAnalysis.computeMeanProjectedResidualsStandardDispersion <- function(data) {

	lix <- grep("rP..x", colnames(data), perl=TRUE, value=FALSE);
	liy <- grep("rP..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	sd <- 0;
	for (j in c(1:number_of_points)) {
		cix = lix[j];
		ciy = liy[j];
		sd <- sd + sd(data[, cix]);
		sd <- sd + sd(data[, ciy]);
	}

	return (sd/(2*number_of_points));
}

##
## Add a label S for point with to small or to large raw S factor
##
## data : The data frame
## sigma : The sigma factor column
## threshold : The threshold in sigma for removing outliers
##
## return the data frame with modified column label
##
lpl.R.dev.faceOutputAnalysis.addLabelSFactorOutlier <- function(data, sigma, threshold) {

	s <- sigma[!is.na(sigma) & as.character(data$label) == "Y"];
	mean <- mean(s);
	sd <- sd(s);

	for (i in c(1:nrow(data))) {
		if (!is.na(sigma[i]) & as.character(data$label[i]) == "Y") {
			if (((sigma[i]-mean)/sd < (-threshold)) | (threshold < (sigma[i]- mean)/sd)) {
				data$label[i]	<- "S";
			}
		}
	}

	return (data);
}

##
## Add a label L to points with large angles with respect to the mean angles of the head
##
## data : The data frame
## min_pitch_angle : The minimal pitch angle (with respect to the mean pitch of the head)
## max_pitch_angle : The maximal pitch angle (with respect to the mean pitch of the head)
## min_yaw_angle : The minimal yaw angle (with respect to the mean yaw of the head)
## max_yaw_angle : The maximal yaw angle (with respect to the mean yaw of the head)
## min_roll_angle : The minimal roll angle (with respect to the mean roll of the head)
## max_roll_angle : The maximal roll angle (with respect to the meanroll of the head)
##
## return the data frame with modified column label
##
lpl.R.dev.faceOutputAnalysis.addLabelForLargeAngle <- function(data, min_pitch_angle, max_pitch_angle, min_yaw_angle, max_yaw_angle, min_roll_angle, max_roll_angle) {

	d <- subset(data, data$label == "Y");
	mpitch <- mean(d$pitch);
	myaw <- mean(d$yaw);
	mroll <- mean(d$roll);

	for (i in c(1:nrow(data))) {
		if (data$label[i] == "Y") {
			if (((data$pitch[i]-mpitch) < min_pitch_angle) | ((data$pitch[i]-mpitch) > max_pitch_angle) | ((data$yaw[i]-myaw) < min_yaw_angle) | ((data$yaw[i]-myaw) > max_yaw_angle) | ((data$roll[i]-mroll) < min_roll_angle) | ((data$roll[i]-mroll) > max_roll_angle)) {
				data$label[i]	<- "L";
			}
		}
	}

	return (data);
}

##
## Add a label P to head position x and y coordinates too large
##
## data : The data frame
## threshold : The threshold in sigma for removing outliers
##
## return the modified data frame
##
lpl.R.dev.faceOutputAnalysis.addLabelOutlierHeadPosition <- function(data, threshold) {

	## The x coordinate of the center point is column 9
	dx <- as.numeric(as.character(data[, 9]));
	xmean <- mean(dx, na.rm=T);
	xsd <- sd(dx, na.rm=T);

	## The x coordinate of the center point is column 10 
	dy <- as.numeric(as.character(data[, 10]));
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

##
## Estimate the parameter of the head model for a given point
##
## data : The data frame
## dec : The coefficients for the direct estimation
## pointindex : The point index
## removeOutliers : Put this flag to TRUE to remove outliers
##
## return the data frame with direct head model parameters
##
lpl.R.dev.faceOutputAnalysis.estimateDirectXYZ <- function(data, dec, pointindex, removeOutliers) {

	pointname <- paste("P", pointindex, sep="");
	column_index_x <- which(colnames(data) == paste(paste("vM", pointindex, sep=""), "x", sep=""));
	column_index_y <- which(colnames(data) == paste(paste("vM", pointindex, sep=""), "y", sep=""));

	nbr_line <- 2*nrow(data);
	value <- numeric(nbr_line);

	for (i in c(1:nrow(data))) {
		value[2*i+1] <- data[i, column_index_x]/data$S[i];
		value[2*i+2] <- data[i, column_index_y]/data$S[i];
	}

	df <- data.frame(dec, value);

	if (removeOutliers) {
		df4 <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(df, 2, 3);
	} else {
		df4 <- df;
	}

	lm <- lm(df4$value ~ df4$cx + df4$cy + df4$cz + 0);

	X <- lm$coefficients[1];
	Y <- lm$coefficients[2];
	Z <- lm$coefficients[3];
	seX <- coef(summary(lm))[, 2][1];
	seY <- coef(summary(lm))[, 2][2];
	seZ <- coef(summary(lm))[, 2][3];

	df <- data.frame(pointname, X, seX, Y, seY, Z, seZ);
	names(df) <- c("point", "X", "seX", "Y", "seY", "Z", "seZ");
	rownames(df) <- NULL;

	return (df);
}

##
## Create the rotation coefficients for the direct estimation
##
## drc : The direct rotation coefficients
##
## return the data frame containing the coefficients
##
lpl.R.dev.faceOutputAnalysis.createCoefficientForDirectEstimation <- function(drc) {

	nbr_line <- 2*nrow(drc);
	cx <- numeric(nbr_line);
	cy <- numeric(nbr_line);
	cz <- numeric(nbr_line);

	for (i in c(1:nrow(drc))) {
		cx[2*i+1] <- drc$cxx[i];
		cy[2*i+1] <- drc$cxy[i];
		cz[2*i+1] <- drc$cxz[i];
		cx[2*i+2] <- drc$cyx[i];
		cy[2*i+2] <- drc$cyy[i];
		cz[2*i+2] <- drc$cyz[i];
	}

	df <- data.frame(cx, cy, cz);
	names(df) <- c("cx", "cy", "cz");
	rownames(df) <- NULL;

	return (df);
}

##
## Create the direct head model for the 11 points of the head
##
## data : The data frame
## drc : The direct rotation coefficients
## removeOutliers : Put this flag to TRUE to remove outliers
##
## return the data frame containing the head model parameters
##
lpl.R.dev.faceOutputAnalysis.createDirectHeadModel <- function(data, drc, removeOutliers) {

	dec <- lpl.R.dev.faceOutputAnalysis.createCoefficientForDirectEstimation(drc); 
	df <- NULL;

	li <- grep("vM..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(li);
	for (i in c(1:number_of_points)) {
		pointindex <- substring(colnames(data)[li[i]], 3,4);
		df <- rbind(df, lpl.R.dev.faceOutputAnalysis.estimateDirectXYZ(data, dec, pointindex, removeOutliers));
		cat(".");
	}
	cat("\n");
	
	return (df);
}

##
## Filter the outliers for a given column
##
## data : The data frame to be filtered
## column_index : The index of the column to be filtered
## alpha : The alpha parameter (in sigma unit) for the outliers rejection 
##
## The filtered data frame
##
lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn <- function(data, column_index, alpha) {

	previous_nol <- nrow(data);
	current_nol <- 0;

	while (previous_nol != current_nol) {
		v <- data[,column_index];
		meanv <- mean(v, na.rm=T);
		sdv <- sd(v, na.rm=T);
		df <- subset(data, (v-meanv)/sdv > -alpha & (v-meanv)/sdv < alpha);
		previous_nol <- nrow(data);
		data <- df
		current_nol <- nrow(data);
	}

	return (df);
}

##
## Merge two filters (logical columns) with the OR operation
##
## f1 : The first filter
## f2 : The second filter
##
## return the new filter
##
lpl.R.dev.faceOutputAnalysis.mergeFilter <- function(f1, f2) {

	f <- logical(length(f1));
	f <- f1 | f2;
	return (f);
}

##
## Create a filter for rejecting outliers (the filter is a column with
## LOGICAL value TRUE for outliers, FALSE otherwise)

## data : The data frame
## column_index : The column index of the value to filter
## alpha : The alpha value (values filtered are (x-mu)/sd < -alpha or 
## (x-mu)/sd > alpha where the mean mu and standard dispersion sd are 
## computed iteratively
##
## return the filtered data frame
##
lpl.R.dev.faceOutputAnalysis.createFilterOutliersForColumn <- function(data, column_index, alpha) {

	filter <- logical(nrow(data));
	filter <- rep(FALSE, nrow(data));

	previous_nol <- nrow(data);
	current_nol <- 0;

	df <- data;

	while (previous_nol != current_nol) {
		v <- ifelse(is.infinite(df[,column_index]), NA, df[,column_index]);
		meanv <- mean(v, na.rm=T);
		sdv <- sd(v, na.rm=T);
		filter <- (data[,column_index]-meanv)/sdv < -alpha | (data[,column_index]-meanv)/sdv > alpha;
		previous_nol <- nrow(df);
		df <- subset(df, (v-meanv)/sdv > -alpha & (v-meanv)/sdv < alpha);
		current_nol <- nrow(df);
	}

	return (filter);
}

##
## Return a two columns one line data frame which converts the time in minutes and seconds 
## from a time in seconds
##
## seconds : The time in seconds
##
## return the minutes and the seconds
##
lpl.R.dev.faceOutputAnalysis.stomns <- function(seconds) {
	
	mn <- numeric(1);
	s <- numeric(1);
	
	mn[1] <- floor(seconds/60);
	s[1] <- seconds - mn*60;

	d <- data.frame(mn, s);
	return(d);
}

##
## Return a String which converts the time in minutes and seconds from a time in seconds
##
## seconds : The time in seconds
##
## return a String with the minutes and the seconds
#"
lpl.R.dev.faceOutputAnalysis.stomnsString <- function(seconds) {
	
	mn <- numeric(1);
	s <- numeric(1);
	
	mn[1] <- floor(seconds/60);
	s[1] <- seconds - mn*60;

	return (paste(paste(round(mn[1], 2), " mn", sep=""), paste(round(s[1], 2), " s", sep="")));
}

##---------------------------------------------------------------------

##
## Create from the data frame with vP.* data, a subset filtered at a 2.5 sigma on the X and Y coordinate of the vector corresponding to pointname
## and add to the data frame the 9 coefficients of rotation in cr
##
filterResidualsOutliers <- function(data, alpha) {

	df <- data;
	for (j in c(1:12)) {
		if (j != 5) {
			column_index_x <- which(colnames(df) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
			column_index_y <- which(colnames(df) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));
			df <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(df, column_index_x, alpha);
			df <- lpl.R.dev.faceOutputAnalysis.filterOutliersForColumn(df, column_index_y, alpha);
		}
	}

	return (df);
}

lpl.R.dev.faceOutputAnalysis.meanDispersionInPixels <- function(software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	topprojectdir <- projectdir;
	projectdir <- paste(projectdir, "/", software, sep="");

	FOLDER_MODEL <- paste(projectdir, "model", sep="/");
        rrsdt <- loadInternalDataFrame(FOLDER_MODEL, "rrsdt.txt");

	return ((mean(rrsdt$sdx) + mean(rrsdt$sdy))/2);
}

##
## Compute the standard deviation by bin
##
computeStandardDeviationByBin <- function(data, axis, value, binfrontier) {

	nbrbin <- length(binfrontier) + 1;
	
	listStandardDeviation <- vector(mode="list", length=23);

	for (j in c(1:11)) {
		listStandardDeviation[[j]] <- numeric(nbrbin);
	}

	for (k in c(1:nbrbin)) {
			
		if (k == 1) {
			ds <- subset(data, value <= binfrontier[k]);
			vs <- subset(value, value <= binfrontier[k]); 
		} else if (k == nbrbin) {
			ds <- subset(data, value > binfrontier[k-1]);
			vs <- subset(value, value > binfrontier[k-1]); 
		} else {
			ds <- subset(data, value <= binfrontier[k] & value > binfrontier[k-1]);
			vs <- subset(value, value <= binfrontier[k] & value > binfrontier[k-1]); 
		}

		listStandardDeviation[[1]][k] <- mean(vs, na.rm=T);

		jc = 1;
		for (j in c(1:12)) {
			if (j != 5) {
				ci = which(colnames(ds) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), axis, sep=""));
				s <- ds[, ci];
				listStandardDeviation[[1+2*(jc-1)+1]][k] <- mean(s, na.rm=T);
				listStandardDeviation[[1+2*(jc-1)+2]][k] <- sd(s, na.rm=T);

				jc <- jc+1;
			}
		}
	}

	return (listStandardDeviation);
}

printResidualsStandardDispersion <- function(data) {

	textx <- "";
	texty <- "";
	jc = 1;
	for (j in c(1:12)) {
		if (j != 5) {
			cix = which(colnames(data) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "x", sep=""));
			ciy = which(colnames(data) == paste(paste("rP", lpl.R.dev.faceOutputAnalysis.getStringNumber(j), sep=""), "y", sep=""));
			textx <- paste(textx, sd(data[, cix]));
			texty <- paste(texty, sd(data[, ciy]));
			jc <- jc+1;
		}
	}
	print(textx);
	print(texty);
}

##
## Create bin boundaries for the column columname (the bins will contain the
## same number of points)
## data : The data frame containing the column
## columname : The name of the column
## nbrbin : The number of bins to create
##
createBinBoundaries <- function(data, columnname, nbrbin) {

	ci <- which(colnames(data) == columnname);
	cdf <- ecdf(data[, ci]);

	binfrontiers <- numeric(nbrbin-1);

	frontierindex <- 1;
	cdffrontiervalue <- frontierindex/nbrbin;

       	for (i in c(1:length(environment(cdf)$y))) {
		if (environment(cdf)$y[i] > cdffrontiervalue) {
			binfrontiers[frontierindex] <- environment(cdf)$x[i];
			frontierindex <- frontierindex + 1;
			cdffrontiervalue <- frontierindex/nbrbin;
		}
	}

	return (binfrontiers);
}



searchMaximalPitch <- function(data, time, dhmt) {

	lix <- grep("vM..x", colnames(data), perl=TRUE, value=FALSE);
	liy <- grep("vM..y", colnames(data), perl=TRUE, value=FALSE);
	number_of_points <- length(lix);

	frame_index = time*25+1;
	pitch_0 <- data$pitch[frame_index];
	yaw_0 <-  data$yaw[frame_index];
	roll_0 <- data$roll[frame_index];

	cat(paste(time, pitch_0, yaw_0, roll_0, "\n"));

	delta_pitch = 1;

	r22y <- 0;
	r23y <- 0;
	r22ymax <- 0;
	r23ymax <- 0;
	pitchmax <- 0;
	dispersionmin <- 1000000000000;

	for (i in c(1:40)) {	

		pitch <- pitch_0 + (i-20)*delta_pitch;
		rc <- lpl.R.dev.faceOutputAnalysis.createRotationCoefficients(pitch, yaw_0, roll_0,"inverse");	
	
		dispersionx <- 0;
		dispersiony <- 0;

		for (j in c(1:number_of_points)) {
		

			xobs <- data[frame_index, lix[j]]/data$S[frame_index];
			yobs <- data[frame_index, liy[j]]/data$S[frame_index];
	
			Xobs = (rc[1,1] - rc[1,3]*rc[3,1]/rc[3,3])*xobs + (rc[1,2] - rc[1,3]*rc[3,2]/rc[3,3])*yobs;
			Yobs = (rc[2,1] - rc[2,3]*rc[3,1]/rc[3,1])*xobs + (rc[2,2] - rc[2,3]*rc[3,2]/rc[3,3])*yobs;

			dispersionx <- dispersionx + (Xobs - dhmt$X[j] + rc[1,3]/rc[3,3]*dhmt$Z[j])^2;
			dispersiony = dispersiony + (Yobs - dhmt$Y[j] + rc[2,3]/rc[3,3]*dhmt$Z[j])^2;

		
			if (j==23) {
				r23y <- Yobs - dhmt$Y[j] + rc[2,3]/rc[3,3]*dhmt$Z[j];
			}
			if (j==22) {
				r22y <- Yobs - dhmt$Y[j] + rc[2,3]/rc[3,3]*dhmt$Z[j];
			}

			
		}
		dispersionx <- sqrt(dispersionx/number_of_points);
		dispersiony = sqrt(dispersiony/number_of_points);

		if (dispersiony < dispersionmin) {
			dispersionmin <- dispersiony;
			r22ymax <- r22y;
			r23ymax <- r23y;
			pitchmax <- pitch;
		}

		##cat(paste(i, pitch, dispersionx, dispersiony, "\n"));
	}
	cat(paste(pitchmax, dispersionmin, r22ymax,r23ymax,  "\n"));
}


