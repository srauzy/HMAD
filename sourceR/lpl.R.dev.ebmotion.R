## lpl.R.dev.ebmotion.R, creator S. Rauzy, LPL, 05/10/2016



lpl.R.dev.ebmotion.insertBadQualityAreaWithZscore <- function(dfres, qf, qt, tmin, tmax, ebmotion) {
	
	bqaopen <- FALSE;
	imin = round((tmin - qf$time[1])*25 + 1, 0);
	imax = round((tmax - qf$time[1])*25 + 1, 0);

	cat(nrow(qf$quality));
	for (i in c(imin:imax)) {
		if (qf$quality[i] < qt) {
			if (!bqaopen) {
				bqaopen <- TRUE;
				bqatmin <- round(qf$time[i], 2);
				mquality <- 0;
				nbrmeasure <- 0;
			}
			mquality <- mquality + qf$quality[i];
			nbrmeasure <- nbrmeasure + 1;
		} else {
			if (bqaopen) {
				bqatmax <- round(qf$time[i], 2);
				duration <- round(bqatmax-bqatmin, 2); 
				dff <- data.frame(bqatmin, bqatmax, duration, NA, NA, (mquality/nbrmeasure), "X", ebmotion);
				names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "class", "ebmotion");
				rownames(dff) <- NULL;
				dfres <- rbind(dfres, dff);
				bqaopen <- FALSE;
			}
		}	
	}

	if (bqaopen) {
		bqatmax <- round(qf$time[i], 2);
		duration <- round(bqatmax-bqatmin, 2); 
		dff <- data.frame(bqatmin, bqatmax, duration, NA, NA, (mquality/nbrmeasure), "X", ebmotion);
		names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "class", "ebmotion");
		rownames(dff) <- NULL;
		dfres <- rbind(dfres, dff);
		bqaopen <- FALSE;
	}
	return (dfres);
}

lpl.R.dev.ebmotion.addBadQualityAreaWithZscore <- function(dlma, qf, qt, ebmotion) {

	tstart <- round(qf$time[1], 2);
	tend <- round(qf$time[nrow(qf)], 2);

	nrow_dlma = nrow(dlma);

	dfres <- NULL;

	if (tstart < round(dlma$tmin[1], 2)) {
		tmin <- tstart;
		tmax <- round(dlma$tmin[1], 2);
		dfres <- lpl.R.dev.ebmotion.insertBadQualityAreaWithZscore(dfres, qf, qt, tmin, tmax, ebmotion); 
	}

	for (i in c(1:nrow_dlma)) {

		dfres <- rbind(dfres, dlma[i, ]);
		
		if (i != nrow_dlma &  round(dlma$tmax[i], 2) < round(dlma$tmin[i+1], 2)) {
			tmin <- round(dlma$tmax[i], 2);
			tmax <- round(dlma$tmin[i+1], 2);
			
			dfres <- lpl.R.dev.ebmotion.insertBadQualityAreaWithZscore(dfres, qf, qt, tmin, tmax, ebmotion); 
		}
	}

	if (tend > round(dlma$tmax[nrow(dlma)], 2)) {
		tmin <- round(dlma$tmax[nrow(dlma)], 2);
		tmax <- tend;
		dfres <- lpl.R.dev.ebmotion.insertBadQualityAreaWithZscore(dfres, qf, qt, tmin, tmax, ebmotion); 
	}

	return (dfres);
}

lpl.R.dev.ebmotion.addNullDetectionAreaWithZscore <- function(dlma, df, ebmotion) {

	tstart <- round(df$time[1], 2);
	tend <- round(df$time[nrow(df)], 2);

	nrow_dlma = nrow(dlma);

	dfres <- NULL;

	if (tstart < dlma$tmin[1]) {
		tmin <- tstart;
		tmax <- round(dlma$tmin[1], 2);
		duration <- round(tmax-tmin, 2); 
		dff <- data.frame(tmin, tmax, duration, NA, NA, NA, "N", ebmotion);
		names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "class", "ebmotion");
		rownames(dff) <- NULL;
		dfres <- rbind(dfres, dff);
	}

	for (i in c(1:nrow_dlma)) {

		dfres <- rbind(dfres, dlma[i, ]);
		
		if (i != nrow_dlma &  round(dlma$tmax[i], 2) < round(dlma$tmin[i+1], 2)) {
			tmin <- round(dlma$tmax[i], 2);
			tmax <- round(dlma$tmin[i+1], 2);
			duration <- round(tmax-tmin, 2); 
			dff <- data.frame(tmin, tmax, duration, NA, NA, NA, "N", ebmotion);
			names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "class", "ebmotion");
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);
		}
	}

	if (tend > round(dlma$tmax[nrow(dlma)], 2)) {
		tmin <- round(dlma$tmax[nrow(dlma)], 2);
		tmax <- tend;
		duration <- round(tmax-tmin, 2); 
		dff <- data.frame(tmin, tmax, duration, NA, NA, NA, "N", ebmotion);
		names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "class", "ebmotion");
		rownames(dff) <- NULL;
		dfres <- rbind(dfres, dff);
	}

	return (dfres);
}

lpl.R.dev.ebmotion.createLabelScoreTable <- function() {

	score <- character(7);
       	score[1] = "A";	
       	score[2] = "B";	
       	score[3] = "C";	
       	score[4] = "D";	
       	score[5] = "E";	
       	score[6] = "N";	
       	score[7] = "X";

	return (data.frame(score));
}

##
## Annotate the media with time intervals encoded with score (from A to E for probability detection of motion),
## N for no motion and X for problematic areas to check by hand. 
## dfwc : The selected wavelet coefficients potentially maxima 
## df1Dqf : The quality table at some reference scale
## qualitythreshold : The quality threshold of the wavelet support
## minimalscale : The minimal scale duration for selecting solution
## ebmotion : The type of eyebrow motion "Raise" or "Frown"
##
lpl.R.dev.ebmotion.automaticAnnotationEBM <- function(dfwc, df1Dqf, qualitythreshold, ebmotion) {

	## Create the score classes
	dfscoredef <- lpl.R.dev.ebmotion.createLabelScoreTable();

	cat(paste("--- Create automatic annotations for", ebmotion,  "function ---\n"));
	cat("Annotate class A, B, C, D and E eyebrow motions...\n");
	d <- lpl.R.dev.ebmotion.computeLocalMaximaByZscore(dfwc, dfscoredef); 
	d$ebmotion <- rep(ebmotion, nrow(d));

	s <- subset(d, d$class == "A");
	cat(paste("Number of A intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "B");
	cat(paste("Number of B intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "C");
	cat(paste("Number of C intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "D");
	cat(paste("Number of D intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "E");
	cat(paste("Number of E intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	

	cat("----------------\n");
	cat("Add bad quality intervals (class X)...\n");
	
	d <- lpl.R.dev.ebmotion.addBadQualityAreaWithZscore(d, df1Dqf, qualitythreshold, ebmotion);
	cat(paste("Quality threshold =", qualitythreshold, "\n"));
	s <- subset(d, d$class == "X");
	cat(paste("Number of X intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));

	cat("Add intervals presumed without eyebrow motion (class N)...\n");
	
	d <- lpl.R.dev.ebmotion.addNullDetectionAreaWithZscore(d, df1Dqf, ebmotion);
	
	s <- subset(d, d$class == "N");
	cat(paste("Number of N intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));
	cat("----------------\n");

	d <- lpl.R.dev.ebmotion.checkAndCorrectIntervalOverlap(d, 1, 2, FALSE);

	return (d);
}

##
## Annotate the media with time intervals encoded with score (from A to E for probability detection of motion),
## N for no motion and X for problematic areas to check by hand. 
## s2df : The 2D wavelet coefficients and quality
## sdwcqmodel : The support quality model
## qt : The quality threshold
## minimalscale : The minimal scale duration for selecting solution
## referencetimescale : The reference scale to build the 1D support
## useZscore : Flag to put to TRUE (because it gives better results)
## ebmotion : The type of eyebrow motion "Raise" or "Frown"
##
automaticAnnotationEBM1 <- function(s2df, sdwcqmodel, qt, minimalscale, referencetimescale, useZscore, ebmotion) {

	score <- character(7);
       	score[1] = "A";	
       	score[2] = "B";	
       	score[3] = "C";	
       	score[4] = "D";	
       	score[5] = "E";	
       	score[6] = "N";	
       	score[7] = "X";

	cat("----------------\n");
	cat("Compute wavelet coefficient thresholds using the quality model...\n");

	cat(paste("Minimal duration parameter for eyebrow motion detection =", (2*minimalscale), "seconds", "\n"));

	wcqt <- sdwcqmodel$intercept[1] + sdwcqmodel$slope[1]*getQualityCodeColumn(s2df$quality);

	cat("Annotate class A, B, C, D and E eyebrow motions...\n");
	s2dfa <- subset(s2df, s2df$wcoef > wcqt & s2df$quality > qt);

	if (useZscore) {
		dlma <- searchLocalMaximaByZscore(s2dfa, sdwcqmodel, minimalscale, ebmotion);
	} else {
		dlma <- searchLocalMaxima(s2dfa, minimalscale);
	}

	class <- character(nrow(dlma));

	for (i in c(1:nrow(dlma))) {
		wcqti <- sdwcqmodel$intercept[1] + sdwcqmodel$slope[1]*getQualityCode(dlma$quality[i]);
		class[i] <- score[getScoreIndex(dlma$wcoef[i], wcqti)];
	}

	d <- data.frame(dlma, class);
	
	s <- subset(d, d$class == "A");
	cat(paste("Number of A intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "B");
	cat(paste("Number of B intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "C");
	cat(paste("Number of C intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "D");
	cat(paste("Number of D intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	
	s <- subset(d, d$class == "E");
	cat(paste("Number of E intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));	

	cat("----------------\n");
	cat("Create 1D quality function...\n");
	qf <- create1DQualityFunction(s2df, referencetimescale);

	cat("----------------\n");
	cat("Add bad quality intervals (class X)...\n");
	if (useZscore) {
		d <- addBadQualityAreaWithZscore(d, qf, qt, ebmotion);
	} else {
		d <- addBadQualityArea(d, qf, qt);
	}
	cat(paste("Quality threshold =", qt, "\n"));
	s <- subset(d, d$class == "X");
	cat(paste("Number of X intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));

	cat("Add intervals presumed without eyebrow motion (class N)...\n");
	if (useZscore) {
		d <- addNullDetectionAreaWithZscore(qf, d, ebmotion);
	} else {
		d <- addNullDetectionArea(qf, d);
	}
	s <- subset(d, d$class == "N");
	cat(paste("Number of N intervals =", nrow(s), ", spanning", round(sum(s$duration), 2), "seconds\n"));
	cat("----------------\n");

	d <- checkAndCorrectIntervalOverlap(d, 1, 2);

	return (d);
}

lpl.R.dev.ebmotion.checkAndCorrectIntervalOverlap <- function(d, tmin_column_index, tmax_column_index, verbose) {

	for (i in c(2:nrow(d))) {
		if (round(d[i,tmin_column_index],2) != round(d[i-1,tmax_column_index], 2)) {
			if (verbose) {
				cat(paste("Problem at line", i, "tmin = ", d[i,tmin_column_index], "and tmax =", d[i-1,tmax_column_index], "\n"));
			}
			if (d$duration[i] > d$duration[i-1]) {
				d[i,tmin_column_index] <-  d[i-1,tmax_column_index];
				d$duration[i] <- round(d$tmax[i] - d$tmin[i], 2);
				if (verbose) {
					cat(paste("tmin changed in ", d[i,tmin_column_index], "\n"));
				}
			} else {
				d[i-1,tmax_column_index] <- d[i,tmin_column_index];
				d$duration[i-1] <- round(d$tmax[i-1] - d$tmin[i-1], 2);
				if (verbose) {
					cat(paste("tmax changed in ", d[i-1,tmax_column_index], "\n"));
				}
			}
		}
	}
	return (d);
}

##
## Return the annotation score index (from 0 to 5) 
## wcoef : The wavelet coefficient
## wcqt : 
##
getScoreIndex <- function(wcoef, wcqt) {

	for (i in c(1:5)) {
	     alpha = 5-(i-1);
	     if (wcoef > alpha*wcqt) {
		     return (i);
	     }
	}
	return (0);
}

##
## Return the annotation score index (from 0 to 5) 
## zscore : The zscore 
##
getScoreIndexForZscore <- function(zscore) {

	for (i in c(1:5)) {
	     alpha = 5-(i-1);
	     if (zscore > alpha) {
		     return (i);
	     }
	}
	return (0);
}

searchMaximalValue <- function(d, v, sigmamin, sigmamax) {

	number_of_points <- nrow(d);
	wcoefmax <- numeric(number_of_points);
	durationmax <- numeric(number_of_points);
	time <- numeric(number_of_points);

	d_step <- floor((sigmamax-sigmamin)*25) + 1;

	for (i in c(1:number_of_points)) {

		wcoefmax[i] <- NA;
		
		time[i] <- d$time[i];
		for (j in c(1:d_step)) {
			sigma <- sigmamin + (j-1)*0.04;
			imin = max(1, (i - (3*sigma)*25));
			imax = min(number_of_points, (i + (3*sigma)*25));
			wcoef <- waveletCoefAtTime(d$time, v, sigma, i, imin, imax)[1];
		        if (is.na(wcoefmax[i]) | wcoef > wcoefmax[i]) {
				wcoefmax[i] <- wcoef;
				durationmax[i] <- sigma;
			}	
		}
		##print(paste(time[i], wcoefmax[i], durationmax[i]));
	}
	
	dd <- data.frame(time, wcoefmax, durationmax);
	return (dd);
}

create1DQualityFunction1 <- function(s2d, scale) {

	d <- subset(s2d, round(s2d$duration, 2) == scale);

	time <- round(d$time, 2);
	quality <- d$quality;

	res <- data.frame(time, quality);
	rownames(res) <- NULL;

	return (res);
}

##
## From a 2D support, create a 1D support where for each time, the quality is 
## the mean quality over scale
## s2d : The 2D support
##
create1DSupport <- function(s2d) {

	l <- levels(factor(s2d$duration));
	ns <- length(l);
	nt <- nrow(s2d)/ns;

	time <- numeric(nt);
	quality <- numeric(nt);

	starting_time = getTimeInHundredthSecond();

	for (i in c(1:nt)) {
		d <- s2d[((i-1)*ns+1):(i*ns),];
		time[i] <- d$time[1];
		quality[i] <- (mean(d$qualityL) + mean(d$qualityC) + mean(d$qualityR))/3;
		if (i == 25*60-1) {
			ed <- estimateDurationInSecond(durationInSecond(starting_time, getTimeInHundredthSecond()), (nt/25/60));
			cat(paste("ESTIMATED COMPUTING DURATION :", lpl.R.dev.faceOutputAnalysis.stomnsString(ed), "\n"));
		}

		if (floor(i/25/60) == (i/25/60)) {
			print(paste("Frame", i, "over", nt));
		}
	}
	return (data.frame(time, quality));
}

##
## Return the quality code (from 0 to 11) from the quality value (from 0 to 1)
## quality : The quality (from 0 to 1)
##
getQualityCode <- function(quality) {
	if (is.na(quality)) {
		return (NA);
	}
	if (quality == 0) {
		return (0);
	}
	for (i in c(1:9)) {
		if (quality > (i-1)/10 & quality <= i/10) {
			return (i);
		}
	}
	if (quality > 9/10 & quality < 1) {
		return (10);
	} else {
		return (11);
	}
}

##
## Return a column with the quality code from the quality column
## qualityColumn : The quality column (values spanning [0, 1] interval
##
getQualityCodeColumn <- function(qualityColumn) {

	res <- numeric(length(qualityColumn));	
	for (i in c(1:length(res))) {
		res[i] <- getQualityCode(qualityColumn[i]);
	}
	return (res);
}


##
## Create the quality model giving the amplitude of noise in function
## of the quality of the support
## sd2 : The 2D support and wavelet coefficients
## qt : The quality threshold for computing the linear model
##
createStandardDeviationQualityModel <- function(qdm, qt) {

	cat("Quality code linear model for wavelet coefficients standard deviation :\n");
      	qct <- getQualityCode(qt);	
	
	qdms <- subset(qdm, qdm$qc > qct);

	qdms$sd <- sqrt(qdms$v/qdms$n - (qdms$m/qdms$n)^2);
	qdms$weight <- sqrt(qdms$n);

	lmw <- lm(qdms$sd ~ qdms$qc, weight=qdms$weight);

	intercept <- lmw$coefficients[[1]]; 	
	slope <- lmw$coefficients[[2]];
	cat(paste("Intercept =", round(intercept, 4), ",slope =", round(slope, 4), "\n"));

	return (data.frame(intercept, slope));		
}

##
## Create the quality model giving the amplitude of noise in function
## of the quality of the support
## sd2 : The 2D support and wavelet coefficients
## qt : The quality threshold for computing the linear model
##
createQualityModel <- function(s2d, qt) {

	cat("Compute standard deviations (in the wavelet space) in function of the quality...\n");  
	sdwcq <- createWaveletCoefficientStandardDeviationByQualityAndScale(s2d);

	##qplot(quality, sdwc, data = sdwcq) + geom_point(size=2, aes(colour = duration));

	cat("Compute quality model for standard deviations...\n");  
	cat(paste("Quality threshold =", qt, "\n"));  
	qm <- createWaveletCoefficientThresholdQualityModel(sdwcq, qt);
	cat(paste(paste("Intercept =", qm$intercept[1], ", slope =", qm$slope[1]), "\n"));  

	return (qm);	
}

##
## Create from the table of the standard deviation of wavelet coefficients
## the linear model for the standard deviation
## sdwcq : The 2D grid in scale and quality the standard deviation of wavelet coefficients
## qt : The quality threshold for computing the linear model
##
createWaveletCoefficientThresholdQualityModel <- function(sdwcq, qt) {

	qct <- getQualityCode(qt);
	sdwcqf <- subset(sdwcq, sdwcq$quality > qct);

	lm <- lm(sdwcqf$sdwc ~ sdwcqf$quality);
	intercept <- lm$coefficients[[1]]; 	
	slope <- lm$coefficients[[2]];

	return (data.frame(intercept, slope));	
}

##
## Compute on a 1D grid in  quality the standard deviation of wavelet coefficients
## s2d : The 2D support and wavelet coefficients
##
createWaveletCoefficientStandardDeviationByQuality <- function(s2d) {

	squality <- vector(mode="list", length=11);
	for (i in c(1:9)) {
		squality[[i]] <- subset(s2d, s2d$quality > (i-1)/10 & s2d$quality <= i/10);
	}
	squality[[10]] <- subset(s2d, s2d$quality > 9/10 & s2d$quality < 1);
	squality[[11]] <- subset(s2d, s2d$quality == 1);

	nq <- length(squality);
	quality <- numeric(nq);
	sdwc <- numeric(nq);
	n <- numeric(nq);

	for (i in c(1:nq)) {
		ss <- squality[[i]];
		sdwc[i] <- sd(ss$wcoef);
		quality[i] <- i;
		n[i] <- nrow(ss);
	}
	return (data.frame(quality, sdwc, n));
}

##
## Compute on a 2D grid in scale and quality the standard deviation of wavelet coefficients
## s2d : The 2D support and wavelet coefficients
##
createWaveletCoefficientStandardDeviationByQualityAndScale <- function(s2d) {

	l <- levels(factor(s2d$duration));
	ns <- length(l);

	squality <- vector(mode="list", length=11);
	for (i in c(1:9)) {
		squality[[i]] <- subset(s2d, s2d$quality > (i-1)/10 & s2d$quality <= i/10);
	}
	squality[[10]] <- subset(s2d, s2d$quality > 9/10 & s2d$quality < 1);
	squality[[11]] <- subset(s2d, s2d$quality == 1);

	nq <- length(squality);
	quality <- numeric(nq*ns);
	duration <- numeric(nq*ns);
	sdwc <- numeric(nq*ns);
	intercept <- numeric(nq*ns);
	n <- numeric(nq*ns);

	cat("----------------\n");
	for (i in c(1:nq)) {
		cat(paste("Quality", i, ": mean quality =", round(mean(squality[[i]]$quality), 3), ", nbr coefficients = ", nrow(squality[[i]]), "\n"));
		for (j in c(1:ns)) {
			ss <- subset(squality[[i]], squality[[i]]$duration == l[j]);
			sdwc[ns*(i-1)+j] <- sd(ss$wcoef);
			duration[ns*(i-1)+j] <- as.numeric(l[j]);
			quality[ns*(i-1)+j] <- i;
			n[ns*(i-1)+j] <- nrow(ss);
		}
	}
	cat("----------------\n");
	return (data.frame(quality, duration, sdwc, n));
}


##
## Return the duration of two time intervals overlap (0 if the
## interval do not overlap)
##
## tmin1 : the left boundary of interval 1
## tmax1 : the right boundary of interval 1
## tmin2 : the left boundary of interval 2
## tmax2 : the right boundary of interval 2
##
## return TRUE if the two intervals overlap
##
intervalOverlap <- function(tmin1, tmax1, tmin2, tmax2) {

	if (tmax1 <= tmin2) return (0);
	if (tmin1 >= tmax2) return (0);
	
	if (tmin1 <= tmin2) {
		if (tmax1 >= tmax2) {
			return (tmax2 - tmin2);
		} else {
			return (tmax1 - tmin2);
		}	
	} else {
		if (tmax2 >= tmax1) {
			return (tmax1 - tmin1);
		} else {
			return (tmax2 - tmin1);
		}
	}
}

##
## Create a 2D support (time and scale) and wavelet coefficients of the
## function v 
## d : The data frame containing the times
## v : The first function for computing the wavelet coefficients
## filter : The filter (value at TRUE are put to NA)
## sigmamin : The minimal scale duration
## sigmamax : The maximal scale duration
## sigmastep : The step in scale
## wDparameter : The parameter tuning the mother wavelet
##
lpl.R.dev.ebmotion.create2DSupport1 <- function(d, v, ev, filter, sigmamin, sigmamax, sigmastep, wDparameter) {

	v[filter == TRUE] <- NA;

	nt <- nrow(d);
	line_number = nt;
	tmin = d$time[1];
	tmax = d$time[nt];
	ns <- (sigmamax-sigmamin)/sigmastep;
	qualityL <- numeric(nt*ns);
	qualityC <- numeric(nt*ns);
	qualityR <- numeric(nt*ns);
	time <- numeric(nt*ns);
	duration <- numeric(nt*ns);
	wcoef <- numeric(nt*ns);
	errm <- numeric(nt*ns);
	errsd <- numeric(nt*ns);

	starting_time = getTimeInHundredthSecond();

	for (j in c(1:nt)) {
		for (i in c(1:ns)) {
			sigma <- sigmamin + (i-1)*sigmastep;
			duration[ns*(j-1) + i] <- sigma;
			time[ns*(j-1) + i] <- d$time[j];
			jmin = round(max(1, (j+1 - ((wDparameter+1)*sigma)*25)), 0);
			jmax = round(min(nt, (j+1 + ((wDparameter+1)*sigma)*25)), 0);
			wc <- lpl.R.dev.wavelets.waveletCoefAtTime(d$time, v, sigma, j, jmin, jmax, wDparameter);
			qualityL[ns*(j-1) + i] <- wc[2]; 
			qualityC[ns*(j-1) + i] <- wc[3]; 
			qualityR[ns*(j-1) + i] <- wc[4];
			wcoef[ns*(j-1) + i] <- wc[1];
			errv <- errorAmplitudeOnInterval(ev, jmin, jmax);
		       	errm[ns*(j-1) + i] <- errv[1];
		       	errsd[ns*(j-1) + i] <- errv[2];
		}
		if (j == 25*20-1) {
			ed <- estimateDurationInSecond(durationInSecond(starting_time, getTimeInHundredthSecond()), (nt/25/20));
			cat(paste("PREDICTED COMPUTING TIME :", lpl.R.dev.faceOutputAnalysis.stomnsString(ed), "\n"));
		}

		if (floor(j/25/20) == (j/25/20)) {
			cat(paste("Frame", j, "over", line_number, "\n"));
		}
	}

	quality <- numeric(nt*ns);
	quality <- (qualityL + qualityC + qualityR)/3;

	return (data.frame(time, duration, wcoef, qualityL, qualityC, qualityR, quality, errm, errsd));
}

##
## Create a 2D support (time and scale) and wavelet coefficients of the
## function v 
## d : The data frame containing the times
## v : The first function for computing the wavelet coefficients
## filter : The filter (value at TRUE are put to NA)
## sigmamin : The minimal scale duration
## sigmamax : The maximal scale duration
## sigmastep : The step in scale
## wDparameter : The parameter tuning the mother wavelet
##
lpl.R.dev.ebmotion.create2DSupport <- function(d, v, filter, sigmamin, sigmamax, sigmastep, wDparameter) {

	v[filter == TRUE] <- NA;

	nt <- nrow(d)-1;
	line_number = nt;
	tmin = d$time[1]+0.02;
	tmax = d$time[nt]+0.02;
	ns <- (sigmamax-sigmamin)/sigmastep;
	qualityL <- numeric(nt*ns);
	qualityC <- numeric(nt*ns);
	qualityR <- numeric(nt*ns);
	time <- numeric(nt*ns);
	duration <- numeric(nt*ns);
	wcoef <- numeric(nt*ns);

	starting_time = getTimeInHundredthSecond();

	for (j in c(1:nt)) {
		for (i in c(1:ns)) {
			sigma <- sigmamin + (i-1)*sigmastep;
			duration[ns*(j-1) + i] <- sigma;
			time[ns*(j-1) + i] <- d$time[j]+0.02;
			jmin = round(max(1, (j + 1 - (wDparameter+1)*sigma/0.08)), 0);
			jmax = round(min(nt, (j + 1 + (wDparameter+1)*sigma/0.08)), 0);
			wc <- lpl.R.dev.wavelets.waveletCoefAtTime(d$time, v, sigma, j, jmin, jmax, wDparameter);
			qualityL[ns*(j-1) + i] <- wc[2]; 
			qualityC[ns*(j-1) + i] <- wc[3]; 
			qualityR[ns*(j-1) + i] <- wc[4];
			wcoef[ns*(j-1) + i] <- wc[1];
		}
		if (j == 25*20-1) {
			ed <- estimateDurationInSecond(durationInSecond(starting_time, getTimeInHundredthSecond()), (nt/25/20));
			cat(paste("PREDICTED COMPUTING TIME :", lpl.R.dev.faceOutputAnalysis.stomnsString(ed), "\n"));
		}

		if (floor(j/25/20) == (j/25/20)) {
			cat(paste("Frame", j, "over", line_number, "\n"));
		}
	}

	quality <- numeric(nt*ns);
	quality <- (qualityL + qualityC + qualityR)/3;

	return (data.frame(time, duration, wcoef, qualityL, qualityC, qualityR, quality));
}

lpl.R.dev.ebmotion.initializeQualityDataModel <- function() {

	qc <- c(0:11);
	n <- rep(0, 12);
	m <- rep(0, 12);
	v <- rep(0, 12);

	qdm <- data.frame(qc, n, m, v);

	return (qdm);
}

lpl.R.dev.ebmotion.fillQualityDataModel <- function(qdm, d) {

	qc <- getQualityCodeColumn(d$quality);
	for (i in c(1:12)) {
		dqc <- subset(d, qc == (i-1));
		qdm$n[i] <- qdm$n[i] + nrow(dqc);
		qdm$m[i] <- qdm$m[i] + sum(dqc$wc);
		qdm$v[i] <- qdm$v[i] + sum(dqc$wc^2);
	}
	return (qdm);
}

##
## Create the quality function of the wavelet transform for some reference scale
##
## d : The data frame containing the times
## filter : The filter (value at TRUE are put to 0
## sigmareference : The reference scale for computing the quality function
## wDparameter : The parameter tuning the mother wavelet
##
## return a two column data fame with time and quality
##
lpl.R.dev.ebmotion.create1DQualityFunction <- function(d, filter, sigmareference, wDparameter) {

	w <- numeric(nrow(d));
	w <- ifelse(filter == TRUE | is.na(filter), 0, 1);
	f <- rep(0, nrow(d));
	nd <- data.frame(f, w);

	isigma <- round(sigmareference/0.08);
	res <- lpl.R.dev.wavelets.iwaveletatscale(nd, isigma, wDparameter);

	time <- d$time;
	quality <- res$q;

	r <- data.frame(time, quality);

	return (r);
}

##
## Create a 2D support (time and scale) and wavelet coefficients of the ## function f 
## d : The data frame containing the times
## f : The function for computing the wavelet coefficients
## filter : The filter (value at TRUE are put to 0)
## sigmamin : The minimal scale duration
## sigmamax : The maximal scale duration
## sigmastep : The step in scale
## wDparameter : The parameter tuning the mother wavelet
##
lpl.R.dev.ebmotion.createWaveletTable <- function(d, f, filter, qt, sigmamin, sigmamax, sigmastep, wDparameter) {

	w <- numeric(nrow(d));
	w <- ifelse(filter == TRUE | is.na(filter), 0, 1);
	f[filter == TRUE | is.na(filter)] <- 0;
	nd <- data.frame(f, w); 

	ns <- round((sigmamax-sigmamin)/sigmastep, 2)+1;
		
	starting_time = getTimeInHundredthSecond();
	rt <- NULL;
	
	qdm <- lpl.R.dev.ebmotion.initializeQualityDataModel();

	for (i in c(1:ns)) {
		sigma <- sigmamin + (i-1)*sigmastep;
		
		isigma <- round(sigma/0.08);
		res <- lpl.R.dev.wavelets.iwaveletatscale(nd, isigma, wDparameter);

		tmin <- d$time+0.02-0.04*isigma;
	        tmax <- d$time+0.02+0.04*isigma;
		wc <- res$wc;
		quality <- res$q;

		r <- data.frame(tmin, tmax, wc, quality);

		qdm <- lpl.R.dev.ebmotion.fillQualityDataModel(qdm, r);

		r <- subset(r, r$wc > 0 & quality > qt);
		rownames(r) <- NULL;
		rt <- rbind(rt, r);

		if (i == 1) {
			ed <- estimateDurationInSecond(durationInSecond(starting_time, getTimeInHundredthSecond()), ns);
			cat(paste("Wavelet coefficient for scale parameters from ", sigmamin, "to ", sigmamax, ":\n"));
			cat(paste("PREDICTED COMPUTING TIME :", lpl.R.dev.faceOutputAnalysis.stomnsString(ed), "\n"));
			cat("s =");
		}
		cat(paste(" ", sigma));	
	}
	cat("\n");	

	sdqm <- createStandardDeviationQualityModel(qdm, qt);
	vsdqm <- sdqm$intercept[1] + sdqm$slope[1]*getQualityCodeColumn(rt$quality);
	rt$zscore <- rt$wc/vsdqm;

	rts <- subset(rt, rt$zscore > 1);

	return(rts);	
}

##
## Create a 2D support (time and scale) and wavelet coefficients of the ## function f 
## d : The data frame containing the times
## f : The function for computing the wavelet coefficients
## filter : The filter (value at TRUE are put to 0)
## sigmamin : The minimal scale duration
## sigmamax : The maximal scale duration
## sigmastep : The step in scale
## wDparameter : The parameter tuning the mother wavelet
##
lpl.R.dev.ebmotion.createWaveletCoefficients <- function(t, f, sigmamin, sigmamax, sigmastep, wDparameter) {

	w <- rep(1: length(f));
	f[is.na(f)] <- 0;
	nd <- data.frame(f, w); 

	ns <- round((sigmamax-sigmamin)/sigmastep, 2)+1;
		
	rt <- NULL;
	
	for (i in c(1:ns)) {
		sigma <- sigmamin + (i-1)*sigmastep;
		
		isigma <- round(sigma/0.08);
		res <- lpl.R.dev.wavelets.iwaveletatscale(nd, isigma, wDparameter);

		tmin <- t+0.02-0.04*isigma;
	        tmax <- t+0.02+0.04*isigma;
		wc <- res$wc;

		r <- data.frame(round((tmax+tmin)/2, 2), round((tmax-tmin), 2), wc);

		rownames(r) <- NULL;
		rt <- rbind(rt, r);
	}
	colnames(rt) <- c("t", "duration", "wcoef");

	return(rt);	
}

##
## Compute the local maxima from the zcore value under the constrains that 
## the resulting intervals have no overlaps
## 
## res : The data frame containing tmin, tmax and zscore for each interval
## dfscoredef : The data frame containing the class definition depending of the zscore value
##
## return the data frame containing the local maxima intervals with no overlap
##
lpl.R.dev.ebmotion.computeLocalMaximaByZscore <- function(res, dfscoredef) {

 	oi <- order(res$zscore, decreasing=T);
	ozscore <- res$zscore[oi];
	otmin <- round(res$tmin[oi], 2);
	otmax <- round(res$tmax[oi], 2);
	odiscard <- logical(length(ozscore));
	d <- data.frame(oi, otmin, otmax, ozscore, odiscard);
	
	k <- 1;
	previousscoreindex <- -1;
	scorecount <- 0;

	while (k < nrow(d)) {

		d$odiscard[(k+1):nrow(d)] <- columnOverlap(d$otmin[k], d$otmax[k], d$otmin[(k+1):nrow(d)], d$otmax[(k+1):nrow(d)]);
		
		##scoreindex <- getScoreIndexForZscore(d$ozscore[k]);
		##if (previousscoreindex == - 1) {
		##	previousscoreindex <- scoreindex;
		##}
		##if (previousscoreindex != scoreindex) {
		##	cat(paste("Number of class", dfscoredef$score[previousscoreindex], "annotation intervals =", scorecount, "\n"));
		##	scorecount <- 1;
		##       previousscoreindex <- scoreindex;	
		##} else {
		##	scorecount <- scorecount + 1;
		##}

		d <- subset(d, d$odiscard == FALSE);
		k <- k +1;
	}

	##scorecount <- scorecount + 1;	
	##cat(paste("Number of class", dfscoredef$score[previousscoreindex], "intervals =", scorecount, "\n"));		
	
	ti <- order(d$otmin);
	dt  <- d[ti, ];
	tmin <- d$otmin[ti];
	tmax <- d$otmax[ti];
	duration <- round(tmax - tmin, 2);
	wcoef <- round(res$wc[d$oi[ti]], 8);
	zscore <- round(d$ozscore[ti], 8);
	class <- dfscoredef$score[sapply(zscore, getScoreIndexForZscore)];
	quality <- round(res$quality[d$oi[ti]], 8);

	res <- data.frame(tmin, tmax, duration, wcoef, zscore, quality, class);

	return (res);
}

##
## Compute the error amplitude mean and dispersion from the error amplitude values on an interval
## values : The error amplitude values
## itimemin : The index of left boundary in time
## itimemax : The index of right boundary in time
##
errorAmplitudeOnInterval <- function(values, itimemin, itimemax) {

	v <- values[itimemin:itimemax];

	mean <- mean(v, na.rm=T);
	sigma <- sd(v, na.rm=T);

	return (c(mean, sigma));
}

##
## Load the annotation table created from the TextGrid by Praat
## directory : The directory of the file
## filename : The file name
##
loadAnnotationTable <- function(directory, filename) {

	## Read the annotation table
	d <- read.table(paste(directory, filename, sep="/"), h=TRUE, sep="\t");

	if (ncol(d) == 4) {
		dd <- data.frame(d$tmin, d$tmax, d[,3]);
	} else {
		dd <- data.frame(d$tmin, d$tmax, d[,4]);
	}
	colnames(dd)[1] = "tmin";
	colnames(dd)[2] = "tmax";
	colnames(dd)[3] = "ebmotion";
	
	return (dd);
}



##
## From the 2D wavelet coefficient table, extract the local maxima in such a way that there
## is no interval overlap in time
## s2d : The wavelet coefficients table
## sigmamin : The minimal threshold in scale for selection a maximum
##
searchLocalMaxima <- function(s2d, sigmamin) {
 
	dfres <- NULL;

	current_maximum_wcoef <- s2d$wcoef[1];
	current_maximum_time <- s2d$time[1];
	current_maximum_duration <- round(s2d$duration[1], 2);
	current_maximum_quality <- s2d$quality[1];

	for (i in c(1:nrow(s2d))) {
		if (s2d$duration[i] > sigmamin) {
			if (s2d$time[i] - s2d$duration[i]  < current_maximum_time + current_maximum_duration) {
				if (s2d$wcoef[i] > current_maximum_wcoef) {
					current_maximum_time <- s2d$time[i];
					current_maximum_duration <- round(s2d$duration[i], 2);
					current_maximum_wcoef <- s2d$wcoef[i];
					current_maximum_quality <- s2d$quality[i];
				}
			} else {
				if (current_maximum_duration > sigmamin) {
					addinterval <- TRUE;
					if (!is.null(dfres)) {
						for (k in c(1:nrow(dfres))) {
							r <- intervalOverlap(dfres$tmin[k], dfres$tmax[k], (current_maximum_time-current_maximum_duration), (current_maximum_time+current_maximum_duration));
							if (r > 0) {
								if (dfres$wcoef[k] < current_maximum_wcoef) {
									dfres$tmin[k] <- round(current_maximum_time-current_maximum_duration, 2);
									dfres$tmax[k] <- round(current_maximum_time+current_maximum_duration, 2);
								       	dfres$duration[k] <- round(2*current_maximum_duration, 2);	
									dfres$wcoef[k] <- current_maximum_wcoef;
									dfres$quality[k] <- current_maximum_quality;
								}
								addinterval <- FALSE;
							}
						}
					}
					if (addinterval) {
						dff <- data.frame(round(current_maximum_time-current_maximum_duration,2), round(current_maximum_time+current_maximum_duration, 2), round(2*current_maximum_duration, 2), current_maximum_wcoef, current_maximum_quality);
						names(dff) <- c("tmin", "tmax", "duration", "wcoef", "quality");
						rownames(dff) <- NULL;
						dfres <- rbind(dfres, dff);
					}
					##print(paste((current_maximum_time-current_maximum_duration), (current_maximum_time+current_maximum_duration), (2*current_maximum_duration), current_maximum_wcoef, current_maximum_quality));
				}
				current_maximum_time <- s2d$time[i];
				current_maximum_duration <- round(s2d$duration[i], 2);
				current_maximum_wcoef <- s2d$wcoef[i];
				current_maximum_quality <- s2d$quality[i];
			}
		}	
	}
	
	return (dfres);
}

##
## From the 2D wavelet coefficient table, extract the local maxima in such a way that there
## is no interval overlap in time
## s2d : The wavelet coefficients table
## sdwcqmodel : The linear model for decrease of standard dispersion wavelet coefficients with quality 
## sigmamin : The minimal threshold in scale for selection a maximum
## ebmotion : The type of eyebrow motion "Raise" or "Frown"
##
searchLocalMaximaByZscore <- function(s2d, sdwcqmodel, sigmamin, ebmotion) {

      	s2d$zscore <- numeric(nrow(s2d));

	for (i in c(1:nrow(s2d))) {
		sq <- sdwcqmodel$intercept[1] + sdwcqmodel$slope[1]*getQualityCode(s2d$quality[i]);
		s2d$zscore[i] <- s2d$wcoef[i]/sq;
	}	
	
	dfres <- NULL;

	maximaltime <- round(s2d$time[nrow(s2d)], 2);
	current_maximum_zscore <- s2d$zscore[1];
	current_maximum_wcoef <- s2d$wcoef[1];
	current_maximum_time <- round(s2d$time[1], 2);
	current_maximum_duration <- round(s2d$duration[1], 2);
	current_maximum_quality <- s2d$quality[1];

	for (i in c(1:nrow(s2d))) {
		if (round(s2d$duration[i], 2) > sigmamin) {
			if (round(s2d$time[i] - s2d$duration[i], 2)  < current_maximum_time + current_maximum_duration) {
				if (s2d$zscore[i] > current_maximum_zscore) {
					current_maximum_time <- round(s2d$time[i], 2);
					current_maximum_duration <- round(s2d$duration[i], 2);
					current_maximum_wcoef <- s2d$wcoef[i];
					current_maximum_zscore <- s2d$zscore[i];
					current_maximum_quality <- s2d$quality[i];
				}
			} else {
				if (current_maximum_duration > sigmamin) {
					addinterval <- TRUE;
					if (!is.null(dfres)) {
						for (k in c(1:nrow(dfres))) {
							r <- intervalOverlap(dfres$tmin[k], dfres$tmax[k], (current_maximum_time-current_maximum_duration), (current_maximum_time+current_maximum_duration));
							if (r > 0) {
								if (dfres$zscore[k] < current_maximum_zscore) {
									amin <- max(0, round(current_maximum_time-current_maximum_duration, 2));
									amax <- min(maximaltime, round(current_maximum_time+current_maximum_duration, 2)); 
									dfres$tmin[k] <- amin;
									dfres$tmax[k] <- amax;
								       	dfres$duration[k] <- round(amax-amin, 2);
									dfres$wcoef[k] <- current_maximum_wcoef;
									dfres$zscore[k] <- current_maximum_zscore;
									dfres$quality[k] <- current_maximum_quality;
								}
								addinterval <- FALSE;
							}
						}
					}
					if (addinterval) {
						amin <- max(0, round(current_maximum_time-current_maximum_duration, 2));
						amax <- min(maximaltime, round(current_maximum_time+current_maximum_duration, 2));
						dff <- data.frame(amin, amax, round(amax-amin, 2), current_maximum_wcoef, current_maximum_zscore, current_maximum_quality, ebmotion);
						names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "ebmotion");
						rownames(dff) <- NULL;
						dfres <- rbind(dfres, dff);
					}
					##print(paste((current_maximum_time-current_maximum_duration), (current_maximum_time+current_maximum_duration), (2*current_maximum_duration), current_maximum_wcoef, current_maximum_quality));
				}
				current_maximum_time <- round(s2d$time[i], 2);
				current_maximum_duration <- round(s2d$duration[i], 2);
				current_maximum_zscore <- s2d$zscore[i];
				current_maximum_wcoef <- s2d$wcoef[i];
				current_maximum_quality <- s2d$quality[i];
			}
		}	
	}

	## Remove multiple lines
	rmindex <- numeric(0);
	for (i in c(2:nrow(dfres))) {
		if (dfres$tmin[i] == dfres$tmin[i-1]) {
			rmindex <- c(rmindex, i);
		}
	}
	if (length(rmindex) != 0) {
		dfres <- dfres[-rmindex, ];
	}
	return (dfres);
}

overlap <- function(tmin1, tmax1, tmin2, tmax2) {
	return (!(tmax2 <= tmin1 | tmax1 <= tmin2));
}

columnOverlap <- function(tmin1, tmax1, ctmin2, ctmax2) {
	return (ifelse(ctmax2 <= tmin1 | tmax1 <= ctmin2, FALSE, TRUE));
}

##
## From the 2D wavelet coefficient table, extract the local maxima in such a way that there
## is no interval overlap in time
## s2d : The wavelet coefficients table
## sdwcqmodel : The linear model for decrease of standard dispersion wavelet coefficients with quality 
## sigmamin : The minimal threshold in scale for selection a maximum
## ebmotion : The type of eyebrow motion "Raise" or "Frown"
##
searchLocalMaximaByZscore1 <- function(s2d, sdwcqmodel, sigmamin, ebmotion) {

	fsmin <- ifelse(s2d$duration < sigmamin, FALSE, TRUE);
	
	s2df <- subset(s2d, fsmin == TRUE);

	##s2df$zscore <- s2df$wcoef/(sdwcqmodel$intercept[1] + sdwcqmodel$slope[1]*getQualityCode(s2df$quality));
	
      	s2df$zscore <- numeric(nrow(s2df));
	for (i in c(1:nrow(s2df))) {
		sq <- (sdwcqmodel$intercept[1] + sdwcqmodel$slope[1]*getQualityCode(s2df$quality[i]));
		s2df$zscore[i] <- s2df$wcoef[i]/sq;
	}	

 	oi <- order(s2d$zscore, decreasing=T);
	otmin <- s2d$time[oi] - s2d$duration[oi]/2
	otmax <- s2d$time[oi] + s2d$duration[oi]/2
	ozscore <- s2d$zscore[oi]
	odiscard <- logical(length(ozscore))
	d <- data.frame(oi, otmin, otmax, ozscore, odiscard);
	
	k <- 1
	while (k < nrow(d)) {
		for (i in c((k+1):nrow(d))) {
			d$odiscard[i] <- overlap(d$otmin[k], d$otmax[k], d$otmin[i], d$otmax[i]);
		}
		print(paste(k, nrow(d)));
		d <- subset(d, d$odiscard == FALSE);
		k <- k +1;
	}
	
	k <- 1
	n <- nrow(d); 
	while (k < n) {
		d$odiscard[(k+1):nrow(d)] <- columnOverlap(d$otmin[k], d$otmax[k], d$otmin[(k+1):n], d$otmax[(k+1):n]);
		print(paste(k, nrow(d)));
		d <- subset(d, d$odiscard == FALSE);
		n <- nrow(d); 
		k <- k +1;
	}
	ti <- order(d$otmin);
	dt  <- d[ti, ];

	maximaltime <- round(s2d$time[nrow(s2d)], 2);
	current_maximum_zscore <- s2d$zscore[1];
	current_maximum_wcoef <- s2d$wcoef[1];
	current_maximum_time <- round(s2d$time[1], 2);
	current_maximum_duration <- round(s2d$duration[1], 2);
	current_maximum_quality <- s2d$quality[1];

	for (i in c(1:nrow(s2d))) {
		if (round(s2d$duration[i], 2) > sigmamin) {
			if (round(s2d$time[i] - s2d$duration[i], 2)  < current_maximum_time + current_maximum_duration) {
				if (s2d$zscore[i] > current_maximum_zscore) {
					current_maximum_time <- round(s2d$time[i], 2);
					current_maximum_duration <- round(s2d$duration[i], 2);
					current_maximum_wcoef <- s2d$wcoef[i];
					current_maximum_zscore <- s2d$zscore[i];
					current_maximum_quality <- s2d$quality[i];
				}
			} else {
				if (current_maximum_duration > sigmamin) {
					addinterval <- TRUE;
					if (!is.null(dfres)) {
						for (k in c(1:nrow(dfres))) {
							r <- intervalOverlap(dfres$tmin[k], dfres$tmax[k], (current_maximum_time-current_maximum_duration), (current_maximum_time+current_maximum_duration));
							if (r > 0) {
								if (dfres$zscore[k] < current_maximum_zscore) {
									amin <- max(0, round(current_maximum_time-current_maximum_duration, 2));
									amax <- min(maximaltime, round(current_maximum_time+current_maximum_duration, 2)); 
									dfres$tmin[k] <- amin;
									dfres$tmax[k] <- amax;
								       	dfres$duration[k] <- round(amax-amin, 2);
									dfres$wcoef[k] <- current_maximum_wcoef;
									dfres$zscore[k] <- current_maximum_zscore;
									dfres$quality[k] <- current_maximum_quality;
								}
								addinterval <- FALSE;
							}
						}
					}
					if (addinterval) {
						amin <- max(0, round(current_maximum_time-current_maximum_duration, 2));
						amax <- min(maximaltime, round(current_maximum_time+current_maximum_duration, 2));
						dff <- data.frame(amin, amax, round(amax-amin, 2), current_maximum_wcoef, current_maximum_zscore, current_maximum_quality, ebmotion);
						names(dff) <- c("tmin", "tmax", "duration", "wcoef", "zscore", "quality", "ebmotion");
						rownames(dff) <- NULL;
						dfres <- rbind(dfres, dff);
					}
					##print(paste((current_maximum_time-current_maximum_duration), (current_maximum_time+current_maximum_duration), (2*current_maximum_duration), current_maximum_wcoef, current_maximum_quality));
				}
				current_maximum_time <- round(s2d$time[i], 2);
				current_maximum_duration <- round(s2d$duration[i], 2);
				current_maximum_zscore <- s2d$zscore[i];
				current_maximum_wcoef <- s2d$wcoef[i];
				current_maximum_quality <- s2d$quality[i];
			}
		}	
	}

	## Remove multiple lines
	rmindex <- numeric(0);
	for (i in c(2:nrow(dfres))) {
		if (dfres$tmin[i] == dfres$tmin[i-1]) {
			rmindex <- c(rmindex, i);
		}
	}
	if (length(rmindex) != 0) {
		dfres <- dfres[-rmindex, ];
	}
	return (dfres);
}

##
## Merge two data frames containing list of local maxima intervals by selecting the higher value when intervals overlap
## lm1 : The first data frame containing list of local maxima intervals
## lm2 : The first data frame containing list of local maxima intervals
##
mergeLocalMaximaByZscore <- function(lm1, lm2) {

      	lm <- rbind(lm1, lm2);
	lm <- lm[order(lm$tmin), ];

	dfres <- NULL;

	dff <- lm[1, ];
	
	dfres <- rbind(dfres, lm[1, ]);

	for (i in c(2:nrow(lm))) {
		
		k <- nrow(dfres);
		r <- intervalOverlap(dfres$tmin[k], dfres$tmax[k], lm$tmin[i], lm$tmax[i]);
		if (r > 0) {
			if (dfres$zscore[k] < lm$zscore[i]) {
				dfres[k, ] <- lm[i, ];
			}
		} else {
			dfres <- rbind(dfres, lm[i, ]);
		}
	}

	rownames(dfres) <- NULL;	

	return (dfres);
}

isLocalMaximum <- function(i, pc, cc, nc, threshold) {

	if (cc[i] <= threshold) return (FALSE);

	if (cc[i] < pc[i+1]) return (FALSE);
	if (cc[i] < cc[i+1]) return (FALSE);
	if (cc[i] < nc[i+1]) return (FALSE);

	if (cc[i] < pc[i]) return (FALSE);
	if (cc[i] < nc[i]) return (FALSE);

	if (cc[i] < pc[i-1]) return (FALSE);
	if (cc[i] < cc[i-1]) return (FALSE);
	if (cc[i] < nc[i-1]) return (FALSE);

	return (TRUE);
}

createXYreferenceByFrame <- function(data, i_frame) {

	d <- subset(data, data$frame == i_frame);

	## Create the 4 vectors
	index <- numeric(0);
	point <- character(0);
	X <- numeric(0);
	Y <- numeric(0);

	point_length = 12;

	for (ipoint in c(1:point_length)) {
	     ## Add for each vector its value for the line ipoint
	     index <- c(index, ipoint);
	     point <- c(point, pointName(ipoint));
	     X <- c(X, d$value[2*ipoint-1]);
	     Y <- c(Y, d$value[2*ipoint]);
	}

	## Add the end create the new data.frame
	dd <- data.frame(index, point, X, Y);

	return (dd);
}


plotResiduals <- function(tref, d, frame_index) {

	tref$Xp  = numeric(12);
	tref$Yp  = numeric(12);
	tref$Zp  = numeric(12);

	for (j in c(1:12)) {
		
		if (j != 5) {
			column_index = which(colnames(d) == paste(pointName(j), "rx", sep=""));
			tref$Xp[j] <- d[frame_index, column_index];
			column_index = which(colnames(d) == paste(pointName(j), "ry", sep=""));
			tref$Yp[j] <-  d[frame_index, column_index];
			column_index = which(colnames(d) == paste(pointName(j), "rz", sep=""));
			tref$Zp[j] <-  d[frame_index, column_index];
		}
	}
	tref$Xp <- tref$Xp + tref$X;
	tref$Yp <- tref$Yp + tref$Y;
	tref$Zp <- tref$Zp + tref$Z;

	qp <- qplot(data=tref, x=X, y=Y, label=index) + geom_point(size=I(1), shape=I(19), colour="blue") + scale_y_reverse(limits = c(45, -15)) + coord_fixed(ratio = 1, expand = TRUE);

	qp <- qp + scale_x_continuous(limits = c(-35, 35));

	qp <- qp + geom_segment(aes(x = X, y = Y, xend = Xp, yend = Yp, colour = "blue"), size=1, arrow = arrow(length = unit(0.1,"cm")), data = tref) + guides(colour=FALSE);

	qp <- qp +  annotate("text", x = -30, y = 22, label = paste("t =", d[frame_index, 2], "s"));
	qp <- qp +  annotate("text", x = -30, y = 26, label = paste("frame :", d[frame_index, 1]));
	##qp <- qp + geom_point(data=tref, aes(x=Xp, y=Yp), size=I(2), shape=I(19), colour="red");
	return (qp);


}

computeNotEnoughDataArea <- function(data, maximal_ratio_of_na_point, window_size) {

	n_frame = nrow(data);
	data$label <- as.character(data$label);

	for (i in c(1:(n_frame-window_size))) {

		na_count = 0;
		for (j in c(i:(i+window_size-1))) {
		    
		    if (data$label[j] == "N" | data$label[j] == "A" | data$label[j] == "C" | data$label[j] == "S") {
			    na_count = na_count+1;
		    }
		}
		if (na_count > (maximal_ratio_of_na_point*window_size)) {
			data$label[i] <- "R";
		} else {
			if (data$label[i] != "N" & data$label[i] != "A" & data$label[i] != "C" & data$label[i] != "S") {
				data$label[i] <- "Y";
			}
		}		
	}

	for (i in c((n_frame-window_size+1):n_frame)) {
		data$label[i] <- data$label[n_frame-window_size];
	}

	data$label <- as.factor(data$label);

	return (data);
}


pointName <- function(ipoint) {

	if (ipoint < 10) {
		point_name = paste("P_0", as.character(ipoint), sep="");
	} else {
		point_name = paste("P_", as.character(ipoint), sep="");
	}
	return (point_name);
}

dimensionName <- function(idimension) {

	if (idimension == 1) {
		dimension_name = "X";
	} else {
		dimension_name = "Y";
	}
	return (dimension_name);
}

hausse <- function(d) { 
	SHx <-  0.5*d$P_02rx - 0.5*d$P_03rx;
	weightx = 1;
	SHy <- -0.5*d$P_01ry - 0.5*d$P_04ry -0.5*d$P_02ry - 0.5*d$P_03ry  + 0.25*d$P_08ry + 0.25*d$P_09ry;
	weighty = 2.5;
	weight = weightx + weighty;
	SH <- SHx + SHy;
	SH <- SH/weight;
 	 return (SH);	
}

createTableOfAnnotatedIntervals <- function(d, label) {

	l = which(d$ebmotion == label);

	index_min = numeric(0);
	index_max = numeric(0);
	ebmotion = character(0);
	

	index_min <- c(index_min, l[1]);
	ebmotion <- c(ebmotion, label);
	for (i in c(2:(length(l)-1))) {
		if (l[i] != l[i-1] + 1) {
			index_max <- c(index_max, l[i-1]);
			index_min <- c(index_min, l[i]);
			ebmotion <- c(ebmotion, label);
		} 
	}
	index_max <- c(index_max, l[length(l)]);
	
 		
	dd <- data.frame(ebmotion, index_min, index_max);

	quality = numeric(0);
	for (i in c(1:nrow(dd))) {
		if ((d$time[index_max[i]] -  d$time[index_min[i]]) == 0) {
			quality = c(quality, 0);
		} else {
			quality = c(quality, 1/((d$time[index_max[i]] -  d$time[index_min[i]])/(index_max[i] -  index_min[i]))/25);
		}
	}

	dd <- data.frame(dd, quality);
	return(dd);
}


createXYmotion <- function(data, lmer) {
	
	## The random effect
	coefrlmer <- ranef(lmer, condVar = TRUE);
	coefnames <- row.names(coefrlmer[[1]]);
	## Get the values and reorder in function of the frame index 
	coefvaluesX <- coefrlmer[[1]]$`(Intercept)`;
	coefvaluesY <- coefrlmer[[1]]$`dimensionY`;

	d <- subset(data, data$point == "P_01" & data$dimension == "X");

	frame <- numeric(nrow(d));
	time <- numeric(nrow(d));
	mX <- numeric(nrow(d));
	mY <- numeric(nrow(d));
	
	for (i in c(1:nrow(d))) {
		mX[i] <- coefvaluesX[which(coefnames == d$frame[i])];
		mY[i] <- coefvaluesY[which(coefnames == d$frame[i])];
		frame[i] <- d$frame[i];
		time[i] <- d$time[i];
	}

	dd <- data.frame(frame, time, mX, mY);
	return (dd);
}

createTableResTimeForAnnotation <- function(d, table, index_annotation, margin) {

	index_min = max(1, (table$index_min[index_annotation] - margin));
	index_max = min(table$index_max[nrow(table)], (table$index_max[index_annotation] + margin));

	dd <- d[index_min:index_max,	];
		

	return(dd);
}

motionAngleMeasurement <- function(d, time_margin) {

	motionangle <- numeric(nrow(d));

	for (i in c(1:nrow(d))) {
		
	   	if (!is.na(d$yaw[i]) & !is.na(d$yaw[i+1])) {
			motionangle[i] = (d$yaw[i+1] - d$yaw[i])^2 + (d$roll[i+1] - d$roll[i])^2 + (d$pitch[i+1] - d$pitch[i])^2;
		} else {
			motionangle[i] <- NA;
		}

	}
	return (motionangle);
}

addMotionAngleMeasurement <- function(d) {

	d$rollprime <- NA;
	d$yawprime <- NA;
	d$pitchprime <- NA;

	for (i in c(1:(nrow(d)-1))) {
		
	   	if (!is.na(d$yaw[i]) & !is.na(d$yaw[i+1])) {
			d$rollprime[i] <- (d$roll[i+1] - d$roll[i])/0.04;
			d$yawprime[i] <- (d$yaw[i+1] - d$yaw[i])/0.04;
			d$pitchprime[i] <- (d$pitch[i+1] - d$pitch[i])/0.04;
		}

	}
	return (d);
}
