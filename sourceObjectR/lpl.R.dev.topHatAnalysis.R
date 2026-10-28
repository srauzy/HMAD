##
## lpl.R.dev.topHatAnalysis.R, creator S. Rauzy, LPL, 12/12/2022
## Various functions for performing a Top Hat wavelet analysis
##

##
## Create and return the data frame containing the wavelet coefficients for the defined time scales space
## times : The time column
## values : The column with the values of the function
## scale_grid_parameters : The data frame containing the scale grid parameters
## verbose : Put this flag at TRUE for printing comments in the R console
## return the data frame containing the wavelet coefficients for the defined time scales space
##
lpl.R.dev.topHatAnalysis.computeTopHatWaveletTransform <- function(times, values, scale_grid_parameters, verbose) {

	sigmamin <- scale_grid_parameters$smin;
	sigmamax <- scale_grid_parameters$smax;
	sigmastep <- scale_grid_parameters$step;

	ns <- (sigmamax-sigmamin)/sigmastep + 1;
	nt <- length(times);
	time <- numeric(nt*ns);
	duration <- numeric(nt*ns);
	tmin <- numeric(nt*ns);
	tmax <- numeric(nt*ns);
	wgof <- numeric(nt*ns);
	wcoef <- numeric(nt*ns);

	if (verbose) cat("Wavelet analysis on scale from", sigmamin, "to", sigmamax, "seconds: ");

	for (i in (1:ns)) {

		timescale <- sigmamin + (i-1)*sigmastep;
		if (verbose) cat(timescale, "  ");

		ftransform <- lpl.R.dev.topHatAnalysis.computeArrayOfWaveletCoefficients(times, values, timescale); 
		for (j in (1:nt)) {
			index <- round(ns*(j-1) + i, 0);
			duration[index] <- timescale;
			time[index] <- times[j];
			tmin[index] <- times[j] - timescale/2;
			tmax[index] <- times[j] + timescale/2;
			wgof[index] <- ftransform[[1]][j];
			wcoef[index] <- ftransform[[2]][j];
		}
		##cat("Scale", i, "completed\n");
	}
	if (verbose) cat("\n");

	## Round the values
	time <- round(time, 2);
	tmin <- round(tmin, 2);
	tmax <- round(tmax, 2);
	duration <- round(duration, 2);
	wgof <- round(wgof, 4);
	wcoef <- round(wcoef, 4);

	## Return the data frame
	return (data.frame(time, duration, tmin, tmax, wgof, wcoef));
}

##
## Compute the Top Hat wavelet coefficient for a given time and a given time scale
## times : The time column
## values : The column with the values of the function
## range : The range of indexes covered by the wavelet at this time scale
## index : The time position index in the times and values column
## a : The time scale in time unit
## gamma : The constant of normalisation gamma (-1 or -1/2)
## return the wavelet coefficient and goodness-of-fit value as a 2 dimensionnal vector
##
lpl.R.dev.topHatAnalysis.computeTopHatWaveletCoefficient <- function(times, values, range, index, a, gamma) {

	## The size of the time array
	n <- length(times);
	## index min and index max (accounting for the boundaries of the arrays)
	imin <- max(1, index - 3*range);
	imax <- min(n, index + 3*range);
	## The effective range of the interval, 2*range if out of the boudaries
	effective_range <- imax - imin;

	sf2 <- 0;
	sw2 <- 0;
	sfw <- 0;

	## Loop on the indexes which matter
	for (i in c(imin:imax)) {
		
		f <- values[i];

		if ((i-imin)/(imax-imin) < 0) {
			w <- 0;
		} else if ((i-imin)/(imax-imin) < 1/3) { 
			w <- -1/sqrt(6);
		} else if ((i-imin)/(imax-imin) < 2/3) {
			w <- 2/sqrt(6);
		} else if ((i-imin)/(imax-imin) < 1) {
			w <- -1/sqrt(6);
		} else {
			w <- 0;
		}	
		## Norms and scalar product
		sf2 <- sf2 + f*f;
		sw2 <- sw2 + w*w;
		sfw <- sfw + f*w;
	}

	## Normalize
	if (sf2 == 0) {
		gof <- 0;
	} else {
		gof <- sfw/sqrt(sf2*sw2);
	}
	## With this normalization a spike sinus function of amplitude A presents a maximum in 2D space of amplitude A
	coefficient <- sfw/effective_range*a^(1+gamma);

	## The result as a couple goodness-of-fit and coefficent
	result <- c(gof, coefficient);

	return (result);
}

lpl.R.dev.topHatAnalysis.topHatWavelet <- function(times, position, scale) {

	timestep <- times[2] - times[1];
	imin <- ((position - 3*scale/2) - times[1])/timestep;
	imax <- ((position + 3*scale/2) - times[1])/timestep;

	f <- numeric(length(times));

	for (i in c(1:length(times))) { 
		if ((i-imin)/(imax-imin) < 0) {
		       f[i] <- 0;	
		} else if ((i-imin)/(imax-imin) < 1/3) { 
			f[i] <- -1/sqrt(6);
		} else if ((i-imin)/(imax-imin) < 2/3) {
			f[i] <- 2/sqrt(6);
		} else if ((i-imin)/(imax-imin) < 1) {
			f[i] <- -1/sqrt(6);
		} else {
			f[i] <- 0;
		}	
	}

	return (f);
}

##
## In the Sinus Spike projection, compute the goodness of fit and the amplitude of the values for a given sinus spike duration
## times The times column
## values : The values of the signal at times
## down : Boolean at TRUE if the basis sinus function starts to go down before going up
## duration : The duration 
## Return a data frame with the wgof and wcoef
##
lpl.R.dev.topHatAnalysis.computeArrayOfWaveletCoefficients <- function(times, values, duration) {

	## Compute the time step from the times array
	timestep <- round(times[2] - times[1], 2);

	## Compute the range (based on the timestep)
	range <- round(duration/2/timestep, 0);

	## The size of the time array
	n <- length(times);
	## The array of predictions in output
	predictions_gof <- numeric(n);
	predictions_coefficient <- numeric(n);

 	for (i in c(1:n)) {
		result <- lpl.R.dev.topHatAnalysis.computeTopHatWaveletCoefficient(times, values, range, i, duration, -1/2);
		predictions_gof[i] <- result[1];
		predictions_coefficient[i] <- result[2];
	}
	predictions <- list(predictions_gof, predictions_coefficient);

	return (predictions);
}

## 
## Store the scale grid parameters as a list
## smin : The minimal time scale
## smax : The maximal time scale
## step : The scale step
## return the scale grid parameters as a list
##
lpl.R.dev.topHatAnalysis.createScaleGridParameters <- function(smin, smax, step) {

	object_result <- vector(mode = "list", length = 3);

	object_result[[1]] <- smin;
	object_result[[2]] <- smax;
	object_result[[3]] <- step;

	names(object_result) <- c("smin", "smax", "step");	
	return (object_result);
}
