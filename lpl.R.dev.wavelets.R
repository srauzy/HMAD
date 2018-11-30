## lpl.R.dev.wavelets.R, creator S. Rauzy, LPL, 15/11/2016

##
## The mother wavelet of integral 0 and scale sigma (the normalisation is in sigma^(-1)
## time : The time on which the wavelet in centered
## sigma : The scale of the wavelet
## wDparameter : For a standard top hat wavelet, D = 2*sigma, but one can tune this parameter to change the area spanned by the negative
## left and right side of the mother wavelet
##
lpl.R.dev.wavelets.motherwavelet <- function(time, sigma, wDparameter) {

	if (time < (-1-wDparameter)*sigma/2) {
		return (0);
	} else if (time < (-sigma/2)) {
		##return ((-sqrt(sigma)/D));
		return ((-1/(wDparameter*sigma)));
		##return ((-sigma/D/sqrt(sigma)));
	} else if (time < sigma/2) {
		##return (sqrt(sigma)/sigma);
		return (1/sigma);
		##return (1/sqrt(sigma));
	} else if (time < (1+wDparameter)*sigma/2) {
		##return ((-sigma/D/sqrt(sigma)));
		##return ((-sqrt(sigma)/D));
		return ((-1/(wDparameter*sigma)));
	} else {
		return (0);
	}
}

##
## The mother wavelet of integral 0 and scale sigma (the normalisation is in sigma^(-1)
## time : The time on which the wavelet in centered
## sigma : The scale of the wavelet
## D : For a standard top hat wavelet, D = 2*sigma, but one can tune this parameter to change the area spanned by the negative
## left and right side of the mother wavelet
##
lpl.R.dev.wavelets.wavelet <- function(si, siD) {

	n <- 2*si + 2*siD;
	w <- numeric(n);
	for (i in c(1:n)) {
		if (i < siD+1) {
			w[i] <- -1/siD;
		} else if (i < siD+2*si+1) {
			w[i] <- 1/si;
		} else {
			w[i] <- -1/siD;
		}
	}
	return (w);
}

##
## The mother wavelet of integral 0 and scale sigma (the normalisation is in sigma^(-1)
## time : The time on which the wavelet in centered
## sigma : The scale of the wavelet
## D : For a standard top hat wavelet, D = 2*sigma, but one can tune this parameter to change the area spanned by the negative
## left and right side of the mother wavelet
##


lpl.R.dev.wavelets.iwaveletatscale <- function(d, si, wDparameter) {

	siD <- wDparameter*si;
	n <- nrow(f);
	##wcvector <- vector(mode="list", length=n);
	wcvector <-numeric(n);

	for (i in c(1:n)) {
		if (i-si-siD+1 < 1 | i+si+siD > n) {
			##wcvector[[i]] <- c(0, -1, -1, -1);
			wcvector[i] <- -1;
		} else {	
			imin = max(1, i-si-siD+1);
			imax <- min(n, i+si+siD);

			dl <- d[imin:(i-si), ]; 
			dc <- d[(i-si+1):(i+si), ]; 
			dr <- d[(i+si+1):imax, ];

			wcvector[i] <- lpl.R.dev.wavelets.iwavelet(dl, dc, dr, si, siD)[1];
			##wcvector[[i]] <- lpl.R.dev.wavelets.iwavelet(dl, dc, dr, si, siD);
		}
	}
	return (wcvector);
}

lpl.R.dev.wavelets.iwaveletatscale1 <- function(d, si, wDparameter) {

	siD <- wDparameter*si;
	n <- nrow(f);
	##wcvector <- vector(mode="list", length=n);
	wcvector <-numeric(n);

	for (i in c(1:n)) {
		if (i-si-siD+1 < 1 | i+si+siD > n) {
			##wcvector[[i]] <- c(0, -1, -1, -1);
			wcvector[i] <- -1;
		} else {	
			imin = max(1, i-si-siD+1);
			imax <- min(n, i+si+siD);

			dl <- d[imin:(i-si), ]; 
			dc <- d[(i-si+1):(i+si), ]; 
			dr <- d[(i+si+1):imax, ];

			wcvector[i] <- lpl.R.dev.wavelets.iwavelet(dl, dc, dr, si, siD)[1];
			##wcvector[[i]] <- lpl.R.dev.wavelets.iwavelet(dl, dc, dr, si, siD);
		}
	}
	return (wcvector);
}

lpl.R.dev.wavelets.iwaveletatscale <- function(d, si, wDparameter) {

	siD <- wDparameter*si;
	n <- nrow(d);
	##wcvector <- vector(mode="list", length=n);
	wc <-numeric(n);
	q <-numeric(n);

	for (i in c(1:n)) {
		if (i-si-siD+1 < 1 | i+si+siD > n) {
			##wcvector[[i]] <- c(0, -1, -1, -1);
			wc[i] <- 0;
			q[i] <- 0;
		} else {	
			imin = max(1, i-si-siD+1);
			imax <- min(n, i+si+siD);

			dl <- d[(i-si-siD+1):(i-si), ]; 
			dc <- d[(i-si+1):(i+si), ]; 
			dr <- d[(i+si+1):(i+si+siD), ];

			res <- lpl.R.dev.wavelets.iwavelet(dl, dc, dr, si, siD);
			wc[i] <- res[1];
			q[i] <- res[2];
			##wcvector[[i]] <- lpl.R.dev.wavelets.iwavelet(dl, dc, dr, si, siD);
		}
	}
	return (data.frame(wc, q));
}

## si= 1; siD=si; 
## si= 1; siD=si; dl <- d[(i-si-siD+1):(i-si), ]; dc <- d[(i-si+1):(i+si), ]; dr <- d[(i+si+1):(i+si+siD), ]
lpl.R.dev.wavelets.iwavelet <- function(dl, dc, dr, si, siD) {

	qwl <- sum(dl$w)/nrow(dl);
	qwc <- sum(dc$w)/nrow(dc);
	qwr <- sum(dr$w)/nrow(dr);

	wc <- 0;
	if (qwl != 0) wc <- wc - sum(dl$f*dl$w)/qwl/siD;
        if (qwc != 0) wc <- wc + sum(dc$f*dc$w)/qwc/si;
        if (qwr != 0) wc <- wc - sum(dr$f*dr$w)/qwr/siD; 

	return (c(wc, (qwl+qwc+qwr)/3));
}	


lpl.R.dev.wavelets.wavelet <- function(si, siD) {

	n <- 2*si + 2*siD;
	w <- numeric(n);
	for (i in c(1:n)) {
		if (i < siD+1) {
			w[i] <- -1/siD;
		} else if (i < siD+2*si+1) {
			w[i] <- 1/si;
		} else {
			w[i] <- -1/siD;
		}
	}
	return (w);
}


topHatMW <- function(time, sigma) {

	if (time < (-3*sigma)) {
		return (0);
	} else if (time < (-sigma)) {
		return (-sqrt(2*sigma));
	} else if (time < sigma) {
		return (2*sqrt(2*sigma));
	} else if (time < (3*sigma)) {
		return (-sqrt(2*sigma));
	} else {
		return (0);
	}
}

waveletcoef <- function(times, values, index, sigma) {

	return (waveletcoefd(times, values, index, sigma, sigma*2));
}

searchClosestIndex <- function(times, time) {

	i = 1;
	while (times[i] < time & i < length(times)) {
		i = i+1;
	}
	return (i);
}

##
## Compute the wavelet coefficients at time indexed by itime in the column times for
## a scale sigma
## times : The column of times
## values : The function to transform (on the times support)
## sigma : The scale 
## itime : The time index in the columns time
## itimemin : The index corresponding to the left boundary of the wavelet to convolve with
## itimemax : The index corresponding to the right boundary of the wavelet to convolve with
##
lpl.R.dev.wavelets.waveletCoefAtTime <- function(times, values, sigma, itime, itimemin, itimemax, wDparameter) {

	coef = 0;
	countL = 0;
	countC = 0;
	countR = 0;
	countLobs = 0;
	countCobs = 0;
	countRobs = 0;

	for (i in c(itimemin:itimemax)) {
	   	
		if ((times[i]-(times[itime]+0.02)) < (-sigma/2)) {
			if (!is.na(values[i])) {
				countLobs = countLobs + 1;
			}
			countL = countL + 1;
		} else if ((times[i]-(times[itime]+0.02)) < (sigma/2)) {
			if (!is.na(values[i])) {
				countCobs = countCobs + 1;
			}
			countC = countC + 1;
		} else {
			if (!is.na(values[i])) {
				countRobs = countRobs + 1;
			}
			countR = countR + 1;
		}
	}

	for (i in c(itimemin:itimemax)) {
	   	if (!is.na(values[i])) {
			if ((times[i]-(times[itime]+0.02)) < (-sigma/2)) {
				coef = coef + countL/countLobs*values[i]*lpl.R.dev.wavelets.motherwavelet(times[i]-times[itime], sigma, wDparameter);
			} else if ((times[i]-(times[itime]+0.02)) < (sigma/2)) {
				coef = coef + countC/countCobs*values[i]*lpl.R.dev.wavelets.motherwavelet(times[i]-times[itime], sigma, wDparameter);
			} else {
				coef = coef + countR/countRobs*values[i]*lpl.R.dev.wavelets.motherwavelet(times[i]-times[itime], sigma, wDparameter);
			}
		}
	}
	if (countL == 0) {
		qL = 0;
	} else {
		qL = countLobs/countL;
	}
	if (countC == 0) {
		qC = 0;
	} else {
		qC = countCobs/countC;
	}
	if (countR == 0) {
		qR = 0;
	} else {
		qR = countRobs/countR;
	}
	coef <- coef*0.04;
	return (c(coef, qL, qC, qR));
}

waveletCoefDAtTime <- function(times, values, time, sigma, D) {

	index = searchClosestIndex(times, time);

	imin = max(1, (index - (sigma+D)*25));
	imax = min(length(times), (index + (sigma+D)*25));
	coef = 0;
	count = 0;
	for (i in c(imin:imax)) {
	   	if (!is.na(values[i])) {
			coef = coef + values[i]*lpl.R.dev.wavelets.motherwavelet(times[i]-times[index], sigma, D);
			count = count + 1;
		}
	}
	return (c(coef, count));
}

waveletcoefd <- function(times, values, index, sigma, D) {

	time = times[index];
	imin = max(1, (index - (sigma+D)*25));
	imax = min(length(times), (index + (sigma+D)*25));
	coef = 0;
	count = 0;
	for (i in c(imin:imax)) {
		
	   	if (!is.na(values[i])) {
			coef = coef + values[i]*lpl.R.dev.wavelets.motherwavelet(times[i]-times[index], sigma, D);
			count = count + 1;
		}
	}
	return (c(coef, count));
}

waveletcoefsd <- function(times, values, sigma, D) {

	coefs <- numeric(length(times));
	for (i in c(1:length(times))) {
		coefs[i] <- waveletcoef(times, values, i, sigma, D)[1];
	}
	return (coefs);
}

waveletcoefs <- function(times, values, sigma) {

	coefs <- numeric(length(times));
	for (i in c(1:length(times))) {
		coefs[i] <- waveletcoef(times, values, i, sigma)[1];
	}
	return (coefs);
}

waveletcoefsAtTime <- function(times, values, time, sigmamin, sigmamax) {

	n_step <- floor((sigmamax-sigmamin)/0.04) + 1;
	coefs <- numeric(n_step);
	for (i in c(1:n_step)) {
		sigma <- sigmamin + (i-1)*0.04;
		coefs[i] <- waveletCoefDAtTime(times, values, time, sigma, 2*sigma)[1];
	}
	return (coefs);
}

waveletcoefs2D <- function(times, values, timemin, timemax, sigmamin, sigmamax) {

	t_step <- floor((timemax-timemin)/0.04) + 1;
	d_step <- floor((sigmamax-sigmamin)/0.04) + 1;
	wcoef <- numeric((t_step*d_step));
	location <- numeric((t_step*d_step));
	duration <- numeric((t_step*d_step));
	for (i in c(1:t_step)) {
		time <- timemin + (i-1)*0.04;
		print(time);
		for (j in c(1:d_step)) {
			sigma <- sigmamin + (j-1)*0.04;

			location[(i-1)*d_step+j] = time;
			duration[(i-1)*d_step+j] = sigma;
			wcoef[(i-1)*d_step+j] <- waveletCoefDAtTime(times, values, time, sigma, 2*sigma)[1];
		}
	}
	dd <- data.frame(location, duration, wcoef);
	return (dd);
}

waveletcoefs2Dts <- function(times, values, timemin, timemax, sigmamin, sigmamax, time_step) {

	t_step <- floor((timemax-timemin)/time_step) + 1;
	d_step <- floor((sigmamax-sigmamin)/time_step) + 1;
	wcoef <- numeric((t_step*d_step));
	location <- numeric((t_step*d_step));
	duration <- numeric((t_step*d_step));
	for (i in c(1:t_step)) {
		time <- timemin + (i-1)*time_step;
		print(time);
		for (j in c(1:d_step)) {
			sigma <- sigmamin + (j-1)*time_step;

			location[(i-1)*d_step+j] = time;
			duration[(i-1)*d_step+j] = sigma;
			wcoef[(i-1)*d_step+j] <- waveletCoefDAtTime(times, values, time, sigma, 2*sigma)[1];
		}
	}
	dd <- data.frame(location, duration, wcoef);
	return (dd);
}
