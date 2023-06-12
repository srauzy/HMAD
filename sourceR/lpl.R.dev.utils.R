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

## lpl.R.dev.utils.R, creator S. Rauzy, LPL, 03/04/2019


##
## Create a data frame containing the values of the boundaries scale of interval of equal size 
## xmin : The minimal boundary of the scale
## xmax : The maximal boundary of the scale
## nbin : The number of bins of the scale
##
## return a data frame containing the boundaries
##
lpl.R.dev.utils.createBinBoundaries <- function(xmin, xmax, nbin) {

	deltabin <- (xmax - xmin)/nbin;
	min <- numeric(nbin);
	max <- numeric(nbin);
	mediane <- numeric(nbin);

	for (i in c(1:nbin)) {
		min[i] <- xmin + (i-1)*deltabin;
		max[i] <- xmin + i*deltabin;
		mediane[i] = min[i] + deltabin/2;
	}
	df <- data.frame(min, max, mediane);
	names(df) <- c("xmin", "xmax", "xmediane");
	rownames(df) <- NULL;

	return (df);
}

##
## Create a data frame containing the values of the boundaries scale 
## support : The support of the values 
## nbin : The number of bins of the scale
##
## return a data frame containing the boundaries
##
lpl.R.dev.utils.createBinBoundariesFromDensitySupport <- function(xmin, xmax, support, nbin) {

	nbrpointbin <- round(length(support)/nbin);
	min <- numeric(nbin);
	max <- numeric(nbin);
	mediane <- numeric(nbin);

	## Sort the support and put it in so
	so <- support[order(support)];
	min[1] <- min(so[1], xmin);

	for (i in c(1:(nbin-1))) {
		max[i] <- so[i*nbrpointbin];
		mediane[i] = so[i*nbrpointbin-round(0.5*nbrpointbin)];
		min[i+1] <- max[i];
		
	}

	mediane[nbin] = so[length(so)-round(0.5*nbrpointbin)];
	max[nbin] <- max(so[length(so)], xmax);

	df <- data.frame(min, max, mediane);
	names(df) <- c("xmin", "xmax", "xmediane");
	rownames(df) <- NULL;

	return (df);
}

##
## Discretize the values of a function following a scale
## boundaries : The data frame containing the boundaries of the scale
## f : The column with the value of the function
##
## return a data frame with the discrete values of the function, and the bin index 
##
lpl.R.dev.utils.discretizeFull <- function(boundaries, f) {
	fdiscrete <- numeric(length(f));
	binindex <- numeric(length(f));
	for (i in c(1:length(f))) {
		binindex[i] <- lpl.R.dev.utils.searchBinIndex(boundaries, f[i]);
		fdiscrete[i] <- boundaries$xmediane[binindex[i]];
	}
	df <- data.frame(f, fdiscrete, binindex);
	names(df) <- c("f", "fd", "bi");
	rownames(df) <- NULL;

	return (df);
}

##
## Discretize the values of a function following a scale
## boundaries : The data frame containing the boundaries of the scale
## f : The column with the value of the function
##
## return a column with the discrete values of the function 
##
lpl.R.dev.utils.discretize <- function(boundaries, f) {
	fdiscrete <- numeric(length(f));
	for (i in c(1:length(f))) {
		binindex <- lpl.R.dev.utils.searchBinIndex(boundaries, f[i]);
		fdiscrete[i] <- boundaries$xmediane[binindex];
	}
	return (fdiscrete);
}

##
## Discretize the values of a function following a scale
## boundaries : The data frame containing the boundaries of the scale
## f : The column with the value of the function
##
## return a column with the discrete values of the function 
##
lpl.R.dev.utils.searchBinIndex <- function(boundaries, value) {

	## Escape with the first bin if the value is not a number
	if (is.na(value)) return (1);

	for (i in c(1:(nrow(boundaries)-1))) {
		if (value < boundaries$xmin[i+1]) return (i);
	}
	return (nrow(boundaries));
}

##
## Compute the number of degrees of freedom for an rXc model (row times columns, r the number of groups and r the number of categories) 
##
## f : The column containing the category of the measurement
## fs : The column containing the group to which the measurment belong to
##
## Return the number of degrees of freedom for an rXc model
##
lpl.R.dev.utils.computeG2dof <- function(f, fs) {

	lfs <- levels(fs);
	lf <- levels(f);

	nr <- length(lfs);
	nc <- length(lf);

	return ((nr-1)*(nc-1));
}

##
## Compute the G2 statistic for an rXc model (row times columns, r the number of groups and r the number of categories) 
##
## f : The column containing the category of the measurement
## fs : The column containing the group to which the measurment belong to
##
## Return the G2 statistic
##
lpl.R.dev.utils.computeG2value <- function(f, fs) {

	## The different groups observed in the data
	lfs <- levels(fs);
	## The differents categories observed in the data
	lf <- levels(f);

	## The number of groups
	nr <- length(lfs);
	## The number of categories
	nc <- length(lf);

	m <- matrix(0, nrow=nr, ncol=nc);
	v <- numeric(nr);
	u <- numeric(nc);

	for (k in c(1:length(f))) {
		## Index of the group
		ii <- which(lfs == fs[k]);
		## Index of the category 
		ij <- which(lf == f[k]);
		## Count for group i and category j
		m[ii, ij] <- m[ii, ij]+1;
		## Count for group i 
		v[ii] <- v[ii]+1;
		## Count for category j 
		u[ij] <- u[ij]+1;
	}

	## Expected number of category j for group i (under the null hypothesis)
	e <- matrix(0, nrow=nr, ncol=nc);
	for (j in c(1:nc)) {
		for (i in c(1:nr)) {
			e[i, j] <- v[i]*u[j]/length(f);
		}
	}

	## the G2 statistic
	G2 <- 0;
	for (j in c(1:nc)) {
		for (i in c(1:nr)) {
			## If m is null m*log(m) tends towards 0
			if (m[i, j] != 0) { 
				G2 <- G2 + m[i, j]*log(m[i, j]/e[i, j]);
			}
		}
	}
	G2 <- 2*G2;

	## Return the value
	return (G2);
}

##
## Compute the entropy of a distribution
##
## f : The column containing the category of the measurement
##
## Return the entropy
##
lpl.R.dev.utils.computeEntropy <- function(f) {

	## The differents categories observed in the data
	lf <- levels(f);

	## The number of  observed categories
	nc <- length(lf);
	v <- numeric(nc);

	## The category density 
	for (k in c(1:length(f))) {
		ij <- which(lf == f[k]);
		v[ij] <- v[ij]+1;
	}
	v <- v/length(f);

	## The formulae for the entropy calculation
	H <- 0;
	for (j in c(1:nc)) {
		if (v[j] != 0) {
		       	H <- H + v[j]*log(v[j]);	
		}
	}
	H <- -H;

	## Return the entropy
	return (H);
}

##
## Compute the mutual entropy for an rXc model (row times columns, r the number of groups and r the number of categories) 
##
## f : The column containing the category of the measurement
## fs : The column containing the group to which the measurment belong to
##
## Return the mutual entropy
##
lpl.R.dev.utils.computeMutalEntropy <- function(f, fs) {

	## The different groups observed in the data
	lfs <- levels(fs);
	## The differents categories observed in the data
	lf <- levels(f);

	## The number of groups
	nr <- length(lfs);
	## The number of categories
	nc <- length(lf);

	m <- matrix(0, nrow=nr, ncol=nc);
	v <- numeric(nr);
	u <- numeric(nc);

	for (k in c(1:length(f))) {
		## Index of the group
		ii <- which(lfs == fs[k]);
		## Index of the category 
		ij <- which(lf == f[k]);
		## Count for group i and category j
		m[ii, ij] <- m[ii, ij]+1;
		## Count for group i 
		v[ii] <- v[ii]+1;
		## Count for category j
		u[ij] <- u[ij]+1;
	}

	## The density of categories for each group, the density of the categories and the density of the groups
	for (ii in c(1:nr)) {
		m[ii, ] <- m[ii, ]/v[ii];
		v[ii] <- v[ii]/length(f);
	}
	for (ij in c(1:nc)) {
		u[ij] <- u[ij]/length(f);
	}
	## The calculation of the entropy without group
	NH <- 0;
	for (j in c(1:nc)) {
		if (u[j] != 0) {
			NH <- NH + u[j]*log(u[j]);
		}
	}

	## The calculation of the splitted entropy
	SH <- 0;
	for (i in c(1:nr)) {
		H <- 0;
		for (j in c(1:nc)) {
			if (m[i, j] != 0) {
				H <- H + m[i, j]*log(m[i, j]);
			}
		}
		## H is the entropy of group i
		SH <- SH + v[i]*H;
	}

	## The calculation of the mutual entropy
	SH <- SH-NH;

	## Return the mutual entropy
	return (SH);
}

##
## Compute the splitted entropy for an rXc model (row times columns, r the number of groups and r the number of categories) 
##
## f : The column containing the category of the measurement
## fs : The column containing the group to which the measurment belong to
##
## Return the splitted entropy
##
lpl.R.dev.utils.computeSplittedEntropy <- function(f, fs) {

	## The different groups observed in the data
	lfs <- levels(fs);
	## The differents categories observed in the data
	lf <- levels(f);

	## The number of groups
	nr <- length(lfs);
	## The number of categories
	nc <- length(lf);

	m <- matrix(0, nrow=nr, ncol=nc);
	v <- numeric(nr);

	for (k in c(1:length(f))) {
		## Index of the group
		ii <- which(lfs == fs[k]);
		## Index of the category 
		ij <- which(lf == f[k]);
		## Count for group i and category j
		m[ii, ij] <- m[ii, ij]+1;
		## Count for group i 
		v[ii] <- v[ii]+1;
	}

	## The density of categories for each group and the density of the groups
	for (ii in c(1:nr)) {
		m[ii, ] <- m[ii, ]/v[ii];
		v[ii] <- v[ii]/length(f);
	}

	## The calculation of the mutual entropy
	SH <- 0;
	for (i in c(1:nr)) {
		H <- 0;
		for (j in c(1:nc)) {
			if (m[i, j] != 0) {
				H <- H + m[i, j]*log(m[i, j]);
			}
		}
		## H is the entropy of group i
		SH <- SH + v[i]*H;
	}
	SH <- -SH;

	## Return the mutual entropy
	return (SH);
}

##
## From a 3 columns data frame with tmin, tmax and the third column the value of a function,
## forms intervals such that the time contiguous intervals with the same value are merged 
## df : The 3 columns data frame
##
## Return the resulting 3 columns data frame 
##
lpl.R.dev.utils.createTimeIntervals <- function(df) {

	nr <- nrow(df);

	tmin <- df$tmin[1];
	tmax <- df$tmax[1];
	value <- df[1, 3];

	ndf <- NULL;
	for (i in c(2:nr)) {
		if (df[i, 3] != value | df$tmin[i] != tmax) {
			nl <- data.frame(tmin, tmax, value);
			ndf <- rbind(ndf, nl);

			tmin <- df$tmin[i];
			value <- df[i, 3];
		}
		tmax <- df$tmax[i];
	}
	nl <- data.frame(tmin, tmax, value);
	ndf <- rbind(ndf, nl);

	names(ndf) <- names(df);
	rownames(ndf) <- NULL;

	return (ndf);
}

##
## Create a data frame sampled on the times column with the values given by
## the time intervals of a data frame
## times : The times column to sample the value on
## df : The data frame containing the intervals and the associated values
##
## Return a data frame with 2 columns, times and associated values

lpl.R.dev.utils.projectTimeIntervals <- function(times, df) {

	values <- numeric(length(times));
	ci <- 1;
	ctmin <- df$tmin[ci];
	ctmax <- df$tmax[ci];
	value <- df[ci, 3];

	for (i in c(1:length(times))) {

		while (ctmax <= times[i] & ci < nrow(df)) {
			ci <- ci + 1;
			ctmin <- df$tmin[ci];
			ctmax <- df$tmax[ci];
			value <- df[ci, 3];
		}
		if (times[i] > ctmax) {
			values[i] <- NA;
		} else if (times[i] < ctmin) {
			values[i] <- NA;
		} else {
			values[i] <- value;
		}
	}
	ndf <- data.frame(times, values);
	names(ndf) <- c("time", names(df)[3]);
	rownames(ndf) <- NULL;

	return (ndf);
}

lpl.R.dev.utils.projectTimeIntervalsCharValues <- function(times, df) {

	values <- character(length(times));
	ci <- 1;
	ctmin <- df$tmin[ci];
	ctmax <- df$tmax[ci];
	value <- as.character(df[ci, 3]);

	for (i in c(1:length(times))) {

		while (ctmax <= times[i] & ci < nrow(df)) {
			ci <- ci + 1;
			ctmin <- df$tmin[ci];
			ctmax <- df$tmax[ci];
			value <- as.character(df[ci, 3]);
		}
		if (times[i] > ctmax) {
			values[i] <- "NA";
		} else if (times[i] < ctmin) {
			values[i] <- "NA";
		} else {
			values[i] <- value;
		}
	}
	ndf <- data.frame(times, values);
	names(ndf) <- c("time", names(df)[3]);
	rownames(ndf) <- NULL;

	return (ndf);
}

##
##
## Return a data frame with 2 columns, times and associated values

smoothCurve <- function(f, nscale) {

	values <- numeric(length(f));

	for (i in c(1:length(f))) {
		n <- 0;
		v <- 0;
		for (j in c(i-nscale, i+nscale)) {
			if (j > 0 & j <= length(f)) {
				n <- n + 1;
				v <- v + f[j];
			}
		}
		if (n == 0) {
			values[i] <- 0;
		} else {
			values[i] <- v/n;
		}
	}
	
	return (values);
}

transformValues <- function(values, transformationdf) {

	output <- character(length(values));
	for (i in c(1:length(values))) {
		j <- which(as.character(transformationdf$input) == as.character(values[i]));
		if (length(j) != 0) {
			output[i] <- as.character(transformationdf$output[j]);
		} else {
			output[i] <- NA;
		}
	}
	df <- data.frame(values, output);
	names(df) <- c("input", "output");
	rownames(df) <- NULL;

	return (df);
}


lpl.R.dev.utils.splitInContinuousIntervals <- function(values) {

	result <- NULL;
	index <- 1;
	v <- values[1];
	for (i in c(2:length(values))) {
		if (!is.na(values[i]) & values[i] != v) {
			result <- rbind(result, c(index, i-1, values[i-1]));
			v <- values[i];
			index <- i;
		}
	}
	result <- rbind(result, c(index, length(values), values[length(values)]));
	return (result);
}

lpl.R.dev.utils.splitCharInContinuousIntervals <- function(values) {

	result <- NULL;
	index <- 1;
	v <- as.character(values[1]);
	for (i in c(2:length(values))) {
		if (!is.na(values[i]) & values[i] != v) {
			result <- rbind(result, c(index, i-1, as.character(values[i-1])));
			v <- values[i];
			index <- i;
		}
	}
	result <- rbind(result, c(index, length(values), as.character(values[length(values)])));
	return (result);
}

transformColumnAtIndex <- function(data, columnIndex, transformationdf) {

	for (i in c(1:nrow(data))) {
		j <- which(as.character(transformationdf$input) == as.character(data[i, columnIndex]));
		if (length(j) != 0) {
			data[i, columnIndex] <- as.character(transformationdf$output[j]);
		} else {
			data[i, columnIndex] <- NA;
		}
		data[i, columnIndex]
	}

	return (data);
}
