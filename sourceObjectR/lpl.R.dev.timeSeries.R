##
## lpl.R.dev.timeSeries.R, creator S. Rauzy, LPL, 12/12/2022
## Various functions for performing time series manipulation 
##

##
## Compute the maximal time intervals (with respect to the value in the column of index column_index)
## and return the data frame of maximal intervals
## df : The data frame with column tmin and tmax defining the time interval and with values to sort
## column_index : The column ndex for the values to sort
## return the data frame of maximal intervals
##
lpl.R.dev.timeSeries.computeMaximalTimeIntervals <- function(df, column_index) {

	## Sort the table by decreasing order on the column column_index
	dfs <- df[order(-df[, column_index]), ];

	## Initialize the loop on the first line
	dfsd <- dfs;
	current_index <- 1;

	## The loop stops when the current line is the last one
 	while (current_index < nrow(dfsd)) {

		dfsd <- lpl.R.dev.timeSeries.removeOverlapIntervals(dfsd, current_index);
		## Go on in the table
		current_index <- current_index + 1;

	}

	## Sort the remaining intervals by increasing time
	dfsd <- dfsd[order(dfsd$tmin), ];

	## Return the result	
	return (dfsd);
}

##
## Remove from the data frame df the time intervals overlapping with df interval of index index
## df The data frame with the tmin and tmax columns
## index The line index of the time interval to consider
## Return the data frame with the overlapping intervals removed
##
lpl.R.dev.timeSeries.removeOverlapIntervals <- function(df, index) {

	tmin <- df$tmin[index];
	tmax <- df$tmax[index];
	remove <- lpl.R.dev.timeSeries.columnOverlap(tmin, tmax, df$tmin, df$tmax);
	remove[index] <- FALSE;

	return (subset(df, !remove));
}

## 
## Return an array of booleans indicating if the intervals [ctim2, ctmax2] overlap with interval [tmin, tmax]
## tmin1 The left boundary of the interval
## tax1 The right boundary of the interval
## ctmin2 The left boundaries of the intervals to check
## ctmax2 The rigth boundaries of the intervals to check
## Return a column of booleans indicating if the intervals [ctim2, ctmax2] overlap with interval [tmin, tmax]
##
lpl.R.dev.timeSeries.columnOverlap <- function(tmin1, tmax1, ctmin2, ctmax2) {

	return (ifelse(round(ctmax2, 2) <= round(tmin1,2) | round(tmax1, 2) <= round(ctmin2, 2), FALSE, TRUE));
}

##
## Interpolate for each time intervals listed in df_intervals the values of a function. 
## For a time interval with boudaries tmin and tmax (of duration delta_t = tmax - tmin), the interpolation starts at
## tmin - delta_t by replacing the values by the mean value between tmin and the current t, between tmin and tmax
## a linear interpolation and from tmax to tmax + delta_t the mean between tmax and the current t
## df : The data fame with the values to interpolate
## column_time_index : The time column index in df
## column_value_index : The values column index in df to interpolate
## df_intervals : The data frame containing the time intervals on which to interpolate
## return the column of interpolated values
##
lpl.R.dev.timeSeries.interpolateValues <- function(df, column_time_index, column_value_index, df_intervals) {

	df$index <- c(1:nrow(df));
	interpolated_values <- df[ , column_value_index];

	## For each interval to interpolate
	for (i in c(1:nrow(df_intervals))) {

		## The tmin and tmax boundaries of the interval
		t_min <- as.numeric(df_intervals$tmin[i]);
		t_max <- as.numeric(df_intervals$tmax[i]);

		## t_min_left is the tmin of the inteval minus the duration of the interval
		## t_max_right is the tmax of the inteval plus the duration of the interval
		duration <- (t_max - t_min);
		t_min_left <- max(0, t_min - duration);
	 	t_max_right <- min(df$time[nrow(df)], t_max + duration);

		## we select the interval before the interpolation area
		df_left <- subset(df, time < t_min & time > t_min_left);
		value_left <- df_left[, column_value_index];
		## That case happens when the first interpolation interval starts at time O for example
		if (nrow(df_left) == 0) {
			df_left <- subset(df, time == t_min_left);
		}

		current_interpolated_values <- numeric(nrow(df_left));
		for (k in c(1:nrow(df_left))) {
			current_interpolated_values[k] <- mean(value_left[1:k]);
			if (is.na(current_interpolated_values[k])) current_interpolated_values[k] <- 0;
		}
		interpolated_values[df_left$index[1]:df_left$index[nrow(df_left)]] <- current_interpolated_values;

		df_right <- subset(df, time < t_max_right & time > t_max);
		if (nrow(df_right) == 0) {
			df_right <- subset(df, time == t_max_right);
		}
		value_right <- df_right[, column_value_index];

		current_interpolated_values <- numeric(nrow(df_right));
		for (k in c(1:nrow(df_right))) {
			value_right <- df_right[, column_value_index];
			current_interpolated_values[k] <- mean(value_right[k:nrow(df_right)]);
			if (is.na(current_interpolated_values[k])) current_interpolated_values[k] <- 0;
		}
		interpolated_values[df_right$index[1]:df_right$index[nrow(df_right)]] <- current_interpolated_values;

		v_left <- interpolated_values[df_left$index[nrow(df_left)]];
		if (is.na(v_left)) v_left <- 0;
		v_right <- interpolated_values[df_right$index[1]];
		if (is.na(v_right)) v_right <- 0;
		t_left <- df$time[df_left$index[nrow(df_left)]];
		t_right <- df$time[df_right$index[1]];

		if (t_left == t_right) {
			slope <- 0;
			intercept <- 0;
		} else {
   	    		slope <- (v_right - v_left)/(t_right - t_left);
			intercept <- v_left - slope*t_left;
		}

		for (k in c((df_left$index[nrow(df_left)]+1):(df_right$index[1]-1))) {
			interpolated_values[k] <- slope*df[k, column_time_index] + intercept;
		}
	}

	return (interpolated_values);
}

##
## From the times and the values create an return the corresponding time intervals
## times : The timeline
## values : The values
## value_column_name : The column header name for the created data frame
## return the corresponding time intervals as a data frame with tmin, tmax and the values column
##
lpl.R.dev.timeSeries.createTimeAnnotationIntervals <- function(times, values, value_column_name) {

	nr <- length(times);
	tmin <- times[1];
	value <- values[1];
	ndf <- NULL;

	for (i in c(2:nr)) {
		if (values[i] != value) {
			ndf <- rbind(ndf, c(tmin, times[i], value));
			tmin <- times[i];
			value <- values[i];
		}
	}
	if (tmin < times[nr]) {
		ndf <- rbind(ndf, c(tmin, times[nr], value));
	}

	ndf <- data.frame(ndf);
	names(ndf) <- c("tmin", "tmax", value_column_name);
	rownames(ndf) <- NULL;

	return (ndf);
}
