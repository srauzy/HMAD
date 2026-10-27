##
## lpl.R.dev.hmm.supervisedLearning.R, creator S. Rauzy, LPL,06/12/2022
##

##
## Create the bin boundaries (of observation) with equals number of observations in each bin (except the last one)
## df : The data frame where to find the observations
## observation_column_index : The column index of the observations in df
## number_of_bins : The desired number of bins
## return a data frame with the name, min and max and mean observation values for each bin
##
lpl.R.dev.hmm.supervisedLearning.createBinBoundaries <- function(df, observation_column_index, number_of_bins) {

	observations <- df[ , observation_column_index];
	number_min_of_data_per_bin <- round(length(observations)/number_of_bins, 0);
	ordered_observations <- observations[order(observations)];
	
	f_min <- numeric(number_of_bins);
	f_max <- numeric(number_of_bins);
	count <-  numeric(number_of_bins);
	f <- numeric(number_of_bins);

	result <- NULL;

	f_min[1] <- -Inf;

	for (i_bin in c(1:(number_of_bins-1))) {
		f_max[i_bin] <- ordered_observations[i_bin*number_min_of_data_per_bin];
		df_bin <- subset(df, observations > f_min[i_bin] &  observations <= f_max[i_bin]);
		count[i_bin] <- nrow(df_bin);
		f[i_bin] <- mean(df_bin[ , observation_column_index]);

		line <- c(f_min[i_bin], f_max[i_bin], f[i_bin]); 
		result <- rbind(result, line);

		f_min[i_bin+1] <- f_max[i_bin]; 
	}

	f_max[number_of_bins] <- Inf;
	df_bin <- subset(df, observations > f_min[number_of_bins] &  observations <= f_max[number_of_bins]);
	f[number_of_bins] <- mean(df_bin[ , observation_column_index]);
	count[number_of_bins] <- nrow(df_bin);
	line <- c(f_min[number_of_bins], f_max[number_of_bins], f[number_of_bins]); 
	result <- rbind(result, line);

	df_result <- data.frame(result);
	line <- c("f_min", "f_max", "f");
	names(df_result) <- line;
	bins <-  as.character(number_of_bins);
	for (i_bin in c(1:(number_of_bins))) {
		bins[i_bin] <- paste("bin_", i_bin, sep=""); 
	}
	rownames(df_result) <- bins;

	return (df_result);
}

##
## Return the bin index corresponding to the value of the variable
## df_bins : The data frame defining the bins
## f_value : The value of the variable
## return the bin index corresponding to the value of the variable
##
lpl.R.dev.hmm.getBinIndex <- function(df_bins, f_value) {

	if (is.na(f_value)) return (1);

	for (i in c(1:nrow(df_bins))) {
		if (f_value < df_bins$f_min[i]) return (i);
	}

	return (nrow(df_bins));
}

##
## Return the column of bin indexes corresponding to the column values of the variable
## df_bins : The data frame defining the bins
## f_values : The column value of the variable
## return the column of bin indexes corresponding to the column values of the variable
##
lpl.R.dev.hmm.getBinIndexes <- function(df_bins, f_values) {

	bin_indexes <- numeric(length(f_values));
	for (i in c(1:length(f_values))) {
		bin_indexes[i] <- lpl.R.dev.hmm.getBinIndex(df_bins, f_values[i]); 
	}

	return (bin_indexes);
}

##
## Create the state emission matrix probability from the emission matrix count. That is not the emission density
## since it is not normalized by the bin size but for a given bin the ratio between the different probabilities
## for each state is correct (and that is what matters in the Viterbi algorithm 
## emission_matrix_count : The emission matrix count
## state_levels : The state levels
## bin_levels : The names of the bins
## return the emission matrix probability (rows are states and columns are bins in observation values) as a data frame
##
lpl.R.dev.hmm.supervisedLearning.createStateEmissionMatrixProbability <- function(emission_matrix_count, state_levels, bin_levels) {

	## For each bin (column), the sum of probabilities on the bins (columns) is 1
	for (j in c(1:ncol(emission_matrix_count))) {
		emission_matrix_count[, j] <- emission_matrix_count[, j]/sum(emission_matrix_count[, j]);
	}

	df_result <- data.frame(emission_matrix_count);
	names(df_result) <- bin_levels; 
	rownames(df_result) <- state_levels;

	return (df_result);
}

##
## Create the emission matrix probability from the emission matrix count. That is not the emission density
## since it is not normalized by the bin size but for a given bin the ratio between the different probabilities
## for each state is correct (and that is what matters in the Viterbi algorithm 
## emission_matrix_count : The emission matrix count
## state_levels : The state levels
## bin_levels : The names of the bins
## return the emission matrix probability (rows are states and columns are bins in observation values) as a data frame
##
lpl.R.dev.hmm.supervisedLearning.createEmissionMatrixProbability <- function(emission_matrix_count, state_levels, bin_levels) {

	## For each state (row), the sum of probabilities on the bins (columns) is 1
	for (i in c(1:nrow(emission_matrix_count))) {
		emission_matrix_count[i, ] <- emission_matrix_count[i, ]/sum(emission_matrix_count[i, ]);
	}

	df_result <- data.frame(emission_matrix_count);
	names(df_result) <- bin_levels;
	rownames(df_result) <- state_levels;

	return (df_result);
}

##
## Create the state levels with state values in the indicated column
## df : The data frame containing the state
## state_column_index : The index of the state column in the table
## return the state levels 
##
lpl.R.dev.hmm.supervisedLearning.createStateLevels <- function(df, state_column_index) {

	states <- df[ , state_column_index];
	state_levels <- levels(factor(states));

	return (state_levels);
}

##
## Create the emission matrix count for the annotations on a continous interval of frames with
## state values in the indicated column
## df : The data frame containing the state and observation column
## observation_column_index : The column index of the observations in df
## state_column_index : The index of the state column in the table
## df_bins : The data frame with the min and max observation values for each bin
## return the emission matrix count (rows are states and columns are bins in observation values)
##
lpl.R.dev.hmm.supervisedLearning.createEmissionMatrixCount <- function(df, observation_column_index, state_column_index, df_bins) {

	states <- df[ , state_column_index];
	state_levels <- levels(factor(states));
	n_states <- length(state_levels);

	n_bins <- nrow(df_bins);

	m_count <- matrix(0, nrow=n_states, ncol=n_bins);

	for (i_state in c(1:n_states)) {
		df_state <- subset(df, states == state_levels[i_state]);
		for (i_bin in c(1:n_bins)) {
			observations <- df_state[ , observation_column_index];
			df_state_bin <- subset(df_state, observations > df_bins$f_min[i_bin] &  observations <= df_bins$f_max[i_bin]);
			m_count[i_state, i_bin] <- nrow(df_state_bin); 
		}
	}

	return (m_count);
}

##
## Create the state transition matrix probability from the state transition matrix count
## m : The state transition matrix count
## state_levels : The state levels
## return the state transition matrix probability (rows are source/previous states and columns target/current states) as a data frame 
##
lpl.R.dev.hmm.supervisedLearning.createStateTransitionMatrixProbability <- function(m, state_levels) {

	## For each previous step (row), the sum of probabilities on the current states (columns) is 1
	for (i in c(1:nrow(m))) {
		m[i, ] <- m[i, ]/sum(m[i, ]);
	}

	df_result <- data.frame(m);
	names(df_result) <- state_levels;
	rownames(df_result) <- state_levels;

	return (df_result);
}

##
## Create the state transition matrix count for the annotations on a continous interval of frames with
## state values in the indicated column
## df_data : The data frame containing the state column
## state_column_index : The index of the state column in the table
## return the state transition matrix count (rows are source/previous states and columns target/current states)
##
lpl.R.dev.hmm.supervisedLearning.createStateTransitionMatrixCount <- function(df_data, state_column_index) {

	## The list of states
	state_levels <- levels(factor(df_data[ , state_column_index]));
	## The number of states
	n_states <- length(state_levels);

	## Create the previous state column by shifting from one frame in future
	## and removing the first entry (which has not previous state)
	previous_state <- c("", df_data[, state_column_index]);
	df_data$previous_state <- previous_state[-nrow(df_data)];
	df <- df_data[-1, ];

	## The count of couples (previous state, current state)
        ## Column for states, row for previous states
	m_count <- matrix(0, nrow=n_states, ncol=n_states);
	for (i in c(1:n_states)) {
		for (j in c(1:n_states)) {
			m_count[i, j] <- nrow(subset(df, previous_state == state_levels[i] & df[ , state_column_index] == state_levels[j]));
		}
	}

	return (m_count);	
}

##
## Create the emission matrix density and cdf for the annotations on a continous interval of frames with
## state values in the indicated column
## df : The data frame containing the state and observation column
## observation_column_index : The column index of the observations in df
## number_of_bins : The desired number of bins
## return the emission matrix count (rows are states and columns are bins in observation values)
##
lpl.R.dev.hmm.supervisedLearning.createEmissionDensity <- function(df, observation_column_index, number_of_bins) {

	observations <- df[ , observation_column_index];
	number_min_of_data_per_bin <- round(length(observations)/number_of_bins, 0);
	ordered_observations <- observations[order(observations)];
	
	f_min <- numeric(number_of_bins);
	f_max <- numeric(number_of_bins);
	cdf <-  numeric(number_of_bins);
	proba <-  numeric(number_of_bins);
	f <- numeric(number_of_bins);

	result <- NULL;

	f_min[1] <- -Inf;

	for (i_bin in c(1:(number_of_bins-1))) {
		f_max[i_bin] <- ordered_observations[i_bin*number_min_of_data_per_bin];
		df_bin <- subset(df, observations > f_min[i_bin] &  observations <= f_max[i_bin]);
		cdf[i_bin] <- i_bin/number_of_bins;
		f[i_bin] <- mean(df_bin[ , observation_column_index]);
	
		line <- c(f_min[i_bin], f_max[i_bin], f[i_bin], cdf[i_bin]); 
		result <- rbind(result, line);
		proba[i_bin] <- 1/(f_max[i_bin] - f_min[i_bin]);
		
		f_min[i_bin+1] <- f_max[i_bin]; 
	}

	## The last bin
	f_max[number_of_bins] <- Inf;
	df_bin <- subset(df, observations > f_min[number_of_bins] &  observations <= f_max[number_of_bins]);
	f[number_of_bins] <- mean(df_bin[ , observation_column_index]);
	cdf[number_of_bins] <- 1;
	line <- c(f_min[number_of_bins], f_max[number_of_bins], f[number_of_bins], cdf[number_of_bins]); 
	result <- rbind(result, line);

	## Interpolation of the probability on the first bin (in log scale for the proba)
	i_position <- 3;
	slope <- (log(proba[i_position]) - log(proba[i_position-1]))/(f[i_position] - f[i_position-1]);
	intercept <- log(proba[i_position]) - slope*f[i_position];
	proba[1] <- exp(slope*f[1] + intercept);

	## Interpolation of the probability on the last bin (in log scale for the proba)
	i_position <- number_of_bins-1;
	slope <- (log(proba[i_position]) - log(proba[i_position-1]))/(f[i_position] - f[i_position-1]);
	intercept <- log(proba[i_position]) - slope*f[i_position];
	proba[number_of_bins] <- exp(slope*f[number_of_bins] + intercept);

        proba <- proba/sum(proba);

	df_result <- data.frame(result, proba);
	line <- c("f_min", "f_max", "f", "cdf", "density");
	names(df_result) <- line;
	rownames(df_result) <- NULL;

	return (df_result);
}


lpl.R.dev.hmm.supervisedLearning.createBinsAndProbabilities <- function(df, state_column_index, observation_colun_index, number_of_bins) {

	state_levels <- levels(factor(df[ , state_column_index]));

	observations <- df[ , observation_column_index];
	number_min_of_data_per_bin <- round(length(observations)/number_of_bins, 0);
	ordered_observations <- observations[order(observations)];
	
	f_min <- numeric(number_of_bins);
	f_max <- numeric(number_of_bins);
	count <-  numeric(number_of_bins);
	f <- numeric(number_of_bins);

	result <- NULL;

	f_min[1] <- -Inf;

	for (i_bin in c(1:(number_of_bins-1))) {
		f_max[i_bin] <- ordered_observations[i_bin*number_min_of_data_per_bin];
		df_bin <- subset(df, observations > f_min[i_bin] &  observations <= f_max[i_bin]);
		count[i_bin] <- nrow(df_bin);
		f[i_bin] <- mean(df_bin[ , observation_column_index]);

		line <- c(f_min[i_bin], f_max[i_bin], f[i_bin]); 
			
		for (i_level in c(1:length(state_levels))) {
			df_bin_level <- subset(df_bin, df_bin[ , state_column_index] == state_levels[i_level]);
			proba_level <- nrow(df_bin_level)/nrow(df_bin);
			line <- c(line, proba_level);
		}

		result <- rbind(result, line);

		f_min[i_bin+1] <- f_max[i_bin]; 
	}

	f_max[number_of_bins] <- Inf;
	df_bin <- subset(df, observations > f_min[number_of_bins] &  observations <= f_max[number_of_bins]);
	f[number_of_bins] <- mean(df_bin[ , observation_column_index]);
	count[number_of_bins] <- nrow(df_bin);
	line <- c(f_min[number_of_bins], f_max[number_of_bins], f[number_of_bins]); 
	for (i_level in c(1:length(state_levels))) {
		df_bin_level <- subset(df_bin, df_bin[ , state_column_index] == state_levels[i_level]);
		proba_level <- nrow(df_bin_level)/nrow(df_bin);
		line <- c(line, proba_level);
	}

	result <- rbind(result, line);

	df_result <- data.frame(result);
	line <- c("f_min", "f_max", "f");
	for (i_level in c(1:length(state_levels))) {
		line <- c(line, paste("state", i_level, sep="_"));
	}

	names(df_result) <- line;

	return (df_result);
}

lpl.R.dev.hmm.supervisedLearning.computeBinCumulativeCountAndProbas <- function(df, state_column_index, state_value, observation_column_index, df_bins) {

	observations <- df[ , observation_column_index];
	
	count <- numeric(nrow(df_bins));	
	
	for (i_bin in c(1:nrow(df_bins))) {
		df_bin_data <- subset(df, observations > df_bins$f_min[i_bin] &  observations <= df_bins$f_max[i_bin]);
		states <- df_bin_data[ , state_column_index];
		df_selected <- subset(df_bin_data, states == state_value);
		count[i_bin] <- nrow(df_selected);
		if (i_bin != 1) {
			count[i_bin] <- count[i_bin] + count[i_bin-1]
		}

	}

	cdf <- count/count[nrow(df_bins)];

	df_result <- data.frame(df_bins, count, cdf);

	return (df_result);
}

lpl.R.dev.hmm.supervisedLearning.computeBinCountAndProbas <- function(df, state_column_index, state_value, observation_column_index, df_bins) {

	observations <- df[ , observation_column_index];
	
	count <- numeric(nrow(df_bins));	
	
	for (i_bin in c(1:nrow(df_bins))) {
		df_bin_data <- subset(df, observations > df_bins$f_min[i_bin] &  observations <= df_bins$f_max[i_bin]);
		states <- df_bin_data[ , state_column_index];
		df_selected <- subset(df_bin_data, states == state_value);
		count[i_bin] <- nrow(df_selected)/(df_bins$f_max[i_bin]-df_bins$f_min[i_bin]);
		if (i_bin == 1 | i_bin == nrow(df_bins)) count[i_bin] <- 0;
	}

	cdf <- count/sum(count);

	df_result <- data.frame(df_bins, count, cdf);

	return (df_result);
}


