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

## lpl.R.dev.smad.R, creator S. Rauzy, LPL, 12/04/2019

##
## Create step by step the automatic smile annotations in the SIS system. 
## Intermediate files are saved, so if the program is stopped when running, the computation will
## start for the next run after the last file created.
##
## software : The software used to created the output
## projectname : The name of the project
## df : The data frame containing the residuals and the column for filtering
## audf : The data frame containing AU values
##
lpl.R.dev.smad.computeSMAD <- function(software, projectname, df, audf) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdir <- paste("projects/", projectname, sep="");
	projectdir <- paste(projectdir, "/", software, sep="");

	## Define the folders	
	FOLDER_MODEL <- paste(projectdir, "models", sep="/");

	smile14dir <- paste(projectdir, "/", "smile14", sep="");
	if (!file.exists(smile14dir)) {
		dir.create(smile14dir);
	}
	FOLDER_TABLES <- paste(smile14dir, "tables", sep="/");
	if (!file.exists(FOLDER_TABLES)) {
		dir.create(FOLDER_TABLES);
	}

	FOLDER_SMILE14_MODEL <- paste(HMAD_DIRECTORY, "/", "models", "/", TRACKING_SOFTWARE, "/", "smile14", sep="");
	if (!exists("stmpdf1")) {
		lpl.R.dev.smad.initializeSmile14Model(FOLDER_SMILE14_MODEL, TRUE);
	}

	qualitythreshold = 0.3;

	## If the result file for automatic annotation of smiles has already been created, just load it
	if (fileExists(FOLDER_TABLES, "smad.txt")) {
		cat("Automatic annotations of smile previously done, loading the resulting file...\n");
		smad <- loadInternalDataFrame(FOLDER_TABLES, "smad.txt");
	} else {
		## A boolean filter function at true if the OpenFace confidence level is low < 0.8, the scale parameter
		## too low or high or the head position angles too far from the frontal pose
		filter <- lpl.dev.openFaceOutputAnalysis.createFilterConfidenceHeadSizeAndLargeAngle(df, 5);
		## Define a smooth function from this filter (of window size 0.4 seconds)
		minimalscale <- 0.4;
		df1Dqf <- lpl.R.dev.ebmotion.create1DQualityFunction(df, filter, minimalscale, 1);
		## Compute the best solution for the automatic smile annotation
		cat("Compute the smile intervals annotation from AU values...\n");
		smadsol <- lpl.R.dev.smad.computeSolution(df, audf, filter, FALSE); 
		## Project the solution
		smadp <- lpl.R.dev.utils.projectTimeIntervalsCharValues(df$time, smadsol);
		filterX <- (df1Dqf$quality < qualitythreshold);
		## Mark the X time intervals to be checked manually because not reliable
		smadp <- lpl.R.dev.smad.markXarea(smadp, 2, filterX);
		## Form the intervals
		sa <- lpl.R.dev.utils.splitCharInContinuousIntervals(smadp$valuef);
		tmin <- 0.04*(as.numeric(as.character(sa[, 1]))-1);
		tmax <- 0.04*as.numeric(as.character(sa[, 2]));
		smile14 <- sa[, 3];
		## Create the automatic annotation data frame (intervals annotated S0 to S4 or X)
		smad <- data.frame(tmin, tmax, smile14);

		## Save the result in a file
		saveInternalDataFrame(smad, FOLDER_TABLES, "smad.txt");
	}

	return (smad);
}

##
## Add to df a column named "valuef" which clones the column at clolumn_index in df replacing
## the annotation by "X" when the value of the filter is true
##
## df : The data frame to work with
## column_index : The column index of the colum in df to work with
## filter : A boolean indicating where to replace the value per "X"
##
## return the data frame df with the valuef column added
##
lpl.R.dev.smad.markXarea <- function(df, column_index, filter) {

	valuef <- character(nrow(df));

	for (i in c(1:nrow(df))) {
		if (filter[i]) {
			valuef[i] <- "X";
		} else {
			valuef[i] <- as.character(df[i, column_index]);
		}
	}
	return (data.frame(df, valuef));
}

##
## Load data and create global variables for the Smile model 
## function_directory The directory containing the files for the Smile model
## verbose Put this flag to true to have comments on console
##
lpl.R.dev.smad.initializeSmile14Model <- function(function_directory, verbose) {

	## The double arrow means that the variable is global (not limited to the scope of the function)
	alphadf <<- loadInternalDataFrame(function_directory, "coeffAlphas.txt");
	if (verbose) cat("Load and create global variables:\n- bindf\n- Variables for step 1 model: stmpdf1, ausmpdf1, svmdf1 and alpha1\n- Variables for step 2 model: stmpdf2, ausmpdf2, svmdf2 and alpha2\n- Variables for step 3 model: stmpdf3, ausmpdf3, svmdf3 and alpha3\n- Variables for step 4 model: stmpdf4, ausmpdf4, svmdf4 and alpha4\n");
	bindf <<- loadInternalDataFrame(function_directory, "binAUVC1.txt");

	stmpdf1 <<- loadInternalDataFrame(function_directory, "smileSTMP1.txt");
	ausmpdf1 <<- loadInternalDataFrame(function_directory, "smileAUSMP1.txt");
	svmdf1 <<- loadInternalDataFrame(function_directory, "coeffAUVC1.txt");
	alpha1 <<- alphadf$alphamax[1];

	stmpdf2 <<- loadInternalDataFrame(function_directory, "smileSTMP2.txt");
	ausmpdf2 <<- loadInternalDataFrame(function_directory, "smileAUSMP2.txt");
	svmdf2 <<- loadInternalDataFrame(function_directory, "coeffAUVC2.txt");
	alpha2 <<- alphadf$alphamax[2];

	stmpdf3 <<- loadInternalDataFrame(function_directory, "smileSTMP3.txt");
	ausmpdf3 <<- loadInternalDataFrame(function_directory, "smileAUSMP3.txt");
	svmdf3 <<- loadInternalDataFrame(function_directory, "coeffAUVC3.txt");
	alpha3 <<- alphadf$alphamax[3];

	stmpdf4 <<- loadInternalDataFrame(function_directory, "smileSTMP4.txt");
	ausmpdf4 <<- loadInternalDataFrame(function_directory, "smileAUSMP4.txt");
	svmdf4 <<- loadInternalDataFrame(function_directory, "coeffAUVC4.txt");
	alpha4 <<- alphadf$alphamax[4];
}

##
## Annotate automatically the smile intervals and return the resulting data frame (need the Smile model to be initialized by a call to lpl.R.dev.smad.initializeSmile14Model(function_directory, verbose))
## df : The data frame containing the HMAD formatted output of the ilandmarks and OPENFACE analysis
## audf : The data frame containing the HMAD formatted output of the OPENFACE analysis for the AU
## filter : The filter which inform whether the AU data are fair 
## verbose Put this flag to true to have comments on console
##
lpl.R.dev.smad.computeSolution <- function(df, audf, filter, verbose) {

	imin <- 1;
	imax <- nrow(df);
	value <- "S";
	isol0 <- data.frame(imin, imax, value);
	intervalsdf <- isol0;

	label <- "S";
	sol1 <- lpl.R.dev.smad.computeBestSolutionOnIntervals(audf, filter, intervalsdf, label, bindf, stmpdf1, ausmpdf1, svmdf1, alpha1, verbose);

	isol1 <- lpl.R.dev.utils.splitInContinuousIntervals(sol1);
	input <- c(1, 2);
	output <- c("S0", "S");
	transformationdf <- data.frame(input, output);
	isol1t <- transformColumnAtIndex(isol1, 3, transformationdf)

	intervalsdf <- isol1t;
	label <- "S";
	sol2 <- lpl.R.dev.smad.computeBestSolutionOnIntervals(audf, filter, intervalsdf, label, bindf, stmpdf2, ausmpdf2, svmdf2, alpha2, verbose);

	isol2 <- lpl.R.dev.utils.splitInContinuousIntervals(sol2);
	input <- c("S0", "1", "2");
	output <- c("S0", "S1", "S");
	transformationdf <- data.frame(input, output);
	isol2t <- transformColumnAtIndex(isol2, 3, transformationdf)

	intervalsdf <- isol2t;
	label <- "S";
	sol3 <- lpl.R.dev.smad.computeBestSolutionOnIntervals(audf, filter, intervalsdf, label, bindf, stmpdf3, ausmpdf3, svmdf3, alpha3, verbose);
	
	isol3 <- lpl.R.dev.utils.splitInContinuousIntervals(sol3);
	input <- c("S0", "S1", "1", "2");
	output <- c("S0", "S1", "S2", "S");
	transformationdf <- data.frame(input, output);
	isol3t <- transformColumnAtIndex(isol3, 3, transformationdf)

	intervalsdf <- isol3t;
	label <- "S";
	sol4 <- lpl.R.dev.smad.computeBestSolutionOnIntervals(audf, filter, intervalsdf, label, bindf, stmpdf4, ausmpdf4, svmdf4, alpha4, verbose);

	isol4 <- lpl.R.dev.utils.splitInContinuousIntervals(sol4);			
	input <- c("S0", "S1", "S2", "1", "2");
	output <- c("S0", "S1", "S2", "S3", "S4");
	transformationdf <- data.frame(input, output);
	isol4t <- transformColumnAtIndex(isol4, 3, transformationdf)

	tmin <- 0.04*(as.numeric(as.character(isol4t[,1]))-1);
	tmax <- 0.04*as.numeric(as.character(isol4t[,2]));
	solution <- data.frame(tmin, tmax, as.character(isol4t[, 3]));
	names(solution) <- c("tmin", "tmax", "value");

	return (solution);
}

##
## Annotate automatically the smile intervals for a given step and return the resulting data frame
## audf : The data frame containing the HMAD formatted output of the OPENFACE analysis for the AU
## filter : For time position where the filter is FALSE, does not use AU values and just make use of the state transition probability of the model
## intervalsdf : The data frame containing the intervals of time where to compute the solution
## label : The label of the value column in the intervalsdf data frame corresponding to the intervals to compute the solution in 
## bindf : The bins data frame for the AU function
## stmpdf : The state transition matrix for the given step
## ausmpdf : The AU state matrix of probability for the given step
## svmdf : The matrix of AU coefficient to define the AU variable for the given step
## alpha : The coefficient alpha of the model for the given step
## verbose Put this flag to true to have comments on console
##
lpl.R.dev.smad.computeBestSolutionOnIntervals <- function(audf, filter, intervalsdf, label, bindf, stmpdf, ausmpdf, svmdf, alpha, verbose) {

	svdf <- lpl.R.dev.smad.createSmile14Variable(svmdf, audf, bindf, verbose);

	sol <- NULL;
	for (i in c(1:nrow(intervalsdf))) {
		if (verbose) {
			cat(paste("tmin = ", intervalsdf[i, 1], ", tmax = ", intervalsdf[i, 2], ", value = ", intervalsdf[i, 3], "\n", sep=""));
		}
		if (!is.na(intervalsdf[i, 3]) & intervalsdf[i, 3] == label) {
			soli <- lpl.R.dev.smad.bestSolution(stmpdf, ausmpdf, alpha, svdf, filter, as.numeric(intervalsdf[i, 1]), as.numeric(intervalsdf[i, 2]), verbose);
		} else {
			soli <- NULL;
			for (j in c(as.numeric(intervalsdf[i, 1]):as.numeric(intervalsdf[i, 2]))) {
				soli <- c(soli, as.character(intervalsdf[i, 3]));
			}
		}
		sol <- c(sol, soli);
	}

	return (sol);	
}

##
## Load the Hidden Markov Model for smile detection with 4 levels
## hmad_directory : The HMAD directory
## tracking_software : The tracking software (OPENFACE)
## verbose : Put to TRUE for verbose mode
##
## Return the data frame containing the probabilities of transition between the 5 states of the HMM for each bin of the AU composite 
## variable
##
lpl.R.dev.smad.loadSmile14HiddenMarkovModel <- function(hmad_directory, tracking_software, verbose) {

	## The folder for the model
	function_directory <- paste(hmad_directory, "/", "models", "/", tracking_software, "/", "smile14", sep="") 

	## Load the file defining the coefficients of the AU for defining the smile variable
	austmdf <- loadInternalDataFrame(function_directory, "smileAUSTM.txt");

	return (austmdf);
}

##
## Load the smile variable model 
## hmad_directory : The HMAD directory
## tracking_software : The tracking software (OPENFACE)
## verbose : Put to TRUE for verbose mode
##
## Return data frame containing the coefficients of the AU for the smile variable model
##
lpl.R.dev.smad.loadSmile14VariableModel <- function(hmad_directory, tracking_software, verbose) {

	## The folder for the model
	function_directory  <- paste(hmad_directory, "/", "models", "/", tracking_software, "/", "smile14", sep="") 

	## Load the file defining the coefficients of the AU for defining the smile variable
	svmdf <- loadInternalDataFrame(function_directory, "smileAUVC.txt");

	return (svmdf);
}

##
## Load the file for bins definition of the smile variable
## hmad_directory : The HMAD directory
## tracking_software : The tracking software (OPENFACE)
## verbose : Put to TRUE for verbose mode
##
## Return the data frame containing the bins definition for the smile composite variable
##
lpl.R.dev.smad.loadSmile14BinsDefinition <- function(hmad_directory, tracking_software, verbose) {

	## The folder for the model
	function_directory <- paste(hmad_directory, "/", "models", "/", tracking_software, "/", "smile14", sep="") 

	## Load the file defining the coefficients of the AU for defining the smile variable
	bindf <- loadInternalDataFrame(function_directory, "binAUVC.txt");

	return (bindf);
}

##
## Create the smile variable from the smile variable model svmdf and the OPENFACE AU data frame audf
## svmdf : The smile variable model
## audf : The OPENFACE AU data frame 
## bindf : The data frame containing the bins definition 
## verbose : Put to TRUE for verbose mode
##
## Return a data frame with the smile variable and its discrete version on 20 bins
##
lpl.R.dev.smad.createSmile14Variable <- function(svmdf, audf, bindf, verbose) {

	## The column indexes of the AU in the audf data frame
	aunames <- c("AU06", "AU07", "AU10", "AU12", "AU14", "AU25", "AU26");
	column_index <- numeric(length(aunames))
	for (k in c(1:length(aunames))) {
		column_index[k] = which(colnames(audf) == paste(aunames[k], "_r", sep=""));
	}

	## The variable for each frame
	v <- numeric(nrow(audf));
	for (k in c(1:length(column_index))) {
		v <- v + svmdf$coeff[k]*audf[, column_index[k]];
	}

	## Discrete version of the variable in a data frame
	vddf <- lpl.R.dev.utils.discretizeFull(bindf, v);

	return (vddf);
}

##
## Create the smile variable from the smile variable model svmdf and the OPENFACE AU data frame audf
## svmdf : The smile variable model
## audf : The OPENFACE AU data frame 
## bindf : The data frame containing the bins definition 
## verbose : Put to TRUE for verbose mode
##
## Return a data frame with the smile variable and its discrete version on 20 bins
##
lpl.R.dev.smad.createSmile14VariableBis <- function(svmdf, audf, bindf, verbose) {

	## The column indexes of the AU in the audf data frame
	aunames <- c("AU06", "AU07", "AU10", "AU12", "AU14", "AU25", "AU26");
	column_index <- numeric(length(aunames))
	for (k in c(1:length(aunames))) {
		column_index[k] = which(colnames(audf) == aunames[k]);
	}

	## The variable for each frame
	v <- numeric(nrow(audf));
	for (k in c(1:length(column_index))) {
		v <- v + svmdf$coeff[k]*audf[, column_index[k]];
	}

	## Discrete version of the variable in a data frame
	vddf <- lpl.R.dev.utils.discretizeFull(bindf, v);

	return (vddf);
}

##
## Return the state probability of transition
## bindf : The data frame containing the bins definition
## austmdf : The data frame containing the state probabilities of transitions
## svdf : The data frame containing the smile variable and its associated bins
##
## Return the Viterbi solution
##
lpl.R.dev.smad.solutionTable <- function(audf, solution) {

	timestep <- round(audf$time[2] - audf$time[1], 2);
	tmin <- round(audf$time, 2);
	tmax <- round(tmin + timestep, 2);

	dfu <- data.frame(tmin, tmax, solution);

	return (lpl.R.dev.utils.createTimeIntervals(df));
}

##
## bindf : The bins data frame for the AU function
## stmpdf : The state transition matrix for the given step
## ausmpdf : The AU state matrix of probability for the given step
## svmdf : The matrix of AU coefficient to define the AU variable for the given step
## alpha : The coefficient alpha of the model for the given step
## verbose Put this flag to true to have comments on console
##
## Return the Viterbi solution
##
lpl.R.dev.smad.bestSolution <- function(stmpdf, ausmpdf, alpha, svdf, filter, itmin, itmax, verbose) {

	## The number of states
	ns <- ncol(stmpdf);
	## The number of frames
	##n <- nrow(svdf);
	n <- itmax - itmin + 1;

	NOT_ACTIVE <- 1;
	NOT_YET_VISITED <- -2;

	max_logproba <- matrix(NOT_ACTIVE, nrow=ns, ncol=n);
	max_index <- matrix(NOT_YET_VISITED, nrow=ns, ncol=n);

	## We start with a uniform distribution of  states (index 1)
	j <- itmin;
	jloop <- 1;
	## Loop on source
	nsa <- 0;
	for (k in c(1:ns)) {
		## Loop on target
		for (i in c(1:ns)) {
			## The probability to transit from source state k to target state i with observation value svdf$bi[j] at position j weighted by the uniform distribution 
			pt <- lpl.R.dev.smad.getProbabilityOfTransition(stmpdf, ausmpdf, alpha, filter[j], svdf$bi[j], k, i)/ns;
			## If the target state is reachable
			if (pt != 0) {
				##cumul_logproba[i, j] <- log(pt);
			       if (max_index[i, jloop] == NOT_YET_VISITED) {	
					max_logproba[i, jloop]  <- log(pt);
					max_index[i, jloop] == k;
			       } else {
				       if (log(pt) > max_logproba[i, jloop]) {
					       max_logproba[i, jloop]  <- log(pt);
					       max_index[i, jloop] == k;
				       }
			       }
			## The target state is not reachable, put its log to 1 in order to identify the case  
			} else {
				max_logproba[i, jloop] <- NOT_ACTIVE;
			}	
			max_index[i, jloop] <- -1;
		}
	}

	if (verbose) cat(paste("Mean value for AU function : ", mean(svdf$f[itmin:itmax]), sep=""));
	## Loop on the time (or position)
	for (j in c((itmin+1):itmax)) {
		jloop <- j - itmin + 1;
		## Loop on source
		nsa <- 0;
		for (k in c(1:ns)) {
			## If the source state is active
			if (max_logproba[k, jloop-1] != NOT_ACTIVE) {
				np0 <- 0;
				## Loop on target
				for (i in c(1:ns)) {
					pt <- lpl.R.dev.smad.getProbabilityOfTransition(stmpdf, ausmpdf, alpha, filter[j], svdf$bi[j], k, i);
					##cat(paste("i", pt, "\n", sep="\t"));
					if (is.na(pt) || pt < 0) {
						cat(paste("warning", pt, "\n", sep=" "));
					}
					## If the target state is reachable
					if (pt != 0) {
						if (max_index[i, jloop] == NOT_YET_VISITED) {
							##cumul_logproba[i, j] <- log(pt) + cumul_logproba[k, j-1];
						} else {
							##cumul_logproba[i, j] <- log(exp(cumul_logproba[i, j]) + exp(log(pt) + cumul_logproba[k, j-1]));
						}

						if (max_index[i, jloop] == NOT_YET_VISITED) {
						       	max_index[i, jloop] <- k;
							max_logproba[i, jloop] <- log(pt) + max_logproba[k, jloop-1];
						} else if ((log(pt) + max_logproba[k, jloop-1]) > max_logproba[i, jloop]) {
							max_index[i, jloop] <- k;
							max_logproba[i, jloop] <- log(pt) + max_logproba[k, jloop-1];
						}
					} else {
						np0 <- np0 + 1;
					}
				}
				if (np0 == 5) {
				cat(paste("j = ", j, " # p = 0 ", np0, "\n", sep=""));
		}
			} else {
				nsa <- nsa + 1;
			}
		}
		if (nsa == ns) {
				cat(paste("j = ", j, " # not active ", nsa, "\n", sep=""));
		}
	}

	j <- itmax;
	jloop <- j - itmin + 1; 
	vmax <- NOT_YET_VISITED;
	imax <- 0;
	for (i in c(1:ns)) {
		if (vmax == NOT_YET_VISITED) {
			vmax <- max_logproba[i, jloop];
			imax <- i;
		} else {
			if (max_logproba[i, jloop] > vmax) {
				vmax <- max_logproba[i, jloop];
				imax <- i;
				
			}
		}
		##cat(paste(max_logproba[i, j], "\n", sep=""));
	}
	##cat(paste("max = ", vmax, " for state ", imax, "\n", sep=""));

	sol <- numeric(n);
	lpmax <- numeric(n);
	j <- itmax;
	jloop <- j - itmin + 1; 
	sol[jloop] <- imax;
	lpmax[jloop] <- max_logproba[sol[jloop],jloop];
	for (j in c((itmax-1):itmin)) {
		jloop <- j - itmin + 1;
		sol[jloop] <- max_index[sol[jloop+1], jloop+1];
		##cat(paste("max = ",j, " state ", sol[j], "\n", sep=""));
		lpmax[jloop] <- max_logproba[sol[jloop], jloop];	
	}
	return (sol);
}

##
## Return the state probability of transition
## bindf : The data frame containing the bins definition
## austmdf : The data frame containing the state probabilities of transitions
## svdf : The data frame containing the smile variable and its associated bins
##
## Return the Viterbi solution
##
lpl.R.dev.smad.bestSolutionFromAUFunction <- function(stmpdf, ausmpdf, alpha, svdf, filter, itmin, itmax, verbose) {

	## The number of states
	ns <- ncol(stmpdf);
	## The number of frames
	##n <- nrow(svdf);
	n <- itmax - itmin + 1;

	NOT_ACTIVE <- 1;
	NOT_YET_VISITED <- -2;

	max_logproba <- matrix(NOT_ACTIVE, nrow=ns, ncol=n);
	max_index <- matrix(NOT_YET_VISITED, nrow=ns, ncol=n);

	## We start with a uniform distribution of  states (index 1)
	j <- itmin;
	jloop <- 1;
	## Loop on source
	nsa <- 0;
	for (k in c(1:ns)) {
		## Loop on target
		for (i in c(1:ns)) {
			## The probability to transit from source state k to target state i with observation value svdf$bi[j] at position j weighted by the uniform distribution 
			pt <- lpl.R.dev.smad.getProbabilityOfTransition(stmpdf, ausmpdf, alpha, filter[j], svdf$bi[j], k, i)/ns;
			## If the target state is reachable
			if (pt != 0) {
				##cumul_logproba[i, j] <- log(pt);
			       if (max_index[i, jloop] == NOT_YET_VISITED) {	
					max_logproba[i, jloop]  <- log(pt);
					max_index[i, jloop] == k;
			       } else {
				       if (log(pt) > max_logproba[i, jloop]) {
					       max_logproba[i, jloop]  <- log(pt);
					       max_index[i, jloop] == k;
				       }
			       }
			## The target state is not reachable, put its log to 1 in order to identify the case  
			} else {
				max_logproba[i, jloop] <- NOT_ACTIVE;
			}	
			max_index[i, jloop] <- -1;
		}
	}

	if (verbose) cat(paste("Mean value for AU function : ", mean(svdf$f[itmin:itmax]), sep=""));
	## Loop on the time (or position)
	for (j in c((itmin+1):itmax)) {
		jloop <- j - itmin + 1;
		## Loop on source
		nsa <- 0;
		for (k in c(1:ns)) {
			## If the source state is active
			if (max_logproba[k, jloop-1] != NOT_ACTIVE) {
				np0 <- 0;
				## Loop on target
				for (i in c(1:ns)) {
					pt <- lpl.R.dev.smad.getProbabilityOfTransition(stmpdf, ausmpdf, alpha, filter[j], svdf$bi[j], k, i);
					##cat(paste("i", pt, "\n", sep="\t"));
					if (is.na(pt) || pt < 0) {
						cat(paste("warning", pt, "\n", sep=" "));
					}
					## If the target state is reachable
					if (pt != 0) {
						if (max_index[i, jloop] == NOT_YET_VISITED) {
							##cumul_logproba[i, j] <- log(pt) + cumul_logproba[k, j-1];
						} else {
							##cumul_logproba[i, j] <- log(exp(cumul_logproba[i, j]) + exp(log(pt) + cumul_logproba[k, j-1]));
						}

						if (max_index[i, jloop] == NOT_YET_VISITED) {
						       	max_index[i, jloop] <- k;
							max_logproba[i, jloop] <- log(pt) + max_logproba[k, jloop-1];
						} else if ((log(pt) + max_logproba[k, jloop-1]) > max_logproba[i, jloop]) {
							max_index[i, jloop] <- k;
							max_logproba[i, jloop] <- log(pt) + max_logproba[k, jloop-1];
						}
					} else {
						np0 <- np0 + 1;
					}
				}
				if (np0 == 5) {
				cat(paste("j = ", j, " # p = 0 ", np0, "\n", sep=""));
		}
			} else {
				nsa <- nsa + 1;
			}
		}
		if (nsa == ns) {
				cat(paste("j = ", j, " # not active ", nsa, "\n", sep=""));
		}
	}

	j <- itmax;
	jloop <- j - itmin + 1; 
	vmax <- NOT_YET_VISITED;
	imax <- 0;
	for (i in c(1:ns)) {
		if (vmax == NOT_YET_VISITED) {
			vmax <- max_logproba[i, jloop];
			imax <- i;
		} else {
			if (max_logproba[i, jloop] > vmax) {
				vmax <- max_logproba[i, jloop];
				imax <- i;
				
			}
		}
		##cat(paste(max_logproba[i, j], "\n", sep=""));
	}
	##cat(paste("max = ", vmax, " for state ", imax, "\n", sep=""));

	sol <- numeric(n);
	lpmax <- numeric(n);
	j <- itmax;
	jloop <- j - itmin + 1; 
	sol[jloop] <- imax;
	lpmax[jloop] <- max_logproba[sol[jloop],jloop];
	for (j in c((itmax-1):itmin)) {
		jloop <- j - itmin + 1;
		sol[jloop] <- max_index[sol[jloop+1], jloop+1];
		##cat(paste("max = ",j, " state ", sol[j], "\n", sep=""));
		lpmax[jloop] <- max_logproba[sol[jloop], jloop];	
	}
	return (sol);
}

##
## Return the state probability of transition
## bindf : The data frame containing the bins definition
## austmdf : The data frame containing the state probabilities of transitions
## ibin : The bin index for the smile variable (if -1, it returns the value ignoring the smile variable value)
## istate : The source state index
## itarget : The targe state index
##
## Return the target state probability of transition 
##
lpl.R.dev.smad.getStateProbabilityOfTransition <- function(bindf, austmdf, ibin, istate, itarget) {

	nb <- nrow(bindf);
	if (ibin == -1) {
		j <- nb+1;
	} else {
		j <- ibin;
	}

	return (austmdf[(istate-1)*(nb+1)+j, itarget+2]);
}

##
## Return the state probability of transition
## bindf : The data frame containing the bins definition
## austmdf : The data frame containing the state probabilities of transitions
## ibin : The bin index for the smile variable (if -1, it returns the value ignoring the smile variable value)
## istate : The source state index
## itarget : The targe state index
##
## Return the target state probability of transition 
##
lpl.R.dev.smad.getProbabilityOfTransition <- function(stmpdf, ausmpdf, alpha, filter, ibin, istate, itarget) {

	if (filter | ibin == -1) {
		p <- stmpdf[istate, itarget];
	} else {
		p <- alpha*stmpdf[istate, itarget] + (1-alpha)*ausmpdf[ibin, itarget];
	}

	return (p);
}

##
## Return the state probability of transition
## bindf : The data frame containing the bins definition
## austmdf : The data frame containing the state probabilities of transitions
## ibin : The bin index for the smile variable (if -1, it returns the value ignoring the smile variable value)
## istate : The source state index
## itarget : The targe state index
##
## Return the target state probability of transition 
##
lpl.R.dev.smad.getProbabilityOfTransitionFromAUFunction <- function(ftdf, filter, ibin, istate, itarget) {

	if (filter | ibin == -1) {
		p <- stmpdf[istate, itarget];
	} else {
		p <- alpha*stmpdf[istate, itarget] + (1-alpha)*ausmpdf[ibin, itarget];
	}

	return (p);
}

lpl.R.dev.smad.discretizeTable <- function(dfx, dfbin, verbose) {
	
	## The discrete AU values as data frame
	dfxd <- NULL; 
	for (i in c(1:ncol(dfx))) {
		if (is.null(dfxd)) {
			dfxd <- data.frame(lpl.R.dev.utils.discretize(dfbin, dfx[,i]));
		} else {
			dfxd <- cbind(dfxd, lpl.R.dev.utils.discretize(dfbin, dfx[,i]));
		}
		colnames(dfxd)[i] <- colnames(dfx)[i];
	}
	rownames(dfxd) <- NULL;

	return (dfxd);
}

##
## Return the state probability of transition
## bindf : The data frame containing the bins definition
## austmdf : The data frame containing the state probabilities of transitions
## svdf : The data frame containing the smile variable and its associated bins
##
## Return the Viterbi solution
##
lpl.R.dev.smad.coverageSolution <- function(austmdf, sol, imin, imax) {

	ns <- ncol(austmdf) - 2;
	n <- nrow(imax-imin);

	NOT_ACTIVE <- 1;
	NOT_YET_VISITED <- -2;

	max_logproba <- matrix(NOT_ACTIVE, nrow=ns, ncol=n);
	max_index <- matrix(NOT_YET_VISITED, nrow=ns, ncol=n);

	## We start with a S0 state (index 1)
	j <- 1;
	## Loop on target
	for (i in c(1:ns)) {
		pt <- lpl.R.dev.smad.getStateProbabilityOfTransition(bindf, austmdf, svdf$bi[j], 1, i);
		## If the target state is reachable
		if (pt != 0) {
			##cumul_logproba[i, j] <- log(pt); 
			max_logproba[i, j]  <- log(pt);
		## The target state is not reachable, put its log to 1 in order to identify the case  
		} else {
			max_logproba[i, j] <- NOT_ACTIVE;
		}	
		max_index[i, j] <- -1;
	}

	## Loop on the time
	for (j in c(2:n)) {
		## Loop on source
		nsa <- 0;
		for (k in c(1:ns)) {
			## If the source state is active
			if (max_logproba[k, j-1] != NOT_ACTIVE) {
				np0 <- 0;
				## Loop on target
				for (i in c(1:ns)) {
					pt <- lpl.R.dev.smad.getStateProbabilityOfTransition(bindf, austmdf, svdf$bi[j], k, i);
					##cat(paste("i", pt, "\n", sep="\t"));
					if (is.na(pt) || pt < 0) {
						cat(paste("warning", pt, "\n", sep=" "));
					}
					## If the target state is reachable
					if (pt != 0) {
						if (max_index[i, j] == NOT_YET_VISITED) {
							##cumul_logproba[i, j] <- log(pt) + cumul_logproba[k, j-1];
						} else {
							##cumul_logproba[i, j] <- log(exp(cumul_logproba[i, j]) + exp(log(pt) + cumul_logproba[k, j-1]));
						}

						if (max_index[i, j] == NOT_YET_VISITED) {
						       	max_index[i, j] <- k;
							max_logproba[i, j] <- log(pt) + max_logproba[k, j-1];
						} else if ((log(pt) + max_logproba[k, j-1]) > max_logproba[i, j]) {
							max_index[i, j] <- k;
							max_logproba[i, j] <- log(pt) + max_logproba[k, j-1];
						}
					} else {
						np0 <- np0 + 1;
					}
				}
				if (np0 == 5) {
				cat(paste("j = ", j, " # p = 0 ", np0, "\n", sep=""));
		}
			} else {
				nsa <- nsa + 1;
			}
		}
		if (nsa == ns) {
				cat(paste("j = ", j, " # not active ", nsa, "\n", sep=""));
		}
	}

	j <- n;
	vmax <- NOT_YET_VISITED;
	imax <- 0;
	for (i in c(1:ns)) {
		if (vmax == NOT_YET_VISITED) {
			vmax <- max_logproba[i, j];
			imax <- i;
		} else {
			if (max_logproba[i, j] > vmax) {
				vmax <- max_logproba[i, j];
				imax <- i;
				
			}
		}
		##cat(paste(max_logproba[i, j], "\n", sep=""));
	}
	cat(paste("max = ", vmax, " for state ", imax, "\n", sep=""));

	sol <- numeric(n);
	lpmax <- numeric(n);
	sol[n] <- imax;
	lpmax[n] <- max_logproba[sol[n], n];
	for (j in c((n-1):1)) {
		sol[j] <- max_index[sol[j+1], j+1];
		##cat(paste("max = ",j, " state ", sol[j], "\n", sep=""));
		lpmax[j] <- max_logproba[sol[j], j];	
	}
	return (sol);
}

##
## Create step by step (first the eyebrow raise automatic annotation and second eyebrow frown one) the automatic
## annotations. Intermediate files are saved, so if the program is stopped when running, the computation will
## start for the next run after the last file created.
##
## software : The software used to created the output
## projectname : The name of the project
## df : The data frame containing the residuals and the column for filtering
## audf : The data frame containing AU values
##
lpl.R.dev.smad.computeAUSMAD <- function(software, projectname, df, audf) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdir <- paste("projects/", projectname, sep="");
	projectdir <- paste(projectdir, "/", software, sep="");

	## Define the folders	
	FOLDER_MODEL <- paste(projectdir, "model", sep="/");

	s4maddir <- paste(projectdir, "/", "S4MAD", sep="");
	if (!file.exists(s4maddir)) {
		dir.create(s4maddir);
	}
	FOLDER_TABLES <- paste(s4maddir, "tables", sep="/");
	if (!file.exists(FOLDER_TABLES)) {
		dir.create(FOLDER_TABLES);
	}

	qualitythreshold = 0.3;

	## If the result file for automatic annotation of eyebrow raise has already been created, just load it
	if (fileExists(FOLDER_TABLES, "smad.txt")) {
		cat("Automatic annotations of smile previously done...\n");
		sbmad <- loadInternalDataFrame(FOLDER_TABLES, "sbmad.txt");
	} else {
		if (fileExists(FOLDER_TABLES, "dfwcsb.txt")) {
			df1Dqf <- loadInternalDataFrame(FOLDER_TABLES, "df1Dqfsb.txt");
			dfwc <- loadInternalDataFrame(FOLDER_TABLES, "dfwcsb.txt");
		} else {
			cat("Create wavelet 2D coefficients for s44 function...\n");

			filter <- lpl.dev.openFaceOutputAnalysis.createFilterForWaveletAnalysis(df);
			## Load the function coefficients model	
			aufm <- lpl.R.dev.ebmad.loadFunctionModel(software, "smile_presence_au_model.txt");
			## Create the function column from AU and coefficients
			sbFunction <- lpl.R.dev.hmad.functionAUColumn(audf, aufm);
			sbFunction <- 	sbFunction - mean(sbFunction, na.rm=T);

			minimalscale <- 0.2;

			cat(paste("Minimal duration parameter for eyebrow motion detection =", minimalscale, "seconds", "\n"));
			dfwc <- lpl.R.dev.ebmotion.createWaveletTable(df, sbFunction, filter, qualitythreshold, minimalscale, 24, 0.4, 1);

			df1Dqf <- lpl.R.dev.ebmotion.create1DQualityFunction(df, filter, minimalscale, 1);

			saveInternalDataFrame(df1Dqf, FOLDER_TABLES, "df1Dqfsb.txt");
			saveInternalDataFrame(dfwc, FOLDER_TABLES, "dfwcsb.txt");
		}

		smad <- lpl.R.dev.ebmotion.automaticAnnotationEBM(dfwc, df1Dqf, qualitythreshold, "Raise");
		## Save the result in a file
		saveInternalDataFrame(smad, FOLDER_TABLES, "smad.txt");
	}

	return (smad);
}

##
## The function depending on the AU
##
## df : The data frame containing the residuals and the column for filtering
## audf : The data frame containing AU values
##
lpl.R.dev.smad.functionAUS4MAD <- function(df, audf) {

	f1 <- 0.7*audf$AU06_r + 0.1*audf$AU07_r + 0.2*audf$AU12_r;
	f2 <- 0.5*audf$AU14_r + 0.5*audf$AU25_r;
	f3 <- 0.7*audf$AU06_r + 0.2*audf$AU12_r + 0.1*audf$AU25_r;
	f4 <- 0.7*audf$AU12_r + 0.3*audf$AU14_r;
      
	result <- numeric();

	for (i in c(1:nrow(df))) {
		if (f1[i] > 1.5) {
			result[i] <- 4;
		} else {
			if (f2[i] > 2) {
				result[i] <- 3;
			} else {
				if (f3[i] > 0.5) {
					result[i] <- 2;
				} else {
					if (f4[i] > 0.3) {
						result[i] <- 1;
					} else {
						result[i] <- 0;
					}
				}
			}
		}
	}
	return (result);
}

##
## Create an html document containing the summary of the analysis
## data : The data table as a data frame
## software : The software used for extracting the video output 
## projectname : The project name
##
lpl.R.dev.ebmad.createHtmlDescription <- function(data, software, projectname) {

	projectdir <- paste("projects/", projectname, sep="");
	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	topprojectdir <- projectdir;
	projectdir <- paste(projectdir, "/", software, sep="");

	FOLDER_RESULT <- paste(projectdir, "result", sep="/"); 

	WORKING_DIRECTORY <- getwd();

	df <- NULL;
	df <- lpl.R.dev.htmldf.openHtmlTag(df);
	df <- lpl.R.dev.htmldf.openHtmlHeader(df, paste(projectname, "ebm result"));
	df <- lpl.R.dev.htmldf.closeHtmlHeader(df);

	df <- lpl.R.dev.htmldf.openHtmlBody(df);

	line <- paste("Project name:", projectname, sep=" ");
	
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("Automatic prediction of eyebrows raising and frowning actions from", software, "output:", sep=" ");
	line <- lpl.R.dev.htmldf.boldHtml(line);
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	globaldf <- paste(projectdir, "EBMAD/tables", sep="/");
	line <- paste("Table of results in",  globaldf, ": ebfmad.txt, ebrmad.txt and ebmad.txt", sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);

	globaldf <- paste(projectdir, "result", sep="/");

	line <- paste("Elan editable output format in",  globaldf, ": ebmadraf.eaf", sep=" ");
	df <- lpl.R.dev.htmldf.addLine(df, line);
	df <- lpl.R.dev.htmldf.addLine(df, "");

	line <- paste("Total duration:",  (nrow(data)/25), "seconds or", (lpl.R.dev.faceOutputAnalysis.stomnsString(nrow(data)/25)), sep=" ");

	df <- lpl.R.dev.htmldf.openHtmlTable(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Class");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Number of<BR>intervals");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Duration<BR>(in seconds)");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "Ratio of<BR>Raise");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- data
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "W");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "X");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "X");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class != "X");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "M");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "A");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "B");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "B");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "C");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "C");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "D");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "D");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
		if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "E");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "E");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "N");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "N");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "");
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "AB");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B" | data$class == "C");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "ABC");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B" | data$class == "C" | data$class == "D");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "ABCD");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	ss <- subset(data, data$class == "A" | data$class == "B" | data$class == "C" | data$class == "D" | data$class == "E");
	df <- lpl.R.dev.htmldf.openHtmlTableLine(df);
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "ABCDE");
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, nrow(ss));
	df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, sum(ss$duration));
	if (nrow(ss) == 0) {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, "NA");
	} else {
		df <- lpl.R.dev.htmldf.addCellToHtmlTableLine(df, round(nrow(subset(ss, ss$ebmotion == "Raise"))/nrow(ss), 2));
	}
	df <- lpl.R.dev.htmldf.closeHtmlTableLine(df);

	df <- lpl.R.dev.htmldf.closeHtmlTag(df);

	descriptionfile <- paste(FOLDER_RESULT, "ebmaddescription.html", sep="/");
	cat(paste("Description of the automatic annotation of eyebrows actions can be found in the file ebmaddescription.html (folder ", FOLDER_RESULT, ")", sep=""));
	cat("\n");

	write.table(df, descriptionfile, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");
}

##
## Get the name of the file (not directory) in the project directory (if all is right, it is the video file) and return it
##
## projectname : The name of the project in the directory "projects"
## extension : The extension of the video file
##
## return the video file name or null in case of problem
##
lpl.R.dev.smad.getVideoFileName<- function(projectname, extension) {

	projectdir <- paste("projects/", projectname, sep="");
	regex <- paste(".*\\.", extension, sep="");
	listfiles <- list.files(projectdir, pattern=regex);
	if (length(listfiles) < 1) {
		return (cat(paste("No video file found in project directory", projectdir, "...")));
	} else if (length(listfiles) > 1) {
		return (cat(paste("More than one candidate for the video file in the project directory", projectdir, "...")));
	}
	return (listfiles[1]);
}

##
## Return the ELAN video Mime type in function of the extension of the video file or null if not found
## extension : The extension of the video file
##
## projectname : The name of the project in the directory "projects"
##
lpl.R.dev.smad.getElanMimeType<- function(projectname, extension) {

	vfn <- lpl.R.dev.smad.getVideoFileName(projectname, extension);
	if (endsWith(vfn, ".mp4")) {
		return ("video/mp4");
	} else if (endsWith(vfn, ".avi")) {
	       return ("video/*");
	} else if (endsWith(vfn, ".mov")) {
	       return ("video/quicktime");
	} else if (endsWith(vfn, ".mpg") || endsWith(vfn, ".mpeg")) {
		return ("video/mpeg");
	}

	cat("Mime Type not found! Please consult the Elan instruction guideline (https://www.mpi.nl/corpus/html/lamus/apa.html) and instantiate manually the VIDEO_MIME_TYPE R variable...");
}	
	
##
## Return the name of the eaf output file name function of the type
##
lpl.R.dev.smad.getDefaultAnnotationFileName<- function() {
	return ("s14mad.eaf");
}

##
## Create the ELAN eaf output file format for editing the automatic annotation of the smile action (in directory
## in projects/projectname/software/result)
##
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.smad.createElanEAFFile<- function(software, projectname, mediafilename, eaffilename, videomimetype) {

	## Define the folders
	projectdir <- paste("projects/", projectname, sep="");

	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	smaddir <- paste(projectdir, "/", "smile14", sep="");
	
	## The folder for the tables
	FOLDER_TABLES <- paste(smaddir, "tables", sep="/");

	if (!fileExists(FOLDER_TABLES, "smad.txt")) {
		cat("Please first run the lpl.R.dev.smad.computeSMAD routine...");
	} else {
		smad <- loadInternalDataFrame(FOLDER_TABLES, "smad.txt");
		lpl.R.dev.smad.cawriteElanEAFFile(smad, software, projectname, mediafilename, eaffilename, videomimetype); 
	}
}

##
## Create the ELAN eaf output file format for editing the automatic annotation of the smile action (in directory
## in projects/projectname/software/result)
##
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.smad.createElanManualAnnotationEAFFile <- function(software, projectname, mediafilename, eaffilename, videomimetype, annotation_file_name) {

	## Define the folders
	projectdir <- paste("projects/", projectname, sep="");

	sm <- loadInternalDataFrame(projectdir, annotation_file_name);
	lpl.R.dev.smad.cawriteElanEAFFile(sm, software, projectname, mediafilename, eaffilename, videomimetype); 
}

##
## Create and write the ELAN eaf output file format for editing the automatic annotation of the smile action (in directory
## in projects/projectname/software/result)
##
## smad : The data frame containing the data
## software : The software used to created the output "OPENFACE" or "INTRAFACE"
## projectname : The name of the project in the directory "projects"
## mediafilename : The media file name in the directory "projects/projectname"
## eaffilename : The name of the eaf output file name in projects/projectname/software/result
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
lpl.R.dev.smad.cawriteElanEAFFile <- function(smad, software, projectname, mediafilename, eaffilename, videomimetype) {

	projectdir <- paste("projects/", projectname, sep="");
	videodir <- projectdir;

	## Convert the string software in uppercases
	software <- toupper(software);
	## The project directory depending on the software used
	projectdir <- paste(projectdir, "/", software, sep="");

	smaddir <- paste(projectdir, "/", "smile14", sep="");

	FOLDER_RESULT <- paste(smaddir, "result", sep="/"); 
	if (!file.exists(FOLDER_RESULT)) {
		dir.create(FOLDER_RESULT);
	}
	WORKING_DIRECTORY <- getwd(); 

	MEDIA_URL <- "file://"; 
	MEDIA_URL <- paste(MEDIA_URL, WORKING_DIRECTORY, sep="/");
	MEDIA_URL <- paste(MEDIA_URL, "projects", sep="/");
	MEDIA_URL <- paste(MEDIA_URL, projectname, sep="/");
	MEDIA_URL <- paste(MEDIA_URL, mediafilename, sep="/");

	RELATIVE_MEDIA_URL <- "..";
	RELATIVE_MEDIA_URL <- paste(RELATIVE_MEDIA_URL, mediafilename, sep="/");

	cat(paste("Creating the Elan eaf annotation file in", FOLDER_RESULT, "...\n"));
	df <- lpl.R.dev.smad.createElanEAFTable(smad, MEDIA_URL, RELATIVE_MEDIA_URL, videomimetype);

	eaffile <- paste(FOLDER_RESULT, eaffilename, sep="/");

	write.table(df, eaffile, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");

	cat(paste("You can now open the annotation file", paste(FOLDER_RESULT, eaffilename, sep="/"), "with the ELAN Software.\n"));
}

##
## Create a one column data frame mimicking the xml ouput of the ELAN eaf format 
##
## da : The data frame containing the automatic annotation smiles action to transform 
## mediaurl : the directory of the media file name
## relativemediaurl : the directory of the relative media file name
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
## return the xml document as a data frame
##
lpl.R.dev.smad.createElanEAFTable <- function(da, mediaurl, relativemediaurl, videomimetype) {

	df <- NULL;

	## Create the xml header
	df <- addXmlHeader(df);

	df <- openAnnotationDocument(df, "SMAD", Sys.time());

	cat("Add smiles (1-4 scale) automatic annotations...\n");
	current_annotation_index = 0;
	cat("Add tier smile14...\n");
       	tiersmotion <- createAnnotationTier(current_annotation_index, da, "smile14", "smile14");

	##current_annotation_index <- tiersmotion[[2]];
	##cat("Add tier ebmqclass...\n");
	##tierebqclass <- createAnnotationTier(current_annotation_index, da, "class", "ebmqclass");

	current_annotation_index <- tiersmotion[[2]];
	lastUsedAnnotationId <- current_annotation_index;

	cat(paste("Number of intervals :", nrow(da), "\n"));
	df <- addHeaderHeader(df, mediaurl, relativemediaurl,  lastUsedAnnotationId, videomimetype); 
	cat("Add time order...\n");
	df <- addTimeOrder(df, da);

	cat("Add tiers...\n");
	df <- rbind(df, tiersmotion[[1]]);
	##df <- rbind(df, tierebqclass[[1]]);

	df <- addConstraints(df);

	## Close the xml header
	df <- closeAnnotationDocument(df);
	colnames(df) <- NULL;
	
       return (df);	
}
