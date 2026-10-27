##
## lpl.R.dev.hmad.bagmadObject.R, creator S. Rauzy, LPL, 12/12/2022
## The bagmad object for blinks and gaze mouvements detection
## Required sources : lpl.R.dev.hmad.bagmadModelObject.R, lpl.R.dev.topHatAnalysis.R, lpl.R.dev.timeSeries.R

##
## Create the bagmad object and return it
## return the bagmad object
##
lpl.R.dev.hmad.bagmadObject.newBagmadObject <- function() {

	object <- vector(mode = "list", length = 6);

	object[[1]] <- "_BAGMAD_OBJECT";

	names(object) <- c("type", "model_object", "df_input", "df_blink_solution", "df_gaze_ouput", "df_gaze_solution");  

	return (object);
}

##
## Create the default bagmad object and return it (and write the various files created by bagmad)
## hmad_directory : The HMAD directory containing the "projects" directory 
## project_name : The project name
## hmad_object : The hmad object containing the imput data
## media_file : The video file to associate with the elan eaf file
## verbose : put this flag to TRUE if you want message in the console
## return the bagmad object completed with the 2 blink and gaze annotation tables and the elan object 
##
lpl.R.dev.hmad.bagmadObject.createDefaultBagmadObject <- function(hmad_directory, project_name, hmad_object, media_file, verbose) {

	## create a bagmad object
	bagmad_object <- lpl.R.dev.hmad.bagmadObject.newBagmadObject();
	## Create the bagmad input (a selection of HMAD output of interest)
	project_directory <- paste(hmad_directory, "projects", project_name, sep="/");
	bagmad_object <- lpl.R.dev.hmad.bagmadObject.loadAndCreateDefaultBagmadModel(hmad_object, hmad_directory, project_directory, verbose);
	bagmad_object <- lpl.R.dev.hmad.bagmadObject.computeBlinkSolution(bagmad_object, verbose);
	bagmad_object <- lpl.R.dev.hmad.bagmadObject.computeGazeSolution(bagmad_object, TRUE, verbose);
	## Write the blink and gaze tables
	lpl.R.dev.hmad.bagmadObject.writeDefaultBagmadOutput(bagmad_object, project_directory, verbose);
	bagmad_result_directory <- paste(project_directory, "OPENFACE", "bagmad", "result", sep="/");
	if (!file.exists(bagmad_result_directory)) {
		dir.create(bagmad_result_directory, recursive=TRUE);
	}
	eaf_file <- paste(bagmad_result_directory, "bagmad.eaf", sep="/");
	## Create the Elan onject
	lpl.R.dev.hmad.bagmadObject.createBagmadElanOutput(bagmad_object, media_file, eaf_file, verbose); 

	return (bagmad_object);
}

##
## Set the bagmad model to the object
## object : The bagmad object
## bagmad_model_object : The bagmad model object
## return the bagmad object completed with bagmad model object
##
lpl.R.dev.hmad.bagmadObject.setBagmadModelObject <- function(object, bagmad_model_object) {

	object[[2]] <- bagmad_model_object;

	return (object);
}

##
## Form the input for the bagmad functions from the full HMAD data
## object : The bagmad object
## hmad_object : The full hmad object
## verbose : put this flag to TRUE if you want message in the console
## return the bagmad object completed with the df_input data frame
##
lpl.R.dev.hmad.bagmadObject.createInput <- function(object, hmad_object, verbose) {

	df_landmark <- hmad_object$df_landmark;
	df_au <- hmad_object$df_au;
	df_gaze <- hmad_object$df_gaze;

	df_input <- data.frame(df_landmark$time, df_landmark$confidence, df_au$AU45_r, df_gaze$gaze_angle_x, df_gaze$gaze_angle_y);
	names(df_input) <- c("time", "confidence", "au_blink", "gaze_angle_x", "gaze_angle_y");

	if (verbose) cat("Create the simplified input for the bagmad object", "...\n");

	object[[3]] <- df_input;

	return (object);
}

##
## Add to the input (the features of interest measured at each frame of the video) the transformed variables
## for the gaze angle in x and y with the center and dispersion given by the local bagmad model
## object : The bagmad object
## return the bagmad object with the df_input data frame modified in order to include transformed variables
##
lpl.R.dev.hmad.bagmadObject.addTransformedGazeAngleVariablesToInput <- function(object) {

	## The barycenter of gaze angle x and dispersion for the eye contact area
	m_x_ec <- object$model_object$local_model$df_parameter$m_x_ec[1];
	sd_x_ec <- object$model_object$local_model$df_parameter$sd_x_ec[1];
	x_reduced <- object$model_object$local_model$x_reduced;

	## The barycenter of gaze angle y and dispersion for the eye contact area
	m_y_ec <- object$model_object$local_model$df_parameter$m_y_ec[1];
	sd_y_ec <- object$model_object$local_model$df_parameter$sd_y_ec[1];
	y_reduced <- object$model_object$local_model$y_reduced;

	## Create the transformed variables
	fn_x <- lpl.R.dev.hmad.bagmadObject.computeTransformedVariable(object$df_input$gaze_angle_x, m_x_ec, sd_x_ec, x_reduced);
	fn_y <- lpl.R.dev.hmad.bagmadObject.computeTransformedVariable(object$df_input$gaze_angle_y, m_y_ec, sd_y_ec, y_reduced);

	## Add them to the data frame and instantiate it again
	df_input <- data.frame(object$df_input, fn_x, fn_y);
	object[[3]] <- df_input;

	return (object);
}

##
## Create the blink and gaze changes annotations file with the Elan format and save it
# bagmad_object : The bagmad object
## media_file : The video file
## eaf_file : he eaf file where to write the result
## verbose : put this flag to TRUE if you want message in the console
##
lpl.R.dev.hmad.bagmadObject.createBagmadElanOutput <- function(bagmad_object, media_file, eaf_file, verbose) {

	if (verbose) cat("Create the Elan eaf output...\n");
	## Create a new elan objet
	elan_object <- lpl.R.dev.elan.elanObject.newElanObject();

	## Set the media file
	elan_object <- lpl.R.dev.hmad.openfaceObject.setMediaFile(elan_object, media_file);
	## Set the author (herein automatic annotation by BAGMAD)
	elan_object <- lpl.R.dev.hmad.openfaceObject.setAuthorName(elan_object, "BAGMAD");

	## Create the blink tier object and add it to the elan object
	eto_blink <- lpl.R.dev.elan.elanObject.newElanTierObject();
	eto_blink <- lpl.R.dev.elan.elanObject.setTierName(eto_blink, "blink");	
	eto_df_blink <- data.frame(bagmad_object$df_blink_solution$tmin, bagmad_object$df_blink_solution$tmax, rep("blink", nrow(bagmad_object$df_blink_solution)));
	names(eto_df_blink) <- c("tmin", "tmax", "label");
	eto_blink <- lpl.R.dev.elan.elanObject.setIntervalDataFrame(eto_blink, eto_df_blink);
	elan_object <- lpl.R.dev.elan.elanObject.addTier(elan_object, eto_blink);

	## Create the gaze tier object and add it to the elan object
	eto_gaze <- lpl.R.dev.elan.elanObject.newElanTierObject();
	eto_gaze <- lpl.R.dev.elan.elanObject.setTierName(eto_gaze, "gaze");
	eto_df_gaze <- data.frame(bagmad_object$df_gaze_solution$tmin, bagmad_object$df_gaze_solution$tmax, bagmad_object$df_gaze_solution$prediction);
	names(eto_df_gaze) <- c("tmin", "tmax", "label");
	eto_gaze <- lpl.R.dev.elan.elanObject.setIntervalDataFrame(eto_gaze, eto_df_gaze);
	elan_object <- lpl.R.dev.elan.elanObject.addTier(elan_object, eto_gaze);

	## Create the XML output corresponding to the eaf format for the two tiers
	df_xml <- lpl.R.elan.elanObject.createElanXmlOutput(elan_object, verbose);

	## Write the eaf XML file 
	lpl.R.elan.elanObject.writeElanEAFFile(df_xml, eaf_file, verbose);

	## Return the elan object
	return (elan_object);	
}

##
## Write the blink and gaze tables in the default directory of the project (project_directory/OPENFACE/bagmad/tables
# bagmad_object : The bagmad object
## project_directory : The directory in projects containing the HMAD files for the project
## verbose : put this flag to TRUE if you want message in the console
##
lpl.R.dev.hmad.bagmadObject.writeDefaultBagmadOutput <- function(bagmad_object, project_directory, verbose) {

	bagmad_tables_directory <- paste(project_directory, "OPENFACE", "bagmad", "tables", sep="/");
	dir.create(bagmad_tables_directory, recursive = TRUE);

	blink_file <-  paste(bagmad_tables_directory, "bmad.txt", sep="/");
	if (verbose) cat("Predicted blinks table saved in", blink_file, "...\n");
	write.table(bagmad_object$df_blink_solution, blink_file, sep="\t", quote = FALSE,  row.names=FALSE, col.names=TRUE, fileEncoding="UTF8");

	gaze_file <-  paste(bagmad_tables_directory, "gmad.txt", sep="/");
	if (verbose) cat("Predicted gaze changes table saved in", gaze_file, "...\n");
	write.table(bagmad_object$df_gaze_solution, gaze_file, sep="\t", quote = FALSE,  row.names=FALSE, col.names=TRUE, fileEncoding="UTF8");
}

##
## Load the general bagmad model and create the local model proper to the proect data
# hmad_object : The hmad object
## hmad_directory : The root HMAD directory containing the models and projects directories
## project_directory : The directory in projects containing the HMAD files for the project
## verbose : put this flag to TRUE if you want message in the console
## return the bagmad object completed with the general and local models instanciated as well as the df_input data frame
##
lpl.R.dev.hmad.bagmadObject.loadAndCreateDefaultBagmadModel <- function(hmad_object, hmad_directory, project_directory, verbose) {

	general_model_directory <- paste(hmad_directory, "models", "OPENFACE", "bagmad", sep="/");
	local_model_directory <- paste(project_directory, "OPENFACE", "model", "bagmad", sep="/");

	object <- lpl.R.dev.hmad.bagmadObject.newBagmadObject();

	object <- lpl.R.dev.hmad.bagmadObject.createInput(object, hmad_object, verbose);
	bagmad_local_model_object <-  lpl.R.dev.hmad.bagmadModelObject.loadOrComputeLocalBagmadModel(local_model_directory, object$df_input, verbose);
        bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadGeneralBagmadModel(general_model_directory, "hmm");

	bagmad_model_object <- lpl.R.dev.hmad.bagmadModelObject.newBagmadModelObject();
	bagmad_model_object <- lpl.R.dev.hmad.bagmadModelObject.setBagmadLocalModelObject(bagmad_model_object, bagmad_local_model_object);
	bagmad_model_object <- lpl.R.dev.hmad.bagmadModelObject.setBagmadGeneralModelObject(bagmad_model_object, bagmad_general_model_object);

	##bagmad_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkModel(bagmad_model_object, "wavelet", general_model_directory); 
	bagmad_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkModel(bagmad_model_object, "hmm", general_model_directory); 

	# Add it to the bagmad object
	object <- lpl.R.dev.hmad.bagmadObject.setBagmadModelObject(object, bagmad_model_object);

	return (object);
}

##
## Add to the input (the features of interest measured at each frame of the video) the logit variable
## object : The bagmad object
## return the bagmad object with the df_input data frame modified in order to include the logit variable
##
lpl.R.dev.hmad.bagmadObject.addLogitVariableToInput <- function(object) {

	df_global_logit_parameters <- object$model_object$global_model$logit_parameters;

	## Create the logit variable
	lv <- lpl.R.dev.hmad.bagmadObject.computeLogitVariable(object$df_input, df_global_logit_parameters);

	## Add the logit varaible to the data frame and instantiate it again
	df_input <- data.frame(object$df_input, lv);
	object[[3]] <- df_input;

	return (object);
}

##
## Create and return the transformed variable (centered, reduced and absolute value)
## variable : The column variable to transform
## center : The center to remove
## dispersion : The dispersion divider
## reduced : flag to put to TRUE if one want to divide by the dispersion
## return the the ntransformed variable (centered, reduced and absolute value)
##
lpl.R.dev.hmad.bagmadObject.computeTransformedVariable <- function(variable, center, dispersion, reduced) {

	transformed_variable <- abs(variable - center);
	if (reduced) transformed_variable <- transformed_variable/dispersion;

	return (transformed_variable);
}

##
## Compute the logit variable for all the rows of the data frame 
## df_data : The data frame with the predictor features involved in the logit formula
## df_global_logit_parameters : The data frame containing the best fit parameters for the logit model
## return a column with the logit variable
##
lpl.R.dev.hmad.bagmadObject.computeLogitVariable <- function(df_data, df_global_logit_parameters) {

	logit_variable <- numeric(nrow(df_data));

	index_column_intercept <- which(names(df_global_logit_parameters) == "Intercept");

	list_of_predictors <- names(df_global_logit_parameters);
	list_of_predictors <- list_of_predictors[-index_column_intercept];

	list_of_predictor_column_indexes <- numeric(length(list_of_predictors));
	list_of_parameter_column_indexes <- numeric(length(list_of_predictors));
	for (j in c(1:length(list_of_predictors))) {
		list_of_predictor_column_indexes[j] <- which(names(df_data) == list_of_predictors[j]);
		list_of_parameter_column_indexes[j] <- which(names(df_global_logit_parameters) == list_of_predictors[j]);
	}

	for (i in c(1:nrow(df_data))) {
		logit_variable[i] <- df_global_logit_parameters[1, index_column_intercept];
	       	for (j in c(1:length(list_of_predictor_column_indexes))) {
			logit_variable[i] <- logit_variable[i] + df_global_logit_parameters[1, list_of_parameter_column_indexes[j]]*df_data[i, list_of_predictor_column_indexes[j]];
		}
	}

	return (logit_variable);
}

##
## Compute the solution
# object : The bagmad object
## hmad_directory : The root HMAD directory containing the models and projects directories
## project_directory : The directory in projects containing the HMAD files for the project
## blink_correction : If TRUE, the au_blink variable is corrected from physiological blink
## verbose : put this flag to TRUE if you want message in the console
## return the bagmad object completed with the general and local models instanciated as well as the df_input data frame
##
lpl.R.dev.hmad.bagmadObject.computeGazeSolution <- function(object, blink_correction, verbose) {

	## Add the transformed angle variables and the logit variable
	if (verbose) cat("Compute the transformed angle variables...\n");
	object <- lpl.R.dev.hmad.bagmadObject.addTransformedGazeAngleVariablesToInput(object);

	## Compute and add the blink correction column to the df_input of the object
	if (verbose) cat("Remove the physiological blinks from the AU blink variable...\n");
	df_input <- object$df_input;
	df_blink_solution <- object$df_blink_solution;
	column_time_index <- which(names(df_input) ==  "time");
	column_value_index <- which(names(df_input) ==  "au_blink");
	au_blink_corrected <- object$au_blink
	if (blink_correction) {
		au_blink_corrected <- lpl.R.dev.hmad.bagmadObject.interpolateSignalWithinBlinkIntervals(df_input, column_time_index, column_value_index, df_blink_solution); 
	}
	object$df_input$au_blink_corrected <- au_blink_corrected;
	## Compute the logit variable corresponding to the logit formula and the data frame df_input
	object <- lpl.R.dev.hmad.bagmadObject.addLogitVariableToInput(object);

	## Compute the bin indexes of the logarithm variables
       	df_lv_bins <- object$model_object$global_model$bins;
	f_bin_indexes <- lpl.R.dev.hmm.getBinIndexes(df_lv_bins, object$df_input$lv);

	## The output solution for each frame, clone the input
	df_output <- object$df_input;
	## The first index and last index
	itmin <- 1;
	itmax <- nrow(df_output);
	if (verbose) cat("Compute the best Viterbi solution...\n");
	bs <- lpl.R.dev.hmad.bagmadObject.bestViterbiSolution(object$model_object$global_model, f_bin_indexes, itmin, itmax, verbose);
       	df_output$prediction <- ifelse(bs == 1, "eye contact", "no eye contact");
	object[[5]] <- df_output;

	times <- df_output$time;
	predictions <- df_output$prediction;
	df_solution <- lpl.R.dev.timeSeries.createTimeAnnotationIntervals(times, predictions, "prediction");
	object[[6]] <- df_solution;

	return (object);
}

##
## From the times and the values create an return the corresponding time intervals
## times : The timeline
## values : The values
## value_column_name : The column header name for the created data frame
## return the corresponding time intervals as a data frame with tmin, tmax and the values column
##
lpl.R.dev.hmad.bagmadObject.createTimeAnnotationIntervals <- function(times, values, value_column_name) {

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

##
## bindf : The bins data frame defining the bins for the logit variable
## stmpdf : The state transition matrix probability  
## alpha : The coefficient alpha of the model for the model type "inertia"
## sempdf : The emission transition matrix of probability
## model : The type of model "inertia" vs "hmm"
## itmin : The time index min to start with
## itmax : The time index max to end with
## verbose Put this flag to true to have comments on console
## return the Viterbi solution as a column with the best state index for each time 
##
lpl.R.dev.hmad.bagmadObject.bestViterbiSolution <- function(global_model, f_bin_indexes, itmin, itmax, verbose) {
	
	model_type <- global_model$model_type;
	alpha <- global_model$alpha;

	if (model_type == "hmm") {
		sempdf <- global_model$observation_emission_transition_matrix;
		stmpdf <- matrix.power(data.matrix(global_model$state_transition_matrix), alpha);
	} else {
		sempdf <- global_model$state_emission_transition_matrix;
		stmpdf <- global_model$state_transition_matrix;
	}

	## The number of states
	ns <- ncol(stmpdf);
	## The number of frames
	n <- itmax - itmin + 1;

	NOT_ACTIVE <- 1;
	NOT_YET_VISITED <- -2;

	max_logproba <- matrix(NOT_ACTIVE, nrow=ns, ncol=n);
	max_index <- matrix(NOT_YET_VISITED, nrow=ns, ncol=n);

	## We start with a uniform distribution of states (index 1)
	## j indexes the time as the index of the current frame 
	j <- itmin;
	jloop <- 1;
	## Loop on source
	nsa <- 0;
	for (k in c(1:ns)) {
		## Loop on target
		for (i in c(1:ns)) {
			## The probability to transit from source state k to target state i with observation value svdf$bi[j] at position j
		        ## weighted by the uniform distribution
			pt <- lpl.R.dev.hmad.bagmadObject.getProbabilityOfTransition(stmpdf, k, i, alpha, sempdf, f_bin_indexes[j], model_type)/ns;
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
					pt <- lpl.R.dev.hmad.bagmadObject.getProbabilityOfTransition(stmpdf, k, i, alpha, sempdf, f_bin_indexes[j], model_type);
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
## Return the target state probability of transition
## isource : The source state index
## itarget : The target state index
## stmpdf : The data frame containing the state transition matrix probability 
## sempdf : The data frame containing the state emission matrix probability 
## ibin : The bin index for the variable (if -1, it returns the value ignoring the smile variable value)
## model_type : The model, "hmm" if HMM formula for the probability or "inertia" for the second model
## return the target state probability of transition 
##
lpl.R.dev.hmad.bagmadObject.getProbabilityOfTransition <- function(stmpdf, isource, itarget, alpha, sempdf, ibin, model_type) {

	if (model_type == "hmm") {
		p <- stmpdf[isource, itarget]*sempdf[itarget, ibin];
	} else {
		p <- alpha*stmpdf[isource, itarget] + (1-alpha)*sempdf[itarget, ibin];
	}
	return (p);
}

##
## Create and return the physiological blinks table
# object : The bagmad object
## verbose : Put this flag at TRUE for printing comments in the R console
## return the bagmad object with the physiological blinks table field instanciated
##
lpl.R.dev.hmad.bagmadObject.computeBlinkSolution <- function(object, verbose) {

	df_input <- object$df_input;
	blink_header_column <- "au_blink";

	blink_model <- object$model_object$blink_model;
	blink_model_type <- blink_model$model_type;
	
	if (verbose) cat("Detect physiological blinks...\n");
	
	if (blink_model_type == "wavelet") {
		reference_scale <- blink_model$blink_detection_model$reference_scale;
		signal_to_noise_ratio <- blink_model$blink_detection_model$signal_to_noise_ratio; 
		df_blink_solution <- lpl.R.dev.hmad.bagmadObject.createPhysiologicalBlinkTable(df_input, blink_header_column, reference_scale, signal_to_noise_ratio, verbose);
	} else {
		df_blink_solution <- lpl.R.dev.hmad.bagmadObject.createPhysiologicalBlinkTableFromHmmModel(df_input, blink_header_column, blink_model$blink_detection_model, verbose);
	}

	object[[4]] <- 	df_blink_solution;
	if (verbose) cat(nrow(df_blink_solution), "physiological blinks in the file\n");
	if (verbose) cat("Mean blink rate per minute =", round(60*nrow(df_blink_solution)/df_input$time[nrow(df_input)], 2), "\n");

	return (object);
}

##
## Create and return the physiological blinks table
## df : The data frame with the blink variable
## blink_header_column : The name of the column with the blink variable
## reference_scale : The scale of reference for computing the noise in the wavelet space
## signal_to_noise_ratio : The signal to noise ratio (noise measured at the reference scale) for selecting physiological blink events
## verbose : Put this flag at TRUE for printing comments in the R console
## return the physiological blinks table
##
lpl.R.dev.hmad.bagmadObject.createPhysiologicalBlinkTableFromHmmModel <- function(df_input, blink_header_column, blink_detection_model, verbose) {

	df_blink_logit_parameters <- blink_detection_model$logit_parameters;
	df_bins <-  blink_detection_model$bins;

	## Create the blink logit variable
	lv <- lpl.R.dev.hmad.bagmadObject.computeLogitVariable(df_input, df_blink_logit_parameters);
	## Create the bin indexes 	
	f_bin_indexes <- lpl.R.dev.hmm.getBinIndexes(df_bins, lv);
	itmin <- 1;
	itmax <- nrow(df_input);

	bs <- lpl.R.dev.hmad.bagmadObject.bestViterbiSolution(blink_detection_model, f_bin_indexes, itmin, itmax, verbose);
	time <- df_input$time;
	df_output <- data.frame(time, bs);
	df_output$prediction <- ifelse(bs == 2, "blink", "no blink");
	df_prediction_intervals <- lpl.R.dev.timeSeries.createTimeAnnotationIntervals(df_output$time, df_output$prediction, "blink");
	df_prediction_intervals$duration <- as.numeric(df_prediction_intervals$tmax) - as.numeric(df_prediction_intervals$tmin);

	df_solution <- subset(df_prediction_intervals, blink == "blink" & duration <= 0.4);

	return (df_solution);
}

##
## Create and return the physiological blinks table
## df : The data frame with the blink variable
## blink_header_column : The name of the column with the blink variable
## reference_scale : The scale of reference for computing the noise in the wavelet space
## signal_to_noise_ratio : The signal to noise ratio (noise measured at the reference scale) for selecting physiological blink events
## verbose : Put this flag at TRUE for printing comments in the R console
## return the physiological blinks table
##
lpl.R.dev.hmad.bagmadObject.createPhysiologicalBlinkTable <- function(df, blink_header_column, reference_scale, signal_to_noise_ratio, verbose) {

	## Define the limits of the scale parameters for the Top Hat wavelet analysis 
	smin <- 0.08
	smax <- 0.4; 
	step <- 0.08;

	## Create the scale grid parameters object for the  wavelet analysis 
	scale_grid_parameters <- lpl.R.dev.topHatAnalysis.createScaleGridParameters(smin, smax, step);

	## The variables
	time <- df$time;
	blink_column_index <- which(names(df) == blink_header_column);
	blink <- df[ , blink_column_index];

	## Compute the Top Hat wavelet coefficients
	dwf <- lpl.R.dev.topHatAnalysis.computeTopHatWaveletTransform(time, blink, scale_grid_parameters, verbose);

	## The coefficients at the given reference scale
	dwf_reference_scale <- subset(dwf, duration == reference_scale);

	## The dispersion of the wavelet coefficients for the reference scale
	sd_reference_scale <- sd(dwf_reference_scale$wcoef);

	## Keep only the intervals with wavelets coefficients greater than the signal to noise ration
	dwf_selected <- subset(dwf, wcoef > signal_to_noise_ratio*sd_reference_scale);

	## Remove the overlapping intervals
	column_value_index <- which(names(dwf_selected) ==  "wcoef");
	df_blink <- lpl.R.dev.timeSeries.computeMaximalTimeIntervals(dwf_selected, column_value_index);

	return (df_blink);
}

##
## Interpolate the signal within (and in the vicinity) of physiological blink intervals and return the interpolated signal
## df : The data frame with the signal
## column_time_index : the column index of the time variable
## column_value_index : the column index of the signal values
## df_blink : The data frame containing the tmin and tmax columns for the blink intervals
## return the interpolated signal
##
lpl.R.dev.hmad.bagmadObject.interpolateSignalWithinBlinkIntervals <- function(df, column_time_index, column_value_index, df_blink) {

	return (lpl.R.dev.timeSeries.interpolateValues(df, column_time_index, column_value_index, df_blink));
}
