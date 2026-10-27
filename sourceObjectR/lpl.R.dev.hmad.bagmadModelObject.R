##
## lpl.R.dev.hmad.bagmadModelObject.R, creator S. Rauzy, LPL, 12/12/2022
## The bagmad object for the model allowing blinks and gaze mouvements detection
##

##
## Create the bagmad model object and return it
## return the bagmad model object
##
lpl.R.dev.hmad.bagmadModelObject.newBagmadModelObject <- function() {

	object <- vector(mode = "list", length = 4);

	object[[1]] <- "_BAGMAD_MODEL_OBJECT";

	names(object) <- c("type", "local_model", "global_model", "blink_model");  

	return (object);
}

##
## Load and return the bagmad general object defined by the 6 specific files in the model directory 
## model_type : The model type, "inertia" or "hmm"
## bagmad_general_model_file : The file where to find the 6 specific files 
## return the general bagmad model object instantiated
##
lpl.R.dev.hmad.bagmadModelObject.loadGeneralBagmadModel <- function(bagmad_general_model_file, model_type) {

	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.newBagmadGeneralModelObject();

	if (verbose) cat("Load the parameters, tables and matrixes of the general model...\n");	
	bagmad_logit_parameters_file <- paste(bagmad_general_model_file, "bagmad_logit_parameters.txt", sep="/");
	bagmad_bins_file <- paste(bagmad_general_model_file, "bagmad_lv_bins.txt", sep="/");
	bagmad_stmp_file <- paste(bagmad_general_model_file, "bagmad_stmp.txt", sep="/");
	bagmad_state_em_file <- paste(bagmad_general_model_file, "bagmad_state_em.txt", sep="/");
	bagmad_observation_em_file <- paste(bagmad_general_model_file, "bagmad_observation_em.txt", sep="/");
	if (model_type == "hmm") {
		bagmad_alpha_file <- paste(bagmad_general_model_file, "bagmad_alpha_hmm.txt", sep="/");
	} else {
		bagmad_alpha_file <- paste(bagmad_general_model_file, "bagmad_alpha_inertia.txt", sep="/");
	}
	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.setModelType(bagmad_general_model_object, model_type);

	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadLogitParametersFile(bagmad_general_model_object, bagmad_logit_parameters_file);
	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadLVBinBoundariesFile(bagmad_general_model_object, bagmad_bins_file);
	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadStateTransitionMatrixProbabilitiesFile(bagmad_general_model_object, bagmad_stmp_file);
	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadStateEmissionMatrixProbabilitiesFile(bagmad_general_model_object, bagmad_state_em_file);
	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadObservationEmissionMatrixProbabilitiesFile(bagmad_general_model_object, bagmad_observation_em_file);
	## That is the alpha which matches the number of predicted gaze intervals with the observed one
	bagmad_general_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadAlphaParameterFile(bagmad_general_model_object, bagmad_alpha_file);

	return (bagmad_general_model_object);
}

##
## Load and compute the bagmad local object (write the file if it does not already exists)
## bagmad_local_model_directory : the directory where to find or to write the local model
## df_bagmad_input : the data frame containing the input used to calculate the local model
## verbose : Put this flag at TRUE for printing comments in the R console
## return the bagmad local model object with the field instantiated
##
lpl.R.dev.hmad.bagmadModelObject.loadOrComputeLocalBagmadModel <- function(bagmad_local_model_directory, df_bagmad_input, verbose) {

	bagmad_local_model_object <- lpl.R.dev.hmad.bagmadModelObject.newBagmadLocalModelObject();

	local_parameters_file <- paste(bagmad_local_model_directory, "gaze_angle_xy_stats.txt", sep="/");
	if (file.exists(local_parameters_file)) {
		if (verbose) cat("Load the gaze angle barycenter in", local_parameters_file, "...\n");
		bagmad_local_model_object <- lpl.R.dev.hmad.bagmadModelObject.loadLocalParameters(bagmad_local_model_object, local_parameters_file);
	} else {
		bagmad_local_model_object <- lpl.R.dev.hmad.bagmadModelObject.computeLocalParameters(bagmad_local_model_object, df_bagmad_input, verbose); 
		dir.create(bagmad_local_model_directory, recursive = TRUE);
		lpl.R.dev.hmad.bagmadModelObject.writeLocalParametersFile(bagmad_local_model_object$df_parameters, local_parameters_file, verbose); 
	}

	return (bagmad_local_model_object);
}

##
## Set the bagmad local model parameters (as a local bagmad model object) in the bagmad model
## object : the bagmad model object 
## bagmad_local_model_object : The local bagmad model object for the locutor
## return the bagmad model object with the local model object instantiated
##
lpl.R.dev.hmad.bagmadModelObject.setBagmadLocalModelObject <- function(object, bagmad_local_model_object) {

	object[[2]] <- bagmad_local_model_object;

	return (object);
}

##
## Set the bagmad general model parameters (as a general bagmad model object) in the bagmad model
## object : the bagmad model object 
## bagmad_general_model_object : The general bagmad model object 
## return the bagmad model object with the general model object instantiated
##
lpl.R.dev.hmad.bagmadModelObject.setBagmadGeneralModelObject <- function(object, bagmad_general_model_object) {

	object[[3]] <- bagmad_general_model_object;

	return (object);
}

##
## Create the bagmad local model object and return it
## return the bagmad local model object
##
lpl.R.dev.hmad.bagmadModelObject.newBagmadLocalModelObject <- function() {

	object <- vector(mode = "list", length = 6);

	object[[1]] <- "_BAGMAD_LOCAL_MODEL_OBJECT";
	## A ratio of 0.93 of the data with zero au_blink which are considered to belong to the eye contact area parameter for the x gaze angle
	object[[2]] <- 0.93;
	## A ratio of 0.78 of the data with zero au_blink which are considered to belong to the eye contact area parameter for the y gaze angle
	object[[3]] <- 0.78;
	## The x centered variable is reduced by the dispersion in gaze angle x (for the au_blink zero population)
	object[[4]] <- TRUE;
	## The y centered variable is not reduced (better result)
	object[[5]] <- FALSE;

	names(object) <- c("type", "eye_contact_area_parameter_x", "eye_contact_area_parameter_y", "x_reduced", "y_reduced", "df_parameters");  

	return (object);
}

##
## Create the bagmad general model object and return it
## return the general local model object
##
lpl.R.dev.hmad.bagmadModelObject.newBagmadGeneralModelObject <- function() {

	object <- vector(mode = "list", length = 8);

	## The type of object
	object[[1]] <- "_BAGMAD_GENERAL_MODEL_OBJECT";
	
	names(object) <- c("type", "logit_parameters", "state_transition_matrix", "bins", "state_emission_transition_matrix", "observation_emission_transition_matrix", "alpha", "model_type");

	return (object);
}

##
## Set the logit parameters of the general bagmad model 
## object : the bagmad general model object
## df_logit_parameters : The data frame containing the logit parameters
## return the general model object with the logit parameters data frame instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setLogitParameters <- function(object, df_logit_parameters) {

	object[[2]] <- df_logit_parameters;
	
	return (object);
}

##
## Load the logit parameters of the general bagmad model 
## object : the bagmad general model object
## df_logit_parameters : The data frame containing the logit parameters
## return the general model object with the logit parameters data frame instanciated
##
lpl.R.dev.hmad.bagmadModelObject.loadLogitParametersFile <- function(object, logit_parameters_file) {

	df_logit_parameters <- read.delim(paste(logit_parameters_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setLogitParameters(object, df_logit_parameters);

	return (object);
}

lpl.R.dev.hmad.bagmadModelObject.newBagmadBlinkModelObject <- function() {

	object <- vector(mode = "list", length = 4);

	object[[1]] <- "_BAGMAD_BLINK_MODEL_OBJECT";

	names(object) <- c("type", "blink_detection_model", "interpolation_parameter", "model_type");  

	return (object);
}

lpl.R.dev.hmad.bagmadModelObject.newBagmadBlinkDetectionWaveletModelObject <- function() {

	object <- vector(mode = "list", length = 7);

	object[[1]] <- "_BAGMAD_BLINK_DETECTION_WAVELET_MODEL_OBJECT";

	names(object) <- c("type", "scale_min", "scale_max", "scale_step", "reference_scale", "signal_to_noise_ratio", "interpolation_parameter");  

	return (object);
}

##
## Create the bagmad general model object and return it
## return the general local model object
##
lpl.R.dev.hmad.bagmadModelObject.newBagmadBlinkDetectionHmmModelObject <- function() {

	object <- vector(mode = "list", length = 8);

	## The type of object
	object[[1]] <- "_BAGMAD_BLINK_DETECTION_HMM_MODEL_OBJECT";
	
	names(object) <- c("type", "logit_parameters", "state_transition_matrix", "bins", "state_emission_transition_matrix", "observation_emission_transition_matrix", "alpha", "model_type");

	return (object);
}

##
## Load the blink parameters, create the blink bagmad model and instanciates it in the bagmad model object
## object : The bagmad model object
## model_type : The type of model for detecting blinks, "wavelet" or "hmm"
## blink_parameters_directory : The directory containing the blink parameters file "blink_parameters.txt"
## return the bagmad model object with the blink model instanciated
##
lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkModel <- function(object, model_type, blink_parameters_directory) {

	blink_model_object <- lpl.R.dev.hmad.bagmadModelObject.newBagmadBlinkModelObject();
	blink_model_object[[4]] <- model_type;

	if (model_type == "wavelet") {
		blink_model_object[[2]] <- lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkDetectionWaveletModel(blink_parameters_directory);
	       	blink_model_object[[3]]	<- blink_model_object$blink_detection_model$interpolation_parameter;
	} else {
		blink_model_object[[2]] <- lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkDetectionHmmModel(blink_parameters_directory);
		blink_model_object[[3]]	<- 1;
	}

	object[[4]] <- blink_model_object;

	return (object);
}

##
## Load the blink detection wavelet model and instanciates it
## blink_parameters_directory : The directory containing the blink parameters file "blink_parameters.txt"
## return the blink detection wavelet model instanciated
##
lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkDetectionWaveletModel <- function(blink_parameters_directory) {

	blink_parameters_file <- paste(blink_parameters_directory, "wavelet_blink_parameters.txt", sep="/");

	df_blink_parameters <- read.delim(paste(blink_parameters_file, sep="/"), fileEncoding="UTF-8");

	blink_wavelet_model_object <- lpl.R.dev.hmad.bagmadModelObject.newBagmadBlinkDetectionWaveletModelObject();
	blink_wavelet_model_object[[2]] <- df_blink_parameters$scale_min[1]; 
	blink_wavelet_model_object[[3]] <- df_blink_parameters$scale_max[1]; 
	blink_wavelet_model_object[[4]] <- df_blink_parameters$scale_step[1]; 
	blink_wavelet_model_object[[5]] <- df_blink_parameters$reference_scale[1]; 
	blink_wavelet_model_object[[6]] <- df_blink_parameters$signal_to_noise_ratio[1]; 
	blink_wavelet_model_object[[7]] <- df_blink_parameters$interpolation_parameter[1]; 

	return (blink_wavelet_model_object);
}

##
## Load the blink detection wavelet model and instanciates it
## blink_parameters_directory : The directory containing the blink parameters file "blink_parameters.txt"
## return the blink detection wavelet model instanciated
##
lpl.R.dev.hmad.bagmadModelObject.loadAndCreateBlinkDetectionHmmModel <- function(blink_parameters_directory) {

	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.newBagmadBlinkDetectionHmmModelObject();

	if (verbose) cat("Load the parameters, tables and matrixes of the blink hmm model...\n");	
	bagmad_logit_parameters_file <- paste(blink_parameters_directory, "hmm_blink_logit_parameters.txt", sep="/");
	bagmad_bins_file <- paste(blink_parameters_directory, "hmm_blink_lv_bins.txt", sep="/");
	bagmad_stmp_file <- paste(blink_parameters_directory, "hmm_blink_stmp.txt", sep="/");
	bagmad_state_em_file <- paste(blink_parameters_directory, "hmm_blink_state_em.txt", sep="/");
	bagmad_observation_em_file <- paste(blink_parameters_directory, "hmm_blink_observation_em.txt", sep="/");
	bagmad_alpha_file <- paste(blink_parameters_directory, "hmm_blink_alpha.txt", sep="/");

	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.loadLogitParametersFile(blink_detection_hmm_object, bagmad_logit_parameters_file);
	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.loadLVBinBoundariesFile(blink_detection_hmm_object, bagmad_bins_file);
	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.loadStateTransitionMatrixProbabilitiesFile(blink_detection_hmm_object, bagmad_stmp_file);
	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.loadStateEmissionMatrixProbabilitiesFile(blink_detection_hmm_object, bagmad_state_em_file);
	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.loadObservationEmissionMatrixProbabilitiesFile(blink_detection_hmm_object, bagmad_observation_em_file);
	blink_detection_hmm_object <- lpl.R.dev.hmad.bagmadModelObject.loadAlphaParameterFile(blink_detection_hmm_object, bagmad_alpha_file);

	blink_detection_hmm_object$model_type <- "hmm";

	return (blink_detection_hmm_object);
}

##
## Set the state transition matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## df_stmp : The data frame containing the the state transition matrix probabilities
## return the general model object with the state transition matrix probabilities data frame instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setStateTransitionMatrixProbabilities <- function(object, df_stmp) {

	object[[3]] <- df_stmp;
	
	return (object);
}

##
## Load the state transition matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## stmp_file : The file containing the matrix
## return the general model object with the state transition matrix probabilities  data frame instanciated
##
##
lpl.R.dev.hmad.bagmadModelObject.loadStateTransitionMatrixProbabilitiesFile <- function(object, stmp_file) {

	df_stmp <- read.delim(paste(stmp_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setStateTransitionMatrixProbabilities(object, df_stmp);

	return (object);
}

##
## Set the bin boundaries of the logit variable for the general bagmad model 
## object : the bagmad general model object
## df_lv_bins : The data frame containing the bin boundaries 
## return the general model object with the bin boundaries  data frame instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setLVBinBoundaries <- function(object, df_lv_bins) {

	object[[4]] <- df_lv_bins;
	
	return (object);
}

##
## Load bin boundaries of the logit variable for the general bagmad model 
## object : the bagmad general model object
## lv_bins_file : The file containing the bin boundaries 
## return the general model object with the bin boundaries  data frame instanciated
##
##
lpl.R.dev.hmad.bagmadModelObject.loadLVBinBoundariesFile <- function(object, lv_bins_file) {

	df_lv_bins <- read.delim(paste(lv_bins_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setLVBinBoundaries(object, df_lv_bins);

	return (object);
}

##
## Set the state emission transition matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## df_em : The data frame containing the the state emission transition matrix probabilities
## return the general model object with the state emission matrix probabilities data frame instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setStateEmissionMatrixProbabilities <- function(object, df_em) {

	object[[5]] <- df_em;
	
	return (object);
}

##
## Set the observation emission transition matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## df_em : The data frame containing the the observation emission transition matrix probabilities
## return the general model object with the observation emission matrix probabilities data frame instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setObservationEmissionMatrixProbabilities <- function(object, df_em) {

	object[[6]] <- df_em;
	
	return (object);
}

##
## Load the state emission matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## em_file : The file containing the matrix
## return the general model object with the state emission matrix probabilities  data frame instanciated
##
##
lpl.R.dev.hmad.bagmadModelObject.loadStateEmissionMatrixProbabilitiesFile <- function(object, em_file) {

	df_em <- read.delim(paste(em_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setStateEmissionMatrixProbabilities(object, df_em);

	return (object);
}

##
## Load the observation emission matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## em_file : The file containing the matrix
## return the general model object with the observation emission matrix probabilities  data frame instanciated
##
##
lpl.R.dev.hmad.bagmadModelObject.loadObservationEmissionMatrixProbabilitiesFile <- function(object, em_file) {

	df_em <- read.delim(paste(em_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setObservationEmissionMatrixProbabilities(object, df_em);

	return (object);
}

##
## Load the observation emission matrix probabilities for the general bagmad model 
## object : the bagmad general model object
## em_file : The file containing the matrix
## return the general model object with the observation emission matrix probabilities  data frame instanciated
##
##
lpl.R.dev.hmad.bagmadModelObject.loadAlphaParameterFile <- function(object, alpha_file) {

	df_alpha <- read.delim(paste(alpha_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setAlphaParameter(object, df_alpha[1, 1]);

	return (object);
}

##
## Set the alpha parameter for the "inertia" model type
## object : the bagmad general model object
## alpha : The alpha parameter
## return the general model object with the alpha parameter instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setAlphaParameter <- function(object, alpha) {

	object[[7]] <- alpha;
	
	return (object);
}

##
## Set the model type, "inertia" or "hmm"
## object : the bagmad general model object
## model_type : The model type
## return the general model object with the model_type parameter instanciated
##
lpl.R.dev.hmad.bagmadModelObject.setModelType <- function(object, model_type) {

	object[[8]] <- model_type;
	
	return (object);
}

##
## Create the bagmad local model parameters (center in x and y gaze angles of the  eye contact area)
## and dispersion around this center). The local model will be used to form transformed variables which
## will enter the logit variable definition
## object : the bagmad local model object 
## df_bagmad_input : The data frame containing the input bagmad data for the locutor
## verbose : Put this flag at TRUE for printing comments in the R console
## return the local bagmad model object with the local parameters data frame instantiated
##
lpl.R.dev.hmad.bagmadModelObject.computeLocalParameters <- function(object, df_bagmad_input, verbose) {

	object[[6]] <- lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenter(object, df_bagmad_input, verbose);

	return (object);
}

##
## Set the bagmad local model parameters (center in x and y gaze angles of the  eye contact area
## and dispersion around this center). The local model is used to form a transformed variables which
## enter the logit variable definition
## object : the bagmad local model object 
## df_local_parameters : The data frame containing the local parameters previously computed
## return the local bagmad model object with the local parameters data frame instantiated
##
lpl.R.dev.hmad.bagmadModelObject.setLocalParameters <- function(object, df_local_parameters) {

	object[[6]] <- df_local_parameters;

	return (object);
}

##
## Save the local parameters data frame in a file
## object : the bagmad local model object 
## df_local_parameters : The data frame containing the local parameters previously computed
## file : The file path to save in
## verbose : Put this flag at TRUE for printing comments in the R console
##
lpl.R.dev.hmad.bagmadModelObject.writeLocalParametersFile <- function(df_local_parameters, local_parameters_file, verbose) {

	if (verbose) cat("Gaze angle barycenter for eye contact area saved in", local_parameters_file, "...\n");
	write.table(df_local_parameters, local_parameters_file, sep="\t", quote = FALSE,  row.names=FALSE, col.names=TRUE, fileEncoding="UTF8");
}

##
## Read the local parameters data frame and load it in the model
## object : the bagmad local model object 
## file : The file path of the local parameters data frame
## return the local bagmad model object with the local parameters data frame instantiated
##
lpl.R.dev.hmad.bagmadModelObject.loadLocalParameters <- function(object, local_parameters_file) {

	df_local_parameters <- read.delim(paste(local_parameters_file, sep="/"), fileEncoding="UTF-8");
	object <- lpl.R.dev.hmad.bagmadModelObject.setLocalParameters(object, df_local_parameters);

	return (object);
}

##
## Compute the center in x and y gaze angles for the  eye contact area and non eye contact area under the assumption
## that the eye contact area is inside the non eye contact area in  x and y gaze angles plane
## object : the bagmad local model object 
## df_input : The data frame containing the input bagmad data for the locutor
## verbose : Put this flag at TRUE for printing comments in the R console
## return a data frame with the barycenter and dispersion for the two eye contact areas 
##
lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenterOld <- function(object, df_input, verbose) {

	df_current <- df_input;
	df_current$distance <- rep(0, nrow(df_current));

	n_limit <- object$eye_contact_area_parameter*nrow(df_input);

	n_step_printed_dot <- 1000;
	n_printed_dot <- round((nrow(df_current) - n_limit)/n_step_printed_dot);
	if (verbose) cat("Compute eye contact barycenter,", n_printed_dot, "dots to run: ");

	result <- NULL;

	while (nrow(df_current) > n_limit) {
		m_x <- mean(df_current$gaze_angle_x);
		m_y <- mean(df_current$gaze_angle_y);
		sd_x <- sd(df_current$gaze_angle_x);
		sd_y <- sd(df_current$gaze_angle_y);

		df_current$distance <- sqrt(((df_current$gaze_angle_x-m_x)/sd_x)^2+((df_current$gaze_angle_y-m_y)/sd_y)^2);
		df_sorted <- df_current[order(df_current$distance), ];
		df_current <- df_sorted[1:(nrow(df_sorted)-1), ];
		if (nrow(df_current)%%n_step_printed_dot == 0) if (verbose) cat(".");
	}

	df_nec <- df_input[-as.numeric(rownames(df_current)), ];
	m_x_nec <- mean(df_nec$gaze_angle_x);
	m_y_nec <- mean(df_nec$gaze_angle_y);
	sd_x_nec <- sd(df_nec$gaze_angle_x);
	sd_y_nec <- sd(df_nec$gaze_angle_y);

	if (verbose) cat("\n");

	df_result <- data.frame(m_x, m_y, sd_x, sd_y, m_x_nec, m_y_nec, sd_x_nec, sd_y_nec);
	names(df_result) <- c("m_x_ec", "m_y_ec", "sd_x_ec", "sd_y_ec", "m_x_nec", "m_y_nec", "sd_x_nec", "sd_y_nec");
	rownames(df_result) <- NULL;

	return (df_result);
}

##
## Compute the center in x and y gaze angles for the  eye contact area under the assumption
## that the eye contact area is inside the non eye contact area in  x and y gaze angles.
## The computation is done by only keeping data with au_blink measurements 
## local_model_object : the bagmad local model object 
## df_input : The data frame containing the input bagmad data for the locutor
## verbose : Put this flag at TRUE for printing comments in the R console
## return a data frame with the barycenter and dispersion for the two eye contact areas 
##
lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenter <- function(local_model_object, df_input, verbose) {

	df_selected_input <- subset(df_input, au_blink == 0);

	if (verbose) cat("Number of zero au_blink measurements =", nrow(df_selected_input), "\n");
	ratio_of_outliers_x <- 1 - local_model_object$eye_contact_area_parameter_x;
	df_x_parameters <- lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenterOnX(df_selected_input, ratio_of_outliers_x, verbose); 
	ratio_of_outliers_y <- 1 - local_model_object$eye_contact_area_parameter_y;
	df_y_parameters <- lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenterOnY(df_selected_input, ratio_of_outliers_y, verbose); 

	df_result <- data.frame(df_x_parameters, df_y_parameters);
	rownames(df_result) <- NULL;

	return (df_result);
}

##
## Compute the center in x gaze angle for the  eye contact area under the assumption that the eye contact area is inside 
## the non eye contact area in on x gaze axis
## df_input : The data frame containing the input bagmad data for the locutor
## ratio_of_outliers : The ratio of outliers to remove 
## verbose : Put this flag at TRUE for printing comments in the R console
## return a data frame with the barycenter and dispersion for the eye contact areas 
##
lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenterOnX <- function(df_input, ratio_of_outliers, verbose) {

	df_current <- df_input;
	df_current$distance <- rep(0, nrow(df_current));

	n_limit <- (1-ratio_of_outliers)*nrow(df_input);

	n_step_printed_dot <- 500;
	n_printed_dot <- round((nrow(df_current) - n_limit)/n_step_printed_dot);
	if (verbose) cat("Compute eye contact barycenter on the x-axis,", n_printed_dot, "dots to run: ");

	result <- NULL;

	while (nrow(df_current) > n_limit) {
		m_x <- mean(df_current$gaze_angle_x);
		df_current$distance <- abs(df_current$gaze_angle_x-m_x);
		df_sorted <- df_current[order(df_current$distance), ];
		df_current <- df_sorted[1:(nrow(df_sorted)-1), ];
		if (nrow(df_current)%%n_step_printed_dot == 0) if (verbose) cat(".");
	}
	sd_x <- sd(df_current$gaze_angle_x);

	df_nec <- df_input[-as.numeric(rownames(df_current)), ];

	if (verbose) cat("\n");

	df_result <- data.frame(m_x, sd_x);
	names(df_result) <- c("m_x_ec", "sd_x_ec");
	rownames(df_result) <- NULL;

	return (df_result);
}

##
## Compute the center in y gaze angle for the  eye contact area under the assumption that the eye contact area is inside 
## the non eye contact area in on x gaze axis
## df_input : The data frame containing the input bagmad data for the locutor
## ratio_of_outliers : The ratio of outliers to remove 
## verbose : Put this flag at TRUE for printing comments in the R console
## return a data frame with the barycenter and dispersion for the eye contact areas 
##
lpl.R.dev.hmad.bagmadModelObject.computeEyeContactBarycenterOnY <- function(df_input, ratio_of_outliers, verbose) {

	df_current <- df_input;
	df_current$distance <- rep(0, nrow(df_current));

	n_limit <- (1 - ratio_of_outliers)*nrow(df_input);

	n_step_printed_dot <- 500;
	n_printed_dot <- round((nrow(df_current) - n_limit)/n_step_printed_dot);
	if (verbose) cat("Compute eye contact barycenter on the y-axis,", n_printed_dot, "dots to run: ");

	result <- NULL;

	while (nrow(df_current) > n_limit) {
		m_y <- mean(df_current$gaze_angle_y);
		df_current$distance <- abs(df_current$gaze_angle_y-m_y);
		df_sorted <- df_current[order(df_current$distance), ];
		df_current <- df_sorted[1:(nrow(df_sorted)-1), ];
		if (nrow(df_current)%%n_step_printed_dot == 0) if (verbose) cat(".");
	}
	sd_y <- sd(df_current$gaze_angle_y);

	df_nec <- df_input[-as.numeric(rownames(df_current)), ];

	if (verbose) cat("\n");

	df_result <- data.frame(m_y, sd_y);
	names(df_result) <- c("m_y_ec", "sd_y_ec");
	rownames(df_result) <- NULL;

	return (df_result);
}

##
## Set the bagmad model to the object
## object : The bagmad object
## bagmad_model_object : The bagmad model object
## return the bagmad object completed with bagmad model object
##
lpl.R.dev.hmad.bagmadModelObject.setBagmadModelObject <- function(object, bagmad_model_object) {

	object[[2]] <- bagmad_model_object;

	return (object);
}

##
## Set the HMAD input to the object
## object : The bagmad object
## hmad_object : The HMAD object containing the data 
## return the bagmad object completed with the HMAD object
##
lpl.R.dev.hmad.bagmadModelObject.setHMADObject <- function(object, hmad_object) {

	object[[3]] <- hmad_object;

	return (object);
}
