##
## lpl.R.dev.hmad.openfaceObject.R
## The functions defining an openface object (input and output)
##

##
## Create the openface object and return it
## return the openface object
##
lpl.R.dev.hmad.openfaceObject.newOpenfaceObject <- function() {

	object <- vector(mode = "list", length = 3);

	object[[1]] <- "_OPENFACE_OBJECT";

	names(object) <- c("type", "openface_csv_file", "df_openface");  

	return (object);
}

##
## Set the OpenFace output csv file
## openface_csv_file : The OpenFace output file
## return the object completed with OpenFace file field
##
lpl.R.dev.hmad.openfaceObject.setOpenFaceCsvFile <- function(object, openface_csv_file) {

	## Instanciate the file variable
	object[[2]] <- openface_csv_file;

	return (object);
}

##
## Read the OpenFace csv file and put in the the field df_openface
## openface_csv_file : The OpenFace output file
## verbose : Put this flag at TRUE for printing comments in the R console
## return the object completed with df_openface field
##
lpl.R.dev.hmad.openfaceObject.readCsvFile <- function(object, verbose) {

	## Read the cvs OpenFace output file
	if (verbose) cat("Read the OpenFace output file in", object$openface_csv_file, "...\n");
	df_openface <- read.table(object$openface_csv_file, h=TRUE, sep=",");
	
	object[[3]] <- df_openface;

	return (object);
}
