##
## lpl.R.dev.elan.elanObject.R, creator S. Rauzy, LPL, 12/12/2022
## The functions defining an elan object (input and output)
##

##
## Create the elan object and return it
## return the elan object
##
lpl.R.dev.elan.elanObject.newElanObject <- function() {

	object <- vector(mode = "list", length = 8);

	object[[1]] <- "_ELAN_OBJECT";

	object[[8]] <- list();

	names(object) <- c("type", "media_file", "media_file_name", "mime_type", "media_url", "relative_media_url", "author", "list_of_tiers");  

	return (object);
}

##
## Set the author name of the Elan document
## object : The elan object 
## author_name : The author name
## return the object completed with author field
##
lpl.R.dev.hmad.openfaceObject.setAuthorName <- function(object, author_name) {

	object[[7]] <- author_name;

	return (object);
}

##
## Set the video media file
## object : The elan object 
## media_file : The media file
## return the object completed with media file field
##
lpl.R.dev.hmad.openfaceObject.setMediaFile <- function(object, media_file) {

	## Instanciate the file variable
	object[[2]] <- media_file;

	## Instanciate the file name variable
	object[[3]] <- basename(media_file);

	## Instanciate the mime type variable
	object[[4]] <- lpl.R.dev.elan.getElanMimeType(media_file);

	## Instanciate the media url variable
	object[[5]] <- lpl.R.dev.elan.getElanMediaUrl(media_file);

	## Instanciate the relative media url variable
	object[[6]] <- lpl.R.dev.elan.getElanRelativeMediaUrl(media_file);

	return (object);
}

##
## Return the ELAN video Mime type in function of the extension of the video file or null if not found
## media_file : The video file
## return the ELAN video Mime type of the video file
##
lpl.R.dev.elan.getElanMimeType <- function(media_file) {

	vfn <- media_file;
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
## Return the ELAN media url
## media_file : The video file
## return the ELAN media url
##
lpl.R.dev.elan.getElanMediaUrl <- function(media_file) {

	return (paste("file://", media_file, sep="/"));
}

##
## Return the ELAN relative media url
## media_file : The video file
## return the ELAN relative media url
##
lpl.R.dev.elan.getElanRelativeMediaUrl <- function(media_file) {

	return (paste("..", basename(media_file), sep="/"));
}

##
## Create the elan tier object and return it
## return the elan tier object
##
lpl.R.dev.elan.elanObject.newElanTierObject <- function() {

	object <- vector(mode = "list", length = 3);

	object[[1]] <- "_ELAN_TIER_OBJECT";

	names(object) <- c("type", "tier_name", "df_intervals");  

	return (object);
}

##
## Set the name of the tier name
## object : The elan tier object 
## tier_name : The tier name
## return the object completed with tier name field
##
lpl.R.dev.elan.elanObject.setTierName <- function(object, tier_name) {

	object[[2]] <- tier_name;

	return (object);
}

##
## Set the data frame of time intervals, 3 columns (tmin, tmax, label)
## object : The elan tier object 
## df_intervals : the data frame of time intervals
## return the object completed with data frame of intervals
##
lpl.R.dev.elan.elanObject.setIntervalDataFrame <- function(object, df_intervals) {

	object[[3]] <- df_intervals;

	return (object);
}

##
## Add a tier object to the elan object
## object : The elan object 
## elan_tier_object : The elan tier object 
## return the object completed with elan tier added
##
lpl.R.dev.elan.elanObject.addTier <- function(object, elan_tier_object) {

	names <- names(object[[8]]);
	object[[8]] <- c(object[[8]], list(elan_tier_object));
	names(object[[8]]) <- c(names, elan_tier_object$tier_name);
	return (object);
}

##
## Create and write the ELAN eaf output file format
## elan_object : The elan object
## eaf_file : The eaf file where to write the result
##
lpl.R.elan.elanObject.writeElanEAFFile <- function(df_xml, eaf_file, verbose) {

	write.table(df_xml, eaf_file, sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8");

	if (verbose) cat(paste("You can now open the annotation file", eaf_file, "with the ELAN Software.\n"));
}

##
## Create a one column data frame mimicking the xml ouput of the ELAN eaf format 
## object : The elan object
## return the xml document as a data frame
##
lpl.R.elan.elanObject.createElanXmlOutput <- function(object, verbose) {

	## Create the global data frame of all tiers with a column melting tmin and tmax
	df_global <- lpl.R.elan.elanObject.createGlobalDataFrame(object, verbose);
	## The number of annotated intervals is half the size of the table
	lastUsedAnnotationId <- nrow(df_global)/2;
	
	## The Xml table
	df_xml <- NULL;

	## Create the Xml header
	df_xml <- lpl.R.elan.elanObject.addXmlHeader(df_xml);

	## Add the open document xml tag
	df_xml <- lpl.R.elan.elanObject.openAnnotationDocument(df_xml, object$author, Sys.time());
	## Add the header xml element
	df_xml <- lpl.R.elan.elanObject.addHeaderHeader(df_xml, object$media_url, object$relative_media_url,  lastUsedAnnotationId, object$video_mime_type);

	## Add the time reference xml element
	if (verbose) cat("Add time order...\n");
	df_xml <- lpl.R.elan.elanObject.addTimeOrder(df_xml, df_global);

	## Add the tiers
	list_of_tiers <- object$list_of_tiers;
	for (i in c(1:length(list_of_tiers))) {
		if (verbose) cat("Add the", list_of_tiers[[i]]$tier_name, "tier ...\n");
		df_xml <- lpl.R.elan.elanObject.addAnnotationTier(df_xml, df_global, list_of_tiers[[i]]$tier_name);
	}

	if (verbose) cat(paste("Number of total annotations :", lastUsedAnnotationId, "\n"));

	## Add the constraints
	df_xml <- lpl.R.elan.elanObject.addConstraints(df_xml);

	## Close the xml header
	df_xml <- lpl.R.elan.elanObject.closeAnnotationDocument(df_xml);
	colnames(df_xml) <- NULL;
	
       return (df_xml);	
}

##
## Add the contraints xml element to the xml data frame and return it
## df_xml : The xml data frame
## return the xml data frame with the element added
##
lpl.R.elan.elanObject.addConstraints <- function(df_xml) {

	line <- "\t<LINGUISTIC_TYPE GRAPHIC_REFERENCES=\"false\" LINGUISTIC_TYPE_ID=\"default\" TIME_ALIGNABLE=\"true\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);
	
	line <- "\t<LOCALE COUNTRY_CODE=\"EN\" LANGUAGE_CODE=\"us\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"Time subdivision of parent annotation's time interval, no time gaps allowed within this interval\" STEREOTYPE=\"Time_Subdivision\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"Symbolic subdivision of a parent annotation. Annotations refering to the same parent are ordered\" STEREOTYPE=\"Symbolic_Subdivision\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"1-1 association with a parent annotation\" STEREOTYPE=\"Symbolic_Association\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"Time alignable annotations within the parent annotation's time interval, gaps are allowed\" STEREOTYPE=\"Included_In\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	return (df_xml);
}

##
## Add the header xml element to the xml data frame and return it
## df_xml : The xml data frame
## mediaurl : The full media url
## relativemediaurl : The relative media url
## lastUsedAnnotationId : The number of annotation intervals in all the tiers
## videomimetype : The video MIME type of the media file
## return the xml data frame with the element added
##
lpl.R.elan.elanObject.addHeaderHeader <- function(df_xml, mediaurl, relativemediaurl, lastUsedAnnotationId, videomimetype) {

	line <- "\t<HEADER MEDIA_FILE=\"";
       	line <- paste(line, "\" TIME_UNITS=\"milliseconds\">", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t\t<MEDIA_DESCRIPTOR MEDIA_URL=\"";
	line <- paste(line, mediaurl, sep="");
	line <- paste(line, "\" MIME_TYPE=\"", sep="");
	line <- paste(line, videomimetype, sep="");
	line <- paste(line, "\" RELATIVE_MEDIA_URL=\"", sep="");
	line <- paste(line, relativemediaurl, sep="");
	line <- paste(line, "\"/>", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t\t<PROPERTY NAME=\"lastUsedAnnotationId\">";
	line <- paste(line, lastUsedAnnotationId, sep="");
	line <- paste(line, "</PROPERTY>", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	line <- "\t</HEADER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	return (df_xml);
}

##
## Add the open document tag in the xml element data frame and return it
## df_xml : The xml data frame
## author : The author of the annotations
## date : The date
## return the xml data frame with the open tag element added
##
lpl.R.elan.elanObject.openAnnotationDocument <- function(df_xml, author, date) {

	line <- "<ANNOTATION_DOCUMENT AUTHOR=\"";
	line <- paste(line, author, sep="");
	line <- paste(line, "\" DATE=\"", sep="");
	line <- paste(line, date, sep="");
	line <- paste(line,  "\" FORMAT=\"3.0\" VERSION=\"3.0\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://www.mpi.nl/tools/elan/EAFv3.0.xsd\">", sep="");
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);
}

##
## Add the close document tag in the xml element data frame and return it
## df_xml : The xml data frame
## return the xml data frame with the close tag element added
##
lpl.R.elan.elanObject.closeAnnotationDocument <- function(df_xml) {

	line <- "</ANNOTATION_DOCUMENT>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;

	df_xml <- rbind(df_xml, dfl);
	
	return (df_xml);
}

##
## Create the xml header
## df_xml : The xml document as a data frame
## return the xml data frame with the header element added
##
lpl.R.elan.elanObject.addXmlHeader <- function(df_xml) {

	line <- "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	return (df_xml);
}

##
## Add the time order xml element to the xml data frame and return it
## df_xml : The xml document as a data frame
## df_global : The global data frame containing all tiers and time information
## return the xml data frame with the header element added
##
lpl.R.elan.elanObject.addTimeOrder <- function(df_xml, df_global) {

	line <- "\t<TIME_ORDER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	## The list of pairs (time slot encoding, time in milliseconds)
	for (i in c(1:nrow(df_global))) {
		line <- "\t\t<TIME_SLOT TIME_SLOT_ID=\"";
		line <- paste(line, df_global$time_slot_encoding[i], sep="");
		line <- paste(line, "\" TIME_VALUE=\"", sep="");
		line <- paste(line, format(round(as.numeric(df_global$time[i])*1000, 0), scientific=FALSE), sep="");
		line <- paste(line, "\"/>", sep="");
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df_xml <- rbind(df_xml, dfl);
	}
	
	line <- "\t</TIME_ORDER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	return (df_xml);
}

##
## Create an return the global data frame containing all tiers and time information
## object : The elan object containing the list of tiers
## return the the global data frame containing all tiers and time information
##
lpl.R.elan.elanObject.createGlobalDataFrame <- function(object, verbose) {

	list_of_tiers <- object$list_of_tiers;
	df_tier_all <- NULL;
	last_annotation_count <- 0;

	if (verbose) cat("Create the global data frame...\n");
	for (i in c(1:length(list_of_tiers))) {
		
		df_tier <- list_of_tiers[[i]]$df_intervals;
		tier_name <- rep(list_of_tiers[[i]]$tier_name, nrow(df_tier));
		annotation_encoding <- c(1:nrow(df_tier));
		annotation_encoding <- annotation_encoding + last_annotation_count;
		annotation_encoding <- paste("a", annotation_encoding, sep="");
		df_tier_all <- rbind(df_tier_all, data.frame(df_tier, tier_name, annotation_encoding));

		last_annotation_count <- nrow(df_tier_all);
	}
	boundary <- rep("tmin", nrow(df_tier_all));
	time <- df_tier_all$tmin;
	df_tier_all_tmin <- data.frame(df_tier_all, boundary, time);
	boundary <- rep("tmax", nrow(df_tier_all));
	time <- df_tier_all$tmax;
	df_tier_all_tmax <- data.frame(df_tier_all, boundary, time);
	df_tier_all_2 <- rbind(df_tier_all_tmin, df_tier_all_tmax);

	df_all_sorted <- df_tier_all_2[order(as.numeric(df_tier_all_2$time)), ];
	time_slot_encoding <- c(1:nrow(df_all_sorted));
	time_slot_encoding <- paste("ts", time_slot_encoding, sep="");
	df_all_sorted <- data.frame(df_all_sorted, time_slot_encoding);
	
	return (df_all_sorted);	
}

##
## Add the elan annotation element for the given tier to the xml data frame and return it
## df_xml : The xml document as a data frame
## df_global : The global data frame containing all tiers and time information
## tier_name : The name of the tier to add
## return the xml data frame with elan annotation tier element added
##
lpl.R.elan.elanObject.addAnnotationTier <- function(df_xml, df_global, tier_name) {

	s_tier_name <- tier_name;
	df_tier <- subset(df_global, tier_name == s_tier_name);
       	
	line <- "\t<TIER DEFAULT_LOCALE=\"us\" LINGUISTIC_TYPE_REF=\"default\" TIER_ID=\"";
	line <- paste(line, tier_name, sep="");
	line <- paste(line, "\">", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	annotation_level <- levels(factor(df_tier$annotation_encoding));

	for (i in c(1:length(annotation_level))) {

		df_interval <- subset(df_tier, annotation_encoding == annotation_level[i]);
		if (!is.na(df_interval$label[1])) {
			line <- "\t\t<ANNOTATION>";
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df_xml <- rbind(df_xml, dfl);

			line <- "\t\t\t<ALIGNABLE_ANNOTATION ANNOTATION_ID=\"";
			line <- paste(line, annotation_level[i], sep="");
			line <- paste(line, "\" TIME_SLOT_REF1=\"", sep="");
			line <- paste(line, df_interval$time_slot_encoding[1], sep="");
			line <- paste(line, "\" TIME_SLOT_REF2=\"", sep="");
			line <- paste(line, df_interval$time_slot_encoding[2], sep="");
			line <- paste(line, "\">", sep="");
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df_xml <- rbind(df_xml, dfl);

			line <- "\t\t\t\t<ANNOTATION_VALUE>";
			line <- paste(line, df_interval$label[1], sep="");
			line <- paste(line, "</ANNOTATION_VALUE>", sep="");
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df_xml <- rbind(df_xml, dfl);

			line <- "\t\t\t</ALIGNABLE_ANNOTATION>";
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df_xml <- rbind(df_xml, dfl);

			line <- "\t\t</ANNOTATION>";
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df_xml <- rbind(df_xml, dfl);
		}
	}
	
	line <- "\t</TIER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df_xml <- rbind(df_xml, dfl);

	return (df_xml);
}
