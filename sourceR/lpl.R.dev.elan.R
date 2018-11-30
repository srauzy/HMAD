## lpl.R.dev.elan.R, creator S. Rauzy, LPL, 15/02/2017

##
## Create a one column data frame mimicking the xml ouput of the ELAN eaf format 
##
## da : The data frame containing the automatic annotation eyebrows action to transform 
## mediaurl : the directory of the media file name
## relativemediaurl : the directory of the relative media file name
## videomimetype : The video mime type, "video/*" for avi, "video/mp4" for mp4, see ELAN instruction guideline
##
## return the xml document as a data frame
##
lpl.R.dev.elan.createElanEAFTable <- function(da, mediaurl, relativemediaurl, videomimetype) {

	df <- NULL;

	## Create the xml header
	df <- addXmlHeader(df);

	df <- openAnnotationDocument(df, "EBMAD", Sys.time());

	cat("Add eyebrow motion automatic annotations...\n");
	current_annotation_index = 0;
	cat("Add tier ebmotion...\n");
       	tierebmotion <- createAnnotationTier(current_annotation_index, da, "ebmotion", "ebmotion");

	current_annotation_index <- tierebmotion[[2]];
	cat("Add tier ebmqclass...\n");
	tierebqclass <- createAnnotationTier(current_annotation_index, da, "class", "ebmqclass");

	current_annotation_index <- tierebqclass[[2]];
	lastUsedAnnotationId <- current_annotation_index;

	cat(paste("Number of intervals :", nrow(da), "\n"));
	df <- addHeaderHeader(df, mediaurl, relativemediaurl,  lastUsedAnnotationId, videomimetype); 
	cat("Add time order...\n");
	df <- addTimeOrder(df, da);

	cat("Add tiers...\n");
	df <- rbind(df, tierebmotion[[1]]);
	df <- rbind(df, tierebqclass[[1]]);

	df <- addConstraints(df);

	## Close the xml header
	df <- closeAnnotationDocument(df);
	colnames(df) <- NULL;
	
       return (df);	
}

##
## Create the xml header
##
## df : The xml document as a data frame
##
## return df
##
addXmlHeader <- function(df) {

	line <- "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

openAnnotationDocument <- function(df, author, date) {

	line <- "<ANNOTATION_DOCUMENT AUTHOR=\"";
	line <- paste(line, author, sep="");
	line <- paste(line, "\" DATE=\"", sep="");
	line <- paste(line, date, sep="");
	line <- paste(line,  "\" FORMAT=\"3.0\" VERSION=\"3.0\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://www.mpi.nl/tools/elan/EAFv3.0.xsd\">", sep="");
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

closeAnnotationDocument <- function(df) {

	line <- "</ANNOTATION_DOCUMENT>";
	
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

addHeaderHeader <- function(df, mediaurl, relativemediaurl, lastUsedAnnotationId, videomimetype) {

	line <- "\t<HEADER MEDIA_FILE=\"";
       	line <- paste(line, "\" TIME_UNITS=\"milliseconds\">", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t\t<MEDIA_DESCRIPTOR MEDIA_URL=\"";
	line <- paste(line, mediaurl, sep="");
	line <- paste(line, "\" MIME_TYPE=\"", sep="");
	line <- paste(line, videomimetype, sep="");
	line <- paste(line, "\" RELATIVE_MEDIA_URL=\"", sep="");
	line <- paste(line, relativemediaurl, sep="");
	line <- paste(line, "\"/>", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t\t<PROPERTY NAME=\"lastUsedAnnotationId\">";
	line <- paste(line, lastUsedAnnotationId, sep="");
	line <- paste(line, "</PROPERTY>", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t</HEADER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	return (df);
}

addTimeOrder <- function(df, da) {

	line <- "\t<TIME_ORDER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	for (i in c(1:nrow(da))) {

		line <- "\t\t<TIME_SLOT TIME_SLOT_ID=\"ts";
		line <- paste(line, (2*(i-1)+1), sep="");
		line <- paste(line, "\" TIME_VALUE=\"", sep="");
		line <- paste(line, format(round(da$tmin[i]*1000, 0), scientific=FALSE), sep="");
		line <- paste(line, "\"/>", sep="");
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);

		line <- "\t\t<TIME_SLOT TIME_SLOT_ID=\"ts";
		line <- paste(line, (2*(i-1)+2), sep="");
		line <- paste(line, "\" TIME_VALUE=\"", sep="");
		line <- paste(line, format(round(da$tmax[i]*1000, 0), scientific=FALSE), sep="");
		line <- paste(line, "\"/>", sep="");
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);
	}
	
	line <- "\t</TIME_ORDER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	return (df);
}

addAnnotationTier <- function(df, da, columnname, tiername) {

	print(colnames(da));
	column_index <- which(colnames(da) == columnname);
	print(paste(columnname, tiername, column_index, "A"));
	value <- da[, column_index];
	value[is.na(value)] <- "";
	print("Here");

	line <- "\t<TIER DEFAULT_LOCALE=\"us\" LINGUISTIC_TYPE_REF=\"default\" TIER_ID=\"";
	line <- paste(line, tiername, sep="");
	line <- paste(line, "\">", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	for (i in c(1:nrow(da))) {

		line <- "\t\t<ANNOTATION>";
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);

		line <- "\t\t\t<ALIGNABLE_ANNOTATION ANNOTATION_ID=\"a";
		line <- paste(line, i, sep="");
		line <- paste(line, "\" TIME_SLOT_REF1=\"ts", sep="");
		line <- paste(line, (2*(i-1)+1), sep="");
		line <- paste(line, "\" TIME_SLOT_REF2=\"ts", sep="");
		line <- paste(line, (2*(i-1)+2), sep="");
		line <- paste(line, "\">", sep="");
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);

		line <- "\t\t\t\t<ANNOTATION_VALUE>";
		line <- paste(line, value[i], sep="");
		line <- paste(line, "</ANNOTATION_VALUE>", sep="");
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);

		line <- "\t\t\t</ALIGNABLE_ANNOTATION>";
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);

		line <- "\t\t</ANNOTATION>";
		dfl <- data.frame(line);
		rownames(dfl) <- NULL;
		df <- rbind(df, dfl);
	}
	
	line <- "\t</TIER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	return (df);
}

createAnnotationTier <- function(current_annotation_index, da, columnname, tiername) {

	column_index <- which(colnames(da) == columnname);
	value <- da[, column_index];

	df <- NULL;

	line <- "\t<TIER DEFAULT_LOCALE=\"us\" LINGUISTIC_TYPE_REF=\"default\" TIER_ID=\"";
	line <- paste(line, tiername, sep="");
	line <- paste(line, "\">", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	for (i in c(1:nrow(da))) {

		
		if (!is.na(value[i])) {
			line <- "\t\t<ANNOTATION>";
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df <- rbind(df, dfl);

			current_annotation_index <- current_annotation_index + 1;
			line <- "\t\t\t<ALIGNABLE_ANNOTATION ANNOTATION_ID=\"a";
			line <- paste(line, current_annotation_index, sep="");
			line <- paste(line, "\" TIME_SLOT_REF1=\"ts", sep="");
			line <- paste(line, (2*(i-1)+1), sep="");
			line <- paste(line, "\" TIME_SLOT_REF2=\"ts", sep="");
			line <- paste(line, (2*(i-1)+2), sep="");
			line <- paste(line, "\">", sep="");
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df <- rbind(df, dfl);

			line <- "\t\t\t\t<ANNOTATION_VALUE>";
			line <- paste(line, value[i], sep="");
			line <- paste(line, "</ANNOTATION_VALUE>", sep="");
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df <- rbind(df, dfl);

			line <- "\t\t\t</ALIGNABLE_ANNOTATION>";
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df <- rbind(df, dfl);

			line <- "\t\t</ANNOTATION>";
			dfl <- data.frame(line);
			rownames(dfl) <- NULL;
			df <- rbind(df, dfl);
		}
	}
	
	line <- "\t</TIER>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	tier <- vector(mode="list", length=2);
	tier[[1]] <- df;
	tier[[2]] <- current_annotation_index;

	return (tier);
}

addConstraints <- function(df) {

	line <- "\t<LINGUISTIC_TYPE GRAPHIC_REFERENCES=\"false\" LINGUISTIC_TYPE_ID=\"default\" TIME_ALIGNABLE=\"true\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
	
	line <- "\t<LOCALE COUNTRY_CODE=\"EN\" LANGUAGE_CODE=\"us\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"Time subdivision of parent annotation's time interval, no time gaps allowed within this interval\" STEREOTYPE=\"Time_Subdivision\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"Symbolic subdivision of a parent annotation. Annotations refering to the same parent are ordered\" STEREOTYPE=\"Symbolic_Subdivision\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"1-1 association with a parent annotation\" STEREOTYPE=\"Symbolic_Association\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t<CONSTRAINT DESCRIPTION=\"Time alignable annotations within the parent annotation's time interval, gaps are allowed\" STEREOTYPE=\"Included_In\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	return (df);
}
