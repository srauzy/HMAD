
createMergedEvaluationTable <- function(software, projectname) {

	return (createMergedEvaluationTable1(software, projectname, "eval"));
}

createMergedEvaluationTable1 <- function(software, projectname, evalprojectname) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdirroot <- paste("projects/", projectname, sep="");
	projectdirectory <- paste(projectdirroot, "/", software, sep="");

	PROJECT_DIRECTORY_TABLES <- paste(projectdirectory, "tables", sep="/");
	if (is.null(evalprojectname)) {
		FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/"); 
	} else {
		FOLDER_EVAL <- paste(projectdirectory, evalprojectname, sep="/"); 
	}
	
	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	fileannotation <- paste(projectdirroot, "ebmotion.txt", sep="/");

	if (!file.exists(fileannotation)) {
		tf <- getTableFiles(projectdirroot);
	        if (is.null(tf)) {
			ebmad <- read.table(paste(PROJECT_DIRECTORY_TABLES, "ebmad.txt", sep="/"), h=TRUE);
			eval <- globalEvaluationWithoutGoldStandard(ebmad);
			saveInternalDataFrame(eval, FOLDER_EVAL, "eval.txt");
			return (cat("No annotation table for the project..."));
		}	
		cat(paste("Create annotation data frame from file :", tf[1], "\n"));
		ebat <- loadAnnotationTable(projectdirroot, tf[1]);
		saveInternalDataFrame(ebat, projectdirroot, "ebmotion.txt");
	} else {
		ebat <- loadInternalDataFrame(projectdirroot, "ebmotion.txt");
		if (nrow(ebat) == 0) {
			ebmad <- read.table(paste(PROJECT_DIRECTORY_TABLES, "ebmad.txt", sep="/"), h=TRUE);
			eval <- globalEvaluationWithoutGoldStandard(ebmad);
			saveInternalDataFrame(eval, FOLDER_EVAL, "eval.txt");
			return (cat("No annotation table for the project..."));
		}
	}
	ebatn <- normalizeEBMlabel(ebat);

	if (nrow(ebatn) != 0) {
	splitfile <- paste(FOLDER_EVAL, "sebmad.txt", sep="/");
	
	if (!file.exists(splitfile)) {
		cat("Split the intervals to match evaluation areas...\n");
		
		ebmad <- read.table(paste(PROJECT_DIRECTORY_TABLES, "ebmad.txt", sep="/"), h=TRUE);
		
		sebmad <- splitResultIntervalForMultipleSolution(ebmad, ebatn);
		saveInternalDataFrame(sebmad, FOLDER_EVAL, "sebmad.txt");
	} else {
		sebmad <- loadInternalDataFrame(FOLDER_EVAL, "sebmad.txt");
	}

	cat("Create evalutation table...\n");
	eval <- globalEvaluation(sebmad, ebatn);

	saveInternalDataFrame(eval, FOLDER_EVAL, "eval.txt");

	evaltable <- evaluationMergedTable(eval);

	saveInternalDataFrame(evaltable, FOLDER_EVAL, "evaltable.txt");
	cat(paste("You can now open the evaluation table file", paste(FOLDER_EVAL, "evaltable.txt", sep="/"), ".\n"));
	}	

	return (evaltable);
}

biasAnalysis <- function(projectname, ebmtype, wDparameter) {

	projectdirectory <- paste("projects", projectname, sep="/");
	PROJECT_DIRECTORY_TABLES <- paste(projectdirectory, "tables", sep="/");
	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");

	if (ebmtype == "Raise") {
		filename <- "evalr.txt";
	} else if (ebmtype == "Frown") {
		filename <- "evalf.txt";
	}

	d <- loadInternalDataFrame(FOLDER_EVAL, filename);

	itdc <- loadInternalDataFrame(PROJECT_DIRECTORY_TABLES, "itdcomplete.txt");
	nt <- nrow(itdc);

	myaw <- numeric(nrow(d));
	mroll <- numeric(nrow(d));
	mpitch <- numeric(nrow(d));
	wyaw <- numeric(nrow(d));
	wroll <- numeric(nrow(d));
	wpitch <- numeric(nrow(d));

	wS <- numeric(nrow(d));
	mS <- numeric(nrow(d));

	meaX <- numeric(nrow(d));
	meaY <- numeric(nrow(d));
	weaX <- numeric(nrow(d));
	weaY <- numeric(nrow(d));

	column_index_S = which(colnames(itdc) == "S");
	column_index_pitch = which(colnames(itdc) == "pitch");
	column_index_yaw = which(colnames(itdc) == "yaw");
	column_index_roll = which(colnames(itdc) == "roll");
	column_index_eaX = which(colnames(itdc) == "eaX");
	column_index_eaY = which(colnames(itdc) == "eaY");
	
	type <- logical(nrow(d));
	 	

	for (i in c(1:nrow(d))) {

		if (d$igs[i]==-1) {
			type[i] <- F;
		} else {
			type[i] <- T;
		}

		if (d$class[i] == "X" | d$class[i] == "N") {
			wS[i] <- NA;
			mS[i] <- NA;
			myaw[i] <- NA;
			mpitch[i] <- NA;
			mroll[i] <- NA;
			wyaw[i] <- NA;
			wpitch[i] <- NA;
			wroll[i] <- NA;
			weaX[i] <- NA;
			weaY[i] <- NA;
			meaX[i] <- NA;
			meaY[i] <- NA;
		} else {

			sigma <- d$duration[i]/2;
			index_min <- round(25*d$tmin[i]+1, 0);
			index_max <- round(25*d$tmax[i]+1, 0);

			mS[i] <- mean(itdc[index_min:index_max, column_index_S], na.rm = T);
			myaw[i] <- mean(itdc[index_min:index_max, column_index_yaw], na.rm = T);
			mpitch[i] <- mean(itdc[index_min:index_max, column_index_pitch], na.rm = T);
			mroll[i] <- mean(itdc[index_min:index_max, column_index_roll], na.rm = T);
			meaX[i] <- mean(itdc[index_min:index_max, column_index_eaX], na.rm = T);
			meaY[i] <- mean(itdc[index_min:index_max, column_index_eaY], na.rm = T);

			jposition <- index_min + floor((index_max-index_min)/2);
			jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
			jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

			wc <- waveletCoefAtTime(itdc$time, itdc[, column_index_S], sigma, jposition, jmin, jmax, wDparameter);
			wS[i] <- wc[1];
			wc <- waveletCoefAtTime(itdc$time, itdc[, column_index_pitch], sigma, jposition, jmin, jmax, wDparameter);
			wpitch[i] <- wc[1];
			wc <- waveletCoefAtTime(itdc$time, itdc[, column_index_yaw], sigma, jposition, jmin, jmax, wDparameter);
			wyaw[i] <- wc[1];
			wc <- waveletCoefAtTime(itdc$time, itdc[, column_index_roll], sigma, jposition, jmin, jmax, wDparameter);
			wroll[i] <- wc[1];
			wc <- waveletCoefAtTime(itdc$time, itdc[, column_index_eaX], sigma, jposition, jmin, jmax, wDparameter);
			weaX[i] <- wc[1];
			wc <- waveletCoefAtTime(itdc$time, itdc[, column_index_eaY], sigma, jposition, jmin, jmax, wDparameter);
			weaY[i] <- wc[1];
		}
	}

	d <- data.frame(d, type, wS, mS, wpitch, mpitch, wyaw, myaw, wroll, mroll, weaX, meaX, weaY, meaY);

	return (d);
}

calibrationAnalysis <- function(projectname, type, qt, typeres) {

	projectdirectory <- paste("projects", projectname, sep="/");
	PROJECT_DIRECTORY_TABLES <- paste(projectdirectory, "tables", sep="/");
	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/"); 
	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	if (type == "Raise") {
		tarname <-  "tarr.txt";
	} else {
		
		tarname <-  "tarf.txt";
	}
	filename <- paste(FOLDER_EVAL, tarname, sep="/");

	if (!file.exists(filename)) {
	
		filename <- paste(FOLDER_EVAL, "itdcan.txt", sep="/");
	
		if (!file.exists(filename)) {
	
			df <- read.table(paste(PROJECT_DIRECTORY_TABLES, "itdcomplete.txt", sep="/"), h=TRUE);

			fileannotation <- paste(FOLDER_EVAL, "ebmotion.txt", sep="/");

			if (!file.exists(fileannotation)) {
				tf <- getTableFiles(projectdirectory);
			       if (is.null(tf)) {
					return (cat("No annotation table for the project..."));
				}	
				cat(paste("Create annotation data frame from file :", tf[1], "\n"));
				ebat <- loadAnnotationTable(projectdirectory, tf[1]);
				saveInternalDataFrame(ebat, FOLDER_EVAL, "ebmotion.txt");
			} else {
				ebat <- loadInternalDataFrame(FOLDER_EVAL, "ebmotion.txt");
			}

			dfa <- addAnnotationToIntrafaceOutput(df, ebat);
			saveInternalDataFrame(dfa, FOLDER_EVAL, "itdcan.txt");
		} else {
			dfa <- loadInternalDataFrame(FOLDER_EVAL, "itdcan.txt");
		}

		if (type == "Raise") {
			w2Ds <- read.table(paste(PROJECT_DIRECTORY_TABLES, "w2Dsr.txt", sep="/"), h=TRUE);
		} else if (type == "Frown") {
			w2Ds <- read.table(paste(PROJECT_DIRECTORY_TABLES, "w2Dsf.txt", sep="/"), h=TRUE);
		}

		itdc <- read.table(paste(PROJECT_DIRECTORY_TABLES, "itdcomplete.txt", sep="/"), h=TRUE);

		##d <- createTableOfAnnotatedIntervalsWithQuality(dfa, w2Ds, type);	

		d <- createTableOfAnnotatedIntervalsForResiduals(dfa, w2Ds, itdc, type);

		saveInternalDataFrame(d,  FOLDER_EVAL, tarname);
	} else {
		d <- loadInternalDataFrame(FOLDER_EVAL, tarname);
	}

	dc <- subset(d, d$quality > qt);
	cat(paste("Number of movements =", nrow(dc)));

	dfres <- NULL;

	kc = 1;
	for (k in c(1:12)) {
		if (k != 5) {
			if (typeres == "raw") {
				namex <- paste(paste("wrP", getStringNumber(k), sep=""), "x", sep="");
				namey <- paste(paste("wrP", getStringNumber(k), sep=""), "y", sep="");
			} else {
				namex <- paste(paste("wrcP", getStringNumber(k), sep=""), "x", sep="");
				namey <- paste(paste("wrcP", getStringNumber(k), sep=""), "y", sep="");
			}
			column_index_x = which(colnames(dc) == namex);
			column_index_y = which(colnames(dc) == namey);
	
			rx <- dc[, column_index_x];
			ry <- dc[, column_index_y];

			dff <- data.frame(namex, mean(rx), sd(rx), (mean(rx)/sd(rx)));
			names(dff) <- c("value", "mean", "sd", "zcore");
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);

			dff <- data.frame(namey, mean(ry), sd(ry), (mean(ry)/sd(ry)));
			names(dff) <- c("value", "mean", "sd", "zcore");
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);

			kc <- kc + 1;
		}

	}
	return(dfres);
}

createMeltedEvaluationTable<- function(software, projectname, meltprojectlist) {
	return (createMeltedEvaluationTable1(software, projectname, meltprojectlist, "eval"))
}
createMeltedEvaluationTable1<- function(software, projectname, meltprojectlist, evaldir) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdirroot <- paste("projects/", projectname, sep="");
	if (!file.exists(projectdirroot)) {
		dir.create(projectdirroot);
	}
	projectdirectory <- paste(projectdirroot, "/", software, sep="");
	if (!file.exists(projectdirectory)) {
		dir.create(projectdirectory);
	}

       	FOLDER_EVAL <- paste(projectdirectory, evaldir, sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	for (i in c(1:length(meltprojectlist))) {
		subprojectname <- meltprojectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		subprojectdirectory <- paste(subprojectdirectory, "/", software, sep="");
		
		subfoldereval <- paste(subprojectdirectory, evaldir, sep="/");
		print(subfoldereval);

       		if (file.exists(paste(subfoldereval, "eval.txt", sep="/"))) {
			
       			eval <- loadInternalDataFrame(subfoldereval, "eval.txt");
			eval$subproject <- subprojectname;
			eval$project <- projectname;
			cat(paste("Melt evaluation of project :", subprojectdirectory, "\n")); 
			evaltot <- rbind(evaltot, eval);
		}
	}
	rownames(evaltot) <- NULL;
	saveInternalDataFrame(evaltot, FOLDER_EVAL, "eval.txt");
	print(nrow(evaltot));
	evaltable <- evaluationMergedTable(evaltot);

	saveInternalDataFrame(evaltable, FOLDER_EVAL, "evaltable.txt");
	cat(paste("You can now open the evaluation table file", paste(FOLDER_EVAL, "evaltable.txt", sep="/"), ".\n"));	

}

createMeltedAnnotatedTable<- function(projectname, meltprojectlist) {

	projectdirectory <- paste("projects", projectname, sep="/");
       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	for (i in c(1:length(meltprojectlist))) {
		subprojectname <- meltprojectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
       		if (file.exists(paste(subfoldereval, "itdcan.txt", sep="/"))) {
       			eval <- loadInternalDataFrame(subfoldereval, "itdcan.txt");
			eval$subproject <- subprojectname;
			eval$project <- projectname;
			cat(paste("Melt annotation of project :", subprojectdirectory, "\n")); 
			evaltot <- rbind(evaltot, eval);
		}
	}
	rownames(evaltot) <- NULL;
	saveInternalDataFrame(evaltot, FOLDER_EVAL, "itdcan.txt");

	return (evaltot);
}

createMeltedAnnotatedTableNormalized <- function(projectname, meltprojectlist) {

	projectdirectory <- paste("projects", projectname, sep="/");
       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	for (i in c(1:length(meltprojectlist))) {
		subprojectname <- meltprojectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
       		if (file.exists(paste(subfoldereval, "itdcan.txt", sep="/"))) {
       			eval <- loadInternalDataFrame(subfoldereval, "itdcan.txt");
			## Create a filter for 3-sigma outliers on the S factor 
			cis = which(colnames(eval) == "S");
			fs <- createFilterOutliersForColumn(eval, cis, 3);
			## Create a filter for 3-sigma outliers on the projection amplitude error variable
			cieay = which(colnames(eval) == "eaY");
			feay <- createFilterOutliersForColumn(eval, cieay, 3);
			## Merge the two filters
			f <- mergeFilter(fs, feay);
			sdrx <- createSDResidual12WithFilter(eval, f, "x");
			sdry <- createSDResidual12WithFilter(eval, f, "y");

			for (i in c(1:12)) {
				cat(paste(sdrx$sr[i], " ", sdry$sr[i]), "\n");
				pointname = paste("P", getStringNumber(i), sep="");
				column_index <- which(colnames(eval) == paste(paste("r", pointname, sep=""), "x", sep=""));
				cat(paste(pointname, column_index , sd(eval[, column_index], na.rm=T), "\n"));
				eval[, column_index] <- eval[, column_index]/sdrx$sr[i];
				cat(paste(sd(eval[, column_index], na.rm=T), "\n"));
				column_index <- which(colnames(eval) == paste(paste("r", pointname, sep=""), "y", sep=""));
				eval[, column_index] <- eval[, column_index]/sdry$sr[i];
			}
			eval$subproject <- subprojectname;
			eval$project <- projectname;
			cat(paste("Melt annotation of project :", subprojectdirectory, "\n")); 
			evaltot <- rbind(evaltot, eval);
		}
	}
	rownames(evaltot) <- NULL;
	saveInternalDataFrame(evaltot, FOLDER_EVAL, "itdcannorm.txt");

	return (evaltot);
}

createMeltedAnnotatedTableNormalizedOpenFace <- function(projectname, meltprojectlist) {

	projectdirectory <- paste("projects", projectname, sep="/");
       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	for (i in c(1:length(meltprojectlist))) {
		subprojectname <- meltprojectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/");
       		if (file.exists(paste(subfoldereval, "itdcan.txt", sep="/"))) {
       			eval <- loadInternalDataFrame(subfoldereval, "itdcan.txt");
			## Create a filter for 3-sigma outliers on the S factor
		        d <- subset(eval, eval$label=="Y");
			tsd <- lpl.R.dev.computeProjectedResidualsStandardDispersion(d);

			sdrx <- tsd$sdx;
			sdry <- tsd$sdy;

			for (i in c(1:68)) {
				cat(paste(sdrx[i], " ", sdry[i]), "\n");
				pointname = paste("P", getStringNumber(i), sep="");
				column_index <- which(colnames(eval) == paste(paste("r", pointname, sep=""), "x", sep=""));
				cat(paste(pointname, column_index , sd(eval[, column_index], na.rm=T), "\n"));
				eval[, column_index] <- eval[, column_index]/sdrx[i];
				cat(paste(sd(eval[, column_index], na.rm=T), "\n"));
				column_index <- which(colnames(eval) == paste(paste("r", pointname, sep=""), "y", sep=""));
				eval[, column_index] <- eval[, column_index]/sdry[i];
			}
			eval$subproject <- subprojectname;
			eval$project <- projectname;
			cat(paste("Melt annotation of project :", subprojectdirectory, "\n")); 
			evaltot <- rbind(evaltot, eval);
		}
	}
	rownames(evaltot) <- NULL;
	saveInternalDataFrame(evaltot, FOLDER_EVAL, "itdcannorm.txt");

	return (evaltot);
}

createMeltedAnnotation <- function(projectname, meltprojectlist) {

	projectdirectory <- paste("projects", projectname, sep="/");
       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	for (i in c(1:length(meltprojectlist))) {
		subprojectname <- meltprojectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
       		if (file.exists(paste(subfoldereval, "ebmotion.txt", sep="/"))) {
       			eval <- loadInternalDataFrame(subfoldereval, "ebmotion.txt");
			if (is.null(eval$duration)) {
			eval$duration <- eval$tmax - eval$tmin;	
			}
			eval$subproject <- subprojectname;
			eval$project <- projectname;
			cat(paste("Melt annotation of project :", subprojectdirectory, "\n")); 
			evaltot <- rbind(evaltot, eval);
		}
	}
	rownames(evaltot) <- NULL;
	saveInternalDataFrame(evaltot, FOLDER_EVAL, "ebmotion.txt");

	return (evaltot);
}


createProjectsEvalutationTable <- function(projectlist) {

	projectdirectory <- paste("projects","ALL", sep="/");
	if (!file.exists(projectdirectory)) {
		dir.create(projectdirectory);
	}

       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	corpus <- character(5*length(projectlist));
	OT <- numeric(5*length(projectlist));
	OF <- numeric(5*length(projectlist));
	class <- numeric(5*length(projectlist));
	tM <- numeric(length(projectlist));
	tC <- numeric(length(projectlist));
	tG <- numeric(length(projectlist));
	tO <- numeric(length(projectlist));
	OFN <- numeric(length(projectlist));
	OTN <- numeric(length(projectlist));
	PTN <- numeric(length(projectlist));
	PFN <- numeric(length(projectlist));
	errm <- numeric(length(projectlist));

	for (i in c(1:length(projectlist))) {
		subprojectname <- projectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
       		if (file.exists(paste(subfoldereval, "eval.txt", sep="/"))) {
			corpus[5*(i-1)+1] <- projectlist[i];
			print(projectlist[i]);
       			eval <- loadInternalDataFrame(subfoldereval, "eval.txt");

			##s <- subset(eval, eval$class != "N" & eval$class != "X" & eval$class != "E" & eval$class != "D");
			s <- subset(eval, eval$class != "N" & eval$class != "X");
		        tC[i] <- sum(s$duration);	
			s <- subset(eval, eval$class != "X");
		        tM[i] <- sum(s$duration);
			#"s <- subset(eval, eval$class != "N" & eval$class != "X" & eval$class != "E"  & eval$class != "D"& eval$igs == -1);
			s <- subset(eval, eval$class != "N" & eval$class != "X" & eval$igs == -1);
			tG[i] <- sum(s$duration);
			errm[i] <- mean(s$errm);
			OFN[i] <- nrow(s);
			##s <- subset(eval, eval$class != "N" & eval$class != "X" & eval$class != "E"  & eval$class != "D"& eval$igs != -1);
			s <- subset(eval, eval$class != "N" & eval$class != "X" & eval$igs != -1);
			OTN[i] <- nrow(s);
			tO[i] <- sum(s$duration);
					
		}
	}

	lm <- lm(formula = OFN ~ tG + 0);
	print(summary(lm));
	PFN <- lm$coefficients[[1]]*tG;
	PTN <- (OFN+OTN)-PFN;

	d <- data.frame(projectlist, tC, tM, tG, tO, OFN, OTN, PFN, PTN, errm);

	return (d);
}

createProjectsEvalutationTable3 <- function(projectlist) {

	projectdirectory <- paste("projects","ALL", sep="/");
	if (!file.exists(projectdirectory)) {
		dir.create(projectdirectory);
	}

       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	corpus <- character(5*length(projectlist));
	OT <- numeric(5*length(projectlist));
	OF <- numeric(5*length(projectlist));
	class <- numeric(5*length(projectlist));
	durationM <- numeric(5*length(projectlist));
	durationX <- numeric(5*length(projectlist));
	tM <- numeric(length(projectlist));
	tC <- numeric(length(projectlist));
	tG <- numeric(length(projectlist));
	tX <- numeric(length(projectlist));
	tP <- numeric(5*length(projectlist));
	OTtotal <- numeric(length(projectlist));
	errm <- numeric(length(5*length(projectlist)));
	md <- numeric(length(5*length(projectlist)));
	rraise <- numeric(length(5*length(projectlist)));
	mzs <- numeric(length(5*length(projectlist)));
	fC <- numeric(length(5*length(projectlist)));
	rraiseT <- numeric(length(5*length(projectlist)));

	for (i in c(1:length(projectlist))) {
		subprojectname <- projectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
       		if (file.exists(paste(subfoldereval, "eval.txt", sep="/"))) {
			corpus[5*(i-1)+1] <- projectlist[i];
			print(projectlist[i]);
       			eval <- loadInternalDataFrame(subfoldereval, "eval.txt");
			evaltable <- loadInternalDataFrame(subfoldereval, "evaltable.txt");
			s <- subset(eval, eval$class == "A" & eval$igs != "-1");
		        OT[5*(i-1)+1] <- nrow(s);	
			s <- subset(eval, eval$class == "A" & eval$igs == "-1");
		        OF[5*(i-1)+1] <- nrow(s);
			sr <- subset(s, s$ebmotion == "Raise");
			rraise[5*(i-1)+1] <- nrow(sr)/OF[5*(i-1)+1];
			class[5*(i-1)+1] <- "A";
			errm[5*(i-1)+1] <- mean(s$errm);
			mzs[5*(i-1)+1] <- mean(s$zscore);
			md[5*(i-1)+1] <- mean(s$duration);
			tP[5*(i-1)+1] <- sum(s$duration);
			corpus[5*(i-1)+2] <- projectlist[i]
			s <- subset(eval, eval$class == "B" & eval$igs != "-1");
		        OT[5*(i-1)+2] <- nrow(s);	
			s <- subset(eval, eval$class == "B" & eval$igs == "-1");
		        OF[5*(i-1)+2] <- nrow(s);
			sr <- subset(s, s$ebmotion == "Raise");
			rraise[5*(i-1)+2] <- nrow(sr)/OF[5*(i-1)+2];
			class[5*(i-1)+2] <- "B";
			errm[5*(i-1)+2] <- mean(s$errm);
			mzs[5*(i-1)+2] <- mean(s$zscore);
			md[5*(i-1)+2] <- mean(s$duration);
			tP[5*(i-1)+2] <- sum(s$duration);
			corpus[5*(i-1)+3] <- projectlist[i]
			s <- subset(eval, eval$class == "C" & eval$igs != "-1");
		        OT[5*(i-1)+3] <- nrow(s);	
			s <- subset(eval, eval$class == "C" & eval$igs == "-1");
		        OF[5*(i-1)+3] <- nrow(s);
			sr <- subset(s, s$ebmotion == "Raise");
			rraise[5*(i-1)+3] <- nrow(sr)/OF[5*(i-1)+3];
			class[5*(i-1)+3] <- "C";
			errm[5*(i-1)+3] <- mean(s$errm);
			mzs[5*(i-1)+3] <- mean(s$zscore);
			md[5*(i-1)+3] <- mean(s$duration);
			tP[5*(i-1)+3] <- sum(s$duration);
			corpus[5*(i-1)+4] <- projectlist[i]
			s <- subset(eval, eval$class == "D" & eval$igs != "-1");
		        OT[5*(i-1)+4] <- nrow(s);	
			s <- subset(eval, eval$class == "D" & eval$igs == "-1");
		        OF[5*(i-1)+4] <- nrow(s);
			sr <- subset(s, s$ebmotion == "Raise");
			rraise[5*(i-1)+4] <- nrow(sr)/OF[5*(i-1)+4];
			class[5*(i-1)+4] <- "D";
			errm[5*(i-1)+4] <- mean(s$errm);
			mzs[5*(i-1)+4] <- mean(s$zscore);
			md[5*(i-1)+4] <- mean(s$duration);
			tP[5*(i-1)+4] <- sum(s$duration);
			corpus[5*(i-1)+5] <- projectlist[i]
			s <- subset(eval, eval$class == "E" & eval$igs != "-1");
		        OT[5*(i-1)+5] <- nrow(s);	
			s <- subset(eval, eval$class == "E" & eval$igs == "-1");
		        OF[5*(i-1)+5] <- nrow(s);
			sr <- subset(s, s$ebmotion == "Raise");
			rraise[5*(i-1)+5] <- nrow(sr)/OF[5*(i-1)+5];
			errm[5*(i-1)+5] <- mean(s$errm);
			mzs[5*(i-1)+5] <- mean(s$zscore);
			tP[5*(i-1)+5] <- sum(s$duration);
			md[5*(i-1)+5] <- mean(s$duration);
			corpus[5*(i-1)+5] <- projectlist[i]
			class[5*(i-1)+5] <- "E";

			s <- subset(eval, eval$class != "N" & eval$class != "X");
		        tC[i] <- sum(s$duration);	
			s <- subset(eval, eval$class != "X");
		        tM[i] <- sum(s$duration);
			s <- subset(eval, eval$class != "N" & eval$class != "X" & eval$igs == -1);
			tG[i] <- sum(s$duration);
			s <- subset(eval, eval$class != "X" & eval$igs == -1);
		        tX[i] <- sum(s$duration);
			s <- subset(eval, eval$igs != -1);

			durationM[5*(i-1)+1] <- tM[i];
			durationM[5*(i-1)+2] <- tM[i];
			durationM[5*(i-1)+3] <- tM[i];
			durationM[5*(i-1)+4] <- tM[i];
			durationM[5*(i-1)+5] <- tM[i];
			durationX[5*(i-1)+1] <- tX[i];
			durationX[5*(i-1)+2] <- tX[i];
			durationX[5*(i-1)+3] <- tX[i];
			durationX[5*(i-1)+4] <- tX[i];
			durationX[5*(i-1)+5] <- tX[i];

			OTtotal[5*(i-1)+1] <- evaltable$Nobs[3];
			OTtotal[5*(i-1)+2] <- evaltable$Nobs[3];
			OTtotal[5*(i-1)+3] <- evaltable$Nobs[3];
			OTtotal[5*(i-1)+4] <- evaltable$Nobs[3];
			OTtotal[5*(i-1)+5] <- evaltable$Nobs[3];

			s <- subset(eval, eval$gsebmotion == "Raise" | eval$gsebmotion == "Frown");
			ss <- subset(s, s$gsebmotion == "Raise");
			rraiseT[5*(i-1)+1] <- nrow(ss)/nrow(s); 
			rraiseT[5*(i-1)+2] <- rraiseT[5*(i-1)+1];
			rraiseT[5*(i-1)+3] <- rraiseT[5*(i-1)+1];
			rraiseT[5*(i-1)+4] <- rraiseT[5*(i-1)+1];
			rraiseT[5*(i-1)+5] <- rraiseT[5*(i-1)+1];


			
 
 
		}
	}
	
	d <- data.frame(corpus, class, OT, OF, errm, durationM, durationX, md, rraise, rraiseT, mzs);

	frequency <- numeric(length(projectlist));
	s <- subset(d, d$class == "A");
	frequency[1] <- sum(s$OF)/sum(tM);
	s <- subset(d, d$class == "B");
	frequency[2] <- sum(s$OF)/sum(tM);
	s <- subset(d, d$class == "C");
	frequency[3] <- sum(s$OF)/sum(tM);
	s <- subset(d, d$class == "D");
	frequency[4] <- sum(s$OF)/sum(tM);
	s <- subset(d, d$class == "E");
	frequency[5] <- sum(s$OF)/sum(tM);

	

	frequencyC <- numeric(length(projectlist));
	s <- subset(d, d$class == "A");
	frequencyC[1] <- sum(s$OF)/sum(tC);
	s <- subset(d, d$class == "B");
	frequencyC[2] <- sum(s$OF)/sum(tC);
	s <- subset(d, d$class == "C");
	frequencyC[3] <- sum(s$OF)/sum(tC);
	s <- subset(d, d$class == "D");
	frequencyC[4] <- sum(s$OF)/sum(tC);
	s <- subset(d, d$class == "E");
	frequencyC[5] <- sum(s$OF)/sum(tC);

	frequencyG <- numeric(length(projectlist));
	s <- subset(d, d$class == "A");
	frequencyG[1] <- sum(s$OF)/sum(tG);
	s <- subset(d, d$class == "B");
	frequencyG[2] <- sum(s$OF)/sum(tG);
	s <- subset(d, d$class == "C");
	frequencyG[3] <- sum(s$OF)/sum(tG);
	s <- subset(d, d$class == "D");
	frequencyG[4] <- sum(s$OF)/sum(tG);
	s <- subset(d, d$class == "E");
	frequencyG[5] <- sum(s$OF)/sum(tG);

	PF <- numeric(5*length(projectlist));
	PFC <- numeric(5*length(projectlist));
	PFG <- numeric(5*length(projectlist));
	PT <- numeric(5*length(projectlist));
	PTC <- numeric(5*length(projectlist));
	PTG <- numeric(5*length(projectlist));
	CPT <- numeric(5*length(projectlist));
	CPF <- numeric(5*length(projectlist));
	COT <- numeric(5*length(projectlist));
 
	for (i in c(1:length(projectlist))) {
		for (j in c(1:5)) {
			fC[5*(i-1)+j] <- frequency[j];
			PF[5*(i-1)+j] <- frequency[j]*tM[i];
			PFC[5*(i-1)+j] <- frequencyC[j]*tC[i];
			PFG[5*(i-1)+j] <- frequencyG[j]*tG[i];
			PT[5*(i-1)+j] <- (OT[5*(i-1)+j]+OF[5*(i-1)+j]) - PF[5*(i-1)+j];
			PTC[5*(i-1)+j] <- (OT[5*(i-1)+j]+OF[5*(i-1)+j]) - PFC[5*(i-1)+j];
			PTG[5*(i-1)+j] <- (OT[5*(i-1)+j]+OF[5*(i-1)+j]) - PFG[5*(i-1)+j];
			if (j == 1) {
				CPT[5*(i-1)+j] <- PT[5*(i-1)+j];
				CPF[5*(i-1)+j] <- PF[5*(i-1)+j];
				COT[5*(i-1)+j] <- OT[5*(i-1)+j];
			} else {
				CPT[5*(i-1)+j] <- PT[5*(i-1)+j] + CPT[5*(i-1)+j-1];
				CPF[5*(i-1)+j] <- PF[5*(i-1)+j] + CPF[5*(i-1)+j-1];
				COT[5*(i-1)+j] <- OT[5*(i-1)+j] + COT[5*(i-1)+j-1];
			}
		}
	}

	d <- data.frame(d, tC, tG, tP, PF, PFC, PFG, PT, PTC, PTG, CPT, CPF, COT, fC, OTtotal);
	return (d);
}

createProjectsEvalutationTable4 <- function(projectlist) {

	projectdirectory <- paste("projects","ALL", sep="/");
	if (!file.exists(projectdirectory)) {
		dir.create(projectdirectory);
	}

       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	corpus <- character(5*length(projectlist));
	class <- character(5*length(projectlist));
	Nobs <- numeric(5*length(projectlist)); 
	cumNobs <- numeric(5*length(projectlist)); 
	Npred <- numeric(5*length(projectlist));
	precision <- numeric(5*length(projectlist));
	recall <- numeric(5*length(projectlist));
	cumrecall <- numeric(5*length(projectlist));
	ebmmatch <- numeric(5*length(projectlist));
	duration <- numeric(5*length(projectlist));
	Ntot <- numeric(5*length(projectlist));

	for (i in c(1:length(projectlist))) {
		subprojectname <- projectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
       		if (file.exists(paste(subfoldereval, "evaltable.txt", sep="/"))) {
			corpus[5*(i-1)+1] <- projectlist[i];
			print(projectlist[i]);
       			eval <- loadInternalDataFrame(subfoldereval, "evaltable.txt");
			for (j in c(1:5)) {
				corpus[5*(i-1)+j] <- projectlist[i];
				class[5*(i-1)+j] <- as.character(eval$class[3+j]);
				Nobs[5*(i-1)+j] <- eval$Nobs[3+j];
				Npred[5*(i-1)+j] <- eval$Npred[3+j];
				precision[5*(i-1)+j] <- eval$precision[3+j];
				recall[5*(i-1)+j] <- eval$recall[3+j];
				if(j == 1) {
					cumNobs[5*(i-1)+j] <- Nobs[5*(i-1)+j];
					cumrecall[5*(i-1)+j] <- recall[5*(i-1)+j];
				} else {
					cumNobs[5*(i-1)+j] <- eval$Nobs[8+j];
					cumrecall[5*(i-1)+j] <- eval$recall[8+j];
				}
				
				ebmmatch[5*(i-1)+j] <- eval$ebmmatch[3+j];
				duration[5*(i-1)+j] <- eval$duration[3+j];
				Ntot[5*(i-1)+j] <- eval$Nobs[3];
			}
 
		}
	}
	
	d <- data.frame(corpus, class, Nobs, cumNobs, Npred, Ntot, precision,recall, cumrecall, ebmmatch, duration);

	rhoobs <- Nobs/Ntot;
        cilNobs <- Nobs - sqrt(Nobs);	
        ciuNobs <- Nobs + sqrt(Nobs);
	cilcumNobs <- cumNobs - sqrt(cumNobs);
	ciucumNobs <- cumNobs + sqrt(cumNobs);
	cilcumrecall <- cilcumNobs/Ntot;
	ciucumrecall <- ciucumNobs/Ntot;
	d <- data.frame(d, rhoobs, cilNobs, ciuNobs, cilcumNobs, ciucumNobs, cilcumrecall, ciucumrecall);
	return (d);
}


createProjectsEvalutationTable1 <- function(projectlist) {

	projectdirectory <- paste("projects","ALL", sep="/");
	if (!file.exists(projectdirectory)) {
		dir.create(projectdirectory);
	}

       	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/");
       	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	evaltot <- NULL;

	corpus <- character(length(projectlist));
	AT <- numeric(length(projectlist));
	BT <- numeric(length(projectlist));
	CT <- numeric(length(projectlist));
	DT <- numeric(length(projectlist));
	ET <- numeric(length(projectlist));
	AF <- numeric(length(projectlist));
	BF <- numeric(length(projectlist));
	CF <- numeric(length(projectlist));
	DF <- numeric(length(projectlist));
	EF <- numeric(length(projectlist));
	tM <- numeric(length(projectlist));
	tN <- numeric(length(projectlist));

	for (i in c(1:length(projectlist))) {
		subprojectname <- projectlist[i];
		subprojectdirectory <- paste("projects", subprojectname, sep="/");
		subfoldereval <- paste(subprojectdirectory, "eval", sep="/"); 
		corpus[i] <- projectlist[i]
       		if (file.exists(paste(subfoldereval, "eval.txt", sep="/"))) {
			print(corpus[i]);
       			eval <- loadInternalDataFrame(subfoldereval, "eval.txt");
			s <- subset(eval, eval$class == "A" & eval$igs != "-1");
		        AT[i] <- nrow(s);	
			s <- subset(eval, eval$class == "A" & eval$igs == "-1");
		        AF[i] <- nrow(s);
			s <- subset(eval, eval$class == "B" & eval$igs != "-1");
		        BT[i] <- nrow(s);	
			s <- subset(eval, eval$class == "B" & eval$igs == "-1");
		        BF[i] <- nrow(s);	
			s <- subset(eval, eval$class == "C" & eval$igs != "-1");
		        CT[i] <- nrow(s);	
			s <- subset(eval, eval$class == "C" & eval$igs == "-1");
		        CF[i] <- nrow(s);	
			s <- subset(eval, eval$class == "D" & eval$igs != "-1");
		        DT[i] <- nrow(s);	
			s <- subset(eval, eval$class == "D" & eval$igs == "-1");
		        DF[i] <- nrow(s);	
			s <- subset(eval, eval$class == "E" & eval$igs != "-1");
		        ET[i] <- nrow(s);	
			s <- subset(eval, eval$class == "E" & eval$igs == "-1");
		        EF[i] <- nrow(s);
			s <- subset(eval, eval$class == "N");
		        tN[i] <- sum(s$duration);	
			s <- subset(eval, eval$class != "X");
		        tM[i] <- sum(s$duration);	

		}
	}
	
	d <- data.frame(corpus, AT, AF, BT, BF, CT, CF, DT, DF, ET, EF, tN, tM);
	df <- data.frame("ALL",sum(AT), sum(AF), sum(BT), sum(BF), sum(CT), sum(CF), sum(DT), sum(DF), sum(ET), sum(EF), sum(tN), sum(tM));
	rownames(df) <- NULL;
	names(df) <- names(d);
	##d <- rbind(d, df);

	df <- data.frame("frequency",sum(AT)/sum(tM), sum(AF)/sum(tM), sum(BT)/sum(tM), sum(BF)/sum(tM), sum(CT)/sum(tM), sum(CF)/sum(tM), sum(DT)/sum(tM), sum(DF)/sum(tM), sum(ET)/sum(tM), sum(EF)/sum(tM), sum(tN)/sum(tM), sum(tM)/sum(tM));
	rownames(df) <- NULL;
	names(df) <- names(d);
	##d <- rbind(d, df);

	AFP <- numeric(length(projectlist)+1);
	BFP <- numeric(length(projectlist)+1);
	CFP <- numeric(length(projectlist)+1);
	DFP <- numeric(length(projectlist)+1);
	EFP <- numeric(length(projectlist)+1);

	AFP <- tM*df[1, 3];
	BFP <- tM*df[1, 5];
	CFP <- tM*df[1, 7];
	DFP <- tM*df[1, 9];
	EFP <- tM*df[1, 11];

	d <- data.frame(d, AFP, BFP, CFP, DFP, EFP);

	return (d);
}

createAnnotatedTable <- function(projectname, overwrite) {

	projectdirectory <- paste("projects", projectname, sep="/");
	PROJECT_DIRECTORY_TABLES <- paste(projectdirectory, "tables", sep="/");
	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/"); 
	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	filename <- paste(FOLDER_EVAL, "itdcan.txt", sep="/");

	if (!file.exists(filename)  | overwrite) {
	
		df <- read.table(paste(PROJECT_DIRECTORY_TABLES, "itdcomplete.txt", sep="/"), h=TRUE);

		fileannotation <- paste(FOLDER_EVAL, "ebmotion.txt", sep="/");

		if (!file.exists(fileannotation)) {
			tf <- getTableFiles(projectdirectory);
		       	if (length(tf) == 0) {
				cat("No annotation table for the project...");
				df$ebmotion <- "Nothing";
				saveInternalDataFrame(df, FOLDER_EVAL, "itdcan.txt");
				return (NULL);
			}	
			cat(paste("Create annotation data frame from file :", tf[1], "\n"));
			ebat <- loadAnnotationTable(projectdirectory, tf[1]);
			saveInternalDataFrame(ebat, FOLDER_EVAL, "ebmotion.txt");
		} else {
			ebat <- loadInternalDataFrame(FOLDER_EVAL, "ebmotion.txt");
		}

		dfa <- addAnnotationToIntrafaceOutput(df, ebat);
		saveInternalDataFrame(dfa, FOLDER_EVAL, "itdcan.txt");
	} else {
		dfa <- loadInternalDataFrame(FOLDER_EVAL, "itdcan.txt");
	}
	return (dfa);
}

createEvaluationTable <- function(software, projectname, type) {

	## Convert the string software in uppercases
	software <- toupper(software);
	projectdirroot <- paste("projects/", projectname, sep="");
	projectdirectory <- paste(projectdirroot, "/", software, sep="");

	PROJECT_DIRECTORY_TABLES <- paste(projectdirectory, "tables", sep="/");
	FOLDER_EVAL <- paste(projectdirectory, "eval", sep="/"); 
	if (!file.exists(FOLDER_EVAL)) {
		dir.create(FOLDER_EVAL);
	}

	filename <- paste(FOLDER_EVAL, "itdcan.txt", sep="/");

	if (!file.exists(filename)) {
	
		df <- read.table(paste(PROJECT_DIRECTORY_TABLES, "itdcomplete.txt", sep="/"), h=TRUE);

		fileannotation <- paste(projectdirroot, "ebmotion.txt", sep="/");

		if (!file.exists(fileannotation)) {
			tf <- getTableFiles(projectdirroot);
		       	if (is.null(tf)) {
				return (cat("No annotation table for the project..."));
			}	
			cat(paste("Create annotation data frame from file :", tf[1], "\n"));
			ebat <- loadAnnotationTable(projectdirroot, tf[1]);
			saveInternalDataFrame(ebat, projectdirroot, "ebmotion.txt");
		} else {
			ebat <- loadInternalDataFrame(projectdirroot, "ebmotion.txt");
		}

		dfa <- addAnnotationToIntrafaceOutput(df, ebat);
		saveInternalDataFrame(dfa, FOLDER_EVAL, "itdcan.txt");
	} else {
		dfa <- loadInternalDataFrame(FOLDER_EVAL, "itdcan.txt");
	}	

	if (type == "Raise") {
		ebmad <- read.table(paste(PROJECT_DIRECTORY_TABLES, "ebrmad.txt", sep="/"), h=TRUE);
	} else if (type == "Frown") {
		ebmad <- read.table(paste(PROJECT_DIRECTORY_TABLES, "ebfmad.txt", sep="/"), h=TRUE);
	}	

	timotion <- createTableOfAnnotatedIntervals(dfa, type);
	if (is.null(timotion)) {
	    return (cat("End of the evalutation..."));
	}

	cat("Split the intervals to match evaluation areas...\n");

	if (type == "Raise") {
		splitfile <- paste(FOLDER_EVAL, "sebmadr.txt", sep="/");
	} else {
		splitfile <- paste(FOLDER_EVAL, "sebmadf.txt", sep="/");
	}

	if (!file.exists(splitfile)) {
		sebmad <- splitResultIntervalForMultipleSolution(ebmad, timotion);
		if (type == "Raise") {
			saveInternalDataFrame(sebmad, FOLDER_EVAL, "sebmadr.txt");
		} else if (type == "Frown") {
			saveInternalDataFrame(sebmad, FOLDER_EVAL, "sebmadf.txt");
		}
	} else {
		if (type == "Raise") {
			sebmad <- loadInternalDataFrame(FOLDER_EVAL, "sebmadr.txt");
		} else if (type == "Frown") {
			sebmad <- loadInternalDataFrame(FOLDER_EVAL, "sebmadf.txt");
		}
	}
	eval <- globalEvaluation(sebmad, timotion);
	dfap <- addEBMPrediction(dfa, sebmad);
	cat("Create evalutation table...\n");
	evaltable <- evaluationTable(eval, dfap);

	if (type == "Raise") {
		saveInternalDataFrame(evaltable, FOLDER_EVAL, "evaltabler.txt");
		cat(paste("You can now open the evaluation table file", paste(FOLDER_EVAL, "evaltabler.txt", sep="/"), ".\n"));	
	} else if (type == "Frown") {
		saveInternalDataFrame(evaltable, FOLDER_EVAL, "evaltablef.txt");
		cat(paste("You can now open the evaluation table file", paste(FOLDER_EVAL, "evaltablef.txt", sep="/"), ".\n"));	
	}

	if (type == "Raise") {
		saveInternalDataFrame(eval, FOLDER_EVAL, "evalr.txt");
	} else if (type == "Frown") {
		saveInternalDataFrame(eval, FOLDER_EVAL, "evalf.txt");
	}
}

##
## Get the list of Praat Table files in the directory "projecdir"
## projectdir : The project directory to be scanned
##
getTableFiles <- function(projectdir) {

	listfiles <- list.files(projectdir, pattern=".*Table");
	cat(paste("In directory", projectdir, ":\n")); 
	if (length(listfiles) < 1) {
		cat(paste("WARNING : Please drop the Praat Table file in your project folder", projectdir, "!\n"));
	        return (NULL);	
	} else if ((length(listfiles) > 1)) {
		for (i in c(1:length(listfiles))) {
			cat(paste("Praat Table file", i, ":", listfiles[i], "\n")); 
		}
	 	print("WARNING : You have more than one Praat Table file in your folder project ! Please select the index of the Praat Table file you want to work on.\n"); 
	}
	return (listfiles);
}

globalEvaluation <- function(res, gs) {

	igs <- numeric(nrow(res));
	igs <- rep(-1, nrow(res));

	r <- res;
	r$li <- c(1:nrow(r));

	g <- gs;
	g$li <- c(1:nrow(g));

	resA <- subset(r, r$class == "A");
	gsA <- g;
	if (nrow(resA) != 0) {
		dA <- evaluationAnnotation(resA, g);
		dfA <- subset(dA, dA$lligs != -1);
		if (nrow(dfA) != 0) {
			gsA <- g[-dfA$lligs, ];
			for (i in c(1:nrow(dA))) {
				igs[dA$li[i]] <- dA$ligs[i];
			}
		}
	}

	resB <- subset(r, r$class == "B");
	gsB <- gsA;
	if (nrow(resB) != 0) {
       		dB <- evaluationAnnotation(resB, gsA);
		dfB <- subset(dB, dB$lligs != -1);
		if (nrow(dfB) != 0) {
			gsB <- gsA[-dfB$lligs, ];
			for (i in c(1:nrow(dB))) {
				igs[dB$li[i]] <- dB$ligs[i];
			}
		}
	}	

	resC <- subset(r, r$class == "C");
	gsC <- gsB;
	if (nrow(resC) != 0) {
		dC <- evaluationAnnotation(resC, gsB);
		dfC <- subset(dC, dC$lligs != -1);
		if (nrow(dfC) != 0) {
			gsC <- gsB[-dfC$lligs, ];	
			for (i in c(1:nrow(dC))) {
				igs[dC$li[i]] <- dC$ligs[i];
			}
		}
	}

	resD <- subset(r, r$class == "D");
	gsD <- gsC;
	if (nrow(resD) != 0) {
       		dD <- evaluationAnnotation(resD, gsC);
		dfD <- subset(dD, dD$lligs != -1);
		if (nrow(dfD) != 0) {
			gsD <- gsC[-dfD$lligs, ];
			for (i in c(1:nrow(dD))) {
				igs[dD$li[i]] <- dD$ligs[i];
			}
		}
	}
	
	resE <- subset(r, r$class == "E");
	gsE <- gsD;
	if (nrow(resE) != 0) {
       		dE <- evaluationAnnotation(resE, gsD);
		dfE <- subset(dE, dE$lligs != -1);
		if (nrow(dfE) != 0) {
			gsE <- gsD[-dfE$lligs, ];
			for (i in c(1:nrow(dE))) {
				igs[dE$li[i]] <- dE$ligs[i];
			}
		}
	}

	resX <- subset(r, r$class == "X");
	gsX <- gsE;
	if (nrow(resX) != 0) {
       		dX <- evaluationAnnotation(resX, gsE);
		dfX <- subset(dX, dX$lligs != -1);
		if (nrow(dfX) != 0) {
			gsX <- gsE[-dfX$lligs, ];
			for (i in c(1:nrow(dX))) {
				igs[dX$li[i]] <- dX$ligs[i];
			}
		}
	}

	resN <- subset(r, r$class == "N");
	if (nrow(resN) != 0) {
       		dN <- evaluationAnnotation(resN, gsX);
		dfN <- subset(dN, dN$lligs != -1);
		if (nrow(dfN) != 0) {
			for (i in c(1:nrow(dN))) {
				igs[dN$li[i]] <- dN$ligs[i];
			}
		}
	}

	gstmin <- numeric(nrow(res));
	gstmax <- numeric(nrow(res));
	gsduration <- numeric(nrow(res));
	gsebmotion <- numeric(nrow(res));

	for (i in c(1:nrow(res))) {
		if (igs[i] == -1) {
			gstmin[i] <- NA;
			gstmax[i] <- NA;
			gsduration[i] <- NA;
			gsebmotion[i] <- NA;
		} else {
			gstmin[i] <- gs$tmin[igs[i]];
			gstmax[i] <- gs$tmax[igs[i]];
			gsduration[i] <- round(gs$tmax[igs[i]]-gs$tmin[igs[i]], 2);
			gsebmotion[i] <- gs$ebmotion[igs[i]];

			
		}
	}

	d <- data.frame(res, igs, gstmin, gstmax, gsduration, gsebmotion);

	return (d);
}

addCorrelation <- function(res, itdc, value, wDparameter) {

	wpitch <- numeric(nrow(res));
	for (i in c(1:nrow(res))) {

		itmin <- round(res$tmin[i]*25, 0);
		itmax <- round(res$tmax[i]*25, 0);

		sigma <- round(res$duration[i]/2, 2);
		jmin = round(max(1, (itmin - ((wDparameter+1)*sigma)*25)), 0);
		jmax = round(min(nrow(itdc), (itmax + ((wDparameter+1)*sigma)*25)), 0);

		jposition <- jmin + floor((jmax-jmin)/2);
		wc <- waveletCoefAtTime(itdc$time, value, sigma, jposition, jmin, jmax, wDparameter);
		wpitch[i] <- wc[1];
	}

	return (data.frame(res, wpitch));
}

globalEvaluationWithoutGoldStandard <- function(res) {

	igs <-  numeric(nrow(res));
	gstmin <- numeric(nrow(res));
	gstmax <- numeric(nrow(res));
	gsduration <- numeric(nrow(res));
	gsebmotion <- numeric(nrow(res));

	for (i in c(1:nrow(res))) {
		igs[i] <- -1;
		gstmin[i] <- NA;
		gstmax[i] <- NA;
		gsduration[i] <- NA;
		gsebmotion[i] <- NA;
	}

	d <- data.frame(res, igs, gstmin, gstmax, gsduration, gsebmotion);

	return (d);
}

evaluationTable <- function(dscore, dfwp) {

	dsobs <- subset(dscore, dscore$igs != -1);
	dfwpx <-  subset(dfwp, dfwp$ebmclass == "X");

	dfres <- NULL;
	dff <- data.frame("W", NA, nrow(dsobs),  NA, NA, nrow(dfwp)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsx <- subset(dscore, dscore$class == "X");
	dsxobs <- subset(dsx, dsx$igs != -1);
	dfwpx <-  subset(dfwp, dfwp$ebmclass == "X");

	dff <- data.frame("X", NA, nrow(dsxobs),  NA, NA, nrow(dfwpx)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsnx <- subset(dscore, dscore$class != "X");
	dsnxobs <- subset(dsnx, dsnx$igs != -1);
	dfwpnx <-  subset(dfwp, dfwp$ebmclass != "X");

	dff <- data.frame("M", NA, nrow(dsnxobs),  NA, NA, nrow(dfwpnx)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	nobstot <- nrow(dsnxobs);

	dfres <- evaluationAnnotationClass(dfres, dscore, dfwp, nrow(dsnxobs), "A");
	dfres <- evaluationAnnotationClass(dfres, dscore, dfwp, nrow(dsnxobs), "B");
	dfres <- evaluationAnnotationClass(dfres, dscore, dfwp, nrow(dsnxobs), "C");
	dfres <- evaluationAnnotationClass(dfres, dscore, dfwp, nrow(dsnxobs), "D");
	dfres <- evaluationAnnotationClass(dfres, dscore, dfwp, nrow(dsnxobs), "E");

	dsnn <- subset(dscore, dscore$class == "N");
	dsnnobs <- subset(dsnn, dsnn$igs != -1);
	dfwpnn <-  subset(dfwp, dfwp$ebmclass == "N");

	dff <- data.frame("N", NA, nrow(dsnnobs),  NA, NA, nrow(dfwpnn)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsc <- subset(dscore, dscore$class == "A" |  dscore$class == "B");
	dscobs <- subset(dsc, dsc$igs != -1);
	dfwpc <-  subset(dfwp, dfwp$ebmclass == "A" | dfwp$ebmclass == "B");
	dff <- data.frame("AB", nrow(dsc), nrow(dscobs),  nrow(dscobs)/nrow(dsc), nrow(dscobs)/nobstot, nrow(dfwpc)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsc <- subset(dscore, dscore$class == "A" |  dscore$class == "B" |  dscore$class == "C");
	dscobs <- subset(dsc, dsc$igs != -1);
	dfwpc <-  subset(dfwp, dfwp$ebmclass == "A" | dfwp$ebmclass == "B" |  dfwp$ebmclass == "C");
	dff <- data.frame("ABC", nrow(dsc), nrow(dscobs),  nrow(dscobs)/nrow(dsc), nrow(dscobs)/nobstot, nrow(dfwpc)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsc <- subset(dscore, dscore$class == "A" |  dscore$class == "B" |  dscore$class == "C" |  dscore$class == "D");
	dscobs <- subset(dsc, dsc$igs != -1);
	dfwpc <-  subset(dfwp, dfwp$ebmclass == "A" | dfwp$ebmclass == "B" |  dfwp$ebmclass == "C" |  dfwp$ebmclass == "D");
	dff <- data.frame("ABCD", nrow(dsc), nrow(dscobs),  nrow(dscobs)/nrow(dsc), nrow(dscobs)/nobstot, nrow(dfwpc)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsc <- subset(dscore, dscore$class == "A" |  dscore$class == "B" |  dscore$class == "C" |  dscore$class == "D"  |  dscore$class == "E");
	dscobs <- subset(dsc, dsc$igs != -1);
	dfwpc <-  subset(dfwp, dfwp$ebmclass == "A" | dfwp$ebmclass == "B" |  dfwp$ebmclass == "C" |  dfwp$ebmclass == "D"  |  dfwp$ebmclass == "E");
	dff <- data.frame("ABCDE", nrow(dsc), nrow(dscobs),  nrow(dscobs)/nrow(dsc), nrow(dscobs)/nobstot, nrow(dfwpc)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);
	
	return (dfres);
}

evaluationAnnotationClass <- function(dfres, dscore, dfwp, nobstot, classname) {

	dsc <- subset(dscore, dscore$class == classname);
	dscobs <- subset(dsc, dsc$igs != -1);
	dfwpc <-  subset(dfwp, dfwp$ebmclass == classname);

	dff <- data.frame(classname, nrow(dsc), nrow(dscobs),  nrow(dscobs)/nrow(dsc), nrow(dscobs)/nobstot, nrow(dfwpc)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	return (dfres);
}

removeMultipleInterval <- function(eval) {

	mii <- -1;
	for (i in c(1:nrow(eval))) {
		if (eval$igs[i] != -1) {
			if (eval$igs[i] == mii) {
				eval$igs[i] <- -1;
			} else {
				mii <- eval$igs[i];
			}
		}
	}

	return (eval);
}


evaluationMergedTable <- function(eval) {

	eval <- removeMultipleInterval(eval);

	dsobs <- subset(eval, eval$igs != -1);

	dfres <- NULL;
	dff <- data.frame("W", NA, nrow(dsobs),  NA, NA, NA, sum(round(eval$duration, 2)));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsx <- subset(eval, eval$class == "X");
	dsxobs <- subset(dsx, dsx$igs != -1);

	dff <- data.frame("X", NA, nrow(dsxobs),  NA, NA, NA, sum(dsx$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsnx <- subset(eval, eval$class != "X");
	dsnxobs <- subset(dsnx, dsnx$igs != -1);

	dff <- data.frame("M", NA, nrow(dsnxobs),  NA, NA, NA, sum(dsnx$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	nobstot <- nrow(dsnxobs);

	classname <- "A";
	dsclass <- subset(eval, eval$class == classname);
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame(classname, nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);
	

	classname <- "B";
	dsclass <- subset(eval, eval$class == classname);
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame(classname, nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	classname <- "C";
	dsclass <- subset(eval, eval$class == classname);
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame(classname, nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	classname <- "D";
	dsclass <- subset(eval, eval$class == classname);
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame(classname, nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	classname <- "E";
	dsclass <- subset(eval, eval$class == classname);
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame(classname, nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	classname <- "N";
	dsclass <- subset(eval, eval$class == classname);
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	
	dff <- data.frame(classname, NA, nrow(dsclassobs),  NA, NA, NA, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsclass <- subset(eval, eval$class == "A" |  eval$class == "B");
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame("AB", nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsclass <- subset(eval, eval$class == "A" |  eval$class == "B" |  eval$class == "C");
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame("ABC", nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsclass <- subset(eval, eval$class == "A" |  eval$class == "B" |  eval$class == "C" |  eval$class == "D");
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame("ABCD", nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	dsclass <- subset(eval, eval$class == "A" |  eval$class == "B" |  eval$class == "C" |  eval$class == "D"  |  eval$class == "E");
	dsclassobs <- subset(dsclass, dsclass$igs != -1);
	dsclassobsebmatch <- subset(dsclassobs, as.character(dsclassobs$ebmotion) == as.character(dsclassobs$gsebmotion));
	if (nrow(dsclassobs)) {
		ebmm <- round(nrow(dsclassobsebmatch)/nrow(dsclassobs), 3);
	} else {
		ebmm <- 0;
	}
	dff <- data.frame("ABCDE", nrow(dsclass), nrow(dsclassobs),  nrow(dsclassobs)/nrow(dsclass), nrow(dsclassobs)/nobstot, ebmm, sum(dsclass$duration));
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "ebmmatch", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	return (dfres);
}

searchMinimalIndex<- function(startingindex, vc, value) {

	for (i in c(startingindex:length(vc))) {
		if (value > vc[i]) {
			return (i-1);
		}
	}
	return (-1);
}

createSplitTable<- function(res, gs) {

	itmin = 1;
	itmax = 1;

	dfres <- NULL;
	for (i in c(1:nrow(gs))) {

		itmin <- searchMinimalIndex(itmin, res$tmin, gs$tmin[i]);
		if (itmin != -1 & res$tmin[itmin] != gs$tmin[i]) {
			dff <- data.frame(itmin,  gs$tmin[i], "tmin");
			names(dff) <- c("index", "tsplit", "class");
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);
		}

		itmax <- searchMinimalIndex(itmax, res$tmax, gs$tmax[i]);
		if (itmax != -1 & res$tmax[itmax] != gs$tmax[i]) {
			dff <- data.frame(itmax,  gs$tmax[i], "tmax");
			names(dff) <- c("index", "tsplit", "class");
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);
		}
		
		return (dfres);
	}
}

## Return the automatic intervals with interval splitted to correspond to boudary of the manual interval
splitResultIntervalForMultipleSolution	<- function(res, gs) {

	d <- computeMultipleSolution(res, gs);
	##write.table(d, "d.txt", sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE); 
	while (!is.null(d)) {
		
		dfres <- NULL;

		index <- d$index[1];
		for (i in c(1:(index-1))) {
			dff <- res[i, ];
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);
		}

		dff <- res[index, ];
		dff$tmin[1] <- d$tmin[1];	
		dff$tmax[1] <- d$tmax[1];	
		dff$duration[1] <- round(dff$tmax[1] - dff$tmin[1], 2);
		rownames(dff) <- NULL;
		dfres <- rbind(dfres, dff);
		dff <- res[index, ];
		dff$tmin[1] <- d$tmax[1];	
		dff$duration[1] <- round(dff$tmax[1] - dff$tmin[1], 2);
		rownames(dff) <- NULL;
		dfres <- rbind(dfres, dff);

		print(paste("Split interval", index, "on", nrow(res)));

		if (index == nrow(res)) {
			return (dfres);
		}

		for (i in c((index+1):nrow(res))) {
			dff <- res[i, ];
			rownames(dff) <- NULL;
			dfres <- rbind(dfres, dff);
		}

		res <- dfres;
		d <- computeMultipleSolution(res, gs);
	}

	return (res);
}


computeMultipleSolution <- function(res, gs) {

	dfres <- NULL;

	## Loop on the automatic intervals
	for (i in c(1:nrow(res))) {
		
		## Min and max of the automatic interval
		restmin <- res$tmin[i];
		restmax <- res$tmax[i];
		## 
		ir <- -1;

		## Loop on the manual intervals
		for (j in c(1:nrow(gs))) {

			## Min and max of the manual interval
			gstmin <- gs$tmin[j];
			gstmax <- gs$tmax[j];

			## If one of manual interval overlap interval i
			if (intervalOverlap(gstmin, gstmax, restmin, restmax) > 0) {
				## If it is the first time
				if (ir == -1) {

					ir <- i;
					tmin <- restmin;
				} else {
					dff <- data.frame(i, tmin, gs$tmax[j-1], res$class[i]);
					tmin <- gstmax;
					
					names(dff) <- c("index", "tmin", "tmax", "class");
					rownames(dff) <- NULL;
					dfres <- rbind(dfres, dff);

					return (dfres);
				}
			}
		}
	}

	return (NULL);
}

evaluationAnnotation <- function(res, gs) {

	ir <- numeric(nrow(res));
	ligs <- numeric(nrow(res));
	lligs <- numeric(nrow(res));

	for (i in c(1:nrow(res))) {
		restmin <- res$tmin[i];
		restmax <- res$tmax[i];
		
		ir[i] <- -1;
		ligs[i] <- -1;
		lligs[i] <- -1;

		if (nrow(gs) != 0) {
		for (j in c(1:nrow(gs))) {
				
			gstmin <- gs$tmin[j];
			gstmax <- gs$tmax[j];
			
			if (intervalOverlap(gstmin, gstmax, restmin, restmax) > 0) {
				if (ir[i] == -1) {
					ir[i] <- i;
					lligs[i] <- j;
					ligs[i] <- gs$li[j];
				} else {
					if ( (abs(gstmin - restmin) + abs(gstmax - restmax)) < (abs(gstmin - res$tmin[ir[i]]) + abs(gstmax - res$tmax[ir[i]])) ) {
					       	ir[i] <- i;
						lligs[i] <- j;
						ligs[i] <- gs$li[j];
					}	
					print(paste("multiple solution : ", j, ir[i]));
				}
				if (i != 1) {
				       	if (ir[i-1] == ir[i]) {
					if ((abs(gstmin - res$tmin[ir[i-1]]) + abs(gstmax - res$tmax[ir[i-1]])) < (abs(gstmin - res$tmin[ir[i]]) + abs(gstmax - res$tmax[ir[i]])) ) {
						ir[i-1] <-  -1;
						ligs[i-1] <- -1;
						lligs[i-1] <- -1;
					} else {
						ir[i] <- -1;
						ligs[i] <- -1;
						lligs[i] <- -1;
					}
					}
				}
			}
		}
		}
	}

	d <- data.frame(res, lligs, ligs);
	return (d);
}

evaluationAnnotationClass <- function(dfres, dscore, dfwp, nobstot, classname) {

	dsc <- subset(dscore, dscore$class == classname);
	dscobs <- subset(dsc, dsc$igs != -1);
	dfwpc <-  subset(dfwp, dfwp$ebmclass == classname);

	dff <- data.frame(classname, nrow(dsc), nrow(dscobs),  nrow(dscobs)/nrow(dsc), nrow(dscobs)/nobstot, nrow(dfwpc)*0.04);
	names(dff) <- c("class", "Npred", "Nobs", "precision", "recall", "duration");
	rownames(dff) <- NULL;
	dfres <- rbind(dfres, dff);

	return (dfres);
}


createTableOfAnnotatedIntervals <- function(d, label) {

	l = which(d$ebmotion == label);

	if (length(l) == 0) {
		cat(paste("There is not", label, "observed eyebrow movement\n")); 
		return (NULL);
	}

	index_min = numeric(0);
	index_max = numeric(0);
	ebmotion = character(0);

	index_min <- c(index_min, l[1]);
	ebmotion <- c(ebmotion, label);
	for (i in c(2:(length(l)-1))) {
		if (l[i] != l[i-1] + 1) {
			index_max <- c(index_max, l[i-1]);
			index_min <- c(index_min, l[i]);
			ebmotion <- c(ebmotion, label);
		} 
	}
	index_max <- c(index_max, l[length(l)]);

	tmin <- numeric(length(index_max));
	tmax <- numeric(length(index_max));
	time <- numeric(length(index_max));
	duration <- numeric(length(index_max));
	
	tmin <- d$time[index_min];
	tmax <- d$time[index_max];
	time <- (tmax+tmin)/2;
	duration <- tmax - tmin;

	dd <- data.frame(ebmotion, index_min, index_max, tmin, tmax, time, duration);

	return(dd);
}

createTableOfAnnotatedIntervalsWithQuality1 <- function(d, s2d, label) {

	l = which(d$ebmotion == label);

	if (length(l) == 0) {
		cat(paste("There is not", label, "observed eyebrow movement\n")); 
		return (NULL);
	}

	index_min = numeric(0);
	index_max = numeric(0);
	ebmotion = character(0);

	index_min <- c(index_min, l[1]);
	ebmotion <- c(ebmotion, label);
	for (i in c(2:(length(l)-1))) {
		if (l[i] != l[i-1] + 1) {
			index_max <- c(index_max, l[i-1]);
			index_min <- c(index_min, l[i]);
			ebmotion <- c(ebmotion, label);
		} 
	}
	index_max <- c(index_max, l[length(l)]);

	tmin <- numeric(length(index_max));
	tmax <- numeric(length(index_max));
	time <- numeric(length(index_max));
	duration <- numeric(length(index_max));
	wcoef <- numeric(length(index_max));
	qualityL <- numeric(length(index_max));
	qualityC <- numeric(length(index_max));
	qualityR <- numeric(length(index_max));
	quality <- numeric(length(index_max));
	errm <- numeric(length(index_max));
	errsd <- numeric(length(index_max));
	
	tmin <- d$time[index_min];
	tmax <- d$time[index_max];
	duration <- tmax - tmin;

	l <- levels(factor(s2d$duration));
	ns <- length(l);

	for (j in c(1:length(index_max))) {
		pdistance = 10000;
		iscale = -1;
		for (i in c(1:ns)) {
			cdistance = abs(as.numeric(l[i]) - (duration[j]/2));
			if (cdistance < pdistance) {
				pdistance <- cdistance;
				iscale <- i;
			}
		}
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);

		time[j] <- s2d$time[ns*(jposition-1)+iscale];
		wcoef[j] <- s2d$wcoef[ns*(jposition-1)+iscale];
		qualityL[j] <- s2d$qualityL[ns*(jposition-1)+iscale];
		qualityC[j] <- s2d$qualityC[ns*(jposition-1)+iscale];
		qualityR[j] <- s2d$qualityR[ns*(jposition-1)+iscale];
		quality[j] <- s2d$quality[ns*(jposition-1)+iscale];
		errm[j] <- s2d$errm[ns*(jposition-1)+iscale];
		errsd[j] <- s2d$errsd[ns*(jposition-1)+iscale];
		if (is.na(quality[j])) {
			quality[j] <- 0;
			errm[j] <- 1;
			errsd[j] <- 1;
		}
	}	

 		
	dd <- data.frame(ebmotion, index_min, index_max, tmin, tmax, time, duration, wcoef, qualityL, qualityC, qualityR, quality, errm, errsd);

	return(dd);
}

createTableOfAnnotatedIntervalsForResiduals <- function(d, s2d, itdc, label, wDparameter) {

	l = which(d$ebmotion == label);

	index_min = numeric(0);
	index_max = numeric(0);
	ebmotion = character(0);

	index_min <- c(index_min, l[1]);
	ebmotion <- c(ebmotion, label);
	for (i in c(2:(length(l)-1))) {
		if (l[i] != l[i-1] + 1) {
			index_max <- c(index_max, l[i-1]);
			index_min <- c(index_min, l[i]);
			ebmotion <- c(ebmotion, label);
		} 
	}
	index_max <- c(index_max, l[length(l)]);

	tmin <- numeric(length(index_max));
	tmax <- numeric(length(index_max));
	time <- numeric(length(index_max));
	duration <- numeric(length(index_max));
	wcoef <- numeric(length(index_max));
	qualityL <- numeric(length(index_max));
	qualityC <- numeric(length(index_max));
	qualityR <- numeric(length(index_max));
	quality <- numeric(length(index_max));
	errm <- numeric(length(index_max));
	errsd <- numeric(length(index_max));
	
	tmin <- d$time[index_min];
	tmax <- d$time[index_max];
	duration <- tmax - tmin;

	l <- levels(factor(s2d$duration));
	ns <- length(l);

	for (j in c(1:length(index_max))) {
		pdistance = 10000;
		iscale = -1;
		for (i in c(1:ns)) {
			cdistance = abs(as.numeric(l[i]) - (duration[j]/2));
			if (cdistance < pdistance) {
				pdistance <- cdistance;
				iscale <- i;
			}
		}
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);

		time[j] <- s2d$time[ns*(jposition-1)+iscale];
		wcoef[j] <- s2d$wcoef[ns*(jposition-1)+iscale];
		qualityL[j] <- s2d$qualityL[ns*(jposition-1)+iscale];
		qualityC[j] <- s2d$qualityC[ns*(jposition-1)+iscale];
		qualityR[j] <- s2d$qualityR[ns*(jposition-1)+iscale];
		quality[j] <- s2d$quality[ns*(jposition-1)+iscale];
		errm[j] <- s2d$errm[ns*(jposition-1)+iscale];
		errsd[j] <- s2d$errsd[ns*(jposition-1)+iscale];
		if (is.na(quality[j])) {
			quality[j] <- 0;
			errm[j] <- 1;
			errsd[j] <- 1;
		}
	}	

	dd <- data.frame(ebmotion, index_min, index_max, tmin, tmax, time, duration, wcoef, qualityL, qualityC, qualityR, quality, errm, errsd, wcoef);

	nt <- nrow(itdc);

	column_index_x = which(colnames(itdc) == "S");
	rx <- itdc[, column_index_x];

	wrcx <- numeric(length(index_max));

 	for (j in c(1:length(index_max))) {
		
		sigma <-duration[j]/2;
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);
		jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
		jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

		wc <- waveletCoefAtTime(itdc$time, rx, sigma, jposition, jmin, jmax, wDparameter);
		wrcx[j] <- wc[1];
	}

	dd <- data.frame(dd, wrcx);
	colnames(dd)[ncol(dd)] <- "wS";

	column_index_x = which(colnames(itdc) == "pitch");
	rx <- itdc[, column_index_x];

	wrcx <- numeric(length(index_max));

 	for (j in c(1:length(index_max))) {
		
		sigma <-duration[j]/2;
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);
		jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
		jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

		wc <- waveletCoefAtTime(itdc$time, rx, sigma, jposition, jmin, jmax, wDparameter);
		wrcx[j] <- wc[1];
	}

	dd <- data.frame(dd, wrcx);
	colnames(dd)[ncol(dd)] <- "wpitch";

	column_index_x = which(colnames(itdc) == "roll");
	rx <- itdc[, column_index_x];

	wrcx <- numeric(length(index_max));

 	for (j in c(1:length(index_max))) {
		
		sigma <-duration[j]/2;
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);
		jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
		jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

		wc <- waveletCoefAtTime(itdc$time, rx, sigma, jposition, jmin, jmax, wDparameter);
		wrcx[j] <- wc[1];
	}

	dd <- data.frame(dd, wrcx);
	colnames(dd)[ncol(dd)] <- "wroll";

	column_index_x = which(colnames(itdc) == "yaw");
	rx <- itdc[, column_index_x];

	wrcx <- numeric(length(index_max));

 	for (j in c(1:length(index_max))) {
		
		sigma <-duration[j]/2;
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);
		jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
		jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

		wc <- waveletCoefAtTime(itdc$time, rx, sigma, jposition, jmin, jmax, wDparameter);
		wrcx[j] <- wc[1];
	}

	dd <- data.frame(dd, wrcx);
	colnames(dd)[ncol(dd)] <- "wyaw";


	listpointX <- vector(mode="list", length=11);
	listpointY <- vector(mode="list", length=11);
	kc = 1;
	for (k in c(1:12)) {
		if (k != 5) {

			column_index_x = which(colnames(itdc) == paste(paste("rP", getStringNumber(k), sep=""), "x", sep=""));
			column_index_y = which(colnames(itdc) == paste(paste("rP", getStringNumber(k), sep=""), "y", sep=""));
	
			rx <- itdc[, column_index_x];
			ry <- itdc[, column_index_y];
	
			wrcx <- numeric(length(index_max));
			wrcy <- numeric(length(index_max));

 			for (j in c(1:length(index_max))) {
		
				sigma <-duration[j]/2;
				t <- time[j];
				jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);
				jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
				jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

				wc <- waveletCoefAtTime(itdc$time, rx, sigma, jposition, jmin, jmax, wDparameter);
				wrcx[j] <- wc[1];
				wc <- waveletCoefAtTime(itdc$time, ry, sigma, jposition, jmin, jmax, wDparameter);
				wrcy[j] <- wc[1];
			}

			listpointX[[kc]] <- wrcx;
			listpointY[[kc]] <- wrcy;
			dd <- data.frame(dd, wrcx, wrcy);
			colnames(dd)[ncol(dd)-1] <- paste(paste("wrP", getStringNumber(k), sep=""), "x", sep="");
			colnames(dd)[ncol(dd)] <-paste(paste("wrP", getStringNumber(k), sep=""), "y", sep="");


			column_index_x = which(colnames(itdc) == paste(paste("rcP", getStringNumber(k), sep=""), "x", sep=""));
			column_index_y = which(colnames(itdc) == paste(paste("rcP", getStringNumber(k), sep=""), "y", sep=""));
	
			rx <- itdc[, column_index_x];
			ry <- itdc[, column_index_y];
	
			wrcx <- numeric(length(index_max));
			wrcy <- numeric(length(index_max));

 			for (j in c(1:length(index_max))) {
		
				sigma <-duration[j]/2;
				t <- time[j];
				jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);
				jmin = round(max(1, (jposition+1 - ((1+wDparameter)*sigma)*25)), 0);
				jmax = round(min(nt, (jposition+1 + ((1+wDparameter)*sigma)*25)), 0);

				wc <- waveletCoefAtTime(itdc$time, rx, sigma, jposition, jmin, jmax, wDparameter);
				wrcx[j] <- wc[1];
				wc <- waveletCoefAtTime(itdc$time, ry, sigma, jposition, jmin, jmax, wDparameter);
				wrcy[j] <- wc[1];
			}

			listpointX[[kc]] <- wrcx;
			listpointY[[kc]] <- wrcy;

			print(paste("P",  getStringNumber(k), sep=""));

			dd <- data.frame(dd, wrcx, wrcy);
			colnames(dd)[ncol(dd)-1] <- paste(paste("wrcP", getStringNumber(k), sep=""), "x", sep="");
			colnames(dd)[ncol(dd)] <-paste(paste("wrcP", getStringNumber(k), sep=""), "y", sep="");

			kc <- kc + 1;
		}

	}

	return(dd);
}

createTableOfMergedAnnotatedIntervals <- function(d, s2d, label) {

	l = which(d$ebmotion == label);

	index_min = numeric(0);
	index_max = numeric(0);
	ebmotion = character(0);

	index_min <- c(index_min, l[1]);
	ebmotion <- c(ebmotion, label);
	for (i in c(2:(length(l)-1))) {
		if (l[i] != l[i-1] + 1) {
			index_max <- c(index_max, l[i-1]);
			index_min <- c(index_min, l[i]);
			ebmotion <- c(ebmotion, label);
		} 
	}
	index_max <- c(index_max, l[length(l)]);

	tmin <- numeric(length(index_max));
	tmax <- numeric(length(index_max));
	duration <- numeric(length(index_max));
	quality <- numeric(length(index_max));
	errm <- numeric(length(index_max));
	errsd <- numeric(length(index_max));
	
	tmin <- d$time[index_min];
	tmax <- d$time[index_max];
	duration <- tmax - tmin;

	l <- levels(factor(s2d$duration));
	ns <- length(l);

	for (j in c(1:length(index_max))) {
		pdistance = 10000;
		iscale = -1;
		for (i in c(1:ns)) {
			cdistance = abs(as.numeric(l[i]) - (duration[j]/2));
			if (cdistance < pdistance) {
				pdistance <- cdistance;
				iscale <- i;
			}
		}
		jposition <- index_min[j] + floor((index_max[j]-index_min[j])/2);

		quality[j] <- s2d$quality[ns*(jposition-1)+iscale];
		errm[j] <- s2d$errm[ns*(jposition-1)+iscale];
		errsd[j] <- s2d$errsd[ns*(jposition-1)+iscale];
		if (is.na(quality[j])) {
			quality[j] <- 0;
			errm[j] <- 1;
			errsd[j] <- 1;
		}
	}	

 		
	dd <- data.frame(ebmotion, index_min, index_max, tmin, tmax, duration, quality, errm, errsd);

	return(dd);
}

##
## da : The data frame containing the eyebrow motion annotation
##
normalizeEBMlabel <- function(da) {

	ebmotion <- character(nrow(da));
	if (nrow(da) != 0) {
	print(nrow(da));
	for (i in c(1:nrow(da))) {
		if (as.character(da$ebmotion[i]) == "Haussement") {
			ebmotion[i] <- "Raise";
		} else {
			ebmotion[i] <- "Frown";
		}
	}
	}
	da$ebmotion <- ebmotion;
	return (da);
}


## Add the annotations of eyebrow motion to the main data frame
addEBMPrediction <- function(d, da) {

	line_number = nrow(d);
	ebmclass = character(line_number);
	current_da_index = 1;
	number_of_annotations = nrow(da);

	for (i in c(1:line_number)) {
		if (da$tmin[current_da_index] > d$time[i]) {
			ebmclass[i] = "N";
		} else if (da$tmax[current_da_index] > d$time[i]) {
			ebmclass[i] <- as.character(da$class[current_da_index]); 
		} else {
			if (current_da_index < number_of_annotations) { 
				current_da_index = current_da_index+1;
			}
			ebmclass[i] = "N";
		}
	}

	dd <- data.frame(d, ebmclass);
	return (dd);
}
## Add the annotations of eyebrow motion to the main data frame
addEBMEvaluation <- function(d, da) {

	line_number = nrow(d);
	ebmdetection = logical(line_number);
	current_da_index = 1;
	number_of_annotations = nrow(da);

	for (i in c(1:line_number)) {
		if (da$tmin[current_da_index] > d$time[i]) {
		} else if (da$tmax[current_da_index] > d$time[i]) {
			if (da$class[current_da_index] != "N" & da$class[current_da_index] != "X") {
				if (da$igs[current_da_index] == -1) {
					ebmdetection[i] <- FALSE;
				} else {
					ebmdetection[i] <- TRUE;
				}
			}
		} else {
			if (current_da_index < number_of_annotations) { 
				current_da_index = current_da_index+1;
			}
		}
	}

	dd <- data.frame(d, ebmdetection);
	return (dd);
}

loadPhonemeTable <- function(directory, filename) {

	## Read the annotation table
	d <- read.table(paste(directory, filename, sep="/"), h=TRUE, sep="\t", comment.char="");

	if (ncol(d) == 4) {
		dd <- data.frame(d$tmin, d$tmax, d[,3]);
	} else {
		dd <- data.frame(d$tmin, d$tmax, d[,4]);
	}
	colnames(dd)[1] = "tmin";
	colnames(dd)[2] = "tmax";
	colnames(dd)[3] = "phoneme";
	
	return (dd);
}

addPhoneme <- function(d, da) {

	line_number = nrow(d);
	phoneme = character(line_number);
	current_da_index = 1;
	number_of_annotations = nrow(da);

	for (i in c(1:line_number)) {
		if (da$tmax[current_da_index] > d$time[i]) {
			phoneme[i] <- as.character(da$phoneme[current_da_index]);
		} else {
			if (current_da_index < number_of_annotations) { 
				current_da_index = current_da_index+1;
			}
			phoneme[i] <- as.character(da$phoneme[current_da_index]);
		}
	}

	dd <- data.frame(d, phoneme);
	return (dd);
}

## Add the annotations of eyebrow motion to the main data frame
addAnnotationToIntrafaceOutput <- function(d, da) {

	line_number = nrow(d);
	ebmotion = character(line_number);
	current_da_index = 1;
	number_of_annotations = nrow(da);

	for (i in c(1:line_number)) {
		if (current_da_index > number_of_annotations || da$tmin[current_da_index] > d$time[i]) {
			ebmotion[i] = "Nothing";
		} else if (da$tmax[current_da_index] > d$time[i]) {
			if (as.character(da$ebmotion[current_da_index]) == "Haussement") {
				ebmotion[i] = "Raise";
			} else {
				ebmotion[i] = "Frown";
			}
		} else {
			if (current_da_index < number_of_annotations) { 
				current_da_index = current_da_index+1;
			}
			ebmotion[i] = "Nothing";
		}
	}

	dd <- data.frame(d, ebmotion);
	return (dd);
}
