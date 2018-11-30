## lpl.R.dev.cputime.R, creator S. Rauzy, LPL, 04/12/2016

getTimeInHundredthSecond <- function() {
	
	t <- as.numeric(Sys.time())*100;
	return (t);
}

durationInSecond <- function(startingtimehs, endingtimehs) {

	return ((endingtimehs-startingtimehs)/100);
}

estimateDurationInSecond <- function(operationduration, nbroperations) {

	return (operationduration*nbroperations);
}
