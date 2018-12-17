##
## lpl.R.dev.cputime.R, creator S. Rauzy, LPL, 04/12/2016
##

## 
## Return the time in hundredth of second
##
getTimeInHundredthSecond <- function() {
	
	t <- as.numeric(Sys.time())*100;
	return (t);
}

## 
## Return the duration in second
##
## startingtimehs : The starting time in hundredth of second
## endingtimehs : The ending time in hundredth of second
##
durationInSecond <- function(startingtimehs, endingtimehs) {

	return ((endingtimehs-startingtimehs)/100);
}

## 
## Return the estimate duration time in second
##
## operationduration : The duration in second of a single operation
## nbroperations : The number of operations
##
estimateDurationInSecond <- function(operationduration, nbroperations) {

	return (operationduration*nbroperations);
}
