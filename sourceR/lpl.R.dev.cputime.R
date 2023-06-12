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
