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

## lpl.R.dev.hmm.R, creator S. Rauzy, LPL, 03/04/2019

##
## Return the state probability of transition 
## austmdf : The data frame containing the state probabilities of transitions
## svdf : The data frame containing the smile variable and its associated bins
##
## Return the Viterbi solution
##
lpl.R.dev.hmm.viterbiSolution <- function(austmdf, svdf) {

	ns <- 5;
	n <- nrow(svdf);

	NOT_ACTIVE <- 1;
	NOT_YET_VISITED <- -2;

	cumul_logproba <- matrix(0, nrow=5, ncol=n);	
	max_logproba <- matrix(NOT_ACTIVE, nrow=5, ncol=n);
	max_index <- matrix(NOT_YET_VISITED, nrow=5, ncol=n);

	## We start with a S0 state (index 1)
	j <- 1;
	## Loop on target
	for (i in c(1:ns)) {
		pt <- lpl.R.dev.smad.getStateProbabilityOfTransition(austmdf, svdf$bi[j], 1, i);
		## If the target state is reachable
		if (pt != 0) {
			cumul_logproba[i, j] <- log(pt); 
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
					pt <- lpl.R.dev.smad.getStateProbabilityOfTransition(austmdf, svdf$bi[j], k, i);
					if (is.na(pt) || pt < 0) {
						cat("warning");
					}
					## If the target state is reachable
					if (pt != 0) {
						if (max_index[i, j] == NOT_YET_VISITED) {
							cumul_logproba[i, j] <- log(pt) + cumul_logproba[k, j-1];
						} else {
							cumul_logproba[i, j] <- log(exp(cumul_logproba[i, j]) + exp(log(pt) + cumul_logproba[k, j-1]));
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
		if (nsa == 5) {
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
		cat(paste(max_logproba[i, j], "\n", sep=""));
	}
	cat(paste("max = ", vmax, " for state ", imax, "\n", sep=""));

	sol <- numeric(n);
	lpmax <- numeric(n);
	sol[n] <- imax;
	lpmax[n] <- max_logproba[sol[n], n];
	for (j in c((n-1):1)) {
		sol[j] <- max_index[sol[j+1], j+1];
		cat(paste("max = ",j, " state ", sol[j], "\n", sep=""));
		lpmax[j] <- max_logproba[sol[j], j];	
	}
	return (max_logproba);
}
