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

## Load the libraries in sourceR
	
## Load the library dealing with cpu time
source("sourceR/lpl.R.dev.cputime.R");
## Load the library to transform Intraface output	
source("sourceR/lpl.R.dev.intraFaceOutputAnalysis.R");
## Load the libraries to transform OpenFace output	
source("sourceR/lpl.R.dev.openFaceOutputAnalysis.R");
source("sourceR/lpl.R.dev.faceOutputAnalysis.R");
## Load the library defining the hmad operations
source("sourceR/lpl.R.dev.hmad.R");
## Load the libraries specific to eyebrow motion detection
source("sourceR/lpl.R.dev.ebmad.R");
source("sourceR/lpl.R.dev.ebmotion.R");
## Load the libraries specific to smile detection
source("sourceR/lpl.R.dev.smad.R");
## Load the libraries specific to blink detection
source("sourceR/lpl.R.dev.bmad.R");
## Load the libraries required during the analysis
source("sourceR/lpl.R.dev.wavelets.R");
source("sourceR/lpl.R.dev.elan.R");
source("sourceR/lpl.R.dev.htmldf.R");
source("sourceR/lpl.R.dev.utils.R");
source("sourceR/lpl.R.dev.hmm.R");

