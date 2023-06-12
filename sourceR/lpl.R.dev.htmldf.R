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

## lpl.R.dev.htmdf.R, creator S. Rauzy, LPL, 15/02/2017

##
## Open the html tag
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.openHtmlTag <- function(df) {

	line <- "<HTML>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Open the html body
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.openHtmlBody <- function(df) {

	line <- "\t<BODY BGCOLOR=\"#FFFFFF\">";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Open the html table
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.openHtmlTable <- function(df) {

	line <- "\t<TABLE HALIGN=\"LEFT\" BGCOLOR=\"#FFFFFF\" CELLPADDING=\"2\" CELLSPACING=\"1\" BORDER=\"1\">";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Open the line of an html table
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.openHtmlTableLine <- function(df) {

	line <- "\t\t<TR VALIGN=\"TOP\">";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Add a cell to the line of an html table
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.addCellToHtmlTableLine <- function(df, celltext) {

	line <- paste("\t\t\t<TD HALIGN=\"LEFT\">", celltext, "</TD>");
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Close the line of an html table
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.closeHtmlTableLine <- function(df) {

	line <- "\t\t</TR>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Close the html table
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.closeHtmlTable <- function(df) {

	line <- "\t</TABLE>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Return the html bold version of the line
##
## line : The line as a String
## return the html bold version of the line
##
lpl.R.dev.htmldf.boldHtml <- function(line) {

	return (paste("<B>", line, "</B>", sep=""));
}


##
## Close the html tag
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.closeHtmlTag <- function(df) {

	line <- "</HTML>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Close the html tag
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.closeHtmlBody <- function(df) {

	line <- "\t</BODY>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}
##
## Create the html header
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.openHtmlHeader <- function(df, title) {

	line <- "\t<HEAD>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- paste("\t\t<TITLE>", title, "</TITLE>", sep="");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);

	line <- "\t\t<META http-equiv=\"Content-Type\" Content=\"text/html; charset=ISO-8859-15\"/>";
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Close the html header
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.closeHtmlHeader <- function(df) {

	line <- "\t</HEAD>";
	dfl <- data.frame(line);

	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}

##
## Close the html tag
##
## df : The html document as a data frame
##
lpl.R.dev.htmldf.addLine <- function(df, line) {

	line <- paste(line, "<BR>");
	dfl <- data.frame(line);
	rownames(dfl) <- NULL;
	df <- rbind(df, dfl);
}


