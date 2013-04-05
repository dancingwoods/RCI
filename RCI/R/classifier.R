# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Classifiers for segmentation of calcium images                                             
# 
# License information
# This file is part of the R package RCI.
# 
# RCI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# RCI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with RCI.  If not, see <http://www.gnu.org/licenses/>.
# 
###############################################################################


InitiateMaskClassifier <- function(path){
# classifier - the classifier object
# db - a SQLite database to store 	 
}


#-
#' Pulls data from a mask database into a data.frame that has one row for
#' each mask with the label, id, and feature values of that mask in the columns
#' 
#' @param db the mask database object
#' @param labeled boolean, should the results be restricted to masks that have been labeled
#' 
#' @return a data.frame with the extracted data
#' 
#' @export
#-
PullData <- function(db, labeled=T){
	features <- dbGetQuery(db$db, "SELECT id, tag FROM features")
	df <- dim(features)
	
	if(nrow(features>0)){
		featurestring <- paste("MAX(CASE WHEN fm.featureid=",
			features[1,1], " THEN fm.fvalue END) as ", features[1,2])
		for(i in 2:nrow(features)){
			featurestring <- paste(featurestring, ", MAX(CASE WHEN fm.featureid=",
				features[i,1], " THEN fm.fvalue END) as '", features[i,2], "' ")
		}
	
		qry <- paste("SELECT m.label, m.id, ",
			featurestring,
			"FROM masks AS m ",
			"JOIN feats_masks AS fm ON m.id = fm.maskid")
		
		if(labeled){
			qry <- paste(qry, "WHERE m.label > 0")
		}
			
		qry <- paste(qry, "GROUP BY m.id")
		return(as.data.frame(dbGetQuery(db$db, qry)))
	}else{
		return(NA)
	}
}


TestData <- function(db){
	
}



Train <- function(maskcf){
	
}