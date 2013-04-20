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


InitiateMaskClassifier <- function(trainingdata){
# classifier - the classifier object
# training - the training data
	td <- trainingdata
	td <- td[,-which(names(td)=="id")]
	td[is.na(td)]=0
	cf <- randomForest(x = td[,-c(which(names(td)=="label"), which(names(td)=="experiment"))], y = as.factor(td[,which(names(td)=="label")]))
	ret <- list("training"=trainingdata, "cf"=cf)
	return(ret)
}

#-
#' @details Predicts an experiment using the classifier and MWIS
#' 
#' @param classifier the classifier
#' @param dn the database to predict
#' @param enforcelabels boolen, should the final segmentation be forced to correspond
#' to the hand labels
#' 
#' @return null, modifies the database
#' 
#' @useDynLib RCI compmasksC
#-
PredictExperiment <- function(classifier, db, enforcelabels=T){
	data <- PullData(db, F)
	data[is.na(data)]=0
	labs <- data[,which(names(data)=="label")]
	data  <- data[, -(which(names(data)=="label"))]
	probs <- predict(classifier$cf, data, type="prob")
	ovlp <- dbGetQuery(db$db, "select maskid1, maskid2 from edges")
	
	# append mask id and labels to probs matrix to keep track when sorting
	probs <- cbind(data[,which(names(data)=="id")], labs, probs)
	probs[is.na(probs)]=0

	# sort by highest non-1 probability
	maxv <- apply(probs[,-(1:3)], 1, max)
	probs <- probs[order(maxv, decreasing=T),]
	
	# change ovlp numbers to refer to the row numbers in probs instead
	#  of ids
	for(i in 1:nrow(probs)){
		ovlp[which(ovlp==probs[i,1])]=i
	}
	
	out <- .C("labelexperimentC",
		as.double(as.vector(probs[,-1])),
		as.integer(nrow(probs)),
		as.integer(ncol(probs)-1),
		as.integer(as.vector(ovlp)),
		as.integer(nrow(ovlp)),
		pvec = as.integer(rep(0, nrow(probs))),
		as.integer(enforcelabels)
	)
	
	# 
	# # results vector to store predicted class
	# pvec  <- rep(NA, nrow(probs))
	# 
	# # iterate through and check the labels (if enforcing labels)
	# if(enforcelabels){
	# 	for(i in 1:nrow(probs)){
	# 		# if not already locked
	# 		if(is.na(pvec[i])){
	# 			# if negative, transfer label TODO rethink this since lots of automatic negatives
	# 			# if(probs[i,2]==1){
	# 			# 		pvec[i] <- 1
	# 			# if has a label, set any overlapping masks to negative
	# 			#}else 
	# 			if(probs[i,2]>1){
	# 				pvec[i] <- probs[i,2]
	# 				oids <- c(ovlp[which(ovlp[,1]==probs[i,1]),2], ovlp[which(ovlp[,2]==probs[i,1]), 1])
	# 				oinds <- which(probs[,1] %in% oids)
	# 				oinds <- oinds[which(oinds>i)]
	#  				pvec[oinds] <- 1
	# 			}
	# 		}
	# 	}
	# }
	# 
	# # Now go through remaining masks and choose any that are segmented (with unknown label)	
	# for(i in 1:nrow(probs)){
	# 	if(is.na(pvec[i])){
	# 		# If haven't already constrained this mask to be negative
	# 		wmax <- which.max(probs[i,-(1:2)])
	# 		pvec[i] <- wmax
	# 		if(wmax!=1){
	# 			# Assign overlapping masks to negative value
	# 			oids <- c(ovlp[which(ovlp[,1]==probs[i,1]),2], ovlp[which(ovlp[,2]==probs[i,1]), 1])
	# 			oinds <- which(probs[,1] %in% oids)
	# 			oinds <- oinds[which(oinds>i)]
	#  			pvec[oinds] <- 1
	# 		}
	# 	}
	# }
	
	addseg <- cbind(out$pvec, probs[,1])
	addseg <- as.data.frame(addseg)
	names(addseg) <- c("segmentation", "id")
	
	dbBeginTransaction(db$db)
	dbGetPreparedQuery(db$db, "UPDATE masks SET segmentation=? where id=?", bind.data=addseg)
	dbCommit(db$db)
	

}


# Checks each database in the given directory and pulls any
# labeled data
PullAllData <- function(directory){
	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
	
	count <- 1
	while(count <= length(dblist)){
		cat("-- ", dblist[count], "\n")
		tmp <- PullData(ConMaskDb(paste(directory, dblist[count], sep="")))
		if(nrow(tmp)>0){
			break
		}
	}
	if(nrow(tmp)==0){
		cat("No labeled data found\n")
		return(NA)
	}else{
		namecol <- as.data.frame(rep(dblist[count], nrow(tmp)))
		names(namecol) <- c("experiment")
		ret <- cbind(namecol, tmp)
	}
	count <-  count+1
	print(names(ret))
	for(i in count:length(dblist)){
		cat("----", dblist[i], "\n")
		tmp <- PullData(ConMaskDb(paste(directory, dblist[i], sep="")))
		if(nrow(tmp)>0){
			if( ncol(tmp) == (ncol(ret)-1) ){
				namecol <- as.data.frame(rep(dblist[i], nrow(tmp)))
				names(namecol) <- c("experiment")
				tmp <- cbind(namecol, tmp)
				print(names(tmp))
		 		ret <- rbind(ret, tmp)
			}
		}
		
	}
	return(ret)
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
	scaleids <- grep(paste("Scale_.*", sep=""), features[,2])
	ids <- 1:nrow(features)
	ids <- ids[-scaleids]
	if(nrow(features>0)){
		featurestring <- paste("MAX(CASE WHEN fm.featureid=",
			features[ids[1],1], " THEN fm.fvalue END) as '", features[ids[1],2], "'")
		for(i in 2:length(ids)){
			featurestring <- paste(featurestring, ", MAX(CASE WHEN fm.featureid=",
				features[ids[i],1], " THEN fm.fvalue END) as '", features[ids[i],2], "' ")
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



