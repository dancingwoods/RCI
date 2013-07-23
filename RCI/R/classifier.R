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

# The commented functions are for evaluating segmentations.  They are specific to the
# evaluations performed for my dissertation, but are here for reference for 
# future work

# MergeEvaluateDataLimited <- function(directory){
# 	filelist <- list.files(directory)
# 	ret  <- matrix(NA, 2*length(filelist), 9)
# 	for(i in 1:length(filelist)){
# 		load(paste(directory, filelist[i], sep=""))
# 		ntrain <- length(strsplit(tail(strsplit(filelist[i], "_")[[1]],1), "db")[[1]])-1
# 		ret[(i*2-1):(i*2),1:8] <- as.matrix(get(filelist[i]))
# 		ret[(i*2-1):(i*2),9] = c(ntrain, ntrain)
# 	}
# 	return(ret)
# }
# 
# SumPerformance <- function(memat, ntrain){
# 	memat <- memat[which(memat[,9]==ntrain),]
# 	data2 <- memat[which(memat[,1]==2),]
# 	data3 <- memat[which(memat[,1]==3),]
# 	return(rbind(apply(data2, 2, sum), apply(data3, 2, sum)))
# }

# # merge the evaluation results in the given directory
# MergeEvaluate <- function(directory){
# 	filelist <- list.files(directory)
# 	load(paste(directory, filelist[1], sep=""))
# 	name <- paste("perf_", strsplit(filelist[1], split=".R")[[1]], sep="")
# 	perf <- get(name)
# 	for(i in 2:length(filelist)){
# 		load(paste(directory, filelist[i], sep=""))
# 		name <- paste("perf_",strsplit(filelist[i], split=".R")[[1]], sep="")
# 		perf <- perf + get(name)
# 	}
# 	return(perf)
# }
# 
# # perfrow is which row of performance (class) is under consideration
# MergeEvaluateROC <- function(directory, perfrow=1){
# 	filelist <- list.files(directory)
# 	#filelist <- filelist[grep("perf", filelist)]
# 	ret <- matrix(NA, length(filelist), 9)
# 	for(i in 1:length(filelist)){
# 		load(paste(directory, filelist[i], sep=""))
# 		name <- paste("perf_",strsplit(filelist[i], split=".Rdata")[[1]], sep="")
# 		perf <- as.matrix(get(name))
# 		ret[i,1:8] <- perf[perfrow,] 
# 		ret[i,9] <- as.double(tail(strsplit(name, "_")[[1]], 1))
# 	}
# 	return(ret)
# }


# EvaluateDataLimited <- function(directory, outdirectory, limit=100){
# 	data <- PullAllData(directory)
# 	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
# 	nexp <- length(dblist)
# 	for(n in 5:(nexp-1)){
# 		cat(n, " experiments\n")
# 		# For each number of training experiments, compute performance
# 		# with all possible sets of training data of that size
# 
# 		# all combinations of training data of this size, each column is
# 		# the names of the selected experiments
# 		combos <- combn(dblist, n)
# 		trainsets <- sample(1:ncol(combos))[1:min(limit, ncol(combos))]
# 		for(i in trainsets){
# 			cat(i, " of ", ncol(combos), " training combos\n")
# 			class <- InitiateMaskClassifier(data[which(data[,1] %in% combos[,i]),])
# 			strclass <- paste(strsplit(combos[,i], split=".sqlite"), collapse="")
# 			
# 			for(d in dblist){
# 				cat(d, " ")
# 				# for each experiment not used to create classifier, test performance
# 				if(!(d %in% combos[,i])){
# 					db <- ConMaskDb(paste(directory, d, sep=""))
# 					PredictExperiment(class, db)
# 					dbstr <- strsplit(d, split=".sqlite")[[1]]
# 					assign(paste("perf_", dbstr, "_", strclass,sep=""), EvaluateSegmentation(db))
# 					#print(get(paste("perf_", dbstr, "_", strclass,sep="")))
# 					save(list=paste("perf_", dbstr, "_", strclass,sep=""), file=paste(outdirectory,"perf_", dbstr, "_", strclass,sep=""))					
# 				}
# 			}
# 	 	}	
# 	 }	
# }

# # Evaluates the performance of each of the experiments in the directory after
# # training on the remaining experiments.  Save the results to the outdirectory.
# # Threshvals specify the minimum necessary threshold for a mask to be selected.
# EvaluateCV <- function(directory, outdirectory, threshvals=NULL){
# 	data <- PullAllData(directory)
# 	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
# 	for(dbname in dblist){
# 		print(dbname)
# 		dbstr <- strsplit(dbname, split=".sqlite")[[1]]
# 		class <- InitiateMaskClassifier(data[-which(data[,1]==dbname),])
# 		db <- ConMaskDb(paste(directory, dbname, sep=""))
# 		if(is.null(threshvals)){
# 			PredictExperiment(class, db)
# 			assign(paste("perf_", dbstr, sep=""), EvaluateSegmentation(db))
# 			print(get(paste("perf_", dbstr,sep="")))
# 			save(list=paste("perf_", dbstr,sep=""), file=paste(outdirectory, dbstr,sep=""))	
# 		}else{
# 			for(t in threshvals){
# 				cat(t, " ")
# 				PredictExperiment(class, db, thresh=t)
# 				assign(paste("perf_", dbstr, "_", t,sep=""), EvaluateSegmentation(db))
# 				print(get(paste("perf_", dbstr, "_", t,sep="")))
# 				save(list=paste("perf_", dbstr, "_", t,sep=""), file=paste(outdirectory, dbstr, "_", t,sep=""))
# 			}
# 		}
# 	}
# }

# EvaluateCVConfidence <- function(directory, outdirectory){
# 	data <- PullAllData(directory)
# 	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
# 	for(dbname in dblist[-1]){
# 		print(dbname)
# 		dbstr <- strsplit(dbname, split=".sqlite")[[1]]
# 		class <- InitiateMaskClassifier(data[-which(data[,1]==dbname),])
# 		db <- ConMaskDb(paste(directory, dbname, sep=""))
# 		PredictExperiment(class, db)
# 		assign(paste("conf_", dbstr, sep=""), EvaluateConfidence(db, class))
# 		print(get(paste("conf_", dbstr,sep="")))
# 		save(list=paste("conf_", dbstr,sep=""), file=paste(outdirectory, dbstr,sep=""))	
# 	}	
# }

#-
#' Creates a mask classifier based on the given training data
#' 
#' @param trainingdata the training data, as returned by PullAllData or PullData
#' 
#' @return a mask classifier
#' 
#' @export
#-
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
#' Predicts an experiment using the classifier and MWIS
#' TODO - this is currently a hack with a heuristic to find cliques.  should really find connected
#' components and solve the MWIS
#' 
#' @param classifier the classifier
#' @param dn the database to predict
#' @param enforcelabels boolen, should the final segmentation be forced to correspond
#' to the hand labels
#' 
#' @return null, modifies the database
#' 
#' @export
#-
PredictExperiment <- function(classifier, db, thresh=NULL){
	data <- PullData(db, F)
	data[is.na(data)]=0
	labs <- data[,which(names(data)=="label")]
	data  <- data[, -(which(names(data)=="label"))]
	probs <- predict(classifier$cf, data, type="prob")
	
	# append mask id and labels to probs matrix to keep track when sorting
	probs <- cbind(data[,which(names(data)=="id")], labs, probs)
	probs[is.na(probs)]=0

	# maxval is assigned class
	maxval <- apply(probs[,-(1:2)], 1, which.max)
	maxp <- apply(probs[,-(1:2)], 1, max)
	if(!is.null(thresh)){
		maxval[which(maxp < thresh)] =1
	}
	
	# find groups for each class
	groups <- rep(0, length(maxval))
	for(val in unique(maxval)){
		if(val!=1){
			# this is an ROI class
			w <- which(maxval==val)
			groups[w] <- GetCliques(db, as.vector(probs[w,1]), max(groups, na.rm=T)+1)
		}
	}
	
	# results vector to store predicted class
	pvec  <- rep(NA, nrow(probs))
	pvec[which(maxval==1)]=1
	if(any(is.na(groups))){
		pvec[which(is.na(groups))]=1
		uniquegroups <- unique(groups[-which(is.na(groups))])
	}else{
		uniquegroups <- unique(groups)
	}
	
	if(0 %in% uniquegroups){
		uniquegroups <- uniquegroups[-which(uniquegroups==0)]
	}
	# for each group, assign the best mask the appropriate class
	for(grp in uniquegroups){
		w <- which(groups==grp)
		class <- maxval[w[1]]
		wbest <- which.max(probs[w,2+class])
		pvec[w] <- 1
		pvec[w[wbest]] <- class
	}
	
	# Now need to check overlap and do something about it.
	# This is pretty ugly right now.  could almost certainly be much better
	
	ovlp <- as.matrix(dbGetQuery(db$db, "select maskid1, maskid2 from edges"))	
	curids <- as.vector(probs[which(pvec>1),1])
	curovlp <- ovlp[which((ovlp[,1] %in% curids) & (ovlp[,2] %in% curids)),,drop=F]
	
	# print(curids[order(curids)])
	# print(curovlp)
	
	# While there are any currently overlapping selected masks
	count <- 1
	while(nrow(curovlp)>0){
		ind1 <- which(probs[,1]==curovlp[1,1])
		ind2 <- which(probs[,1]==curovlp[1,2])
	
		grp1 <- groups[ind1]
		grp2 <- groups[ind2]
		
		class1 <- pvec[ind1]
		class2 <- pvec[ind2]
		
		grpinds1 <- which(groups==grp1)
		grpinds2 <- which(groups==grp2)
		
		bestval <- 0
		bestind1 <- NA
		bestind2 <- NA
		
		# For each mask in group 1
		for(i in grpinds1){
			qry <- paste("select maskid1, maskid2 from edges where maskid1=",
							probs[i,1], " OR maskid2 = ", probs[i,1])
			restrict <- as.matrix(dbGetQuery(db$db, qry))
	
			# These are the indices of the masks which overlap with the mask
			# under consideration
			connectioninds <- which(probs[,1] %in% unique(restrict))
	
			# tmp is the vector of indices in group 2 which overlap with the mask
			# under consideration
			tmp <- grpinds2 %in% connectioninds
	
			# if there are any overlapping masks
			if(any(tmp)){
				# grpinds2tmp are the ones that don't overlap
				grpinds2tmp <- grpinds2[-which(tmp)]
			}else{
				# all masks in group2 
				grpinds2tmp <- grpinds2
			}
			
			# For each of the masks in group 2 which don't overlap with the mask under consideration
			for(j in grpinds2tmp){
				# if the sum of the 2 non-overlapping masks is better than the current best pair
				#  keep track of these two indices
				if(maxp[i] + maxp[j] > bestval){
					bestind1 <- i
					bestind2 <- j
				}
			
			}
		}
		
		
		# Remove the originally overlapping masks
		pvec[ind1] <- 1
		pvec[ind2] <- 1
		
		# Assign the best mask the class labels
		pvec[bestind1] <- class1
		pvec[bestind2] <- class2
		
		# Get the mask IDs of the currently selected masks
		curids <- as.vector(probs[which(pvec>1),1])
		
		# Update curovlp 
		curovlp <- ovlp[which((ovlp[,1] %in% curids) & (ovlp[,2] %in% curids)),,drop=F]
		
		#if in infinite loop - TODO fix so that this doesn't happen
		count <- count+1
		if(count>100){
			# remove one of each pair of overlapping masks (the one with the lower value)
			for(ri in 1:nrow(curovlp)){
				ind1 <- which(probs[,1]==curovlp[1,1])
				ind2 <- which(probs[,1]==curovlp[1,2])
				class1 <- pvec[ind1]
				class2 <- pvec[ind2]
				
				if(maxp[ind1]>maxp[ind2]){
					pvec[ind2] <- 1 
				}else{
					pvec[ind1] <- 1
				}
			}
			break
		}
		
	}

	addseg <- cbind(pvec, probs[,1])
	addseg <- as.data.frame(addseg)
	names(addseg) <- c("segmentation", "id")
	
	dbBeginTransaction(db$db)
	dbGetPreparedQuery(db$db, "UPDATE masks SET segmentation=? where id=?", bind.data=addseg)
	dbCommit(db$db)
}

#-
#' Takes the current labels and segmentation in the given database and evaluates perfomance
#' of the segmenter against the hand labels.
#' For each group of labels cells, a false negative if no mask overlaps, false positive if segmented
#' mask overlaps with no group, marginal result if segmented cell overlaps with but is not of of the 
#' labeled masks.
#' 
#' @param db the database for which to evaluate the segmentation.  should have both
#' a segmentation and hand labels.
#' 
#' @return a data.frame with a row for each class and columns specifying performance
#' @export
#-
EvaluateSegmentation <- function(db){
	labels <- dbGetQuery(db$db, "select distinct label from masks")[,1]
	labels <- labels[which(labels>1)]	
	labelstring=paste(labels, collapse=", ")
	
	# 1 - class id 
	# 2 - labeled groups
	# 3 - correct
	# 4 - correct mislabeled
	# 5 - marginal
	# 6 - marginal mislabeled
	# 7 - missed
	# 8 - new 
	ret <- matrix(0, length(labels), 8)
	ret[,1] <- labels

	# data is 3 column: id, label, group
	data <- dbGetQuery(db$db, paste("select id, label from masks where label in (", labelstring, ")"))
	data <- cbind(data, rep(0, nrow(data)))

	# Assign clique group numbers
	for(l in labels){
		mval <- max(data[,3], na.rm=T)
		w <- which(data[,2]==l)
		data[w,3] <- GetCliques(db, data[w,1], mval+1)
	}

	# seg is a matrix giving the id and segmentation for all segmented data
	seg <- dbGetQuery(db$db, paste("select id, segmentation from masks where segmentation in (", labelstring, ")"))

	# assign number of annotated ROIs
	for(l in labels){
		ret[which(ret[,1]==l),2] <- length(unique(data[which(data[,2]==l),3]))
	}
	
	# for each annotated group, assign to appropriate class and remove mask from seg (if it exists)
	# check for exact matches first so that won't get removed as a marginal match for a different ROI
	rois <- unique(data[,3], na.rm=T)
	for(roi in rois){
		# the masks that were marked as ok for this ROI
		maskids <- data[which(data[,3]==roi),1]
		
		# the label for the ROI
		lab <- data[which(data[,1]==maskids[1]), 2]
		labrow <- which(ret[,1]==lab)

		w = which(seg[,1] %in% maskids)
		if(length(w)==1){
			# there was one matching mask.  this is either correct or (mis)correct
			if(seg[w,2]==lab){
				#correct
				ret[labrow,3] = ret[labrow,3]+1
				seg <- seg[-w,]
			}else{
				#mis(correct)
				ret[labrow,4] = ret[labrow,4]+1
				seg <- seg[-w,]
			}
			rois <- rois[-which(rois==roi)]
		}
	}
	
	while(length(rois>=1)){
		for(roi in rois){
			# there was no matching mask - check for overlapping mask
			maskids <- data[which(data[,3]==roi),1]
			
			# the label for the ROI
			lab <- data[which(data[,1]==maskids[1]), 2]
			labrow <- which(ret[,1]==lab)
			# this is either a missed ROI or a marginal ROI
			# check if it is a marginal ROI by finding all masks in seg which overlap with
			#  any of the masks for the ROI
			ovlp <- as.matrix(dbGetQuery(db$db, paste("select maskid1, maskid2 from edges where (maskid1 in (",
				paste(maskids, collapse=", "), ") and maskid2 in (", paste(seg[,1], collapse=", "), 
				")) or (maskid1 in (", paste(seg[,1], collapse=", "), ") and maskid2 in (", 
				paste(maskids, collapse=", "),"))")))
					
			#print(ovlp)
			# if there is any overlap, this ROI has at least one marginal match
			if(nrow(ovlp)>=1){
				# which of the segmented masks overlap with this ROI
				segovlprows <- which(seg[,1] %in% as.vector(ovlp))
				if(length(segovlprows)==1){
					# there was one marginal mask
					if(seg[segovlprows,2]==lab){
						ret[labrow,5] <- ret[labrow,5] + 1
						seg <- seg[-segovlprows,]
					}else{
						ret[labrow,6] <- ret[labrow,6] + 1
						seg <- seg[-segovlprows,]						
					}
					rois <- rois[-which(rois==roi)]
				}else{
					
					# there were multiple masks overlapping with this ROI.
					if(any(seg[segovlprows,2]==lab)){
						# at least one matched the ROI label.  count the match and remove the ROI
						ret[labrow,5] <- ret[labrow,5] + 1
						rois <- rois[-which(rois==roi)]
						seg <- seg[-segovlprows[which(seg[segovlprows,2]==lab)[1]],]
					}else{
						# no correct label overlapping.  count a (mis)marginal and remove the ROI
						# and the first overlap
						# NOTE: this is fragile.  might be an error with very close ROIs of same type
						ret[labrow,6] <- ret[labrow,6] + 1
						rois <- rois[-which(rois==roi)]
						seg <- seg[-segovlprows[1],]
					}
				}
			}else{
				
				# this is a missed ROI
				ret[labrow, 7] = ret[labrow, 7] + 1
				rois <- rois[-which(rois==roi)]
			}
		}
	}

	# Any remaining masks in seg are false positives
	for(l in labels){
		labrow <- which(ret[,1]==l)
		ret[labrow, 8] = ret[labrow, 8] + length(which(seg[,2]==l))
	}
	
	ret <- as.data.frame(ret)
	names(ret) <- c("class", "labeled", "correct", "(mis)correct", "marginal", "(mis)marginal", "missed", "new")
	return(ret[order(ret[,1]),])
}

#-
#' Computes several confidence measures on the segmentation in a database these include
#' - number of masks returned by the classifier
#' - min, max, and mean probability assigned by the classifier
#' 
#' @param db the database for which to evaluate the segmentation confidence. 
#' @param class the classifier used to generate the segmentation
#' 
#' @return a data.frame with a row for each ROI and columns specifying confidence measures
#' @export
#-
EvaluateConfidence <- function(db, class){
	labels <- dbGetQuery(db$db, "select distinct label from masks")[,1]
	labels <- labels[which(labels>1)]	
	labelstring=paste(labels, collapse=", ")
	
	
	# for each ROI and false positive, return a matrix
	# 1 - the class (1,2,3)
	# 2 - status (1=correct, 2=marginal, 3=false pos)
	# 3 - number of classified masks
	# 4 - mean probability
	# 5 - min probability
	# 6 - max probability

	ret <- matrix(0, 1, 6)

	
	####### Getting info on classifier results
	data <- PullData(db, F)
	data[is.na(data)]=0
	labs <- data[,which(names(data)=="label")]
	data  <- data[, -(which(names(data)=="label"))]
	probs <- predict(class$cf, data, type="prob")
	
	# append mask id and labels to probs matrix to keep track when sorting
	probs <- cbind(data[,which(names(data)=="id")], labs, probs)
	probs[is.na(probs)]=0

	# maxval is assigned class
	maxval <- apply(probs[,-(1:2)], 1, which.max)
	maxp <- apply(probs[,-(1:2)], 1, max)
	
	# find groups for each class
	groups <- rep(0, length(maxval))
	for(val in unique(maxval)){
		if(val!=1){
			# this is an ROI class
			w <- which(maxval==val)
			groups[w] <- GetCliques(db, as.vector(probs[w,1]), max(groups, na.rm=T)+1)
		}
	}
	
	######## PROBS
	# probs is now (id, label, class, probability, group)
	probs  <- cbind(probs[,1:2], maxval, maxp, groups)
	
	
	############# Getting info on labeling
	
	# data is 3 column: id, label, group
	data <- dbGetQuery(db$db, paste("select id, label from masks where label in (", labelstring, ")"))
	data <- cbind(data, rep(0, nrow(data)))

	# Assign clique group numbers
	for(l in labels){
		mval <- max(data[,3], na.rm=T)
		w <- which(data[,2]==l)
		data[w,3] <- GetCliques(db, data[w,1], mval+1)
	}

	######### SEG
	# seg is a matrix giving the id and segmentation for all segmented data
	seg <- dbGetQuery(db$db, paste("select id, segmentation from masks where segmentation in (", labelstring, ")"))

	
	# for each annotated group, assign to appropriate class and remove mask from seg (if it exists)
	# check for exact matches first so that won't get removed as a marginal match for a different ROI
	rois <- unique(data[,3], na.rm=T)
	for(roi in rois){
		# the masks that were marked as ok for this ROI
		maskids <- data[which(data[,3]==roi),1]
		
		# the label for the ROI
		lab <- data[which(data[,1]==maskids[1]), 2]

		w = which(seg[,1] %in% maskids)
		if(length(w)==1){
			# there was one matching mask.  this is correct 
			classgrp <- probs[which(probs[,1]==seg[w,1]),5]
			whichclass <- which(probs[,5]==classgrp)
			ret <- rbind(ret, c(lab, 1, length(whichclass), mean(probs[whichclass, 4]), min(probs[whichclass, 4]), max(probs[whichclass, 4])))
			seg <- seg[-w,]
			rois <- rois[-which(rois==roi)]
		}
	}
	
	while(length(rois>=1)){
		for(roi in rois){
			# there was no matching mask - check for overlapping mask
			maskids <- data[which(data[,3]==roi),1]
			
			# the label for the ROI
			lab <- data[which(data[,1]==maskids[1]), 2]
			labrow <- which(ret[,1]==lab)
			# this is either a missed ROI or a marginal ROI
			# check if it is a marginal ROI by finding all masks in seg which overlap with
			#  any of the masks for the ROI
			ovlp <- as.matrix(dbGetQuery(db$db, paste("select maskid1, maskid2 from edges where (maskid1 in (",
				paste(maskids, collapse=", "), ") and maskid2 in (", paste(seg[,1], collapse=", "), 
				")) or (maskid1 in (", paste(seg[,1], collapse=", "), ") and maskid2 in (", 
				paste(maskids, collapse=", "),"))")))
					
			#print(ovlp)
			# if there is any overlap, this ROI has at least one marginal match
			if(nrow(ovlp)>=1){
				# which of the segmented masks overlap with this ROI
				segovlprows <- which(seg[,1] %in% as.vector(ovlp))
				if(length(segovlprows)==1){
					# there was one marginal mask
					maskid <- seg[segovlprows, 1]
					classgrp <- probs[which(probs[,1]==maskid),5]
					whichclass <- which(probs[,5]==classgrp)
					ret <- rbind(ret, c(seg[segovlprows,2], 2, length(whichclass), mean(probs[whichclass, 4]), min(probs[whichclass, 4]), max(probs[whichclass, 4])))					
					seg <- seg[-segovlprows,]						
					rois <- rois[-which(rois==roi)]
				}else{					
					# at least one matched the ROI label.  count the match and remove the ROI
					maskid <- seg[segovlprows[1],1]
					classgrp <- probs[which(probs[,1]==maskid),5]
					whichclass <- which(probs[,5]==classgrp)
					ret <- rbind(ret, c(seg[segovlprows[1],2], 2, length(whichclass), mean(probs[whichclass, 4]), min(probs[whichclass, 4]), max(probs[whichclass, 4])))					
					rois <- rois[-which(rois==roi)]
					seg <- seg[-segovlprows[which(seg[segovlprows,2]==lab)[1]],]
				}
			}else{
				# Missed ROI -- not doing anything with this now
				rois <- rois[-which(rois==roi)]
			}
		}
	}

	# Any remaining masks in seg are false positives
	for(rnum in 1:nrow(seg)){
		maskid <- seg[rnum, 1]
		classgrp <- probs[which(probs[,1]==maskid),5]
		whichclass <- which(probs[,5]==classgrp)
		ret <- rbind(ret, c(lab, 3, length(whichclass), mean(probs[whichclass, 4]), min(probs[whichclass, 4]), max(probs[whichclass, 4])))					
	}
	
	ret <- as.data.frame(ret[-1,])
	names(ret) <- c("class", "label", "nmasks", "mean prob", "min prob", "max prob")
	return(ret)
}


#-
#' Assigns each of the specified masks to a clique to use when
#' solving the MWIS segmentation
#' 
#' @param db the mask database
#' @param ids the vector of ids of the masks to assign to cliques
#' @param minid the minimum value of the clique ids to return
#' 
#' @return a vector giving the clique id for each of the specified masks
#-
GetCliques <- function(db, ids, minid=1){
	# returns a vector where groups each have a unique integer id >0
	# returns in the same order as the ids are passed into the function
	# minid is the starting id assigned to groups.  will increment from there.

	# note: there's a problem where 4 masks with the overlap characteristics 1-2-3-4 get assigned clique
	# IDs 1-2-2-3 which results in an infinite loop in the segmentation algorithm.  should be assigned
	# IDs 1-1-2-2.  Want to encourage smaller number of cliques?  Maximin clique size?
	# current hack - remove any element that is in more than one clique, because these shouldn't be chosen anyway.
	#    set the return value to 0 for these elements

	idsstring <- paste(ids, collapse=", ")
	qry <- paste("select maskid1, maskid2 from edges where maskid1 IN (",idsstring, ") AND maskid2 IN (", idsstring,")")
	cons <- dbGetQuery(db$db, qry)

	cons <- as.matrix(cons)
	cons <- apply(cons, c(1,2), toString)
	graph <- graph.edgelist(cons, directed=F)


	# cliques is list of maximal cliques with the graph vertex ids
	# we need to convert these back to mask ids
	cliques <- maximal.cliques(graph)
	vertexids <- as.integer(V(graph)$name)
	
	ret <- rep(NA, length(ids))
	
	for(i in 1:length(ids)){
		# For each mask, if it is in one group, give it that number.  If it 
		# is in multiple groups, keep it at NA
		id = ids[i]
		if(!(id %in% vertexids)){
			ret[i]=0
		}else{		
			clid <- which(vertexids==id)-1
			clq <- which(sapply(cliques, function(list, val){return(val %in% list)}, val=clid))  
			if(length(clq)==1){
				ret[i]=clq
			}
		}
	}

	if(all(is.na(ret))){
		mval <- 2
	}else{
		mval <- max(ret, na.rm=T)
		if(any(ret[-which(is.na(ret))]==0)){
			ret[which(ret==0)] = seq(mval+1, mval+length(which(ret==0)))
		}
	}
	
	ret= ret+minid
	return(ret)
}

# GetConnectedGroups <- function(db, ids, minid=1){
# 	idsstring <- paste(ids, collapse=", ")
# 	qry <- paste("select maskid1, maskid2 from edges where maskid1 IN (",idsstring, ") AND maskid2 IN (", idsstring,")")
# 	cons <- dbGetQuery(db$db, qry)
# 
# 	cons <- as.matrix(cons)
# 	cons <- apply(cons, c(1,2), toString)
# 	graph <- graph.edgelist(cons, directed=F)
# 
# 	members <- clusters(graph)$membership + minid
# 	vertexids <- as.integer(V(graph)$name)
# 	o <- order(ids)
# 	o2  <- order(vertexids)
# 	
# 	return((members[o2])[o])
# 	
# }

#-
#' Checks each database in the given directory and pulls any
#' labeled data
#' 
#' @param directory the directory to pull data from
#' 
#' @return a data.frame with the extracted data
#' 
#' @export
#-
PullAllData <- function(directory){
	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
	
	count <- 1
	# Get first experiment with data
	while(count <= length(dblist)){
		tmp <- PullData(ConMaskDb(paste(directory, dblist[count], sep="")))
		if(nrow(tmp)>0){
			break
		}
		count <- count+1
	}
	if(nrow(tmp)==0){
		# If never found any experiment with data
		stop("No data found.")
	}else{
		# Else add the relevant experiment name
		namecol <- as.data.frame(rep(dblist[count], nrow(tmp)))
		names(namecol) <- c("experiment")
		ret <- cbind(namecol, tmp)
	}
	# Iterate through the remaining experiments and add any data
	count <-  count+1
	for(i in count:length(dblist)){
		tmp <- PullData(ConMaskDb(paste(directory, dblist[i], sep="")))
		if(nrow(tmp)>0){
			if( ncol(tmp) == (ncol(ret)-1) ){
				namecol <- as.data.frame(rep(dblist[i], nrow(tmp)))
				names(namecol) <- c("experiment")
				tmp <- cbind(namecol, tmp)
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
PullData <- function(db, labeled=TRUE, group=TRUE){
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
		stop("No data found in the supplied database.")
	}
}

#-
#' Get a matrix giving the segmentation stored in a mask database
#' 
#' @param db the mask database
#' @param classids a vector giving the ids of the classes to include in the returned segmentation
#' @param val if "id", puts the mask ids in the ROI in the returned matrix, otherwise uses the class id
#' 
#' @return a data.frame with the extracted data
#' 
#' @export
#-
GetSegmentation <- function(db, classids=c(2,3), val="id"){
	maskids <- dbGetQuery(db$db, paste("select id, segmentation from masks where segmentation in (",paste(classids, collapse=", "),")"))
	dims <- as.integer(dbGetQuery(db$db, "select nx, ny from experiment"))
	ret <- matrix(NA, dims[1], dims[2])
	for(id in maskids[,1]){
		if(val=="id"){
			ret[GetMask(db, id)]=id
		}else{
			ret[GetMask(db, id)]=maskids[which(maskids[,1]==id),2]
		}
	}
	return(ret)
}
