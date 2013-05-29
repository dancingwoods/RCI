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

# tmp function for evaluation
EvaluateBulk <- function(directory, source=T){
	# for each evaluation object in the target directory, add a row to the matrix
	# where the columns are
	# correct, marginal, missed, new (summed astrocytes and neurons)

	if(!source){
		filelist <- list.files(directory)
		load(paste(directory, filelist[1], sep=""))
		name <- paste(strsplit(filelist[1], split=".R")[[1]], sep="")
		perf <- get(name)
		for(i in 2:length(filelist)){
			load(paste(directory, filelist[i], sep=""))
			name <- paste(strsplit(filelist[i], split=".R")[[1]], sep="")
			perf <- perf + get(name)
			# ret[i,5] <- perf[2,2]-perf[2,3]-perf[2,5] 
			# ret[i,6] <- perf[2,5]
			# ret[i,7] <- perf[2,3]
			# ret[i,8] <- perf[2,4]
		}
		return(perf)
	}else{
		filelist <- list.files(directory)
		load(paste(directory, filelist[1], sep=""))
		name <- paste(strsplit(filelist[1], split=".R")[[1]], sep="")
		ret <- matrix(NA, length(filelist), 4)
		for(i in 1:length(filelist)){
			load(paste(directory, filelist[i], sep=""))
			name <- paste("perf_", strsplit(filelist[i], split=".Rdata")[[1]], sep="")
			perf <- get(name)
			ret[i,1] <- perf[1,2]-perf[1,3]-perf[1,5] 
			ret[i,2] <- perf[1,5]
			ret[i,3] <- perf[1,3]
			ret[i,4] <- perf[1,4] 
		}
		return(ret)
	}
}

EvaluateROC <- function(directory, perfrow=1){
	filelist <- list.files(directory)
	filelist <- filelist[grep("perf", filelist)]
	ret <- matrix(NA, length(filelist), 5)
	for(i in 1:length(filelist)){
		load(paste(directory, filelist[i], sep=""))
		name <- paste(strsplit(filelist[i], split=".Rdata")[[1]], sep="")
		perf <- get(name)
		ret[i,1] <- perf[perfrow,2]-perf[perfrow,3]-perf[perfrow,5] 
		ret[i,2] <- perf[perfrow,5]
		ret[i,3] <- perf[perfrow,3]
		ret[i,4] <- perf[perfrow,4] 
		ret[i,5] <- as.double(tail(strsplit(name, "_")[[1]], 1))
	}
	return(ret)
}

EvaluateDataLimited <- function(directory, outdirectory){
	data <- PullAllData(directory)
	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
	nexp <- length(dblist)
	for(n in 1:(nexp-1)){
		cat(n, " experiments\n")
		# For each number of training experiments, compute performance
		# with all possible sets of training data of that size

		# all combinations of training data of this size, each column is
		# the names of the selected experiments
		combos <- combn(dblist, n)
		for(i in 1:ncol(combos)){
			cat(i, " of ", ncol(combos), " training combos\n")
			class <- InitiateMaskClassifier(data[which(data[,1] %in% combos[,i]),])
			strclass <- paste(strsplit(combos[,i], split=".sqlite"), collapse="")
			
			for(d in dblist){
				cat(d, " ")
				# for each experiment not used to create classifier, test performance
				if(!(d %in% combos[,i])){
					db <- ConMaskDb(paste(directory, d, sep=""))
					PredictExperiment(class, db)
					dbstr <- strsplit(d, split=".sqlite")[[1]]
					assign(paste("perf_", dbstr, "_", strclass,sep=""), EvaluateSegmentation(db))
					save(list=paste("perf_", dbstr, "_", strclass,sep=""), file=paste(outdirectory,"perf_", dbstr, "_", strclass,sep=""))					
				}
			}
		}	
	}	
}

EvaluateCV <- function(directory, threshvals=seq(0.4,1,0.005), sources=F){
	data <- PullAllData(directory)
	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
	for(dbname in dblist){
		print(dbname)
		dbstr <- strsplit(dbname, split=".sqlite")[[1]]
		if(sources){
			scales <- dbGetQuery(db$db, "select id, tag from features")
			scales <- scales[grep("Scale", scales[,2]),]
			for(i in 1:nrow(scales)){
				assign(paste(dbstr, "_", scales[i,2], sep=""), EvaluateSegmentation(db, scales[i,1]))
				save(list=paste(dbstr, "_", scales[i,2], sep=""), file=paste("performance/sourcelevels/", dbstr, "_", scales[i,2],sep=""))
			}
		}else{
			class <- InitiateMaskClassifier(data[-which(data[,1]==dbname),])
			db <- ConMaskDb(paste(directory, dbname, sep=""))
			for(t in threshvals){
				cat(t, " ")
				PredictExperiment(class, db, thresh=t)
				assign(paste("perf_", dbstr, "_", t,sep=""), EvaluateSegmentation(db))
				save(list=paste("perf_", dbstr, "_", t,sep=""), file=paste("../performance/perf_", dbstr, "_", t,sep=""))
			}
		}	
	}
}




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
PredictExperiment <- function(classifier, db, thresh=NULL, enforcelabels=F){
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
			groups[w] <- GetCliques(db, as.vector(probs[w,1]), max(groups)+1)
		}
	}

	# results vector to store predicted class
	pvec  <- rep(NA, nrow(probs))
	pvec[which(maxval==1)]=1
	pvec[which(is.na(groups))]=1
	uniquegroups <- unique(groups[-which(is.na(groups))])

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
	while(nrow(curovlp)>0){
		print(curovlp)
		ind1 <- which(probs[,1]==curovlp[1,1])
		ind2 <- which(probs[,1]==curovlp[1,2])

		grp1 <- groups[ind1]
		grp2 <- groups[ind2]
		cat(grp1, grp2,"\n")
		
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
				if(probs[i,2+class1] + probs[j,2+class2] > bestval){
					bestind1 <- i
					bestind2 <- j
				}
			
			}
		}
		
		cat(probs[bestind1,1], probs[bestind2,1], "\n")
		
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
	}

	
	addseg <- cbind(pvec, probs[,1])
	addseg <- as.data.frame(addseg)
	names(addseg) <- c("segmentation", "id")
	
	dbBeginTransaction(db$db)
	dbGetPreparedQuery(db$db, "UPDATE masks SET segmentation=? where id=?", bind.data=addseg)
	dbCommit(db$db)
}

# Takes the current labels and segmentation in the given database and evaluates perfomance
# of the segmenter against the hand labels.
# For each group of labels cells, a false negative if no mask overlaps, false positive if segmented
#  mask overlaps with no group, marginal result if segmented cell overlaps with but is not of of the 
#  labeled masks. 
# Returns a matrix with a row for each class
# each row is 4 columns (class id, no. label groups, false negatives, false positives, marginal results)
EvaluateSegmentation <- function(db, sourceids=NULL){
	labels <- dbGetQuery(db$db, "select distinct label from masks")[,1]
	labels <- labels[which(labels>1)]
	
	if(!is.null(sourceids)){
		labelstring=paste(labels, collapse=", ")
		labels=c(1)
	}
	
	ret <- matrix(0, length(labels), 5)
	count <- 1
	
	for(label in labels){
		ret[count,1] <- label
		
		if(is.null(sourceids)){
			ids <- dbGetQuery(db$db, paste("select id from masks where label=", label))[,1]
		}else{
			ids <- dbGetQuery(db$db, paste("select id from masks where label in (", labelstring, ")"))[,1]
		}
		
		groups <- GetCliques(db, ids)
		labdata <- cbind(ids, groups)
		
		# number of cells hand labeled with this class - number of groups + number of singletons
		ret[count,2] <- length(unique(groups[which(groups>0)]))+length(which(groups==0))

		# Get ID of segmented masks - these will each be in their own group due to 
		# non-overlap being enforced during segmentation.
		if(is.null(sourceids)){
			idseg <- dbGetQuery(db$db, paste("select id from masks where segmentation=", label))[,1]
		}else{
			idseg <- dbGetQuery(db$db, 
				paste("select masks.id from masks join feats_masks on masks.id=feats_masks.maskid",
				"where feats_masks.featureid in (", paste(sourceids, collapse=", "), ")"))[,1]
		}
		
		
		# segmented masks which are also labeled with the class. these are correct
		matchedidseg <- which(idseg %in% ids)
		matchedgrps <- groups[which(ids %in% idseg)]
		
		# possible border case: what if there are two masks segmented that both are part of one labeled
		#  group -- no, this can't happen because a group is a clique
		# Could have a segmented mask, plus one that partially overlaps the group (but not the segmented
		#  mask).  This currently will be a false positive, not a marginal
		
		# remove any groups from the labeled data which have been matched
		if(length(matchedgrps)>0){
			labdata <- labdata[-which(groups %in% matchedgrps),,drop=F]
		}
		# remove the segmented masks which have been matched
		if(length(matchedidseg)>0){
			idseg <- idseg[-matchedidseg]
		}
		
		# groups that are left are either false negatives or marginal matches
		# idsegs that are left are either false positives or marginal matches
		# find marginal matches
		if(nrow(labdata)>0){
			
			for(id in idseg){
				
				idstring <- paste(labdata[,1], collapse=", ")
						
				ovlp <- dbGetQuery(db$db, paste("select maskid1, maskid2 from edges where (maskid1=",id,
				"and maskid2 in (", idstring, ")) or (maskid2=",id,"and maskid1 in (", idstring,"))"))
			
				if(nrow(ovlp)==0){
					# false positive
					ret[count,4] <- ret[count,4]+1
				}else{
	 				# marginal match
					ret[count,5] = ret[count,5]+1
					# remove marginal match groups 
					mgroups <- labdata[which(labdata[,1] %in% as.vector(as.matrix((ovlp)))),2]
					labdata <- labdata[-which(labdata[,2] %in% mgroups),]
					if(length(labdata)==2){
						labdata <- matrix(labdata, 1, 2)
					}
				}
			}
			# false negatives are those groups which were not partially or fully matched
			ret[count,3] <- length(unique(labdata[,2]))

		}else{
			ret[count,4] <- length(idseg)
		}

		# increment class index
		count <- count+1
	}
	ret <- as.data.frame(ret)
	names(ret) <- c("class", "labeled", "falseneg", "falsepos", "marginal")
	return(ret[order(ret[,1]),])
}


# returns a vector where groups each have a unique integer id >0
# returns in the same order as the ids are passed into the function
# minid is the starting id assigned to groups.  will increment from there.

# note: there's a problem where 4 masks with the overlap characteristics 1-2-3-4 get assigned clique
# IDs 1-2-2-3 which results in an infinite loop in the segmentation algorithm.  should be assigned
# IDs 1-1-2-2.  Want to encourage smaller number of cliques?  Maximin clique size?
# current hack - remove any element that is in more than one clique, because these shouldn't be chosen anyway.
#    set the return value to 0 for these elements
GetCliques <- function(db, ids, minid=1){
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
	
	print(vertexids[order(vertexids)])
	print(ids[order(ids)])
	
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
	}
	ret[which(ret==0)] = seq(mval+1, mval+length(which(ret==0)))

	ret= ret+minid
	return(ret)
}

# Checks each database in the given directory and pulls any
# labeled data
PullAllData <- function(directory){
	dblist <- list.files(directory)[grep(".*sqlite", list.files(directory))]
	
	count <- 1
	# Get first experiment with data
	while(count <= length(dblist)){
		tmp <- PullData(ConMaskDb(paste(directory, dblist[count], sep="")))
		if(nrow(tmp)>0){
			break
		}
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



