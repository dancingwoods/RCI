# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Segmentation of calcium images                                              
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



##########################################
## SQLite databases for mask management ##
##########################################

#-
#' Connects to an experiment's mask database
#' 
#' @param path the path to the SQLite database to connect to
#' 
#' @return a connection object as returned by dbConnect in the DBI package
#' 
#' @export
#-
ConMaskDb <- function(path){
	db <- dbConnect(dbDriver("SQLite"), dbname=path)
	ret = list("db"=db )
	class(ret) <- "MaskDb"
	return(ret)
}

#-
#' INTERNAL
#' Creates an empty mask database with the appropriate tables
#' 
#' @param db the database object for which to create the mask tables
#' 
#' @return NULL
#' 
#' @export
#-
MaskDbSetup <- function(db, calexp, tag, dt=NA, description=NA, datafile=NA){
	# 
	# TABLE: masks - constrained to be unique on mask
	#  id - primary key int
	#  mask - sparse representation of the mask (index of 1 pixels)
	#  label - integer, hand label.  0/NA=none, 1=not cell, 2=neuron, 3=astrocyte
	#  segementation - output of segmenter, 0/NA=none, 1=not cell, 2=neuron, 3=astrocyte
	# 
	# TABLE: experiment
	#  tag - text, shorthand name of experiment (ie 't782')
	#  description - text, description of experiment in long form giving parameters
	#  nx - int, x dimension of images
	#  ny - int, y dimension of images
	#  nt - int number of images
	#  nc - int number of channels
	#  dt - 1/sampling rate of experiment
	#  filename - text filename of actual data for this experiment
	# 
	# TABLE: features
	#  id - primary key int
	#  tag - text short name for feature
	#  desc - long form description of feature
	# 
	# TABLE: feats_masks (linking table between masks and features)
	#   id - primary key int
	#   maskid - primary key of mask
	#   featureid - primary key of feature
	#   fvalue - real value of the feature for this mask
	# 
	# TABLE: edges (has entry for masks which overlap and therefore can't both be 1)
	#   id - primary key int
	#   maskid1 - id of first mask
	#   maskid2 - id of second mask
	# 
	# TABLE: images
	#   id - primary key int
	#   tag - name of summary image
	#   image - text string of image
	
	if(is.null(tag)){
		cat("error: no tag provided\n")
		return()
	}
	if(is.null(calexp)){
		cat("error: no calexp object provided\n")
		return()
	}

	dbGetQuery(db$db, "CREATE TABLE masks (id INTEGER PRIMARY KEY, mask TEXT NOT NULL, label INTEGER, segmentation INTEGER, UNIQUE(mask))")
	dbGetQuery(db$db, "CREATE TABLE experiment (tag TEXT NOT NULL, description TEXT, nx INTEGER NOT NULL, ny INTEGER NOT NULL, nt INTEGER,nc INTEGER, dt REAL, filename TEXT, UNIQUE(tag))")
	dbGetQuery(db$db, "CREATE TABLE features (id INTEGER PRIMARY KEY, tag TEXT NOT NULL, description TEXT, UNIQUE(tag))")
	dbGetQuery(db$db, "CREATE TABLE feats_masks (id INTEGER PRIMARY KEY, maskid INTEGER NOT NULL, featureid INTEGER NOT NULL, fvalue REAL NOT NULL, UNIQUE(maskid, featureid))")
	dbGetQuery(db$db, "CREATE TABLE edges (id INTEGER PRIMARY KEY, maskid1 INTEGER NOT NULL, maskid2 INTEGER NOT NULL)")
	dbGetQuery(db$db, "CREATE TABLE images (id INTEGER PRIMARY KEY, tag TEXT NOT NULL, image TEXT NOT NULL)")
	if(is.na(dt)){
		dt <- 0
	}
	if(is.na(description)){
		description <- "(unknown)"
	}
	if(is.na(datafile)){
		datafile <- "(unknown)"
	}

	d <- dim(calexp$data)
	qry <- paste("INSERT INTO experiment (tag, description, nx, ny, nt, nc, dt, filename) values ('",tag,"', '",description,"', ", 
		d[3],", ",d[4],", ", d[2], ", ", d[1], ", ", dt,", '", datafile, "')", sep="")
	dbGetQuery(db$db, qry)
}

#-
#' Summary function to quickly see the details/statistics of a mask database object
#' 
#' @param db the mask database option
#' @param flag type of summary to give
#' 
#' @return NULL
#' 
#' @export
#-
summary.MaskDb <- function(db, flag="default"){
	nmasks <- as.integer(dbGetQuery(db$db, "SELECT count(*) from masks"))

	features <- dbGetQuery(db$db, "SELECT id, tag FROM features")

	# Get information about the sources
	sourceinds <- grep("Source_", features[,2])
	sources <- as.vector(sapply(features[sourceinds,2], substr, start=8, stop=100))
		
	# Keep track of features and IDs associated with sources
	sourcefeatures <- features[sourceinds,]
	
	# Information about the remaining non-source features
	features <- features[-sourceinds,]

	# How many masks from each source
	sourcecounts <- rep(NA, length(sources))
	# How many scales for each source
	scalecounts <- rep(NA, length(sources))
	for(i in 1:length(sources)){
		fid <- dbGetQuery(db$db, paste("SELECT id FROM features WHERE tag='Source_", sources[i], "'", sep=""))
		sourcecounts[i] <- nrow(dbGetQuery(db$db, paste("SELECT fvalue FROM feats_masks WHERE featureid=",fid , sep="")))
		scaleids <- grep(paste("Scale_",sources[i], "_[0-9].*", sep=""), features[,2])
		scalecounts[i] <- length(scaleids)
		features <- features[-scaleids,]
	}

	cat("\nExperiment:", dbGetQuery(db$db, "SELECT tag FROM experiment")[[1]], " \n\n")
	cat("There are", nmasks, "unique masks in the database.\n\n")

	cat("Mask sources:\n")
	for(i in 1:length(sources)){
		cat(sources[i], "  (",sourcecounts[i]," unique masks, ",scalecounts[i], " scales)\n", sep="")
	}
	cat("\n\n")

	# For each source, get information about the scales
	cat("Features:\n")
	if(length(features)==0){
		cat("(no features yet)\n")
	}else{
		for(i in 1:nrow(features)){
			cat(features[i,2], "( ID:",features[i,1],"; ", as.integer(dbGetQuery(db$db, paste("SELECT count(*) from feats_masks where featureid=", features[i,1], sep="")))," masks) \n")
		}
	}
	cat("\n")

	
}




#######################################################################
############### Functions to add masks to a database ##################
#######################################################################

#-
#' Converts a matrix with unique positive integers on each mask into a sparse mask list
#' 
#' @param maskmat the matrix with unique integers on each mask
#' 
#' @param return a sparsemasks object with a list of masks in the $masks element and the 
#' dimensions of the original image in the $dims element
#' 
#' @export
#-
MatrixToSparseMasks <- function(maskmat){
	ids <- unique(as.vector(maskmat))
	ids <- ids[which(ids>0)]
	masklist <- sapply(ids, function(id, img){return(which(img==id))}, img=maskmat)
	ret <- list("masks"=masklist, "dims"=dim(maskmat))
	class(ret) <- c("sparsemasks")
	return(ret)
}

#-
#' Generates a set of masks using the Laplacian of Gaussian technique
#' for the given scale and kernel size
#' 
#' @param image the image to use to generate masks
#' @param scale the scale of the kernel to use
#' @param kside the size of the kernel
#' @param sparse boolean, should the function return sparse masks instead of a matrix for plotting
#' 
#' @return a matrix with unique integers at mask locations and 0 in the background, or a a sparse masks object
#' 
#' @export
#-
LoGMasks <- function(image, scale, ksize=15, sparse=T){
	log.image <- ConvolveImage(image, LoGKernel(ksize,scale))
	region <- matrix(0, nrow(log.image), ncol(log.image))
	region[which(log.image<0)]=1
	masks <- AssignToPeaks(region, image, restrict=T) 
	if(!sparse){
		return(masks)
	}else{
		return(MatrixToSparseMasks(masks))
	}
}

#-
#' Generates a set of masks using thresholding of the sliding histogram equalized
#' version of an image
#' 
#' @param image the image to use to generate masks
#' @param thresh the threshold to use (pixels above thresh in equalized image are found)
#' @param radius the radius for the window used for equalization
#' @param fullmax the maximum value possible (for equalization)
#' @param sparse boolean, should the function return sparse masks instead of a matrix for plotting
#' 
#' @return a matrix with unique integers at mask locations and 0 in the background, or a a sparse masks object
#' 
#' @export
#-
EqualThreshMasks <- function(image, thresh, radius=8, fullmax=4096, sparse=T){
	eq.image <- SlidingHistEqualC(image, radius, fullmax)
	region <- matrix(0, nrow(eq.image), ncol(eq.image))
	region[which(eq.image>thresh)] <- 1
	masks <- AssignToPeaks(region, image, restrict=T)
	if(!sparse){
		return(masks)
	}else{
		return(MatrixToSparseMasks(masks))
	}
}


# TODO: there's a warning about closing with pending results set - but seems to work.  Should find
# the error sometime.
# TODO: This could almost cerainly be sped up by clever SQLite manipulations
#-
#' Generate masks according to the given method and add them to the database
#' 
#' @param db the mask database object
#' @param calexp the calcium experiment data object
#' @param method what method should be used to generate masks to add. 'LoG' Laplacian of
#' Gaussian, EqThresh thresholding of equalized image
#' @param channel which data channel to use for computing the masks
#' @param scales which smoothing scales to use for method (LoG)
#' @param invert boolean, should the mean image be inverted to find dark regions instead of bright?
#' 
#' @export
#-
AddMasks <- function(db, calexp, method, channel=2, scales=NULL, invert=F){
	mimg <- apply(calexp$data[channel,,,], c(2,3), mean)
	if(invert){
		mimg <- -1*mimg
		mimg <- mimg+min(mimg)
	}
	masks <- c()
	
	if(is.null(scales)){
		if(method=="LoG"){
			scales <- LogSeq(0,5,50)
		}else if(method=="EqThresh"){
			scales <- seq(2500, 4000, length.out=100)
		}
	}

	if(method=="LoG"){
		# Laplacian of Gaussian mask finding procedure
		masks <- as.vector(sapply(LoGMasks(mimg, scales[1])$masks, toString))
		scalevec <- rep(scales[1], length(masks))
		if(length(scales)>1){
			for(i in 2:length(scales)){
				newmasks <- as.vector(sapply(LoGMasks(mimg, scales[i])$masks, toString))
				masks <- c(masks, newmasks)
				scalevec <- c(scalevec, rep(scales[i], length(newmasks)))
			}
		}
	}else if(method=="EqThresh"){
		# Thresholding of sliding histogram equalized mean image
		masks <- as.vector(sapply(EqualThreshMasks(mimg, scales[1])$masks, toString))
		scalevec <- rep(scales[1], length(masks))
		if(length(scales)>1){
			for(i in 2:length(scales)){
				newmasks <- as.vector(sapply(EqualThreshMasks(mimg, scales[i])$masks, toString))
				masks <- c(masks, newmasks)
				scalevec <- c(scalevec, rep(scales[i], length(newmasks)))
			}
		}
	}else{
		cat("Invalid method selection\n")
		return()
	}
	
	######
	# Add the vector of (string) sparse masks to the database
	cat("Adding ", length(masks), " masks (", length(unique(masks))," unique)\n")
	methodtag <- paste(method, "_ch", channel, sep="")
	if(invert){
		methodtag <- paste(methodtag, "_inverse", sep="")
	}
	for(i in 1:length(masks)){
		if(i%%1000==0){
			cat(i, " ")
		}
		# Get the sourceid, or add the source if it doesn't exist
		qry <- paste("SELECT id FROM features WHERE tag='Source_",methodtag,"'", sep="")
		res <- dbGetQuery(db$db, qry)
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			qry <- paste("INSERT INTO features (tag) values ('Source_",methodtag,"')", sep="")
			dbGetQuery(db$db, qry)
			qry <- paste("SELECT id FROM features WHERE tag='Source_",methodtag,"'", sep="")
			res <- dbGetQuery(db$db, qry)
		}
		sourceid <- res[1,1]

		# Get the sourceid, or add the SCALE if it doesn't exist
		qry <- paste("SELECT id FROM features WHERE tag='Scale_", methodtag,"_", scalevec[i], "'", sep="")
		res <- dbGetQuery(db$db, qry)
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			qry <- paste("INSERT INTO features (tag) values ('Scale_", methodtag,"_", scalevec[i], "')", sep="")
			dbGetQuery(db$db, qry)
			qry <- paste("SELECT id FROM features WHERE tag='Scale_", methodtag,"_", scalevec[i], "'", sep="")
			res <- dbGetQuery(db$db, qry)
		}
		scaleid <- res[1,1]

		
		# Try adding the mask.  Duplicate masks won't be added
		tryCatch(dbGetQuery(db$db, paste("INSERT into masks (mask) values ('", masks[i], "')", sep="")), error=function(er){})

		# Find the maskid
		qry <- paste("SELECT id FROM masks WHERE mask = '", masks[i], "'", sep="")
		res <- dbGetQuery(db$db, qry)
		maskid <- res[1,1]


		# If there is already a link between this SOURCE and this mask, increment it by 1
		qry <- paste("SELECT id FROM feats_masks WHERE maskid = ", maskid, " AND featureid = ", sourceid, sep="")
		res <- dbGetQuery(db$db, qry)[[1]]

		if(length(res)==0){
			# There was no link (new source for this mask)
			dbGetQuery(db$db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (", maskid, ", ", sourceid, ", 1)", sep=""))
		}else{
			dbGetQuery(db$db, paste("UPDATE feats_masks SET fvalue=(fvalue + 1) WHERE id=", res[1], sep=""))
		}


		# If there is already a link between this SCALE and this mask, increment it by 1
		qry <- paste("SELECT id FROM feats_masks WHERE maskid = ", maskid, " AND featureid = ", scaleid, sep="")
		res <- dbGetQuery(db$db, qry)[[1]]

		if(length(res)==0){
			# There was no link (new scale for this mask)
			dbGetQuery(db$db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (", maskid, ", ", scaleid, ", 1)", sep=""))
		}else{
			dbGetQuery(db$db, paste("UPDATE feats_masks SET fvalue=(fvalue + 1) WHERE id=", res[1], sep=""))
		}

	}
	cat("\n")
}

AddImages <- function(db, calexp){
	mimg1 <- apply(calexp$data[1,,,], c(2,3), mean)
	mimg2 <- apply(calexp$data[2,,,], c(2,3), mean)
	mimg1eq <- SlidingHistEqualC(mimg1, 8)
	mimg2eq <- SlidingHistEqualC(mimg2, 8)

	dbGetQuery(db$db, paste("INSERT into images (tag, image) values ('mimg1', '", paste(as.vector(mimg1), collapse=", "),"')"))
	dbGetQuery(db$db, paste("INSERT into images (tag, image) values ('mimg1eq', '", paste(as.vector(mimg1eq), collapse=", "),"')"))
	dbGetQuery(db$db, paste("INSERT into images (tag, image) values ('mimg2', '", paste(as.vector(mimg2), collapse=", "),"')"))
	dbGetQuery(db$db, paste("INSERT into images (tag, image) values ('mimg2eq', '", paste(as.vector(mimg2eq), collapse=", "),"')"))

}


#####################################################################################
################ Functions to manipulate masks (extract, get features etc)  #########
#####################################################################################

GetImage <- function(db, tag){
	qry = paste("SELECT image from images where tag='", tag, "'", sep="")
	imgstr <- dbGetQuery(db$db, qry)[1,1]
	dims <- dbGetQuery(db$db, "SELECT nx, ny from experiment")
	vec <- as.double(strsplit(imgstr, split=", ")[[1]])
	return(matrix(vec, dims[1,1], dims[1,2]))
}

#-
#' Sets the label field for a particular mask in a mask database
#' 
#' @param db a database connection
#' @param id the id of the mask to label
#' @param label the label to assign to the mask (0=unknown, 1=cell, 2=not cell)
#' 
#' @return NULL
#' 
#' @export
#-
SetMaskLabel <- function(db, id, label){
	dbGetQuery(db$db, paste("UPDATE masks SET label=", label, " WHERE id= ", id, sep=""))
}

#-
#' Return the requested mask from the specified database
#' 
#' @param db a database connection
#' @param id the id of the mask to return
#' @param format "sparse" for a sparse mask in vector form, "matrix" for a matrix mask
#' 
#' @return either a vector giving the indices of the requested mask or a matrix version of
#' the mask
#' 
#' @export
#-
GetMask <- function(db, id, sparse=T){
	rawd <- dbGetQuery(db$db, paste("select mask from masks where id=", id, sep=""))
	smask <- as.integer(unlist(strsplit(rawd[1,1], split=" ")))
 	if(sparse){
		return(smask)
	}else{
		dims <- dbGetQuery(db, "select nx, ny from experiment")
		ret <- matrix(0, dims[1,1], dims[1,2])
		ret[smask] <- 1
		return(ret)		
	}
}

#-
#' Remove a mask and its feature links from the database
#' 
#' @param db the mask database object
#' @param id the id of the mask to remove
#' 
#' @export
#-
RemoveMask <- function(db, id){
	dbGetQuery(db$db, paste("DELETE FROM masks WHERE id=",id, sep=""))
	dbGetQuery(db$db, paste("DELETE FROM feats_masks WHERE maskid=",id, sep=""))
}

#-
#' Removes masks from the database whose size is less than min pixels or greater
#' than max pixels
#' 
#' @param db the mask database object
#' @param min the minimum mask size (pixels) to retain
#' @param max the maximum mask size (pixels) to retain
#' 
#' @export
#-
RestrictMaskSize <- function(db, minsize=NA, maxsize=NA){
	
	res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='Size'")
	if(length(res[[1]])==0){
		cat("Need to compute size feature first\n")
		return()
	}
	sizeid <- res[1,1]
	masks <- dbGetQuery(db$db, paste("SELECT maskid, fvalue FROM feats_masks where featureid=", sizeid,sep=""))
	if(!is.na(minsize)){
		w = which(masks[,2]<minsize)
		if(length(w)>0){
			for(i in 1:length(w)){
				RemoveMask(db, masks[w[i],1])
			}
			masks <- masks[-w,]
		}
	}
	if(!is.na(maxsize)){
		w = which(masks[,2]>maxsize)
		if(length(w)>0){
			for(i in 1:length(w)){
				RemoveMask(db, masks[w[i],1])
			}
			masks <- masks[-w,]
		}
	}
}

#-
#' Computed features of masks that are currently in the mask database.
#' This is sort of hack - should be more general and allow users to
#' specify feature functions
#' 
#' @param db the mask database object
#' @param calexp the calcium experiment object
#' @param feature string indicating which feature to compute
#' 
#' @details Current feature options: 
#' "size" - mask size in pixels
#' 
#' @export
#-
ComputeMaskFeatures <- function(db, calexp, feature){

	# Features computed and then added to the database with a prepared query
	masks <- dbGetQuery(db$db, "SELECT id, mask FROM masks")
	nr <- dbGetQuery(db$db, "SELECT ny FROM experiment")[1,1]
	nc <- dbGetQuery(db$db, "SELECT nx FROM experiment")[1,1]

	LimitMasks <- function(masks, fid){
		# given a set of mask ids, return those ids which already have a value for the specified feature
		havefeat <- dbGetQuery(db$db, paste("SELECT masks.id FROM ",
			"masks LEFT JOIN feats_masks ON masks.id=feats_masks.maskid",
			" WHERE feats_masks.featureid=", fid, sep=""))
		# Remove the masks which already have the feature
		w <- which(masks[,1] %in% havefeat)
		if(length(w)>0){
			masks <- masks[-w,]
		}
		return(masks)
	}

	## STRUCTURAL FEATURES (SHAPE-BASED)

	# SIZE
	if(feature=="size"){
		# Select masks from database that don't already have a size feature
		# Get the feature id for size
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='Size'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('Size')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='Size'")
			fid <- res[1,1]
		}else{
			fid <- res[1,1]
			masks <- LimitMasks(masks, fid)
		}
		
		# Get the size for the remaining masks
		GetSize <- function(string){
			return(length(strsplit(string, ", ")[[1]]))
		}
		size <- as.vector(sapply(masks[,2], GetSize))
		featdf <- cbind(rep(fid, length(size)), masks[,1], size)
		featdf <- as.data.frame(featdf)
		colnames(featdf) <- c("featureid", "maskid", "fvalue")
	}

	# CONVEX HULL PERCENTAGE	
	else if(feature=="chull"){
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='CHull'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('CHull')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='CHull'")
			fid <- res[1,1]
		}else{
			fid <- res[1,1]
			masks <- LimitMasks(masks, fid)
		}
		
		# Get the convex hull fraction for the masks
		GetHullPercentage <- function(string, nr, nc){
			vec <- as.integer(strsplit(string, split=", ")[[1]])
			hullsize <- length(MaskHull(vec, nr, nc))
			return(length(vec)/hullsize)
		}
		chullp <- as.vector(sapply(masks[,2], GetHullPercentage, nr=nr, nc=nc))
		featdf <- cbind(rep(fid, length(chullp)), masks[,1], chullp)
		featdf <- as.data.frame(featdf)
		colnames(featdf) <- c("featureid", "maskid", "fvalue")
	}
	
	# NUMBER OF HOLES
	else if(feature=="holes"){
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='Holes'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('Holes')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='Holes'")
			fid <- res[1,1]
		}else{
			fid <- res[1,1]
			masks <- LimitMasks(masks, fid)
		}
		
		# Get the convex hull fraction for the masks
		GetHoles <- function(string, nr, nc){
			vec <- as.integer(strsplit(string, split=", ")[[1]])
			return(CountHolesC(vec, nr, nc))
		}
		holes <- as.vector(sapply(masks[,2], GetHoles, nr=nr, nc=nc))
		featdf <- cbind(rep(fid, length(holes)), masks[,1], holes)
		featdf <- as.data.frame(featdf)
		colnames(featdf) <- c("featureid", "maskid", "fvalue")
	}
	
	# ACTIVITY FEATURES (FLUORESCENCE DATA BASED)
	
	# Means - for each channel, and for the equalized version
	else if(feature=="means"){
		GetMean <- function(string, mimg){
			vec <- as.integer(strsplit(string, split=", ")[[1]])
			return(mean(mimg[vec]))
		}
		
		mimg1 <- apply(calexp$data[1,,,], c(2,3), mean)
		mimg2 <- apply(calexp$data[2,,,], c(2,3), mean)
		mimg1eq <- SlidingHistEqualC(mimg1, 8, fullmax=4096)
		mimg2eq <- SlidingHistEqualC(mimg2, 8, fullmax=4096)
		
		# Removes any mask that has _any_ of the mean feautures already
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh1'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('MeanCh1')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh1'")
			mc1id <- res[1,1]
		}else{
			mc1id <- res[1,1]
			masks <- LimitMasks(masks, mc1id)
		}
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh2'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('MeanCh2')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh2'")
			mc2id <- res[1,1]
		}else{
			mc2id <- res[1,1]
			masks <- LimitMasks(masks, mc2id)
		}
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh1Eq'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('MeanCh1Eq')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh1Eq'")
			mc1eqid <- res[1,1]
		}else{
			mc1eqid <- res[1,1]
			masks <- LimitMasks(masks, mc1eqid)
		}
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh2Eq'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('MeanCh2Eq')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='MeanCh2Eq'")
			mc2eqid <- res[1,1]
		}else{
			mc2eqid <- res[1,1]
			masks <- LimitMasks(masks, mc2eqid)
		}
		
		mch1 <- as.vector(sapply(masks[,2], GetMean, mimg=mimg1))
		mch2 <- as.vector(sapply(masks[,2], GetMean, mimg=mimg2))
		mch1eq <- as.vector(sapply(masks[,2], GetMean, mimg=mimg1eq))
		mch2eq <- as.vector(sapply(masks[,2], GetMean, mimg=mimg2eq))

		nmasks <- nrow(masks)
		print(nmasks)
		print(c(mc1id, mc2id, mc1eqid, mc2eqid))
		featdf <- cbind(c(rep(mc1id, nmasks), rep(mc2id, nmasks), rep(mc1eqid, nmasks), rep(mc2eqid, nmasks)),
			rep(masks[,1], 4), c(mch1, mch2, mch1eq, mch2eq))
		featdf <- as.data.frame(featdf)
		colnames(featdf) <- c("featureid", "maskid", "fvalue")
	}
	
	else if(feature=="correlation"){
		GetCorrelation <- function(string, datamat){
			vec <- as.integer(strsplit(string, split=", ")[[1]])
			cormat <- cor(datamat[,vec])
			diag(cormat) <- NA
			return(mean(cormat, na.rm=T))
		}
		
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='CorCh1'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('CorCh1')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='CorCh1'")
			corch1id <- res[1,1]
		}else{
			corch1id <- res[1,1]
			masks <- LimitMasks(masks, corch1id)
		}
		res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='CorCh2'")
		if(length(res[[1]])==0){
			# Add the feature if it's not already in the database
			dbGetQuery(db$db, "INSERT INTO features (tag) values ('CorCh2')")
			res <- dbGetQuery(db$db, "SELECT id FROM features WHERE tag='CorCh2'")
			corch2id <- res[1,1]
		}else{
			corch2id <- res[1,1]
			masks <- LimitMasks(masks, corch2id)
		}
		
		ntime <- dim(calexp$data)[2]
		npixels <- prod(dim(calexp$data)[3:4])
		datamat1 <- matrix(calexp$data[1,,,], ntime, npixels)
		datamat2 <- matrix(calexp$data[2,,,], ntime, npixels)
		
		corch1 <- as.vector(sapply(masks[,2], GetCorrelation, datamat=datamat1))
		corch2 <- as.vector(sapply(masks[,2], GetCorrelation, datamat=datamat2))
		
		nmasks <- nrow(masks)
		featdf <- cbind(c(rep(corch1id, nmasks), rep(corch2id, nmasks)),
			rep(masks[,1], 2), c(corch1, corch2))
		featdf <- as.data.frame(featdf)
		colnames(featdf) <- c("featureid", "maskid", "fvalue")
		
	}
	else{
		cat("Error: Unrecognized feature name.\n")
		return()
	}
	
	# We have the featdf data.frame and the featid, now add the computed values to the database
	dbBeginTransaction(db$db)
	dbGetPreparedQuery(db$db, paste("INSERT INTO feats_masks (featureid, maskid, fvalue) VALUES (?,?,?)", sep=""), bind.data=featdf)
	dbCommit(db$db)
		
}


#-
#' INTERNAL
#' Compute the overlap matrix between a set of masks
#' 
#' @details Computes the overlap matrix between a set of masks using C
#' code for efficiency.
#' 
#' @param masklist A list of sparse masks.  Each element
#' of this list is a vector whose first element is the negative id of the mask
#' and whose other elements are the sorted indices of the mask pixels.
#' 
#' @return a matrix whose elements are 0 or 1 giving the overlap relationships
#' between the masks.  The masks are sorted in the matrix in the same order as they
#' are given in masklist.
#' 
#' @useDynLib RCI compmasksC
#-
CompMasksC <- function(masklist){
	nmasks <- length(masklist)

	# Use the C function
	out <- .C("compmasksC",
		as.integer(c(unlist(masklist), 0)),
		as.integer(nmasks),
		ret = as.logical(rep(FALSE, nmasks*nmasks))
	)
	return(matrix(out$ret, nmasks, nmasks))
}


#-
#' INTERNAL
#' Adds the overlap edges to a mask database
#' 
#' @param db the mask database object
#' 
#' @return NULL
#' 
#' @export
#-
AddMaskConnections <- function(db){

	# Remove any existing connections in the database
	dbGetQuery(db$db, "delete from edges")
	
	masks <- dbGetQuery(db$db, "select id, mask from masks order by id")
	ToVector <- function(row){
		return(c(-strtoi(row[1]), as.integer(strsplit(row[2], split=", ")[[1]])))
	}
	masklist <- apply(masks, 1, ToVector)
	ids <- strtoi(masks[,1])
	
	if(length(masklist) <= 46340){
		# few enough masks to compute overlap simultaneously.  need to rethink implementation
		# for scale
		cmat <- CompMasksC(masklist)
	
		# Make data.frame of connections to add
		cons <- arrayInd(which(cmat), dim(cmat))
		cons <- cons[which(cons[,1]>cons[,2]),]
		cons[,1] <- ids[cons[,1]]
		cons[,2] <- ids[cons[,2]]
		cons <- as.data.frame(cons)
		names(cons) <- c("maskid1", "maskid2")
		# Add edges from the data.frame using a prepared query
		dbBeginTransaction(db$db)
		dbGetPreparedQuery(db$db, "INSERT INTO edges (maskid1, maskid2) VALUES (?,?)", bind.data=cons)
		dbCommit(db$db)
	}else{
		for(i in 1:(length(masklist)-1)){
			id <- ids[i]
			tmpids <- ids[-(1:i)]
			cvec <- CompMaskC(masklist[[i]], masklist[-(1:i)])
			w <- which(cvec==1)
			cons <- cbind(rep(id, length(w)), tmpids[w])
			cons <- as.data.frame(cons)
			names(cons) <- c("maskid1", "maskid2")
			dbBeginTransaction(db$db)
			dbGetPreparedQuery(db$db, "INSERT INTO edges (maskid1, maskid2) VALUES (?,?)", bind.data=cons)
			dbCommit(db$db)
		}
	}
	
}


#-
#' Returns a list of the masks in a database
#' 
#' @param db a database connection
#' 
#' @return a list of vectors, each vector specifying a mask.  The first element of each mask
#' vector is the negative index of the mask.  The remaining elements of each vector are the
#' indices of the mask pixels.
#' @export
#-
GetSparseMasks <- function(db){
	ret <- dbGetQuery(db$db, "select id, mask from masks order by id")
	convertmask <- function(row){
		return(c(-strtoi(row[1]), as.integer(strsplit(row[2], split=", ")[[1]])))
	}
	return(apply(ret, 1, convertmask))
}


# ###########################
# ## Manipulation of masks ##
# ###########################

#-
#' INTERNAL
#' Compute the overlap of a single mask with a list of masks
#' 
#' @details Computes the overlap matrix between a mask and a list of other masks
#' using C code for efficiency
#' 
#' @param mask A single mask formatted as a vector of mask indices with or without
#' the negative id as the first element of the vector.
#' @param masklist A list of sparse masks as returned by GetMasks.  Each element
#' of this list is a vector whose first element is the negative id of the mask
#' and whose other elements are the sorted indices of the mask pixels.
#' 
#' @return a vector whose elements are 0 or 1 giving the overlap relationships
#' between the masks.  The values are sorted in the vector in the same order as they
#' are given in masklist.
#' 
#' @useDynLib RCI compmaskC
#-
CompMaskC <- function(mask, masklist){
	nmasks <- length(masklist)
	# Use the C function
	out <- .C("compmaskC",
		as.integer(c(mask[which(mask>0)], 0)),
		as.integer(c(unlist(masklist), 0)),
		as.integer(nmasks),
		ret = as.logical(rep(FALSE, nmasks))
	)
	return(out$ret)
}


#-
#' INTERNAL
#' Inverts a mask matrix so that the mask region is turned to backgroun
#' and vice versa
#' 
#' @param mask the mask matrix to invert, with NA in the background
#' 
#' @return a matrix with the inverted mask
#-
InvertMask <- function(mask){
	mat <- matrix(1, nrow(mask), ncol(mask))
	mat[which(mask==1)]=NA
	return(mat)
}


#-
#' INTERNAL
#' Converts a sparse mask to a matrix mask
#' 
#' @param sm the sparse representation of the mask (a vector whose positive values are the
#' indices of the mask pixels)
#' @param ny the number of rows of the matrix mask
#' @param nx the number of columbs of the matrix mask
#' @param background the value to put in the non-mask pixels of the matrix
#' 
#' @return A matrix of dimension (ny, nx) with 1's in the mask pixels and background
#' elsewhere
#-
SparseToMatrix <- function(sm, ny=126, nx=126, background=NA){
	mat <- matrix(background, ny, nx)
	mat[sm[which(sm>0)]] <- 1
	return(mat)
}


#-
#' INTERNAL
#' Converts a matrix to a sparse mask
#' 
#' @param mat the matrix to convert
#' 
#' @return A list of indices representing the mask
#-
MatrixToSparse <- function(mat){
	mat[which(is.na(mat))]=0
	return(which(mat>0))
}


#-
#' INTERNAL
#' Selects the masks from the given list that are contained in a region
#' 
#' @details Given a list of masks (with negative ids as first element) and a matrix with a mask specifying a
#' region, returns the masks in the masklist that are completely contained in the given region.
#' 
#' @param framemat a matrix of the same size as the masks in masklist with non-NA pixels specifying
#' the region in which to find masks
#' @param masklist a list of sparse masks (vectors where the first element
#' is the negative id of the mask and the other elements are the mask indices)
#-
GetInnerMasks <- function(framemat, masklist){
	tmp <- CompMaskC(MatrixToSparse(InvertMask(framemat)), masklist)
	return(masklist[which(!tmp)])
}
