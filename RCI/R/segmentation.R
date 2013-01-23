#############################################################
# Bronwyn Woods                                             #
# 2013                                                      #
#                                                           #
# Functions to handle segmentation of calcium images data.  #
# Currently the segmentation is only cell vs not cell.      #
#############################################################

#-
#' INTERNAL
#' Computes the sliding window histogram equalization of a matrix
#' 
#' @details Uses C code in slidinghistequalC.c
#' 
#' @param mat the matrix to equalize
#' @param radius the radius of the sliding window (total window size is a square window
#' with sides 2*radius+1)
#' @param fullmax the maximum value in the equalized image
#' 
#' @return The equalized matrix
#' 
#' @useDynLib RCI slidinghistequalC
#-
SlidingHistEqualC <- function(mat, radius, fullmax=4096){
	# Interfaces with the C function for sliding window histogram equalization
    d = dim(mat)
	
    # Use the C function
	out <- .C("slidinghistequalC",
		mat = as.double(mat),
		as.integer(d),
		as.integer(radius),
		as.integer(fullmax)
	)
		
	ret = matrix(out$mat, d[1], d[2])
	return(ret)
}

#-
#' INTERNAL
#' Computed the histogram equalization of a matrix.
#' 
#' @details Uses C code in histequalC.c
#' 
#' @param mat the matrix to equalize
#' @param fullmax the range to equalize to
#' 
#' @return the equalized matrix
#' @useDynLib RCI histequalC
#-
HistEqualC <- function(mat, fullmax=4096){
    d <- dim(mat)

    # Use the C function
	out <- .C("histequalC",
		as.double(mat),
		as.integer(d[1]),
		as.integer(d[2]),
		as.integer(fullmax),
		res = as.double(mat)
	)

	ret <- matrix(out$res, d[1], d[2])
	return(ret)
}

#-
#' INTERNAL
#' Finds the local maxima in an image
#' 
#' @details Uses C code in localmaxC.c
#' 
#' @param img the image in which to find the local maxima
#' @param min if this is set to 1, find local minima instead
#' 
#' @return matrix with 1 at the maxima and NA elsewhere
#' @useDynLib RCI localmaxC
#-
SimpleModesC <- function(img, min=0){
	# Interfaces with the C function for local max finding (local min if min==0)
	# the map flag determines whether background should be 0 or NA
    d = dim(img)
	
    # Use the C function
	out <- .C("localmaxC",
		mat = as.double(img),
		as.integer(d),
		as.integer(min)
	)
		
	ret <- matrix(out$mat, d[1], d[2])
	ret[which(ret==0)] <- NA
	return(ret)
}

#-
#' INTERNAL
#' Perform hill climbing on a matrix starting from a given point and returning the
#' local maxima that is reached.
#' 
#' @details Uses C code in hillclimbC.c
#' 
#' @param y Starting row
#' @param x Starting column
#' @param img The matrix on which to performt the hillclimbing
#' 
#' @return a vector of 2 numbers giving the coordinates of the peak found by hillclimbing
#' @useDynLib RCI hillclimbC
#-
HillClimbC <- function(y,x,img){
    d <- dim(img)
	
	out <- .C("hillclimbC",
		as.double(img),
		as.integer(d),
		x = as.integer(x),
		y = as.integer(y)
	)
		
	ret <- c(out$x, out$y)

	return(ret)
    
}


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
	return(db)
}


#######################################################################
############### Functions to set up and add to a database #############
#######################################################################

#-
#' INTERNAL
#' Creates an empty mask database with the appropriate tables
#' 
#' @param db the database object for which to create the mask tables
#' 
#' @return NULL
#-
DbSetup <- function(db){
	# TABLE: masks - constrained to be unique on (mask, expid)
	#  id - primary key int
	#  mask - text representation of the mask (sequence of 0 and 1)
	#  expid - integer what experiment (dataset) does this mask apply to
	#  label - integer, hand label.  0/NA=none, 1=cell, 2=not cell
	#  segementation - output of segmenter, 0/NA=not cell, 1=cell
	# 
	# TABLE: experiment
	#  id - primary key int
	#  tag - text, shorthand name of experiment (ie 't782')
	#  desc - text, description of experiment in long form giving parameters
	#  nx - int, x dimension of images
	#  ny - int y dimension of images
	#  nt - int number of images
	#  nc - int number of channels
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

	dbGetQuery(db, "CREATE TABLE masks (id INTEGER PRIMARY KEY, mask TEXT NOT NULL, label INTEGER, segmentation INTEGER, UNIQUE(mask))")
	dbGetQuery(db, "CREATE TABLE experiment (id INTEGER PRIMARY KEY, tag TEXT NOT NULL, description TEXT, nx INTEGER NOT NULL, ny INTEGER NOT NULL, nt INTEGER, nc INTEGER, filename TEXT, UNIQUE(tag))")
	dbGetQuery(db, "CREATE TABLE features (id INTEGER PRIMARY KEY, tag TEXT NOT NULL, description TEXT, UNIQUE(tag))")
	dbGetQuery(db, "CREATE TABLE feats_masks (id INTEGER PRIMARY KEY, maskid INTEGER NOT NULL,  featureid INTEGER NOT NULL, fvalue REAL NOT NULL, UNIQUE(maskid, featureid))")
	dbGetQuery(db, "CREATE TABLE edges (id INTEGER PRIMARY KEY, maskid1 INTEGER NOT NULL,  maskid2 INTEGER NOT NULL)")
}

#-
#' INTERNAL
#' Adds the overlap edges to a mask database
#' 
#' @param db a database connection object
#' @param cmat a matrix of 0/1 values giving the locations of edges between masks
#' @param ids a vector giving the ids of the masks in cmat (in order)
#' 
#' @return NULL
#-
AddConMat <- function(db, cmat, ids){
	# Make data.frame of connections to add
	cons <- arrayInd(which(cmat), dim(cmat))
	cons <- cons[which(cons[,1]>cons[,2]),]
	cons[,1] <- ids[cons[,1]]
	cons[,2] <- ids[cons[,2]]
	cons <- as.data.frame(cons)
	names(cons) <- c("maskid1", "maskid2")
	# Add edges from the data.frame using a prepared query
	dbBeginTransaction(db)
	dbGetPreparedQuery(db, "INSERT INTO edges (maskid1, maskid2) VALUES (?,?)", bind.data=cons)
	dbCommit(db)
}


#####################################################################################
################ Functions to manipulate masks (extract, get features etc)  #########
#####################################################################################

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
	dbGetQuery(db, paste("UPDATE masks SET label=", label, " WHERE id= ", id, sep=""))
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
GetMask <- function(db, id, format="sparse"){
	rawd <- dbGetQuery(db, paste("select mask, expid from masks where id=", id, sep=""))
	smask <- as.integer(unlist(strsplit(rawd[1,1], split=" ")))
	expid <- rawd[1,2]
	if(format=="sparse"){
		return(smask)
	}
	dims <- dbGetQuery(db, paste("select nx, ny from experiment where id= ", expid, sep=""))
	ret <- matrix(NA, dims[1,1], dims[1,2])
	ret[smask] <- 1
	return(ret)
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
GetMasks <- function(db){
	ret <- dbGetQuery(db, "select id, mask from masks order by id")
	convertmask <- function(row){
		return(c(-strtoi(row[1]), as.integer(unlist(strsplit(row[2], split=" ")))))
	}
	return(apply(ret, 1, convertmask))
}

#-
#' Remove a mask from a mask database
#' 
#' @details Removes a mask, as well as any associated features and edges
#' 
#' @param db a database connection
#' @param maskid the ID of the mask to remove
#' 
#' @return NULL
#' @export
#-
RemoveMask <- function(db, maskid){
	dbGetQuery(db, paste("DELETE from masks where id=",maskid, sep=""))
	dbGetQuery(db, paste("DELETE from feats_masks where maskid=",maskid, sep=""))
	dbGetQuery(db, paste("DELETE from edges where maskid1=",maskid," or maskid2=", maskid, sep=""))
}

#-
#' INTERNAL
#' Computes the features related to just the shape of masks, adding them to the database
#'
#' @details Computes features of all masks in the database or a list of masks specified by
#' id.  The features computed are currently:\cr
#' npixels - the number of pixels in a mask \cr
#' nholes - the number of non-mask pixels that are surrounded by at least 3 mask pixels\cr
#' bboxratio - the ratio of the area of the mask's bounding box and the number of pixels\cr
#' in the mask
#' hullratio - the ratio of the area of the mask's convex hull and the number of pixels
#' in the mask\cr
#' 
#' @param db a database connection
#' @param mids an optional vector of the mask ids for which to extract features
#' 
#' @return NULL
#-
# TODO: allow specification of which features to compute
GetShapeFeatures <- function(db, mids=NULL){
	# If mids is null, compute the features for all masks in the database
	if(is.null(mids)){
		mids <- dbGetQuery(db, "SELECT id from masks")[,1]
	}
	
	d <- dbGetQuery(db, "SELECT nx, ny from experiment")
	nrow <- d[1,2]
	ncol <- d[1,1]

	# Add the shape features to the database if they don't already exist
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('npixels', 'The number of pixels that are covered by this mask')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('nholes', 'The number of holes in this mask, defined as pixels that are not in the mask surrounded on 3 or more cardinal sides by pixels in the mask')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('bboxratio', 'The ratio of the number of pixels to the size of the rectangular bounding box')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('hullratio', 'The ratio of the number of pixels to the size of the convex hull of the mask')"))

	# Get the ids for the shape feature
	npixelsid <- dbGetQuery(db, "SELECT id from features where tag='npixels'")
	nholesid <- dbGetQuery(db, "SELECT id from features where tag='nholes'")
	bboxratioid <- dbGetQuery(db, "SELECT id from features where tag='bboxratio'")
	hullratioid <- dbGetQuery(db, "SELECT id from features where tag='hullratio'")

	mask <- matrix(NA, nrow, ncol)
	for(id in mids){
		mask <- GetMask(db, id)
		mask[which(is.na(mask))] <- FALSE
		mask[which(mask==1)] <- TRUE

		# Number of pixels
		numPixels <- length(which(!is.na(mask)))
	
		# Number of holes (defined as empty spaces with at least three map pixels on the
		# non-diagonal sides)
		numHoles <- countholesC(mask)
	
		# Percent of bounding box filled
		pixelInds <- arrayInd(which(!is.na(mask)), dim(mask))
		bboxarea <- (diff(range(pixelInds[,1]))+1) * (diff(range(pixelInds[,2]))+1)
		bBoxRatio <- numPixels/bboxarea
	
		# Percent of convex hull filled
		hull <- MaskHull(mask)
		hullRatio <- numPixels/(sum(hull, na.rm=T))
		foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", npixelsid, ", ", numPixels,")", sep="")))
		foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", nholesid, ", ", numHoles,")", sep="")))
		foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", bboxratioid, ", ", bBoxRatio,")", sep="")))
		foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", hullratioid, ", ", hullRatio,")", sep="")))		
	}

}

#-
#' INTERNAL
#' Computes the features related to the data under a mask, adding them to the database
#'
#' @details Computes features of all masks in the database or a list of masks specified by
#' id.  The features computed are currently:\cr
#' var1 - the variance of the pixel means for channel 1 \cr
#' var2 - the variance of the pixel means for channel 2\cr
#' var1eq - the variance of the pixel means for the equalized version of channel 1\cr
#' var2eq - the variance of the pixel means for the equalized version of channel 2\cr
#' mean1eq - the mean of the pixel means for the equalized version of channel 1\cr
#' mean2eq - the mean of the pixel means for the equalized version of channel 2\cr
#' cor2 - the mean pixel-pixel correlation between the map pixels in channel 2\cr
#' cor2min - the min pixel-pixel correlation between the map pixels in channel 2\cr
#' cor2max - the max pixel-pixel correlation between the map pixels in channel 2\cr
#' 
#' @param db a database connection
#' @param data the data array for this experiment
#' @param cormat the pixel-pixel corrlations for channel 2 for this data
#' 
#' @return NULL
#-
# TODO: allow specification of which features to compute
# TODO: reformat so that cormat isn't required to be passed in (too memory intensive)
GetDataFeatures <- function(db, data, cormat){
	mids <- dbGetQuery(db, "SELECT id from masks")[,1]
	d <- dbGetQuery(db, "SELECT nx, ny from experiment")
	nrow <- d[1,2]
	ncol <- d[1,1]
	
	#### Setting up the data
	# time x pixels
	data1 <- wrap(data[1,,,], list("time"=c(1), "pixels"=c(2,3)))
	data2 <- wrap(data[2,,,], list("time"=c(1), "pixels"=c(2,3)))
	# Just for computing equalized
	meandata2 <- apply(data[2,,,], c(2,3), mean)
	meandata1 <- apply(data[1,,,], c(2,3), mean)
	
	# equalized, 1xpixels
	meandata2eq <- wrap(array(slidinghistequalC(meandata2,8), c(1,dim(meandata2))), list("time"=c(1), "pixels"=c(2,3)))
	meandata1eq <- wrap(array(slidinghistequalC(meandata1,8), c(1,dim(meandata1))), list("time"=c(1), "pixels"=c(2,3)))
	# pixel means in the write form
	meandata1 <- wrap(array(meandata1, c(1, dim(meandata1))), list("time"=c(1), "pixels"=c(2,3)))
	meandata2 <- wrap(array(meandata2, c(1, dim(meandata2))), list("time"=c(1), "pixels"=c(2,3)))
	gc()
	mask <- rep(FALSE, nrow*ncol)

	#### Setting up the database
	# Add the shape features to the database if they don't already exist
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('var1', 'The variance of the pixel means in channel 1')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('var2', 'The variance of the pixel means in channel 2')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('var1eq', 'The variance of the pixel means in equalized channel 1')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('var2eq', 'The variance of the pixel means in equalized channel 2')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('mean1eq', 'The mean of the pixel means in equalized channel 1')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('mean2eq', 'The mean of the pixel means in equalized channel 2')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('cor2', 'The correlation of the pixel time series in channel 2')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('cor2min', 'The correlation of the pixel time series in channel 2')"))
	try(dbGetQuery(db, "INSERT INTO features (tag, description) values ('cor2max', 'The correlation of the pixel time series in channel 2')"))

	# Get the ids for the shape feature
	var1id <- dbGetQuery(db, "SELECT id from features where tag='var1'")	
	var2id <- dbGetQuery(db, "SELECT id from features where tag='var2'")	
	var1eqid <- dbGetQuery(db, "SELECT id from features where tag='var1eq'")	
	var2eqid <- dbGetQuery(db, "SELECT id from features where tag='var2eq'")	
	mean1eqid <- dbGetQuery(db, "SELECT id from features where tag='mean1eq'")	
	mean2eqid <- dbGetQuery(db, "SELECT id from features where tag='mean2eq'")	
	cor2id <- dbGetQuery(db, "SELECT id from features where tag='cor2'")	
	cor2minid <- dbGetQuery(db, "SELECT id from features where tag='cor2min'")	
	cor2maxid <- dbGetQuery(db, "SELECT id from features where tag='cor2max'")	

	#print("about to compute cormat")
	#cormat = cor(data2)
	#diag(cormat)=NA
	#print("finished cormat")
	print(length(mids))
	for(id in mids[1:length(mids)]){
		if(id%%100==0){
			cat(id, " ")
		}
		mask <- GetMask(db, id)

		# Only compute properties if mask size > 1 pixel and <300 pixels
		if(length(which(mask==1))>1 && length(which(mask==1))<300){

			var1 <- var(meandata1[1,which(mask==1)])
			var2 <- var(meandata2[1,which(mask==1)])
			var1eq <- var(meandata1eq[1,which(mask==1)])
			var2eq <- var(meandata2eq[1,which(mask==1)])
		
			meaneq1 <- mean(meandata1eq[1,which(mask==1)])
			meaneq2 <- mean(meandata2eq[1,which(mask==1)])
		
			smcormat <- cormat[which(mask==1), which(mask==1)]
			diag(smcormat) <- NA
			cor2mean <- mean(smcormat, na.rm=T)
			cor2min <- min(smcormat, na.rm=T)
			cor2max <- max(smcormat, na.rm=T)
						
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", cor2id, ", ", cor2mean,")", sep="")))		
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", cor2minid, ", ", cor2min,")", sep="")))		
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", cor2maxid, ", ", cor2max,")", sep="")))		
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", var1id, ", ", var1,")", sep="")))
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", var2id, ", ", var2,")", sep="")))
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", var1eqid, ", ", var1eq,")", sep="")))
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", var2eqid, ", ", var2eq,")", sep="")))		
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", mean1eqid, ", ", meaneq1,")", sep="")))		
			foo <- try(dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (",id, ", ", mean2eqid, ", ", meaneq2,")", sep="")))		
		}
	}	
}

#-
#' INTERNAL
#' Computes the convex hull of a mask
#' 
#' @details FIXME: there's the issue that maphull(maphull(x))!=maphull(x), but using this anyway
#' 
#' @param mask the mask for which to find the convex hull.  Background pixels should be NA
#' 
#' @return a matrix with 1's on the convex hull of the mask and NA in the background
#-
MaskHull <- function(mask){
	pts <- arrayInd(which(!is.na(mask)), dim(mask))
	if(nrow(pts)==0){
		return(mask)
	}
	# Compute the points on the center of each edge of each pixel
	pts2 <- matrix(0, nrow(pts)*4,2)
	pts2[,1] <- c(pts[,1]-0.5, pts[,1]+0.5, pts[,1], pts[,1])
	pts2[,2] <- c(pts[,2], pts[,2], pts[,2]-0.5, pts[,2]+0.5)
	pts2 <- pts2[!duplicated(pts2),]
	
	# Convex hull of corner points
	ptspoly <- pts2[chull(pts2),]
	np <- nrow(mask)*ncol(mask)
	testpts <- matrix(0, np,2)
	testpts[,1] <- rep(1:nrow(mask), ncol(mask))
	testpts[,2] <- sort(rep(1:ncol(mask), nrow(mask)))
	
	# test whether each point center is within the convex hull
	# uses in.out from package mgcv
	flags <- in.out(ptspoly, testpts)
	
	ret <- matrix(NA, nrow(mask), ncol(mask))
	ret[flags] <- 1
	return(ret)
}

#-
#' INTERNAL
#' Counts the number of pixels not in a mask that are surrounded by at least 3 mask pixels
#' 
#' @details Uses C code from the file countholesC.c
#' 
#' @param mask the mask in which to count holes.  NA or 0 in the background.
#' 
#' @return an integer giving the number of holes in the mask
#' @useDynLib RCI countholesC
#-
CountHolesC <- function(mask){
    d = dim(mask)
	count = 0;
	mask[which(is.na(mask))]=0
	mask[which(mask>0)]=1
    # Use the C function
	out <- .C("countholes",
		mat = as.double(mask),
		as.integer(d),
		ct = as.integer(count)
	)
		
	ret = out$ct
	return(ret)
}


############################################################################
############# Functions to add masks to the database #######################
############################################################################

#TODO: add experiment details to database
#TODO: infrastructure for creating masks in the first place (clean up tree stuff)

#-
#' Add a mask to a database
#' 
#' @details Adds the given mask to the database.  If the mask is already in the database,
#' increments the count for the source of the mask (or adds a new count for a new source)
#' 
#' @param db a database connection object
#' @param mask a matrix giving the mask to add to the database (T/F, 0/1, or NA/1)
#' @param source a string giving the tag for the source of the mask
#' 
#' @return NULL
#-
DbAddMask <- function(db, mask, source){
	# db is the database to connect to
	# mask is the matrix representing the mask (TRUE/FALSE or 0/1 or NA/1)
	# source is the tag for the source the mask is from (string)
	# exp is the tag for the experiment the mask is from
	
	# This adds the mask to the database.
	#  If there is no identical mask for this experiment, adds this mask.  If there is an identical mask, does not add it.
	#  Either way, adds the source if it does not exist (as a feature called source_SOURCETAG), and updates the count in the matching between this source and this mask

	d <- dim(mask)
	if(mode(mask)=="logical"){
		mask[which(mask==FALSE)]=0
		mask[which(mask==TRUE)]=1
	}else{
		mask[which(is.na(mask))]=0
	}

	# Create a sparse mask vector
	mask <- which(as.vector(mask)==1)
	mask <- paste(mask, collapse=" ")
		
	# Get the sourceid, or add the source if it doesn't exist
	qry <- paste("SELECT id FROM features WHERE tag='source_", source,"'", sep="")
	res <- dbGetQuery(db, qry)
	if(length(res[[1]])==0){
		# Add the feature if it's not already in the database
		dbGetQuery(db, paste("INSERT INTO features (tag) values ('source_", source,"')", sep=""))
		res <- dbGetQuery(db, qry)
	}
	sourceid <- res[1,1]
	
	# Try adding the mask with the given source and experiment.  Duplicate masks won't be added
	tryCatch(dbGetQuery(db, paste("INSERT into masks (mask) values ('", mask, "')", sep="")), error=function(er){})
	
	# Find the maskid
	qry <- paste("SELECT id FROM masks WHERE mask = '", mask, "'", sep="")
	res <- dbGetQuery(db, qry)
	maskid <- res[1,1]
	
	# If there is already a link between this source and this mask, increment it by 1
	qry <- paste("SELECT id FROM feats_masks WHERE maskid = ", maskid, " AND featureid = ", sourceid, sep="")
	res <- dbGetQuery(db, qry)[[1]]
	if(length(res)==0){
		# There was no link (new source for this mask)
		dbGetQuery(db, paste("INSERT INTO feats_masks (maskid, featureid, fvalue) values (", maskid, ", ", sourceid, ", 1)", sep=""))
	}else{
		dbGetQuery(db, paste("UPDATE feats_masks SET fvalue=(fvalue + 1) WHERE id=", res[1], sep=""))
	}
}

###########################
## Manipulation of masks ##
###########################

#-
#' INTERNAL
#' Compute the overlap matrix between a set of masks
#' 
#' @details Computes the overlap matrix between a set of masks using C
#' code for efficiency.
#' 
#' @param masklist A list of sparse masks as returned by GetMasks.  Each element
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
	dyn.load("~/Dropbox/thesis/Rcalc/CSource/compmasksC.so")
	nmasks <- length(masklist)

	# Use the C function
	out <- .C("compmasksC",
		as.integer(c(unlist(masklist), 0)),
		as.integer(nmasks),
		ret = as.logical(rep(FALSE, nmasks*nmasks))
	)
   	dyn.unload("~/Dropbox/thesis/Rcalc/CSource/compmasksC.so")
	return(matrix(out$ret, nmasks, nmasks))
}

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
SparseToMatrix <- function(sm, ny=128, nx=128, background=NA){
	mat <- matrix(background, ny, nx)
	mat[sm[which(sm>0)]] <- 1
	return(mat)
}

#-
#' INTERNAL
#' Converts a matrix mask into a sparse mask.  Assumes that the non-mask pixels
#' of the matrix are NA.
#' 
#' @param mat The mask as a matrix with NA in non-mask pixels
#' 
#' @return a vector of indices of the mask pixels
#-
MatrixToSparse <- function(mat){
	return(which(!is.na(as.vector(mat))))
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
#' Selects the masks from the given list that are contained in a region
#' 
#' @details Given a list of masks as returned by GetMasks and a matrix with a mask specifying a
#' region, returns the masks in the masklist that are completely contained in the given region.
#' 
#' @param framemat a matrix of the same size as the masks in masklist with non-NA pixels specifying
#' the region in which to find masks
#' @param masklist a list, as returned by GetMasks, of sparse masks (vectors where the first element
#' is the negative id of the mask and the other elements are the mask indices)
#-
GetInnerMasks <- function(framemat, masklist){
	tmp <- CompMaskC(MatrixToSparse(InvertMask(framemat)), masklist)
	return(masklist[which(!tmp)])
}
