################################################################
# Bronwyn Woods                                                #
# 2013                                                         #
#                                                              #
# Functions to inport data, and also for basic manipulation    #
# of data that isn't specific to some task.                    #
#                                                              #
################################################################

#-
#' Convert a folder of text images to a calexp data object 
#'
#' @details This function Converts a directory of csv text files into a calexp data object in R. 
#' Assumes that the images are individual csv text files and that they are alphabetically 
#' in order by channel and then by time index.  The directory must contain only 
#' these csv image files. Each image must have the same dimensions, and there must
#' be the same number of images for each channel.
#'
#' @param name a short name to identify this experiment
#' @param imgdir a string giving the directory path for the directory containing the csv images
#' @param nchans the number of channels that exist in the data
#' @return an object of class calexp
#'  \item{name}{the name passed in as an argument to this function}
#'  \item{data}{an array containing the image data, with dimensions nchans-nrows-ncols}
#' @export
#-
CreateCalExpFromCSV <- function(name, imgdir, nchans=2){
	
	file.list <- list.files(path=imgdir)
	nframes <- length(file.list)/nchans

	frame <- read.csv(paste(imgdir, file.list[1], sep=""), header=F, sep="\t")
	frame.dim <- dim(frame)

	data = array(NA, dim=c(nchans, nframes, frame.dim))

	for(j in 1:nchans){
		for(i in 1:nframes){
			data[j,i,,] <- as.matrix(read.csv(paste(imgdir, file.list[(nframes*(j-1))+i], sep=""), header=F, sep="\t"))
		}
	}
	
	ret <- list("name"=name, "data"=data)
	attr(ret, 'class') <- 'calexp'
	return(ret)
}

#-
#' INTERNAL
#' Converts an image matrix to a matrix with coordinates and values in the columns
#' 
#' @param img the matrix to convert
#' 
#' @return A matrix of size npixels-by-3. The first coordinate is the row, the 
#' the column and the third the intensity.
#-
ImageToCoordMat <- function(img){
	# Takes an image and returns a matrix with 3 columns: x coord, y coord, value
	ymax <- nrow(img)
	xmax <- ncol(img)
	col1 <- sort(rep(1:ymax, xmax))
	col2 <- rep(1:xmax, ymax)
	col3 <- as.vector(img)
	return(cbind(col1, col2, col3))
}

#-
#' INTERNAL
#' Filters a vector by frequency using a butterworth filter
#' 
#' @param vec the vector to filter
#' @param low the lower value of the filter
#' @param high the higher value of the filter
#' @param order the order of the butterworth filter
#' @param dt the time (in seconds) of one datapoint.  1/frequency in hz
#' @param type the type of filter, defaults to "BP" bandpass filter.
#' Can also choose other filters offered by the butfilt function
#' 
#' @return the filtered vector
#-
FilterVector <- function(vec, low, high, order=8, dt = 1/1000, type="BP"){
	# Filters a vector using forward-backward butterworth filter
	vec <- butfilt(vec, fl=low, fh=high, npoles=order, type=type, deltat=dt)
	vec <- rev(butfilt(rev(vec), fl=low, fh=high, npoles=order, type=type, deltat=dt))
	return(vec)
}

