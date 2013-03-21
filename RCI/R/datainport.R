# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Functions to read in data                                    
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


#-
#' Convert a folder of text images to a calexp data object 
#'
#' @details This function Converts a directory of text files into a calexp data object in R. 
#' Assumes that the images are individual text files and that they are alphabetically 
#' in order by channel and then by time index.  The directory must contain only 
#' these image files. Each image must have the same dimensions, and there must
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
CreateCalExpFromText <- function(name, imgdir, nchans=2){
	
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


