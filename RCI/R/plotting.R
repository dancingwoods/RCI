################################################################
# Bronwyn Woods                                                #
# 2013                                                         #
#                                                              #
# Functions for plotting things.                               #
#                                                              #
################################################################


#-
#' Plots an image of the given matrix with the origin in the upper left
#' 
#' @param img the image matrix to plot
#' @param col a list of colors to use for plotting, defaults to grey
#' @param ... additional graphing parameters
#' 
#' @return NULL
#' 
#' @export
#-
Image <-function(img, col=grey(seq(0,1,0.001)), ...){
	n = nrow(img)
	image(t(img[n:1,]), col=col, axes=F, ...)
}

#-
#' Plots a mask over an already plotted image
#' 
#' @details Given a mask as either a matrix of logicals or a matrix with 1's on the mask,
#' over-plot a semi-transparent colored region on an already plotted image.
#' 
#' @param mask the specification of the mask
#' @param rgb a vector of length 3 giving the color of the mask in RGB (defaults to random)
#' @param alpha the alpha transparency value of the mask (between 0 and 1)
#' @param ... additional graphing parameters
#' 
#' @return NULL
#' 
#' @export
#-
AddMask <- function(mask, rgb=runif(3), alpha=0.5, ...){
	if(mode(mask)=="logical"){
		mask[which(mask==TRUE)]=1
		mask[which(mask==FALSE)]=NA
	}
	mask[which(mask==0)]=NA
	n = nrow(mask)
	r = rgb[1]
	g = rgb[2]
	b = rgb[3]
	image(t(mask[n:1,]), col=rgb(r,g,b,alpha=alpha), add=T, ...)
}