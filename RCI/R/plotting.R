# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Specialized plotting for RCI                                                
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

#-
#' Plots sets of masks over an already plotted image
#' 
#' @details Given a matrix with unique integers for each mask set, overplot each 
#' mask set in a different color (randomly chosen)
#' 
#' @param mask the specification of the mask, unique values for each mask set, and 0 or NA in background
#' @param alpha the alpha transparency value of the mask (between 0 and 1)
#' @param ... additional graphing parameters
#' 
#' @return NULL
#' 
#' @export
#-
AddMaskSet <- function(mask, alpha=0.5, ...){
	mask[which(mask==0)]=NA
	uids <- unique(as.vector(mask))
	uids = uids[which(!is.na(uids))]
	nids <- length(uids)
	cvec = rgb(runif(nids), runif(nids), runif(nids), alpha=alpha)	
	n = nrow(mask)
	for(i in 1:nids){
		submask = mask
		submask[which(submask!=uids[i])]=NA
		image(t(submask[n:1,]), col=cvec[i], add=T, ...)
	}
}

#-
#' Plots a multitaper spectral estimate created by MultiTaperSpectrum
#' 
#' @param spect the multitaper spectrum object
#' @param maglog should the magnitide be plotted on the log scale
#' @param minfreq the minimum frequency to plot
#' @param maxfreq the maximum frequency to plot (NULL for Nyquist frequency)
#' @param ... other graphical parameters
#' 
#' @return NULL
#' 
#' @export
#-
plot.MTSpectrum <- function(spect, maglog=TRUE, minfreq=0, maxfreq=NULL,...){
	mini <- min(which(spect$freq>=minfreq))
	if(is.null(maxfreq)){
		maxi <- length(spect$freq)
	}else{
		maxi <- max(which(spect$freq<=maxfreq))
	}
	
	d = spect$spec[mini:maxi]
	if(maglog){
		d = log(d)
	}
	
	plot(spect$freq[mini:maxi], d, type="l", ...)
}
