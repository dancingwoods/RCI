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

PlotPerformance <- function(perf, ...){
	
	if(nrow(perf)==2){
		plotmat <- matrix(NA, 3,4)
		plotmat[,1] = c(perf[1,2]-perf[1,3]-perf[1,5], perf[1,5], perf[1,3])
		plotmat[,2] = c(perf[1,4], 0, 0)
		plotmat[,3] = c(perf[2,2]-perf[2,3]-perf[2,5], perf[2,5], perf[2,3])
		plotmat[,4] = c(perf[2,4], 0, 0)
		print(plotmat)
		barplot(plotmat[,c(1,3)], col=c("black", grey(0.75), grey(1)), 
			names.arg=c("Neurons", "Astrocytes"),  cex.names=1.5, ...)
		barplot(plotmat[,c(2,4)], col=rgb(0.5, 0.7, 0.5), add=T, width=0.5, space=c(1.5, 1.4))
		legend("topright", cex=1.5, bty="n", fill=c("black", grey(0.75), "white", rgb(0.5, 0.7, 0.5)), 
			legend=c("Correct", "Marginal", "Missed", "New"))
	
		text(0.5, cex=1.5,plotmat[1,1]/2, paste(toString(round(plotmat[1,1]/sum(plotmat[,1])*100)), "%"), col="white")
		text(0.5, cex=1.5,plotmat[2,1]/2 +plotmat[1,1], paste(toString(round(plotmat[2,1]/sum(plotmat[,1])*100)), "%"))
		text(0.5, cex=1.5,plotmat[3,1]/2 + sum(plotmat[1:2,1]), paste(toString(round(plotmat[3,1]/sum(plotmat[,1])*100)), "%"))

		text(1.7, cex=1.5,plotmat[1,3]/2, paste(toString(round(plotmat[1,3]/sum(plotmat[,3])*100)), "%"), col="white")
		text(1.7, cex=1.5,plotmat[2,3]/2 +plotmat[1,3], paste(toString(round(plotmat[2,3]/sum(plotmat[,3])*100)), "%"))
		text(1.7, cex=1.5,plotmat[3,3]/2 + sum(plotmat[1:2,3]), paste(toString(round(plotmat[3,3]/sum(plotmat[,3])*100)), "%"))
	}
	if(nrow(perf)==1){
		plotmat <- matrix(NA, 3,2)
		plotmat[,1] = c(perf[1,2]-sum(perf[1,c(3,5)]), perf[1,5], perf[1,3])
		plotmat[,2] = c(perf[1,4], 0, 0)

		totalcells <- sum(plotmat[,1])
		maxy = min(max(totalcells, perf[1,4]), totalcells*1.5)

		barplot(as.matrix(plotmat[,1], 3, 1), col=c("black", grey(0.75), grey(1)), 
			names.arg=c("Cells"),  ylim=c(0,maxy), cex.names=1.5, xlim=c(0,1.5), ...)
		barplot(plotmat[,2], col=rgb(0.5, 0.7, 0.5, alpha=0.75), add=T, width=0.5, space=c(1.7))
	
		text(0.5, cex=1.5,plotmat[1,1]/2, paste(toString(round(plotmat[1,1]/sum(plotmat[,1])*100)), "%"), col="white")
		text(0.5, cex=1.5,plotmat[2,1]/2 +plotmat[1,1], paste(toString(round(plotmat[2,1]/sum(plotmat[,1])*100)), "%"))
		text(0.5, cex=1.5,plotmat[3,1]/2 + sum(plotmat[1:2,1]), paste(toString(round(plotmat[3,1]/sum(plotmat[,1])*100)), "%"))
		
		text(1.1, cex=1.5, maxy-50, toString(plotmat[1,2]))
		
	}
}

ImageDb <- function(db, imagetag="mimg2"){
	Image(GetImage(db, imagetag))
}

PlotMaskSetByID <- function(db, ids, rgb=NULL){
	dims <- as.integer(dbGetQuery(db$db, "select nx, ny from experiment"))
	mask <- matrix(NA, dims[2], dims[1])
	for(id in ids){
		mask[GetMask(db,id)]=id
	}
	if(is.null(rgb)){
		PlotMaskSet(mask)
	}else{
		PlotMask(mask, rgb=rgb)
	}
}

PlotSegmentation <- function(db, classid, rgb=NULL){
	maskids <- dbGetQuery(db$db, paste("select id from masks where segmentation=",classid))[,1]
	PlotMaskSetByID(db, maskids, rgb)
}

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
PlotMask <- function(mask, rgb=runif(3), alpha=0.5, ...){
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
PlotMaskSet <- function(mask, alpha=0.5, rgb=NULL,...){
	mask[which(mask==0)]=NA
	uids <- unique(as.vector(mask))
	uids <- uids[which(!is.na(uids))]
	nids <- length(uids)
	
	if(is.null(rgb)){	
		cvec <- rgb(runif(nids), runif(nids), runif(nids), alpha=alpha)	
	}else{
		cvec <- rgb(rep(rgb[1], nids), rep(rgb[2], nids), rep(rgb[3], nids), alpha=alpha) 
	}
	
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
