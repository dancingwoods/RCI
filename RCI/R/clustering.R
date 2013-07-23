# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Functions to perform clustering of cells                           
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
################################################################################

#-
#' Return the corrected Rand index for the two given clusterings
#' 
#' @param clust a vector giving the first clustering
#' @param clust2 a vector giving the second clustering
#' 
#' @return the corrected Rand index value
#' 
#' @export
#-
ClustDistance <- function(clust, clust2){
	return(cluster.stats( as.dist(matrix(runif(length(clust)), length(clust), length(clust))), clust, alt.clustering=clust2)$corrected.rand)
}

#-
#' Return array specifying the correlation matrix for a sliding
#' window of the data
#' 
#' @param seriesmat a matrix with the calcium traces on the columns
#' @param window the size of the sliding window
#' 
#' @return an array with the first two dimensions giving the correlation matrices and the
#' third dimension indicating the start time of the window
#' 
#' @export
#-
CorByTime <- function(seriesmat, window=500){

	WindowCor <- function(start, mat, len){
		return(cor(mat[start:(start+len),]))
	}
	
	nwindow <- dim(seriesmat)[1]-window
	nseries <- dim(seriesmat)[2]
	return(array(sapply(1:nwindow, WindowCor, mat=seriesmat, len=window), c(nseries, nseries, nwindow)))

}


#-
#' Compute circular phase distance between two phases
#' 
#' @param phase1 a phase value in [-pi,pi]
#' @param phase2 a phase value in [-pi,pi]
#' 
#' @return the circular distance between the two phases
#' 
#' @export
#-
PhaseDist <- function(phase1, phase2){
	tmp <- max(phase1, phase2) - min(phase1, phase2)
	if(tmp>pi){
		tmp <- 2*pi-tmp
	}
	return(tmp)
}

#-
#' Get phase distance matrix 
#' 
#' @param mat a matrix with the calcium traces on the columns
#' @param low the lower bound of the frequency to consider
#' @param high the upper bound of the frequency to consider
#' @param dt the sampling rate of the calcium traces
#' 
#' @return a matrix of phase distances between the calcium traces
#' 
#' @export
#-
PhaseDistMat <- function(mat, low=0.789, high=0.791, dt=0.1247232){
	phs <- apply(mat, 2, GetPhase, low=low, high=high, dt=dt)
	ln <- length(phs)
	ret <- matrix(NA, ln, ln)
	for(i in 1:ln){
		for(j in 1:ln){
			ret[i,j] <- PhaseDist(phs[i], phs[j])
		}
	}
	return(ret)
}


#-
#' Cluster segmented ROIs based on correlation or phase distance using k-means
#' 
#' @param calexp the calexp object with the data
#' @param mask a mask identifying the cells to be clustered.  Each unique non-zero/NA value
#' in the mask indicates a cell to be clustered.
#' @param k the number of clusters to find
#' @param criteria the criteria to use for clustering -- 'cor' (correlation)
#' 'phase' (phase of frequency specified in freq)
#' @param freq the frequency band to use to extract the phase for phase-clustering
#' 
#' @export
#-
ClusterCells <- function(calexp, mask, k, criteria="cor", freq=c(0.78,0.81), dt=0.1247232){
	ser <- GetSeries(calexp, mask)
	
	if(criteria=="cor"){
		cmat <- cor(ser)
		clust <- pam(1-cmat, k, diss=T)
	}else if(criteria=="phase"){
		serf <- apply(ser, 2, GetPhase, low=freq[1], high=freq[2], dt=dt)
		PhaseDif=function(i, vec){
			ret <- vec - vec[i]
			w <- which(ret > pi)
			ret[w] <- 2*pi - ret[w]
			return(abs(ret))
		}
		dmat <- sapply(1:length(serf), PhaseDif, vec=serf)
		clust <- pam(dmat, k, diss=T)
	}

	clmask <- matrix(NA, nrow(mask), ncol(mask))
	for(i in 1:nrow(clust$clusinfo)){
		wid <- as.integer(names(which(clust$clustering == i)))
		for(id in wid){
			clmask[which(mask==id)]=i
		}
	}
	return(list("series"=ser, "clusters"=clust, "mask"=clmask))
}

#-
#' Get average time series of the given class in the segmentation stored in a mask database
#' 
#' @param db the database
#' @param calexp the calcium experiment with the data to use to extract the series
#' @param classids a vector of ids specifying which types of ROI to extract traces for
#' @param chan the channel to use for activity traces (defaults to 2)
#' @export
#-
GetAllSeries <- function(db, calexp, classids, chan=2){
	# Note - no check that the data in calexp actually matches the db (size etc)
	
	GetSeries <- function(smaskstring){
		return(apply(apply(arrayInd(as.integer(strsplit(smaskstring, split=", ")[[1]]), 
		       dim(calexp$data)[3:4]), 1, function(vec){return(calexp$data[chan,,vec[1], vec[2]])}), 1, mean))
	}
	
	masks <- dbGetQuery(db$db, paste("select id, mask, segmentation from masks where segmentation in (",paste(classids, collapse=", "),")"))
	series <- sapply(masks[,2], GetSeries)
	colnames(series) <- masks[,1]
	celltype <- masks[,3]
	o <- order(masks[,1])
	series <- series[,o]
	celltype <- celltype[o]
	o <- order(celltype)
	series <- series[,o]
	celltype <- celltype[o]
	return(list("data"=series, "type"=celltype, "maskid"=colnames(series)))
}

#-
#' Return average time series for each cell in a mask
#' 
#' @param calexp the calexp with the data
#' @param mask the mask identifying cells.  Each unique non-zero/NA value in the mask
#' indicates a cell to be clustered.
#' @param channel the channel to get the cell traces from
#' @export
#-
GetSeries <- function(mask, calexp, channel=2){

	# Get unique non-zero, non-NA values
	uvs <- unique(as.vector(mask))
	uvs <- uvs[which(!is.na(uvs))]
	uvs <- uvs[which(!(uvs==0))]
	
	d <- dim(calexp$data)
	data <- matrix(calexp$data[channel,,,], d[2], d[3]*d[4])
	retval <- matrix(NA, d[2], length(uvs))
	for(u in 1:length(uvs)){
		w = which(mask==uvs[u])
		retval[,u] = apply(data[,w], 1, mean)
	}
	colnames(retval)=uvs
	return(retval)
}