

#-
#' Clustering
#' 
#' @param calexp the calexp object with the data
#' @param mask a mask identifying the cells to be clustered.  Each unique non-zero/NA value
#' in the mask indicates a cell to be clustered.
#' @param k the number of clusters to find
#' @param criteria the criteria to use for clustering -- 'cor' (correlation)
#' 'phase' (phase of frequency specified in freq)
#' @param freq the frequency band to use to extract the phase for phase-clustering
#-
ClusterCells <- function(calexp, mask, k, criteria="cor", freq=c(0.78,0.81), dt=0.1247232){
	ser <- GetSeries(calexp, mask)
	
	if(criteria=="cor"){
		cmat <- cor(ser)
		clust <- pam(1-cmat, k, diss=T)
	}else if(criteria=="phase"){
		serf <- apply(ser, FilterVector, 1, low=freq[1], high=freq[2], dt=dt)
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
#' See correlation with clusters
#' 
#' @param calexp the calexp object
#' @param clusters the cluster object as returned from ClusterCells
#-
ClusterCorrelation <- function(calexp, clusters){
	clmeans <- matrix()
}


#-
#' Return average time series for each cell in a mask
#' 
#' @param calexp the calexp with the data
#' @param mask the mask identifying cells.  Each unique non-zero/NA value in the mask
#' indicates a cell to be clustered.
#' @param channel the channel to get the cell traces from
#-
GetSeries <- function(calexp, mask, channel=2){

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