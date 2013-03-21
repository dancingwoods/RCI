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

GetPhase <- function(vec, low, high, dt=1/1000){
	vec <- FilterVector(vec, low, high, dt=dt)
	vecf <- fft(vec)
	freqlabs <-  seq(0, 0.5/dt, length.out=ceiling(length(vec)/2))
}

#-
#' Uses multi-taper methods to etimate a spectrum for the given vector
#-
MultiTaperSpectrum <- function(vec, dt=0.1247232, dif=T){
	if(dif){
		l = length(vec)
		vec = vec[2:l]-vec[1:(l-1)]
	}
	sp <- spec.mtm(ts(vec[which(!is.na(vec))], deltat=dt), plot=F)
	class(sp) <- c("MTSpectrum")
	return(sp)
}


PeriodSpectrum <- function(vec, dt=0.1247232){
	freq <- list("data"=fft(vec)[1:ceiling(length(vec)/2)])
	freq$freqlabs <-  seq(0, 0.5/dt, length.out=ceiling(length(vec)/2))
	class(freq) <- c("PeriodSpectrum")
	return(freq)
}

