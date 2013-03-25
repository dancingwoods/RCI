# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Functions to process time series (filtering, spectral estimation)                      
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
	vecf <- fft(vec)[1:ceiling(length(vec)/2)]
	freqlabs <-  seq(0, 0.5/dt, length.out=ceiling(length(vec)/2))
	w <- which(freqlabs>low & freqlabs<high)
	return(Arg(vecf[w[ceiling(length(w)/2)]]))
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

