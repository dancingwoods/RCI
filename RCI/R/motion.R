# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Functions to perform motion correction techniques to calicum 
# imaging data.                                                
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



#TODO Fix this to use any of the motion estimation methods
#-
#' Removes in-plane motion effects using rigid body alignment of the image frames 
#' 
#' @details Registers the images in a calexp object by rigid body image alignment of
#' the images in a particular channel to the reference image given.  Initial translation
#' parameters are estimated using Phase-Only correlation.  The parameters are then 
#' optimized using Nelder-Mead optimization of the mean squared error between the images.
#' 
#' @param calexp a calexp object with a \$data field
#' @param refimg a reference image to use for alignment.  Should be the
#' same size as the images in calexp\$data
#' @param channel the channel to use for alignment (typically the structual channel)
#' @param upsamp the upsampling factor (this gives the sup-pixel precision of 1/upsamp)
#' 
#' @return a calexp object with a \$registration field.  The \$data in the returned object has
#' been registered.  The \$registration field records the details of the estimated shifts.
#'  \item{refimg}{the reference image used}
#'  \item{mpars}{the estimated shifts. This is a matrix of size nframes-by-2}
#' 
#' @export
#-
RegisterCalExp <-function(calexp, refimg, channel=1, bigsize=c(256,256)){
	# Array to hold the registered data
	dmcor <- array(NA, dim=dim(calexp$data))
	
	cat("Calculating image shifts\n")
	
	# Calculate movement parameters
	# Negative 1 is since will be moving img1 instead of img2
	mpars <- -1*t(apply(calexp$data[1,,,], 1, OptimShift, img2=refimg, bigsize=bigsize))
	
	# Shift data
	nchannels <- dim(calexp$data)[1]
	nframes <- dim(calexp$data)[2]
	
	cat("Shifting images\n")
	
	for(i in 1:nframes){
		for(j in 1:nchannels){
			dmcor[j,i,,] <- ShiftFFT(calexp$data[j,i,,], mpars[i,])
		}
	}
	
	cat("Completed registration\n")
		
	registrationinfo <- list("referenceimage"=refimg, "mpars"=mpars)
	ret = list("name"=calexp$name, "registration"=registrationinfo, "data"=dmcor)
	attr(ret, 'class') <- 'calexp'
	return(ret)
}


#-
#' INTERNAL
#' Computes the rigid body motion alignment parameters by optimizing some error function
#' comparing the two images. (uses optimization routines in the neldermead package)
#' 
#' @param img1 the reference image
#' @param img2 the image to align
#' @param taper boolean, should the images be tapered before aligning
#' @param bigsize the size of the array in which to embed the tapered images (defaults to next power of 2)
#' @param error the error function to use. Can be "mse" for mean squared error, mae" for mean absolute error, or "cor" for correlation.
#' @param startval a length 3 vector giving the inital values for the optimization (xshift, yshift, theta)
#' @param pocstart should the POC method be used to initialize the start values
#' 
#' @return a vector of length 3 giving the translation and rotation estimates
#' 
#' @export
#-
OptimShift <- function(img1, img2, taper=TRUE, error="mse", startval=c(0.1,0.1,0), 	pocstart=TRUE, bigsize=NULL){
		
	if(pocstart){
		startval[1:2] <- FFTPhaseCor(img1, img2, upsamp=3, subpixel="poc", taper=taper, cortaper=taper)
	}
	
	if(is.null(bigsize)){
		bigsize <- c(2^(ceiling(log(nrow(img1)+1,2))), 2^(ceiling(log(ncol(img1)+1, 2))))
	}
	if(taper){
		img1 <- EmbedAndTaperImage(img1-mean(img1), floor(min(dim(img1))/4), size=bigsize)
		img2 <- EmbedAndTaperImage(img2-mean(img2), floor(min(dim(img2))/4), size=bigsize)
	}else{
		img1 <- EmbedImage(img1-mean(img1), size=bigsize)
		img2 <- EmbedImage(img2-mean(img2), size=bigsize)
	}
	
	optimf <- function(x, index, fmsfundata){
		
		img2.corrected <- Re(fft(ShiftFFT(fmsfundata$fft2, x, fdomain=T), inverse=T))/prod(dim(fmsfundata$fft2))
		if(fmsfundata$error=="mse"){
			ret <- (mean((img2.corrected - fmsfundata$img1)^2))
		}else if(fmsfundata$error=="cor"){
			ret <- (-1*cor(as.vector(img2.corrected), as.vector(fmsfundata$img1)))
		}else{
			ret <- (mean(abs(img2.corrected - fmsfundata$img1)))
		}
		
		return(list(f = ret, index = index, this = list(costfargument = fmsfundata)))
	}
	
	fft2 <- fft(img2)
	fmsfundata <- list(img1=img1,error=error, fft2=fft2)
	attr(fmsfundata, "type") <- "T_FARGS"
	nm <- neldermead.new()
	nm <- neldermead.configure(nm, "-numberofvariables", 3)
	nm <- neldermead.configure(nm, "-function", optimf)
	nm <- neldermead.configure(nm, "-x0", transpose(startval))
	nm <- neldermead.configure(nm, "-costfargument", fmsfundata)
	nm <- neldermead.configure(nm, "-method", 'variable')
	nm <- neldermead.configure(nm, "-maxiter", 50)
	nm <- neldermead.configure(nm, "-tolxmethod", TRUE)
	nm <- neldermead.configure(nm, "-tolxabsolute", 0.001)
	nm <- neldermead.search(nm)
	return(as.vector(neldermead.get(nm, "-xopt")))
}

#-
#' Uses optimization of an objective function to compute the best alignment translation
#' between two images (uses optimization routines in the neldermead package)
#' 
#' @param img1 the reference image
#' @param img2 the image to align
#' @param taper should the images be tapered before the rotation is computed (hanning window)
#' @param error objective function to be used - "mse" mean squared error, "mae" mean absolute error, "cor" correlation
#' @param startval the inital estimate of the shift parameters
#' 
#' @return a real valued vector of length 2, giving estimates of x and y translation
#' 
#' @export
#-
OptimTranslate <- function(img1, img2, taper=TRUE, error="mse", startval=c(0.1,0.1), bigsize=NULL){

	if(is.null(bigsize)){
		bigsize <- c(2^(ceiling(log(nrow(img1)+1,2))), 2^(ceiling(log(ncol(img1)+1, 2))))
	}
	if(taper){
		img1 <- EmbedAndTaperImage(img1-mean(img1), floor(min(dim(img1))/4), size=bigsize)
		img2 <- EmbedAndTaperImage(img2-mean(img2), floor(min(dim(img2))/4), size=bigsize)
	}else{
		img1 <- EmbedImage(img1-mean(img1), size=bigsize)
		img2 <- EmbedImage(img2-mean(img2), size=bigsize)
	}

	optimf <- function(pars){
		img2.corrected <- Re(fft(TranslateFFT(fmsfundata$fft2, x[1], x[2], fdomain=T), inverse=T))/prod(dim(fmsfundata$fft2))
		if(error=="mse"){
			return(sum((img2.corrected - img1)^2))
		}else if(error=="cor"){
			return(-1*cor(as.vector(img2.corrected), as.vector(img1)))
		}else if(error=="mae"){
			return(sum(abs(img2.corrected - img1)))
		}else{
			print("Error: invalid error function specified.")
			return(NA)
		}
	}

	fft2 <- fft(img2)
	fmsfundata <- list(img1=img1,error=error, fft2=fft2)
	attr(fmsfundata, "type") <- "T_FARGS"
	nm <- neldermead.new()
	nm <- neldermead.configure(nm, "-numberofvariables", 2)
	nm <- neldermead.configure(nm, "-function", optimf)
	nm <- neldermead.configure(nm, "-x0", transpose(startval))
	nm <- neldermead.configure(nm, "-costfargument", fmsfundata)
	nm <- neldermead.configure(nm, "-method", 'variable')
	nm <- neldermead.configure(nm, "-maxiter", 50)
	nm <- neldermead.configure(nm, "-tolxmethod", TRUE)
	nm <- neldermead.configure(nm, "-tolxabsolute", 0.001)
	nm <- neldermead.search(nm)
	return(as.vector(neldermead.get(nm, "-xopt")))
}

#-
#' INTERNAL
#' Uses optimization of an objective function to compute the best alignment rotation
#' between two images
#' 
#' @param img1 the reference image
#' @param img2 the image to align
#' @param taper should the images be tapered before the rotation is computed (hanning window)
#' @param error objective function to be used - "mse" mean squared error, "mae" mean absolute error, "cor" correlation
#' @param searchrange the range of rotations to search over in the optimization
#' 
#' @return a real valued estimate of the optimal alignment rotation
#' 
#' @export
#-
OptimRotate <- function(img1, img2, taper=TRUE, error="mse", searchrange=c(-0.1, 0.1)){
	if(taper){
		img1 = EmbedAndTaperImage(img1-mean(img1), min(dim(img1))/4, 20)
		img2 = EmbedAndTaperImage(img2-mean(img2), min(dim(img2))/4, 20)
	}else{
		img1 <- EmbedImage(img1-mean(img1), 20)
		img2 <- EmbedImage(img2-mean(img2), 20)
	}

	optimf <- function(theta){
		img2.corrected <- RotateFFT(img2, theta)
		if(error=="mse"){
			return(sum((img2.corrected - img1)^2))
		}else if(error=="cor"){
			return(-1*cor(as.vector(img2.corrected), as.vector(img1)))
		}else if(error=="mae"){
			return(sum(abs(img2.corrected - img1)))
		}else{
			print("Error: invalid error function specified.")
			return(NA)
		}
	}
	
	offset <- optimize(optimf, searchrange)$minimum
	return(offset)
}


#-
#' INTERNAL 
#' Computes sub-pixel shifts values using phase correlation (FFT implementation)
#' 
#' @details Computes the sub-pixel shifts by computing the upsampled
#' phase correlation between the two images and finding the maximum.  If
#' the parameter gausfit is TRUE, then a gaussian is fit around the peak
#' of the phase correlation function to get additional sub-pixel shift information.
#' This is on top of any upsampling
#' 
#' @param img1 matrix giving the first image (the reference)
#' @param img2 matrix giving the second image (to be shifted)
#' @param upsamp the factor by which the fft matrix should be expanded
#' @param taper boolean, should the images be tapered before alignment
#' @param cortaper boolean, should the normalized cross-spectrum be tapered before being (inverse) transformed
#' @param subpixel 'none' for no additional subpixel fitting, 'gauss' for Gaussian fit, 'poc' for poc function fitting
#' @param subrad the radius of the submatrix used to compute the subpixel fits
#' 
#' @return a vector of length 2 giving the magnitude of the estimted x and y shift
#' returns NA in the case of improper input
#' 
#' @export
#-
FFTPhaseCor <- function(img1, img2, upsamp=2, taper=TRUE, cortaper=TRUE, subpixel="gauss", subrad=3){
	
	# check that img1 and img2 are the same dimensions
	if(any(dim(img1) != dim(img2))){
		cat("Error: Images must have the same dimensions.\n")
		return(NA)
	}
	
	if(taper){
		img1 = EmbedAndTaperImage(img1, floor(min(dim(img1))/4), border=20)
		img2 = EmbedAndTaperImage(img2, floor(min(dim(img2))/4), border=20)
	}else{
		img1 <- EmbedImage(img1, border=20)
		img2 <- EmbedImage(img2, border=20)
	}
	
	# Calculate the ffts of the images
	fft1 <- fft(img1)
	fft2 <- fft(img2)

	# Calculate prod of FFT and conj(FFT)
	ccmat <- fft1*Conj(fft2)
	
	# Normalize to isolate the phase shift
	ccmat <- ccmat/Mod(ccmat)


	# Taper the crosscorrelation matrix cc
	if(cortaper){
		ccmat <- ReorderFFT(EmbedAndTaperImage(ReorderFFT(ccmat), floor(min(dim(ccmat))/4), border=0), inverse=T)
	}

	# Embed in larger matrix (by factor of 'upsample')
	d <- dim(img1)
	ccup <- matrix(0, d[1]*upsamp, d[2]*upsamp)
	dup <- dim(ccup)
	ccup[(floor(upsamp/2)*d[1]+1):(floor(upsamp/2)*d[1]+d[1]), (floor(upsamp/2)*d[2]+1):(floor(upsamp/2)*d[2]+d[2])] <- ReorderFFT(ccmat)

	ccup <- ReorderFFT(ccup, inverse=T)

	# Take inverse transform, keep real part
	ccup <- ReorderFFT(Re(fft(ccup, inverse=T))/prod(dim(ccup)))
	# Find maximum
	trans <- arrayInd(which.max(ccup), dim(ccup))

	if(subpixel=="gauss" || subpixel=="poc"){
	
		#create sub-matrix around peak
		#fit this sub-matrix with some function
		#refine estimate of peak with fitted peak

		submat <- ccup[(trans[1]-subrad):(trans[1]+subrad), (trans[2]-subrad):(trans[2]+subrad)]
		submat <- submat/sum(submat)
		xinds = (-subrad):(subrad)
		yinds = (-subrad):(subrad)
			
		gmat <- function(pars){
			#pars = c(y,x,var)
			if(pars[3]<0){return(matrix(0, 2*subrad+1, 2*subrad+1))}
			ret = matrix(1, 2*subrad+1, 2*subrad+1)
			ret = ret*dnorm(yinds,pars[1], pars[3])
			ret = t(t(ret)*dnorm(yinds,pars[2], pars[3]))
			return(ret)
		}
	
		# Using the function from Takita2003
		pocmat <- function(pars){
			#pars = c(y,x)
		
			#hack right now - if either par 0, return matrix of zeros
			if(pars[1]==0 || pars[2]==0){
				return(matrix(0, 2*subrad+1, 2*subrad+1))
			}
			loc <- function(y,x,pars){
				ret  <- (pars[3]/(dup[2]*dup[1])) * (sin(pi*(y+pars[1])) * 
						sin(pi*(x+pars[2]))) / (sin(pi/dup[1]*(y+pars[1])) * 
						sin(pi/dup[2]*(x+pars[2])))
			}
			ret = matrix(NA, 2*subrad+1, 2*subrad+1)
			for(x in xinds){
				for(y in yinds){
					ret[y+subrad+1,x+subrad+1]=loc(y,x,pars)
				}
			}
			return(ret)
		}
	
		optimf <- function(pars){
			if(subpixel=="gauss"){
				return(sum((gmat(pars)-submat)^2))
			}else{
				return(sum((pocmat(pars)-submat)^2))
			}
		}
	
		if(subpixel=="gauss"){
			offset <- optim(c(0,0,1), optimf)$par
		}else{
			offset <- -1*optim(c(0.1,0.1,1), optimf)$par
		}
		trans <- trans+offset[1:2]
	}

	# Return translation coordinates
	return(-1*rev((trans-c(nrow(ccup)/2,ncol(ccup)/2))/upsamp))

}


#-
#' Performs intensity correction on the given calcium experiment
#' 
#' @param calexp the data to be corrected is in the $data element of this calexp object
#' @param cortype the type of correction to perform.  `ar` for autoregressive filter
#' @param order the order of the model to fit (for ar type)
#' @param naclip should NAs produced at the beginning of the experiment be clipped off (by AR model, for instance)
#' 
#' @export
#-
IntensityCorrection <- function(calexp, cortype="ar", order=25, naclip=T){
	ARFilter <- function(vec, order){
		# Takes the first difference of the given vector, fits an AR
		# model of the given order, and returns the residuals of the 
		# model plus the vector mean
		vecmean <- mean(vec)
		vec <- vec-vecmean
		dv <- vec[-1]-vec[1:(length(vec)-1)]
		return(c(NA, ar(dv, F, order)$resid+vecmean))
	}
	
	RegMPars <- function(vec, mpars){
		#mpars should be a matrix of size length(vec):mpars
		# if mpars is longer than vec, the first bit of mpars will be ignored
		if(nrow(mpars)>length(vec)){
			d = nrow(mpars)-length(vec)
			mpars = mpars[-(1:d),]
		}
		mod = lm(vec ~ ., data = as.data.frame(mpars))
		return(mod$residuals+mean(vec, na.rm=T))
	}
	
	
	chs <- dim(calexp$data)[1]
	for(i in 1:chs){
		if(cortype=="ar"){
			calexp$data[i,,,] <- apply(calexp$data[i,,,], c(2,3), ARFilter, order=order)
		}else if(cortype=="regmpars"){
			calexp$data[i,,,] <- apply(calexp$data[i,,,], c(2,3), RegMPars, mpars=calexp$registration$mpars)			
		}
	}
	
	# Clip NA at beginning from AR filtering
	if(naclip){
		w = which(is.na(calexp$data[1,,1,1]))
		if(length(w)>0){
			calexp$data = calexp$data[,-(1:max(w)),,]
			calexp$naclip = max(w)
		}
	}
	
	calexp$intensitycor <- c(calexp$intensitycor, cortype)
	return(calexp)
}















