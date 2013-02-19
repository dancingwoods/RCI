################################################################
# Bronwyn Woods                                                #
# 2013                                                         #
#                                                              #
# Functions to perform motion correction techniques to calicum #
# imaging data.                                                #
#                                                              #
################################################################

#TODO Fix this to use any of the motion estimation methods, not just phase correlation
#-
#' Removes in-plane motion effects using rigid body alignment of the image frames 
#' 
#' @details Registers the images in a calexp object by rigid body image alignment of
#' the images in a particular channel to the reference image given.  The shifts for each
#' image are estimated by maximum cross-correlation. Sub-pixel shifts are achieved using upsampling
#' by the factor given (providing accuracy of 1/upsamp pixels).  Only one channel is used to 
#' register the data.  This should be the channel with the cleanest spatial information for the
#' best results
#' 
#' @param calexp a calexp object with a \$data field
#' @param refimg a reference image to use for alignment.  Should be the
#' same size as the images in calexp\$data
#' @param channel the channel to use for alignment
#' @param upsamp the upsampling factor (this gives the sup-pixel precision of 1/upsamp)
#' 
#' @return a calexp object with a \$registration field.  The \$data in the returned object has
#' been registered.  The \$registration field records the details of the estimated shifts.
#'  \item{upsamp}{the usampling factor used}
#'  \item{refimg}{the reference image used}
#'  \item{mpars}{the estimated shifts. This is a matrix of size nframes-by-2}
#' 
#' @export
#-
RegisterCalExp <-function(calexp, refimg, channel=1, upsamp=2){
	# Array to hold the registered data
	dmcor <- array(NA, dim=dim(calexp$data))
	
	cat("Calculating image shifts\n")
	
	# Calculate movement parameters
	# Negative 1 is since will be moving img1 instead of img2
	mpars <- -1*t(apply(calexp$data[1,,,], 1, FFTPhaseCor, img2=refimg, upsamp=upsamp))
	
	# Shift data
	nchannels <- dim(calexp$data)[1]
	nframes <- dim(calexp$data)[2]
	
	cat("Shifting images\n")
	
	for(i in 1:nframes){
		for(j in 1:nchannels){
			dmcor[j,i,,] <- TranslateFFT(calexp$data[j,i,,], mpars[i,1], mpars[i,2])
		}
	}
	
	cat("Completed registration\n")
		
	registrationinfo <- list("upsamp"=upsamp, "referenceimage"=refimg, "mpars"=mpars)
	ret = list("name"=calexp$name, "registration"=registrationinfo, "data"=dmcor)
	attr(ret, 'class') <- 'calexp'
	return(ret)
}


#-
#' INTERNAL
#' Computes the rigid body motion alignment parameters by optimizing some error function
#' comparing the two images.
#' 
#' @param img1 the reference image
#' @param img2 the image to align
#' @param taper boolean, should the images be tapered before aligning
#' @param error the error function to use. Can be "mse" for mean squared error or "mae" for mean absolute error
#' @param startval a length 3 vector giving the inital values for the optimization
#' @param pocstart should the POC method be used to initialize the start values
#' 
#' @return a vector of length 3 giving the translation and rotation estimates
#-
OptimShift <- function(img1, img2, taper=TRUE, error="mse", startval=c(0.1,0.1,0), 
					pocstart=TRUE){
	if(pocstart){
		startval[1:2] <- FFTPhaseCor(img1, img2, upsamp=3, subpixel="poc", taper=taper, cortaper=taper)
	}
	
	bigsize <- c(2^(ceiling(log(nrow(img1)+1,2))), 2^(ceiling(log(ncol(img1)+1, 2))))
	if(taper){
		img1 <- EmbedAndTaperImage(img1-mean(img1), min(dim(img1))/4, size=bigsize)
		img2 <- EmbedAndTaperImage(img2-mean(img2), min(dim(img2))/4, size=bigsize)
	}else{
		img1 <- EmbedImage(img1-mean(img1), size=bigsize)
		img2 <- EmbedImage(img2-mean(img2), size=bigsize)
	}
	
	optimf <- function(x, index, fmsfundata){
		
		img2.corrected <- ShiftFFT(fmsfundata$img2, x)
		if(fmsfundata$error=="mse"){
			ret <- (mean((img2.corrected - fmsfundata$img1)^2))
		}else if(fmsfundata$error=="cor"){
			ret <- (-1*cor(as.vector(img2.corrected), as.vector(fmsfundata$img1)))
		}else{
			ret <- (mean(abs(img2.corrected - fmsfundata$img1)))
		}
		
		return(list(f = ret, index = index, this = list(costfargument = fmsfundata)))
	}
	fmsfundata <- list(img1=img1, img2=img2, error=error)
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


	#return(optim(startval, optimf, method=optimtype, control=list("maxit"=2)))
#	offset <- $par[1:3]
#	return(offset)
}

#-
#' Uses optimization of an objective function to compute the best alignment translation
#' between two images
#' 
#' @param img1 the reference image
#' @param img2 the image to align
#' @param taper should the images be tapered before the rotation is computed (hanning window)
#' @param error objective function to be used - "mse" mean squared error, "mae" mean absolute error, "cor" correlation
#' @param optimtype the optimization method to use, as given to optim()
#' 
#' @return a real valued vector of length 2, giving estimates of x and y translation
#-
OptimTranslate <- function(img1, img2, taper=TRUE, error="mse", optimtype="Nelder-Mead"){

	if(taper){
		img1 = EmbedAndTaperImage(img1-mean(img1), min(dim(img1))/4, 20)
		img2 = EmbedAndTaperImage(img2-mean(img2), min(dim(img2))/4, 20)
	}else{
		img1 <- EmbedImage(img1-mean(img1), 20)
		img2 <- EmbedImage(img2-mean(img2), 20)
	}
	
	optimf <- function(pars){
		img2.corrected <- TranslateFFT(img2, pars[1], pars[2])
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

	offset <- optim(c(0,0), optimf, method=optimtype)$par[1:2]
	return(offset)
}

#-
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

#TODO: sometimes max ends up on edge and submat gets messed up - only when registration fails miserably though
#TODO: gaussian filter (or FFT or something) to lowpass filter first - are high frequencies causing problems??
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
#-
FFTPhaseCor <- function(img1, img2, upsamp=2, taper=TRUE, cortaper=TRUE, subpixel="gauss", subrad=3){
	
	# check that img1 and img2 are the same dimensions
	if(any(dim(img1) != dim(img2))){
		cat("Error: Images must have the same dimensions.\n")
		return(NA)
	}
	
	if(taper){
		img1 = EmbedAndTaperImage(img1, min(dim(img1))/4, 20)
		img2 = EmbedAndTaperImage(img2, min(dim(img2))/4, 20)
	}else{
		img1 <- EmbedImage(img1, 20)
		img2 <- EmbedImage(img2, 20)
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
		ccmat <- ReorderFFT(EmbedAndTaperImage(ReorderFFT(ccmat), min(dim(ccmat))/4, 0), inverse=T)
	}

	# Embed in larger matrix (by factor of 'upsample')
	d <- dim(img1)
	ccup <- matrix(0, d[1]*upsamp, d[2]*upsamp)
	dup <- dim(ccup)
	ccup[(floor(upsamp/2)*d[1]+1):(floor(upsamp/2)*d[1]+d[1]), (floor(upsamp/2)*d[2]+1):(floor(upsamp/2)*d[2]+d[2])] <- ReorderFFT(ccmat)

	ccup <- ReorderFFT(ccup, inverse=T)

	# Take inverse transform, keep real part
	ccup <- ReorderFFT(Re(fft(ccup, inverse=T)))
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
				ret  <- (pars[3]/(dup[2]*dup[1])) * (sin(pi*(y+pars[1])) * sin(pi*(x+pars[2]))) / (sin(pi/dup[1]*(y+pars[1])) * sin(pi/dup[2]*(x+pars[2])))
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
			offset <- optim(c(0,0,1), optimf)$par[1:2]
		}else{
			offset <- -1*optim(c(0.1,0.1,1), optimf)$par[1:2]
		}
		trans <- trans+offset
	
	}

	# Return translation coordinates
	return(-1*rev((trans-c(nrow(ccup)/2,ncol(ccup)/2))/upsamp))

}


