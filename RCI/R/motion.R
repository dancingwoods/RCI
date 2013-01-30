################################################################
# Bronwyn Woods                                                #
# 2013                                                         #
#                                                              #
# Functions to perform motion correction techniques to calicum #
# imaging data.                                                #
#                                                              #
################################################################


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
#' @param imagetaper the type of window to use to taper each image
#' @param cortaperwidth number of pixels to taper the data on the edges of the normalized cross-spectrum
#' @param subpixel 'none' for no additional subpixel fitting, 'gauss' for Gaussian fit
#' 
#' @return a vector of length 2 giving the magnitude of the estimted x and y shift
#' returns NA in the case of improper input
#-
FFTPhaseCor <- function(img1, img2, upsamp=1, imagetaper="hanning", cortaperwidth=0, subpixel="none"){
	# check that img1 and img2 are the same dimensions
	if(any(dim(img1) != dim(img2))){
		cat("Error: Images must have the same dimensions.\n")
		return(NA)
	}
	
	TaperImage <- function(img, taper){
		if(taper=="hanning"){
			row = 0.5 * (1 - cos((2*pi*1:ncol(img))/(ncol(img)-1)))
			col = 0.5 * (1 - cos((2*pi*1:nrow(img))/(nrow(img)-1)))
		}else{
			row = rep(1, ncol(img))
			col = rep(1, nrow(img))
		}
		ret = matrix(row, nrow(img), ncol(img), byrow=T)
		ret = ret*col
		return(ret*img)
	}
	img1 = TaperImage(img1, imagetaper)
	img2 = TaperImage(img2, imagetaper)
	
	# Calculate prod of FFT and conj(FFT)
	cc <- fft(img1)*Conj(fft(img2))
	# Normalize to isolate the phase shift
	cc <- cc/Mod(cc)
	
	# Taper the crosscorrelation matrix cc
	if(cortaperwidth>0){
		taper <- cortaperwidth
		tapmat <- matrix(1, nrow(img1), ncol(img1))
		l <- 1:(taper*2)
		lt <- 0.5*(1 - cos(2*pi*l/(taper*2-1)))
		tapmat[1:taper,] <- tapmat[1:taper,]*lt[1:taper]
		tapmat[(nrow(tapmat)-taper+1):nrow(tapmat),] <- tapmat[(nrow(tapmat)-taper+1):nrow(tapmat),]*(matrix(lt[(taper+1):(taper*2)], taper, nrow(tapmat)))
		tapmat[,1:taper] <- tapmat[,1:taper]*t(matrix(lt[1:taper], taper, ncol(tapmat)))
		tapmat[,(ncol(tapmat)-taper+1):ncol(tapmat)] <- tapmat[,(ncol(tapmat)-taper+1):ncol(tapmat)]*t(matrix(lt[(taper+1):(taper*2)], taper, ncol(tapmat)))
		cc <- ReorderFFT(ReorderFFT(cc)*tapmat, T)	
	}
	
	# Embed in larger matrix (by factor of 'upsample')
	d <- dim(img1)
	ccup <- matrix(0, d[1]*upsamp, d[2]*upsamp)
	dup <- dim(ccup)
	ccup[(floor(upsamp/2)*d[1]+1):(floor(upsamp/2)*d[1]+d[1]), (floor(upsamp/2)*d[2]+1):(floor(upsamp/2)*d[2]+d[2])] <- ReorderFFT(cc)
	
	ccup <- ReorderFFT(ccup, inverse=T)
	
	# Take inverse transform, keep modulus
	ccup <- ReorderFFT(Mod(fft(ccup, inverse=T)))
		
	# Find maximum
	trans <- arrayInd(which.max(ccup), dim(ccup))
		
	if(subpixel=="gauss" || subpixel=="poc"){
		
		#create sub-matrix around peak
		#fit this sub-matrix with some function
		#refine estimate of peak with fitted peak

		subrad <- 3*upsamp
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
			#print(optim(c(0.1,0.1,1), optimf))
			offset <- -1*optim(c(0.1,0.1,1), optimf)$par[1:2]
		}
	#	print(trans)
	#	print(offset)
		trans <- trans+offset
		
	}
	
	# Return translation coordinates
	return(-1*rev((trans-c(nrow(ccup)/2,ncol(ccup)/2))/upsamp))
}


#-
#' INTERNAL
#' Reorders the matrix returned by fft 
#' 
#' @details Reorders the matrix returned by the R function fft.  The R function returns
#' the coefficients from low-to-high-to-low frequencies in both dimensions.  The reordering
#' puts the low frequencies in the center of the matrix so that the coefficients go
#' from high-to-low-to-high in each dimension
#' 
#' @param mat a matrix of values to reorder
#' @param inverse if true, takes reordered matrix and returns to order expected by fft.
#' if false, takes matrix from fft and reorders it
#' 
#' @return the reordered matrix
#-
ReorderFFT <- function(mat, inverse=F){
	# Reorders matrix from fft to the logical order
	# if inverse=T, takes intuitive matrix and reorders for fft
	
	if(inverse){
		x <- floor(ncol(mat)/2)+1
		y <- floor(nrow(mat)/2)+1
		return(RotateImg(mat, x, y))		
	}else{
		x <- ceiling(ncol(mat)/2)-1
		y <- ceiling(nrow(mat)/2)-1
		return(RotateImg(mat, x, y))
	}
}

#-
#' INTERNAL
#' Rotates an image by a given number of integer rows and columns
#' 
#' @param mat the matrix to rotate
#' @param x the number of columns to rotate
#' @param y the number of rows to rotate
#' 
#' @return the rotated matrix
#-
RotateImg <- function(mat, x, y){
	nr <- nrow(mat)
	nc <- ncol(mat)
	ret <- matrix(0, nr, nc)
	if(y>0){
		ret[(y+1):nr,] <- mat[1:(nr-y), ]
		ret[1:y,] <- mat[(nr-y+1):nr,]
	}else{
		ret <- mat
	}
	if(x>0){
		mat[,(x+1):nc] <- ret[,1:(nc-x)]
		mat[,1:x] <- ret[,(nc-x+1):nc]
	}else{
		mat <- ret
	}
	return(mat)
}

#-
#' INTERNAL
#' Shifts an image by the given (fractional pixel) amounts
#' 
#' @details Uses the shift theorem to shift the given image by transforming to the
#' Fourier domain.  The shift can be sub-pixel, resulting in Fourier interpolation.
#' 
#' @param img the image (matrix) to shift
#' @param xshift the amount to shift the in x dimension (columns)
#' @param yshift the amount to shift in the y dimension (rows)
#' 
#' @return the shifted image (matrix)
#-
TranslateFFT <- function(img, xshift, yshift){
	# Shifts an image by the specified amounts using FFT

	myshift <- function(vec, amt){
		#shifts a vector by a specified amount using fft
		l <- length(vec)

		# indices of fourier coefficients.  Remeber 0 is on the left, then positive numbers then negative
		#  for example: 0, 1, 2, 3, 4, -5, -4, -3, -2, -1
		myind <- c(0:(ceiling(l/2)-1),-(floor(l/2)):-1)

		# take the fft, multiply by shift
		foo <- fft(vec) * complex(real=cos(2*pi*myind*amt/l), im=sin(2*pi*myind*amt/l))

		# return the inverse fft
		return(fft(foo, inverse=TRUE)/l)
	}

	img_shiftx <- t(apply(img, 1, myshift, amt=xshift))
	img_shiftxy <- apply(img_shiftx, 2, myshift, amt=yshift)
	return(Mod(img_shiftxy))
}

