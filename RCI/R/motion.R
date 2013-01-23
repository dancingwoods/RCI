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
	mpars <- -1*t(apply(calexp$data[1,,,], 1, FFTXCor, img2=refimg, upsamp=upsamp))
	
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
#' Computes sub-pixel shifts values using FFT
#' 
#' @details Computes the sub-pixel shifts by computing the upsampled
#' cross-correlation between the two images and finding the maximum.
#' Computes the upsampled cross-correlation by embedding the product of 
#' FT(img1)* and FFT(img2) in a larger matrix of 0's determined by the 
#' upsampling factor.
#' 
#' @param img1 matrix giving the first image (the reference)
#' @param img2 matrix giving the second image (to be shifted)
#' @param upsamp the factor by which the fft matrix should be expanded
#' @param taper number of pixels to taper the data on the edges of the image
#' 
#' @return a vector of length 2 giving the magnitude of the estimted x and y shift
#' returns NA in the case of improper input
#-
FFTXCor <- function(img1, img2, upsamp=1, taper=0){
	# check that img1 and img2 are the same dimensions
	if(any(dim(img1) != dim(img2))){
		cat("Error: Images must have the same dimensions.\n")
		return(NA)
	}
	# Calculate prod of FFT and conj(FFT)
	cc <- fft(img1)*Conj(fft(img2))
	
	# Taper the crosscorrelation matrix cc
	if(taper>0){
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
	
	ccup <- ReorderFFT(ccup)
	
	# Take inverse transform, keep modulus
	ccup <- Mod(fft(ccup, inverse=T))
		
	# Find maximum
	trans <- arrayInd(which.max(ccup), dim(ccup))
	
	# Translate location of maximum into translation coordinates
	if(trans[1]>(dup[1]/2)){
		ytrans <- (dup[1]-trans[1]+1)/upsamp
	}else{
		ytrans <- -1*(trans[1]-1)/upsamp
	}
	
	if(trans[2]>(dup[1]/2)){
		xtrans <- (dup[2]-trans[2]+1)/upsamp
	}else{
		xtrans <- -1*(trans[2]-1)/upsamp
	}
	
	# Return translation coordinates
	return(c(xtrans, ytrans))
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

