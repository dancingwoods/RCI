# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: Functions to perform image processing tasks                           
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
#' Finds the extrema in an image.
#' 
#' @param image the image matrix
#' @param maxima boolean, should this function find maxima (default).  If false, finds minima
#' 
#' @return a matrix with 1 at maxima (or minima) and 0 elsewhere
#' 
#' @useDynLib RCI localmaxC
#' @export
#-
GetExtrema <- function(image, maxima=T){	
    d <- dim(image)	
    # Use the C function
	out <- .C("localmaxC",
		mat = as.double(image),
		as.integer(d),
		as.integer(!maxima)
	)			
	ret <- matrix(out$mat, d[1], d[2])
	return(ret)
}


#-
#' Assigns the non-zero pixels of 'region' to one of the maxima of the image by hillclimbing
#' on image
#' 
#' @param region a matrix with 1 in the regions to be assigned and 0 elsewhere
#' @param image the image matrix
#' @param restrict boolean.  should the hill-climbing be restricted to a path entirely within region
#' 
#' @return a matrix with unique integers in the pixels of region corresponding to each local maxima
#' 
#' @useDynLib RCI labelbumpsC
#' @export
#-
AssignToPeaks <- function(region, image, restrict=T){
	peaks <- GetExtrema(image)
	w <- which(peaks==1)
	inds = 1:length(w)
	peaks[w] <- inds
    d <- dim(image)

	if(restrict){
		# Set non-bump pixels of image to 0 to require hill climbing to
		# only occur within bump regions
		image[which(region==0)]=0
	}
	
    # Use the C function
	out <- .C("labelbumpsC",
		ret = as.integer(region),
		as.integer(peaks),
		as.double(image),
		as.integer(d)
	)
		
	ret <- matrix(out$ret, d[1], d[2])
	return(ret)
    
}


#-
#' Convolves an image with the given kernel matrix
#' 
#' @details Uses Fourier methods to convolve the given image with the given kernel
#' 
#' @param image a matrix with the image
#' @param kernel a matrix with the kernel (should be smaller than the image)
#' @param circular boolean, should the convolution be circular (default) or should the
#' image be padded with zeros to prevent circular convolution
#' 
#' @return a matrix of the same size as image with the convolved image
#' @export
#-
ConvolveImage <- function(image, kernel, circular=T){
	PadKernel <- function(kernel, size){
		mat <- matrix(0, size[1], size[2])
		hcenter <- floor(dim(kernel)[2]/2+1)
		hremain <- ncol(kernel)-hcenter
		vcenter <- floor(dim(kernel)[1]/2+1)
		vremain <- nrow(kernel)-vcenter

		mat[1:vcenter, 1:hcenter] <- kernel[(vremain+1):nrow(kernel), (hremain+1):ncol(kernel)]
		mat[(size[1]-vremain+1):size[1], (size[2]-hremain+1):size[2]] <- kernel[1:(vremain), 1:(hremain)]
		mat[(size[1]-vremain+1):size[1], 1:hcenter] <- kernel[1:(vremain), (hremain+1):ncol(kernel)]
		mat[1:vcenter, (size[2]-hremain+1):size[2]] <- kernel[(vremain+1):nrow(kernel), 1:(hremain)]
		return(mat)
	}
	
	dk <- dim(kernel)
	di <- dim(image)
	ipad <- ceiling(max(dk))/2
	if(!circular){
		image <- EmbedImage(image, border=ipad)
	}
	kernel <- PadKernel(kernel, size=dim(image))
	imaget <- fft(image)
	kernelt <- fft(kernel)
	return(ClipImage(Re(fft(imaget*kernelt, inverse=T))/prod(dim(imaget)), size=di))
}

#-
#' Returns a Laplacian of Gaussian kernel
#' 
#' @param kdim the dimension of the (square) kernel to generate
#' @param sigma the standard deviation of the gaussian smoother
#' 
#' @return a matrix giving the LoG kernel
#' 
#' @export
#- 
LoGKernel <- function(kdim, sigma){
	kmat <- matrix(NA, kdim, kdim)
	for(i in 1:kdim){
		for(j in 1:kdim){
			x <- i-ceiling(kdim/2)
			y <- j-ceiling(kdim/2)
			kmat[i,j] <- -1*(1/(pi*sigma^4)) * ( 1 - (x^2+y^2)/(2*sigma^2) ) * exp(-1*(x^2+y^2)/(2*sigma^2))
		}
	}
	kmat <- kmat-(sum(kmat)/kdim^2)
	return(kmat)
}

#-
#' Embeds an image matrix in a larger matrix with a border of 0's
#' 
#' @param image the image to embed
#' @param border the width of the border to add around the edges
#' @param size the resulting size of the image - this must be bigger than the dimensions of image
#' 
#' @return a matrix of size (nrow+2*border) by (ncol+ 2*border)
#' 
#' @details uses the size argument if given, else uses the border argument, else returns the original image
#' 
#' @export
#-
EmbedImage <- function(image, border=NULL, size=NULL){
	if(!is.null(size)){
		ret <- matrix(0, size[1], size[2])
		rb <- floor((size[1]-nrow(image))/2)
		cb <- floor((size[2]-ncol(image))/2)
		ret[(rb+1):(rb+nrow(image)), (cb+1):(cb+ncol(image))] <- image
	}else if(!is.null(border)){
		ret <- matrix(0, nrow(image)+2*border, ncol(image)+2*border)
		ret[(border+1):(border+nrow(image)), (border+1):(border+nrow(image))] <- image
	}else{
		ret <- image
	}
	
	return(ret)
}

#-
#' Clips a border from around an image matrix
#' 
#' @param image the image matrix to clip
#' @param border the size of the border to clip.  Must be less than half the image size
#' @param size the resulting size of the image.
#' 
#' @return a matrix with the center (nrow-2*border) by (ncol-2*border) pixels of the image
#' 
#' @details uses the size argument if given, else uses the border argument, else returns the original image
#' 
#' @export
#-
ClipImage <- function(image, border=NULL, size=NULL){
	if(!is.null(size)){
		nr <- size[1]
		nc <- size[2]
		rb <- floor((nrow(image)-nr)/2)
		cb <- floor((ncol(image)-nc)/2)
	}else if(!is.null(border)){
		nr <- nrow(image)-2*border
		nc <- ncol(image)-2*border
		rb <- border
		cb <- border
	}else{
		return(image)
	}
	return(image[(rb+1):(rb+nr), (cb+1):(cb+nc)])
}

#-
#' Embeds an image in a larger matrix of 0's and tapers the image edges using a Hanning window
#' 
#' @param img the image to embed and taper
#' @param taperamt the width of the taper on the edges of the image.  Must be less than or equal to half the image width
#' @param border the width of the border of 0's to add
#' 
#' @return an image that has been embedded and tapered
#' 
#' @details uses size if given, else uses border, else doesn't embed
#' 
#' @export
#-
EmbedAndTaperImage <- function(img, taperamt, size=NULL, border=NULL){
	trow <- 0.5 * (1 - cos((2*pi*1:(taperamt*2))/((taperamt*2)-1)))
	tcol <- 0.5 * (1 - cos((2*pi*1:(taperamt*2))/((taperamt*2)-1)))
	rb <- 0
	cb <- 0
	if(!is.null(border)){
		rb <- border
		cb <- border
	}
	if(!is.null(size)){
		rb <- floor((size[1]-nrow(img))/2)
		cb <- floor((size[2]-ncol(img))/2)
	}
	if(rb>0){
		rbr <- rep(0,rb)
	}else{
		rbr <- c()
	}
	if(cb>0){
		cbr <- rep(0,cb)
	}else{
		cbr <- c()
	}
	trow <- c(rbr, trow[1:taperamt], rep(1, ncol(img)-(taperamt*2)), 
				trow[(taperamt+1):(2*taperamt)], rbr)
	tcol <- c(cbr, tcol[1:taperamt], rep(1, nrow(img)-(taperamt*2)), 
				tcol[(taperamt+1):(2*taperamt)], cbr)

	# Make larger matrix
	ret <- matrix(trow, length(tcol), length(trow), byrow=T)
	ret <- ret*tcol

	ret[(rb+1):(rb+nrow(img)),(cb+1):(cb+ncol(img))] <- 
		ret[(rb+1):(rb+nrow(img)),(cb+1):(cb+ncol(img))]*img

	return(ret)

}

#-
#' Rotates an image by the given angle using a sequence of Fourier domain shears as
#' described in Eddy 1996.
#' 
#' @param img the image to rotate
#' @param theta the angle to rotate the image
#' @param fdomain is the image given already in the Fourier domain?  It will be returned
#' in the same domain as given (passing in the Fourier domain is helpful to reduce 
#' superfluous transforms if performing additional operations in the Fourier domain).
#' 
#' @export
#-
RotateFFT <- function(img, theta, fdomain=FALSE){	
	# Note that the shearing approach works for rotations less than pi/2.  But rotations
	# that are multiples of pi/2 can be done exactly by just transposing the matrix.  So
	# do that first, then use FFT methods for the remaining rotation.
	
	#Remove full rotations
	if(theta>0){
		theta = theta %% (2*pi)
	}else{
		theta = -1*(-theta %% (2*pi))
	}
	
	# Remove half rotations
	if(theta >= pi){
		theta = theta-pi
		img = img[nrow(img):1,ncol(img):1]
	}
	if(theta <= (-pi)){
		theta = theta+pi
		img = img[nrow(img):1,ncol(img):1]
	}
	
	# Remove quarter rotations
	if(theta >= pi/2){
		theta = theta - pi/2
		img = t(img)[1:nrow(img), ncol(img):1]
	}
	if(theta <= -pi/2){
		theta = theta + pi/2
		img = t(img)[nrow(img):1, 1:ncol(img)]
	}
	
	# Shear to rotate the rest of the way
	
	ShearRow <- function(row, img, shvec){
		return(ShiftFFTVector(img[row,], shvec[row]))
	}

	theta = -1*theta

	# Shear 1 - rows by -tan(theta/2)
	nr <- nrow(img)
	rsec = 1:nr - (nr+1)/2
	shvec <- rsec * (-1*tan(theta/2))
	
	if(!fdomain){
		# If not in Fourier domain, perform the fft on the rows
		img <- t(apply(img, 1, fft))
	}else{
		# If in the Fourier domain, undo transform on columns
		img <- apply(img, 2, fft, inverse=T)/nrow(img)
	}
	img <- t(sapply(1:nr, ShearRow, img=img, shvec=shvec))

	# Regardless of input, inverse transform on rows
	# (note that the two applies cancel themselves out in terms of 
	#  transposing the image)
	img <- apply(img, 1, fft, inverse=T)
	img <- apply(img, 1, Re)/ncol(img)
	
	# Shear 2 - cols by sin(theta)
	img <- apply(img, 2, fft)
	nc <- ncol(img)
	csec <- 1:nc - (nc+1)/2
	shvec <- csec * sin(theta) 
	img <- sapply(1:nc, ShearRow, img=t(img), shvec=shvec)
	
	# inverse transform on cols
	img <- apply(img, 2, fft, inverse=T)
	img <- apply(img, 2, Re)/nrow(img)

	# Shear 3 - rows by -tan(theta/2)	
	img <- t(apply(img, 1, fft))
	shvec <- rsec * -1*tan(theta/2)
	img <- t(sapply(1:nr, ShearRow, img=img, shvec=shvec))


	if(!fdomain){
		# inverse transform on rows
		img <- apply(img, 1, fft, inverse=T)
		img <- apply(img, 1, Re)/ncol(img)
	}else{
		# transform on columns
		return(apply(img, 2, fft))
	}
	
}

#-
#' Shifts an image by the given (fractional pixel) amounts
#' 
#' @details Uses the shift theorem to shift the given image by transforming to the
#' Fourier domain.  The shift can be sub-pixel, resulting in Fourier interpolation.
#' 
#' @param img the image (matrix) to shift
#' @param xshift the amount to shift the in x dimension (columns)
#' @param yshift the amount to shift in the y dimension (rows)
#' @param fdomain is the image given in the Fourier domain?  It will be returned
#' in the same domain as given (passing in the Fourier domain is helpful to reduce 
#' superfluous transforms if performing additional operations in the Fourier domain). 
#' 
#' @return the shifted image (matrix)
#' 
#' @export
#-
TranslateFFT <- function(img, xshift, yshift, fdomain=FALSE){
	if(!fdomain){
		img <- fft(img)
	}
	
	# Shifts an image by the specified amounts using FFT
	img <- t(apply(img, 1, ShiftFFTVector, amt=xshift))
	img <- apply(img, 2, ShiftFFTVector, amt=yshift)
	
	if(!fdomain){
		return(Re(fft(img, inverse=T))/prod(dim(img)))
	}else{
		return(img)
	}
}

#-
#' Shifts an image by the given amount, both translation and rotation
#' 
#' @param img the image to shift
#' @param pars a length-3 vector giving (x-translation, y-translation, rotation angle)
#' @param fdomain is the image given in the Fourier domain?  It will be returned
#' in the same domain as given (passing in the Fourier domain is helpful to reduce 
#' superfluous transforms if performing additional operations in the Fourier domain).
#' @param rotatefirst should rotation be performed before translation
#' 
#' @return the shifted image
#' 
#' @details Uses RotateFFT and TranslateFFT to compute result
#' 
#' @export
#-
ShiftFFT <- function(img, pars, fdomain=FALSE, rotatefirst=FALSE){
	if(!fdomain){
		img <- fft(img)
	}

	if(rotatefirst){
		if(fdomain){
			return(RotateFFT(TranslateFFT(img, pars[1], pars[2], fdomain=TRUE), pars[3], fdomain=TRUE))
		}else{
			return(Re(fft(RotateFFT(TranslateFFT(img, pars[1], pars[2],fdomain=TRUE), pars[3], 
							fdomain=TRUE), inverse=TRUE))/prod(dim(img)))		
		}		
	}else{
		if(fdomain){
			return(TranslateFFT(RotateFFT(img, pars[3], fdomain=TRUE), pars[1], pars[2], fdomain=TRUE))
		}else{
			return(Re(fft(TranslateFFT(RotateFFT(img, pars[3], fdomain=TRUE), pars[1], pars[2], 
							fdomain=TRUE), inverse=TRUE))/prod(dim(img)))		
		}
	}
}

#-
#' INTERNAL
#' Shifts a vector by the specified amount using FFT phase shift,
#' but assuming the Fourier transform has already been performed.
#' 
#' @param vec the vector to shift
#' @param amt the amount to shift
#' 
#' @return the circularly shifted vector
#' 
#' @export
#-
ShiftFFTVector <- function(vec, amt){
	#shifts a vector by a specified amount using fft
	l <- length(vec)

	# indices of fourier coefficients.  Remeber 0 is on the left, then positive numbers then negative
	#  for example: 0, 1, 2, 3, 4, -5, -4, -3, -2, -1
	if(l%%2==1){
		myind <- c(0:(ceiling(l/2)-1),-(floor(l/2)):-1)
	}else{
		myind <- c(0:(ceiling((l-1)/2)-1),0,-(floor((l-1)/2)):-1)
	}

	# multiply by shift
	return(vec * complex(real=cos(2*pi*myind*amt/l), im=sin(2*pi*myind*amt/l)))
}

#-
#' INTERNAL
#' Shifts a vector by the specified amount using FFT
#' 
#' @param vec the vector to shift
#' @param amt the amount to shift
#' 
#' @return the circularly shifted vector
#' 
#' @export
#-
ShiftVector <- function(vec, amt){
	#shifts a vector by a specified amount using fft
	l <- length(vec)

	# indices of fourier coefficients.  Remeber 0 is on the left, then positive numbers then negative
	#  for example: 0, 1, 2, 3, 4, -4, -3, -2, -1  FOR ODD LENGTH
	#               0, 1, 2, 3, 4, 0, -4, -3, -2, -1 FOR EVEN LENGTH
	if(l%%2==1){
		myind <- c(0:(ceiling(l/2)-1),-(floor(l/2)):-1)
	}else{
		myind <- c(0:(ceiling((l-1)/2)-1),0,-(floor((l-1)/2)):-1)
	}

	# take the fft, multiply by shift
	foo <- fft(vec) * complex(real=cos(2*pi*myind*amt/l), im=sin(2*pi*myind*amt/l))

	# return the inverse fft
	return(Re(fft(foo, inverse=TRUE)/l))
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
#' 
#' @export
#-
ReorderFFT <- function(mat, inverse=F){
	# Reorders matrix from fft to the logical order
	# if inverse=T, takes intuitive matrix and reorders for fft

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


