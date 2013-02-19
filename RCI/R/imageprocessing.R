
# Requires a kernel of the same size as the image
# Zero padding allows for non-circular convolution
ConvolveImage <- function(image, kernel){
	
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
#' @export
#-
EmbedImage <- function(image, border=NULL, size=NULL){
	if(!is.null(border)){
		ret <- matrix(0, nrow(image)+2*border, ncol(image)+2*border)
		ret[(border+1):(border+nrow(image)), (border+1):(border+nrow(image))] <- image
	}
	if(!is.null(size)){
		ret <- matrix(0, size[1], size[2])
		rb <- floor((size[1]-nrow(image))/2)
		cb <- floor((size[2]-ncol(image))/2)
		ret[(rb+1):(rb+nrow(image)), (cb+1):(cb+ncol(image))] <- image
	}
	return(ret)
}

#-
#' Clips a border from around an image matrix
#' 
#' @param image the image matrix to clip
#' @param border the size of the border to clip.  Must be less than half the image size
#' 
#' @return a matrix with the center (nrow-2*border) by (ncol-2*border) pixels of the image
#' 
#' @export
#-
ClipImage <- function(image, border=NULL, size=NULL){
	if(!is.null(border)){
		nr <- nrow(image)-2*border
		nc <- ncol(image)-2*border
		rb <- border
		cb <- border
	}
	if(!is.null(size)){
		nr <- size[1]
		nc <- size[2]
		rb <- floor((nrow(image)-nr)/2)
		cb <- floor((ncol(image)-nc)/2)
	}
	return(image[(rb+1):(rb+nr), (cb+1):(cb+nc)])
}

#-
#' Embeds an image in a larger matrix of 0's and tapers the image edges using a Hanning window
#' 
#' @param img the image to embed and taper
#' @param taperamt the width of the taper on the edges of the image.  Must be less than or equal to half the image width
#' @param borderamt the width of the border of 0's to add
#' 
#' @export
#-
EmbedAndTaperImage <- function(img, taperamt, borderamt=NULL, size=NULL){
	trow <- 0.5 * (1 - cos((2*pi*1:(taperamt*2))/((taperamt*2)-1)))
	tcol <- 0.5 * (1 - cos((2*pi*1:(taperamt*2))/((taperamt*2)-1)))
	rb <- 0
	cb <- 0
	if(!is.null(borderamt)){
		rb <- borderamt
		cb <- borderamt
	}
	if(!is.null(size)){
		rb <- floor((size[1]-nrow(img))/2)
		cb <- floor((size[2]-ncol(img))/2)
	}
	
	trow <- c(rep(0, rb), trow[1:taperamt], rep(1, ncol(img)-(taperamt*2)), 
				trow[(taperamt+1):(2*taperamt)], rep(0, rb))
	tcol <- c(rep(0, cb), tcol[1:taperamt], rep(1, nrow(img)-(taperamt*2)), 
				tcol[(taperamt+1):(2*taperamt)], rep(0, cb))
	
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
#' @param domain is the image given in the image or Fourier domain?  It will be returned
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
#' @param domain is the image given in the image or Fourier domain?  It will be returned
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
#' Shifts a vector by the specified amount using FFT,
#' but assuming the transform has already been performed.
#' 
#' @param vec the vector to shift
#' @param amt the amount to shift
#' 
#' @return the circularly shifted vector
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


