#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>


/******* This function taken from the following source ************/
/* Efficient spatial convolution code                                    */
/* (c) Stefan Gustavson, ISY-LiTH, Linkoping, Sweden, 1993               */
/* email: stefang@isy.liu.se                                             */
/* This code may be used and altered freely for non-commercial purposes. */
/* For commercial use, contact me. I'm cheap, but I want some credit.    */

(int*) convolve(int* kernel, int ksize, int* image, int isize, int gain){
/*
   kernel: pointer to the kernel array
   ksize : side of the square kernel
   image : pointer to the image array
   isize : side of the square image
   gain  : normalization factor for the kernel
*/

    int *in,*out,*inbuf,*outbuf;
    int bsize = isize+ksize-1;
    int pad = ksize/2;
    int ix,iy;
    int *offset;
    int *ilower,*iupper,*rlower;
    register int a, *img,*krn,*res,*offs,*kupper; /* heavily used */

    in=image;
    
    /* Pad the input image with zeros */
    inbuf = (int*)malloc(bsize*bsize*sizeof(int));
    for(ix=0;ix<bsize*bsize;ix++)
        inbuf[ix]=0;
    /* The padding could be done with another constant value, or with
       mirrored replicas of the image data to make a circular convolution,
       compatible with FFT filtering methods */

    /* Copy the image data into the input buffer */
    for(iy=0;iy<isize;iy++)
        for(ix=0;ix<isize;ix++)
	    inbuf[(iy+pad)*bsize+ix+pad] = in[iy*isize+ix];

    /* Here you could dispose of the input array */

    /* Build a lookup table to make an efficient convolution loop */
    offset = (int*)malloc(ksize*ksize*sizeof(int));
    for (ix=0;ix<ksize;ix++)
	for(iy=0;iy<ksize;iy++)
	    offset[iy*ksize+ix] = (iy-ksize/2)*bsize+ix-ksize/2;
    
    /* Set up the output buffer and some loop control constants */
    outbuf = (int*)malloc(bsize*bsize*sizeof(int));
    kupper = &kernel[ksize*ksize-1];
    iupper = &inbuf[(bsize-ksize+ksize/2+1)*bsize-ksize+ksize/2+1];
    ilower = &inbuf[ksize/2*bsize+ksize/2];
    rlower = &outbuf[ksize/2*bsize+ksize/2];
    
    /* Perform the actual convolution */
    for(img=ilower, res=rlower; img<iupper; img++)
	{
	    a = 0;
	    for(offs=offset, krn=kernel; krn<=kupper; )
		a += img[*offs++] * (*krn++);  /* I just love this line */
	    *res++ = a / gain;
	}

    /* Here you could dispose of the inbuf and offset arrays */
    
    /* Clip the output image to the original size */
    out = (int*)malloc(isize*isize*sizeof(int));
    for(iy=0;iy<isize;iy++)
        for(ix=0;ix<isize;ix++)
	    out[iy*isize+ix] = outbuf[(iy+pad)*bsize+ix+pad];
    
    /* Here you could dispose of the outbuf array */

    return(out);
}

/* Bronwyn Woods */
/* Function to interface with R */
void ConvolveImageC(double *img, int *kernel, int *imgsize, int *kernesize){

}
