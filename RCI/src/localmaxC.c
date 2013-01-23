#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>


int localmax(double **mat, int i, int j, int nrow, int ncol){
	
	double ref;
	// Check each of the 8 surrounding pixels.  If any are greater
	// than pixel of interest, return 0, else return 1
	ref = mat[i][j];
	if(j>0 && mat[i][j-1]>ref){
		return(0);
	}else if(i>0 && j>0 && mat[i-1][j-1]>ref){
		return(0);
	}else if(i>0 && mat[i-1][j]>ref){
		return(0);
	}else if(i>0 && j<(ncol-1) && mat[i-1][j+1] > ref){
		return(0);
	}else if(j<(ncol-1) && mat[i][j+1]>ref){
		return(0);
	}else if(j<(ncol-1) && i<(nrow-1) && mat[i+1][j+1]>ref){
		return(0);
	}else if(i<(nrow-1) && mat[i+1][j]>ref){
		return(0);
	}else if(i<(nrow-1) && j>0 && mat[i+1][j-1]>ref){
		return(0);
	}else{	
		return(1);
	}
}

// Interfaces with the R functions in mode.R
// mat is the data (and returned data), dim is the image dimensions, min determines whether it's
//  local max (0) or local min (1)
void localmaxC(double *mat, int *dim, int *min){

	double **img, **ret;
	int nrow, ncol, i, j, k, minflag;
	
	nrow = dim[0];
	ncol = dim[1];
	
	minflag = *min;
	
	
 	//allocate image array
	img = (double **)malloc(nrow * sizeof(double *));
	ret = (double **)malloc(nrow * sizeof(double *));
	for(i=0; i<nrow; i++){
		img[i] = (double *)malloc(ncol * sizeof(double));
		ret[i] = (double *)malloc(ncol * sizeof(double));
	}
	
	// Assign values to image array, initialize return array to 0
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			img[j][i] = mat[k];
			ret[j][i] = 0;
			k++;
		}
	}
	
	// For each pixel, check all neighbors and place a 1 in the return array if
	//  that pixel is a local max (min)
	
	//if min==1, negate the data first.
	if(minflag==1){
		for(i=0; i<nrow; i++){
			for(j=0; j<ncol; j++){
				img[i][j] = img[i][j]*(-1);
			}
		}
	}

	// check each pixel for being local max
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			ret[i][j] = localmax(img, i, j, nrow, ncol);
		}
	}
	

	// put return values back in mat to pass them back to R
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			mat[k] = ret[j][i];
			k++;
		}
	}
	
	// Free stuff
	for(i=0; i<nrow; i++){
		free(img[i]);
		free(ret[i]);
	}
	free(img);
	free(ret);
}

