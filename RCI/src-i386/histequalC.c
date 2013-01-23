#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


void histequalC(double *matrix, int *nrows, int *ncols, int *maxval, double *result){
	int cdf, cdfmin, i, j, k, l, count, ncol, nrow, fullmax;
	double **mat, **resmat, val, min;
	
	nrow=*nrows;
	ncol=*ncols;
	fullmax = *maxval;
	
	setDouble2D(&mat, matrix, nrow, ncol);
	setDouble2D(&resmat, matrix, nrow, ncol);
	
	//compute mincdf
	min = mat[0][0];
	count = 0;
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			if(mat[i][j]!=-0.12345){
				if(mat[i][j] == min) count++;
				if(mat[i][j] < min){
					min = mat[i][j];
					count = 1;
				} 
			}
		}
	}
	cdfmin = count;
	
	//compute normalized value for each pixel
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			
			//compute CDF for point
			val = mat[i][j];
			cdf = 0;
			for(i=0; i<nrow; i++){
				for(j=0; j<ncol; j++){
					if(mat[i][j] <= val) cdf++; 
				}
			}
			
			resmat[i][j] = rint( 1.0*(cdf - cdfmin)/(nrow*ncol - cdfmin) * (fullmax-1) );
		}
	}
	// put return value in return vector
	returnDouble2D(mat, result, nrow, ncol);
	freeDouble2D(mat, nrow);
	freeDouble2D(resmat, nrow);
}


