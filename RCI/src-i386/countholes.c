#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


// Interfaces with the R functions
// mat is the data (and returned data), dim is the image dimensions
void countholes(double *mat, int *dim, int *count){

	double **img;
	int nrow, ncol, ct, subct;
	
	nrow = dim[0];
	ncol = dim[1];
	setDouble2D(&img, mat, nrow, ncol);

	ct=0;
	for(int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			subct = 0;
			if(img[i][j]==0){
				if(j>0){
					subct+=img[i][j-1];
				}
				//if(i>0 && j>0){
				//	subct+=img[i-1][j-1];
				//}
				if(i>0){
					subct+=img[i-1][j];
				}
				//if(i>0 && j<(ncol-1)){
				//	subct+=img[i-1][j+1];
				//}
				if(j<(ncol-1)){
					subct+=img[i][j+1];
				}
				//if(j<(ncol-1) && i<(nrow-1)){
				//	subct+=img[i+1][j+1];
				//}
				if(i<(nrow-1)){
					subct+=img[i+1][j];
				}
				//if(i<(nrow-1) && j>0){
				//	subct+=img[i+1][j-1];
				//}
				if(subct>=3){ct+=1;}
			}
		}
	}
	*count=ct;
	// Free stuff
	freeDouble2D(img, nrow);
}
