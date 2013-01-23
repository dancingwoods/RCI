#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


// Interfaces with the R functions
// bumps is 0/1 matrix, peaks has non-zero labels at peaks, mat is image for hill climbing
void labelbumpsC(int *bumpmat, int *peakmat, double *mat, int *dim){

	double **img;
	int **bumps, **peaks;
	int nrow, ncol, i, j, k, x, y, curx, cury;
	
	nrow = dim[0];
	ncol = dim[1];
	
	setDouble2D(&img, mat, nrow, ncol);
	setInt2D(&bumps, bumpmat, nrow, ncol);
	setInt2D(&peaks, peakmat, nrow, ncol);
	
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			if(bumps[i][j]>0){
				// Start at each x and y.  hclimb until indices don't change
				curx = x = i;
				cury = y = j;
				hclimb(img, &x, &y, nrow, ncol);
				while(curx!=x || cury!=y){
					curx = x;
					cury = y;
					hclimb(img, &x, &y, nrow, ncol);
				}
				bumps[i][j] = peaks[x][y];
			
			}						
		}
	}
	returnInt2D(bumps, bumpmat, nrow, ncol);
	freeDouble2D(img, nrow);
	freeInt2D(bumps, nrow);
	freeInt2D(peaks, nrow);

}

