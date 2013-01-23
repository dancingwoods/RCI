#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


// Interfaces with the R functions
// mat is the data (and returned data), dim is the image dimensions
void hillclimbC(double *mat, int *dim, int *x, int*y){

	double **img;
	int nrow, ncol, i, j, k, curx, cury;
	
	nrow = dim[0];
	ncol = dim[1];
	
	*x = *x-1;
	*y = *y-1;
	
	setDouble2D(&img, mat, nrow, ncol);
	
	// Start at given x and y.  hclimb until indices don't change
	curx = *x;
	cury = *y;
	hclimb(img, x, y, nrow, ncol);
	while(curx!=*x || cury!=*y){
		curx = *x;
		cury = *y;
		hclimb(img, x, y, nrow, ncol);
	}
	
	*x = *x+1;
	*y = *y+1;
	
	// Free stuff
	freeDouble2D(img, nrow);
}

