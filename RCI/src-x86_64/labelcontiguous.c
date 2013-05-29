#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"

void LabelNeighbors(int **bumpmat, int row, int col, int val, int nrow, int ncol){
	if(row>0 && bumpmat[row-1][col]==1){
		bumpmat[row-1][col]=val;
		LabelNeighbors(bumpmat, row-1, col, val, nrow, ncol);
	}
	if(col>0 && bumpmat[row][col-1]==1){
		bumpmat[row][col-1]=val;
		LabelNeighbors(bumpmat, row, col-1, val, nrow, ncol);
	}
	if(row<(nrow-1) && bumpmat[row+1][col]==1){
		bumpmat[row+1][col]=val;
		LabelNeighbors(bumpmat, row+1, col, val, nrow, ncol);
	}
	if(col<(ncol-1) && bumpmat[row][col+1]==1){
		bumpmat[row][col-1]=val;
		LabelNeighbors(bumpmat, row, col+1, val, nrow, ncol);
	}
}

// Interfaces with the R functions
// bumps is 0/1 matrix
void labelcontiguousC(int *bumpmat, int *dim){

	int **bumps;
	int nrow, ncol, i, j, val;
	
	nrow = dim[0];
	ncol = dim[1];
	
	setInt2D(&bumps, bumpmat, nrow, ncol);
	
	val=2;
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			if(bumps[i][j]==1){
				bumps[i][j] = val;
				LabelNeighbors(bumps, i, j, val, nrow, ncol);
				val++;
			}						
		}
	}
	
	returnInt2D(bumps, bumpmat, nrow, ncol);
	freeInt2D(bumps, nrow);
}

