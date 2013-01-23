#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


// Interfaces with the R functions
// Takes a matrix with nmask rows and npixel+nsources columns.  Columns 1:npixel are binary indicating if that
// pixel is in the mask.  Columns npixel+1:npixel+nsources are binary indicating the source of the mask.  
// Returns a matrix that is the same size as the input matrix, but merges identical masks and adds the source columns (so that 
// the source column in the result gives the number of times each mask occurs in each source).  removed rows are set to all -1.
void mergemasksC(int *maskmat, int *nmasks, int *nsources, int *npixels){

	int **masks;
	int cmasks, csources, cpixels, flag;
	
	cmasks = nmasks[0];
	csources = nsources[0];
	cpixels = npixels[0];
	
	setInt2D(&masks, maskmat, cmasks, cpixels+csources);
	
	for(int j=1; j<cmasks; j++){
		// For all masks starting with the second
		for(int i=0; i<j; i++){
			// Compare with each previous mask
			flag=1;
			for(int p=0; p<cpixels; p++){
				// For each pixel, compare and break if different
				if(masks[i][p]!=masks[j][p]){
					flag=0;
					break;
				}
			}
			if(flag){
				// If no pixels are different...
				// Add the source pixels to the earlier mask (i)
				for(int s=cpixels; s<(cpixels+csources); s++){
					masks[i][s] = masks[i][s]+masks[j][s];
				}
				// Set everything in row j to -1
				for(int a=0; a<(cpixels+csources); a++){
					masks[j][a]=-1;
				}
			}
		}
	}


	returnInt2D(masks, maskmat, cmasks, cpixels+csources);
	freeInt2D(masks, cmasks);
}

