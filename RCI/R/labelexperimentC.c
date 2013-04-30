#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"

// NOT REALLY DONE - MIGHT NOT WANT TO

/*
probs is a nrow-by-ncol matrix with columns 
(hand label, probs for each class with negative class first)

ovlp is the overlap matrix, 2 columns with nrowovlp rows.  each number is the relevant row
of the probs matrix

segmentation is a vector for the return value, the column index of the chosen class

enforce indicates whether the hand labels should be enforced
*/
void labelexperimentC(double *probs, int *nrow, int *ncol, int *overlap, int *nrowovlp, 
						int *segmentation, int *enforce){
	double **data, **ovlp;
	int force, nr, nc, nro, flag;
	int *tmpvec;

	nc = *ncol;
	nr = *nrow;
	force = *enforce;
	nro = *nrowovlp;

	// Copy data into data array
	data = (double **)malloc(nc*sizeof(double *));
	for(int i=0; i<nc; i++){
		data[i] = (double *)malloc(nr*sizeof(double));
	}
	for(int i=0; i<nc; i++){
		for(int j=0; j<nr; j++){
			data[i][j] = probs[i*nr + j];
		}
	}	
	// Copy data into overlap array
	ovlp = (integer **)malloc(2*sizeof(int *));
	for(int i=0; i<nro; i++){
		ovlp[i] = (int *)malloc(nro*sizeof(int));
	}
	for(int i=0; i<2; i++){
		for(int j=0; j<nro; j++){
			ovlp[i][j] = overlap[i+nro + j]-1;
		}
	}

	// Initialize tmpvec
	tmpvec = (int *)malloc(nr*sizeof(int));
	for(int i=0; i<nr; i++){
		tmpvec[i]=0;
	}

	for(int curind=0; curind<nrow; curind++){
		// if the seg has not already been set
		if(seg[curind]==0){
			
			// find the set of masks that are connected to the mask under consideration
			flag=1;
			tmpvec[curind]=2;
			// until no more masks are added to set
			while(flag==1){
				flag=0;
				for(i=0; i<nr; i++){
					// if this mask was added on the last pass
					if(tmpvec[i]==2){
					
						// Set tmpvec indices which overlap with the curent mask
						for(int ovind=0; ovind<nro; ovind++){
							if(ovlp[1][ovind] == i && tmpvec[ovlp[2][ovind]]==0){
								tmpvec[ovlp[2][ovind]] = 2;
								flag=1;
							}else if(ovlp[2][ovind] == i && tmpvec[ovlp[1][ovind]]==0){
								tmpvec[ovlp[1][ovind]] = 2;
								flag=1;
							}
						}
						tmpvec[i]=1;
					}
				}				
			}			
			//tmpvec now has 1's in the set of masks to be considered
		
		
			// Find the best set of those marked with 1's in tmpvec
			
			
		}
	}


}

