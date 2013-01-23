#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


/* Interfaces with R - Compares all masks in masks to all other masks (computes the complete overlap matrix)
 * 
 * masks is a vecor of integers where each mask (of unknown length) is preceeded by the negative of 
 *   the mask ID and the vector is terminated by a "0".  Assumes that the indices are sorted from smallest
 *   to largest within each mask
 * nmasks is the number of masks
 * conmat is the matrix to put the results.  It is assumed that this matrix is nmasks-by-nmasks and that
 *   the masks are in the same order as in the masks input (IDs are not checked against anything)
 */
void compmasksC(int *masks, int *nummasks, int *conmat){

	int nmasks, i, cur, s1, s2, e1, e2, flag;
	int *mstarts;
	int res;
		
	nmasks = *nummasks;
	mstarts = (int *)malloc((nmasks+1)*sizeof(int));
	
	// Find the indices of the masks starts in masks.  The last
	//  entry in mstarts is actually the index of the 0 signalling the end of the input
	i=0;
	cur=0;
	flag=1;
	while(flag){
		if(masks[i]<=0){
			mstarts[cur]=i;
			cur++;
		}
		if(masks[i]==0){
			flag=0;
		}
		i++;
	}
	
	// For each pair of masks, check overlap
	for(int i=0; i<(nmasks-1); i++){

		for(int j=i+1; j<nmasks; j++){	
			
			// Find the first and last indices of each mask
			s1 = mstarts[i]+1;
			s2 = mstarts[j]+1;
			e1 = mstarts[i+1]-1;
			e2 = mstarts[j+1]-1;
			//Rprintf("i range: %d - %d , j range: %d - %d  (end is %d)\n", s1, e1, s2, e2, mstarts[nmasks]);

			// Start with the assumption that they don't overlap
			res = 0;
			// Narrow the ranges until there is no overlap or a match is found
			while(1){			

				// if there is no overlap, break and set overlap to 0
				if(masks[s1]>masks[e2] || masks[s2]>masks[e1]){

					break;
				}

				// if any of the endpoints match, break and set overlap to 1
				if(masks[s1]==masks[s2] || masks[s1]==masks[e2] || masks[e1]==masks[s2] || masks[e1]==masks[e2]){
					res=1;

					break;
				}
				
				// if the ranges overlap but there are no matches of the endpoint
				if(masks[s1]<masks[s2]){
					// if the first element of 1 is smaller, search for the first element of 1
					// that is larger or equal to the first element of 2
					s1 = mybinsearch(masks, s1, e1, masks[s2]);
				}else{
					s2 = mybinsearch(masks, s2, e2, masks[s1]);					
				}
				
			}

			// Set the appropriate elements of the conmat

			conmat[j*nmasks+i]=res;
			conmat[i*nmasks+j]=res;
		}
	}
	
	free(mstarts);

}

