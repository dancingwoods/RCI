#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


/* Interfaces with R - Compares one mask to a set of masks and returns the overlap vector
 * 
 * newmask is the mask to compare to the others.  Assumed to start with the first element of the mask and end with "0".
 *  assumes a sorted vector of indices for the mask
 * masks is a vecor of integers where each mask (of unknown length) is preceeded by the negative of 
 *   the mask ID and the vector is terminated by a "0".  Assumes that the indices are sorted from smallest
 *   to largest within each mask
 * nmasks is the number of masks to be compared (not counting the new mask)
 * conmat is the vector to put the results.  It is assumed that this vector is length nmasks and that
 *   the masks are in the same order as in the masks input (IDs are not checked against anything)
 */
void compmaskC(int *newmask, int *masks, int *nummasks, int *convec){

	int nmasks, i, cur, s1,e1, flag, snew, enew;
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
	
	// the start of the new mask is the beginning of the vector
	snew = 0;
	i=0;
	// find the length of the new mask
	while(newmask[i]!=0){
		i++;
	}
	enew = i-1;

	// For each mask, check overlap with the new mask
	for(int i=0; i<nmasks; i++){
		// Find the first and last indices of each mask
		s1 = mstarts[i]+1;
		e1 = mstarts[i+1]-1;
		snew=0;

		//Rprintf("Mask %d: %d - %d, %d - %d\n", i, masks[s1], masks[e1], newmask[snew], newmask[enew]);

		// Start with the assumption that they don't overlap
		res = 0;
		// Narrow the ranges until there is no overlap or a match is found
		while(1){
			//Rprintf("  %d - %d, %d - %d\n", masks[s1], masks[e1], newmask[snew], newmask[enew]);
			// if there is no overlap, break and set overlap to 0
			if(masks[s1]>newmask[enew] || newmask[snew]>masks[e1]){
				//Rprintf("no match\n");
				break;
			}

			// if any of the endpoints match, break and set overlap to 1
			if(masks[s1]==newmask[snew] || masks[s1]==newmask[enew] || masks[e1]==newmask[snew] || masks[e1]==newmask[enew]){
				res=1;
				//Rprintf("match\n");
				break;
			}
			
			// if the ranges overlap but there are no matches of the endpoint
			if(masks[s1]<newmask[snew]){
				// if the first element of 1 is smaller, search for the first element of 1
				// that is larger or equal to the first element of 2
				s1 = mybinsearch(masks, s1, e1, newmask[snew]);
			}else{
				snew = mybinsearch(newmask, snew, enew, masks[s1]);					
			}
			
			//Rprintf("  %d - %d, %d - %d\n", masks[s1], masks[e1], newmask[snew], newmask[enew]);
			//break;
		}

		// Set the appropriate elements of the conmat

		convec[i]=res;
	}
	
	free(mstarts);

}

