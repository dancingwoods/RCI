#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>
#include "mylib.h"


// Interfaces with the R functions
void maskcorC(int *maskcount, int *data, int *nt, double *cors){

	int mlen, count, n;
	double meancor, mincor, maxcor, cor, sumx, sumy, sumx2, sumy2, sumxy, x, y;
	double *sums, *sums2, **sumsxy;
	
	count = maskcount[0];
	n = nt[0];
	
	meancor = 0;
	maxcor=-1;
	mincor=1;
	
	sums = (double *)malloc(count * sizeof(double));
	sums2 = (double *)malloc(count * sizeof(double));
	sumsxy = (double **)malloc(count*sizeof(double*));
	for(int i=0; i<count; i++){
		sumsxy[i] = (double *)malloc(count*sizeof(double));
	}
	
	// Compute sum x, x2
	for(int i=0; i<count; i++){
		sums[i]=0;
		sums2[i]=0;
		for(int j=0; j<n; j++){
			sums[i]+=data[i*n+j];
			sums[i]+=data[i*n+j]*data[i*n+j];
		}
	}
	
	// Compute sum xy
	for(int i=0; i<count-1; i++){
		for(int j=i+1; j<count; j++){
			sumsxy[i][j]=0;
			for(int k=0; k<n; k++){
				sumsxy[i][j]+=data[i*n+k]*data[j*n+k];
			}
		}
	}
	
	mincor=1;
	maxcor=-1;
	meancor=0;
	
	for(int i=0; i<count-1; i++){
		for(int j=i+1; j<count; j++){
			cor = (n*sumsxy[i][j] - sums[i]*sums[j])/(sqrt(n*sums2[i] - sums[i]*sums[i])*sqrt(n*sums2[j] - sums[j]*sums[j]));
			Rprintf("%f  ", cor);
			meancor+=cor;
			if(cor<mincor){mincor=cor;}
			if(cor>maxcor){maxcor=cor;}
		}
	}	
	meancor=meancor/count;
	// 
	// // For each pair of pixels, if they are both in the mask, compute correlation
	// // Update meancor (sum), mincor and maxcor
	// for(int i=0; i<(mlen-1); i++){
	// 	if(mask[i]==1){
	// 		sumx = sumx2 = 0;
	// 		for(int k=0; k<nt[0]; k++){
	// 			x = data[i*nt[0]+k];
	// 			sumx+=x;
	// 			sumx2+= x*x;
	// 		}
	// 	}
	// 	for(int j=i+1; j<mlen; j++){
	// 		//Rprintf("%d %d; ", mask[i], mask[j]);
	// 		if(mask[i]==1 && mask[j]==1){
	// 			 sumy = sumy2 = sumxy = 0;
	// 			for(int k=0; k<nt[0]; k++){
	// 				//Rprintf("%d (%d) %d (%d); ", i*nt[0]+k, i, j*nt[0]+k, j);
	// 				x = data[i*nt[0]+k];
	// 				y = data[j*nt[0]+k];
	// 				sumy+=y;
	// 				sumy2+= y*y;
	// 				sumxy+= x*y;
	// 			}
	// 			cor = (nt[0]*sumxy - sumx*sumy)/(sqrt(nt[0]*sumx2 - sumx*sumx)*sqrt(nt[0]*sumy2 - sumy*sumy));
	// 			if(cor<mincor){
	// 				mincor = cor;
	// 			}
	// 			if(cor>maxcor){
	// 				maxcor = cor;
	// 			}
	// 			meancor+=cor;
	// 			count+=1;
	// 		}
	// 	}
	// }
	// meancor = meancor/count;

	free(sums);
	free(sums2);
	for(int i=0; i<count; i++){
		free(sumsxy[i]);
	}
	free(sumsxy);
	cors[0] = meancor;
	cors[1] = mincor;
	cors[2] = maxcor;
	
}

