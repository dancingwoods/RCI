#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>

double histequal(double **mat, int nrow, int ncol, int fullmax){
	int cdf, cdfmin, i, j, k, l, count;
	double val, min;
		
	//compute CDF for center point
	val = mat[nrow/2][ncol/2];
	cdf = 0;
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			if(mat[i][j] <= val && mat[i][j]!=-0.12345) cdf++; 
		}
	}
	
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
	
	//compute normalized value for center pixel
	return(rint( 1.0*(cdf - cdfmin)/(nrow*ncol - cdfmin) * (fullmax-1) ));
}

// interfaces with the R functions in features.R
void slidinghistequalC(double *mat, int *dim, int *radius, int *fullmax){

	double **img, **ret, **subar;
	int nrow, ncol, i, j, k, l, fm, rad, px, py;
	
	nrow = dim[0];
	ncol = dim[1];
	fm = *fullmax;
	rad = *radius;
	
 	//allocate image array
	img = (double **)malloc(nrow * sizeof(double *));
	ret = (double **)malloc(nrow * sizeof(double *));
	subar = (double **)malloc((rad*2+1)*sizeof(double *));
	for(i=0; i<nrow; i++){
		img[i] = (double *)malloc(ncol * sizeof(double));
		ret[i] = (double *)malloc(ncol * sizeof(double));
	}
	for(i=0; i<(rad*2+1); i++){
		subar[i] = (double *)malloc((rad*2+1) * sizeof(double));
	}
	
	//assign values to image array, initialize return array to 0
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			img[j][i] = mat[k];
			ret[j][i] = 0;
			k++;
		}
	}
	
	// for each non-edge pixel, calculate normalized value based on
	// the histogram normalization function
	// truncate the window for the edges
//	for(i=rad; i<(nrow-rad); i++){
//		for(j=rad; j<(nrow-rad); j++){
	for(i=0; i<nrow; i++){
		for(j=0; j<nrow; j++){
			
			for(k=0; k<(rad*2+1); k++){
				for(l=0; l<(rad*2+1); l++){
					px = i-rad+k;
					py = j-rad+l;
					if(px>=0 && px<nrow && py>=0 && py<ncol){
						subar[k][l] = img[i-rad+k][j-rad+l];
					}else{
						// Marker that this value shouldn't be considered
						subar[k][l] = -0.12345;
					}
				}
			}
			ret[i][j] = histequal(subar, rad*2+1, rad*2+1, fm);
		}
	}
	
	// place return data in mat
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			mat[k] = ret[j][i];
			k++;
		}
	}

	for(i=0; i<nrow; i++){
		free(img[i]);
		free(ret[i]);
	}
	for(i=0; i<(rad*2+1); i++){
		free(subar[i]);
	}
	free(img);
	free(ret);
	free(subar);
}

