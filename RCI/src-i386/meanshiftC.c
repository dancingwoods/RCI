#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <R.h> 
#include <Rmath.h>

// Fraction of bandwidth considered convergence
int MAXITERATE = 1000;

// Helper function to copy vectors
void copyval(double *vfrom, double *vto, int L){
	for(int i=0; i<L; i++){
		vto[i] = vfrom[i];
	}
}


// Compute L2 norm of H(x-xi)
// Note that H should be the bandwidth matrix ^ -1/2
double compnorm(double *x, double *xi, double **H, int F, double *tmp, double *tmp2){
	double ret;
	
	
	for(int i=0; i<F; i++){
		tmp[i] = (x[i] - xi[i]);
	}
	
	for(int i=0; i<F; i++){
		tmp2[i] = 0;
		for(int j=0; j<F; j++){
			tmp2[i] += H[i][j] * tmp[j];
		}
	}
	
	ret = 0;
	for(int i=0; i<F; i++){
		ret+= tmp2[i] * tmp2[i];
	}
	ret = sqrt(ret);
	
	return(ret);
}

/* Computes the mean shift vector
 * INPUT:
 *   double *shiftv - the vector to put the answer
 *   double *center - the current center for the kernel
 *   double **fmat - the feature space data
 *   double *hvec - the vector of bandwidths (H ^ -1/2)
 *   int F - the number of features
 *   int P - the number of pixels
 */
void mshiftvec(double *shiftv, double *center, double **fmat, double **hmat, int F, int P){
	double *topv, bottom, gofnorm, norm, *tmp, *tmp2;
	
	topv = (double *)malloc(F * sizeof(double *));
	
	//tmp is x-xi
	tmp = (double *)malloc(F * sizeof(double *));
	//tmp2 multiplies by H
	tmp2 = (double *)malloc(F * sizeof(double *));
	
		
	bottom=0;
	for(int i=0; i<F; i++){
		topv[i] = 0;
	}
	for(int i=0; i<P; i++){
		//Using normal kernel
		gofnorm = exp(-1 * compnorm(center, fmat[i], hmat, F, tmp, tmp2));
		// //Rprintf("Weight: %f Compnorm: %f  -- center = %f %f -- fmat = %f %f \n", gofnorm, norm, center[0], center[1],
		// //		fmat[i][0], fmat[i][1]);
		for(int j=0; j<F; j++){
			topv[j] += fmat[i][j] * gofnorm;
		}
		bottom += gofnorm;
	}
	
	for(int i=0; i<F; i++){
		shiftv[i] = topv[i]/bottom - center[i];
	}
	
	free(topv);
	free(tmp);
	free(tmp2);
	
}

/* Function for actually performing the meanshift algorithm 
 *   For each pixel, perform mean shift alg to find convergence point
 * 
 */
void meanshift(double **fmat, double **cmat, double **hmat, double *thresh, int F, int P, int maxit){
	double *shiftv, foo;
	int indic;
	int count;
	
	shiftv = (double *)malloc(F * sizeof(double *));
		
	for(int i=0; i<P; i++){
		if(i%100 == 0){
			Rprintf("Pixel number %d out of %d\n", i, P);
		}
		// Initialize the convergence point to the starting point
		copyval(fmat[i], cmat[i], F);
		indic = 1;
		count = 0;
		//Rprintf("Iteration ");
		while(indic){
			count++;
			//Rprintf(" %d ", count);
			mshiftvec(shiftv, cmat[i], fmat, hmat, F, P);
			
			//Rprintf("Center: %f %f %f Shift: %f %f %f\n", cmat[i][0], cmat[i][1], cmat[i][2], shiftv[0],
			//	shiftv[1], shiftv[2]);
			
			indic=0;
			for(int j=0; j<F; j++){
				//foo = fabs(shiftv[j]);
				//Rprintf("Test: %f > %f\n", foo, thresh[j] );
				if(fabs(shiftv[j]) > thresh[j]){
					// If any coordinate changed too much, keep iterating
					indic=1;
				}
				cmat[i][j] += shiftv[j];
			}
			if(count > maxit){
				indic=0;
				for(int j=0; j<F; j++){
					// set all coordinates to arbitrary value to indicate failure to converge
					cmat[i][j] = 0.12345;
				}
			}
		}
		//Rprintf("Converged in %d steps: %f %f %f\n", count, cmat[i][0], cmat[i][1], cmat[i][2]);
	}
	
	free(shiftv);
}


/* Interface with the R function in meanshift.R 
 *
 * INPUT:
 *   fmat - npixels x nfeatures matrix giving the feature values for each pixels (column-wise)
 *   cmat - empty matrix for returning the convergence point for each pixel in feature space (column-wise)
 *   hmatrix - bandwidth matrix, raised to the -1/2
 *   nfeatures - number of features
 *   npixels - number of pixels
 * 
 * The work of the algorithm is done by the meanshift function
 */
void meanshiftC(double *features, double *converged, double *hmatrix, double *thresh, int *nfeatures, 
				int *npixels, int *maxiterations){
	double **fmat, **cmat, **hmat;
	int F, P, M;
	

	F = *nfeatures;
	P = *npixels;
	M = *maxiterations;
	
	//Make the matrices in the right shapes
	//Row-major matrices -- indexed as matrix[pixel][feature]
	fmat = (double **)malloc(P * sizeof(double *));
	for(int i=0; i<P; i++){
		fmat[i] = (double *)malloc(F * sizeof(double));
	}
	cmat = (double **)malloc(P * sizeof(double *));
	for(int i=0; i<P; i++){
		cmat[i] = (double *)malloc(F * sizeof(double));
	}
	
	// hmat is F by F
	hmat = (double **)malloc(F * sizeof(double *));
	for(int i=0; i<F; i++){
		hmat[i] = (double *)malloc(F * sizeof(double));
	}
	
	// Initialize feat to input, 
	for(int i=0; i<P; i++){
		for(int j=0; j<F; j++){
			//Rprintf("%d\n", j*P+i);
			fmat[i][j] = features[j*P + i];
			cmat[i][j] = 0;
		}
	}	
	for(int i=0; i<F; i++){
		for(int j=0; j<F; j++){
			hmat[i][j] = hmatrix[j*F + i];
		}
	}
	
	// Do the actual computation
	meanshift(fmat, cmat, hmat, thresh, F, P, M);
	
	
	// put values from conv back into *cmat to return to R
	for(int i=0; i<P; i++){
		for(int j=0; j<F; j++){
			converged[j*P + i] = cmat[i][j];
		}
	}
	
	//Free the stuff allocated here
	for(int i=0; i<P; i++){
		free(fmat[i]);
		free(cmat[i]);
	}
	for(int i=0; i<F; i++){
		free(hmat[i]);
	}
	free(hmat);
	free(fmat);
	free(cmat);
}




