#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>


/* Binary search.  Returns the index of the match to val.  If no match found,
 *  returns the index of the first value larger than val.  If no value larger than
 *  val, returns the last index
 *
 * v : pointer to the vector in question
 * start : index of first value to consider
 * end : index of the last value to consider
 * val : the value being searched for
 */
int mybinsearch(int *v, int start, int end, int val){
	int cind;
	
	
	// if there is no early termination, the range will eventually be just 2 numbers and 
	// the algorithm will stall
	if(end == (start+1)){
		if(v[start]<val){
			return(end);
		}else{
			return(start);
		}
	}
		
	// find the center index of the given range
	cind = (start+end)/2;

	// compare the center value in v to val
	if(v[cind]<val){
		// if the center value is still smaller, recursively search the top half of the range
		return(mybinsearch(v, cind, end, val));
	}else if(v[cind]>val){
		// note that the cind-1 index must be valid.  otherwise the array is only two ints long and the 
		// check at the beginning of the function will catch it
		if(v[cind-1]<val){
			// check if this is the first value larger than val, and if so return
			return(cind);
		}else{
			// if this is not the first value larger than val, recursively search the bottom half of the range
			return(mybinsearch(v, start, cind, val));
		}
	}else{
		return(cind);
	}
	
	
}


void setDouble2D(double ***imgpointer, double *Rmat, int nrow, int ncol){
	double **img;
	int i, j, k;
	
 	//allocate image array
	img = (double **)malloc(nrow * sizeof(double *));
	for(i=0; i<nrow; i++){
		img[i] = (double *)malloc(ncol * sizeof(double));
	}
	
	// Assign values to image array
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			img[j][i] = Rmat[k];
			k++;
		}
	}	
	*imgpointer = img;
}



void setInt2D(int ***imgpointer, int *Rmat, int nrow, int ncol){
	int **img;
	int i, j, k;
	
 	//allocate image array
	img = (int **)malloc(nrow * sizeof(int *));
	for(i=0; i<nrow; i++){
		img[i] = (int *)malloc(ncol * sizeof(int));
	}
	
	// Assign values to image array
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			img[j][i] = Rmat[k];
			k++;
		}
	}	
	*imgpointer = img;
}

// mallocs a new matrix and copies the values of mat into it
void copyDouble2D(double ***copy, double **mat, int nrow, int ncol){
	double **cp;
	int i,j;
	
	cp = (double **)malloc(nrow * sizeof(double *));
	for(i=0; i<nrow; i++){
		cp[i] = (double *)malloc(ncol * sizeof(double));
	}
	
	for(i=0; i<nrow; i++){
		for(j=0; j<ncol; j++){
			cp[i][j] = mat[i][j];
		}
	}	
	
	*copy = cp;
}


void returnDouble2D(double **retmat, double *retvec, int nrow, int ncol){
	int k, j, i;
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			retvec[k] = retmat[j][i];
			k++;
		}
	}	
}

void returnInt2D(int **retmat, int *retvec, int nrow, int ncol){
	int k, j, i;
	k=0;
	for(i=0; i<ncol; i++){
		for(j=0; j<nrow; j++){
			retvec[k] = retmat[j][i];
			k++;
		}
	}	
}




void freeDouble2D(double **img, int nrow){
	for(int i=0; i<nrow; i++){
		free(img[i]);
	}
	free(img);
}

void freeInt2D(int **img, int nrow){
	for(int i=0; i<nrow; i++){
		free(img[i]);
	}
	free(img);
}


void hclimb(double **mat, int *i, int *j, int nrow, int ncol){
	
	double ref;
	int curx, cury;
	
	// Check each of the 8 surrounding pixels.  set i and j to
	// be the index of the largest value 
	curx = *i;
	cury = *j;

	ref = mat[*i][*j];

	if(*j>0 && mat[*i][*j-1]>ref){
		curx = *i;
		cury = *j-1;
		ref = mat[curx][cury];
	}
	if(*i>0 && *j>0 && mat[*i-1][*j-1]>ref){
		curx = *i-1;
		cury = *j-1;
		ref = mat[curx][cury];
	}
	if(*i>0 && mat[*i-1][*j]>ref){
		curx = *i-1;
		cury = *j;
		ref = mat[curx][cury];
	}
	if(*i>0 && *j<(ncol-1) && mat[*i-1][*j+1] > ref){
		curx = *i-1;
		cury = *j+1;
		ref = mat[curx][cury];
	}
	if(*j<(ncol-1) && mat[*i][*j+1]>ref){
		curx = *i;
		cury = *j+1;
		ref = mat[curx][cury];
	}
	if(*j<(ncol-1) && *i<(nrow-1) && mat[*i+1][*j+1]>ref){
		curx = *i+1;
		cury = *j+1;
		ref = mat[curx][cury];
	}
	if(*i<(nrow-1) && mat[*i+1][*j]>ref){
		curx = *i+1;
		cury = *j;
		ref = mat[curx][cury];
	}
	if(*i<(nrow-1) && *j>0 && mat[*i+1][*j-1]>ref){
		curx = *i+1;
		cury = *j-1;
		ref = mat[curx][cury];
	}
	
	*i = curx;
	*j = cury;
}
