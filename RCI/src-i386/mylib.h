#include <stdlib.h>

void hclimb(double **mat, int *i, int *j, int nrow, int ncol);
void setDouble2D(double ***imgpointer, double *Rmat, int nrow, int ncol);
void returnDouble2D(double **retmat, double *retvec, int nrow, int ncol);
void freeDouble2D(double **img, int nrow);
void freeInt2D(int **img, int nrow);
void copyDouble2D(double ***copy, double **mat, int nrow, int ncol);
void setInt2D(int ***imgpointer, int *Rmat, int nrow, int ncol);
void returnInt2D(int **retmat, int *retvec, int nrow, int ncol);
int mybinsearch(int *v, int start, int end, int val);