#ifndef MATRIX_H
#define MATRIX_H

#include "Commondata.h"

int  bginv(double *a, int m, int n, double *aa, double eps1, double *u, double* v, int ka);
int  bmuav(double *a, int m, int  n, double* u, double* v, double eps1, int ka);
void sss(double *fg, double* cs);
void ppp(double*a, double*e, double* s, double* v, int m, int n);

void mult(double *m1,double *m2,double *result,int i_1,int j_12,int j_2);
int  invers_matrix(double *m1,int n);
void transpose(double *m1,double *m2,int m,int n);

bool MatrixTranspose(int *m1, int row1, int col1, int *m2);
bool MatrixTranspose(float *m1, int row1, int col1, float *m2);
bool MatrixInverse(float *m1, int row1, int col1);

bool MatrixMulti(int *m1, int row1, int col1, float *m2, int row2, int col2, float *m3);
bool MatrixMulti(float *m1, int row1, int col1, int *m2, int row2, int col2, float *m3);
bool MatrixMulti(float *m1, int row1, int col1, float *m2, int row2, int col2, float *m3);


void CalPCATransParam(float *lambda, 
					  float **trans_vector,
					  float **pTrainFeat, 
					  int nTrainSam, int ndim, 
					  int PCA_ndim, float *total_avg);

#endif