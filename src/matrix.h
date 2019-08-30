#ifndef MATRIX_H
#define MATRIX_H

/*********************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>

#define TINY 1.0e-20
#define E 2.718281828
#define PI 3.141593
#define PI2 6.283185

#ifndef SWAP
#define SWAP(a, b) {swap=(a); (a)=(b), (b)=swap;}
#endif

#ifndef MAX
  #define MAX(X,Y) ((X)>(Y)?(X):(Y))
#endif

#ifndef MIN
  #define MIN(X,Y) ((X)<(Y)?(X):(Y))
#endif

#ifndef FALSE
  #define FALSE 0
#endif

/********************************************************************/

void matmult(double   **A, double   **B, double   **C, int n);
void matrixop(double   **A, double   *x, double   *y, int n, int p);
void matrixcp(double   **A, double   **B, int n, int p);
void veccp(double   *x, double   *y, int p);
void vecadd(double   *x, double   *y, double   *z, int p);
void vecsub(double   *x, double   *y, double   *z, int p);
void vecmult(double   *x, double   *y, int p, double   **M);
void transpose(double   **A, double   **AT, int n, int p);
void printMatrix(double   **M, int n, int p);
void printVector(double *v, int n);

void ludcmp(double **a,int n,int *indx,double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void invert_and_get_det(double **m, double **y, double *val, int N);
void invert(double **m, double **y, int N);


#ifndef ROTATE
  #define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#endif

void jacobi(double **a, int n, double d[], double **v, int *nrot);
void eigsrt(double  *d, double  **v, int n);
void eigen_analysis(double **m, double **ev, double *el, int N);

	
#endif

