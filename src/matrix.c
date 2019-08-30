/********************************************************************

	matrix.c

	Contains functions for doing linear algebra

	Author:         Fiona Evans
    Date:           June 2000
	Last updated    July 2002

	Functions:
		matmult		- matrix multiplication
		matrixop	- multiplies a matrix by a vector
		matrixcp	- copies a matrix
		veccp		- copies a vector
		vecadd		- adds two vectors
		vecsub		- subtracts a vector from another
		vecmult		- multiplies two vectors (elementwise)
		transpose	- transpose of a matrix
		printMatrix	- prints to stderr
		printVector	- prints to stderr

		eigen, eigsrt   - eigenvalues and vectors
		invert_and_get_det
		invert


    Modifications:


*********************************************************************/
#include "include.h"

/* General operations */

/* matmult multiplies sqaure A  by B  to produce C ) */

void matmult(A, B, C, n)
double   **A, **B, **C;
int      n;
{
        int    i, j, k;
	for (i=0;i<n;i++) for (j=0;j<n;j++) C[i][j] = 0;

        for (i=0;i<n;i++) 
		 for (j=0;j<n;j++) 
                        for (k=0;k<n;k++) 
                                C[i][j] += A[i][k]*B[k][j];
                                
}

/* matrixop multiplies A (n*p) by x (n) to produce y (p) */

void matrixop(A, x, y, n, p)
        double   **A, *x, *y;
        int      n, p;

{
        int    i, j;

        for (i=0;i<n;i++) {
                for (j=0;j<p;j++) {
                        y[j] += A[i][j] * x[i];
                        }
                }
}

/* copies matrix A into matrix B (both n by p) */
void matrixcp(A, B, n, p)
	double **A, **B;
	int 	n, p;
{
	int	i, j;
        for (i=0;i<n;i++) {
                for (j=0;j<p;j++) {
                        B[i][j] = A[i][j];
			}
		}
}

/* copies vector x into y (both size p) */
void veccp(x, y, p)
        double *x, *y;
        int p;
{
        int     i;

        for (i=0; i<p; i++) y[i] = x[i];
}


/* adds two vectors of size p */

void vecadd(x, y, z, p)
        double *x, *y, *z;
        int     p;
{
        int     i;

        for (i=0; i<p; i++) z[i] = x[i] + y[i];
}

/* subtracts two vectors of size p */

void vecsub(x, y, z, p)
        double *x, *y, *z;
        int     p;
{
        int     i;

        for (i=0; i<p; i++) z[i] = x[i] - y[i];
}

/* multiplies x by the transpose of y to get the matrix M */

void vecmult(x, y, p, M)
        double *x, *y;
        int     p;
        double **M;
{

        int i, j;

        for (i=0; i<p; i++) {
                for (j=0; j<p; j++) {
                        M[i][j] = x[i] * y[j];
                        }
                }
}

void transpose(A, AT, n, p)
        double **A, **AT;
        int    n, p;
{
        int i, j;


        for (i=0;i<p;i++) {
                for (j=0;j<n;j++) {
                        AT[i][j] = A[j][i];
                        }
                }
}

void printMatrix(M, n, p)
        double **M;
        int     n, p;
{
        int i, j;

        for (i=0; i<n; i++) {
                for (j=0; j<p; j++) {
                        fprintf(stderr, "\t%lf ", M[i][j]);
                        }
                fprintf(stderr, "\n");
                }
}

void printVector(v, n)
        double *v;
        int     n;
{
        int     i;

        for (i=0; i<n; i++) fprintf(stderr, "\t%lf ", v[i]);
        fprintf(stderr, "\n");
}

/************************************************************************/

void ludcmp(double **a,int n,int *indx,double *d)
   /*   Performs LU decomposition of a matrix        */

{
        int i,imax=0,j,k;
        double big,dum,sum,temp;
        double *vv;

	vv = alloc_vector(n);

        vv--;
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0)  {
			/* fprintf(stderr, "Singular matrix in ludcmp\n");
			exit(1); */
			*d = 0;
			for (j=1;j<=n;j++) for (i=1;i<=n;i++) a[i][j]=0;	
			return;
			}
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                              big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        vv++;
        free(vv);
}



void lubksb(double **a, int n, int *indx, double b[])
    /*  Performs back substitution      */

{
        int i,ii=0,ip,j;
        double sum;

        for (i=1;i<=n;i++)
        {
                ip=indx[i];
                sum=b[ip];
                b[ip]=b[i];
                if (ii)
                        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
                else if (sum) ii=i;
                b[i]=sum;
        }
        for (i=n;i>=1;i--)
        {
                sum=b[i];
                for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
                b[i]=sum/a[i][i];
        }
}

/* m is unchanged */

void invert_and_get_det(m, y, val, N)
	double **m;		/* input matrix */
	double **y;		/* inverse */
	double *val;		/* determinant */
	int	N;		/* matrix dimension */

{
    int     i,j,*index;
    double  d,*col;
    double  **cp;

    /* copy input so as not to overwrite */
    cp = alloc_matrix(N, N);
    for (i=0; i<N; i++) for (j=0;j<N; j++) cp[i][j] = m[i][j];

    index = (int *)malloc(N*sizeof(int));
    index--;
    col = alloc_vector(N);
    col--;

    for(i=0;i<N;i++) { /* change offsets so that indices start at 1 */
        cp[i]--;
        y[i]--;
    }
    cp--;
    y--;

    ludcmp(cp,N,index,&d);
    for(j=1;j<=N;j++)
        d = d * cp[j][j];
    *val = d;

    if (d != 0) {	/* if the determinant is zero, then the matrix
			   cannot be inverted */
    	for(j=1;j<=N;j++) {
        	for(i=1;i<=N;i++) col[i] = 0.0;
        	col[j] = 1.0;
        	lubksb(cp,N,index,col);
        	for(i=1;i<=N;i++)
            		y[i][j] = col[i];
    	}
    	col++;
    	index++;
    }
    else {
	fprintf(stderr, "Singular matrix\n");
	for(j=1;j<=N;j++) for(i=1;i<=N;i++) y[i][j] = 0.0;
	}

    /* change offsets back to zero */
    for(i=1;i<=N;i++) {
	cp[i]++;
	y[i]++;
    }
    y++;
    cp++;

    free(index);
    free_vector(col);
    free_matrix(cp, N);
}

void invert(m, y, N)
        double **m;             /* input matrix */
        double **y;             /* inverse */
        int     N;              /* matrix dimension */

{
        double  val;
        invert_and_get_det(m, y, &val, N);
}



/* code taken nearly directly from numerical recipes */
/* unaltered except for change from float to double */

void jacobi(a,n,d,v,nrot)
double **a,d[],**v;
int n,*nrot;
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*nr_vector();
	void nr_error(),nr_free_vector();

	b=(double *) malloc (n*sizeof(double));
	z=(double *) malloc (n*sizeof(double));
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free(z);
			free(b);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr, "Too many iterations in routine JACOBI");
}


void eigsrt(d, v, n)
double  *d;
double  **v;
int     n;
{
        int     i, j, k;
        double  p;

        for (i=1; i<=n; i++) {
                p = d[k=i];
                for (j=i+1; j<=n; j++) if (d[j] >= p) p = d[k=j];
                if (k != i) {
                        d[k] = d[i];
                d[i] = p;
                for (j=1; j<=n; j++) {
                        p = v[j][i];
                        v[j][i] = v[j][k];
                        v[j][k] = p;
                        }
                }
        }
}

/************************************************************************/
void eigen_analysis(m, ev, el, N)
        double **m;             /* input matrix */
        double **ev;            /* eigenvectors */
        double *el;            /* eigenvalues */
        int     N;              /* matrix dimension */

{
    int     i,j;
    int	    nrot;
    double  **cp;


    /* Copy input so as not to over-write */
    cp = (double **) malloc (N * sizeof (double *));
    for (i=0; i<N; i++) cp[i] = (double *) malloc (N *sizeof(double));
    for (i=0; i<N; i++) for (j=0;j<N; j++) cp[i][j] = m[i][j];

    /* Change offsets so that indices start at 1 */
    el--;
    for(i=0;i<N;i++)  {
        cp[i]--;
        ev[i]--;
        }
    cp--;
    ev--;

    jacobi(cp, N, el, ev, &nrot);
    eigsrt(el, ev, N);

	for (i=0; i<N; i++) free(cp[i]);
	free(cp);
}

	


