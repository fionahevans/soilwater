/********************************************************************

        memory.c

	      Contains functions for allocating / deallocating memory

        Author:         Fiona Evans
        Date:           July 2002


*********************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>

#include "soilwater.h"


double *alloc_vector(int n)
{
        double *v;

        v = (double *) malloc(n*sizeof(double));
        if (!v) {
                fprintf(stderr, "Memory allocation failure, exitting...\n");
                exit(0);
                }

        return v;
}

double **alloc_matrix(int n, int p)
{
        int i;
        double **m;

        /* allocate pointers to rows */
        m =(double **) malloc(n*sizeof(double*));
        if (!m) {
                fprintf(stderr, "Memory allocation failure, exitting...\n");
                exit(0);
                }

        /* allocate rows and set pointers to them */
        for (i=0; i<n; i++) {
                m[i] = (double *)malloc(p*sizeof(double));
                if (!m[i]) {
                        fprintf(stderr, "Memory allocation failure, exitting...\n");
                        exit(0);
                        }
                }

        /* return pointer to array of pointers to rows */
        return m;
}

double ***alloc_array(int n, int p, int q)
{
	int i, j;
	double ***a;

        a =(double ***) malloc(n*sizeof(double**));
        if (!a) {
                fprintf(stderr, "Memory allocation failure, exitting...\n");
                exit(0);
                }

        for (i=0; i<n; i++) {
                a[i] = (double **)malloc(p*sizeof(double *));
                if (!a[i]) {
                        fprintf(stderr, "Memory allocation failure, exitting...\n");
                        exit(0);
                        }
                }

	for (i=0; i<n; i++) for (j=0; j<p; j++) {
		a[i][j] =  (double *)malloc(q*sizeof(double));
		if (!a[i][j]) {
                        fprintf(stderr, "Memory allocation failure, exitting...\n");
                        exit(0);
                        }
                }

        return a;
}


void free_vector(double *v)
{
        free(v);
}

void free_matrix(double **m, int n)
{
        int i;
        for (i=0; i<n; i++) free (m[i]);
        free(m);
}

void free_array(double ***a, int n, int p)
{
	int i, j;
	for (i=0; i<n; i++)  for (j=0; j<p; j++) free (a[i][j]);
	for (i=0; i<n; i++)  free (a[i]);
	free(a);
}


