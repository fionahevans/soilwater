#ifndef MEMORY_H
#define MEMORY_H

double *alloc_vector(int n);
double **alloc_matrix(int n, int p);
double ***alloc_array(int n, int p, int q);
void free_vector(double *v);
void free_matrix(double **m, int n);
void free_array(double ***a, int n, int p);

#endif
