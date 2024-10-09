/**
 * @file clat.h
 * @author Arata Sato (23lc002y@rikkyo.ac.jp)
 * @brief 
 * @version 0.1
 * @date 2024-10-08
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef CLAT_H
#define CLAT_H

#include <stdio.h>

typedef struct _lattice {
    int nrows;
    int ncols;
    long **basis;
    long double **basis_star;
    long double *B;
    long double **mu;
} lattice;

void print(lattice b, char* flag);
lattice random_lattice(int nrows, int ncols);
int *coef2lat(int* v, lattice b);
lattice GSO(lattice b);
long double vol(lattice b);
lattice LLL(lattice b, const float d);
lattice DeepLLL(lattice b, const float d);
lattice PotLLL(lattice b, const float d);
lattice BKZ(lattice b, const int beta, const float d, const int lp);
int *ENUM(long double** mu, long double* B, const int n, const double R);
int *enumerate(long double **mu, long double *B, const int n);

#endif
