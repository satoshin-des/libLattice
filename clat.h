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
#include <stdbool.h>

typedef struct _lattice {
    int nrows;
    int ncols;
    long **basis;
    long double **basis_star;
    long double *B;
    long double **mu;
} lattice;

// for sub-routine
bool isZero(const int* v, const int n);
double dot_dbl_dbl(long double *x, long double *y, const int n);
double dot_int_dbl(long *x, long double *y, const int n);
long dot_int_int(long *x, long *y, const int n);
lattice _GSO(lattice b);
void SizeReduce(long **b, long double **mu, const int i, const int j, const int m);
lattice _LLL(lattice b, const float d, const int start, const int end);
double _ENUM(long double** mu, long double* B, const double R, int* v, int qq, const int start, const int end);
double _enumerate(long double **mu, long double *B, int *v, const int qq, const int start, const int end);


void print(lattice b, char* flag);
void all_information(lattice b);
lattice Lattice(long** b, const int n, const int m);
lattice random_lattice(int nrows, int ncols);
int *coef2lat(int* v, lattice b);
lattice GSO(lattice b);
long double vol(lattice b);
void Lagrange(lattice b);
lattice LLL(lattice b, const float d);
lattice DeepLLL(lattice b, const float d);
lattice PotLLL(lattice b, const float d);
lattice NanchatteBKZ(lattice b, const int beta, const float d);
void ENUM(long double** mu, long double* B, const int n, const double R, int* v);
void enumerate(long double **mu, long double *B, const int n, int* v);
void Babai(lattice b, long double* w, long* v);

#endif
