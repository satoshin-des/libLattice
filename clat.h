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
lattice GSO(lattice b);
lattice LLL(lattice b, const double d);

#endif
