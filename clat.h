#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct _lattice {
    int nrows;
    int ncols;
    long **basis;
    long double **basis_star;
    long double *B;
    long double **mu;
} lattice;


lattice random_lattice(int nrows, int ncols){
    lattice b;
    srand((unsigned int)time(NULL));
    b.nrows = nrows;
    b.ncols = ncols;
    b.basis = (long **)malloc(nrows * sizeof(long *));
    for(int i = 0, j; i < nrows; ++i){
        b.basis[i] = (long *)malloc(ncols * sizeof(long));
        for(j = 0; j < ncols; ++j)
            b.basis[i][j] = rand();
    }
}