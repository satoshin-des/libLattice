#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "clat.h"


/**
 * @brief Computes inner product value of vectors ``x`` and ``y``.
 * 
 * @param x vector
 * @param y vector
 * @param n size of vector
 * @return double 
 */
double dot_dbl_dbl(long double *x, long double *y, const int n){
    double s = 0.0;
    for(int i = 0; i < n; ++i) s += x[i] * y[i];
    return s;
}
/**
 * @brief Computes inner product value of vectors ``x`` and ``y``.
 * 
 * @param x vector
 * @param y vector
 * @param n size of vector
 * @return double 
 */
double dot_int_dbl(long *x, long double *y, const int n){
    double s = 0.0;
    for(int i = 0; i < n; ++i) s += y[i] * x[i];
    return s;
}

/// @brief The function ```print``` prints basis 
/// @param b 
/// @param flag 
void print(lattice b, char* flag){
    if(b.nrows == 0 || b.ncols == 0){
        puts("Lattice initialization error: nrows or ncols is undefined.");
        exit(0);
    }
    if(strcmp(flag, "basis") == 0){
        for(int i = 0, j; i < b.nrows; ++i){
            for(j = 0; j < b.ncols; ++j) printf("%ld ", b.basis[i][j]);
            puts("");
        }
    }else if(strcmp(flag, "GSO") == 0){
        puts("GSO-vectors' information:");
        for(int i = 0, j; i < b.nrows; ++i){
            for(j = 0; j < b.ncols; ++j) printf("%Lf ", b.basis_star[i][j]);
            puts("");
        }
        puts("GSO-coefficient's information:");
        for(int i = 0, j; i < b.nrows; ++i){
            for(j = 0; j < b.nrows; ++j) printf("%Lf ", b.mu[i][j]);
            puts("");
        }
    }else{
        printf("Flag erroe: flag %s is not defined.\n", flag);
    }
}

/// @brief The function ```random_lattice``` generates random ```nrows```-dimensional lattice basis.
/// @param nrows Rank of lattice
/// @param ncols 
/// @return 
lattice random_lattice(int nrows, int ncols){
    lattice b;
    srand((unsigned int)time(NULL));
    b.nrows = nrows;
    b.ncols = ncols;
    b.basis = (long **)malloc(nrows * sizeof(long *));
    for(int i = 0, j; i < nrows; ++i){
        b.basis[i] = (long *)malloc(ncols * sizeof(long));
        for(j = 0; j < ncols; ++j){
            b.basis[i][i] = 1;
            b.basis[i][0] = rand() % 100;
        }
    }
    return b;
}


/**
 * @brief The function ```GSO``` computes GSO-informations of lattice basis ```b.basis```.
 * 
 * @param b lattice
 * @return lattice computed GSO-informations.
 */
lattice GSO(lattice b){
    int i, j, k;
    b.mu = (long double **)malloc(b.nrows * sizeof(long double *));
    b.basis_star = (long double **)malloc(b.nrows * sizeof(long double *));
    b.B = (long double *)malloc(b.nrows * sizeof(long double));
    for(i = 0; i < b.nrows; ++i){
        b.basis_star[i] = (long double *)malloc(b.ncols * sizeof(long double));
        b.mu[i] = (long double *)malloc(b.nrows * sizeof(long double));
    }

    for(i = 0; i < b.nrows; ++i){
        b.mu[i][i] = 1.0;
        for(j = 0; j < b.nrows; ++j) b.basis_star[i][j] = b.basis[i][j];
        for(j = 0; j < i; ++j){
            b.mu[i][j] = dot_int_dbl(b.basis[i], b.basis_star[j], b.ncols) / dot_dbl_dbl(b.basis_star[j], b.basis_star[j], b.ncols);
            for(k = 0; k < b.ncols; ++k) b.basis_star[i][k] -= b.mu[i][j] * b.basis_star[j][k];
        }
        b.B[i] = dot_dbl_dbl(b.basis_star[i], b.basis_star[i], b.ncols);
    }
    return b;
}


void SizeReduce(long **b, long double **mu, const int i, const int j, const int m){
    int k;
    if(mu[i][j] > 0.5 || mu[i][j] < -0.5){
        const int q = round(mu[i][j]);
        for(k = 0; k < m; ++k) b[i][k] -= q * b[j][k];
        for(k = 0; k <= j; ++k) mu[i][k] -= mu[j][k] * q;
    }
}


/**
 * @brief 
 * 
 * @param b 
 * @param d 
 * @param n 
 * @param m 
 */
lattice LLL(lattice b, const double d){
    double nu, BB, C, t;
    b = GSO(b);
    
    for(int k = 1, tmp, i, j, h; k < b.nrows;){
        h = k - 1;
        for(j = h; j > -1; --j) SizeReduce(b.basis, b.mu, k, j, b.ncols);

        if(k > 0 && b.B[k] < (d - b.mu[k][h] * b.mu[k][h]) * b.B[h]){
            for(i = 0; i < b.ncols; ++i){tmp = b.basis[h][i]; b.basis[h][i] = b.basis[k][i]; b.basis[k][i] = tmp;}
            
            nu = b.mu[k][h]; BB = b.B[k] + nu * nu * b.B[h]; C = 1.0 / BB;
            b.mu[k][h] = nu * b.B[h] * C; b.B[k] *= b.B[h] * C; b.B[h] = BB;

            for(i = 0; i < h; ++i){
                t = b.mu[h][i]; b.mu[h][i] = b.mu[k][i]; b.mu[k][i] = t;
            }
            for(i = k + 1; i < b.nrows; ++i){
                t = b.mu[i][k]; b.mu[i][k] = b.mu[i][h] - nu * t;
                b.mu[i][h] = t + b.mu[k][h] * b.mu[i][k];
            }
            
            k = h;
        }else ++k;
    }
    return b;
}
