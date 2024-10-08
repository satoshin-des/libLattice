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

/**
 * @brief Prints lattice information.
 * 
 * @param b 
 * @param flag 
 */
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
        printf("Flag error: flag %s is not defined.\n", flag);
    }
}

/**
 * @brief Generates random ```nrows```-dimensional lattice basis.
 * 
 * @param nrows 
 * @param ncols 
 * @return lattice 
 */
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
    b = GSO(b);
    return b;
}


int *coef2lat(int* v, lattice b){
    int *x = (int *)malloc(b.nrows * sizeof(int));
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


/**
 * @brief Computes volume (or determinat) of lattice
 * 
 * @param b lattice basis
 * @return long double 
 */
long double vol(lattice b){
    long double p = 1.0;
    b = GSO(b);
    for(int i = 0; i < b.nrows; ++i) p *= b.B[i];
    return sqrt(p);
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
 * @brief LLL-reduces the lattice basis ``b``.
 * 
 * @param b lattice basis
 * @param d reduction parameter
 * @return lattice
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

/**
 * @brief Enumerates a lattice vector whose squared norm is shorter than ``R``.
 * 
 * @param mu 
 * @param B 
 * @param n 
 * @param R 
 * @return int* 
 */
int *ENUM(long double** mu, long double* B, const int n, const double R) {
    int n1 = n + 1;
    int i, j, *r, *w, *v;
    double tmp;
    long double *c, *rho, **sigma;

    r = (int *)malloc(n1 * sizeof(int));
    w = (int *)malloc(n * sizeof(int));
    v = (int *)malloc(n * sizeof(int));
    c = (long double *)malloc(n * sizeof(long double));
    rho = (long double *)malloc(n1 * sizeof(long double));
    sigma = (long double **)malloc(n1 * sizeof(long double *));
    for (i = 0; i < n; ++i){
        r[i] = i;
        w[i] = v[i] = 0;
        c[i] = rho[i] = 0;
        sigma[i] = (long double *)malloc(n * sizeof(long double));
        for(j = 0; j < n; ++j) sigma[i][j] = 0;
    }
    v[0] = 1;
    r[n] = 0;
    rho[n] = 0;
    sigma[n] = (long double *)malloc(n * sizeof(long double));
    for(j = 0; j < n; ++j) sigma[n][j] = 0;

    for (int k = 0, last_nonzero = 0; ;) {
        printf("%d\n", k);
        tmp = c[k] - v[k]; tmp *= tmp;
        rho[k] = rho[k + 1] + tmp * B[k];
        if (rho[k] <= R) {
            if (k == 0){
                free(r); free(w); free(c); free(rho); free(sigma);
                return v;
            }else{
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k]: r[k + 1]);
                for (i = r[k]; i > k; --i) sigma[i][k] = sigma[i + 1][k] + mu[i][k] * v[i];
                c[k] = -sigma[k + 1][k];
                v[k] = round(c[k]);
                w[k] = 1; // 小さいやつが見つかったら、変分を元に戻す
            }
        }else{
            ++k;
            if (k == n) {
                free(r); free(w); free(c); free(rho); free(sigma);
                free(v); v = NULL;
                return v;
            }else{
                r[k] = k;
                if (k >= last_nonzero) {
                    last_nonzero = k;
                    ++v[k];
                }
                else {
                    if(v[k] > c[k]) v[k] -= w[k]; else v[k] += w[k];
                    ++w[k];
                }
            }
        }
    }
}


int *enumerate(long double **mu, long double *B, const int n) {
    int i, *enum_v, *pre_enum_v;
    enum_v = (int *)malloc(n * sizeof(int));
    pre_enum_v = (int *)malloc(n * sizeof(int));
    for(i = 0; i < n; ++i) enum_v[i] = pre_enum_v[i] = 0;
    for (double R = B[0];; R *= 0.99) {
        for(i = 0; i < n; ++i) pre_enum_v[i] = enum_v[i];
        enum_v = ENUM(mu, B, n, R);
        if (enum_v == NULL){
            free(enum_v);
            return pre_enum_v;
        }
    }
}