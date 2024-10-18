#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "clat.h"

/**
 * @brief Tests if the vector is zero or not.
 * 
 * @param v vector
 * @param n size of vector
 * @return true 
 * @return false 
 */
bool isZero(const int* v, const int n){
    for(int i = 0; i < n; ++i) if(v[i] != 0) return false;
    return true;
}

bool IsZero(const long* v, const int n){
    for(int i = 0; i < n; ++i) if(v[i] != 0) return false;
    return true;
}

/**
 * @brief Computes inner product value of vectors ``x`` and ``y``.
 * 
 * @param x vector
 * @param y vector
 * @param n size of vectors
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
 * @param n size of vectors
 * @return double 
 */
double dot_int_dbl(long *x, long double *y, const int n){
    double s = 0.0;
    for(int i = 0; i < n; ++i) s += y[i] * x[i];
    return s;
}
/**
 * @brief Computes inner product value of vectors ``x`` and ``y``.
 * 
 * @param x vector
 * @param y vector
 * @param n size of vectors
 * @return long 
 */
long dot_int_int(long *x, long *y, const int n){
    long s = 0;
    for(int i = 0; i < n; ++i) s += y[i] * x[i];
    return s;
}

lattice _GSO(lattice b){
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
 * @brief Size-reduces lattice basis and updates GSO-information.
 * 
 * @param b lattice basis
 * @param mu GSO-coefficient matrix
 * @param i 
 * @param j 
 * @param m 
 */
void SizeReduce(long **b, long double **mu, const int i, const int j, const int m){
    if(mu[i][j] > 0.5 || mu[i][j] < -0.5){
        int k;
        const int q = round(mu[i][j]);
        for(k = 0; k < m; ++k) b[i][k] -= q * b[j][k];
        for(k = 0; k <= j; ++k) mu[i][k] -= mu[j][k] * q;
    }
}

lattice _LLL(lattice b, const float d, const int start, const int end){
    double nu, BB, C, t;
    int n = end - start;
    b = GSO(b);
    
    for(int k = 1, tmp, i, j, h; k < n;){
        h = k - 1;
        for(j = h; j > -1; --j) SizeReduce(b.basis, b.mu, k, j, b.ncols);

        if(k > 0 && b.B[k] < (d - b.mu[k][h] * b.mu[k][h]) * b.B[h]){
            for(i = 0; i < b.ncols; ++i){tmp = b.basis[h][i]; b.basis[h][i] = b.basis[k][i]; b.basis[k][i] = tmp;}
            
            nu = b.mu[k][h]; BB = b.B[k] + nu * nu * b.B[h]; C = 1.0 / BB;
            b.mu[k][h] = nu * b.B[h] * C; b.B[k] *= b.B[h] * C; b.B[h] = BB;

            for(i = 0; i < h; ++i){
                t = b.mu[h][i]; b.mu[h][i] = b.mu[k][i]; b.mu[k][i] = t;
            }
            for(i = k + 1; i < n; ++i){
                t = b.mu[i][k]; b.mu[i][k] = b.mu[i][h] - nu * t;
                b.mu[i][h] = t + b.mu[k][h] * b.mu[i][k];
            }
            
            k = h;
        }else ++k;
    }
    return b;
}

double _ENUM(long double** mu, long double* B, const double R, int* v, int qq, const int start, const int end) {
    const int n = end - start;
    const int n1 = n + 1;
    int i, j, *r, *w;
    double tmp;
    long double *c, *rho, **sigma;

    r = (int *)malloc(n1 * sizeof(int));
    w = (int *)malloc(n * sizeof(int));
    c = (long double *)malloc(n * sizeof(long double));
    rho = (long double *)malloc(n1 * sizeof(long double));
    sigma = (long double **)malloc(n1 * sizeof(long double *));
    for (i = 0; i < n; ++i){
        r[i] = i;
        c[i] = rho[i] = w[i] = v[i] = 0;
        sigma[i] = (long double *)malloc(n * sizeof(long double));
        for(j = 0; j < n; ++j) sigma[i][j] = 0;
    }
    v[0] = 1;
    rho[n] = r[n] = 0;
    sigma[n] = (long double *)malloc(n * sizeof(long double));
    for(j = 0; j < n; ++j) sigma[n][j] = 0;

    for (int k = 0, last_nonzero = 0; ;) {
        tmp = c[k] - v[k]; tmp *= tmp;
        rho[k] = rho[k + 1] + tmp * B[k + start];
        if (rho[k] < R || fabs(rho[k] - R) < 0.001) {
            if (k == 0){
                free(r); free(w); free(c);
                for(i = 0; i < n1; ++i) free(sigma[i]);
                free(sigma);
                tmp = rho[qq];
                free(rho);
                return tmp;
            }else{
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k]: r[k + 1]);
                for (i = r[k]; i > k; --i) sigma[i][k] = sigma[i + 1][k] + mu[i + start][k + start] * v[i];
                c[k] = -sigma[k + 1][k];
                v[k] = round(c[k]);
                w[k] = 1;
            }
        }else{
            ++k;
            if (k == n) {
                free(r); free(w); free(c);
                for(i = 0; i < n1; ++i) free(sigma[i]);
                free(sigma);
                for(i = 0; i < n; ++i) v[i] = 0;
                tmp = rho[qq];
                free(rho);
                return tmp;
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

/**
 * @brief 
 * 
 * @param mu 
 * @param B 
 * @param n 
 * @return int* 
 */
double _enumerate(long double **mu, long double *B, int *v, const int qq, const int start, const int end) {
    int i, *enum_v, *pre_enum_v, n = end - start;
    enum_v = (int *)malloc(n * sizeof(int));
    pre_enum_v = (int *)malloc(n * sizeof(int));
    for(i = 0; i < n; ++i) enum_v[i] = pre_enum_v[i] = 0;
    for (double R = B[start], pre_rho = 0, rho = 0;; R *= 0.99){
        for(i = 0; i < n; ++i) pre_enum_v[i] = enum_v[i];
        pre_rho = rho;
        rho = _ENUM(mu, B, R, enum_v, qq, start, end);
        if (isZero(enum_v, n)){
            for(i = 0; i < n; ++i) v[i] = pre_enum_v[i];
            free(pre_enum_v); free(enum_v);
            return pre_rho;
        }
    }
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

void all_information(lattice b){
    printf("nrows = %d, ncols = %d\n", b.nrows, b.ncols);
    puts("lattice basis matrix(spanned by row vectors) =");
    for(int i = 0, j; i < b.nrows; ++i){
        for(j = 0; j < b.ncols; ++j) printf("%ld ", b.basis[i][j]);
        puts("");
    }
}

lattice Lattice(long** b, const int n, const int m){
    lattice bb;
    srand((unsigned int)time(NULL));
    bb.nrows = n;
    bb.ncols = m;
    bb.basis = (long **)malloc(n * sizeof(long *));
    for(int i = 0, j; i < n; ++i){
        bb.basis[i] = (long *)malloc(m * sizeof(long));
        for(j = 0; j < m; ++j) bb.basis[i][j] = b[i][j];
    }
    bb = _GSO(bb);
    if(fabs(vol(bb)) < 0.999){
        puts("Lattice initialization Error: input matrix is not regular. Please input lattice basis.");
        exit(1);
    }
    return bb;
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
        b.basis[i][i] = 1;
        b.basis[i][0] = rand() % 900 + 100;
    }
    b = _GSO(b);
    return b;
}

/**
 * @brief Converts coefficient vector to lattice vector
 * 
 * @param v 
 * @param b 
 * @return int* 
 */
int *coef2lat(int* v, lattice b){
    int *x = (int *)malloc(b.ncols * sizeof(int)), i, j;
    for(i = 0; i < b.ncols; ++i){
        x[i] = 0;
        for(j = 0; j < b.nrows; ++j) x[i] += v[j] * b.basis[j][i];
    }
    return x;
}

/**
 * @brief 
 * 
 * @param v 
 * @param b 
 * @return long* 
 */
long *Coef2Lat(long* v, lattice b){
    long *x = (long *)malloc(b.ncols * sizeof(long)), i, j;
    for(i = 0; i < b.ncols; ++i){
        x[i] = 0;
        for(j = 0; j < b.nrows; ++j) x[i] += v[j] * b.basis[j][i];
    }
    return x;
}

/**
 * @brief The function ```GSO``` computes GSO-informations of lattice basis ```b.basis```.
 * 
 * @param b lattice
 * @return lattice computed GSO-informations.
 */
lattice GSO(lattice b){
    int i, j, k;
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
    for(int i = 0; i < b.nrows; ++i) p *= b.B[i];
    return sqrt(p);
}

/**
 * @brief Computes Lagrange(J.-L. Lagrange 1773)-reduced lattice basis
 * 
 * @param b 
 */
void Lagrange(lattice b){
    int q, i;
    long *v = (long *)malloc(b.ncols * sizeof(long));

    if(dot_int_int(b.basis[0], b.basis[0], b.ncols) > dot_int_int(b.basis[1], b.basis[1], b.ncols))
        for(i = 0; i < b.ncols; ++i){v[i] = b.basis[0][i]; b.basis[0][i] = b.basis[1][i]; b.basis[1][i] = v[i];}

    do{
        q = -round((double)dot_int_int(b.basis[0], b.basis[1], b.ncols) / dot_int_int(b.basis[0], b.basis[0], b.ncols));
        for(i = 0; i < b.ncols; ++i) v[i] = b.basis[1][i] + q * b.basis[0][i];
        for(i = 0; i < b.ncols; ++i){
            b.basis[1][i] = b.basis[0][i]; b.basis[0][i] = v[i];
        }
    }while(dot_int_int(b.basis[0], b.basis[0], b.ncols) < dot_int_int(b.basis[1], b.basis[1], b.ncols));
    
    for(i = 0; i < b.ncols; ++i){
        v[i] = b.basis[0][i]; b.basis[0][i] = b.basis[1][i]; b.basis[1][i] = v[i];
    }
}

/**
 * @brief Computes LLL(A. K. Lenstra, H. W. Lenstra, and L. Lovasz 1982)-reduced lattice basis.
 * 
 * @param b lattice basis
 * @param d reduction parameter
 * @return lattice
 */
lattice LLL(lattice b, const float d){
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
 * @brief Computes DeepLLL(C. P. Schnorr and M. Euchner 1994)-reduced lattice basis
 * 
 * @param b 
 * @param d 
 * @return lattice 
 */
lattice DeepLLL(lattice b, const float d){
    b = GSO(b);

	double C;
	for(int k = 1, j, i, t, l; k < b.nrows;){
		for(j = k - 1; j >= 0; --j) SizeReduce(b.basis, b.mu, k, j, b.ncols);

		C = dot_int_int(b.basis[k], b.basis[k], b.ncols);
		for(i = 0; i < k;){
			if(C >= d * b.B[i]){
				C -= b.mu[k][i] * b.mu[k][i] * b.B[i];
				++i;
			}else{
                for(l = 0; l < b.ncols; ++l){
				    t = b.basis[k][l];
				    for(j = k; j > i; --j) b.basis[j][l] = b.basis[j - 1][l];
				    b.basis[i][l] = t;
                }
                b = GSO(b);
				k = fmax(i - 1, 0);
			}
		}
		++k;
	}
    return b;
}

/**
 * @brief Computes PotLLL(F. Fontein, M. Schneider, and U. Wagner 2018)-reduced latttice basis.
 * 
 * @param b lattice basis
 * @param d reduction parameter
 * @return lattice 
 */
lattice PotLLL(lattice b, const float d){
    b = GSO(b);

    double P, P_min, S;
    long *t = (long *)malloc(b.ncols * sizeof(long));
    for(int l = 0, j, i, k, h; l < b.nrows;){
        h = l - 1;

        for(j = h; j > -1; --j) SizeReduce(b.basis, b.mu, l, j, b.ncols);

        P = P_min = 1.0; k = 0;
        for(j = h; j >= 0; --j){
            S = 0.0;
            for(i = j; i < l; ++i) S += b.mu[l][i] * b.mu[l][i] * b.B[i];
            P *= (b.B[l] + S) / b.B[j];
            if(P < P_min){k = j; P_min = P;}
        }
        
        if(d > P_min){
            t = b.basis[l];
            for(j = l; j > k; --j) b.basis[j] = b.basis[j - 1];
            b.basis[k] = t;

            b = GSO(b);
            l = k;
        }else ++l;
    }
    free(t);
    return b;
}

/**
 * @brief Computes BKZ-like reduced lattice basis. The algorithm in this function is not BKZ algorithm, but it behaves like BKZ algorithm.
 * 
 * @param b lattice basis
 * @param beta block-size
 * @param d reduction parameter
 * @return lattice 
 */
lattice NanchatteBKZ(lattice b, const int beta, const float d) {
    const int n1 = b.nrows - 1;
    int *v, *w;
    double s;
    w = (int *)malloc(b.nrows * sizeof(int));
    v = (int *)malloc(b.ncols * sizeof(int));

    b = LLL(b, d);

    for (int z = 0, j, t, BKZTour = 0, i, k = 0, h, lk1, l, m; z < n1;) {
        if(BKZTour >= 10) break;
        //printf("z = %d\n", z);

        if (k == n1){k = 0; ++BKZTour;} ++k;
        l = fmin(k + beta - 1, b.nrows); h = fmin(l + 1, b.nrows);
        lk1 = l - k + 1;

        s = _enumerate(b.mu, b.B, w, k - 1, k - 1, l);

        if (b.B[k - 1] > s && ! isZero(w, lk1)) {
            z = 0;

            for(i = 0; i < b.ncols; ++i){
                v[i] = 0;
                for(j = 0; j < lk1; ++j) v[i] += w[j] * b.basis[j + k - 1][i];
            }

            for(i = lk1 - 1; i >= 0; --i) if(w[i] != 0){m = i; break;}

            for(i = 0; i < b.ncols; ++i) b.basis[m + k - 1][i] = v[i];
            for(i = 0; i < b.ncols; ++i){
				t = b.basis[m + k - 1][i];
				for(j = m + k - 1; j > k - 1; --j) b.basis[j][i] = b.basis[j - 1][i];
				b.basis[k - 1][i] = t;
            }

            b = _LLL(b, d, 0, h);
            b = GSO(b);
        }else{
            ++z;

            b = _LLL(b, d, 0, h);
            b = GSO(b);
        }
    }
    free(v); free(w);
    return b;
}

/**
 * @brief Enumerates a lattice vector whose squared norm is shorter than ``R``.
 * 
 * @param mu GSO-coefficient matrix
 * @param B squared norms of GSO-vectors
 * @param n rank of lattice vector
 * @param R radius of enumeration
 * @param v integral coefficients vector of the lattice vector whose squared norm is shorter than ``R``
 */
void ENUM(long double** mu, long double* B, const int n, const double R, int* v) {
    int n1 = n + 1;
    int i, j, *r, *w;
    double tmp;
    long double *c, *rho, **sigma;

    r = (int *)malloc(n1 * sizeof(int));
    w = (int *)malloc(n * sizeof(int));
    c = (long double *)malloc(n * sizeof(long double));
    rho = (long double *)malloc(n1 * sizeof(long double));
    sigma = (long double **)malloc(n1 * sizeof(long double *));
    for (i = 0; i < n; ++i){
        r[i] = i;
        c[i] = rho[i] = w[i] = v[i] = 0;
        sigma[i] = (long double *)malloc(n * sizeof(long double));
        for(j = 0; j < n; ++j) sigma[i][j] = 0;
    }
    v[0] = 1;
    rho[n] = r[n] = 0;
    sigma[n] = (long double *)malloc(n * sizeof(long double));
    for(j = 0; j < n; ++j) sigma[n][j] = 0;

    for (int k = 0, last_nonzero = 0; ;) {
        tmp = c[k] - v[k]; tmp *= tmp;
        rho[k] = rho[k + 1] + tmp * B[k];
        if (rho[k] <= R) {
            if (k == 0){
                free(r); free(w); free(c); free(rho);
                for(i = 0; i < n1; ++i) free(sigma[i]);
                free(sigma);
                break;
            }else{
                --k;
                r[k] = (r[k] > r[k + 1] ? r[k]: r[k + 1]);
                for (i = r[k]; i > k; --i) sigma[i][k] = sigma[i + 1][k] + mu[i][k] * v[i];
                c[k] = -sigma[k + 1][k];
                v[k] = round(c[k]);
                w[k] = 1;
            }
        }else{
            ++k;
            if (k == n) {
                free(r); free(w); free(c); free(rho);
                for(i = 0; i < n1; ++i) free(sigma[i]);
                free(sigma);
                for(i = 0; i < n; ++i) v[i] = 0;
                break;
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

/**
 * @brief Enumerates the shortest lattice vectors.
 * 
 * @param mu GSO-coefficients matrix
 * @param B squared norms of GSO-vectors
 * @param n rank of lattice
 * @param v the shortest lattice vector
 */
void enumerate_SV(long double **mu, long double *B, const int n, int* v) {
    int i, *enum_v, *pre_enum_v;
    enum_v = (int *)malloc(n * sizeof(int));
    pre_enum_v = (int *)malloc(n * sizeof(int));
    for(i = 0; i < n; ++i) enum_v[i] = pre_enum_v[i] = 0;
    for (double R = B[0];; R *= 0.99) {
        for(i = 0; i < n; ++i) pre_enum_v[i] = enum_v[i];
        ENUM(mu, B, n, R, enum_v);
        if (isZero(enum_v, n)){
            for(i = 0; i < n; ++i) v[i] = pre_enum_v[i];
            free(enum_v); free(pre_enum_v);
            break;
        }
    }
}

/**
 * @brief Solves approx-CVP by Babai's nearest plane algorithm.
 * 
 * @param b lattice basis
 * @param w target vector
 * @param v close vector to ``w``
 */
void Babai(lattice b, long double* w, long* v){
	int i, j;
    long double *t, c;
    t = (long double *)malloc(b.ncols * sizeof(long double));
    for(i = 0; i < b.ncols; ++i) t[i] = w[i];

	for(i = b.nrows - 1; i >= 0; --i){
		c = round(dot_dbl_dbl(t, b.basis_star[i], b.ncols) / dot_dbl_dbl(b.basis_star[i], b.basis_star[i], b.ncols));
		for(j = 0; j < b.ncols; ++j) t[j] -= c * b.basis_star[i][j];
	}
	for(i = 0; i < b.ncols; ++i) v[i] = w[i] - t[i];
    free(t);
}
