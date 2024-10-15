#include <stdio.h>
#include <stdlib.h>

#include "clat.h"

#define RANK 60

int main(){
    int r, *v;
    long double *w;
    long *y;
    y = (long *)malloc(RANK * sizeof(long));
    v = (int *)malloc(RANK * sizeof(int));
    w = (long double *)malloc(RANK * sizeof(long double));
    lattice b = random_lattice(RANK, RANK);

    // Generates random lattice basis
    puts("random basis");
    print(b, "basis");

    // Computes the shortest vector on the lattive basis
    puts("\nThe shortest vector");
    enumerate(b.mu, b.B, RANK, v);
    int* x = coef2lat(v, b);
    for(int i = 0; i < RANK; ++i) printf("%d ", x[i]);
    puts("");

    // LLL
    puts("\nLLL-reduced basis");
    b = LLL(b, 0.99);
    print(b, "basis");

    // DeepLLL
    puts("\nDeepLLL-reduced basis");
    b = DeepLLL(b, 0.99);
    print(b, "basis");

    // NanchatteBKZ
    puts("\n\"Nanchatte\"BKZ-reduced basis");
    b = NanchatteBKZ(b, 40, 0.99);
    print(b, "basis");

    // Babai's nearest plane algorithm
    puts("\nTarget vector");
    for(int i = 0; i < RANK; ++i){
        w[i] = rand();
        printf("%ld ", (long)w[i]);
    }
    Babai(b, w, y);
    puts("\n\nBabai's vector");
    for(int i = 0; i < RANK; ++i) printf("%ld ", y[i]);

    puts("\n\nDeference vector");
    for(int i = 0; i < RANK; ++i) printf("%ld ", y[i] - (long)w[i]);
    puts("");
    return 0;
}
