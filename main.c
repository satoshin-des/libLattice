#include <stdio.h>
#include <stdlib.h>

#include "clat.h"

#define RANK 60

int main(){
    int r, *v; v = (int *)malloc(RANK * sizeof(int));
    lattice b = random_lattice(RANK, RANK);
    puts("random basis");
    print(b, "basis");

    puts("\nThe shortest vector");
    enumerate(b.mu, b.B, RANK, v);
    int* x = coef2lat(v, b);
    for(int i = 0; i < RANK; ++i) printf("%d ", x[i]);
    puts("");

    puts("\nLLL-reduced basis");
    b = LLL(b, 0.99);
    print(b, "basis");

    puts("\n\"Nanchatte\"BKZ-reduced basis");
    //b = DeepLLL(b, 0.99);
    b = NanchatteBKZ(b, 40, 0.99);
    //all_information(b);
    print(b, "basis");
    return 0;
}
