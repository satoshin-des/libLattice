#include <stdio.h>

#include "clat.h"

#define RANK 10

int main(){
    int r;
    lattice b = random_lattice(10, 10);
    puts("random basis");
    print(b, "basis");

    puts("\nThe shortest vector");
    int* v = enumerate(b.mu, b.B, 10);
    int* x = coef2lat(v, b);
    for(int i = 0; i < 10; ++i) printf("%d ", x[i]);
    puts("");

    puts("\nDeepLLL-reduced basis");
    //b = DeepLLL(b, 0.99);
    b = BKZ(b, 5, 0.99, 10);
    //printf("%d\n", b.ncols);
    all_information(b);
    //print(b, "basis");
    return 0;
}
