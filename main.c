#include <stdio.h>

#include "clat.h"

int main(){
    lattice b = random_lattice(10, 10);
    puts("random basis");
    print(b, "basis");

    puts("\nThe shortest vector");
    int* v = enumerate(b.mu, b.B, 10);
    int* x = coef2lat(v, b);
    for(int i = 0; i < 10; ++i) printf("%d ", x[i]);
    puts("");

    puts("\nDeepLLL-reduced basis");
    b = BKZ(b, 5, 0.99, 10);
    print(b, "basis");
    return 0;
}
