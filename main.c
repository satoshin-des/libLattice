#include <stdio.h>

#include "clat.h"

int main(){
    lattice b = random_lattice(3, 3);
    print(b, "basis");
    int* v = enumerate(b.mu, b.B, 3);
    for(int i = 0; i < 3; ++i) printf("%d ", v[i]);
    puts("");
    return 0;
}