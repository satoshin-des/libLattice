#include <stdio.h>

#include "clat.h"

int main(){
    lattice b = random_lattice(3, 3);
    print(b, "basis");
    b = GSO(b);
    print(b, "GSO");
    b = LLL(b, 0.99);
    print(b, "basis");
    return 0;
}