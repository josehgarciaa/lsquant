#include <cmath>
#include <stdio.h>

extern "C" {

struct complex_struct {
    double real;
    double imag;
};

double compute_norm(struct complex_struct* complex_array, int length) {
    double norm = 0;
    printf("The dimension of the array is %d\n",length );

    for(int i=0; i<length; i++) {
        norm += std::sqrt(std::pow(complex_array[i].real, 2) + std::pow(complex_array[i].imag, 2));
    }
    return norm;
}

} // extern "C"