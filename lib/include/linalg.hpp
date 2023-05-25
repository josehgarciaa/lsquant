#include <cblas.h>

extern "C" {
std::complex<double> 
dot(const std::complex<double>* x,const  std::complex<double>* y, int size) 
{
    std::complex<double>
    norm = cblas_ddot(size,
                      reinterpret_cast<const double*>(x), 1,
                      reinterpret_cast<const double*>(y), 1);
    return norm;
}}