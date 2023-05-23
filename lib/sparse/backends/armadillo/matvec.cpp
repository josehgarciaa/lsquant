#include <armadillo>
#include <complex>
#include <cstdio>
extern "C" {
struct complex_struct {
    double real;
    double imag;
};}

extern "C" {
void matvec_mult(complex_struct* values, int* col_indices, int* row_ptr, int nrows, complex_struct* vector, complex_struct* result) {
    arma::cx_vec v(nrows);
    for (int i = 0; i < nrows; i++) 
        v[i] = std::complex<double>(vector[i].real, vector[i].imag);

    arma::sp_cx_mat A(nrows, nrows);
    for (int i = 0; i < nrows; i++) 
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) 
            A(i, col_indices[j]) = std::complex<double>(values[i].real, values[i].imag);

    arma::cx_vec r = A * v;
    for (int i = 0; i < nrows; i++) {
        result[i].real = r[i].real();
        result[i].imag = r[i].imag();
    }
}

} // extern "C"

extern "C" {
void sparse_dense_matmult(complex_struct* values, int* col_indices, int* row_ptr, int nrows, complex_struct* matrix, int mcols, complex_struct* result) {
    arma::cx_mat arma_matrix(nrows, mcols);
    for (int i = 0; i < nrows*mcols; i++) 
        arma_matrix[i] = std::complex<double>(matrix[i].real, matrix[i].imag);

    arma::sp_cx_mat arma_spmat(nrows, nrows);
    for (int i = 0; i < nrows; i++) 
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) 
            arma_spmat(i, col_indices[j]) = std::complex<double>(values[i].real, values[i].imag);

    arma::cx_mat arma_result = arma_spmat * arma_matrix;
    for (int i = 0; i < nrows; i++)
        for (int j = 0; j < mcols; j++) {
            result[i].real = arma_result(i,j).real();
            result[i].imag = arma_result(i,j).imag();
            }
}

} // extern "C"