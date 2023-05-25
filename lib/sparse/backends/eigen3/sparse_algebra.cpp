#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef std::complex<double> Complex;
typedef Eigen::SparseMatrix< Complex ,Eigen::RowMajor > EigSpMatrix;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic> EigMatrix;
typedef Eigen::Map<Eigen::Matrix<Complex, Eigen::Dynamic, 1>> EigMap;
typedef Eigen::VectorX<Complex> EigVector;

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
#include "ctypes_utils.hpp"



extern "C"{
EigSpMatrix 
convert_to_eigen(const ctype_csr* x)
{
    printf("%d %d ", x->dim, x-> nnz);
    EigSpMatrix X(x->dim, x->dim);
    X.reserve(x->nnz);
    for (int i = 0; i < x->dim; ++i)
        for (int j = x->row_ptr[i]; j < x->row_ptr[i+1]; ++j)
        X.insert(i, x->col_indices[j]) = x->values[j];
    X.makeCompressed();
    return X;
}};


extern "C"{
void  chebyshev_density( int num_moms, 
                         const ctype_csr* ham,
                         const ctype_csr* op, 
                         const ctype_vector* T0,
                         ctype_vector* cheb_moms)
{
    const int dim = ham->dim;
    
    const EigSpMatrix eigen_ham = convert_to_eigen(ham);
    const EigSpMatrix eigen_op = convert_to_eigen(op);
   
    // Convert the three Chebyshev vectors into a matrix
    EigVector Tm0(dim);
    EigVector Tm1(dim);

    // Reinterpret T0 values as an eigen Vector by reference
    // and pass it to the columns of the cheb_vecs 
    // for iterations
    EigMap eigenT0(T0->values, dim);
    Tm0 = eigenT0;
    Tm1 = eigen_ham * eigenT0;

    for (int m=0; m < 2; m++)
    {
        printf("\n" );
        for (int i = 0 ; i < dim; i++)
            printf("%f %f \n", (Tm0.data())[i].real(), (Tm0.data())[i].imag() );

        cheb_moms->values[m] = dot(T0->values,Tm0.data(),dim);
//        Tm0 = 2*eigen_ham * Tm1 - Tm0;
        Tm1 = eigen_ham * Tm0;
        Tm0.swap(Tm1);
    }
}};
