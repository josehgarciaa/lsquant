#include <iostream>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "ctypes_utils.hpp"

#include <cblas.h>

extern "C" 
std::complex<double> 
dot(std::complex<double>* x, std::complex<double>* y, int size) 
{
    std::complex<double>
    norm = cblas_ddot(size,
                      reinterpret_cast<double*>(x), 1,
                      reinterpret_cast<double*>(y), 1);
    return norm;
}

typedef std::complex<double> Complex;
typedef Eigen::SparseMatrix< Complex ,Eigen::RowMajor > EigSpMatrix;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic> EigMatrix;


extern "C" 
void 
ctype_to_eigenspmat(ctype_util::csr x, EigSpMatrix X)
{
    X = EigSpMatrix(x.dim, x.dim);
    X.reserve(x.nnz);
    for (int i = 0; i < x.dim; ++i)
        for (int j = x.row_ptr[i]; j < x.row_ptr[i+1]; ++j)
        X.insert(i, x.col_indices[j]) = Complex(x.values[j].real,x.values[j].imag);
    X.finalize();
    //free(x.row_ptr); free(x.col_indices); free(x.values);
};

extern "C" 
void chebyshev_density( int num_moms, 
                        ctype_util::csr ham,
                        ctype_util::csr op, 
                        ctype_util::complex* cheb_vecs)
{

    EigSpMatrix eigen_ham;
    ctype_to_eigenspmat(ham, eigen_ham);
   
    EigSpMatrix eigen_op;
    ctype_to_eigenspmat(op, eigen_op);

    // Convert the three Chebyshev vectors into a matrix
    const int nvec = 3;
    EigMatrix eig_cheb_vecs(ham.dim,nvec);
    for (int i = 0; i < ham.dim; ++i)
        for (int j = 0; j < nvec; ++j) 
        {
            ctype_util::complex x = cheb_vecs[j * ham.dim + i];
            eig_cheb_vecs(i, j) = Complex(x.real,x.imag);
        }

    // Perform sparse-dense matrix multiplication
    EigMatrix result(ham.dim,nvec);
    for (int m=0; m < num_moms; m++)
    {
        Complex* Tm = eig_cheb_vecs.data();
        //for (int n=0; n < nvec; n++)
            //Complex mom = dot(&Tm[n*ham.dim],&Tm[n*ham.dim] , ham.dim);

        // result = eigen_ham * eig_cheb_vecs;
        // eig_cheb_vecs = result;
        
    }
}
