#include <iostream>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "ctypes_utils.hpp"

#include <cblas.h>

extern "C" 
std::complex<double> 
dot(const std::complex<double>* x,const  std::complex<double>* y, int size) 
{
    std::complex<double>
    norm = cblas_ddot(size,
                      reinterpret_cast<const double*>(x), 1,
                      reinterpret_cast<const double*>(y), 1);
    return norm;
}

typedef std::complex<double> Complex;
typedef Eigen::SparseMatrix< Complex ,Eigen::RowMajor > EigSpMatrix;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic> EigMatrix;


extern "C" 
void 
ctype_to_eigenspmat(ctype_util::csr x, EigSpMatrix& X, bool release_memory=false)
{
    X = EigSpMatrix(x.dim, x.dim);
    X.reserve(x.nnz);
    for (int i = 0; i < x.dim; ++i)
        for (int j = x.row_ptr[i]; j < x.row_ptr[i+1]; ++j)
        X.insert(i, x.col_indices[j]) = Complex(x.values[j].real,x.values[j].imag);
    X.finalize();
    if(release_memory)
    {
        printf("Realing memory\n");
        free(x.row_ptr); 
        free(x.col_indices); 
        free(x.values);
    }
};

void 
ctype_to_eigenmat(int nrows, int ncols, 
                  ctype_util::complex* x, EigMatrix& X, 
                  bool release_memory=false)
{
    // Convert the three Chebyshev vectors into a matrix
    X = EigMatrix(nrows,ncols);
    for (int i = 0; i < nrows; ++i)
        for (int j = 0; j < ncols; ++j) 
        {
            auto idx = i * ncols + j;
            X(i, j) = Complex(x[idx].real,x[idx].imag);
        }
    if(release_memory)
    {
        printf("Realing memory\n");
        free(x); 
    }
};


extern "C" 
void  chebyshev_density( int num_moms, 
                        ctype_util::csr ham,
                        ctype_util::csr op, 
                        ctype_util::complex* cheb_vecs)
{
    const int dim = ham.dim;

    EigSpMatrix eigen_ham;
    ctype_to_eigenspmat(ham, eigen_ham);
   
    EigSpMatrix eigen_op;
    ctype_to_eigenspmat(op, eigen_op);

    // Convert the three Chebyshev vectors into a matrix
    const int nvec = 3;
    EigMatrix eig_cheb_vecs;
    ctype_to_eigenmat(dim, nvec, cheb_vecs,eig_cheb_vecs);

    // Perform sparse-dense matrix multiplication
    EigMatrix result(dim,nvec);

    const int eff_num_moms = int(num_moms/2)*2; 
    auto cheb_momemts = new std::vector<Complex>(eff_num_moms);

    for (int m=0; 2*m < num_moms; m++)
    {
        const Complex* Tm0 = eig_cheb_vecs.data()+0*dim;
        const Complex* Tm1 = eig_cheb_vecs.data()+1*dim;
        (*cheb_momemts)[2*m+0] = dot(Tm0,Tm0,dim);
        (*cheb_momemts)[2*m+1] = dot(Tm0,Tm1,dim);
        result = eigen_ham * eig_cheb_vecs;
        eig_cheb_vecs = result;
    }
}
