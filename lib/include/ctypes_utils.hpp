
namespace ctype_util{

extern "C" 
struct complex {
    double real;
    double imag;
};

extern "C" 
struct csr {
    int dim;
    int nnz;
    int* row_ptr;
    int* col_indices;
    complex* values; 
};

};