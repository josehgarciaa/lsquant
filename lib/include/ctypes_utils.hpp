

extern "C"{
    struct ctype_csr {
    int dim;
    int nnz;
    int* row_ptr;
    int* col_indices;
    std::complex<double>* values; 
};}

extern "C"{
    struct ctype_matrix {
    int nrows;
    int ncols;
    std::complex<double>* values; 
};}

extern "C"{
    struct ctype_vector {
    int dim;
    std::complex<double>* values; 
};}
