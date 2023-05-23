import numpy as np
import scipy.sparse as sp
from ctypes import CDLL, c_int, c_double, POINTER, Structure
from numpy.ctypeslib import as_ctypes

from numpy2ctype import copy_complex_nparray, copy_nparray, ComplexStruct, CSRStruct


# Load the shared library into ctypes
lib = CDLL("lib/eigen3_spalg.so")

lib.chebyshev_density.argtypes = [ c_int, CSRStruct, CSRStruct, POINTER(ComplexStruct)]
lib.chebyshev_density.restype = None


# Create a sparse matrix and a complex vector
Ham= sp.csr_matrix([[0.0+2.j, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.+2.j]], dtype=complex)
cheb_vecs  = np.eye(3, dtype=complex).flatten()


#Convert numpy array to c_types
nnz   = c_int(Ham.nnz)
dim = c_int(Ham.shape[0])
col_indices = copy_nparray(Ham.indices)
row_ptr = copy_nparray(Ham.indptr)
values = copy_complex_nparray(Ham.data)
c_cheb_vecs = copy_complex_nparray(cheb_vecs)


csr_tuple = CSRStruct(dim, nnz, col_indices, row_ptr, values )

ntimes = c_int(1)
# Perform the matrix-vector multiplication
lib.chebyshev_density(ntimes, csr_tuple, csr_tuple, c_cheb_vecs)

# Convert result back to numpy array
cheb_vecs = np.array([complex(r.real, r.imag) for r in c_cheb_vecs])
print(f"The result of the matrix-vector multiplication is: {cheb_vecs}")
