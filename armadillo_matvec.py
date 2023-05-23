import numpy as np
import scipy.sparse as sp
from ctypes import CDLL, c_int, c_double, POINTER, Structure
from numpy.ctypeslib import as_ctypes

from numpy2ctype import copy_complex_nparray, copy_nparray, Complex


# Load the shared library into ctypes
lib = CDLL("lib/armadillo_matvec.so")

lib.sparse_dense_mult.argtypes = [c_int, c_int, c_int, 
                            POINTER(Complex), POINTER(c_int), POINTER(c_int), 
                            c_int, c_int,
                            POINTER(Complex), POINTER(Complex)]

lib.matvec_mult.restype = None


lib.sparse_dense_matmult.argtypes = [POINTER(Complex), POINTER(c_int), POINTER(c_int), c_int, POINTER(Complex),c_int, POINTER(Complex)]
lib.sparse_dense_matmult.restype = None

# Create a sparse matrix and a complex vector
A = sp.csr_matrix([[0.0+2.j, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.+2.j]], dtype=complex)
v = np.array([+2.j, 1.0j, 1.0+0.0j,], dtype=complex)

# Create ctypes-compatible data
#values = np.array(A.data, dtype=c_double)
col_indices = copy_nparray(A.indices)
row_ptr = copy_nparray(A.indptr)
values = copy_complex_nparray(A.data)
vector = copy_complex_nparray(v)
result = copy_complex_nparray(0.0*v)
nrows = c_int(A.shape[0])

# Perform the matrix-vector multiplication
lib.matvec_mult(values, col_indices, row_ptr, nrows, vector, result)
# Convert result back to numpy array
result_np = np.array([complex(r.real, r.imag) for r in result])
print(f"The result of the matrix-vector multiplication is: {result_np}")
lib.sparse_dense_matmult(values, col_indices, row_ptr, nrows, vector,1, result)

# Convert result back to numpy array
result_np = np.array([complex(r.real, r.imag) for r in result])

print(f"The result of the matrix-vector multiplication is: {result_np}")