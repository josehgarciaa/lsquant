import numpy as np
import scipy.sparse as sp
from ctypes import CDLL, c_int, c_double, POINTER, byref
from numpy2ctype import make_ndpointer, CSRStruct, MatrixStruct, VectorStruct

# Load the shared library into ctypes
lib = CDLL("lib/eigen3_spalg.so")

lib.chebyshev_density.argtypes = [ c_int, POINTER(CSRStruct), POINTER(CSRStruct), POINTER(VectorStruct)]
lib.chebyshev_density.restype =  None

# Create a sparse matrix and a complex vector
Ham= sp.csr_matrix([[1.0+np.sqrt(2)*1j, 1.0, 2.0], [0.0, np.pi, 0.787686], [.541560, 1j + .0, 1.+2.j]], dtype=np.complex128)
cheb_vecs  = np.eye(3, dtype=complex).flatten()

#Convert numpy array to c_types
csr = CSRStruct();
csr.dim = Ham.shape[0]
csr.nnz = Ham.nnz
csr.pindices= make_ndpointer(Ham.indices, dtype=np.int32)
csr.pindptr = make_ndpointer(Ham.indptr, dtype=np.int32)
csr.data    = make_ndpointer(Ham.data, dtype=np.complex128)

#We create a structure for the number of moments
cheb_vec0 = np.zeros(csr.dim).astype(np.complex128);
cheb_vec0[0] = 1.5
cheb_vec0[1] = 3
c_cheb_vec0 = VectorStruct()
c_cheb_vec0.dim= Ham.shape[0]
c_cheb_vec0.data = make_ndpointer( cheb_vec0, dtype=np.complex128 )

#Create the array that will storage the chebyshev array
nmom=10
cheb_moms=np.zeros(nmom, dtype=np.complex128)
c_cheb_moms = VectorStruct()
c_cheb_moms.dim= nmom
c_cheb_moms.data = make_ndpointer( cheb_moms, dtype=np.complex128 ) 

# Perform the matrix-vector multiplication
lib.chebyshev_density(nmom,  byref(csr),  byref(csr),  byref(c_cheb_vec0), byref(c_cheb_moms))
print(cheb_moms)
