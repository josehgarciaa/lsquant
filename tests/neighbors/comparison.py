import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer
from numpy2ctype import make_ndpointer


# Load the shared library into ctypes
lib = ctypes.CDLL("../../lib/libneighbors.so")

# define the numpy array types that will be used as arguments to the function
array_1d_bool = ndpointer(ctypes.c_bool, flags="C_CONTIGUOUS")
array_1d_double = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
array_2d_double = ndpointer(ctypes.c_double, ndim=2, flags="C_CONTIGUOUS")
array_2d_int = ndpointer(ctypes.c_int, ndim=2, flags="C_CONTIGUOUS")

# set argument types
lib.compute_neighbor_list.argtypes = [ctypes.c_int, 
                                      array_1d_bool,
                                      array_1d_double,array_1d_double,array_1d_double,
                                      ctypes.c_int, array_1d_double,
                                      ctypes.c_double]
# set return type
lib.compute_neighbor_list.restype = None

# Now you can call your function
pbc = np.array([True, False, True], dtype=bool)
lat0 = np.array([10.0, 0.0, 0.0], dtype=float)
lat1 = np.array([0.0, 10.0, 0.0], dtype=float)
lat2 = np.array([0.0, 0.0, 10.0], dtype=float)
positions = np.array([[1.1, 2.0,30], [3.0, 4.0,10]], dtype=float)
cutoff_radius = 1.0
neighbor_list = np.empty((2, 2), dtype=int)

lib.compute_neighbor_list(3,pbc, lat0, lat1,lat2,len(positions), positions.flatten(), cutoff_radius )


"""

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
"""