from ctypes import CDLL, c_int, c_double, POINTER, Structure
import numpy as np

import numpy as np
from ctypes import c_int, c_float, c_double, POINTER
from numpy.ctypeslib import as_ctypes, ndpointer

def make_ndpointer(array, dtype):
    array = np.ascontiguousarray(array, dtype=dtype)
    return array.ctypes.data_as(ndpointer(dtype=dtype, ndim=1, flags='C_CONTIGUOUS, ALIGNED'))
class CSRStruct(Structure):
    _fields_ = [
                ("dim", c_int),
                ("nnz", c_int),
                ("pindices",ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS, ALIGNED')),
                ("pindptr", ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS, ALIGNED')),
                ("data", ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS, ALIGNED'))
                 ]
class MatrixStruct(Structure):
    _fields_ = [
                ("nrows", c_int),
                ("ncols", c_int),
                ("data", ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS, ALIGNED'))
                 ]
        
class VectorStruct(Structure):
    _fields_ = [
                ("dim", c_int),
                ("data", ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS, ALIGNED'))
                 ]
