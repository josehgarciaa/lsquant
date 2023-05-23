from ctypes import CDLL, c_int, c_double, POINTER, Structure
import numpy as np

import numpy as np
from ctypes import c_int, c_float, c_double, POINTER
from numpy.ctypeslib import as_ctypes

# Define the mapping from numpy data types to ctypes data types
dtype_dict = {
    np.dtype(np.int32).name: c_int,
    "int": c_int,
    np.dtype(np.float32).name: c_float,
    np.dtype(np.float64).name: c_double
}

class ComplexStruct(Structure):
    _fields_ = [('real', c_double), ('imag', c_double)]

class CSRStruct(Structure):
    _fields_ = [
        ("dim", c_int),
        ("nnz", c_int),
        ("pindices", POINTER(c_int)),
        ("pindptr", POINTER(c_int)),
        ("pvalues",POINTER(ComplexStruct))
    ]



def copy_nparray(arr):
    """
    Convert a numpy array to a ctypes array.
    """
    if arr.dtype.name not  in  dtype_dict.keys():
        raise ValueError(f"Unsupported data type: {arr.dtype}")
    ctype = dtype_dict[arr.dtype.name]
    return arr.ctypes.data_as(POINTER(ctype))

def copy_complex_nparray(arr):
    arr_complex = np.empty((arr.shape[0],), dtype=[('real', np.float64), ('imag', np.float64)])
    arr_complex['real'] = arr.real
    arr_complex['imag'] = arr.imag
    #The Indicates that an array of size A_complex.shape[0]
    # for data type complex should be created
    return  (ComplexStruct * arr_complex.shape[0]).from_buffer(arr_complex)