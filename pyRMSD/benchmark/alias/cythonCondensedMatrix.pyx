#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
"""
Created on 12/07/2012

@author: victor
"""
import numpy
cimport numpy
import math
from libc.stdlib cimport malloc, free
cimport cython

numpy.import_array()

cdef inline int condensedsubscript(int row_length_minus_one,int i,int j) :
#    return i * (row_length-1) - i - ((( i-1)*i) /2 ) + j
    return i*(row_length_minus_one) - i - ((( i-1)*i) >> 1 ) + j

@cython.final
cdef class CythonCondensedMatrixWithDiagonal(object):
    cdef double *inner_data
    cdef public int row_length, size
    
    def __cinit__(self, data):
        cdef int i
        self.size = len(data)
        self.row_length = int((1 + math.sqrt(1+8*self.size))/2)
        # We add some more room for diagonal
        self.size += self.row_length 
        # Alloc. memory
        self.inner_data = <double *> malloc(self.size * sizeof(double))
        # Now we'll copy the data as if the origin was another smaller
        # condensed matrix, so we'll create one
        cond =  CythonCondensedMatrix(data)
        
        for i in range(self.row_length):
            for j in range(i,self.row_length):
                pos = condensedsubscript(self.row_length,i,j)
                self.inner_data[pos] = cond[i,j]
        
    
    def __init__(self, data):
        pass
            
    def __dealloc__(self):
        free(self.inner_data)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, tuple index):
        return self.getitem(index[0],index[1])
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double getitem(self,int i,int j):
        if i < j:
            return self.inner_data[i*(self.row_length) - i - ((( i-1)*i) >> 1 ) + j]
        else:
            return self.inner_data[j*(self.row_length) - j - ((( j-1)*j) >> 1 ) + i]
   
    cpdef get_data(self):
        cdef numpy.npy_intp shape[1]
        shape[0] = <numpy.npy_intp> self.size
        ndarray = numpy.PyArray_SimpleNewFromData(1, shape, numpy.NPY_DOUBLE, self.inner_data)
        return ndarray
        
@cython.final
cdef class CythonCondensedMatrix(object):
    
    cdef double *inner_data
    cdef public int row_length, row_length_minus_one, size
    
    def __cinit__(self, data):
        cdef int i
        self.size = len(data)
        self.inner_data = <double *> malloc(self.size * sizeof(double))
        for i in range(self.size):
            self.inner_data[i] = data[i]
        self.row_length = int((1 + math.sqrt(1+8*self.size))/2)
        self.row_length_minus_one = self.row_length - 1
    
    def __init__(self, data):
        pass
            
    def __dealloc__(self):
        free(self.inner_data)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, tuple index):
        return self.getitem(index[0],index[1])

        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef inline double getitem(self, int i, int j):    
        if i == j:
            return 0.0 
        else:
            if i < j:
                return self.inner_data[ i*(self.row_length_minus_one) - i - ((( i-1)*i) >> 1 ) + j-1]
            else:
                return self.inner_data[j*(self.row_length_minus_one) - j - ((( j-1)*j) >> 1 ) + i-1]
   
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __setitem__(self, tuple index, double item):
        cdef int real_pos,i,j #@DuplicatedSignature
        i = index[0]
        j = index[1]
        
        if i != j:
            if i < j:
                real_pos = condensedsubscript(self.row_length_minus_one,i,j-1)
            else:
                real_pos = condensedsubscript(self.row_length_minus_one,j,i-1)
            self.inner_data[real_pos] = item

    cpdef get_data(self):
        cdef numpy.npy_intp shape[1]
        shape[0] = <numpy.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = numpy.PyArray_SimpleNewFromData(1, shape, numpy.NPY_DOUBLE, self.inner_data)
        return ndarray

    def get_number_of_rows(self):
        return self.row_length