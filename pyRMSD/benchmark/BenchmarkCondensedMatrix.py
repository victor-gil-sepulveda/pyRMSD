"""
Created on 07/08/2012

@author: victor
"""
from __future__ import print_function
import time
import random
from pyRMSD.condensedMatrix import CondensedMatrix
import pyRMSD.benchmark.alias.condensedMatrix as PythonCondensedMatrix
import pyRMSD.benchmark.alias.CythonCondensedMatrix as CythonCondensedMatrixes #@UnresolvedImport

if __name__ == '__main__':
    
    print("Creating data...")
    row_size = 10000
    matrix_elem_size = row_size*(row_size-1)/2
    contents = random.sample(xrange(matrix_elem_size+1),matrix_elem_size)
    matrixes = {}
    matrixes["Python ConndensedMatrix"] = PythonCondensedMatrix.CondensedMatrix(contents)
    matrixes["C CondensedMatrix"] = CondensedMatrix(contents)
    matrixes["CythonCondensedMatrix"] = CythonCondensedMatrixes.CythonCondensedMatrix(contents)
    matrixes["CythonCondensedMatrixWithDiagonal"] = CythonCondensedMatrixes.CythonCondensedMatrixWithDiagonal(contents)
    
    print("========================")
    print("Serial access Benchmark")
    print("========================")
    
    for matrix_type in matrixes:
        irange = range(row_size)
        jrange = range(row_size)
        
        time_start = time.time()
        for i in irange:
            for j in jrange:
                data = matrixes[matrix_type][i,j]
        time_end = time.time()
        print(matrix_type+" access benchmark took %.3f s"%(time_end-time_start))
       
    print("========================")
    print("Random Access Benchmark")
    print("========================")
    random_i = range(row_size)
    random.shuffle(random_i)
    random_j = range(row_size)
    random.shuffle(random_j)
    
    for matrix_type in matrixes:
        time_start = time.time()
        for i in random_i:
            for j in random_j:
                data =  matrixes[matrix_type][i,j]
        time_end = time.time()
        print("Random access benchmark took %.3f s for %s"%(time_end-time_start, matrix_type))
