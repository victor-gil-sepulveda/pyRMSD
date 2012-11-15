'''
Created on 07/08/2012

@author: victor
'''
import time
import random
from pyRMSD.condensedMatrix import CondensedMatrix
from pyRMSD.cython.CythonCondensedMatrix import CythonCondensedMatrix
from pyRMSD.cython.CythonCondensedMatrix import CythonCondensedMatrixWithDiagonal

if __name__ == '__main__':
    
    print "Creating data..."
    row_size = 1000
    matrix_elem_size = row_size*(row_size-1)/2
    contents = random.sample(xrange(matrix_elem_size+1),matrix_elem_size)
    matrix = CythonCondensedMatrixWithDiagonal(contents)
    matrix2 = CondensedMatrix(contents)
    matrix3= CythonCondensedMatrix(contents)
    
    
    print "========================"
    print "Serial access Benchmark"
    print "========================"
    
    irange = range(row_size)
    jrange = range(row_size)
    
    time_start = time.time()
    for i in irange:
        for j in jrange:
            data = matrix[i,j]
    time_end = time.time()
    print "Serial access benchmark took %.3f s"%(time_end-time_start)
    
    time_start = time.time()
    for i in irange:
        for j in jrange:
            data = matrix2[i,j]
    time_end = time.time()
    print "Serial access benchmark took %.3f s"%(time_end-time_start)
    
    time_start = time.time()
    for i in irange:
        for j in jrange:
            data = matrix3[i,j]
    time_end = time.time()
    print "Serial access benchmark took %.3f s"%(time_end-time_start)
    
    print "========================"
    print "Random Access Benchmark"
    print "========================"
    random_i = range(row_size)
    random.shuffle(random_i)
    random_j = range(row_size)
    random.shuffle(random_j)
    
    time_start = time.time()
    for i in random_i:
        for j in random_j:
            data = matrix[i,j]
    time_end = time.time()
    print "Random access benchmark took %.3f s"%(time_end-time_start)
    
    time_start = time.time()
    for i in random_i:
        for j in random_j:
            data = matrix2[i,j]
    time_end = time.time()
    print "Random access benchmark took %.3f s"%(time_end-time_start)
    
    time_start = time.time()
    for i in random_i:
        for j in random_j:
            data = matrix3[i,j]
    time_end = time.time()
    print "Random access benchmark took %.3f s"%(time_end-time_start)
    
    print "========================"
    print "Validation"
    print "========================"
    for i in irange:
        for j in jrange:
            if not matrix3[i,j] == matrix2[i,j]:
                print "Error, matrices not equal when i= ",i, " and j=",j
                exit()
            if not matrix3[i,j] == matrix[i,j]:
                print "Error, matrices not equal when i= ",i, " and j=",j
                exit()
    print "OK"
    
    m2_data = matrix2.get_data()
    m3_data = matrix3.get_data()
    
    for i in range(matrix_elem_size):
        if(m3_data[i]!=m2_data[i]):
            print "Error, matrices not equal when i= ",i," ", m2_data[i]," ", m3_data[i]
            exit()
    print "OK"