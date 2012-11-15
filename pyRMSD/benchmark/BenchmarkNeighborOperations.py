'''
Created on 10/08/2012

@author: victor
'''
import time
import random
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyRMSD.cython.CythonCondensedMatrix import CythonCondensedMatrix
from pyRMSD.cython.CythonCondensedMatrix import CythonCondensedMatrixWithDiagonal
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproclust.algorithms.gromos.gromosAlgorithmTools import get_neighbors_for_node

if __name__ == '__main__':
    
    print "Creating data..."
    row_size = 10000
    matrix_elem_size = row_size*(row_size-1)/2
    contents = numpy.array(random.sample(xrange(matrix_elem_size+1),matrix_elem_size))
    float_contents = contents / float(max(contents))
    del contents
    matrix = CondensedMatrix(float_contents)
    matrix2 = CondensedDistanceMatrix(float_contents)
    remaining_nodes = range(row_size)

    print "======================================"
    print "'get_neighbors_for_node' benchmark"
    print "======================================"
    time_start = time.time()
    neighbors1 = matrix.get_neighbors_for_node(1,remaining_nodes, 0.5)
    time_end = time.time()
    print len(neighbors1)," neighbors found."
    print "Neighbor search for fast matrix took %.3fs"%(time_end-time_start)
    
    time_start = time.time()
    neighbors2 = matrix2.get_neighbors_for_node(1,remaining_nodes, 0.5)
    time_end = time.time()
    print len(neighbors2)," neighbors found."
    print "Neighbor search for slow matrix took %.3fs"%(time_end-time_start)
    print "======================================"
    print "'get_neighbors_for_node' validation"
    print "======================================"
    try:
        numpy.testing.assert_array_equal(neighbors1,neighbors2)
    except Exception,message:
        print message
        print "KO"
    
    print "================================================"
    print "'choose_node_with_higher_cardinality' benchmark"
    print "================================================"
    time_start = time.time()
    neighbors1 = matrix.choose_node_with_higher_cardinality(remaining_nodes, 0.5)
    time_end = time.time()
    print "Neighbor cardinality for fast matrix took %.3fs"%(time_end-time_start)
    
    time_start = time.time()
    neighbors2 = matrix2.choose_node_with_higher_cardinality(remaining_nodes, 0.5)
    time_end = time.time()
    print "Neighbor cardinality for slow matrix took %.3fs"%(time_end-time_start)
    print "================================================"
    print "'choose_node_with_higher_cardinality' validation"
    print "================================================"
    if(neighbors1 == neighbors2):
        print "OK"
    else:
        print "KO"