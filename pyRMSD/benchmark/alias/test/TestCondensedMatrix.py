'''
Created on 30/01/2012

@author: victor
'''
import unittest
import scipy.spatial.distance as distance
import cStringIO
import random
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix, load_condensed_matrix, calc_number_of_rows,complete_to_condensed,zero_condensed
from pyproclust.matrix.completeMatrix import CompleteDistanceMatrix
import numpy as np

class Test(unittest.TestCase):

    def test_equal(self):
        cm1 = CondensedDistanceMatrix([1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9])
        cm2 = CondensedDistanceMatrix([1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9])
        cm3 = CondensedDistanceMatrix([6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4])
        cm4 = CondensedDistanceMatrix([6,7,8,9,0,1,2,3])
        
        self.assertEqual(cm1 == cm2, True)
        self.assertEqual(cm1 == cm3, False)
        self.assertEqual(cm1 == cm4, False)
        self.assertEqual(cm2 == cm3, False)
        self.assertEqual(cm2 == cm4, False)
        self.assertEqual(cm3 == cm4, False)

    def test_compare_condensed_matrixes(self):
        cm1 = CondensedDistanceMatrix([1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9])
        cm2 = CondensedDistanceMatrix([6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4])
        cm3 = CondensedDistanceMatrix([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
        cm4 = CondensedDistanceMatrix([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
        result_1 = cm1.compare_with(cm2) 
        result_2 = cm1.compare_with(cm3) 
        result_3 = cm3.compare_with(cm4,1.,2.)
        result_4 = cm3.compare_with(cm4,1.,1.)
        self.assertEqual(result_1, (5.0, 0.0))
        self.assertEqual(result_2, (3.8421052631578947, 2.6008734948643863))
        self.assertEqual(result_3, (0., 0.))
        self.assertEqual(result_4, (0.5, 0.))

    def test_get_number_of_rows(self):
        random.seed()
        for i in range(100): #@UnusedVariable
            rows = random.randint(1,1000)
            number_of_elements =  (rows *(rows-1)) / 2
            calculated_rows = calc_number_of_rows(number_of_elements)
            self.assertEqual(rows,calculated_rows)
    
    def test_normalize_condensed_matrix(self):
        condensed = CondensedDistanceMatrix([ 1., 4.5, 8.5, 7.2, 4.5, 7.8, 6.7, 3.6,2.2, 2.])
        expected = CondensedDistanceMatrix([0.0, 0.47, 1.0, 0.83, 0.47, 0.91, 0.76, 0.35, 0.16, 0.13])
        minmax = condensed.get_minimum_and_maximum()
        condensed.normalize(minmax[0], minmax[1])
        for i in range(len(condensed.get_data())):
            self.assertAlmostEqual(condensed.get_data()[i],expected.get_data()[i],2)
       
    def test_data_sharing(self):
        mylist = [ 1., 4.5, 8.5, 7.2, 4.5, 7.8, 6.7, 3.6,2.2, 2.]        
        myarray = np.array([ 1., 4.5, 8.5, 7.2, 4.5, 7.8, 6.7, 3.6,2.2, 2.])
        mylistaarray = np.array(mylist)
        condensed1 = CondensedDistanceMatrix(mylist)
        condensed2 = CondensedDistanceMatrix(myarray)
        condensed3 = CondensedDistanceMatrix(mylistaarray)
        
        mylist[5] = 0.
        self.assertEqual(False, mylist[5] == condensed1.get_data()[5])
        myarray[5] = 0.
        self.assertEqual(False, myarray[5] == condensed2.get_data()[5])
        mylistaarray[5] = 0.
        self.assertEqual(False, mylistaarray[5] == condensed3.get_data()[5])
        
        mycontents = condensed3.get_data()
        mycontents[5] = 0.
        self.assertEqual(True, mycontents[5] == condensed3.get_data()[5] and\
                         condensed3.get_data()[5] == 0.)

    def test_gen_condensed_matrix(self):
        obs = [(1,1),(2,1),(4,5),(7,7),(5,7)]
        ## distance matrix
        distance_matrix =  CompleteDistanceMatrix(distance.cdist(obs,obs))
        
        ## lower distance matrix (wo diagonal)
        expected_distance_condensed = CondensedDistanceMatrix(distance.pdist(obs))
        distance_condensed = complete_to_condensed(distance_matrix)
        self.assertEqual(True,distance_condensed ==  expected_distance_condensed)
        
    def test_validate_dimensions(self):
        condensed_matrix_1 = CondensedDistanceMatrix([ 1., 4.5, 8.5, 7.2, 4.5, 7.8, 6.7, 3.6,2.2, 2.])
        self.assertEqual(True,condensed_matrix_1._CondensedDistanceMatrix__validate_dimensions())
        condensed_matrix_2 = CondensedDistanceMatrix([ 1., 4.5, 8.5, 7.2, 4.5, 7.8, 6.7, 3.6])
        self.assertEqual(False,condensed_matrix_2._CondensedDistanceMatrix__validate_dimensions())
        
    def test_minmax_condensed(self):
        condensed_matrix = CondensedDistanceMatrix([ 1., 
                                                    4.5, 8.5, 
                                                    7.2, 4.5, 7.8, 
                                                    6.7, 3.6,2.2, 2.0])
        expected = (1,8.5)
        self.assertEqual(condensed_matrix.get_minimum_and_maximum(),expected)

    def test_save_condensed_matrix(self):
        # with final spaces!
        expected_matrix_string =  """1.0 4.5 7.2 6.7 
8.5 4.5 3.6 
7.8 2.2 
2.0 
"""
        condensed_matrix = CondensedDistanceMatrix([1.0, 4.5, 7.2, 6.7, 
                                                         8.5, 4.5, 3.6, 
                                                              7.8, 2.2, 
                                                                   2.0]) 
        output = cStringIO.StringIO()
        condensed_matrix.save(output)
        self.assertEqual(expected_matrix_string,output.getvalue())
        
    def test_load_condensed_matrix(self):
        matrix_string =  """1.0 
4.5 8.5 
7.2 4.5 7.8 
6.7 3.6 2.2 2.0
"""
        expected_matrix = CondensedDistanceMatrix([ 1., 4.5, 8.5, 7.2, 4.5, 7.8, 6.7, 3.6,2.2, 2.])
        input = cStringIO.StringIO(matrix_string)
        loaded_matrix =  load_condensed_matrix(input)
        for i in range(len(expected_matrix.get_data())):
            self.assertAlmostEqual(expected_matrix.get_data()[i],\
                                   loaded_matrix.get_data()[i],3)

    def test_item_access(self):
        condensed_matrix_1 = CondensedDistanceMatrix([1.0, 4.5,7.2, 
                                                        8.5, 4.5, 
                                                             7.8])
        condensed_matrix_2 = CondensedDistanceMatrix([.0]*6)
        
        complete_matrix = [[0.0, 1.0, 4.5, 7.2],
                           [1.0, 0.0, 8.5, 4.5],
                           [4.5, 8.5, 0.0, 7.8],
                           [7.2, 4.5, 7.8, 0.0]]
        
        row_len = condensed_matrix_1.row_length
        
        for i in range(row_len):
            for j in range(row_len):
                condensed_matrix_2[i,j] = complete_matrix[i][j]
        
        ## The access for a complete and a condensed matrix is exactly the same
        for i in range(row_len):
            for j in range(row_len):
                self.assertEquals(condensed_matrix_1[i,j],complete_matrix[i][j])
        
        ## And we can build a condensed matrix as a complete matrix
        self.assertItemsEqual(condensed_matrix_1.get_data(), condensed_matrix_2.get_data())
    
    def test_zero_condensed(self):
        row_len = 5
        zeroed_condensed = zero_condensed(row_len)
        self.assertEqual(row_len,zeroed_condensed.row_length)
        for i in range(row_len):
            for j in range(row_len):
                self.assertEquals(zeroed_condensed[i,j],0.)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()