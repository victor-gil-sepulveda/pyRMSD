'''
Created on 07/08/2012

@author: victor
'''
import unittest

from pyproct.matrix.cython.CythonCondensedMatrix import CythonCondensedMatrix #@UnresolvedImport
from pyproct.matrix.condensedMatrix import CondensedDistanceMatrix

class Test(unittest.TestCase):
    

    def test_item_get(self):
        condensed_matrix_1 = CondensedDistanceMatrix([1.0, 4.5,7.2, 
                                                        8.5, 4.5, 
                                                             7.8])
        condensed_matrix_2 = CythonCondensedMatrix([.0]*6)
        
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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_item_get']
    unittest.main()