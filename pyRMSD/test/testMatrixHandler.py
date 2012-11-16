'''
Created on 20/09/2012

@author: victor
'''
import unittest
from pyRMSD.matrixHandler import MatrixHandler
from pyRMSD.condensedMatrix import CondensedMatrix 
import numpy
import os

class Test(unittest.TestCase):


    def test_write_and_load(self):
        mh = MatrixHandler(None)
        data = range(1000)
        matrix = CondensedMatrix(data)
        mh.distance_matrix = matrix
        mh.saveMatrix("matrix")
        
        mh2 = MatrixHandler(None)
        mh2.loadMatrix("matrix")
        recovered_data = mh2.distance_matrix.get_data()
        
        numpy.testing.assert_array_equal(mh.distance_matrix.get_data(), data)
        numpy.testing.assert_array_equal(recovered_data, data) 
        
        # Clean it!
        os.system("rm matrix.bin")
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_write_and_load']
    unittest.main()