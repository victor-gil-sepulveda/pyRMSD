'''
Created on 15/08/2012

@author: victor
'''
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix

class Test(unittest.TestCase):


    def test_Laplacian(self):
        
        expected_Laplacian = numpy.array([[ 1., -0.05154639, -0.23195877, -0.37113402, -0.34536082],
                                          [-0.05681818,  1.        , -0.48295455, -0.25568182, -0.20454545],
                                          [-0.19565217, -0.36956521,  1.        , -0.33913044, -0.09565217],
                                          [-0.33488371, -0.20930233, -0.36279071,  1.        , -0.09302326],
                                          [-0.46206896, -0.24827586, -0.15172414, -0.13793104,  1.        ]])
        matrix = CondensedMatrix([1.0, 4.5, 7.2, 6.7, 
                                       8.5, 4.5, 3.6, 
                                            7.8, 2.2, 
                                                 2.0]) 
        laplacian =  matrix.calculate_rw_laplacian()
        
        numpy.testing.assert_array_almost_equal(expected_Laplacian, laplacian, 5)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()