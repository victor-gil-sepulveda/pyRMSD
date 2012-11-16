'''
Created on 15/11/2012

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
import random 
import numpy
import scipy.stats

class Test(unittest.TestCase):
    
    def setUp(self):
        num_elems = 50*49/2;
        self.contents = random.sample(xrange(num_elems+1),num_elems)
        self.condensedMatrix = CondensedMatrix(self.contents)
        
    def test_mean(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateMean(), numpy.mean(self.contents))#, delta = 1)
    
    def test_variance(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateVariance(), numpy.var(self.contents))#, delta = 1)
    
    def test_skewness(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateSkewness(), scipy.stats.skew(self.contents))#d, delta = 1)
    
    def test_kurtosis(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateKurtosis(), scipy.stats.kurtosis(self.contents))#, delta = 1)
    
    def test_max(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateMax(), numpy.max(self.contents))#, delta = 1)

    def test_min(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateMin(), numpy.min(self.contents))#, delta = 1)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_mean']
    unittest.main()