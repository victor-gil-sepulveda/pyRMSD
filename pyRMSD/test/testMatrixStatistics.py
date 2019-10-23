"""
Created on 15/11/2012

@author: victor
"""
from __future__ import division
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
import random 
import numpy
import sys

if sys.version_info[0] >= 3:
    xrange = range

class testMatrixStatistics(unittest.TestCase):
    
    def setUp(self):
        random.seed(12345)
        num_elems = 50*49//2;
        self.contents = random.sample(xrange(num_elems+1),num_elems)
        self.condensedMatrix = CondensedMatrix(self.contents)
        
    def test_mean(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateMean(), numpy.mean(self.contents))#, delta = 1)
    
    def test_variance(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateVariance(), numpy.var(self.contents))#, delta = 1)
    
    def test_skewness(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateSkewness(), -0.000733274016748)#scipy.stats.skew(self.contents))#d, delta = 1)
    
    def test_kurtosis(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateKurtosis(), -1.20119577613)#scipy.stats.kurtosis(self.contents))#, delta = 1)
    
    def test_max(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateMax(), numpy.max(self.contents))#, delta = 1)

    def test_min(self):
        self.assertAlmostEquals(self.condensedMatrix.calculateMin(), numpy.min(self.contents))#, delta = 1)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_mean']
    unittest.main()
