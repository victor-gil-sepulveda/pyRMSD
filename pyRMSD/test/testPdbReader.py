'''
Created on 29/11/2012

@author: victor
'''
import unittest
import numpy
from pyRMSD.utils.proteinReading import Reader

class testPdbReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mini_CA_coords = numpy.reshape(numpy.loadtxt("data/mini_ca_coords"), (3,9,3))
        cls.mini_all_coords = numpy.reshape(numpy.loadtxt("data/mini_all_coords"),(3,148,3))

    def test_read_one(self):
        reader = Reader().readThisFile('data/amber_mini.pdb')
        coordinates = reader.read()
        self.assertEqual(148, reader.numberOfAtoms)
        self.assertEqual(3, reader.numberOfFrames)
        numpy.testing.assert_almost_equal(testPdbReader.mini_all_coords, coordinates, 12)
    
    def test_read_multiple(self):
        reader = Reader().readThisFile('data/amber_mini.pdb').andThisOtherFile('data/amber_mini.pdb')
        coordinates = reader.read()
        self.assertEqual(148, reader.numberOfAtoms)
        self.assertEqual(6, reader.numberOfFrames)
        coord_shape = coordinates.shape
        self.assertEqual(coord_shape[0]*coord_shape[1]*coord_shape[2],148*6*3)
        numpy.testing.assert_almost_equal(numpy.reshape(self.mini_all_coords,3*148*3), numpy.reshape(coordinates,148*6*3)[0:148*3*3], 10)
        numpy.testing.assert_almost_equal(numpy.reshape(self.mini_all_coords,3*148*3), numpy.reshape(coordinates,148*6*3)[148*3*3:], 10)
     
    def test_read_one_only_CA(self):
        reader = Reader().readThisFile('data/amber_mini.pdb').gettingOnlyCAs()
        coordinates = reader.read() 
        self.assertEqual(9, reader.numberOfAtoms)
        self.assertEqual(3, reader.numberOfFrames)
        numpy.testing.assert_almost_equal(self.mini_CA_coords, coordinates, 10)     
     
    def test_read_multiple_only_CA(self):
        reader = Reader().readThisFile('data/amber_mini.pdb').andThisOtherFile('data/amber_mini.pdb').gettingOnlyCAs()
        coordinates = reader.read()
        self.assertEqual(9, reader.numberOfAtoms)
        self.assertEqual(6, reader.numberOfFrames)
        coord_shape = coordinates.shape
        self.assertEqual(coord_shape[0]*coord_shape[1]*coord_shape[2],9*6*3)
        numpy.testing.assert_almost_equal(numpy.reshape(self.mini_CA_coords,3*9*3), numpy.reshape(coordinates,6*9*3)[0:9*3*3], 10)
        numpy.testing.assert_almost_equal(numpy.reshape(self.mini_CA_coords,3*9*3), numpy.reshape(coordinates,6*9*3)[9*3*3:], 10)
         
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()