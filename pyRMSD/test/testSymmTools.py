'''
Created on 29/07/2013

@author: victor
'''
import unittest
from symmTools import  symm_group_permutator, swap_atoms
import numpy


class Test(unittest.TestCase):


    def test_symm_group_permutations(self):
        permutations = []
        symm_group_permutator([ ((5,7),(4,6)), ((3,),(4,))],[], permutations)
        self.assertItemsEqual([[], [((3,), (4,))], [((5, 7), (4, 6))], [((5, 7), (4, 6)), ((3,), (4,))]], permutations)
    
    def test_swap_atoms(self):
        coordsets = numpy.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
        swap_atoms(coordsets, 0, 2)
        numpy.testing.assert_array_equal([[7, 8, 9],[4, 5, 6],[1, 2, 3],[10, 11, 12]], coordsets)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_symm_group_permutations']
    unittest.main()