'''
Created on 15/02/2012

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix

class TestMatrixNeighbours(unittest.TestCase):

    def test_get_neighbors_for_node(self):
        condensed_matrix = CondensedMatrix( [0.2, 1.0, 0.3,  .0,
                                                  0.5, 0.6, 0.7,
                                                       0.9, 0.8,
                                                            0.4])
        
        self.assertItemsEqual([0,2],condensed_matrix.get_neighbors_for_node(1,range(5),0.5))

    
    def test_choose_node_with_higher_cardinality(self):
        condensed_matrix = CondensedMatrix([1., 4.5, 8.5, 7.2, 
                                                4.5, 7.8, 6.7, 
                                                     3.6, 2.2, 
                                                          2.0])
        nodes = range(5)
        self.assertEqual( condensed_matrix.choose_node_with_higher_cardinality( nodes, 4.),\
                          (2,2))
        
        nodes = range(1,5)
        self.assertEqual( condensed_matrix.choose_node_with_higher_cardinality( nodes, 4.),\
                          (2,2))
        
        nodes = [2,0,3]
        self.assertEqual( condensed_matrix.choose_node_with_higher_cardinality( nodes, 4.),\
                         (2,1))
        
        nodes = range(5)
        self.assertEqual( condensed_matrix.choose_node_with_higher_cardinality( nodes, 7.),\
                          (2,4))
         
        nodes = [2,0,3]
        self.assertEqual( condensed_matrix.choose_node_with_higher_cardinality(  nodes, 7.),\
                         (2,2))
        
    def test_get_cluster_for_node(self):
        condensed_matrix = CondensedMatrix([1., 4.5, 8.5, 7.2, 
                                                4.5, 7.8, 6.7, 
                                                     3.6, 2.2, 
                                                          2.0])
        row_len = 5
        nodes_left = range(row_len)
        
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 0, nodes_left, 4.)),1)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 1, nodes_left, 4.)),1)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 2, nodes_left, 4.)),2)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 3, nodes_left, 4.)),2)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 4, nodes_left, 4.)),2)
        
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 0, nodes_left, 7.)),2)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 1, nodes_left, 7.)),3)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 2, nodes_left, 7.)),4)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 3, nodes_left, 7.)),2)
        self.assertEqual(len(condensed_matrix.get_neighbors_for_node( 4, nodes_left, 7.)),3)
        
        neig_4_3 = [2,4]
        neig = condensed_matrix.get_neighbors_for_node( 3, nodes_left, 4.)
        for i in range(len(neig_4_3)):
            self.assertEqual(neig[i],neig_4_3[i])
        
        neig_7_2 = [0,1,3,4]
        neig = condensed_matrix.get_neighbors_for_node( 2, nodes_left, 7.)
        for i in range(len(neig_7_2)):
            self.assertEqual(neig[i],neig_7_2[i])
            
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()