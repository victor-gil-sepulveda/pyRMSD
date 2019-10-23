"""
Created on 12/03/2012

@author: victor
"""
import math
import numpy as np
import numpy
from pyRMSD.benchmark.alias.neighbourOps import get_neighbors_for_node,\
    choose_node_with_higher_cardinality


def calc_number_of_rows(num_elements):
    """
    Calculates the dimension of the square matrix being represented by this one.
    """
    return int((1 + math.sqrt(1+8*num_elements))/2)

class CondensedMatrix(object):
    """
    Stores a squared distance matrix in 'condensed form' (like the resulting
    matrix from scipy.spatial.distance.pdist .
    """

    def __init__(self, contents=[],validate = False):
        """
        Sets the contents and the number of elements and calculates the number the row length
        that the original squared matrix had.
        It also checks that the matrix origin was a square matrix.
        """
        self.contents = np.array(contents)
        self.number_of_elements = len(self.contents)
        self.row_length =  calc_number_of_rows(self.number_of_elements)
        # Check if the condensed matrix array has the expected number of elements.
        # Abort if not
        if validate and not self.__validate_dimensions():
            print "[Error in CondensedDistanceMatrix::__init__]: this condensed matrix is not a square matrix."
            exit(0)

    def __validate_dimensions(self):
        """
        Checks if the condensed matrix array has the expected number of elements
        """
        expected_num_elements =  self.row_length * (self.row_length - 1) / 2
        return expected_num_elements == self.number_of_elements

    def get_data(self):
        """
        Returns the real contents array of the condensed matrix.
        """
        return self.contents

    def get_number_of_rows(self):
        """
        Calculates the number of rows of the original squared matrix from a condensed
        matrix' number of elements
        """
        return calc_number_of_rows(self.number_of_elements)


    def __getitem__(self, index):
        """
        Calculates the index of an element in a condensed matrix given a position of
        the real square matrix it represents.
        """
        try:
            if index[0] == index[1]:
                return 0.0
            else:
                if index[0] < index[1]:
                    real_pos = self.__condensedsubscript(index[0],index[1]-1)
                else:
                    real_pos = self.__condensedsubscript(index[1],index[0]-1)
                return self.contents[real_pos]
        except Exception as msg:
            print ("[Error in CondensedDistanceMatrix::__getitem__]",msg)
            exit(0)

    def __setitem__(self, index, item):
        try:
            if index[0] == index[1]:
                return
            else:
                if index[0] < index[1]:
                    real_pos = self.__condensedsubscript(index[0],index[1]-1)
                else:
                    real_pos = self.__condensedsubscript(index[1],index[0]-1)
            self.contents[real_pos] = item
        except Exception as msg:
            print ("[Error in CondensedDistanceMatrix::__setitem__]",msg)
            exit(0)

    def __condensedsubscript(self,i,j):
        """
        Calculates the index of an element in a condensed matrix given a position of
        the real square matrix it represents.
        """

        return i*(self.row_length-1) - i - ( (( i-1)*i)/2) + j

    def get_neighbors_for_node(self,node,nodes_left,cutoff):
        return get_neighbors_for_node(self,node,nodes_left,cutoff)

    def choose_node_with_higher_cardinality(self,nodes,cutoff):
        return choose_node_with_higher_cardinality(self,nodes,cutoff)


