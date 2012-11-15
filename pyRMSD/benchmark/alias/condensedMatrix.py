'''
Created on 12/03/2012

@author: victor
'''
import math
import numpy as np
from pyproclust.matrix.matrixCommon import generate_row
from pyproclust.algorithms.gromos.gromosAlgorithmTools import get_neighbors_for_node,\
    choose_node_with_higher_cardinality

def load_condensed_matrix(file_handler):
    """
    Loads a condensed matrix from an opened file handler.
    The file must have the upper triangle (without diagonals) of a distance matrix
    where each line is one row of this upper matrix.
    """
    m_contents= []
    for l in file_handler:
        row = []
        generate_row(l, row)
        m_contents.extend(row)
    return CondensedDistanceMatrix(m_contents)

def calc_number_of_rows(num_elements):
    """
    Calculates the dimension of the square matrix being represented by this one.
    """
    return int((1 + math.sqrt(1+8*num_elements))/2)

def complete_to_condensed(complete_matrix):
    """
    Generates a condensed distance matrix from a squared distance matrix as the 
    pdist function (in scipy.spatial.distance) does, which is the lower matrix part without 
    the diagonal.
    """
    column_len = len(complete_matrix.get_data())
    row_len = complete_matrix.get_row_dimension()
    matrix_data = complete_matrix.get_data()
    condensed_data = []
    for i in range(0,column_len):
        for j in range(i+1,row_len): 
            condensed_data.append(matrix_data[i][j]) 
    return CondensedDistanceMatrix(condensed_data) 

def zero_condensed(row_length):
    """
    Creates a condensed matrix representing a row_len x row_len matrix, filled with 0s.
    """
    num_elements = row_length * (row_length - 1) / 2
    return CondensedDistanceMatrix([0]*num_elements)

class CondensedDistanceMatrix(object):
    '''
    Stores a squared distance matrix in 'condensed form' (like the resulting
    matrix from scipy.spatial.distance.pdist .
    '''

    def __init__(self, contents=[],validate = False):
        '''
        Sets the contents and the number of elements and calculates the number the row length 
        that the original squared matrix had.
        It also checks that the matrix origin was a square matrix.
        '''
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
            
    def save(self, file_handler,low_precision = False):
        """
        Writes the condensed matrix into a file.
        """
        acc = self.row_length - 1
        i = 0
        while i < self.number_of_elements:
            for j in range(acc): #@UnusedVariable
                if low_precision:
                    file_handler.write("""%.4f """%(self.contents[i]))
                else:
                    file_handler.write(str(self.contents[i])+" ")
                i = i + 1
            file_handler.write("\n")
            acc = acc - 1

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
    
    def get_minimum_and_maximum(self):
        """
        Returns a tuple with the minimum and maximum values of a distance matrix
        in condensed form.
        """
        # In the condensed matrix there's no 0 diagonal
        # it should be impossible to have a value lesser than 0 but...
        return(min(self.contents),max(self.contents)) 
    
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
        except Exception,msg:
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
        except Exception,msg:
            print ("[Error in CondensedDistanceMatrix::__setitem__]",msg)
            exit(0)

    def __condensedsubscript(self,i,j):
        """
        Calculates the index of an element in a condensed matrix given a position of
        the real square matrix it represents.
        """
        
        return i*(self.row_length-1) - i - ( (( i-1)*i)/2) + j

    def normalize(self,minval,maxval):
        """
        Converts all the values of this matrix into a number in the range [0,1].
        It modifies the matrix!!
        """
        n_range = float(maxval-minval)
        self.contents = self.contents - minval
        self.contents = self.contents / n_range 
    
    def compare_with(self, this_other_one, mult_A = 1., mult_B = 1.):
        """
        Given two condensed matrixes (It's not going to check if they are or not) calculates
        the mean difference and the std deviation of them.
        """
        difference = abs(mult_A*self.get_data() - mult_B*this_other_one.get_data())
        return np.mean(difference), np.std(difference)
    
    def __eq__(self,other):
        """
        Equality operator.
        """
        if len(self.get_data()) != len(other.get_data()):
            return False
        else:
            for i in range(len(self.get_data())):
                if(self.get_data()[i] != other.get_data()[i]):
                        return False 
        return True

    def distance_distribution(self,element=None,granularity=50):
        """
        Calculates the distance distribution of this matrix with a granularity.
        """
        if element == None:
            return zip(np.histogram(a = self.get_data(),bins = granularity)[1],np.histogram(a = self.get_data(),bins = granularity)[0])
        else:
            distances = []
            for i in range(self.row_length):
                distances.append(self[element,i])
            return zip(np.histogram(a = distances,bins = granularity)[1],np.histogram(a = distances,bins = granularity)[0])
    
    
    def mean_distance_distribution(self):
        """
        Calculates the mean distance to all other elements for every element.
        """
        mean_distances = [0]*self.row_length
        for i in range(self.row_length-1):
            for j in range(i+1,self.row_length):
                mean_distances[i] += self[i,j]
                mean_distances[j] += self[i,j]
        
        for i in range(self.row_length):
            mean_distances[i] = mean_distances[i]/self.row_length
        
        mean_distances.sort()
        return mean_distances
        
    def percent_of_elements_within_cutoff_per_element(self,cutoff):
        print "TESTEAME PORFA (neighbors_per_element)"
        neighbors = [0]*self.row_length
        for i in range(self.row_length-1):
            for j in range(i+1,self.row_length):
                if self[i,j]<cutoff:
                    neighbors[i] += 1
                    neighbors[j] += 1
        for i in range(self.row_length):
            neighbors[i] = neighbors[i]/float(self.row_length)
        
        neighbors.sort()
        return neighbors
    
    def mean_distance_per_element(self,cutoff):
        print "TESTEAME PORFA (density_per_element)"
        densities = [0]*self.row_length
        neighbors = [1]*self.row_length
        for i in range(self.row_length-1):
            for j in range(i+1,self.row_length):
                if self[i,j]<cutoff:
                    densities[i] += self[i,j]
                    densities[j] += self[i,j]
                    neighbors[i] += 1
                    neighbors[j] += 1
                    
        for i in range(self.row_length):
            densities[i] =  densities[i]/neighbors[i]
            
        densities.sort()
        return densities
    
    def get_neighbors_for_node(self,node,nodes_left,cutoff):
        return get_neighbors_for_node(self,node,nodes_left,cutoff)
    
    def choose_node_with_higher_cardinality(self,nodes,cutoff):
        return choose_node_with_higher_cardinality(self,nodes,cutoff)
        
     
        