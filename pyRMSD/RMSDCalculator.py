import numpy.linalg
import pyRMSD.calculators
from pyRMSD.utils.proteinReading import flattenCoords
from pyRMSD.availableCalculators import availableCalculators
import numpy

class RMSDCalculator(object):
    def __init__(self, coordsets, calculatorType):
        """
        Class constructor
        
        @param coordsets: An array containing the used coordinates of each conformation. It has the following form:
            coordsets: [Conformation 1, Conformation 2, ..., Conformation N]  
            Conformation: [Atom 1, Atom 2,..., Atom M]  
            Atom: [x,y,z]
            
            This coordinates will be used for both structural superposition and RMSD calculation if the 'calculation 
            coordinates' parameter is not defined.
        
        @param calculatorType: One of the calculators returned by 'availableCalculators()'. i.e. OMP_CALCULATOR
         
        @author: vgil
        @date: 26/11/2012
        """
        self.calculation_coordinates = None

        if not calculatorType in availableCalculators():
            print "Calculator ", calculatorType, " is not an available calculator."
            raise KeyError
        else:
            self.fitting_coordinates = coordsets
            self.calculatorType = calculatorType
            self.number_of_conformations = len(coordsets)
            self.number_of_fitting_atoms = len(coordsets[0])
            self.__threads_per_block = 32
            self.__blocks_per_grid = 8
            self.__number_of_threads = 8
    
    def setCalculationCoordinates(self, coordsets):
        """
        Defines the coordinates used for RMSD calculation (will be 'fitting_coordinates' if this function is not used).
        
        @param coordsets:
        
        @author: vgil
        @date: 27/03/2013
        """
        self.calculation_coordinates = coordsets
        self.number_of_calculation_atoms = len(coordsets[0])
    
    def __getInternalCalculator(self, fitting_coordinates, calculation_coordinates):
        """
        Encapsulates the creation of an RMSD calculator for internal usage.
        
        @return: An RMSD calculator.
        
        @author: vgil
        @date: 27/03/2013
        """
        calculator = RMSDCalculator(fitting_coordinates, self.calculatorType)
        if (calculation_coordinates is not None):
            calculator.setCalculationCoordinates(calculation_coordinates)
        return calculator
    
    def pairwise(self, first_conformation_number, second_conformation_number):
        """
        Calculates the rmsd of two conformations.
        
        @param first_conformation_number: Id of the reference conformation.
        
        @param second_conformation_number: Id of the conformation that will be superposed into the reference conformation.
        
        @return: The RMSD value.
        
        @author: vgil
        @date: 26/11/2012
        """
        
        # Working with fitting coordinates
        first_coords = self.fitting_coordinates[first_conformation_number]
        second_coords = self.fitting_coordinates[second_conformation_number]
        tmp_coordsets = numpy.array([first_coords,second_coords])
        
        # Then with calculation coordinates (if available)
        tmp_calculation_coordsets = None
        if ( self.calculation_coordinates is not None):
            first_calculation_coords = self.calculation_coordinates[first_conformation_number]
            second_calculation_coords = self.calculation_coordinates[second_conformation_number]
            tmp_calculation_coordsets = numpy.array([first_calculation_coords,second_calculation_coords])
        
        return self.__getInternalCalculator(tmp_coordsets, tmp_calculation_coordsets).oneVsFollowing(0)[0]
    
    def oneVsTheOthers(self, conformation_number):
        """
        Calculates the RMSD between a reference conformation and all other conformations in the set.
        
        @param conformation_number: The id of the reference structure.
        
        @return: A numpy array of RMSD values.
        
        @author: vgil
        @date: 26/11/2012
        """
        # Rearrange fitting coordinates
        previous_coords = self.fitting_coordinates[:conformation_number]
        following_coords = self.fitting_coordinates[conformation_number+1:]
        rearranged_coords_list = [self.fitting_coordinates[conformation_number]]
                
        for coords in previous_coords:
            rearranged_coords_list.append(coords)
        
        for coords in following_coords:
            rearranged_coords_list.append(coords)
        
        rearranged_coords = numpy.array(rearranged_coords_list)
        
        # Rearrange calculation coordinates
        rearranged_calculation_coords = None
        if ( not self.calculation_coordinates is None):
            previous_coords = self.calculation_coordinates[:conformation_number]
            following_coords = self.calculation_coordinates[conformation_number+1:]
            rearranged_calculation_coords_list = [self.calculation_coordinates[conformation_number]]
    
            for coords in previous_coords:
                rearranged_coords_list.append(coords)
            
            for coords in following_coords:
                rearranged_coords_list.append(coords)
            
            rearranged_calculation_coords = numpy.array(rearranged_calculation_coords_list)
        
        return self.__getInternalCalculator(rearranged_coords, rearranged_calculation_coords).oneVsFollowing(0)
    
    def oneVsFollowing(self, conformation_number):
        """
        Calculates the RMSD between a reference conformation and all other conformations with an id 
        greater than it.
        
        @param conformation_number: The id of the reference structure.
        
        @return: A numpy array of RMSD values.
        
        @author: vgil
        @date: 26/11/2012
        """     
        np_coords_fit = flattenCoords(self.fitting_coordinates)
        if (self.calculation_coordinates is None):
            return pyRMSD.calculators.oneVsFollowing(
                                                     availableCalculators()[self.calculatorType], 
                                                     np_coords_fit, self.number_of_fitting_atoms, 
                                                     numpy.array([]), 0,
                                                     conformation_number, self.number_of_conformations,
                                                     self.__number_of_threads, self.__threads_per_block, self.__blocks_per_grid)
        else:
            np_coords_calc = flattenCoords(self.calculation_coordinates)
            return pyRMSD.calculators.oneVsFollowing(
                                                     availableCalculators()[self.calculatorType], 
                                                     np_coords_fit, self.number_of_fitting_atoms, 
                                                     np_coords_calc, self.number_of_calculation_atoms,
                                                     conformation_number, self.number_of_conformations,
                                                     self.__number_of_threads, self.__threads_per_block, self.__blocks_per_grid)
        
    def pairwiseRMSDMatrix(self):
        """
        Calculates the pairwise RMSD matrix for all conformations in the coordinates set.
        
        @return: A numpy array with the upper triangle of the matrix, in row major format.
        
        @author: vgil
        @date: 26/11/2012
        """
        np_coords = flattenCoords(self.fitting_coordinates)
        if (self.calculation_coordinates is None):
            return pyRMSD.calculators.calculateRMSDCondensedMatrix(
                                                                   availableCalculators()[self.calculatorType], 
                                                                   np_coords, self.number_of_fitting_atoms, 
                                                                   numpy.array([]), 0,
                                                                   self.number_of_conformations,
                                                                   self.__number_of_threads, self.__blocks_per_grid, self.__blocks_per_grid)
        else:
            np_coords_calc = flattenCoords(self.calculation_coordinates)
            return pyRMSD.calculators.calculateRMSDCondensedMatrix(
                                                                   availableCalculators()[self.calculatorType], 
                                                                   np_coords, self.number_of_fitting_atoms, 
                                                                   np_coords_calc, self.number_of_calculation_atoms,
                                                                   self.number_of_conformations,
                                                                   self.__number_of_threads, self.__blocks_per_grid, self.__blocks_per_grid)
            
    def setNumberOfOpenMPThreads(self, number_of_threads):
        """
        Sets the number of threads to be used by an OpenMP calculator.
        
        @param number_of_threads: The number of threads to be used by OpenMP runtime.
        
        @author: vgil
        @date: 26/11/2012
        """
        if ("OMP" in self.calculatorType):
            self.__number_of_threads = number_of_threads
        else:
            print "Cannot set any OpenMP related parameter using this calculator: ", self.calculatorType
            raise KeyError
    
    def setCUDAKernelThreadsPerBlock(self, number_of_threads, number_of_blocks):
        """
        Sets the number of threads per block and blocks per grid in CUDA calculators.
        
        @param number_of_threads: Number of threads per block to be used when launching CUDA Kernels.
        
        @param number_of_blocks: Number of blocks per grid to be used when launching CUDA Kernels.
        
        @author: vgil
        @date: 26/11/2012
        """
        if ("CUDA" in self.calculatorType):
            self.__threads_per_block = number_of_threads
            self.__blocks_per_grid = number_of_blocks
        else:
            print "Cannot set any CUDA related parameter using this calculator: ", self.calculatorType
            raise KeyError
