import pyRMSD.calculators
from pyRMSD.availableCalculators import availableCalculators
import numpy

class RMSDCalculator(object):
    """
    Using the calculator will sometimes modify your input coordinates (shape and value). If the user wants to preserve the 
    input coordinates, the best way would be to save them to disk and recover them after the operation (i.e. numpy.save or 
    similar functions).
    """
    def __init__(self, calculatorType, 
                 fittingCoordsets, 
                 calculationCoordsets = None, 
                 symmetryGroups = []):
        """
        Class constructor.
        
        @param calculatorType: One of the calculators returned by 'availableCalculators()'. i.e. KABSCH_OMP_CALCULATOR
        
        @param fittingCoordsets: An array containing the used coordinates of each conformation. It has the following form:
            coordsets: [Conformation 1, Conformation 2, ..., Conformation N]
            Conformation: [Atom 1, Atom 2,..., Atom M]
            Atom: [x,y,z]
            
            This coordinates will be used for both structural superposition and RMSD calculation if the 'calculation 
            coordinates' parameter is not defined.
            
            The array type is constrained to be a numpy.array object with dtype = numpy.float64 (dtype will not be always double
            by default).
            
            Input coordinates are modified after each operation (usually centered and superposed into reference conformation).
        
        @param calculationCoordsets: An array containing the coordinates used to calculate the RMSD. Must have the same structure 
            than 'fittingCoordinates'.
        
        @param symmetryGroups: List of symmetry groups. Each symmetry group is a 2-Tuple of n-tuples defining interchangeable positions
            of current calculation coordinates (which will be the fitting coordinates, or the calculation coordinates if defined). 
            For instance, given a calculator with this fitting coordinates (without calculation coordinates):
                 a,  b,  c,  d,  e , f ;  for conf1
                 a1, b1, c3, d4, e5, f6 ; for conf2
                 a2, b2, c3, d4, e5, f6 ; for conf3
            
            and a symmetry group ((1,3),(2,4)). The final RMSDs would be the minimum RMSDs of:
                a,  b,  c,  d,  e , f ;  for conf1 (original)
                a,  c,  b,  e,  d , f ;  for conf1 (applying symmetry group)
            
            with:
                a1, b1, c3, d4, e5, f6 ; for conf2
                a2, b2, c3, d4, e5, f6 ; for conf3
                
            This is a low level description (specially because is structure agnostic) of the type of symmetries that can be found in 
            some ligands, i.e. in rotating benzene groups. It can also be used in symmetries of bigger selections though. 
            
        @author: vgil
        @date: 26/11/2012
        """
        if not calculatorType in availableCalculators():
            print "Calculator ", calculatorType, " is not an available calculator."
            raise ValueError
        else:
            self.fitting_coordinates = fittingCoordsets
            self.calculator_type = calculatorType
            self.number_of_conformations = self.fitting_coordinates.shape[0]
            self.number_of_fitting_atoms = self.fitting_coordinates.shape[1]
            
            self.calculation_coordinates = calculationCoordsets
            if self.calculation_coordinates is not None:
                self.calculation_coordinates = calculationCoordsets
                if self.number_of_conformations != self.calculation_coordinates.shape[0]:
                    print "Calculation coordinates must hold the same number of conformations than fitting coordinates."
                    raise ValueError
                self.number_of_calculation_atoms = self.calculation_coordinates.shape[1]
            
            # Default values for openMP and CUDA flags
            self.__threads_per_block = 32
            self.__blocks_per_grid = 8
            self.__number_of_threads = 8
    
    def pairwise(self, first_conformation_number, second_conformation_number, get_superposed_coordinates = False):
        """
        Calculates the rmsd of two conformations. As it must create a copy of the conformations, it does not modify the
        input coordinates. However input coordinates constantness is never guaranteed.
        
        @param first_conformation_number: Id of the reference conformation.
        
        @param second_conformation_number: Id of the conformation that will be superposed into the reference conformation.
        
        @param get_superposed_coordinates: If true, the function will also return the superposed coordinates.
    
        @return: The RMSD value or a tuple consisting of the RMSD of both conformations, the superposed fitting coordinates (copy), and the
        superposed calculation coordinates (copy, this last only if calculation coordinates were defined). Superposed fitting and 
        calculation resulting coordinates is a 2 conformation array.
        
        @author: vgil
        @date: 26/11/2012
        """
        
        # Working with fitting coordinates
        first_coords = self.fitting_coordinates[first_conformation_number]
        second_coords = self.fitting_coordinates[second_conformation_number]
        tmp_coordsets = numpy.copy(numpy.array([first_coords,second_coords]))
        
        # Then with calculation coordinates (if available)
        tmp_calculation_coordsets = None
        if  self.calculation_coordinates is not None:
            first_calculation_coords = self.calculation_coordinates[first_conformation_number]
            second_calculation_coords = self.calculation_coordinates[second_conformation_number]
            tmp_calculation_coordsets = numpy.array([first_calculation_coords,second_calculation_coords])
        
        rmsd = RMSDCalculator(self.calculator_type, tmp_coordsets, tmp_calculation_coordsets).oneVsFollowing(0)[0]
        
        if get_superposed_coordinates:
            if tmp_calculation_coordsets is None:
                return rmsd, tmp_coordsets
            else:
                return rmsd, tmp_coordsets, tmp_calculation_coordsets
        else:
            return rmsd
    
    def oneVsTheOthers(self, conformation_number, get_superposed_coordinates = False):
        """
        Calculates the RMSD between a reference conformation and all other conformations in the set.
        
        @param conformation_number: The id of the reference structure.
        
        @return: The array of RMSD values or a tuple consisting of the RMSD array, a superposed copy of fitting coordinates, and a 
        superposed copy of calculation coordinates (this last only if calculation coordinates were defined). 
        
        @author: vgil
        @date: 26/11/2012
        """
        # Rearrange fitting coordinates
        previous_coords = self.fitting_coordinates[:conformation_number]
        following_coords = self.fitting_coordinates[conformation_number+1:]
        rearranged_coords_list = [numpy.copy(self.fitting_coordinates[conformation_number])]
                
        for coords in previous_coords:
            rearranged_coords_list.append(numpy.copy(coords))
        
        for coords in following_coords:
            rearranged_coords_list.append(numpy.copy(coords))
        
        rearranged_coords = numpy.array(rearranged_coords_list)
        
        # Rearrange calculation coordinates
        rearranged_calculation_coords = None
        if self.calculation_coordinates is not None:
            previous_coords = self.calculation_coordinates[:conformation_number]
            following_coords = self.calculation_coordinates[conformation_number+1:]
            rearranged_calculation_coords_list = [self.calculation_coordinates[conformation_number]]
    
            for coords in previous_coords:
                rearranged_calculation_coords_list.append(coords)
            
            for coords in following_coords:
                rearranged_calculation_coords_list.append(coords)
            
            rearranged_calculation_coords = numpy.array(rearranged_calculation_coords_list)
        
        rmsd_array= RMSDCalculator(self.calculator_type, rearranged_coords, rearranged_calculation_coords).oneVsFollowing(0)
    
        if get_superposed_coordinates:
            if rearranged_calculation_coords is None:
                return rmsd_array, rearranged_coords
            else:
                return rmsd_array, rearranged_coords, rearranged_calculation_coords
        else:
            return rmsd_array
    
    def oneVsFollowing(self, conformation_number):
        """
        Calculates the RMSD between a reference conformation and all other conformations with an id greater than it.
        
        @param conformation_number: The id of the reference structure.
        
        @return: A numpy array of RMSD values.
        
        @author: vgil
        @date: 26/11/2012
        """
        rmsds = None
        
        # This kind of reshaping must do it locally (without copy)
        np_coords_fit = numpy.reshape(self.fitting_coordinates, self.number_of_conformations*self.number_of_fitting_atoms*3)
        
        if (self.calculation_coordinates is None):
            rmsds =  pyRMSD.calculators.oneVsFollowing(
                             availableCalculators()[self.calculator_type], 
                             np_coords_fit, 
                             self.number_of_fitting_atoms, 
                             numpy.array([]), 
                             0,
                             conformation_number, 
                             self.number_of_conformations,
                             self.__number_of_threads, 
                             self.__threads_per_block, 
                             self.__blocks_per_grid)
        else:
            np_coords_calc = numpy.reshape(self.calculation_coordinates, self.number_of_conformations*self.number_of_calculation_atoms*3)
            rmsds =  pyRMSD.calculators.oneVsFollowing(
                             availableCalculators()[self.calculator_type], 
                             np_coords_fit, 
                             self.number_of_fitting_atoms, 
                             np_coords_calc, 
                             self.number_of_calculation_atoms,
                             conformation_number, 
                             self.number_of_conformations,
                             self.__number_of_threads, 
                             self.__threads_per_block, 
                             self.__blocks_per_grid)
        return rmsds
        
    def pairwiseRMSDMatrix(self):
        """
        Calculates the pairwise RMSD matrix for all conformations in the coordinates set.
        
        @return: A numpy array with the upper triangle of the matrix, in row major format.
        
        @author: vgil
        @date: 26/11/2012
        """
        np_coords = numpy.reshape(self.fitting_coordinates,self.number_of_conformations*self.number_of_fitting_atoms*3)
        
        if self.calculation_coordinates is None:
            return pyRMSD.calculators.calculateRMSDCondensedMatrix(availableCalculators()[self.calculator_type], 
                                                                   np_coords, 
                                                                   self.number_of_fitting_atoms, 
                                                                   numpy.array([]), 0,
                                                                   self.number_of_conformations,
                                                                   self.__number_of_threads, 
                                                                   self.__threads_per_block, 
                                                                   self.__blocks_per_grid)
        else:
            np_coords_calc = numpy.reshape(self.calculation_coordinates,self.number_of_conformations*self.number_of_calculation_atoms*3)
            return pyRMSD.calculators.calculateRMSDCondensedMatrix(availableCalculators()[self.calculator_type], 
                                                                   np_coords, 
                                                                   self.number_of_fitting_atoms, 
                                                                   np_coords_calc, 
                                                                   self.number_of_calculation_atoms,
                                                                   self.number_of_conformations,
                                                                   self.__number_of_threads, 
                                                                   self.__threads_per_block, 
                                                                   self.__blocks_per_grid)
    
    def iterativeSuperposition(self):
        """
        Calculates an iterative superposition of a set of conformations. When using this function,
        the input coordinates are changed (it wouldn't have too much sense otherwise). 
        In this case calculation coordinates are not used as such, but as a secondary coordinates set
        that is rotate along with the primary fitting set (which allows to use different selections 
        to do the iterative fit and any other calculation).
        
        @author: vgil
        @date: 27/03/2013
        """
        np_coords = numpy.reshape(self.fitting_coordinates,self.number_of_conformations*self.number_of_fitting_atoms*3)
        if self.calculation_coordinates is None:
            pyRMSD.calculators.iterativeSuperposition(availableCalculators()[self.calculator_type], 
                               np_coords, self.number_of_fitting_atoms, 
                               numpy.array([]), 
                               0,
                               self.number_of_conformations,
                               self.__number_of_threads, 
                               self.__threads_per_block, 
                               self.__blocks_per_grid)
        else:
            np_coords_calc = numpy.reshape(self.calculation_coordinates,self.number_of_conformations*self.number_of_calculation_atoms*3)
            pyRMSD.calculators.iterativeSuperposition(availableCalculators()[self.calculator_type], 
                               np_coords, 
                               self.number_of_fitting_atoms, 
                               np_coords_calc, 
                               self.number_of_calculation_atoms,
                               self.number_of_conformations,
                               self.__number_of_threads, 
                               self.__threads_per_block, 
                               self.__blocks_per_grid)
    
    def setNumberOfOpenMPThreads(self, number_of_threads):
        """
        Sets the number of threads to be used by an OpenMP calculator.
        
        @param number_of_threads: The number of threads to be used by OpenMP runtime.
        
        @author: vgil
        @date: 26/11/2012
        """
        if ("OMP" in self.calculator_type):
            self.__number_of_threads = number_of_threads
        else:
            print "Cannot set any OpenMP related parameter using this calculator: ", self.calculator_type
            raise KeyError
    
    def setCUDAKernelThreadsPerBlock(self, number_of_threads, number_of_blocks):
        """
        Sets the number of threads per block and blocks per grid in CUDA calculators.
        
        @param number_of_threads: Number of threads per block to be used when launching CUDA Kernels.
        
        @param number_of_blocks: Number of blocks per grid to be used when launching CUDA Kernels.
        
        @author: vgil
        @date: 26/11/2012
        """
        if ("CUDA" in self.calculator_type):
            self.__threads_per_block = number_of_threads
            self.__blocks_per_grid = number_of_blocks
        else:
            print "Cannot set any CUDA related parameter using this calculator: ", self.calculator_type
            raise KeyError
