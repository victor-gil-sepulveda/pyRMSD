"""
Created on 16/11/2012

@author: victor
"""
from __future__ import print_function
import unittest
import pyRMSD.RMSDCalculator
import numpy

class TestRMSDCalculators(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        
        cls.expected_rmsd = {0:[ 0.60179184,  0.70814575,  0.88785042,  0.92862096,  0.69024252,  0.59267699,  0.66596155,  0.81180133,  0.87438831,  1.00465129],
                              1:[ 0.61473279,  0.82416178,  0.96955624,  0.71842781,  0.5359385,   0.68621908,  0.90540226,  0.83185205,  0.96145774],
                              2:[ 1.02156795,  1.16059055,  0.80778577,  0.72752425,  0.80478222,  0.98594799,  1.04869932,  1.01149253],
                              3:[ 0.69628994,  1.04059251,  0.77859792,  0.74962628,  0.73856698,  0.70444404,  0.92168545]}
        
        cls.expected_serial_matrix = [0.60179184,0.70814575,0.88785042,0.92862096,0.69024252,0.59267699,
                                       0.66596155,0.81180133,0.87438831,1.00465129,0.61473279,0.82416178,
                                       0.96955624,0.71842781,0.5359385, 0.68621908,0.90540226,0.83185205,
                                       0.96145774,1.02156795,1.16059055,0.80778577,0.72752425,0.80478222,
                                       0.98594799,1.04869932,1.01149253,0.69628994,1.04059251,0.77859792,
                                       0.74962628,0.73856698,0.70444404,0.92168545,1.08217543,0.86196576,
                                       0.89731473,0.96848922,0.84721509,1.13748551,0.64892912,0.87248355,
                                       1.00029474,1.01622641,1.10694473,0.68347196,0.83819283,0.7589582,
                                       0.93694602,0.76944618,0.82288799,0.91196003,0.75938856,0.68278426,
                                       0.76302383]
    def setUp(self):
#         print("In method", self._testMethodName)
        # Each test gets a fresh coordinates set, as they will modify the coordinates
        # QCP is specially sensitive to variations in input coordinates and results can vary
        self.coordsets_mini =  numpy.load("data/coordsets_mini.npy")
        self.coordsets =  numpy.load("data/coordsets.npy")
        self.number_of_conformations = self.coordsets.shape[0]
        self.number_of_atoms = self.coordsets.shape[1]
    
    def test_serial_omp_pairwise(self):
        """
        Calculates all matrix elements one by one with the pairwise operation.
        """
        expected_rmsd_data = [0.22677106513739653, 0.44598234794295144, 0.37817804816455303]
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_OMP_CALCULATOR",self.coordsets_mini)
        rmsd = []
        for i in range(len(self.coordsets_mini)):
            for j in range(i+1, len(self.coordsets_mini)):
                rmsd.append(calculator.pairwise(i, j))
        numpy.testing.assert_array_almost_equal(rmsd, expected_rmsd_data,8)
      
    @unittest.skipIf(not "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_cuda_pairwise(self):
        """
        Calculates all matrix elements one by one with the pairwise operation.
        """
        expected_rmsd_data = [0.22677106513739653, 0.44598234794295144, 0.37817804816455303]
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_CUDA_CALCULATOR",self.coordsets_mini)
        rmsd = []
        for i in range(len(self.coordsets_mini)):
            for j in range(i+1, len(self.coordsets_mini)):
                rmsd.append(calculator.pairwise(i, j))
        numpy.testing.assert_array_almost_equal(rmsd, expected_rmsd_data,4)
      
    def test_one_vs_others_serial_omp(self):
        """
        Calculates the reference vs the others with the OpenMP functions.
        """
        expected = [0.88785042, 0.82416178, 1.02156795, 0.69628994, 1.04059251, 0.77859792, 0.74962628, 0.73856698, 0.70444404, 0.92168545]
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_OMP_CALCULATOR", self.coordsets).oneVsTheOthers(3)
        numpy.testing.assert_array_almost_equal(rmsd, expected,8)
          
    def test_kabsch_serial(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("KABSCH_SERIAL_CALCULATOR",self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
      
    def test_kabsch_OpenMP(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("KABSCH_OMP_CALCULATOR", self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
       
    def test_serial(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", self.coordsets )
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
       
    def test_OpenMP(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_OMP_CALCULATOR", self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
               
    def test_theobald_serial(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_SERIAL_CALCULATOR", self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
               
    def test_theobald_OpenMP(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
       
    @unittest.skipIf(not "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_CUDA(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_CUDA_CALCULATOR", self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],4)
       
    @unittest.skipIf(not "QCP_CUDA_MEM_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_experimental_CUDA(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator( "QCP_CUDA_MEM_CALCULATOR", self.coordsets)
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],4)
               
    def test_serial_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", self.coordsets)
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
       
    def test_openmp_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_OMP_CALCULATOR", self.coordsets)
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
       
    def test_theobald_serial_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_SERIAL_CALCULATOR", self.coordsets)
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
       
    def test_theobald_OpenMP_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", self.coordsets)
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
       
    @unittest.skipIf(not "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_cuda_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_CUDA_CALCULATOR", self.coordsets)
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix, 4)
       
    @unittest.skipIf(not "QCP_CUDA_MEM_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_cuda_experimental_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_CUDA_MEM_CALCULATOR", self.coordsets)
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix, 4)
       
    def test_coordinates_change(self):
        """
        Tests if 
        """
        number_of_coordsets = 5;
        number_of_atoms = 3239;
        number_of_CAs = 224;
   
        not_aligned_CA = numpy.reshape(numpy.loadtxt("data/ligand_mini_CAs"), (number_of_coordsets,number_of_CAs,3))
           
        reference_copy = numpy.copy(not_aligned_CA[0])
        target_copy = numpy.copy(not_aligned_CA[4])
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_OMP_CALCULATOR", not_aligned_CA)
        rmsd = calculator.pairwise(0, 4, get_superposed_coordinates=False)
        self.assertAlmostEqual(rmsd, 0.3383320758562839, 12)
  
        # After the calculations, the target and reference coordinates have not changed, because 
        # a copy was performed ...
        numpy.testing.assert_almost_equal(not_aligned_CA[0],reference_copy,16)
        numpy.testing.assert_almost_equal(not_aligned_CA[4],target_copy,16)
      
      
    # Helper functions
    def read_coords_file(self, path):
        file_handler = open(path, "r")
        data = []
        for line in file_handler:
            if line.strip() != "":
                parts = line.split()
                data.append([float(x) for x in parts])
        file_handler.close()
        float_shape = data.pop(0)
        shape = (int(float_shape[0]),int(float_shape[1]),int(float_shape[2]))
        return numpy.array(data).reshape(shape)
      
    def read_rmsd_file(self, path):
        file_handler = open(path, "r")
        data = []
        for line in file_handler:
            if line.strip() != "":
                data.append(float(line.strip()))
        file_handler.close()
        return numpy.array(data)
    
    def test_fitting_symmetry(self):
        coords = self.read_coords_file( "../../src/calculators/test/data/Symmetry/OneVsAllFit/Native_0_ligand.pdb.aligned_BEN.coords")
        expected_rmsds = self.read_rmsd_file("../../src/calculators/test/data/Symmetry/OneVsAllFit/minimum.rmsd")
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          fittingCoordsets=coords,
                                                          fitSymmetryGroups=[
                                                                                [[0,2],[3,5]],
                                                                                [[7,8]]
                                                                             ])
        rmsds = calculator.oneVsFollowing(0)
        numpy.testing.assert_almost_equal(expected_rmsds, rmsds, 12)
    
    def test_calc_symmetry(self):
        
        CA_coords_native = self.read_coords_file( "../../src/calculators/test/data/Symmetry/Models/Natives/Native_3.CA.coords")
        CA_coords        = numpy.concatenate((CA_coords_native,
                                             self.read_coords_file( "../../src/calculators/test/data/Symmetry/Models/Trajectory/traj_testset.CA.coords")))
        
        lig_coords_native = self.read_coords_file( "../../src/calculators/test/data/Symmetry/Models/Natives/Native_3.ligand.coords")
        lig_coords        = numpy.concatenate((lig_coords_native,
                                              self.read_coords_file( "../../src/calculators/test/data/Symmetry/Models/Trajectory/traj_testset.ligand.coords")))
        
        
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          fittingCoordsets= numpy.array(CA_coords),
                                                          calculationCoordsets= numpy.array(lig_coords), 
                                                          calcSymmetryGroups=[
                                                                                [[0,2],[3,5]],
                                                                                [[7,8]]
                                                                             ])
        expected_rmsds = [0.705712355495, 0.683974525459, 0.685348934979, 0.687110565349, 0.689813612304, 0.694005214697, 0.694778845639, 0.697774242113]
        rmsds = calculator.oneVsFollowing(0)
        numpy.testing.assert_almost_equal(expected_rmsds, rmsds, 12)
        
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          fittingCoordsets= numpy.array(CA_coords),
                                                          calculationCoordsets= numpy.array(lig_coords), 
                                                          calcSymmetryGroups=[
                                                                                [[0,2],[3,5]],
                                                                                [[7,8]]
                                                                             ])
        
        rmsds = calculator.pairwiseRMSDMatrix()
        print(rmsds)
        
        
    # ------------------------------------------------------
    # Replication of some C tests (in  src/calculators/test)
    # ------------------------------------------------------
     
    def test_C_test_replication__superposition_with_fit_and_calc(self):
        # Loading initial conformations
        initial_CA_coords_for_fitting = self.read_coords_file(          "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords")
        initial_ligand_coords_for_calculating = self.read_coords_file(  "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords")
        # Loading final (superposed) conformations
        superposed_CA_coords = self.read_coords_file(       "../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_CA.coords")
        superposed_ligand_coords = self.read_coords_file(   "../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.coords")
        # Center
        centers = [ numpy.mean(conf,0) for conf in superposed_CA_coords]
        superposed_CA_coords = [ conf - centers[i] for i, conf in enumerate(superposed_CA_coords)]
        superposed_ligand_coords = [ conf - centers[i] for i, conf in enumerate(superposed_ligand_coords)]
          
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          initial_CA_coords_for_fitting, 
                                                          initial_ligand_coords_for_calculating)
        rmsds = calculator.oneVsFollowing(0)
          
        # Coordinates have been changed and have the expected_value
        numpy.testing.assert_almost_equal(superposed_CA_coords, initial_CA_coords_for_fitting, 12)
        numpy.testing.assert_almost_equal(superposed_ligand_coords, initial_ligand_coords_for_calculating,12)
          
        # RMSDS have the expected value
        expected_rmsds = self.read_rmsd_file("../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.rmsd")
        numpy.testing.assert_almost_equal(expected_rmsds[1:],rmsds,12)
      
    def test_C_test_replication__superposition_with_fit_and_calc_using_oneVsTheOthers(self):
        ### Results must be the same ones as in the test_C_test_replication__superposition_with_fit_and_calc,
        ### but coordinates change is handled differently
          
        # Loading initial conformations
        initial_CA_coords_for_fitting = self.read_coords_file(          "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords")
        initial_ligand_coords_for_calculating = self.read_coords_file(  "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords")
        initial_CA_coords_for_fitting_copy = self.read_coords_file(          "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords")
        initial_ligand_coords_for_calculating_copy = self.read_coords_file(  "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords")
          
        # Loading final (superposed) conformations
        superposed_CA_coords = self.read_coords_file(       "../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_CA.coords")
        superposed_ligand_coords = self.read_coords_file(   "../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.coords")
        # Center
        centers = [ numpy.mean(conf,0) for conf in superposed_CA_coords]
        superposed_CA_coords = [ conf - centers[i] for i, conf in enumerate(superposed_CA_coords)]
        superposed_ligand_coords = [ conf - centers[i] for i, conf in enumerate(superposed_ligand_coords)]
          
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          initial_CA_coords_for_fitting, 
                                                          initial_ligand_coords_for_calculating)
          
        # We expect that input coordinates are not changed, and to receive a copy of superposed coordinates.
        rmsds, superposed_fit, superposed_calc = calculator.oneVsTheOthers(0,get_superposed_coordinates = True)
          
        #Initial coordinates were not changed
        numpy.testing.assert_almost_equal(initial_CA_coords_for_fitting_copy, initial_CA_coords_for_fitting, 16)
        numpy.testing.assert_almost_equal(initial_ligand_coords_for_calculating_copy, initial_ligand_coords_for_calculating, 16)
                                            
        # Coordinates have been changed and have the expected_value
        numpy.testing.assert_almost_equal(superposed_CA_coords, superposed_fit, 12)
        numpy.testing.assert_almost_equal(superposed_ligand_coords, superposed_calc,12)
          
        # RMSDS have the expected value
        expected_rmsds = self.read_rmsd_file("../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.rmsd")
        numpy.testing.assert_almost_equal(expected_rmsds[1:],rmsds,12)
          
    def test_C_test_replication__superposition_with_fit_and_calc_using_pairwise(self):
        ### Results must be the same ones as in the test_C_test_replication__superposition_with_fit_and_calc,
        ### but again coordinates change
          
        # Loading initial conformations
        initial_CA_coords_for_fitting = self.read_coords_file(          "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords")
        initial_ligand_coords_for_calculating = self.read_coords_file(  "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords")
        initial_CA_coords_for_fitting_copy = self.read_coords_file(          "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords")
        initial_ligand_coords_for_calculating_copy = self.read_coords_file(  "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords")
          
        # Loading final (superposed) conformations
        superposed_CA_coords = self.read_coords_file(       "../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_CA.coords")
        superposed_ligand_coords = self.read_coords_file(   "../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.coords")
        # Center
        centers = [ numpy.mean(conf,0) for conf in superposed_CA_coords]
        superposed_CA_coords = [ conf - centers[i] for i, conf in enumerate(superposed_CA_coords)]
        superposed_ligand_coords = [ conf - centers[i] for i, conf in enumerate(superposed_ligand_coords)]
          
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          initial_CA_coords_for_fitting, 
                                                          initial_ligand_coords_for_calculating)
          
        # We expect that input coordinates are not changed, and to receive a copy of superposed coordinates.
        rmsd, superposed_fit, superposed_calc = calculator.pairwise(0, 1, get_superposed_coordinates = True)
          
        #Initial coordinates were not changed
        numpy.testing.assert_almost_equal(initial_CA_coords_for_fitting_copy, initial_CA_coords_for_fitting, 16)
        numpy.testing.assert_almost_equal(initial_ligand_coords_for_calculating_copy, initial_ligand_coords_for_calculating, 16)
                                            
        # Coordinates have been changed and have the expected_value
        numpy.testing.assert_almost_equal(superposed_CA_coords[0:2], superposed_fit, 12)
        numpy.testing.assert_almost_equal(superposed_ligand_coords[0:2], superposed_calc,12)
          
        # RMSDS have the expected value
        expected_rmsds = self.read_rmsd_file("../../src/calculators/test/data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.rmsd")
        numpy.testing.assert_almost_equal(expected_rmsds[1], rmsd, 12)
  
    def test_C_test_replication__iterative_superposition_with_fit_and_calc_rotation(self):
        # Loading initial conformations
        initial_fitting_coords = self.read_coords_file( "../../src/calculators/test/data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.initial_all.coords")
        initial_rotating_coords = self.read_coords_file( "../../src/calculators/test/data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.initial_BEN.coords" )
          
        final_fitting_coords = self.read_coords_file("../../src/calculators/test/data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.iterposed_all.coords")
        final_rotating_coords = self.read_coords_file("../../src/calculators/test/data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.iterposed_BEN.coords")
        # Center
        centers = [ numpy.mean(conf,0) for conf in final_fitting_coords]
        final_fitting_coords = [ conf - centers[i] for i, conf in enumerate(final_fitting_coords)]
        final_rotating_coords = [ conf - centers[i] for i, conf in enumerate(final_rotating_coords)]
          
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          initial_fitting_coords, 
                                                          initial_rotating_coords)
          
        calculator.iterativeSuperposition()
          
        # Coordinates MUST be modified
        numpy.testing.assert_almost_equal(final_fitting_coords, initial_fitting_coords, 9) # Decrease of precision also noted in C version
        numpy.testing.assert_almost_equal(final_rotating_coords, initial_rotating_coords,9)
      
    def test_C_test_replication__matrix_with_fit_and_calculation_coordinates(self):
        # Loading initial conformations
        initial_fitting_coords = self.read_coords_file( "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords")
        initial_rotating_coords = self.read_coords_file( "../../src/calculators/test/data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords")
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", 
                                                          initial_fitting_coords, 
                                                          initial_rotating_coords)
          
        expected_rmsds = self.read_rmsd_file("../../src/calculators/test/data/Matrix_Fit_CA_Calc_BEN/prot_plus_ligand_offset_very_different.CA.rmsd_matrix")
        numpy.testing.assert_almost_equal(expected_rmsds, calculator.pairwiseRMSDMatrix(), 12)
      
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
