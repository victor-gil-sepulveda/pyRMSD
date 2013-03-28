'''
Created on 16/11/2012

@author: victor
'''
import unittest
import pyRMSD.utils.proteinReading
import pyRMSD.RMSDCalculator
import numpy
from prody import *

class TestRMSDCalculators(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        
        self.expected_rmsd = {0:[ 0.60179184,  0.70814575,  0.88785042,  0.92862096,  0.69024252,  0.59267699,  0.66596155,  0.81180133,  0.87438831,  1.00465129],
                              1:[ 0.61473279,  0.82416178,  0.96955624,  0.71842781,  0.5359385,   0.68621908,  0.90540226,  0.83185205,  0.96145774],
                              2:[ 1.02156795,  1.16059055,  0.80778577,  0.72752425,  0.80478222,  0.98594799,  1.04869932,  1.01149253],
                              3:[ 0.69628994,  1.04059251,  0.77859792,  0.74962628,  0.73856698,  0.70444404,  0.92168545]}
        
        self.expected_serial_matrix = [0.60179184,0.70814575,0.88785042,0.92862096,0.69024252,0.59267699,
                                       0.66596155,0.81180133,0.87438831,1.00465129,0.61473279,0.82416178,
                                       0.96955624,0.71842781,0.5359385, 0.68621908,0.90540226,0.83185205,
                                       0.96145774,1.02156795,1.16059055,0.80778577,0.72752425,0.80478222,
                                       0.98594799,1.04869932,1.01149253,0.69628994,1.04059251,0.77859792,
                                       0.74962628,0.73856698,0.70444404,0.92168545,1.08217543,0.86196576,
                                       0.89731473,0.96848922,0.84721509,1.13748551,0.64892912,0.87248355,
                                       1.00029474,1.01622641,1.10694473,0.68347196,0.83819283,0.7589582,
                                       0.93694602,0.76944618,0.82288799,0.91196003,0.75938856,0.68278426,
                                       0.76302383]
        
        reader = pyRMSD.utils.proteinReading.Reader("PRODY_READER")
        reader.readThisFile('data/amber_mini.pdb').gettingOnlyCAs()
        self.coordsets_mini =  reader.read()
        
        reader = pyRMSD.utils.proteinReading.Reader("PRODY_READER")
        reader.readThisFile('data/amber_short.pdb').gettingOnlyCAs()
        self.coordsets =  reader.read()
        self.number_of_conformations = reader.numberOfFrames
        self.number_of_atoms = reader.numberOfAtoms
    
    def test_serial_omp_pairwise(self):
        """
        Calculates all matrix elements one by one with the pairwise operation.
        """
        expected_rmsd_data = [0.22677106513739653, 0.44598234794295144, 0.37817804816455303]
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets_mini, "QTRFIT_OMP_CALCULATOR")
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
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets_mini, "QCP_CUDA_CALCULATOR")
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
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QTRFIT_OMP_CALCULATOR").oneVsTheOthers(3)
        numpy.testing.assert_array_almost_equal(rmsd, expected,8)
        
    def test_kabsch_serial(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "KABSCH_SERIAL_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
    
    def test_kabsch_OpenMP(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "KABSCH_OMP_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
    
    def test_serial(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QTRFIT_SERIAL_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
    
    def test_OpenMP(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QTRFIT_OMP_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
            
    def test_theobald_serial(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QCP_SERIAL_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
            
    def test_theobald_OpenMP(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QCP_OMP_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],8)
    
    @unittest.skipIf(not "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_CUDA(self):
        """
        Calculates the whole pairwise matrix by calculating each of the rows of the matrix.
        """
        for conf_num in self.expected_rmsd:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QCP_CUDA_CALCULATOR")
            rmsd = calculator.oneVsFollowing(conf_num)
            numpy.testing.assert_array_almost_equal(rmsd, self.expected_rmsd[conf_num],4)
            
    def test_serial_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QTRFIT_SERIAL_CALCULATOR")
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
    
    def test_openmp_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QTRFIT_OMP_CALCULATOR")
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
    
    def test_theobald_serial_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QCP_SERIAL_CALCULATOR")
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
    
    def test_theobald_OpenMP_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QCP_OMP_CALCULATOR")
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix,8)
    
    @unittest.skipIf(not "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators(),"CUDA calculator not available")
    def test_theobald_cuda_matrix_generation(self):
        """
        Calculates the whole matrix.
        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QCP_CUDA_CALCULATOR")
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, self.expected_serial_matrix, 4)
    
    def test_coordinates_change(self):
        """
        Tests if 
        """
        number_of_coordsets = 5;
        number_of_atoms = 3239;
        number_of_CAs = 224;
        

#      load_vector(not_aligned_CA, "../../src/calculators/test/data/ligand_mini_CAs");
        not_aligned_CA = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_CAs"),(number_of_coordsets,number_of_CAs,3))
        
#       load_vector(aligned_coordinates, "data/ligand_mini_all_aligned");
        aligned_coordinates = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_all_aligned"),(number_of_coordsets,number_of_atoms,3))

        reference_copy = numpy.copy(not_aligned_CA[0])
        target_copy = numpy.copy(not_aligned_CA[4])
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(not_aligned_CA, "QTRFIT_OMP_CALCULATOR", False)
        self.assertAlmostEqual(calculator.pairwise(0,4), 0.3383320758562839, 12)
        
        # The RMSD result must be the same
        calculator.modify_coordinates = True
        rmsd, coords =  calculator.pairwise(0,4)
        self.assertAlmostEqual(rmsd, 0.3383320758562839, 12)

        # After the calculations, the target and reference coordinates have not changed ...
        numpy.testing.assert_almost_equal(not_aligned_CA[0],reference_copy,16)
        numpy.testing.assert_almost_equal(not_aligned_CA[4],target_copy,16)

        # .. but returned reference has been modified (centered)
        reference_coords = not_aligned_CA[0]
        reference_geom_center = reference_coords.mean(0)
        reference_at_origin = (reference_coords - reference_geom_center)    
        numpy.testing.assert_almost_equal(reference_at_origin, coords[0], 12)


    ################################################
    # PART NOT NECESSARY WITH REGRESSION TEST      #
    ################################################
#         # Returned target has also been modified (centered and superposed)
#         pdb = parsePDB('../../src/calculators/test/data/ligand_mini_CA_aligned_to_frame_0.pdb')
#         print pdb.numAtoms()
#         pdb.addCoordset(aligned_coordinates[4])
#         alphas = pdb.select("name CA")
#         pdb_target_coords= alphas.getCoordsets()[5]
#         pdb_target_geom_center = pdb_target_coords.mean(0)
#         pdb_target_at_origin = (pdb_target_coords - pdb_target_geom_center)
#         numpy.testing.assert_almost_equal(pdb_target_at_origin, coords[1], 2)
#         
#         
#         pdb = parsePDB('../../src/calculators/test/data/ligand_mini.pdb')
#         alphas = pdb.select("name CA")
#         pdb_target_coords= alphas.getCoordsets()[4]
#         numpy.testing.assert_almost_equal(pdb_target_coords, target_copy, 2)
#         
#         # Proof that the superposition was done
#         self.assertRaises( AssertionError, numpy.testing.assert_almost_equal, *( pdb_target_coords, coords[1], 2) )

        regression_coords1 = numpy.loadtxt("./data/ligand_mini_CA_4_aligned_to_0")
        numpy.testing.assert_almost_equal(regression_coords1, coords[1], 12)
    
    ########################
    # PORTS OF C TESTS     #
    ########################
    def test_superposition_with_different_fit_and_calc_coordsets(self):
        number_of_coordsets = 5;
        number_of_atoms = 3239;
        number_of_CAs = 224;
        

#      load_vector(not_aligned_CA, "../../src/calculators/test/data/ligand_mini_CAs");
        not_aligned_CA = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_CAs"),(number_of_coordsets,number_of_CAs,3))
        
#      load_vector(not_aligned_coordinates, "data/ligand_mini_all");
        not_aligned_coordinates = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_all"),(number_of_coordsets,number_of_atoms,3))
            
#     // RMSD of all atoms using CA for superposition
#     RMSDCalculator* calculator1 = RMSDCalculatorFactory::createCalculator(
#                                                             type,
#                                                             number_of_coordsets,
#                                                             number_of_CAs,
#                                                             &(not_aligned_CA[0]),
#                                                             // Setting calculation coordinates
#                                                             number_of_atoms,
#                                                             &(not_aligned_coordinates[0]));

        calculator1 = pyRMSD.RMSDCalculator.RMSDCalculator(not_aligned_CA, "QCP_OMP_CALCULATOR")
        calculator1.setCalculationCoordinates(not_aligned_coordinates)
        numpy.testing.assert_almost_equal(
                              [1.864003731005552, 2.076760850428891, 3.596135117728627, 2.182685209336899],
                              calculator1.oneVsFollowing(0),
                              12)
        
#     // RMSD of CA using CA for superposition (default behavior)
#     RMSDCalculator* calculator2 = RMSDCalculatorFactory::createCalculator(
#                                                                 type,
#                                                                 number_of_coordsets,
#                                                                 number_of_CAs,
#                                                                 &(not_aligned_CA[0]));

        calculator2 = pyRMSD.RMSDCalculator.RMSDCalculator(not_aligned_CA, "QCP_OMP_CALCULATOR")
        numpy.testing.assert_almost_equal(
                              [0.767947519172927, 0.8838644164683896, 0.4177715823462121, 0.3383320758562839],
                              calculator2.oneVsFollowing(0),
                              12)

#     // RMSD  of CA using CA for superposition (using the same selection and RMSD subsets)
#     RMSDCalculator* calculator3 = RMSDCalculatorFactory::createCalculator(
#                                                                 type,
#                                                                 number_of_coordsets,
#                                                                 number_of_CAs,
#                                                                 &(not_aligned_CA[0]),
#                                                                 // Setting calculation coordinates
#                                                                 number_of_CAs,
#                                                                 &(not_aligned_CA[0]));
        calculator3 = pyRMSD.RMSDCalculator.RMSDCalculator(not_aligned_CA, "QTRFIT_OMP_CALCULATOR")
        calculator3.setCalculationCoordinates(not_aligned_CA)
        numpy.testing.assert_almost_equal(
                              [0.767947519172927, 0.8838644164683896, 0.4177715823462121, 0.3383320758562839],
                              calculator3.oneVsFollowing(0),
                              12)
        
    def test_iterative_superposition_with_equal_calc_and_fit_sets(self):
        number_of_coordsets = 5
        number_of_atoms = 3239
        
        not_aligned_coordinates = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_all"),(number_of_coordsets,number_of_atoms,3))

        iterposed_coordinates = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_iterposed_all"),(number_of_coordsets,number_of_atoms,3))
        
#     RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
#             type,
#             number_of_coordsets,
#             number_of_atoms,
#             &(not_aligned_coordinates[0]));

        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(not_aligned_coordinates, "QTRFIT_OMP_CALCULATOR")
        
        calculator.iterativeSuperposition()
        
        numpy.testing.assert_almost_equal(
                              iterposed_coordinates,
                              not_aligned_coordinates,
                              6)

    def test_iterative_superposition_with_different_calc_and_fit_sets(self):

        number_of_coordsets = 5
        number_of_atoms = 3239
        number_of_CAs = 224
        
        not_aligned_CA = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_CAs"),(number_of_coordsets,number_of_CAs,3))
        
        not_iterposed_coordinates = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_all"),(number_of_coordsets,number_of_atoms,3))
        
        iterposed_coordinates = numpy.reshape(numpy.loadtxt("../../src/calculators/test/data/ligand_mini_iterposed_with_cas_all_atom"),(number_of_coordsets,number_of_atoms,3))
        
#     RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
#                 type,
#                 number_of_coordsets,
#                 number_of_CAs,
#                 &(not_aligned_CA[0]),
#                 // "Calculation" coordinates
#                 number_of_atoms,
#                 &(not_iterposed_coordinates[0]));
        
        # Iterposes CA atoms using all atoms (weird but useful as a test)
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(not_aligned_CA, "QTRFIT_OMP_CALCULATOR")
        
        calculator.setCalculationCoordinates(not_iterposed_coordinates)
        
        calculator.iterativeSuperposition()
        
        numpy.testing.assert_almost_equal(
                              iterposed_coordinates,
                              not_iterposed_coordinates,
                              6)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()