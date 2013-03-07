'''
Created on 16/11/2012

@author: victor
'''
import unittest
import pyRMSD.utils.proteinReading
import pyRMSD.RMSDCalculator
import numpy

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
    

    def test_python_pairwise(self):
        """
        Calculates all matrix elements one by one with the pairwise operation.
        """
        expected_rmsd_data = [0.22677106513739653, 0.44598234794295144, 0.37817804816455303]
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets_mini, "KABSCH_PYTHON_CALCULATOR")
        rmsd = []
        for i in range(len(self.coordsets_mini)):
            for j in range(i+1, len(self.coordsets_mini)):
                rmsd.append(calculator.pairwise(i, j))
        numpy.testing.assert_array_almost_equal(rmsd, expected_rmsd_data,8)
        
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
    
    def test_one_vs_others_python(self):
        """
        Calculates the reference vs the others with the python functions.
        """
        expected = [0.88785042, 0.82416178, 1.02156795, 0.69628994, 1.04059251, 0.77859792, 0.74962628, 0.73856698, 0.70444404, 0.92168545]
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "KABSCH_PYTHON_CALCULATOR").oneVsTheOthers(3)
        numpy.testing.assert_array_almost_equal(rmsd, expected,8)
        
    def test_one_vs_others_serial_omp(self):
        """
        Calculates the reference vs the others with the OpenMP functions.
        """
        expected = [0.88785042, 0.82416178, 1.02156795, 0.69628994, 1.04059251, 0.77859792, 0.74962628, 0.73856698, 0.70444404, 0.92168545]
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets, "QTRFIT_OMP_CALCULATOR").oneVsTheOthers(3)
        numpy.testing.assert_array_almost_equal(rmsd, expected,8)
        
    def test_mini_python_serial(self):
        """
        Calculates the whole pairwise matrix.
        """
        expected_rmsd_data = [0.22677106513739653, 0.44598234794295144, 0.37817804816455303]
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(self.coordsets_mini, "KABSCH_PYTHON_CALCULATOR")
        rmsd = calculator.pairwiseRMSDMatrix()
        numpy.testing.assert_array_almost_equal(rmsd, expected_rmsd_data,4)
    
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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()