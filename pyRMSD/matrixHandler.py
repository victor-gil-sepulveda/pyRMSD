'''
Created on 19/09/2012

@author: victor
'''
import pyRMSD.RMSDCalculator
import pickle
from pyRMSD.condensedMatrix import CondensedMatrix #@UnresolvedImport
from pyRMSD.utils.proteinReading import Reader

class MatrixHandler(object):
    
    def __init__(self,statistics_folder=None):
        self.distance_matrix = None
        self.statistics_folder = statistics_folder

    def createMatrix(self, pdb_coordsets, calculator = "THEOBALD_SERIAL_OMP_CALCULATOR"):
        print "Calculating matrix..."
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(pdb_coordsets, calculator).pairwiseRMSDMatrix()
        self.distance_matrix = CondensedMatrix(rmsd)
        self.distance_matrix.recalculateStatistics()
        self.__save_statistics()
        print " Done\n"
    
    def createMatrixReadingOnlyCAs(self, pdb_file, reader_type = "LITE_READER", calculator = "THEOBALD_SERIAL_OMP_CALCULATOR"):
        reader = Reader(reader_type).readThisFile(pdb_file).gettingOnlyCAs()
        self.__createMatrixWithReader(reader, calculator)
    
    def createMatrixWithReader(self, pdb_file, reader_type = "LITE_READER", calculator = "THEOBALD_SERIAL_OMP_CALCULATOR"):
        reader = Reader(reader_type).readThisFile(pdb_file)
        self.__createMatrixWithReader(reader, calculator)
    
    def __createMatrixWithReader(self, reader, calculator = "THEOBALD_SERIAL_OMP_CALCULATOR"):
        pdb_coordsets = reader.read()
        print "Calculating matrix..."
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(pdb_coordsets, calculator).pairwiseRMSDMatrix()
        self.distance_matrix = CondensedMatrix(rmsd)
        self.distance_matrix.recalculateStatistics()
        self.__save_statistics()
        print " Done\n"
        
    def saveMatrix(self, matrix_file_without_extension):
        print "Writing matrix data (in "+matrix_file_without_extension+".bin) ..."
        output_handler = open(matrix_file_without_extension+".bin",'w')
        pickle.dump(self.distance_matrix.get_data(),output_handler)
        output_handler.close()
        print " Done\n"
    
    def loadMatrix(self, matrix_file_without_extension):
        print "Loading matrix data from "+matrix_file_without_extension+".bin ..."
        input_file_handler = open(matrix_file_without_extension+".bin","r")
        rmsd_data = pickle.load(input_file_handler)
        input_file_handler.close()
        self.distance_matrix = CondensedMatrix(list(rmsd_data))
        self.distance_matrix.recalculateStatistics()
        self.__save_statistics()
        print " Done\n"
        
    def __save_statistics(self):
        if self.statistics_folder!=None:
            print "SAVING"
            file_handler = open(self.statistics_folder+"/"+"statistics.txt","w")
            file_handler.write( "--------------\n")
            file_handler.write( "Matrix values\n") 
            file_handler.write( "-------------\n")
            file_handler.write( "Minimum:      %.5f\n"%(self.distance_matrix.calculateMax()))
            file_handler.write( "Maximum:      %.5f\n"%(self.distance_matrix.calculateMin()))
            file_handler.write( "Mean:         %.5f\n"%(self.distance_matrix.calculateMean()))
            file_handler.write( "Std. Dev.:    %.5f\n"%(self.distance_matrix.calculateVariance()))
            file_handler.write( "Skewness:     %.5f\n"%(self.distance_matrix.calculateSkewness()))
            file_handler.write( "Kurtosis:     %.5f\n"%(self.distance_matrix.calculateKurtosis()))
            file_handler.close()