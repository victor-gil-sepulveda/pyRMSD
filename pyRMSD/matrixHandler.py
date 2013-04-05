'''
Created on 19/09/2012

@author: victor
'''
import json
import numpy
import pyRMSD.RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix #@UnresolvedImport
from pyRMSD.utils.proteinReading import Reader

class MatrixHandler(object):
    
    def __init__(self, statistics_folder=None):
        self.distance_matrix = None
        self.statistics_folder = statistics_folder
        
    def getMatrix(self):
        return self.distance_matrix

    def createMatrix(self, pdb_coordsets, calculator = "QCP_OMP_CALCULATOR"):
        print "Calculating matrix..."
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(pdb_coordsets, calculator).pairwiseRMSDMatrix()
        self.distance_matrix = CondensedMatrix(rmsd)
        self.distance_matrix.recalculateStatistics()
        self.__save_statistics()
        print " Done\n"
        return self.distance_matrix
    
    def createMatrixReadingOnlyCAs(self, pdb_file, reader_type = "LITE_READER", calculator = "QCP_OMP_CALCULATOR"):
        reader = Reader(reader_type).readThisFile(pdb_file).gettingOnlyCAs()
        return self.__createMatrixWithReader(reader, calculator)
    
    def createMatrixWithReader(self, pdb_file, reader_type = "LITE_READER", calculator = "QCP_OMP_CALCULATOR"):
        reader = Reader(reader_type).readThisFile(pdb_file)
        return self.__createMatrixWithReader(reader, calculator)
    
    def __createMatrixWithReader(self, reader, calculator = "QCP_OMP_CALCULATOR"):
        pdb_coordsets = reader.read()
        print "Calculating matrix..."
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(pdb_coordsets, calculator).pairwiseRMSDMatrix()
        self.distance_matrix = CondensedMatrix(rmsd)
        self.distance_matrix.recalculateStatistics()
        MatrixHandler.save_statistics(self.statistics_folder, self.distance_matrix)
        print " Done\n"
        return self.distance_matrix
        
    def saveMatrix(self, matrix_file_without_extension):
        print "Writing matrix data (in "+matrix_file_without_extension+".bin) ..."
        MatrixHandler.save_matrix(matrix_file_without_extension, self.distance_matrix)
        print " Done\n"
    
    def loadMatrix(self, matrix_file_without_extension):
        print "Loading matrix data from "+matrix_file_without_extension+".bin ..."
        self.distance_matrix = MatrixHandler.load_matrix(matrix_file_without_extension)
        self.distance_matrix.recalculateStatistics()
        MatrixHandler.save_statistics(self.statistics_folder, self.distance_matrix)
        print " Done\n"
    
    @classmethod
    def save_matrix(cls, matrix_file_without_extension, distance_matrix):
        numpy.save(matrix_file_without_extension+".bin",distance_matrix.get_data())
    
    @classmethod
    def load_matrix(cls, matrix_file_without_extension):
        return CondensedMatrix(list(numpy.load(matrix_file_without_extension+".bin")))
    
    @classmethod
    def save_statistics(cls, statistics_folder , distance_matrix):
        if statistics_folder!=None:
            stats_dic = {}
            stats_dic["Minimum"] =  distance_matrix.calculateMax() 
            stats_dic["Maximum"] =  distance_matrix.calculateMin()
            stats_dic["Mean"] =     distance_matrix.calculateMean()
            stats_dic["Std. Dev."] =distance_matrix.calculateVariance()
            stats_dic["Skewness"] = distance_matrix.calculateSkewness()
            stats_dic["Kurtosis"] = distance_matrix.calculateKurtosis()
            open( statistics_folder+"/"+"statistics.json","w").write(json.dumps(stats_dic,indent=4, separators=(',', ': ')))
