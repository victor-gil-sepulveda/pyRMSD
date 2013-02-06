'''
Created on 19/09/2012

@author: victor
'''
import pyRMSD.RMSDCalculator
import pickle
from pyRMSD.condensedMatrix import CondensedMatrix #@UnresolvedImport
from pyRMSD.utils.proteinReading import Reader
import json

class MatrixHandler(object):
    
    def __init__(self,statistics_folder=None):
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
        self.__save_statistics()
        print " Done\n"
        return self.distance_matrix
        
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
            stats_dic = {}
            stats_dic["Minimum"] =  self.distance_matrix.calculateMax() 
            stats_dic["Maximum"] =  self.distance_matrix.calculateMin()
            stats_dic["Mean"] =     self.distance_matrix.calculateMean()
            stats_dic["Std. Dev."] =self.distance_matrix.calculateVariance()
            stats_dic["Skewness"] = self.distance_matrix.calculateSkewness()
            stats_dic["Kurtosis"] = self.distance_matrix.calculateKurtosis()
            open( self.statistics_folder+"/"+"statistics.json","w").write(json.dumps(stats_dic,indent=4, separators=(',', ': ')))
