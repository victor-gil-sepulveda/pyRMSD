"""
Created on 19/09/2012

@author: victor
"""
from __future__ import print_function
import json
import numpy
import os.path
import pyRMSD.RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix #@UnresolvedImport
import math

class MatrixHandler(object):

    def __init__(self, statistics_folder=None):
        """
        Class constructor.

        @param statistics_folder: Path to the folder that will contain the statistics summary. If None,
        no statistics will be saved.
        """
        self.distance_matrix = None
        self.statistics_folder = statistics_folder

    def getMatrix(self):
        """
        @return: The inner condensed matrix.
        """
        return self.distance_matrix

    def createMatrix(self, pdb_coordsets, calculator = "QCP_OMP_CALCULATOR"):
        """
        Generates a condensed matrix given a set of coordinates (usually representing a trajectory) and
        a calculator definition.

        @param pdb_coordsets: Input coordinates in Prody format.
        @param calculator: One of the possible available calculator types.

        @return: A condensed matrix with the pairwise RMSD matrix.
        """
        print("Calculating matrix ...")
        rmsd = pyRMSD.RMSDCalculator.RMSDCalculator(calculator, pdb_coordsets).pairwiseRMSDMatrix()
        self.distance_matrix = CondensedMatrix(rmsd)
        MatrixHandler.save_statistics(self.statistics_folder,self.distance_matrix)
        print(" Done")
        return self.distance_matrix

    def saveMatrix(self, matrix_file_without_extension):
        """
        Saves the internal matrix to a file. It adds the extension '.bin' to that name.

        @param matrix_file_without_extension: Is the matrix file name without any
        extension ('.bin' will be added).
        """
        print("Writing matrix data (in "+matrix_file_without_extension+".bin) ...", end=' ')
        MatrixHandler.save_matrix(matrix_file_without_extension, self.distance_matrix)
        print(" Done")

    def loadMatrix(self, matrix_file_without_extension):
        """
        Loads a condensed matrix from disk.

        @param matrix_file_without_extension: Is the matrix file name without any
        extension ('.bin' will be added).
        """
        print("Loading matrix data from "+matrix_file_without_extension+".bin ...", end=' ')
        self.distance_matrix = MatrixHandler.load_matrix(matrix_file_without_extension)
        print(" Done")

    @classmethod
    def save_matrix(cls, matrix_file_without_extension, distance_matrix):
        """
        Saves the internal matrix to a file. It adds the extension '.npy' to that name.

        @param matrix_file_without_extension: Is the matrix file name without any
        extension ('.npy' will be added).
        @param distance_matrix: The matrix to be saved.
        """
        numpy.save(matrix_file_without_extension,distance_matrix.get_data())

    @classmethod
    def load_matrix(cls, matrix_file_without_extension):
        """
        Loads a condensed matrix from disk. <ATT> Without the cast to list, the matrix
        is not loaded properly.</ATT>

        @param matrix_file_without_extension: Is the matrix file name without any
        extension ('.bin' will be added).

        @return: The loaded condensed matrix.
        """
        data = numpy.load(matrix_file_without_extension+".npy").astype(numpy.float64, copy=False)
        return CondensedMatrix(list(data))

    @classmethod
    def save_statistics(cls, statistics_folder, distance_matrix):
        """
        Saves a file with some statistics of the internal matrix in json format.

        @param statistics_folder: Folder where the file 'statistics.json' will be stored.
        @param distance_matrix: The distance matrix from which the statistics are calculated.
        """
        if statistics_folder is not None:
            print("Calculating statistics ...")
            stats_dic = {}
            stats_dic["Minimum"] =  distance_matrix.calculateMin()
            stats_dic["Maximum"] =  distance_matrix.calculateMax()
            stats_dic["Mean"] =     distance_matrix.calculateMean()
            stats_dic["Std. Dev."] = math.sqrt(distance_matrix.calculateVariance())
            stats_dic["Skewness"] = distance_matrix.calculateSkewness()
            stats_dic["Kurtosis"] = distance_matrix.calculateKurtosis()
            open( os.path.join(statistics_folder,"statistics.json"),"w").write(json.dumps(stats_dic,indent=4, separators=(',', ': ')))
            return os.path.join(statistics_folder,"statistics.json")
        return None
