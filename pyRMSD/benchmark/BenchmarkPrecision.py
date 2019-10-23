"""
Created on 22/02/2013

@author: victor
"""

from __future__ import print_function
import pyRMSD.RMSDCalculator
import numpy
import math

NUMBER_OF_THREADS = 4

if __name__ == '__main__':
    # Reading coords

    coordsets = numpy.load("data/tmp_amber_long.npy")
    number_of_atoms = coordsets.shape[1]
    number_of_conformations = coordsets.shape[0]
    
    
    rmsds = {}
    calculators = {"KABSCH":'KABSCH_OMP_CALCULATOR',
                   "QTRFIT":"QTRFIT_OMP_CALCULATOR",
                   "QCP":"QCP_OMP_CALCULATOR",
                   "QCP CUDA":"QCP_CUDA_CALCULATOR"}
    
    for calculator_key in calculators:
        tmp_coords = numpy.copy(coordsets)
        calculator_type = calculators[calculator_key]
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(tmp_coords, calculator_type)
        if "OMP" in calculator_type:
            calculator.setNumberOfOpenMPThreads(NUMBER_OF_THREADS)
        if "CUDA" in calculator_type:
            calculator.setCUDAKernelThreadsPerBlock(2,16)
        rmsds[calculator_key] = calculator.oneVsFollowing(0)
    
    #---------------#
    rms = {}
    calculator_names = rmsds.keys()
    for calculator_name_i in calculator_names:
        print("* ",calculator_name_i)
        rms[calculator_name_i] = {}
        for calculator_name_j in calculator_names:
            print("\t ",calculator_name_j)
            rmsd_diff = rmsds[calculator_name_i] - rmsds[calculator_name_j]
            rms[calculator_name_i][calculator_name_j] = math.sqrt(numpy.sum(rmsd_diff**2))
    #---------------#
    
    handler = open("root_mean_square","w")
    for calculator_name in calculator_names:
        handler.write("%s "%calculator_name)
    handler.write("\n")
    for calculator_name_i in calculator_names:
        handler.write("%s "%calculator_name_i)
        for calculator_name_j in calculator_names:
            handler.write("%.03e "%(rms[calculator_name_i][calculator_name_j]))
        handler.write("\n") 
