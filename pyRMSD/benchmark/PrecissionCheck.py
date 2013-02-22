'''
Created on 22/02/2013

@author: victor
'''

import pyRMSD.RMSDCalculator
import time
import bz2
from pyRMSD.utils.proteinReading import Reader
import numpy

if __name__ == '__main__':
    # Reading coords
    print "Loading file..."
    t1 = time.time()
    print "\tUncompressing..."
    open("tmp_amber_long.pdb","w").write(bz2.BZ2File("data/amber_long.pdb.tar.bz2").read())
    print "\tLoading..."
    reader = Reader("PRODY_READER").readThisFile('tmp_amber_long.pdb').gettingOnlyCAs()
    coordsets = reader.read() 
    number_of_atoms = reader.numberOfAtoms
    number_of_conformations = reader.numberOfFrames

    rmsds = {}
    #---------------#
    rmsds["KABSCH"] = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "KABSCH_PYTHON_CALCULATOR").oneVsFollowing(0)
    #---------------# 
    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "QTRFIT_OMP_CALCULATOR")
    calculator.setNumberOfOpenMPThreads(4)
    rmsds["QTRFIT"] = calculator.oneVsFollowing(0)
    #---------------#
    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "QCP_OMP_CALCULATOR")
    calculator.setNumberOfOpenMPThreads(4)
    rmsds["QCP"] = calculator.oneVsFollowing(0)
    #---------------#
    mean_diffs = {}
    std_diffs = {}
    calculator_names = rmsds.keys()
    for calculator_name_i in calculator_names:
        print "i ",calculator_name_i
        mean_diffs[calculator_name_i] = {}
        std_diffs[calculator_name_i] = {}
        for calculator_name_j in calculator_names:
            print "\tj ",calculator_name_j
            rmsd_diff = abs(rmsds[calculator_name_i] - rmsds[calculator_name_j]) 
            mean_diffs[calculator_name_i][calculator_name_j] = numpy.mean(rmsd_diff, dtype=numpy.float64)
            std_diffs[calculator_name_i][calculator_name_j] = numpy.std(rmsd_diff, dtype=numpy.float64)
    #---------------#
    
    handler = open("mean_and_std","w")
    handler.write("Mean\n")
    for calculator_name in calculator_names:
        handler.write("%s "%calculator_name)
        handler.write("\n")
    for calculator_name_i in calculator_names:
        handler.write("%s "%calculator_name_i)
        for calculator_name_j in calculator_names:
            handler.write("%e "%mean_diffs[calculator_name_i][calculator_name_j])
        handler.write("\n") 
    handler.write("Std.\n")
    for calculator_name in calculator_names:
        handler.write("%s "%calculator_name)
        handler.write("\n")
    for calculator_name_i in calculator_names:
        handler.write("%s "%calculator_name_i)
        for calculator_name_j in calculator_names:
            handler.write("%e "%std_diffs[calculator_name_i][calculator_name_j])
        handler.write("\n") 
    handler.close()    