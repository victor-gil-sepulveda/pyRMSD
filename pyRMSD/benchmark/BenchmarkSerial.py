"""
Created on 14/11/2012

@author: victor
"""
from __future__ import print_function
import pyRMSD.RMSDCalculator
import time
import numpy
import sys
#With Serial and  amber_5k.pdb  it took:  7.08920001984 [ 0.0 ]
if __name__ == '__main__':

    files = [
            ("amber_5k.pdb",10),
            ("amber_10k.pdb",10),
            ("amber_15k.pdb",10),
            ("amber_20k.pdb",8),
            ("amber_25k.pdb",7),
            ("amber_30k.pdb",5),
            ("amber_35k.pdb",3)]
    
    for pdb_file, tries in files:
        print("Reading file ", "data/"+pdb_file,"...", end=' ')
        sys.stdout.flush()
        coordsets = numpy.load("data/%s.npy"%pdb_file.split(".")[0])
        print("OK")
        sys.stdout.flush()
        number_of_conformations = coordsets.shape[0]
        number_of_atoms = coordsets.shape[1]
        
        times = []
        for i in range(tries):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType="QCP_SERIAL_CALCULATOR", fittingCoordsets=coordsets)
            t1 = time.time()
            rmsd = calculator.pairwiseRMSDMatrix()
            t2 = time.time()
            del rmsd
            times.append(t2-t1)
        print("With Serial and ",pdb_file, " it took: ",numpy.mean(times),"[", numpy.std(times), "]")
        sys.stdout.flush()

#28
# 8.831276
