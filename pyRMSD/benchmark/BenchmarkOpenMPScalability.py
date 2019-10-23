"""
Created on 14/11/2012

@author: victor
"""
from __future__ import print_function
import pyRMSD.RMSDCalculator
import time
import numpy
import sys
#With OpenMP and amber_5k.pdb  it took:  2.5554060936 [ 0.0 ]
if __name__ == '__main__':

    coordsets = numpy.load("data/amber_30k.npy")
    number_of_conformations = coordsets.shape[0]
    number_of_atoms = coordsets.shape[1]
    print("Coordsets read")
    sys.stdout.flush()
    calculator_types = ["KABSCH_OMP_CALCULATOR","QTRFIT_OMP_CALCULATOR","QCP_OMP_CALCULATOR"]
    for calculator_type in calculator_types:
        for num_threads in range(1,7):
            times = []
            for i in range(5):
                calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType=calculator_type,
                                                                  fittingCoordsets=coordsets)
                calculator.setNumberOfOpenMPThreads(6)
                t1 = time.time()
                rmsd = calculator.pairwiseRMSDMatrix()
                t2 = time.time()
                del rmsd
                times.append(t2-t1)
            print(calculator_type, num_threads, numpy.mean(times),"[", numpy.std(times), "]")
            sys.stdout.flush()

