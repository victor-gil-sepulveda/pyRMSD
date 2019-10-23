"""
Created on 14/11/2012

@author: victor
"""
from __future__ import print_function
import pyRMSD.RMSDCalculator
import time
import numpy
import sys

#With CUDA and  amber_5k.pdb  it took:  18.7449800968
#With CUDA and  amber_5k.pdb  it took:  1.81559896469 [ 0.0 ]
if __name__ == '__main__':
    using_cuda = "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators()
    
    if not using_cuda:
        print("Build it using --cuda.")
        exit()
    
    files = ["amber_5k.pdb","amber_10k.pdb","amber_15k.pdb","amber_20k.pdb","amber_25k.pdb","amber_30k.pdb","amber_35k.pdb"]  
    for pdb_file in files:
        print("Reading file data/"+pdb_file,"...", end=' ')
        sys.stdout.flush()
        coordsets = numpy.load("data/%s.npy"%pdb_file.split(".")[0])
        print("OK")
        number_of_conformations = coordsets.shape[0]
        number_of_atoms = coordsets.shape[1]
        
        times = []
        for i in range(20):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType="QCP_CUDA_CALCULATOR", 
                                                              fittingCoordsets=coordsets)
            calculator.setCUDAKernelThreadsPerBlock(2, 16)
#             calculator.setCUDAKernelThreadsPerBlock(128, 64)
#             calculator.setCUDAKernelThreadsPerBlock(256, 256)
            t1 = time.time()
            rmsd = calculator.pairwiseRMSDMatrix()
            t2 = time.time()
            del rmsd
            times.append(t2-t1)
            print(t2-t1)
            sys.stdout.flush()
        print("With CUDA and ",pdb_file, " it took: ",numpy.mean(times),"[", numpy.std(times), "]")
        sys.stdout.flush()

