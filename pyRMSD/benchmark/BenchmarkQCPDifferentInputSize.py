"""
Created on 22/02/2013

@author: victor
"""
from __future__ import print_function
import pyRMSD.RMSDCalculator
import time
import numpy
import sys


def grow_coordsets(coordsets, times):
    new_coordsets = []
    for coordset in coordsets:
        new_coordsets.append(coordset*times)
    return numpy.array(new_coordsets)

def add_coordsets_copy(coordsets,original_size):
    new_coordsets = []
    for coordset in coordsets:
        new_coordsets.append(numpy.append(coordset,coordset[:original_size], axis = 0))
    return numpy.array(new_coordsets)
if __name__ == '__main__':
    
    
    #------------------
    # CUDA
    #------------------
    #Best in Minotauro (NVIDIA M2090): 128, 64
    #Best with Quadro FX 580: 2, 16
    coordsets = numpy.load("data/tmp_amber_long.npy")
    number_of_atoms = coordsets.shape[1]
    number_of_conformations = coordsets.shape[0]
    print("Coordinates read (%d models, %d atoms)"%(number_of_conformations, number_of_atoms))
    sys.stdout.flush()
    original_size = coordsets.shape[1]
    for times in range(10):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType="QCP_CUDA_CALCULATOR", fittingCoordsets=coordsets)
        calculator.setCUDAKernelThreadsPerBlock(2, 16)
        t1 = time.time()
        rmsd = calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        del rmsd
        print("With CUDA and num. atoms %d it took: %fs"%(coordsets.shape[1],t2-t1))
        sys.stdout.flush()
        coordsets = add_coordsets_copy(coordsets, original_size)
    #------------------
    # OpenMP
    #------------------
    #Best in Minotauro (NVIDIA M2090): 6 threads
    #Best with Quadro FX 580: 4 threads
    coordsets = numpy.load("data/tmp_amber_long.npy")
    number_of_atoms = coordsets.shape[1]
    number_of_conformations = coordsets.shape[0]
    original_size = coordsets.shape[1]
    print("Coordinates read (%d models, %d atoms)"%(number_of_conformations, number_of_atoms))
    sys.stdout.flush()
    for times in range(10):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator( calculatorType="QCP_OMP_CALCULATOR", fittingCoordsets=coordsets)
        calculator.setNumberOfOpenMPThreads(4)
        t1 = time.time()
        rmsd = calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        del rmsd
        print("With OpenMP and num. atoms %d it took: %fs"%(coordsets.shape[1],t2-t1))
        sys.stdout.flush()
        coordsets = add_coordsets_copy(coordsets, original_size)
