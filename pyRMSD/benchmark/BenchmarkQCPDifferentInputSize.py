'''
Created on 22/02/2013

@author: victor
'''
import pyRMSD.RMSDCalculator
import time
import bz2
from pyRMSD.utils.proteinReading import Reader
import numpy
import os

if __name__ == '__main__':
    # Reading coords
    print "Loading file..."
    t1 = time.time()
    print "\tUncompressing..."
    open("tmp_amber_long.pdb","w").write(bz2.BZ2File("data/amber_long.pdb.tar.bz2").read())
    print "\tLoading..."
    reader = Reader().readThisFile('tmp_amber_long.pdb').gettingOnlyCAs()
    coordsets = reader.read() 
    number_of_atoms = reader.numberOfAtoms
    
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
    #------------------
    # CUDA
    #------------------
    #Best in Minotauro (NVIDIA M2090): 128, 64
    #Best with Quadro FX 580: 2, 16
    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "QCP_CUDA_CALCULATOR")
    original_size = coordsets.shape[1]
    calculator.setCUDAKernelThreadsPerBlock(2, 16)
    t1 = time.time()
    rmsd = calculator.pairwiseRMSDMatrix()
    t2 = time.time()
    print "With CUDA and size %d it took: %fs"%(coordsets.shape[1]*3,t2-t1)
    
    for times in range(10):
        coordsets = add_coordsets_copy(coordsets, original_size)
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "QCP_CUDA_CALCULATOR")
        calculator.setCUDAKernelThreadsPerBlock(2, 16)
        t1 = time.time()
        rmsd = calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        del rmsd
        print "With CUDA and size %d it took: %fs"%(coordsets.shape[1]*3,t2-t1)

    #------------------
    # OpenMP
    #------------------
    #Best in Minotauro (NVIDIA M2090): 6 threads
    #Best with Quadro FX 580: 4 threads
    reader = Reader().readThisFile('tmp_amber_long.pdb').gettingOnlyCAs()
    coordsets = reader.read() 
    number_of_atoms = reader.numberOfAtoms
    os.system("rm tmp_amber_long.pdb")
    for times in range(10):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "QCP_OMP_CALCULATOR")
        calculator.setNumberOfOpenMPThreads(6)
        t1 = time.time()
        rmsd = calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        del rmsd
        print "With OpenMP and size %d it took: %fs"%(coordsets.shape[1]*3,t2-t1)
