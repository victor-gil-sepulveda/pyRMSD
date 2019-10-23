"""
Created on 14/11/2012

@author: victor
"""
from __future__ import print_function
import pyRMSD.RMSDCalculator
import time
import os
import bz2
from pyRMSD.utils.proteinReading import Reader

if __name__ == '__main__':
    using_cuda = "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators()
    
    if not using_cuda:
        print("Build it using --cuda.")
        exit()
        
    ######################
    # BENCHMARKING
    ######################
    # Reading coords
    print("Loading file...")
    t1 = time.time()
    print("\tUncompressing...")
    open("tmp_amber_long.pdb","w").write(bz2.BZ2File("data/amber_long.pdb.tar.bz2").read())
    print("\tLoading...")
    reader = Reader().readThisFile('tmp_amber_long.pdb').gettingOnlyCAs()
    coordsets = reader.read() 
    number_of_atoms = reader.numberOfAtoms
    number_of_conformations = reader.numberOfFrames
    os.system("rm tmp_amber_long.pdb")
    print("\tDeleting temporary file")
    t2 = time.time()
    print('Loading took %0.3f s' % (t2-t1))
    
    for number_of_threads in range(1,13):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType="QCP_OMP_CALCULATOR", 
                                                          fittingCoordsets=coordsets)
        calculator.setNumberOfOpenMPThreads(number_of_threads)
        t1 = time.time()
        rmsd = calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        del rmsd
        print('OpenMP calculation took %0.3fs with number of threads %d' % (t2-t1, number_of_threads))
         
    # Generate test parameters
    max_n_threads = 512
    max_n_blocks = 512 
    threads = []
    blocks = []
    tmp_thread = 2
    while tmp_thread<=max_n_threads:
        threads.append(tmp_thread)
        tmp_thread *= 2
    tmp_blocks = 8
    while tmp_blocks<=max_n_blocks:
        blocks.append(tmp_blocks)
        tmp_blocks *= 2
    
    # Do test
    for number_of_threads in threads:
        for number_of_blocks in blocks:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType="QCP_CUDA_CALCULATOR",
                                                              fittingCoordsets=coordsets)
            calculator.setCUDAKernelThreadsPerBlock(number_of_threads, number_of_blocks)
            t1 = time.time()
            rmsd = calculator.pairwiseRMSDMatrix()
            t2 = time.time()
            del rmsd
            print('Calculating took %0.3fs with CUDA and numthreads:%d and numblocks:%d' % (t2-t1, number_of_threads, number_of_blocks))

        #Best in Minotauro (NVIDIA M2090): 128, 64
        #Best with Quadro FX 580: 2, 16
