'''
Created on 14/11/2012

@author: victor
'''
import pyRMSD.RMSDCalculator
import pyRMSD.utils.proteinReading
import time
import os
import bz2

if __name__ == '__main__':
    using_cuda = "THEOBALD_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators()
    
    if not using_cuda:
        print "Build it using --cuda."
        exit()
        
    ######################
    # BENCHMARKING
    ######################
    print "Loading file..."
    t1 = time.time()
    print "\tUncompressing..."
    open("tmp_amber_long.pdb","w").write(bz2.BZ2File("data/amber_long.pdb.tar.bz2").read())
    print "\tLoading..."
    coordsets,number_of_conformations,number_of_atoms = pyRMSD.utils.proteinReading.getCoordsetsFromPDB('tmp_amber_long.pdb')
    os.system("rm tmp_amber_long.pdb")
    print "\tDeleting temporary file"
    np_coords = pyRMSD.utils.proteinReading.flattenCoords(coordsets)
    t2 = time.time()
    print 'Loading took %0.3f s' % (t2-t1)
    
    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "THEOBALD_SERIAL_OMP_CALCULATOR")
    for number_of_threads in range(1,9):
        calculator.setNumberOfOpenMPThreads(number_of_threads)
        t1 = time.time()
        calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        print 'Calculating took %0.3f s' % (t2-t1)
        
    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "THEOBALD_CUDA_CALCULATOR")
    
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
        
    for number_of_threads in threads:
        for number_of_blocks in blocks:
            #calculator.setCUDAKernelThreadsPerBlock(number_of_threads, number_of_blocks)
            t1 = time.time()
            #calculator.pairwiseRMSDMatrix()
            t2 = time.time()
#            print 'Calculating took %0.3f s' % (t2-t1)
            print number_of_threads, number_of_blocks