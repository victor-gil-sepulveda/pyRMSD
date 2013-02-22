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
    
    
    
    def grow_coordsets(coordsets, times):
        new_coordsets = []
        for coordset in coordsets:
            new_coordsets.append(coordset*times)
        return numpy.array(new_coordsets)
    
    for times in [1,2,3,4,5,6]:
        new_coordsets = grow_coordsets(coordsets,times)
        new_number_of_atoms = len(new_coordsets[0])
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(new_coordsets, "QCP_CUDA_CALCULATOR")
        calculator.setCUDAKernelThreadsPerBlock(8, 8)
        t1 = time.time()
        rmsd = calculator.pairwiseRMSDMatrix()
        t2 = time.time()
        del rmsd
        print "With CUDA and size %d it took: %fs"%(new_number_of_atoms*3,t2-t1)
    
#     print grow_coordsets([
#                           [
#                            [1,2,3],
#                            [4,5,6],
#                            [7,8,9]
#                           ],
#                           [
#                            [10,11,12],
#                            [13,14,15],
#                            [16,17,18]]
#                           ],2)
