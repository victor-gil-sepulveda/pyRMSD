"""
Created on 14/11/2012

@author: victor
"""
from __future__ import print_function
import collections
import numpy.testing
import pyRMSD.RMSDCalculator
import time
import os
import bz2
from pyRMSD.utils.proteinReading import Reader
if __name__ == '__main__':
    using_cuda = "QCP_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators()
    
    ######################
    # BENCHMARKING
    ######################
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
    
    #######################
    # CALCULATION
    #######################
    types = pyRMSD.RMSDCalculator.availableCalculators()
    
    DEFAULT_PRECISSION = 1e-8
    precissions = collections.defaultdict(lambda:DEFAULT_PRECISSION)
    
    if using_cuda:
        precissions["QCP_CUDA_CALCULATOR"] = 1e-4
    
    rmsds = {}
    times = {}
    golden = "QTRFIT_OMP_CALCULATOR"
    
    for CALC_TYPE in types:
        print("Calculating RMSD with ", CALC_TYPE)
        t1 = time.time()
        rmsds[CALC_TYPE] = pyRMSD.RMSDCalculator.RMSDCalculator(calculatorType=CALC_TYPE, 
                                                                fittingCoordsets=coordsets).pairwiseRMSDMatrix()
        t2 = time.time()
        times[CALC_TYPE] = t2-t1
        print("\tRmsd array generated. ", len(rmsds[CALC_TYPE]), " elements.")
        print('\tMatrix generation with %s took %0.3fs' %(CALC_TYPE, t2-t1))
    
    ######################
    # VALIDATION
    #####################
    for CALC_TYPE in types:
        if CALC_TYPE != golden:
            print("Comparing results (%s vs %s)..."%(golden, CALC_TYPE))
            t1 = time.time()
            numpy.testing.assert_array_almost_equal(rmsds[golden], rmsds[CALC_TYPE],precissions[CALC_TYPE])
            t2 = time.time()
            print('Comparison took %0.3fs' % (t2-t1))
    
    for CALC_TYPE in types:
        print(CALC_TYPE,": ",times[CALC_TYPE])
