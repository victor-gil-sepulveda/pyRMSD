'''
Created on 14/11/2012

@author: victor
'''
import collections
import numpy.testing
import pyRMSD.RMSDCalculator
import pyRMSD.utils.proteinReading
import time
import os
import bz2
import pyRMSD.pdbReader
if __name__ == '__main__':
    using_cuda = "THEOBALD_CUDA_CALCULATOR" in pyRMSD.RMSDCalculator.availableCalculators()
    
    ######################
    # BENCHMARKING
    ######################
    print "Loading file..."
    t1 = time.time()
    print "\tUncompressing..."
    open("tmp_amber_long.pdb","w").write(bz2.BZ2File("data/amber_long.pdb.tar.bz2").read())
    print "\tLoading..."
    coordsets,number_of_conformations,number_of_atoms = pyRMSD.utils.proteinReading.getCoordsetsFromPDB('tmp_amber_long.pdb')
    pyRMSD.pdbReader.readPDB('tmp_amber_long.pdb'," CA ")
    os.system("rm tmp_amber_long.pdb")
    print "\tDeleting temporary file"
    np_coords = pyRMSD.utils.proteinReading.flattenCoords(coordsets)
    t2 = time.time()
    print 'Loading took %0.3f s' % (t2-t1)
    
    #######################
    # CALCULATION
    #######################
    types = pyRMSD.RMSDCalculator.availableCalculators()
    
    DEFAULT_PRECISSION = 1e-8
    precissions = collections.defaultdict(lambda:DEFAULT_PRECISSION)
    
    if using_cuda:
        types.append("THEOBALD_CUDA_CALCULATOR")
        precissions["THEOBALD_CUDA_CALCULATOR"] = 1e-4
    
    rmsds = {}
    times = {}
    golden = "OMP_CALCULATOR"
    
    for CALC_TYPE in types:
        print "Calculating RMSD with ", CALC_TYPE
        t1 = time.time()
        pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, CALC_TYPE)
        rmsds[CALC_TYPE] = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, CALC_TYPE).pairwiseRMSDMatrix()
        t2 = time.time()
        times[CALC_TYPE] = t2-t1
        print "\tRmsd array generated. ", len(rmsds[CALC_TYPE]), " elements."
        print '\tMatrix generation with %s took %0.3fs' %(CALC_TYPE, t2-t1)
    
    ######################
    # VALIDATION
    #####################
    for CALC_TYPE in types:
        if CALC_TYPE != golden:
            print "Comparing results (%s vs %s)..."%(golden, CALC_TYPE)
            t1 = time.time()
            numpy.testing.assert_array_almost_equal(rmsds[golden], rmsds[CALC_TYPE],precissions[CALC_TYPE])
            t2 = time.time()
            print 'Comparison took %0.3fs' % (t2-t1)
    
    for CALC_TYPE in types:
        print CALC_TYPE,": ",times[CALC_TYPE]