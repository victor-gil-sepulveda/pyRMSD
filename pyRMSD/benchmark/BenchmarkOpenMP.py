'''
Created on 14/11/2012

@author: victor
'''
import pyRMSD.RMSDCalculator
import time
import numpy
import sys
#With OpenMP and amber_5k.pdb  it took:  2.5554060936 [ 0.0 ]
if __name__ == '__main__':

    files = [
            "amber_5k.pdb",
            "amber_10k.pdb",
            "amber_15k.pdb",
            "amber_20k.pdb",
            "amber_25k.pdb",
            "amber_30k.pdb",
            "amber_35k.pdb"
             ]  

    for pdb_file in files:
        print "Reading file ", "data/"+pdb_file,"...",
        sys.stdout.flush()
        coordsets = numpy.load("data/%s.npy"%pdb_file.split(".")[0])
        print "OK"
        number_of_conformations = coordsets.shape[0]
        number_of_atoms = coordsets.shape[1]
        
        times = []
        for i in range(1):
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, "QCP_OMP_CALCULATOR")
            calculator.setNumberOfOpenMPThreads(6)
            t1 = time.time()
            rmsd = calculator.pairwiseRMSDMatrix()
            t2 = time.time()
            del rmsd
            times.append(t2-t1)
        print "With OpenMP and",pdb_file, " it took: ",numpy.mean(times),"[", numpy.std(times), "]"
        sys.stdout.flush()

