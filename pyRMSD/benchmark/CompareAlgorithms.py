'''
Created on 14/11/2012

@author: victor
'''
import pyRMSD.RMSDCalculator
import time
from pyRMSD.utils.proteinReading import Reader
import numpy
import sys

if __name__ == '__main__':

	filename = "data/amber_30k.pdb" 
	print "Reading file ", filename
	sys.stdout.flush()
	reader = Reader().readThisFile(filename).gettingOnlyCAs()
	coordsets = reader.read()
	number_of_atoms = reader.numberOfAtoms
	number_of_conformations = reader.numberOfFrames

	for calculator_type in ["QCP_SERIAL_CALCULATOR","QTRFIT_SERIAL_CALCULATOR","KABSCH_SERIAL_CALCULATOR",
						"QCP_OMP_CALCULATOR","QTRFIT_OMP_CALCULATOR","KABSCH_OMP_CALCULATOR"]:
		times = []
		for i in range(20):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type)
			t1 = time.time()
			rmsd = calculator.oneVsFollowing(0)
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print "Using ",calculator_type, "with",filename, "and without modifying coordinates it took: ",numpy.mean(times),"[", numpy.std(times), "]"
	
		times = []
		for i in range(20):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type, modifyCoordinates=True)
			t1 = time.time()
			rmsd = calculator.oneVsFollowing(0)
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print "Using ",calculator_type, "with",filename, "and modifying coordinates it took: ",numpy.mean(times),"[", numpy.std(times), "]"
		times = []
		for i in range(20):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type, modifyCoordinates=True)
			t1 = time.time()
			rmsd = calculator.iterativeSuperposition()
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print "Iterative superposition using",calculator_type, "with",filename, "and modifying coordinates it took: ",numpy.mean(times),"[", numpy.std(times), "]"



