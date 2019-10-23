"""
Created on 14/11/2012

@author: victor
"""
from __future__ import print_function
import pyRMSD.RMSDCalculator
import time
import numpy
import sys

if __name__ == '__main__':

	print("Reading file data/amber_30k.npy ...", end=' ')
	sys.stdout.flush()
	coordsets = numpy.load("data/amber_30k.npy")
	print("OK")
	number_of_conformations = coordsets.shape[0]
	number_of_atoms = coordsets.shape[1]

	for calculator_type in ["QCP_SERIAL_CALCULATOR","QTRFIT_SERIAL_CALCULATOR","KABSCH_SERIAL_CALCULATOR",
						"QCP_OMP_CALCULATOR","QTRFIT_OMP_CALCULATOR","KABSCH_OMP_CALCULATOR"]:
		times = []
		for i in range(10):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type)
			t1 = time.time()
			rmsd = calculator.oneVsFollowing(0)
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print("Using ",calculator_type, " without modifying coordinates it took: ",numpy.mean(times),"[", numpy.std(times), "]")
	
		times = []
		for i in range(10):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type, modifyCoordinates=True)
			t1 = time.time()
			rmsd = calculator.oneVsFollowing(0)
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print("Using ",calculator_type, " and modifying coordinates it took: ",numpy.mean(times),"[", numpy.std(times), "]")
		times = []
		for i in range(10):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type, modifyCoordinates=True)
			t1 = time.time()
			rmsd = calculator.iterativeSuperposition()
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print("Iterative superposition using",calculator_type, " and modifying coordinates it took: ",numpy.mean(times),"[", numpy.std(times), "]")
		for i in range(10):
			calculator = pyRMSD.RMSDCalculator.RMSDCalculator(coordsets, calculator_type, modifyCoordinates=True)
			t1 = time.time()
			rmsd = calculator.pairwiseRMSDMatrix()
			t2 = time.time()
			del rmsd
			times.append(t2-t1)
		print("Matrix using",calculator_type, " took", numpy.mean(times),"[", numpy.std(times), "]")



