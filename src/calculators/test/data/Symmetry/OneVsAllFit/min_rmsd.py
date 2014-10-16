import os
import numpy

files = os.listdir(".") 
files_to_process = []
for filename in files:
    if ".rmsd" in filename and not "minimum" in filename:
        files_to_process.append(filename)

values = []

for a_file in files_to_process:
    values.append(numpy.loadtxt(a_file))

min_rmsds = []
for rmsds in numpy.array(values).T:
    min_rmsds.append(min(rmsds))

print min_rmsds
numpy.savetxt("minimum.rmsds",min_rmsds)
