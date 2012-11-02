from math import fabs
import time
import pyRMSD.RMSD as python_pure_functions
import pyRMSD.utils

#################
# Assertion
################
def rmsds_are_equal(rmsd1,rmsd2, precission):
    for i in range(max(len(rmsd1),len(rmsd2))):
        if fabs(rmsd1[i]- rmsd2[i]) > precission:
            print "Different RMSDs"
            print rmsd1
            print rmsd2
            return False;
    return True

#################
# Checking
################
def checkRMSDs(rmsd1,rmsd2,precission = 0.00001):
    if rmsds_are_equal(rmsd1, rmsd2, precission):
        print "OK"
    else:
        print "\t KK    KK     OOOOOO  "
        print "\t KK   KK     OO     O "
        print "\t KK  KK      OO     O "
        print "\t KK KK       OO     O "
        print "\t KKKK        OO     O "
        print "\t KKK         OO     O "
        print "\t KK K        OO     O "
        print "\t KK KK       OO     O "
        print "\t KK  KK      OO     O "
        print "\t KK   KK     OO     O "
        print "\t KK    KK    OO    OO "
        print "\t KK    KK     OOOOO   "
        exit()

#######################
# EXPECTED
######################
expected_rmsd = {0:[ 0.60179184,  0.70814575,  0.88785042,  0.92862096,  0.69024252,  0.59267699,  0.66596155,  0.81180133,  0.87438831,  1.00465129],
                 1:[ 0.61473279,  0.82416178,  0.96955624,  0.71842781,  0.5359385,   0.68621908,  0.90540226,  0.83185205,  0.96145774],
                 2:[ 1.02156795,  1.16059055,  0.80778577,  0.72752425,  0.80478222,  0.98594799,  1.04869932,  1.01149253],
                 3:[ 0.69628994,  1.04059251,  0.77859792,  0.74962628,  0.73856698,  0.70444404,  0.92168545]}

#######################
# MINI TEST
######################
#######################
# READING COORDINATES
######################
coordsets,number_of_conformations,number_of_atoms = pyRMSD.utils.getCoorsetsFromPDB('amber_mini.pdb')
np_coords = pyRMSD.utils.flattenCoords(coordsets)
########################
## TESTING PYTHON FUNCTIONS
#######################
expected_rmsd_data = [0.22677106513739653, 0.44598234794295144, 0.37817804816455303]
print "Testing python serial:"
rmsd =  python_pure_functions.__calculateRMSDCondensedMatrix(coordsets)
checkRMSDs(rmsd, expected_rmsd_data)
#######################
# TESTING CUDA
######################
#print "Testing CUDA serial:"
#rmsd =  python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"CUDA_CALCULATOR")
#checkRMSDs(rmsd, expected_rmsd_data, 0.0001)


#######################
# MEDIUM TEST, ONE VS OTHERS
######################
#######################
# READING COORDINATES
######################
coordsets,number_of_conformations,number_of_atoms = pyRMSD.utils.getCoorsetsFromPDB('amber_short.pdb')
np_coords = pyRMSD.utils.flattenCoords(coordsets)
#######################
# TEST SERIAL
######################
print "Testing serial:"
for conf_num in expected_rmsd:
    rmsd = python_pure_functions.oneVsTheOthers(conf_num,coordsets,"SERIAL_CALCULATOR")
    print "Conf. ",conf_num, " ", 
    checkRMSDs(rmsd, expected_rmsd[conf_num])
#######################
# TEST CUDA
######################
#print "Testing Cuda:"
#for conf_num in expected_rmsd:
#    rmsd = python_pure_functions.oneVsTheOthers(conf_num,coordsets,"CUDA_CALCULATOR")
#    print "Conf. ",conf_num, " ", 
#    checkRMSDs(rmsd, expected_rmsd[conf_num])
#######################
# TEST OPENMP
######################
print "Testing OpenMP:"
for conf_num in expected_rmsd:
    rmsd = python_pure_functions.oneVsTheOthers(conf_num,coordsets,"OMP_CALCULATOR")
    print "Conf. ",conf_num, " ", 
    checkRMSDs(rmsd, expected_rmsd[conf_num])

#######################
# MEDIUM TEST, WHOLE MATRIX
###################### 
expected_serial_matrix = [  0.60179184,0.70814575,0.88785042,0.92862096,0.69024252,0.59267699,
                            0.66596155,0.81180133,0.87438831,1.00465129,0.61473279,0.82416178,
                            0.96955624,0.71842781,0.5359385, 0.68621908,0.90540226,0.83185205,
                            0.96145774,1.02156795,1.16059055,0.80778577,0.72752425,0.80478222,
                            0.98594799,1.04869932,1.01149253,0.69628994,1.04059251,0.77859792,
                            0.74962628,0.73856698,0.70444404,0.92168545,1.08217543,0.86196576,
                            0.89731473,0.96848922,0.84721509,1.13748551,0.64892912,0.87248355,
                            1.00029474,1.01622641,1.10694473,0.68347196,0.83819283,0.7589582,
                            0.93694602,0.76944618,0.82288799,0.91196003,0.75938856,0.68278426,
                            0.76302383]

#######################
# TEST SERIAL MATRIX
######################
print "Testing Serial Matrix Generation:"
rmsd = python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"SERIAL_CALCULATOR")
checkRMSDs(rmsd, expected_serial_matrix)
#######################
# TEST CUDA MATRIX
######################
#print "Testing CUDA Matrix Generation:"
#rmsd = python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"CUDA_CALCULATOR")
#checkRMSDs(rmsd, expected_serial_matrix)
#######################
# TEST OMP MATRIX
######################
print "Testing OpenMP Matrix Generation:"
rmsd = python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"OMP_CALCULATOR")
checkRMSDs(rmsd, expected_serial_matrix)


#######################
# BENCHMARKING
######################
print "Loading file..."
t1 = time.time()
coordsets,number_of_conformations,number_of_atoms = pyRMSD.utils.getCoorsetsFromPDB('amber_long.pdb')
np_coords = pyRMSD.utils.flattenCoords(coordsets)
t2 = time.time()
print 'Loading took %0.3f s' % (t2-t1)

#######################
# CALCULATION
######################
t1 = time.time()
rmsd_cuda = python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"CUDA_CALCULATOR")
t2 = time.time()
print "Rmsd array generated. ",len(rmsd_cuda), " elements."
print 'Matrix generation with CUDA took %0.3fs' % (t2-t1)

t1 = time.time()
rmsd_serial = python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"SERIAL_CALCULATOR")
t2 = time.time()
print "Rmsd array generated. ",len(rmsd_serial), " elements."
print 'Matrix generation with Serial took %0.3fs' % (t2-t1)

t1 = time.time()
rmsd_omp = python_pure_functions.calculateRMSDCondensedMatrix(coordsets,"OMP_CALCULATOR")
t2 = time.time()
print "Rmsd array generated. ",len(rmsd_serial), " elements."
print 'Matrix generation with OpenMP took %0.3fs' % (t2-t1)

## !! this one takes tooooo long to finish
#t1 = time.time()
#rmsd_py = func.calculateRMSDs(coordsets)
#t2 = time.time()
#print "Rmsd array generated. ",len(rmsd_py), " elements."
#print 'Matrix generation with Pure Python took %0.3fs' % (t2-t1)

#######################
# VALIDATION
######################

print "Comparing results (Serial vs CUDA)..."
t1 = time.time()
checkRMSDs(rmsd_serial, rmsd_cuda)
t2 = time.time()
print 'Comparison took %0.3fs' % (t2-t1)

print "Comparing results (Serial vs OpenMP)..."
t1 = time.time()
checkRMSDs(rmsd_serial, rmsd_omp)
t2 = time.time()
print 'Comparison took %0.3fs' % (t2-t1)

