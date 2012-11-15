'''
Created on 14/11/2012

@author: victor
'''
from math import fabs

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
def checkRMSDs(rmsd1,rmsd2,precission = 1e-8):
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