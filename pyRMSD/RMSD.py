import numpy.linalg
import numpy
import pyRMSD_cfuncs
from utils import flattenCoords

def oneVsTheOthers(target,coordsets,calcType = "PYTHON_CALCULATOR"):
    calculators = {"PYTHON_CALCULATOR":-1,"SERIAL_CALCULATOR":0,"OMP_CALCULATOR":1,"CUDA_CALCULATOR":2}
    if not calcType in calculators:
        print "Calculator ",calcType, " is not a calculator."
        return []
    
    # TODO: Check for availability of the device ie. CUDA
    number_of_conformations = len(coordsets)
    number_of_atoms = len(coordsets[0])
    conf_num = target
    np_coords = flattenCoords(coordsets)
    
    if calcType == "PYTHON_CALCULATOR":
        rmsd = []
        target = coordsets[conf_num]
        __oneVsTheOthers(target,coordsets[conf_num+1:],rmsd)
        return rmsd
    else:
        return pyRMSD_cfuncs.oneVsTheOthers(calculators[calcType],np_coords, number_of_atoms, conf_num, number_of_conformations)
    
    return []

def calculateRMSDCondensedMatrix(coordsets,calcType = "PYTHON_CALCULATOR"):
    calculators = {"PYTHON_CALCULATOR":-1,"SERIAL_CALCULATOR":0,"OMP_CALCULATOR":1,"CUDA_CALCULATOR":2}
    if not calcType in calculators:
        print "Calculator ",calcType, " is not a calculator."
        return []
    
    # TODO: Check for availability of the device ie. CUDA
    number_of_conformations = len(coordsets)
    number_of_atoms = len(coordsets[0])
    np_coords = flattenCoords(coordsets)
    if calcType == "PYTHON_CALCULATOR":
        return __calculateRMSDCondensedMatrix(coordsets)
    else:
        return pyRMSD_cfuncs.calculateRMSDCondensedMatrix(calculators[calcType],np_coords, number_of_atoms, number_of_conformations)
    
    return []



def __calculateRMSDCondensedMatrix(coordsets):
    rmsd_data = []
    for i, coords in enumerate(coordsets[:-1]):
        __oneVsTheOthers(coords,coordsets[i+1:],rmsd_data)
    return rmsd_data

def __oneVsTheOthers(target,coordsets,rmsd_data):
    """
    Modified Prody's code for superposition. This code in python is very very slow anyway.
    """
    divByN = 1.0 / target.shape[0]
    tar_com = target.mean(0)
    tar_org = (target - tar_com)
    mob_org = numpy.zeros(tar_org.shape, dtype=coordsets.dtype)
    tar_org = tar_org.T

    for i, mob in enumerate(coordsets):    
        mob_com = mob.mean(0)        
        matrix = numpy.dot(tar_org, numpy.subtract(mob, mob_com, mob_org))
            
        u, s, v = numpy.linalg.svd(matrix)
        Id = numpy.array([ [1, 0, 0], [0, 1, 0], [0, 0, numpy.sign(numpy.linalg.det(matrix))] ])
        rotation = numpy.dot(v.T, numpy.dot(Id, u.T))

        coordsets[i] = numpy.dot(mob_org, rotation) 
        numpy.add(coordsets[i], tar_com, coordsets[i]) 
        rmsd_data.append(numpy.sqrt(((coordsets[i]-target) ** 2).sum() * divByN))