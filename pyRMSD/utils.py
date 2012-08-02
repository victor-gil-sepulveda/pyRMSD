from prody.ensemble import PDBEnsemble
from prody.proteins import parsePDB
import numpy

def getCoorsetsFromPDB(pdb_path):
    pdb = parsePDB(pdb_path)
    ca_atoms = pdb.select('calpha')
    coordsets = ca_atoms.getCoordsets()
    number_of_conformations = len(coordsets)
    number_of_atoms = len(ca_atoms)
    return coordsets, number_of_conformations, number_of_atoms

def flattenCoords(coordsets):    
    ca_coords_list = []
    for coords in coordsets:
        for coord_triplet in coords:
            ca_coords_list.append(coord_triplet[0]) 
            ca_coords_list.append(coord_triplet[1]) 
            ca_coords_list.append(coord_triplet[2]) 
    return numpy.array(ca_coords_list)

def coalesceCoords(coordsets, number_of_conformations, number_of_atoms):
    ca_coords_list = [0]*(number_of_conformations*number_of_atoms*3)
    confId = 0
    for coords in coordsets:
        atomId = 0
        for coord_triplet in coords:
            ca_coords_list[atomId*number_of_conformations*3+0*number_of_conformations+confId] = coord_triplet[0]
            coordinate = coord_triplet[1]
            ca_coords_list[atomId*number_of_conformations*3+1*number_of_conformations+confId] = coord_triplet[1]
            ca_coords_list[atomId*number_of_conformations*3+2*number_of_conformations+confId] = coord_triplet[2]
            atomId += 1
        confId += 1
    return numpy.array(ca_coords_list)
