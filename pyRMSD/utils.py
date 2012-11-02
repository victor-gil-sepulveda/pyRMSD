from prody.ensemble import PDBEnsemble
from prody.proteins import parsePDB
import numpy

def getCoordsetsFromPDB(pdb_path,selection_string=""):
    pdb = parsePDB(pdb_path, subset='calpha')
    print "PDB parsed (",pdb_path,")"
    if selection_string != "":
        selection = pdb.select(selection_string)
        pdb = selection.copy()
    coordsets = pdb.getCoordsets()
    number_of_conformations = len(coordsets)
    number_of_atoms = len(pdb)
    return coordsets, number_of_conformations, number_of_atoms

def getCoordsetsFromTwoPDBs(pdb_path1,pdb_path2,selection_string=""):
    coords1,num_conf1,num_atoms1 = getCoordsetsFromPDB(pdb_path1,selection_string)
    coords2,num_conf2,num_atoms2 = getCoordsetsFromPDB(pdb_path2,selection_string)
    complete_coordset =  numpy.append(coords1,coords2,axis=0)   
    if(num_atoms1 != num_atoms2):
        print "[WARNING] pdb1 and pdb2 have different number of CA atoms (",num_atoms1,",",num_atoms2,")"
    return complete_coordset, num_conf1+num_conf2, num_atoms1+num_atoms2

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
            ca_coords_list[atomId*number_of_conformations*3+1*number_of_conformations+confId] = coord_triplet[1]
            ca_coords_list[atomId*number_of_conformations*3+2*number_of_conformations+confId] = coord_triplet[2]
            atomId += 1
        confId += 1
    return numpy.array(ca_coords_list)
