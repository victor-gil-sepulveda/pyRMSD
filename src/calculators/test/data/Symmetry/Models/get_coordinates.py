import sys
import prody 
import numpy
import os.path

pdb_data = prody.parsePDB(sys.argv[1])
pdb_trajectory = prody.PDBEnsemble("aligned_CA")
pdb_name = os.path.splitext(sys.argv[1])[0]

prot = pdb_data.select("name CA not resname CA")
pdb_trajectory.setAtoms(prot)
pdb_trajectory.addCoordset(prot.getCoordsets())
pdb_trajectory.setCoords(prot.getCoordsets()[0])
prody.writePDB(pdb_name+"_CA.pdb", pdb_trajectory)
with file(pdb_name+".CA.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)
        
pdb_trajectory = prody.PDBEnsemble("")
lig = pdb_data.select("resname BEN not name H1 H2 H3 H4 H5 H6 H7 HN1 HN2 ")
pdb_trajectory.setAtoms(lig)
pdb_trajectory.addCoordset(lig.getCoordsets())
pdb_trajectory.setCoords(lig.getCoordsets()[0])
prody.writePDB(pdb_name+"_ligand.pdb", pdb_trajectory)
with file(pdb_name+".ligand.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

