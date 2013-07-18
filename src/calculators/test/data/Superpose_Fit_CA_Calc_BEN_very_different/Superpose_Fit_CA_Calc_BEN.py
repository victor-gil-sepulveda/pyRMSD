import prody
import numpy

pdb_data = prody.parsePDB("../Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.pdb")
pdb_trajectory = prody.PDBEnsemble("aligned_CA")

prot = pdb_data.select("name CA not resname CA")
pdb_trajectory.setCoords(pdb_data.getCoordsets()[0])
pdb_trajectory.addCoordset(pdb_data.getCoordsets())
pdb_trajectory.setAtoms(prot)
pdb_trajectory.superpose()

prody.writePDB("prot_plus_ligand_similar.aligned_CA.pdb", pdb_trajectory)
with file("prot_plus_ligand_similar.aligned_CA.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

lig = pdb_data.select("resname BEN not element H")
pdb_trajectory.setAtoms(lig)
rmsds =  pdb_trajectory.getRMSDs()

prody.writePDB("prot_plus_ligand_similar.aligned_BEN.pdb", pdb_trajectory)
with file("prot_plus_ligand_similar.aligned_BEN.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

print rmsds

numpy.savetxt("prot_plus_ligand_similar.aligned_BEN.rmsd", rmsds)
