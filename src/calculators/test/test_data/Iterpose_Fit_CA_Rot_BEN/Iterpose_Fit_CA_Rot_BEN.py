import prody
import numpy

pdb_data = prody.parsePDB("../Models/prot_stretching/stretching_trajectory_offset_ligand.pdb")
pdb_trajectory = prody.PDBEnsemble("iterposed_CA")

prot = pdb_data.select("name CA")
pdb_trajectory.setCoords(pdb_data.getCoordsets()[0])
pdb_trajectory.addCoordset(pdb_data.getCoordsets())
pdb_trajectory.setAtoms(prot)
pdb_trajectory.iterpose()

prody.writePDB("stretching_trajectory_offset_ligand.iterposed_CA.pdb", pdb_trajectory)
with file("stretching_trajectory_offset_ligand.iterposed_CA.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)


lig = pdb_data.select("resname BEN not element H")
pdb_trajectory.setAtoms(lig)
rmsds =  pdb_trajectory.getRMSDs()

prody.writePDB("stretching_trajectory_offset_ligand.rot_BEN.pdb", pdb_trajectory)
with file("stretching_trajectory_offset_ligand.rot_BEN.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

print rmsds

numpy.savetxt("prot_plus_ligand_similar.rot_BEN.rmsd", rmsds)
