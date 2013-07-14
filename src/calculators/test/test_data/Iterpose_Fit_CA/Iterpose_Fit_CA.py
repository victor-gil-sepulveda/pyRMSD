import prody
import numpy

pdb_data = prody.parsePDB("../Models/prot_stretching/stretching_trajectory_offset_ligand.pdb")
pdb_trajectory = prody.PDBEnsemble("iterposed_CA")

# Write the initial coordsets
prot = pdb_data.select("name CA")
prody.writePDB("stretching_trajectory_offset_ligand.iterposed_all.pdb", prot)
with file("stretching_trajectory_offset_ligand.initial_CA.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%prot.getCoordsets().shape)
    for coordset in prot.getCoordsets():
        numpy.savetxt(outfile, coordset)

# We only want to work with CAs. If we use the 'all coordinates+atom selection" trick
# Prody will still use all coordinates for iterative superposition
pdb_trajectory.setCoords(prot.getCoordsets()[0])
pdb_trajectory.addCoordset(prot.getCoordsets())
pdb_trajectory.setAtoms(prot)
pdb_trajectory.iterpose()

prody.writePDB("stretching_trajectory_offset_ligand.iterposed_CA.pdb", pdb_trajectory)
with file("stretching_trajectory_offset_ligand.iterposed_CA.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

