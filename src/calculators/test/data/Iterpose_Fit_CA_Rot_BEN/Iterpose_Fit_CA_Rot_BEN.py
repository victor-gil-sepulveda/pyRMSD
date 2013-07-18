import prody
import numpy

pdb_data = prody.parsePDB("../Models/prot_stretching/stretching_trajectory_offset_ligand.pdb")
pdb_trajectory = prody.PDBEnsemble("iterposed_CA")

# Write the initial coordsets
all = pdb_data.select("all")
lig = pdb_data.select("resname BEN")
with file("stretching_trajectory_offset_ligand.initial_all.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%all.getCoordsets().shape)
    for coordset in all.getCoordsets():
        numpy.savetxt(outfile, coordset)
with file("stretching_trajectory_offset_ligand.initial_BEN.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%lig.getCoordsets().shape)
    for coordset in lig.getCoordsets():
        numpy.savetxt(outfile, coordset)

pdb_trajectory.setCoords(pdb_data.getCoordsets()[0])
pdb_trajectory.addCoordset(pdb_data.getCoordsets())
pdb_trajectory.setAtoms(all)
pdb_trajectory.iterpose()

prody.writePDB("stretching_trajectory_offset_ligand.iterposed_all.pdb", pdb_trajectory)
with file("stretching_trajectory_offset_ligand.iterposed_all.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

pdb_trajectory.setAtoms(lig) 
with file("stretching_trajectory_offset_ligand.iterposed_BEN.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)
        
rmsds =  pdb_trajectory.getRMSDs()

print rmsds

numpy.savetxt("stretching_trajectory_offset_ligand.iterposed_BEN.rmsd", rmsds)

