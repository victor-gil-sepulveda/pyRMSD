import prody
import numpy

pdb_data = prody.parsePDB("not_aligned_offset_prot_plus_ligand.pdb")
pdb_trajectory = prody.PDBEnsemble("aligned_CA")

prot = pdb_data.select("name CA not resname CA")
pdb_trajectory.setAtoms(prot)
pdb_trajectory.addCoordset(prot.getCoordsets())
pdb_trajectory.setCoords(prot.getCoordsets()[0])
pdb_trajectory.superpose()
rmsds =  pdb_trajectory.getRMSDs()

prody.writePDB("prot_plus_ligand.aligned_CA.pdb", pdb_trajectory)
with file("prot_plus_ligand.aligned_CA.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)
numpy.savetxt("prot_plus_ligand.aligned_CA.rmsd", rmsds)
print rmsds
