import prody
import numpy

pdb_data = prody.parsePDB("../Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.pdb")
pdb_trajectory = prody.PDBEnsemble("aligned_CA")

rmsd_total = []

prot = pdb_data.select("name CA not resname CA")
pdb_trajectory.setCoords(pdb_data.getCoordsets()[0])
pdb_trajectory.addCoordset(pdb_data.getCoordsets())
pdb_trajectory.setAtoms(prot)

for i in range(0, pdb_trajectory.numConfs()):
    pdb_trajectory.setCoords(pdb_data.getCoordsets()[i])
    pdb_trajectory.superpose()
    rmsds = pdb_trajectory.getRMSDs()
    print rmsds
    rmsd_total.extend(rmsds[i+1:pdb_trajectory.numConfs()])

print rmsd_total
numpy.savetxt("prot_plus_ligand_offset_very_different.CA.rmsd_matrix", rmsd_total)