import prody
import numpy
import sys
import os.path

pdb_data = prody.parsePDB(sys.argv[1])
pdb_name = os.path.basename(sys.argv[1])
pdb_trajectory = prody.PDBEnsemble("aligned_BEN")

lig = pdb_data.select("resname BEN not name H1 H2 H3 H4 H5 H6 H7 HN1 HN2")
pdb_trajectory.setCoords(pdb_data.getCoordsets()[0])
pdb_trajectory.addCoordset(pdb_data.getCoordsets())
pdb_trajectory.setAtoms(lig)
pdb_trajectory.superpose()
rmsds =  pdb_trajectory.getRMSDs()

prody.writePDB(pdb_name+".aligned_BEN.pdb", pdb_trajectory)
with file(pdb_name+".aligned_BEN.coords", 'w') as outfile:
    outfile.write("%d %d %d\n"%pdb_trajectory.getCoordsets().shape)
    for coordset in pdb_trajectory.getCoordsets():
        numpy.savetxt(outfile, coordset)

print rmsds

numpy.savetxt(pdb_name+".aligned_BEN.rmsd", rmsds)
