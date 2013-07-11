import prody
import numpy
import sys

offsets = numpy.array([[[12, 40, 2]]*3239,
           [[5, -5, 8]]*3239,
           [[-13, -14, 7]]*3239,
           [[13, 14, 23]]*3239,
           [[-8,-4,-2]]*3239,
           [[6,89,2]]*3239])

pdb_data = prody.parsePDB(sys.argv[1])
pdb_trajectory = prody.PDBEnsemble(sys.argv[1])
pdb_trajectory.setAtoms(pdb_data)
pdb_trajectory.addCoordset(pdb_data.getCoordsets()+offsets[0:len(pdb_data.getCoordsets())])
pdb_trajectory.setCoords(pdb_data.getCoordsets()[0])
prody.writePDB(sys.argv[2], pdb_trajectory)
