#ifndef _KERNEL_SERIAL_OMP_H_
#define _KERNEL_SERIAL_OMP_H_

#define floating_point_type double

namespace ThRMSDSerialOmpKernel{
	void calcRMSDOfOneVsOthers( floating_point_type* all_coordinates, int base_conformation_id, int other_conformations_starting_id,
											int number_of_conformations, int number_of_atoms, floating_point_type* rmsd);

}
#endif
